/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "qgen_resp.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <map>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/coulombintegrals/coulombintegrals.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/textreader.h"

#include "poldata.h"
#include "regression.h"
#include "units.h"

namespace alexandria
{

QgenResp::QgenResp(const QgenResp *src)
{
    ChargeType_         = src->ChargeType_;      
    watoms_             = src->watoms_;                  
    qtot_               = src->qtot_;                    
    qshell_             = src->qshell_;                  
    rms_                = src->rms_;                     
    rrms_               = src->rrms_;                    
    cosangle_           = src->cosangle_;                
    for(int m = 0; m < DIM; m++)
    {
        origin_[m] = src->origin_[m];  
        space_[m]  = src->space_[m];
        nxyz_[m]   = src->nxyz_[m];  
    }  
    uniqueQ_            = src->uniqueQ_;                 
    fitQ_               = src->fitQ_;                    
    nAtom_              = src->nAtom_;                   
    dzatoms_            = src->dzatoms_;                 
    stoichiometry_      = src->stoichiometry_;           
    ep_                 = src->ep_;                      
    symmetricAtoms_     = src->symmetricAtoms_;          
}

void QgenResp::updateAtomCoords(const gmx::HostVector<gmx::RVec> &x)
{
    for (int i = 0; i < nAtom_; i++)
    {
        x_[i] = x[i];
    }
}

void QgenResp::updateAtomCharges(t_atoms  *atoms)
{
    GMX_RELEASE_ASSERT(nAtom_ == atoms->nr,
                       "Inconsistency between number of resp atoms and topology atoms");

    for (int i = 0; i < nAtom_; i++)
    {
        q_[i] = atoms->atom[i].q;
    }
}

void QgenResp::setAtomInfo(t_atoms                          *atoms,
                           const alexandria::Poldata        *pd,
                           const gmx::HostVector<gmx::RVec> &x,
                           const int                         qtotal)
{
    nAtom_   = atoms->nr;
    qtot_    = qtotal;
    qshell_  = 0;
    x_       = x;
    auto zzz = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    nFixed_ = 0;
    for (int i = 0; i < atoms->nr; i++)
    {
        auto atype = pd->findAtype(*atoms->atomtype[i]);
        auto ztype = atype->id(InteractionType::CHARGEDISTRIBUTION);
        q_.push_back(atype->charge());
        row_.push_back(atype->row());
        zeta_.push_back(zzz.findParameterTypeConst(ztype, "zeta").value());
        
        mutable_.push_back(atype->mutability() != Mutability::Fixed);
        ptype_.push_back(atoms->atom[i].ptype);
        if (atype->mutability() == Mutability::Fixed)
        {
            nFixed_++;
            qshell_ += q_[i];
        }
    }
}

void QgenResp::summary(FILE *fp)
{
    if (nullptr != fp)
    {
        fprintf(fp, "There are %d atoms for (R)ESP fitting.\n", nAtom_);
        for (int i = 0; (i < nAtom_); i++)
        {
            fprintf(fp, " %d", symmetricAtoms_[i]);
        }
        fprintf(fp, "\n");
    }
}

void QgenResp::setAtomSymmetry(const std::vector<int> &symmetricAtoms)
{
    std::string msg = gmx::formatString("Please pass me a correct symmetric atoms vector. There are %d symmetricAtoms and %d RespAtoms",
                                        static_cast<int>(symmetricAtoms.size()),
                                        nAtom_);
    GMX_RELEASE_ASSERT(symmetricAtoms.size() == 0 ||
                       symmetricAtoms.size() == static_cast<size_t>(nAtom_),
                       msg.c_str());

    if (symmetricAtoms.size() == 0)
    {
        for (int i = 0; i < nAtom_; i++)
        {
            symmetricAtoms_.push_back(i);
        }
    }
    else
    {
        symmetricAtoms_ = symmetricAtoms;
    }
    uniqueQ_ = 0;
    fitQ_    = 0;
    for (int i = 0; i < nAtom_; i++)
    {
        if (mutable_[i])
        {
            fitQ_ += 1;
            if (symmetricAtoms_[i] == static_cast<int>(i))
            {
                uniqueQ_ += 1;
            }
        }
    }
    GMX_RELEASE_ASSERT(fitQ_ > 0, "No charges to fit");
}

void QgenResp::writeHisto(const std::string      &fn,
                          const std::string      &title,
                          const gmx_output_env_t *oenv)
{
    FILE       *fp;
    gmx_stats_t gs;
    real       *x, *y;
    int         nbin = 100;

    if (0 == fn.size())
    {
        return;
    }
    gs = gmx_stats_init();
    for (size_t i = 0; (i < nEsp()); i++)
    {
        gmx_stats_add_point(gs, i, convertFromGromacs(ep_[i].vCalc(), "Hartree/e"), 0, 0);
    }

    gmx_stats_make_histogram(gs, 0, &nbin, ehistoY, 1, &x, &y);

    fp = xvgropen(fn.c_str(), title.c_str(), "Pot (1/a.u.)", "()", oenv);
    for (int i = 0; (i < nbin); i++)
    {
        fprintf(fp, "%10g  %10g\n", x[i], y[i]);
    }
    free(x);
    free(y);
    fclose(fp);
    gmx_stats_free(gs);
}

void QgenResp::copyGrid(QgenResp &src)
{
    int m;

    for (m = 0; (m < DIM); m++)
    {
        origin_[m] = src.origin_[m];
        space_[m]  = src.space_[m];
        nxyz_[m]   = src.nxyz_[m];
    }
    int nesp = src.nEsp();
    ep_.clear();
    for (m = 0; (m < nesp); m++)
    {
        ep_.push_back(src.ep_[m]);
    }
}

void QgenResp::addEspPoint(double x, double y,
                           double z, double V)
{
    gmx::RVec rv(x, y, z);
    ep_.push_back(EspPoint(rv, V));
}

real QgenResp::myWeight(int iatom) const
{
    if (iatom < nAtom_)
    {
        return watoms_;
    }
    else
    {
        return 1.0;
    }
}

void QgenResp::plotLsq(const gmx_output_env_t *oenv,
                       const char             *ESPcorr)
{
    real        x, y;
    const char *leg = "Alexandria";
    gmx_stats_t lsq = gmx_stats_init();
    std::string potUnit("Hartree/e");
    for (size_t i = 0; i < nEsp(); i++)
    {
        gmx_stats_add_point(lsq,
                            convertFromGromacs(ep_[i].v(), potUnit),
                            convertFromGromacs(ep_[i].vCalc(), potUnit),
                            0, 0);
    }
    FILE *fp = xvgropen(ESPcorr, "Electrostatic Potential (Hartree/e)", "QM", "Calc", oenv);
    xvgr_legend(fp, 1, &leg, oenv);
    xvgr_line_props(fp, 0, elNone, ecBlack, oenv);
    fprintf(fp, "@ s%d symbol %d\n", 0, 1);
    fprintf(fp, "@type xy\n");
    while (gmx_stats_get_point(lsq, &x, &y, nullptr, nullptr, 0) == estatsOK)
    {
        fprintf(fp, "%10g  %10g\n", x, y);
    }
    fprintf(fp, "&\n");
    gmx_stats_free(lsq);
}

void QgenResp::regularizeCharges()
{
    double qtot   = 0;
    for (int ii = 0; ii < nAtom_; ii++)
    {
        if (mutable_[ii])
        {
            qtot += q_[ii];
        }
    }
    double dq = (qtot_ - (qtot + qshell_))/(nAtom_-nFixed_);
    if (debug)
    {
        fprintf(debug, "Found qtot %g, should be %d. Subtracting %g from each atom.\n",
                qtot, qtot_, dq);
    }
    for (int ii = 0; ii < nAtom_; ii++)
    {
        if (mutable_[ii])
        {
            q_[ii] += dq;
        }
    }
}

void QgenResp::calcRms()
{
    double pot2 = 0, sum2 = 0, ip = 0, calc2 = 0;
    for (size_t i = 0; (i < nEsp()); i++)
    {
        double vesp  = ep_[i].v();
        double vcalc = ep_[i].vCalc();
        double diff  = vesp - vcalc;
        if (debug && (i < static_cast<size_t>(4*nAtom_)))
        {
            fprintf(debug, "ESP %zu QM: %g FIT: %g DIFF: %g\n",
                    i, ep_[i].v(), ep_[i].vCalc(), diff);
        }
        sum2  += gmx::square(diff);
        pot2  += gmx::square(vesp);
        ip    += vesp*vcalc;
        calc2 += gmx::square(vcalc);
    }
    if (nEsp() > 0)
    {
        rms_     = convertFromGromacs(sqrt(sum2/nEsp()), "Hartree/e");
    }
    else
    {
        rms_     = 0;
    }
    rrms_     = sqrt(sum2/pot2);
    double denom = std::sqrt(pot2*calc2);
    if (denom > 0)
    {
        cosangle_ = ip/denom;
    }
    else
    {
        cosangle_ = 0.0;
    }
}

real QgenResp::getRms(real *rrms, real *cosangle)
{
    calcRms();
    *rrms     = rrms_;
    *cosangle = cosangle_;
    return rms_;
}

static double calcJ(ChargeType  chargeType,
                    rvec        espx,
                    rvec        rax,
                    double      zeta,
                    int         watoms,
                    int         row)
{
    rvec   dx;
    double r    = 0;
    double eTot = 0;

    rvec_sub(espx, rax, dx);
    r = norm(dx);
    if (zeta <= 0)
    {
        chargeType = ChargeType::Point;
    }
    if (watoms == 0 && r == 0)
    {
        gmx_fatal(FARGS, "Zero distance between the atom and the grid.");
    }
    if (ChargeType::Gaussian == chargeType)
    {
        eTot = Nuclear_GG(r, zeta);
    }
    else if (ChargeType::Slater == chargeType)
    {
        eTot = Nuclear_SS(r, row, zeta);
    }
    else
    {
        if (r == 0)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("r = 0 when computing Coulomb potential in %s, %d.", __FILE__, __LINE__)));
        }
        eTot = (1.0/r);
    }
    return ONE_4PI_EPS0*eTot;
}

void QgenResp::calcPot(double epsilonr)
{
    double scale_factor = 1.0/epsilonr;
    
    // Loop over ESP points
    for (size_t i = 0; i < nEsp(); i++)
    {
        double vv    = 0;
        double qtot  = 0;
        auto   espx  = ep_[i].esp();
        
        // Loop over RESP atoms
        for (int j = 0; j < nAtom_; j++)
        {
            auto epot = calcJ(ChargeType_, espx, x_[j], zeta_[j],
                              watoms_, row_[j]);
            qtot += q_[j];
            vv += (scale_factor*q_[j]*epot);
        }
        if (debug)
        {
            fprintf(debug, "CalcESP[%d]: vv = %8.5f qtot = %g Vqm = %8.5f\n",
                    static_cast<int>(i), vv, qtot, ep_[i].v());
        }
        ep_[i].setVCalc(vv);
    }
}

void QgenResp::optimizeCharges(double epsilonr)
{
    /*
       Increase number of rows for the symmetric atoms. E.g.
       if we know that atoms 2, 3 and 4 have the same charge we
       add two equation q2 - q3 = 0 and q2 - q4 = 0.
       One extra row is needed to reproduce the total charge.
     */
    int                   nrow     = nEsp() + 1 + fitQ_ - uniqueQ_;
    int                   factor   = nEsp()*50;
    int                   ncolumn  = fitQ_;
    MatrixWrapper         lhs(ncolumn, nrow);
    std::vector<double>   rhs;
    double                scale_factor = 1.0/std::sqrt(epsilonr);
    
    GMX_RELEASE_ASSERT(nEsp() > static_cast<size_t>(fitQ_), gmx::formatString("WARNING: Only %zu ESP points for %d atoms. Cannot generate charges.", nEsp(), nAtom_).c_str());
    
    // Algorithm as described in Ghahremanpour et al., JCTC 14 (2018)
    // pp. 5553-5566.
    for (size_t j = 0; j < nEsp(); j++)
    {
        // Initialize the right  hand side with the QM potential
        rhs.push_back(ep_[j].v());
        // Get grid point coordinates
        auto espx  = ep_[j].esp();
        int i = 0;
        for (int ii = 0; ii < nAtom_; ii++)
        {
            // Compute potential due to this partice at grid
            // position j
            auto pot = scale_factor*calcJ(ChargeType_, espx, x_[ii],
                                          zeta_[ii], watoms_, row_[ii]);
            // If this is a mutable atom
            // TODO check with symmetric charges and virtual sites
            if (mutable_[ii])
            {
                GMX_RELEASE_ASSERT(i < fitQ_, "trying to change wrong atom");
                lhs.set(i, j, lhs.get(i, j) + pot);
                i++;
            }
            else
            {
                // Typically for a shell, reduce right hand side by the
                // product of the calcJ function and the charge on  the
                // shell. Note that this is not written explicitly in the
                // JCTC paper (the charge is not mentioned).
                rhs[j] -= pot*q_[ii];
            }
        }
    }
    // Now fix the total charge
    int i = 0;
    for (int ii = 0; ii < nAtom_; ii++)
    {
        if (mutable_[ii])
        {
            lhs.set(i, nEsp(), factor);
            i++;
        }
    }
    rhs.push_back(factor * (qtot_ - qshell_)); // Add the total charge
    
    if (debug)
    {
        for (int j = 0; j < 4*nAtom_; j++)
        {
            auto espx = ep_[j].esp();
            fprintf(debug, "ESP[%d] espx = %g espy = %g espz = %g V= %g  rhs=%g\n",
                    j, espx[XX], espx[YY], espx[ZZ], ep_[j].v(), rhs[j]);
        }
    }

    // Add the equations to ascertain symmetric charges
    // We store the index of the mutable cores in ii1.
    std::map<int, int> ii1;
    int                i1   = 0;
    for (int i = 0; i < nAtom_; i++)
    {
        if (mutable_[i])
        {
            ii1.insert({i, i1});
            i1++;
        }
    }
    int j1 = rhs.size();
    for (int i = 0; i < nAtom_; i++)
    {
        if (symmetricAtoms_[i] < i)
        {
            lhs.set(ii1[i], j1, factor);
            lhs.set(ii1[symmetricAtoms_[i]], j1, -factor);
            rhs.push_back(0);
            j1++;
        }

    }
    GMX_RELEASE_ASSERT(j1 == static_cast<int>(rhs.size()), "Inconsistency adding equations for symmetric charges");
    GMX_RELEASE_ASSERT(j1 == nrow, "Something fishy adding equations for symmetric charges");
    if (debug)
    {
        fprintf(debug, "ncolumn = %d nrow = %d point = %zu nfixed = %d nUnique = %d\n",
                ncolumn, nrow, nEsp(), fitQ_, uniqueQ_);
        for (int i = 0; i < nrow; i++)
        {
            fprintf(debug, "ROW: %d", i);
            for (int j = 0; j < ncolumn; j++)
            {
                fprintf(debug, "  %8g", lhs.get(j, i));
            }
            fprintf(debug, "  RHS: %8g\n", rhs[i]);
        }
        fprintf(debug, "QCore in the r.h.s:%2g\n", rhs[nrow-1]);
        fprintf(debug, "Qtot:%2d\n",   qtot_);
        fprintf(debug, "QShell:%g\n", qshell_);
    }
    // Fit the charge
    std::vector<double> q;
    q.resize(ncolumn, -123.0);
    lhs.solve(rhs, &q);
    if (debug)
    {
        fprintf(debug, "Fitted Charges from optimizeCharges:\n");
    }
    i = 0;
    for (int ii = 0; ii < nAtom_; ii++)
    {
        if (mutable_[ii])
        {
            q_[ii] = q[i];
            if (debug)
            {
                fprintf(debug, "q[%d] = %0.3f\n", i, q_[ii]);
            }
            i++;
        }
    }
    regularizeCharges();
}

void QgenResp::updateZeta(t_atoms *atoms, const Poldata *pd)
{
    auto    fs   = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    for (int i = 0; i < nAtom_; i++)
    {
        auto atype = pd->findAtype(*(atoms->atomtype[i]));
        auto myid  = atype->id(InteractionType::CHARGEDISTRIBUTION);
        auto eep   = fs.findParametersConst(myid);
        zeta_[i]   = eep.find("zeta")->second.value();
    }
}

} // namespace
