/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
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

#include "act/alexandria/topology.h"
#include "act/basics/msg_handler.h"
#include "act/coulombintegrals/gaussian_integrals.h"
#include "act/coulombintegrals/slater_integrals.h"
#include "act/forcefield/forcefield.h"
#include "act/utility/regression.h"
#include "act/utility/units.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/textreader.h"

namespace alexandria
{

void QgenResp::updateAtomCoords(const std::vector<gmx::RVec> &x)
{
    x_.resize(nAtom_);
    for (size_t i = 0; i < nAtom_; i++)
    {
        x_[i] = x[i];
    }
}

void QgenResp::updateAtomCharges(const std::vector<ActAtom> &atoms)
{
    if (nAtom_ != atoms.size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Inconsistency between number of resp atoms %zu and topology atoms %lu", nAtom_, atoms.size()).c_str()));
    }
    qshell_  = 0;
    for (size_t i = 0; i < nAtom_; i++)
    {
        q_[i] = atoms[i].charge();
        if (!mutable_[i])
        {
            q_[i]    = charge_[i]->value();
            qshell_ += q_[i];
        }
    }
}

void QgenResp::setAtomInfo(MsgHandler                   *msg_handler,
                           const std::vector<ActAtom>   &atoms,
                           alexandria::ForceField       *pd,
                           const int                     qtotal)
{
    if (nAtom_ != 0)
    {
        msg_handler->msg(ACTStatus::Warning,
                         "setAtomInfo called more than once. Ignoring second time.");
        return;
    }
    GMX_RELEASE_ASSERT(nAtom_ == 0 || (nAtom_ == atoms.size()),
                       gmx::formatString("Changing the number of atoms from %lu to %lu", nAtom_, atoms.size()).c_str());
    nAtom_   = atoms.size();
    qtot_    = qtotal;
    qshell_  = 0;
    auto qt       = pd->findForcesConst(InteractionType::ELECTROSTATICS);
    auto eqtModel = potentialToChargeDistributionType(qt.potential());
    bool haveZeta = eqtModel != ChargeDistributionType::Point;
    nFixed_ = 0;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        auto atype = pd->findParticleType(atoms[i].ffType());
        auto ztype = atype->interactionTypeToIdentifier(InteractionType::ELECTROSTATICS);
        auto qparm = atype->parameter("charge");
        q_.push_back(qparm->value());
        charge_.push_back(qparm);
        row_.push_back(atype->row());
        if (haveZeta)
        {
            zeta_.push_back(qt.findParameterTypeConst(ztype, "zeta").value());
        }
        else
        {
            zeta_.push_back(0.0);
        }
        
        mutable_.push_back(qparm->mutability() == Mutability::ACM);
        if (qparm->mutability() != Mutability::ACM)
        {
            nFixed_++;
            qshell_ += q_[i];
        }
    }
}

void QgenResp::summary(MsgHandler *msg_handler)
{
    if (msg_handler)
    {
        msg_handler->msg(ACTStatus::Info,
                         gmx::formatString("There are %zu atoms for (R)ESP fitting.", nAtom_));
        std::string msg;
        for (size_t i = 0; (i < nAtom_); i++)
        {
            msg += gmx::formatString(" %zu", symmetricAtoms_[i]);
        }
        msg_handler->msg(ACTStatus::Info, msg);
    }
}

void QgenResp::setAtomSymmetry(const std::vector<int> &symmetricAtoms)
{
    std::string msg = gmx::formatString("Please pass me a correct symmetric atoms vector. There are %zu symmetricAtoms and %zu RespAtoms",
                                        symmetricAtoms.size(),
                                        nAtom_);
    GMX_RELEASE_ASSERT(symmetricAtoms.size() == 0 ||
                       symmetricAtoms.size() == nAtom_,
                       msg.c_str());

    if (symmetricAtoms.size() == 0)
    {
        for (size_t i = 0; i < nAtom_; i++)
        {
            symmetricAtoms_.push_back(i);
        }
    }
    else
    {
        for(auto s : symmetricAtoms)
        {
            symmetricAtoms_.push_back(s);
        }
    }
    uniqueQ_ = 0;
    fitQ_    = 0;
    for (size_t i = 0; i < nAtom_; i++)
    {
        if (mutable_[i])
        {
            fitQ_ += 1;
            if (i - symmetricAtoms_[i] == 0)
            {
                uniqueQ_ += 1;
            }
        }
    }
}

void QgenResp::writeHisto(const std::string      &fn,
                          const std::string      &title,
                          const gmx_output_env_t *oenv)
{
    FILE       *fp;
    gmx_stats   gs;
    std::vector<double> x, y;
    int         nbin = 100;

    if (0 == fn.size())
    {
        return;
    }
    for (size_t i = 0; (i < nEsp()); i++)
    {
        gs.add_point(i, convertFromGromacs(ep_[i].vCalc(), "Hartree/e"), 0, 0);
    }

    gs.make_histogram(0, &nbin, eHisto::Y, 1, &x, &y);

    fp = xvgropen(fn.c_str(), title.c_str(), "Pot (1/a.u.)", "()", oenv);
    for (int i = 0; (i < nbin); i++)
    {
        fprintf(fp, "%10g  %10g\n", x[i], y[i]);
    }
    xvgrclose(fp);
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

void QgenResp::plotLsq(const gmx_output_env_t *oenv,
                       const char             *ESPcorr)
{
    real        x, y;
    const char *leg = "Alexandria";
    gmx_stats   lsq;
    std::string potUnit("Hartree/e");
    for (size_t i = 0; i < nEsp(); i++)
    {
        lsq.add_point(convertFromGromacs(ep_[i].v(), potUnit),
                      convertFromGromacs(ep_[i].vCalc(), potUnit),
                      0, 0);
    }
    FILE *fp = xvgropen(ESPcorr, "Electrostatic Potential (Hartree/e)", "QM", "Calc", oenv);
    xvgr_legend(fp, 1, &leg, oenv);
    xvgr_line_props(fp, 0, elNone, ecBlack, oenv);
    fprintf(fp, "@ s%d symbol %d\n", 0, 1);
    fprintf(fp, "@type xy\n");
    while (lsq.get_point(&x, &y, nullptr, nullptr, 0) == eStats::OK)
    {
        fprintf(fp, "%10g  %10g\n", x, y);
    }
    fprintf(fp, "&\n");
    xvgrclose(fp);
}

void QgenResp::regularizeCharges(MsgHandler *msg_handler)
{
    double qtot   = 0;
    for (size_t ii = 0; ii < nAtom_; ii++)
    {
        if (mutable_[ii])
        {
            qtot += q_[ii];
        }
    }
    double dq = (qtot_ - (qtot + qshell_))/(nAtom_-nFixed_);
    // If the algorithm works well, dq should be negligable.
    if (msg_handler && fabs(dq) > 0.001)
    {
        msg_handler->msg(ACTStatus::Warning,
                         gmx::formatString("Found qtot %g, should be %d. Subtracting %g from each atom.", qtot, qtot_, dq));
    }
    for (size_t ii = 0; ii < nAtom_; ii++)
    {
        if (mutable_[ii])
        {
            q_[ii] += dq;
        }
    }
}

void QgenResp::calcStatistics(MsgHandler *msg_handler)
{
    double pot2 = 0, sum2 = 0, ip = 0, calc2 = 0;
    mae_ = 0;
    mse_ = 0;
    for (size_t i = 0; (i < nEsp()); i++)
    {
        double vesp  = ep_[i].v();
        double vcalc = ep_[i].vCalc();
        double diff  = vesp - vcalc;
        if (msg_handler && msg_handler->debug() && (i < 4*nAtom_))
        {
            msg_handler->writeDebug(gmx::formatString("ESP %zu QM: %g FIT: %g DIFF: %g\n",
                                                      i, ep_[i].v(), ep_[i].vCalc(), diff));
        }
        mae_  += std::fabs(diff);
        mse_  += diff;
        sum2  += gmx::square(diff);
        pot2  += gmx::square(vesp);
        ip    += vesp*vcalc;
        calc2 += gmx::square(vcalc);
    }
    if (nEsp() > 0)
    {
        mae_ /= nEsp();
        mse_ /= nEsp();
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

real QgenResp::getStatistics(MsgHandler *msg_handler,
                             real       *rrms,
                             real       *cosangle,
                             real       *mae,
                             real       *mse)
{
    calcStatistics(msg_handler);
    *rrms     = rrms_;
    *cosangle = cosangle_;
    *mse      = mse_;
    *mae      = mae_;
    return rms_;
}

static double calcJ(ChargeDistributionType chargeType,
                    rvec                   espx,
                    rvec                   rax,
                    double                 zeta,
                    int                    row)
{
    rvec   dx;
    double r    = 0;
    double eTot = 0;

    rvec_sub(espx, rax, dx);
    r = norm(dx);
    if (zeta <= 0)
    {
        chargeType = ChargeDistributionType::Point;
    }
    if (chargeType == ChargeDistributionType::Point && r == 0)
    {
        gmx_fatal(FARGS, "Zero distance between the atom and the grid. Zeta = %g", zeta);
    }
    if (ChargeDistributionType::Gaussian == chargeType)
    {
        eTot = Nuclear_GG(r, zeta);
    }
    else if (ChargeDistributionType::Slater == chargeType)
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

void QgenResp::calcPot(MsgHandler *msg_handler,
                       double      epsilonr)
{
    double scale_factor = 1.0/epsilonr;
    
    // Loop over ESP points
    for (size_t i = 0; i < nEsp(); i++)
    {
        double vv    = 0;
        double qtot  = 0;
        auto   espx  = ep_[i].esp();
        
        // Loop over RESP atoms
        for (size_t j = 0; j < nAtom_; j++)
        {
            auto epot = calcJ(ChargeDistributionType_, espx, x_[j], zeta_[j],
                              row_[j]);
            qtot += q_[j];
            vv += (scale_factor*q_[j]*epot);
        }
        if (msg_handler && msg_handler->debug())
        {
            msg_handler->writeDebug(gmx::formatString("CalcESP[%lu]: vv = %8.5f qtot = %g Vqm = %8.5f\n",
                                                      i, vv, qtot, ep_[i].v()));
        }
        ep_[i].setVCalc(vv);
    }
}

void QgenResp::optimizeCharges(MsgHandler *msg_handler,
                               double      epsilonr)
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
    if (ncolumn == 0)
    {
        msg_handler->msg(ACTStatus::Warning, "No charges to fit with RESP. Check your input.");
        return;
    }
    MatrixWrapper         lhs(ncolumn, nrow);
    std::vector<double>   rhs;
    double                scale_factor = 1.0/std::sqrt(epsilonr);

    if (nEsp() < fitQ_)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("WARNING: Only %zu ESP points for %zu atoms. Cannot generate charges.", nEsp(), nAtom_).c_str()));
    }
    
    // Algorithm as described in Ghahremanpour et al., JCTC 14 (2018)
    // pp. 5553-5566.
    for (size_t j = 0; j < nEsp(); j++)
    {
        // Initialize the right  hand side with the QM potential
        rhs.push_back(ep_[j].v());
        // Get grid point coordinates
        auto espx  = ep_[j].esp();
        size_t i = 0;
        for (size_t ii = 0; ii < nAtom_; ii++)
        {
            // Compute potential due to this partice at grid
            // position j
            auto pot = scale_factor*calcJ(ChargeDistributionType_, espx, x_[ii],
                                          zeta_[ii], row_[ii]);
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
    size_t i = 0;
    for (size_t ii = 0; ii < nAtom_; ii++)
    {
        if (mutable_[ii])
        {
            lhs.set(i, nEsp(), factor);
            i++;
        }
    }
    rhs.push_back(factor * (qtot_ - qshell_)); // Add the total charge
    
    if (msg_handler && msg_handler->debug())
    {
        for (size_t j = 0; j < 4*nAtom_; j++)
        {
            auto espx = ep_[j].esp();
            msg_handler->writeDebug(gmx::formatString("ESP[%zu] espx = %g espy = %g espz = %g V= %g  rhs=%g\n",
                                                      j, espx[XX], espx[YY], espx[ZZ], ep_[j].v(), rhs[j]));
        }
    }

    // Add the equations to ascertain symmetric charges
    // We store the index of the mutable cores in ii1.
    std::map<int, int> ii1;
    int                i1   = 0;
    for (size_t i = 0; i < nAtom_; i++)
    {
        if (mutable_[i])
        {
            ii1.insert({i, i1});
            i1++;
        }
    }
    size_t j1 = rhs.size();
    for (size_t i = 0; i < nAtom_; i++)
    {
        if (symmetricAtoms_[i] < i)
        {
            lhs.set(ii1[i], j1, factor);
            lhs.set(ii1[symmetricAtoms_[i]], j1, -factor);
            rhs.push_back(0);
            j1++;
        }

    }
    GMX_RELEASE_ASSERT(j1 == rhs.size(), "Inconsistency adding equations for symmetric charges");
    GMX_RELEASE_ASSERT(j1 - nrow == 0, gmx::formatString("Something fishy adding equations for symmetric charges. j1 = %zu, nrow = %d", j1, nrow).c_str());
    if (msg_handler && msg_handler->debug())
    {
        msg_handler->writeDebug(gmx::formatString("ncolumn = %d nrow = %d point = %zu nfixed = %zu nUnique = %d\n",
                                                  ncolumn, nrow, nEsp(), fitQ_, uniqueQ_));
        for (int i = 0; i < nrow; i++)
        {
            std::string msg = gmx::formatString("ROW: %d", i);
            for (int j = 0; j < ncolumn; j++)
            {
                msg += gmx::formatString("  %8g", lhs.get(j, i));
            }
            msg += gmx::formatString("  RHS: %8g\n", rhs[i]);
            msg_handler->writeDebug(msg);
        }
        msg_handler->writeDebug(gmx::formatString("QCore in the r.h.s:%2g\n", rhs[nrow-1]));
        msg_handler->writeDebug(gmx::formatString("Qtot:%2d\n",   qtot_));
        msg_handler->writeDebug(gmx::formatString("QShell:%g\n", qshell_));
    }
    // Fit the charge
    std::vector<double> q;
    q.resize(ncolumn, -123.0);
    lhs.solve(rhs, &q);
    if (msg_handler && msg_handler->debug())
    {
        msg_handler->writeDebug("Fitted Charges from optimizeCharges:\n");
    }
    i = 0;
    for (size_t ii = 0; ii < nAtom_; ii++)
    {
        if (mutable_[ii])
        {
            q_[ii] = q[i];
            if (msg_handler && msg_handler->debug())
            {
                msg_handler->writeDebug(gmx::formatString("q[%zu] = %0.3f\n", i, q_[ii]));
            }
            i++;
        }
    }
    regularizeCharges(msg_handler);
}

void QgenResp::updateZeta(const std::vector<ActAtom> &atoms,
                          const ForceField              *pd)
{
    auto    fs   = pd->findForcesConst(InteractionType::ELECTROSTATICS);
    for (size_t i = 0; i < nAtom_; i++)
    {
        auto atype = pd->findParticleType(atoms[i].ffType());
        auto myid  = atype->interactionTypeToIdentifier(InteractionType::ELECTROSTATICS);
        auto eep   = fs.findParametersConst(myid);
        zeta_[i]   = eep.find("zeta")->second.value();
    }
}

} // namespace
