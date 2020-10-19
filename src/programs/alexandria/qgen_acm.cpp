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

#include "qgen_acm.h"

#include <cctype>

#include "gromacs/coulombintegrals/coulombintegrals.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"

#include "molprop.h"
#include "poldata.h"
#include "regression.h"


namespace alexandria
{

QgenAcm::QgenAcm(const Poldata *pd,
                 t_atoms       *atoms,
                 int            qtotal)
{
    int          j;
    bool         bSupport = true;

    ChargeType_     = pd->chargeType();
    bWarned_        = false;
    bAllocSave_     = false;
    bHaveShell_     = pd->polarizable();
    eQGEN_          = eQGEN_OK;
    chieq_          = 0;
    Jcs_            = 0;
    Jss_            = 0;
    rms_            = 0;
    natom_          = 0;
    qtotal_         = qtotal;
    
    for (int i = 0; i < atoms->nr; i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            coreIndex_.push_back(i);
            natom_++;
        }
        else if (atoms->atom[i].ptype == eptShell)
        {
            shellIndex_.push_back(i);
        }
    }
    chi0_.resize(natom_, 0);
    rhs_.resize(natom_ + 1, 0);
    id_.resize(natom_);
    atomnr_.resize(natom_, 0);
    row_.resize(natom_);
    Jcc_.resize(natom_ + 1);
    zeta_.resize(natom_);
    j00_.resize(natom_, 0);
    q0_.resize(natom_, 0);
    q_.resize(natom_ + 1);   
    x_.resize(atoms->nr);

    auto fs = pd->findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
    for (int i = j = 0; (i < atoms->nr) && bSupport; i++)
    {
        // TODO This code assuses that shells are to be kept fixed
        // without checking the Poldata.
        if (atoms->atom[i].ptype == eptAtom)
        {
            auto atype = pd->findAtype(*atoms->atomtype[i]);
            Jcc_[j].resize(natom_ + 1, 0);
            int atm = atoms->atom[i].atomnumber;
            if (atm == 0)
            {
                gmx_fatal(FARGS, "Don't know atomic number for %s %s",
                          *(atoms->resinfo[i].name),
                          *(atoms->atomname[j]));
            }
            id_[j] = atype->id(eitELECTRONEGATIVITYEQUALIZATION);
            if (!fs.parameterExists(id_[j]))
            {
                fprintf(stderr, "No charge distribution support for atom %s\n",
                        id_[j].id().c_str());
                bSupport = false;
            }
            if (bSupport)
            {
                auto eem   = fs.findParametersConst(id_[j]);
                atomnr_[j] = atm;
                q0_[j]     = eem.find("ref_charge")->second.value();
                // TODO Fix this business
                q_[j]    = eem.find("charge")->second.value();
                row_[j]  = eem.find("row")->second.value();
                    
                if (ChargeType_ == eqtSlater) 
                { 
                    if (row_[j] > SLATER_MAX)
                    {
                        if (debug)
                        {
                            fprintf(debug, "Cannot handle Slaters higher than %d for atom %s %s\n",
                                    SLATER_MAX,
                                    *(atoms->resinfo[i].name),
                                    *(atoms->atomname[j]));
                        }
                        row_[j] = SLATER_MAX;
                    }
                }
                j++;
            }
        }
    }
    // Call routine to set the Chi, J00 and Zeta.
    updateParameters(pd);
}

void QgenAcm::updateParameters(const Poldata *pd)
{
    auto fs = pd->findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
    for (auto i = 0; i < natom_; i++)
    {
        auto ei  = fs.findParametersConst(id_[i]);
        chi0_[i] = ei.find("chi")->second.value();
        j00_[i]  = ei.find("jaa")->second.value();
        zeta_[i] = ei.find("zeta")->second.value();
    }
}

double QgenAcm::getQ(int atom)
{
    if ((0 <= atom) && (atom < natom_))
    {
        return q_[atom];
    }
    return 0;
}

int QgenAcm::getRow( int atom)
{
    if ((0 <= atom) && (atom < natom_))
    {
        return row_[atom];
    }
    return 0;

}

double QgenAcm::getZeta(int atom)
{
    if ((0 <= atom) && (atom < natom_))
    {
        return zeta_[atom];
    }
    return 0;
}

static double Coulomb_PP(double r)
{
    return 1/r;
}

void QgenAcm::dump(FILE *fp, t_atoms *atoms)
{
    if (fp && eQGEN_ == eQGEN_OK)
    {
        auto  i  = 0, j = 0;
        rvec  mu = { 0, 0, 0 };

        fprintf(fp, "Jcc_ matrix:\n");
        for (i = 0; i < natom_; i++)
        {
            for (j = 0; j <= i; j++)
            {
                fprintf(fp, "  %6.2f", Jcc_[i][j]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");

        fprintf(fp, "RHS_ matrix:\n");
        for (i = 0; i < natom_; i++)
        {
            fprintf(fp, "%2d RHS_ = %10g\n", i, rhs_[i]);
        }
        fprintf(fp, "\n");

        fprintf(fp, "                                          Core                 Shell\n");
        fprintf(fp, "Res  Atom   Nr       J0    chi0  row   q   zeta        row    q     zeta\n");
        for (i = j = 0; i < atoms->nr; i++)
        {
            for (auto m = 0; m < DIM; m++)
            {
                mu[m] += atoms->atom[i].q * x_[i][m] * ENM2DEBYE;
            }
            if (atoms->atom[i].ptype == eptAtom)
            {
                fprintf(fp, "%4s %4s%5d %8g %8g",
                        *(atoms->resinfo[atoms->atom[i].resind].name),
                        *(atoms->atomname[i]), i+1, j00_[j], chi0_[j]);
                fprintf(fp, " %3d %8.5f %8.4f\n", row_[j], q_[j], zeta_[j]);
                j++;
            }
        }
        fprintf(fp, "\n");
        fprintf(fp, "chieq = %10g\n|mu| = %8.3f ( %8.3f  %8.3f  %8.3f )\n\n",
                chieq_, norm(mu), mu[XX], mu[YY], mu[ZZ]);
    }
}

const char *QgenAcm::message() const
{
    switch (eQGEN_)
    {
        case eQGEN_OK:
            return "Charge generation finished correctly";
        case eQGEN_NOTCONVERGED:
            return "Charge generation did not converge.";
        case eQGEN_NOSUPPORT:
            return "No charge generation support for (some of) the atomtypes.";
        case eQGEN_ERROR:
        default:
            return "Unknown status %d in charge generation";
    }
    return nullptr;
}

double QgenAcm::calcSij(int i, int j)
{
    double dist, dism, Sij = 1.0;
    rvec   dx;
    int    l, m, tag;

    rvec_sub(x_[i], x_[j], dx);
    dist = norm(dx);
    if ((dist < 0.118) && (atomnr_[i] != 1) && (atomnr_[j] != 1))
    {
        Sij = Sij*1.64;
    }
    else if ((dist < 0.122) && (atomnr_[i] != 1) && (atomnr_[j] != 1))
    {
        if ((atomnr_[i] != 8) && (atomnr_[j] != 8))
        {
            Sij = Sij*2.23;
        }
        else
        {
            Sij = Sij*2.23;
        }
    }
    else if (dist < 0.125)
    {
        tag = 0;
        if ((atomnr_[i] == 6) && (atomnr_[j] == 8))
        {
            tag = i;
        }
        else if ((atomnr_[i] == 8) && (atomnr_[j] == 6))
        {
            tag = j;
        }
        if (tag != 0)
        {
            printf("found CO\n");
            for (l = 0; (l < natom_); l++)
            {
                if (atomnr_[l] == 1)
                {
                    printf("found H\n");
                    dism = 0.0;
                    for (m = 0; (m < DIM); m++)
                    {
                        dism = dism+gmx::square(x_[tag][m] - x_[l][m]);
                    }

                    printf("dist: %8.3f\n", sqrt(dism));
                    if (sqrt(dism) < 0.105)
                    {
                        printf("dist %5d %5d %5s  %5s %8.3f\n",
                               i, l, id_[tag].id().c_str(), id_[l].id().c_str(), sqrt(dism));
                        Sij = Sij*1.605;
                    }
                }
            }
        }
    }
    else if ((atomnr_[i] == 6) && (atomnr_[j] == 8))
    {
        Sij = Sij*1.03;
    }
    else if (((atomnr_[j] == 6) && (atomnr_[i] == 7) && (dist < 0.15)) ||
             ((atomnr_[i] == 6) && (atomnr_[j] == 7) && (dist < 0.15)))
    {
        if (atomnr_[i] == 6)
        {
            tag = i;
        }
        else
        {
            tag = j;
        }
        for (l = 0; (l < natom_); l++)
        {
            if (atomnr_[l] == 8)
            {
                printf("found Oxy\n");
                dism = 0.0;
                for (m = 0; (m < DIM); m++)
                {
                    dism = dism+gmx::square(x_[tag][m] - x_[l][m]);
                }
                if (sqrt(dism) < 0.130)
                {
                    printf("found peptide bond\n");
                    Sij = Sij*0.66;
                }
                else
                {
                    Sij = Sij*1.1;
                }
            }
        }
    }
    return Sij;
}

double QgenAcm::calcJ(rvec    xI,
                      rvec    xJ,
                      double  zetaI,
                      double  zetaJ,
                      int     rowI,
                      int     rowJ,
                      double  epsilonr)
{
    rvec   dx;
    double r    = 0;
    double eTot = 0;
    rvec_sub(xI, xJ, dx);
    r = norm(dx);
    
    if (r == 0 && ChargeType_ == eqtPoint)
    {
        gmx_fatal(FARGS, "Zero distance between atoms!\n");
    }
    switch (ChargeType_)
    {
    case eqtPoint:
        {
            eTot = Coulomb_PP(r);
            break;
        }
    case eqtGaussian:
        {
            eTot = Coulomb_GG(r, zetaI, zetaJ);
            break;
        }
    case eqtSlater:
        {
            eTot = Coulomb_SS(r, rowI, rowJ, zetaI, zetaJ);
            break;
        }
    default:
        {
            gmx_fatal(FARGS, "Unsupported model %s in calc_jcc",
                      chargeTypeName(ChargeType_).c_str());
        }
    }
    return (ONE_4PI_EPS0*eTot)/(epsilonr*ELECTRONVOLT);
}

void QgenAcm::calcJcc(t_atoms *atoms, double epsilonr,
                      bool bYang, bool bRappe)
{
    auto Jcc = 0.0;
    auto i   = 0;
    for (auto n = 0; n < atoms->nr; n++)
    {
        if (atoms->atom[n].ptype == eptAtom)
        {
            auto j = i;
            for (auto m = n; m < atoms->nr; m++)
            {
                if (atoms->atom[m].ptype == eptAtom)
                {
                    if (n != m)
                    {
                        j++;
                        Jcc = calcJ(x_[n],
                                    x_[m],
                                    zeta_[i],
                                    zeta_[j],
                                    row_[i],
                                    row_[j],
                                    epsilonr);
                        if (bYang)
                        {
                            Jcc *= calcSij(i, j);
                        }
                        Jcc_[i][j] = Jcc_[j][i] = (0.5 * Jcc);
                    }
                    else
                    {
                        auto j0 = j00_[i];
                        if ((bYang || bRappe) && (atomnr_[i] == 1))
                        {
                            auto zetaH = 1.0698;
                            j0 = (1+q_[i]/zetaH)*j0; /* Eqn. 21 in Rappe1991. */
                            if (j0 < 0 && !bWarned_)
                            {
                                bWarned_ = true;
                            }
                            Jcc_[i][i] = (j0 > 0) ? j0 : 0;
                        }
                        else
                        {
                            Jcc_[i][i] = (j0 > 0) ? j0 : 0;
                        }
                    }
                }
            }
            i++;
        }
    }
}

void QgenAcm::calcJcs(t_atoms *atoms,
                      int      core_ndx_gromacs,
                      int      core_ndx_eem,
                      double   epsilonr)
{
    Jcs_ = 0;
    if (atoms->atom[core_ndx_gromacs].ptype == eptAtom)
    {
        auto itsShell = core_ndx_gromacs + 1;
        for (auto i = 0; i < natom_; i++)
        {
            auto shell_ndx_gromacs = shellIndex_[i];
            if (atoms->atom[shell_ndx_gromacs].ptype == eptShell)
            {
                if (shell_ndx_gromacs != itsShell)
                {
                    Jcs_ += (q_[i]*calcJ(x_[core_ndx_gromacs],
                                         x_[shell_ndx_gromacs],
                                         zeta_[core_ndx_eem],
                                         zeta_[i],
                                         row_[core_ndx_eem],
                                         row_[i],
                                         epsilonr));  
                    if (debug)
                    {
                        fprintf(debug, "core_ndx: %d shell_ndx: %d shell_charge: %0.1f\n", 
                                core_ndx_gromacs, shell_ndx_gromacs, q_[i]);
                    }
                }
            }
            else
            {
                gmx_fatal(FARGS, "atom %d must be eptShell, but it is not\n", shell_ndx_gromacs);
            }
        }
        Jcs_ *= 0.5;
        if (debug)
        {
            fprintf(debug, "Jcs_:%0.3f\n", Jcs_);
        }
    }
    else
    {
        gmx_fatal(FARGS, "atom %d must be eptAtom, but it is not\n", core_ndx_gromacs);
    }
}

void QgenAcm::calcRhs(t_atoms *atoms, double epsilonr)
{
    auto   qcore     = 0.0;
    auto   qshell    = 0.0;

    for (auto i = 0; i < natom_; i++)
    {
        rhs_[i]   = 0;
        rhs_[i]  -= chi0_[i];                              // Electronegativity
        rhs_[i]  += j00_[i]*q0_[i];                        // Hardness * q0
        if (bHaveShell_)
        {
            calcJcs(atoms, coreIndex_[i], i, epsilonr);
            rhs_[i]   -= j00_[i]*q_[i]; // Hardness * qs
            rhs_[i]   -= Jcs_;                             // Core-Shell interaction
            qshell    += q_[i];
        }
    }
    qcore        = qtotal_ - qshell;
    rhs_[natom_] = qcore;
}

void QgenAcm::copyChargesToAtoms(t_atoms *atoms)
{
    if (bHaveShell_)
    {
        auto j = 0;
        for (auto i = j = 0; i < atoms->nr; i++)
        {
            if (atoms->atom[i].ptype == eptAtom)
            {
                atoms->atom[i].q = atoms->atom[i].qB = q_[j];
            }
            else if (atoms->atom[i].ptype == eptShell)
            {
                atoms->atom[i].q = atoms->atom[i].qB = q_[j];
                j++;
            }
        }
    }
    else
    {
        for (auto i = 0; i < atoms->nr; i++)
        {
            atoms->atom[i].q = atoms->atom[i].qB = q_[i];
        }
    }
}

void QgenAcm::updatePositions(gmx::HostVector<gmx::RVec> x,
                              t_atoms                   *atoms)
{
    for (auto i = 0; i < atoms->nr; i++)
    {
        copy_rvec(x[i], x_[i]);
    }
}

void QgenAcm::checkSupport(const Poldata *pd)
{
    bool bSupport = true;
    auto fs = pd->findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
    for (auto i = 0; i < natom_; i++)
    {
        if (!fs.parameterExists(id_[i]))
        {
            fprintf(stderr, "No charge generation support for atom %s, model %s\n",
                    id_[i].id().c_str(),
                    chargeTypeName(pd->chargeType()).c_str());
            bSupport = false;
        }
    }
    if (bSupport)
    {
        eQGEN_ = eQGEN_OK;
    }
    else
    {
        eQGEN_ = eQGEN_NOSUPPORT;
    }
}

void QgenAcm::solveEEM(FILE *fp)
{
    std::vector<double> q;

    q.resize(natom_+1, 0.0);
    MatrixWrapper lhs(natom_+1, natom_+1);

    for (int i = 0; i < natom_; i++)
    {
        for (int j = 0; j < natom_; j++)
        {
            lhs.set(i, j, Jcc_[i][j]);
        }
        lhs.set(i, i, Jcc_[i][i]);
        lhs.set(i, natom_, 1);
        lhs.set(natom_, i, -1);
    }
    lhs.set(natom_, natom_, 0);

    lhs.solve(rhs_, &q);

    double qtot    = 0;
    for (int i = 0; i < natom_; i++)
    {
        q_[i] = q[i];
        qtot += q_[i];
    }
    chieq_  = q_[natom_] = q[natom_];

    if (fp && (fabs(qtot - qtotal_) > 1e-2))
    {
        fprintf(fp, "qtot = %g, it should be %g\n", qtot, qtotal_);
    }
}

void QgenAcm::solveSQE(FILE                    *fp,
                       const Poldata           *pd,
                       const std::vector<Bond> &bonds)
{
    std::vector<double> rhs;
    
    int nbonds = bonds.size();
    MatrixWrapper lhs(nbonds, nbonds);
    rhs.resize(nbonds, 0.0);

    auto fs = pd->findForcesConst(eitBONDCORRECTIONS);
    for (int bij = 0; bij < nbonds; bij++)
    {
        auto ai    = bonds[bij].getAi()-1;
        auto aj    = bonds[bij].getAj()-1;
        Identifier bccId({id_[ai].id(), id_[aj].id()}, CanSwap::No);
        double hardness = fs.findParameterTypeConst(bccId, "hardness").value();
        double deltachi = fs.findParameterTypeConst(bccId, "electronegativity").value();

        for (int bkl = 0; bkl < nbonds; bkl++)
        {
            auto ak  = bonds[bkl].getAi()-1;
            auto al  = bonds[bkl].getAj()-1;
            
            double J = 2*(Jcc_[ai][ak] - Jcc_[ai][al] - Jcc_[aj][ak] + Jcc_[aj][al]);
            if (bij == bkl)
            {
                J += hardness;
            }
            lhs.set(bij, bkl, J);
        }
        rhs[bij] = chi0_[aj] - chi0_[ai];
        rhs[bij] -= 2*deltachi;
        
        if (fp)
        {
            fprintf(fp, "(");
            for(int i = 0; i < nbonds; i++)
            {
                fprintf(fp, " %6.2f", lhs.get(bij, i));
            }
            fprintf(fp, " ) p[%2d,%2d] = %7.2f\n", ai, aj, rhs[bij]);
        }
    }
    std::vector<double> pij, q;
    pij.resize(nbonds, 0.0);
    q.resize(natom_, 0.0);
    lhs.solve(rhs, &pij);
    for (int bij = 0; bij < nbonds; bij++)
    {
        auto ai  = bonds[bij].getAi()-1;
        auto aj  = bonds[bij].getAj()-1;
        q[ai] += pij[bij];
        q[aj] -= pij[bij];
    }
    for (int i = 0; i < natom_; i++)
    {
        q_[i] = q[i];
    }
    chieq_  = 0;

    double qtot    = 0;
    if (bHaveShell_)
    {
        for (int i = 0; i < natom_; i++)
        {
            qtot += q_[i];
        }
    }
    else
    {
        for (int i = 0; i < natom_; i++)
        {
            qtot += q_[i];
        }
    }
    if (fp && (fabs(qtot - qtotal_) > 1e-2))
    {
        fprintf(fp, "qtot = %g, it should be %g\n", qtot, qtotal_);
    }
}

int QgenAcm::generateCharges(FILE                      *fp,
                             const std::string         &molname,
                             const Poldata             *pd,
                             t_atoms                   *atoms,
                             gmx::HostVector<gmx::RVec> x,
                             const std::vector<Bond>   &bonds)
{
    if (fp)
    {
        fprintf(fp, "Generating charges for %s using %s algorithm\n",
                molname.c_str(), chargeTypeName(pd->chargeType()).c_str());
    }
    checkSupport(pd);
    if (eQGEN_OK == eQGEN_)
    {
        updateParameters(pd);
        updatePositions(x, atoms);
        calcJcc(atoms, pd->getEpsilonR(), pd->yang(), pd->rappe());
        if (pd->interactionPresent(eitBONDCORRECTIONS))
        {
            solveSQE(fp, pd, bonds);
        }
        else
        {
            calcRhs(atoms, pd->getEpsilonR());
            solveEEM(fp);
        }
        copyChargesToAtoms(atoms);
        dump(fp, atoms);
    }
    return eQGEN_;
}
}
