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
    auto qt     = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    ChargeType_ = name2ChargeType(qt.optionValue("chargetype"));
    bHaveShell_ = pd->polarizable();
    eQGEN_      = eQgen::OK;
    natom_      = atoms->nr;
    qtotal_     = qtotal;
    if (debug)
    {
        fprintf(debug, "QgenACM: qtotal = %g\n", qtotal_);
    }
    auto eem = pd->findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
    for (int i = 0; i < atoms->nr; i++)
    {
        auto atype = pd->findParticleType(*atoms->atomtype[i]);
        atomnr_.push_back(atype->atomnumber());
        auto qparam = atype->parameter("charge");
        if (qparam.mutability() != Mutability::Fixed)
        {
            auto acmtype = atype->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            if (acmtype.id().empty())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("No ACM information for %s",
                                                               atype->id().id().c_str()).c_str()));
            }
            jaa_.push_back(eem.findParameterTypeConst(acmtype, "jaa").value());
            auto myrow = std::min(atype->row(), SLATER_MAX);
            row_.push_back(myrow);
            chi0_.push_back(eem.findParameterTypeConst(acmtype, "chi").value());
            nonFixed_.push_back(i);
            nfToGromacs_.insert({ i, nonFixed_.size()-1 });
            q_.push_back(0.0);
            // If this particle has a shell, it is assumed to be the next particle
            // TODO: check the polarization array instead.
            if (atype->hasInteractionType(InteractionType::POLARIZATION) &&
                i < atoms->nr-1 &&
                atoms->atom[i+1].ptype == eptShell)
            {
                myShell_.insert({ i, i+1 });
            }
        }
        else
        {
            fixed_.push_back(i);
            q_.push_back(atoms->atom[i].q);
            jaa_.push_back(0.0);
            chi0_.push_back(0.0);
            row_.push_back(0);
        }
        auto qtype = atype->interactionTypeToIdentifier(InteractionType::CHARGEDISTRIBUTION);
        zeta_.push_back(qt.findParameterTypeConst(qtype, "zeta").value());
        qdist_id_.push_back(qtype);
        auto acmtype = atype->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
        acm_id_.push_back(acmtype);
    }
    rhs_.resize(nonFixed_.size() + 1, 0);
    Jcc_.resize(nonFixed_.size() + 1);
    x_.resize(atoms->nr);

    for (size_t ii = 0; ii < nonFixed_.size()+1; ii++)
    {
        Jcc_[ii].resize(nonFixed_.size() + 1, 0);
    }
}

void QgenAcm::updateParameters(const Poldata *pd)
{
    // Zeta must exist for atoms and shells
    auto qt = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    for (size_t i = 0; i < qdist_id_.size(); i++)
    {
        if (!qt.parameterExists(qdist_id_[i]))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find %s", 
                                                           qdist_id_[i].id().c_str()).c_str()));
        }
        zeta_[i] = qt.findParameterTypeConst(qdist_id_[i], "zeta").value();
    }
    // ACM typically for atoms only
    auto fs = pd->findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
    for (size_t i = 0; i < acm_id_.size(); i++)
    {
        if (!acm_id_[i].id().empty())
        {
            if (!fs.parameterExists(acm_id_[i]))
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find %s", 
                                                               acm_id_[i].id().c_str()).c_str()));
            }
            chi0_[i] = fs.findParameterTypeConst(acm_id_[i], "chi").value();
            jaa_[i]  = fs.findParameterTypeConst(acm_id_[i], "jaa").value();
        }
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
    if (fp && eQGEN_ == eQgen::OK)
    {
        rvec  mu = { 0, 0, 0 };

        fprintf(fp, "Jcc_ matrix:\n");
        for (size_t i = 0; i < nonFixed_.size(); i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                fprintf(fp, "  %6.2f", Jcc_[i][j]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");

        fprintf(fp, "RHS_ matrix:\n");
        for (size_t i = 0; i < nonFixed_.size(); i++)
        {
            fprintf(fp, "%2zu RHS_ = %10g\n", i, rhs_[i]);
        }
        fprintf(fp, "\n");

        fprintf(fp, "                                          Core                 Shell\n");
        fprintf(fp, "Res  Atom   Nr       J0    chi0  row   q   zeta        row    q     zeta\n");
        for (int i = 0; i < atoms->nr; i++)
        {
            for (auto m = 0; m < DIM; m++)
            {
                mu[m] += atoms->atom[i].q * x_[i][m] * ENM2DEBYE;
            }
            fprintf(fp, "%4s %4s%5d %8g %8g",
                    *(atoms->resinfo[atoms->atom[i].resind].name),
                    *(atoms->atomname[i]), i+1, jaa_[i], chi0_[i]);
            fprintf(fp, " %3d %8.5f %8.4f\n", row_[i], q_[i], zeta_[i]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "|mu| = %8.3f ( %8.3f  %8.3f  %8.3f )\n\n",
                norm(mu), mu[XX], mu[YY], mu[ZZ]);
    }
}

const char *QgenAcm::message() const
{
    switch (eQGEN_)
    {
        case eQgen::OK:
            return "Charge generation finished correctly";
        case eQgen::NOTCONVERGED:
            return "Charge generation did not converge.";
        case eQgen::NOSUPPORT:
            return "No charge generation support for (some of) the atomtypes.";
        case eQgen::ERROR:
            return "Unknown error in charge generation";
    }
    return nullptr;
}

double QgenAcm::calcSij(int i, int j)
{
    double dist, dism, Sij = 1.0;
    rvec   dx;
    int    l, m, tag;

    if (atomnr_[i] == 0 || atomnr_[j] == 0)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot compute Sij for shells").c_str()));
    }
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
                               i, l, qdist_id_[tag].id().c_str(), qdist_id_[l].id().c_str(), sqrt(dism));
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
    
    if (r == 0 && ChargeType_ == ChargeType::Point)
    {
        gmx_fatal(FARGS, "Zero distance between atoms!\n");
    }
    switch (ChargeType_)
    {
    case ChargeType::Point:
        {
            eTot = Coulomb_PP(r);
            break;
        }
    case ChargeType::Gaussian:
        {
            eTot = Coulomb_GG(r, zetaI, zetaJ);
            break;
        }
    case ChargeType::Slater:
        {
            eTot = Coulomb_SS(r, rowI, rowJ, zetaI, zetaJ);
            break;
        }
    }
    return (ONE_4PI_EPS0*eTot)/(epsilonr*ELECTRONVOLT);
}

void QgenAcm::calcJcc(double epsilonr,
                      bool   bYang, 
                      bool   bRappe)
{
    auto Jcc = 0.0;
    for (size_t n = 0; n < nonFixed_.size(); n++)
    {
        auto i = nonFixed_[n];
        for (size_t m = n; m < nonFixed_.size(); m++)
        {
            auto j = nonFixed_[m];
            if (j != i)
            {
                Jcc = calcJ(x_[i],
                            x_[j],
                            zeta_[i],
                            zeta_[j],
                            row_[i],
                            row_[j],
                            epsilonr);
                if (bYang)
                {
                    Jcc *= calcSij(i, j);
                }
                Jcc_[n][m] = Jcc_[m][n] = (0.5 * Jcc);
            }
            else
            {
                auto j0 = jaa_[j];
                if ((bYang || bRappe) && (atomnr_[i] == 1))
                {
                    auto zetaH = 1.0698;
                    j0 = (1+q_[i]/zetaH)*j0; /* Eqn. 21 in Rappe1991. */
                    if (j0 < 0 && !bWarned_)
                    {
                        bWarned_ = true;
                    }
                }
                Jcc_[n][n] = (j0 > 0) ? j0 : 0;
            }
        }
    }
}

double QgenAcm::calcJcs(int    core_ndx,
                        double epsilonr)
{
    double Jcs     = 0;
    int    shellId = -1;
    if (myShell_.find(core_ndx) != myShell_.end())
    {
        shellId =  myShell_.find(core_ndx)->second;
    }
    for (size_t i = 0; i < fixed_.size(); i++)
    {
        auto shell_ndx = fixed_[i];
        // Exclude interaction between a particle and its shell
        // TODO: check exclusion array instead
        if (shell_ndx !=  shellId)
        {
            Jcs += (q_[shell_ndx]*calcJ(x_[core_ndx],
                                        x_[shell_ndx],
                                        zeta_[core_ndx],
                                        zeta_[shell_ndx],
                                        row_[core_ndx],
                                        row_[shell_ndx],
                                        epsilonr));  
            if (debug)
            {
                fprintf(debug, "core_ndx: %d shell_ndx: %d shell_charge: %0.1f\n", 
                        core_ndx, shell_ndx, q_[shell_ndx]);
            }
        }
    }
    Jcs *= 0.5;
    if (debug)
    {
        fprintf(debug, "Jcs:%0.3f\n", Jcs);
    }
    return Jcs;
}

void QgenAcm::calcRhs(double epsilonr)
{
    // Compute total fixed charge
    double qfixed = 0.0;
    for (size_t i = 0; i < fixed_.size(); i++)
    {
        qfixed += q_[fixed_[i]];
    }
    for (size_t i = 0; i < nonFixed_.size(); i++)
    {
        auto nfi = nonFixed_[i];
        // Electronegativity
        rhs_[i]  = -chi0_[nfi]; 
        // TODO Check this stuff: Hardness * q0
        // rhs_[i] += jaa_[nfi]*q0_[i];
        if (bHaveShell_)
        {
            // Core-Shell interaction
            rhs_[i] -= calcJcs(nfi, epsilonr);
            // Hardness of atom multiplied by charge of shell
            if (myShell_.find(nfi) != myShell_.end())
            {
                auto shellId = myShell_.find(nfi)->second;
                rhs_[i]     -= jaa_[nfi]*q_[shellId];
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "QgenACM: qtotal_ = %g qfixed = %g\n", qtotal_, qfixed);
    }
    rhs_[nonFixed_.size()] = qtotal_ - qfixed;
}

void QgenAcm::copyChargesToAtoms(t_atoms *atoms)
{
    for (auto i = 0; i < atoms->nr; i++)
    {
        atoms->atom[i].q = atoms->atom[i].qB = q_[i];
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
    auto fs = pd->findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
    for (size_t i = 0; i < acm_id_.size(); i++)
    {
        if (!acm_id_[i].id().empty())
        {
            if (!fs.parameterExists(acm_id_[i]))
            {
                fprintf(stderr, "No %s charge generation support for atom %s.\n",
                        chargeTypeName(ChargeType_).c_str(),
                        acm_id_[i].id().c_str());
                bSupport = false;
            }
        }
    }
    if (bSupport)
    {
        eQGEN_ = eQgen::OK;
    }
    else
    {
        eQGEN_ = eQgen::NOSUPPORT;
    }
}

void QgenAcm::solveEEM(FILE *fp)
{
    std::vector<double> q;
    size_t              nelem = nonFixed_.size();

    q.resize(nelem+1, 0.0);
    MatrixWrapper lhs(nelem+1, nelem+1);

    for (size_t i = 0; i < nelem; i++)
    {
        for (size_t j = 0; j < nelem; j++)
        {
            lhs.set(i, j, Jcc_[i][j]);
        }
        lhs.set(i, i, Jcc_[i][i]);
        lhs.set(i, nelem, 1);
        lhs.set(nelem, i, -1);
    }
    lhs.set(nelem, nelem, 0);

    lhs.solve(rhs_, &q);

    double qtot    = 0;
    for (size_t i = 0; i < nelem; i++)
    {
        q_[nonFixed_[i]] = q[i];
        qtot += q[i];
    }

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
    
    auto itype = InteractionType::BONDCORRECTIONS;
    auto fs    = pd->findForcesConst(itype);
    for (int bij = 0; bij < nbonds; bij++)
    {
        auto ai      = bonds[bij].getAi()-1;
        auto aj      = bonds[bij].getAj()-1;
        auto canSwap = fs.canSwap();
        bool swapped = false;
        // The bonds use the original numbering, to get to the ACM data
        // we have to map to numbers including shells.
        Identifier bccId({acm_id_[nonFixed_[ai]].id(), 
                acm_id_[nonFixed_[aj]].id()},
            bonds[bij].getBondOrder(), canSwap);
        if (!fs.parameterExists(bccId))
        {
            if (CanSwap::Yes == canSwap)
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find identifier %s in list %s", bccId.id().c_str(), interactionTypeToString(itype).c_str()).c_str()));
            }
            else
            {
                bccId   = Identifier({acm_id_[nonFixed_[aj]].id(),
                        acm_id_[nonFixed_[ai]].id()},
                    bonds[bij].getBondOrder(), canSwap);
                swapped = true;
            }
        }
        double hardness = fs.findParameterTypeConst(bccId, "hardness").value();
        double deltachi = fs.findParameterTypeConst(bccId, "electronegativity").value();
        if (swapped)
        {
            deltachi *= -1;
        }

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
        rhs[bij] = chi0_[nonFixed_[aj]] - chi0_[nonFixed_[ai]];
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
    q.resize(nonFixed_.size(), 0.0);
    lhs.solve(rhs, &pij);
    for (int bij = 0; bij < nbonds; bij++)
    {
        auto ai  = bonds[bij].getAi()-1;
        auto aj  = bonds[bij].getAj()-1;
        q[ai] += pij[bij];
        q[aj] -= pij[bij];
    }
    // Compute total fixed charge
    double qfixed = 0.0;
    for (size_t i = 0; i < fixed_.size(); i++)
    {
        qfixed += q_[fixed_[i]];
    }
    for (size_t i = 0; i < nonFixed_.size(); i++)
    {
        q[i] += (qtotal_-qfixed)/nonFixed_.size();
        q_[nonFixed_[i]] = q[i];
    }

    double qtot    = 0;
    for (int i = 0; i < natom_; i++)
    {
        qtot += q_[i];
    }
    if (fp && (fabs(qtot - qtotal_) > 1e-2))
    {
        fprintf(fp, "qtot = %g, it should be %g\n", qtot, qtotal_);
    }
}

eQgen QgenAcm::generateCharges(FILE                      *fp,
                               const std::string         &molname,
                               const Poldata             *pd,
                               t_atoms                   *atoms,
                               gmx::HostVector<gmx::RVec> x,
                               const std::vector<Bond>   &bonds)
{
    if (fp)
    {
        fprintf(fp, "Generating charges for %s using %s algorithm\n",
                molname.c_str(), 
                chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str());
    }
    checkSupport(pd);
    if (eQgen::OK == eQGEN_)
    {
        updateParameters(pd);
        updatePositions(x, atoms);
        calcJcc(pd->getEpsilonR(), pd->yang(), pd->rappe());
        if (pd->interactionPresent(InteractionType::BONDCORRECTIONS))
        {
            solveSQE(fp, pd, bonds);
        }
        else
        {
            calcRhs(pd->getEpsilonR());
            solveEEM(fp);
        }
        copyChargesToAtoms(atoms);
        dump(fp, atoms);
    }
    return eQGEN_;
}
}
