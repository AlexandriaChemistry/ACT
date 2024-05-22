/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

#include "qgen_acm.h"

#include <cctype>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"

#include "act/alexandria/topology.h"
#include "act/coulombintegrals/gaussian_integrals.h"
#include "act/coulombintegrals/slater_integrals.h"
#include "act/molprop/molprop.h"
#include "act/forcefield/forcefield.h"
#include "act/utility/regression.h"

namespace alexandria
{

QgenAcm::QgenAcm(const ForceField           *pd,
                 const std::vector<ActAtom> &atoms,
                 int                         qtotal)
{
    auto qt     = pd->findForcesConst(InteractionType::ELECTROSTATICS);
    ChargeType_ = potentialToChargeType(qt.potential());
    bHaveShell_ = pd->polarizable();
    eQGEN_      = eQgen::OK;
    natom_      = atoms.size();
    qtotal_     = qtotal;
    if (!ffOption(*pd, InteractionType::ELECTROSTATICS, 
                  "epsilonr", &epsilonr_))
    {
        epsilonr_ = 1;
    }
    
    auto eem = pd->findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
    for (size_t i = 0; i < atoms.size(); i++)
    {
        auto atype = pd->findParticleType(atoms[i].ffType());
        atomicNumber_.push_back(atype->atomnumber());
        auto qparam = atype->parameterConst("charge");
        if (qparam.mutability() == Mutability::ACM)
        {
            auto acmtype = atype->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            if (acmtype.id().empty())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("No ACM information for %s",
                                                               atype->id().id().c_str()).c_str()));
            }
            eta_.push_back(eem.findParameterTypeConst(acmtype, "eta").value());
            auto myrow = std::min(atype->row(), SLATER_MAX);
            row_.push_back(myrow);
            chi0_.push_back(eem.findParameterTypeConst(acmtype, "chi").value());
            nonFixed_.push_back(i);
            q_.push_back(0.0);
        }
        else
        {
            fixed_.push_back(i);
            q_.push_back(atoms[i].charge());
            eta_.push_back(0.0);
            chi0_.push_back(0.0);
            row_.push_back(0);
        }
        if (qt.potential() != Potential::COULOMB_POINT)
        {
            auto qtype = atype->interactionTypeToIdentifier(InteractionType::ELECTROSTATICS);
            zeta_.push_back(qt.findParameterTypeConst(qtype, "zeta").value());
            qdist_id_.push_back(qtype);
        }
        else
        {
            zeta_.push_back(0.0);
        }
        auto acmtp = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
        if (atype->hasInteractionType(acmtp))
        {
            auto acmtype = atype->interactionTypeToIdentifier(acmtp);
            acm_id_.push_back(acmtype);
        }
        else
        {
            acm_id_.push_back(Identifier(""));
        }
    }
    rhs_.resize(nonFixed_.size() + 1, 0);
    Jcc_.resize(nonFixed_.size() + 1, {0});
    x_.resize(atoms.size(), {0, 0, 0});

    for (size_t ii = 0; ii < nonFixed_.size()+1; ii++)
    {
        Jcc_[ii].resize(nonFixed_.size() + 1, 0);
    }
}

void QgenAcm::updateParameters(const ForceField           *pd,
                               const std::vector<ActAtom> &atoms)
{
    auto qt = pd->findForcesConst(InteractionType::ELECTROSTATICS);
    if (qt.potential() != Potential::COULOMB_POINT)
    {
        // Zeta must exist for atoms and shells in this case
        for (size_t i = 0; i < qdist_id_.size(); i++)
        {
            if (!qt.parameterExists(qdist_id_[i]))
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find %s", 
                                                               qdist_id_[i].id().c_str()).c_str()));
            }
            zeta_[i] = qt.findParameterTypeConst(qdist_id_[i], "zeta").value();
        }
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
            eta_[i]  = fs.findParameterTypeConst(acm_id_[i], "eta").value();
        }
    }
    // Update fixed charges
    for (auto &i : fixed_)
    {
        // TODO if fixed charges are in the atoms struct we do not need to look them up
        q_[i] = pd->findParticleType(atoms[i].ffType())->paramValue("charge");
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

void QgenAcm::dump(FILE                       *fp,
                   const std::vector<ActAtom> *atoms) const
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
        for (size_t i = 0; i < atoms->size(); i++)
        {
            for (auto m = 0; m < DIM; m++)
            {
                mu[m] += (*atoms)[i].charge() * x_[i][m] * ENM2DEBYE;
            }
            fprintf(fp, "%4s %4s%5lu %8g %8g",
                    "",
                    (*atoms)[i].name().c_str(), i+1, eta_[i], chi0_[i]);
            fprintf(fp, " %3d %8.5f %8.4f\n", row_[i], q_[i], zeta_[i]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "|mu| = %8.3f ( %8.3f  %8.3f  %8.3f )\n\n",
                norm(mu), mu[XX], mu[YY], mu[ZZ]);
    }
}

const char *QgenAcm::status() const
{
    std::map<eQgen, const char *> stat = {
        { eQgen::OK,           "Charge generation finished correctly" },
        { eQgen::NOTCONVERGED, "Charge generation did not converge." },
        { eQgen::NOSUPPORT,    "No charge generation support for (some of) the atomtypes." },
        { eQgen::MATRIXSOLVER, "A problem occurred solving a matrix equation." },
        { eQgen::INCORRECTMOL, "Molecule information inconsistent, e.g. missing bonds" },
        { eQgen::ERROR,        "Unknown error in charge generation" }
    };
    return stat[eQGEN_];
}

double QgenAcm::calcSij(int i, int j)
{
    double dist, dism, Sij = 1.0;
    rvec   dx;
    int    l, m, tag;

    if (atomicNumber_[i] == 0 || atomicNumber_[j] == 0)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot compute Sij for shells").c_str()));
    }
    rvec_sub(x_[i], x_[j], dx);
    dist = norm(dx);
    if ((dist < 0.118) && (atomicNumber_[i] != 1) && (atomicNumber_[j] != 1))
    {
        Sij = Sij*1.64;
    }
    else if ((dist < 0.122) && (atomicNumber_[i] != 1) && (atomicNumber_[j] != 1))
    {
        if ((atomicNumber_[i] != 8) && (atomicNumber_[j] != 8))
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
        if ((atomicNumber_[i] == 6) && (atomicNumber_[j] == 8))
        {
            tag = i;
        }
        else if ((atomicNumber_[i] == 8) && (atomicNumber_[j] == 6))
        {
            tag = j;
        }
        if (tag != 0)
        {
            printf("found CO\n");
            for (l = 0; (l < natom_); l++)
            {
                if (atomicNumber_[l] == 1)
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
    else if ((atomicNumber_[i] == 6) && (atomicNumber_[j] == 8))
    {
        Sij = Sij*1.03;
    }
    else if (((atomicNumber_[j] == 6) && (atomicNumber_[i] == 7) && (dist < 0.15)) ||
             ((atomicNumber_[i] == 6) && (atomicNumber_[j] == 7) && (dist < 0.15)))
    {
        if (atomicNumber_[i] == 6)
        {
            tag = i;
        }
        else
        {
            tag = j;
        }
        for (l = 0; (l < natom_); l++)
        {
            if (atomicNumber_[l] == 8)
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
                auto j0 = eta_[j];
                if ((bYang || bRappe) && (atomicNumber_[i] == 1))
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

double QgenAcm::calcJcs(int                         core_ndx,
                        const std::vector<ActAtom> &atoms,
                        double                      epsilonr)
{
    double Jcs     = 0;
    for (size_t i = 0; i < fixed_.size(); i++)
    {
        auto fixed_ndx = fixed_[i];
        // Exclude interaction between a particle and its own shell
        const auto myshells = atoms[core_ndx].shells();
        if (std::find(myshells.begin(), myshells.end(), fixed_ndx) 
            == myshells.end())
        {
            Jcs += (q_[fixed_ndx]*calcJ(x_[core_ndx],
                                        x_[fixed_ndx],
                                        zeta_[core_ndx],
                                        zeta_[fixed_ndx],
                                        row_[core_ndx],
                                        row_[fixed_ndx],
                                        epsilonr));  
            if (debug)
            {
                fprintf(debug, "calcJcs: core_ndx: %d fixed_ndx: %d shell_charge: %0.1f\n", 
                        core_ndx, fixed_ndx, q_[fixed_ndx]);
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "Jcs:%0.3f\n", Jcs);
    }
    return Jcs;
}

void QgenAcm::calcRhs(const std::vector<ActAtom> &atoms,
                      double                      epsilonr)
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
        // TODO Check this stuff: Delta_Eta * q0
        // rhs_[i] += eta_[nfi]*q0_[i];
        if (bHaveShell_)
        {
            // Core-Shell interaction
            rhs_[i] -= 0.5*calcJcs(nfi, atoms, epsilonr);
            // Hardness of atom multiplied by charge of shell(s)
            for(auto shellId : atoms[nfi].shells())
            {
                rhs_[i] -= eta_[nfi]*q_[shellId];
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "QgenACM: qtotal_ = %d qfixed = %g\n", qtotal_, qfixed);
    }
    rhs_[nonFixed_.size()] = qtotal_ - qfixed;
}

void QgenAcm::copyChargesToAtoms(std::vector<ActAtom> *atoms)
{
    for (size_t i = 0; i < atoms->size(); i++)
    {
        (*atoms)[i].setCharge(q_[i]);
    }
}

void QgenAcm::updatePositions(const std::vector<gmx::RVec> &x)
{
    if (x.size() - x_.size() != 0)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Arrays not equally long. New %zu Old %zu. Help!", x.size(), x_.size()).c_str()));
    }
    for (size_t i = 0; i < x.size(); i++)
    {
        copy_rvec(x[i], x_[i]);
    }
}

void QgenAcm::checkSupport(const ForceField *pd)
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

    int info = lhs.solve(rhs_, &q);

    if (info == 0)
    {
        double qtot    = 0;
        for (size_t i = 0; i < nelem; i++)
        {
            q_[nonFixed_[i]] = q[i];
            qtot += q[i];
        }
        
        if (fp && (fabs(qtot - qtotal_) > 1e-2))
        {
            fprintf(fp, "qtot = %g, it should be %d rhs[%zu] = %g\n",
                    qtot, qtotal_, nelem, rhs_[nelem]);
        }
    }
    else
    {
        eQGEN_ = eQgen::MATRIXSOLVER;
    }    
}

void QgenAcm::getBccParams(const ForceField  *pd,
                           int                ai,
                           int                aj,
                           double             bondorder,
                           double            *delta_chi,
                           double            *delta_eta)
{
    auto itype    = InteractionType::BONDCORRECTIONS;
    auto fs       = pd->findForcesConst(itype);
    bool swapped = false;
    int  nf      = static_cast<int>(nonFixed_.size());
    if (ai >= nf || aj >= nf)
    {
        fprintf(stderr, "Atom index %d or %d out of range (max %zu). Check your input.\n", ai, aj, nonFixed_.size());
        eQGEN_ = eQgen::INCORRECTMOL;
        return;
    }
    Identifier bccId({acm_id_[nonFixed_[ai]].id(), 
            acm_id_[nonFixed_[aj]].id()}, { bondorder }, fs.canSwap());
    if (!fs.parameterExists(bccId))
    {
        if (CanSwap::Yes == fs.canSwap())
        {
            eQGEN_ = eQgen::NOSUPPORT;
            return;
        }
        else
        {
            bccId   = Identifier({acm_id_[nonFixed_[aj]].id(),
                    acm_id_[nonFixed_[ai]].id()},
                { bondorder }, fs.canSwap());
            swapped = true;
        }
    }
    if (!fs.parameterExists(bccId))
    {
        eQGEN_ = eQgen::NOSUPPORT;
        return;
    }
    *delta_chi = fs.findParameterTypeConst(bccId, "delta_chi").value();
    if (swapped)
    {
        *delta_chi *= -1;
    }
    *delta_eta = fs.findParameterTypeConst(bccId, "delta_eta").value();
}

void QgenAcm::solveSQE(FILE                    *fp,
                      const ForceField        *pd,
                      const std::vector<Bond> &bonds)
{
    std::vector<double> myq(nonFixed_.size(), 0.0);
    int info   = 0;
    int nbonds = bonds.size();
    if (nbonds > 0)
    {
        MatrixWrapper lhs(nbonds, nbonds);
        std::vector<double> rhs(nbonds, 0.0);

        std::vector<double> delta_chis;
        // First fill the matrix
        for (int bij = 0; bij < nbonds; bij++)
        {
            auto ai      = bonds[bij].aI();
            auto aj      = bonds[bij].aJ();
            
            double delta_chi = 0, delta_eta = 0;
            getBccParams(pd, ai, aj, bonds[bij].bondOrder(),
                         &delta_chi, &delta_eta);
            if (eQgen::OK != eQGEN_)
            {
                return;
            }
            if (debug)
            {
                fprintf(debug, "Delta_chi %10g Delta_Eta %10g\n",
                        delta_chi, delta_eta);
            }
            // The bonds use the original numbering, to get to the ACM data
            // we have to map to numbers including shells.
            delta_chis.push_back(delta_chi);
            
            for (int bkl = 0; bkl < nbonds; bkl++)
            {
                auto ak  = bonds[bkl].aI();
                auto al  = bonds[bkl].aJ();
            
                double J = (Jcc_[ai][ak] - Jcc_[ai][al] - Jcc_[aj][ak] + Jcc_[aj][al]);
                if (bij == bkl)
                {
                    J += delta_eta;
                }
                lhs.set(bij, bkl, J);
            }
        }
        if (debug)
        {
            fprintf(debug, "Jsqe\n");
            for(int i = 0; i < nbonds; i++)
            {
                for(int j = 0; j <= i; j++)
                {
                    fprintf(debug, " %7.2f", lhs.get(i, j));
                }
                fprintf(debug, "\n");
            }
        }
        // Now fill the right hand side
        std::vector<double> chi_corr;
        chi_corr.resize(nonFixed_.size(), 0.0);
        for(size_t i = 0; i < nonFixed_.size(); i++)
        {
            int  ai  = static_cast<int>(i);
            auto nfi = nonFixed_[ai];
            chi_corr[i] += chi0_[nfi];
            for (int bij = 0; bij < nbonds; bij++)
            {
                if (ai == bonds[bij].aI())
                {
                    chi_corr[i] += delta_chis[bij];
                }
                else if (ai == bonds[bij].aJ())
                {
                    chi_corr[i] -= delta_chis[bij];
                }
            }
            double qcorr = (1.0*qtotal_)/nonFixed_.size();
            for(size_t l = 0; l < nonFixed_.size(); l++)
            {
                chi_corr[i] += Jcc_[i][l]*qcorr;
            }
            // TODO Check this! Only to be done when there are shells!
            //if (false && nonFixed_.size() < static_cast<size_t>(natom_))
            //{
            //   chi_corr[i] += eta_[nfi] * q_[myShell_.find(nfi)->second];
            //   chi_corr[i] += 0.5*calcJcs(nfi, epsilonr_);
            //}
        }
        for (int bij = 0; bij < nbonds; bij++)
        {
            auto ai    = bonds[bij].aI();
            auto aj    = bonds[bij].aJ();
            rhs[bij]   = chi_corr[aj] - chi_corr[ai];
        }
        
        std::vector<double> pij(nbonds, 0.0); 
        info = lhs.solve(rhs, &pij);
        if (info > 0)
        {
            eQGEN_ = eQgen::MATRIXSOLVER;
            return;
        }
        if (fp)
        {
            fprintf(fp, "rhs: ");
            for(auto &p: rhs)
            {
                fprintf(fp, " %8g", p);
            }
            fprintf(fp, "\n");
            fprintf(fp, "pij: ");
            for(auto &p: pij)
            {
                fprintf(fp, " %8g", p);
            }
            fprintf(fp, "\n");
        }
        for (int bij = 0; bij < nbonds; bij++)
        {
            auto ai  = bonds[bij].aI();
            auto aj  = bonds[bij].aJ();
            myq[ai] += pij[bij];
            myq[aj] -= pij[bij];
        }
    }
    double qfixed = qtotal_;
    for (size_t i = 0; i < fixed_.size(); i++)
    {
        qfixed -= q_[fixed_[i]];
    }
    for (size_t i = 0; i < nonFixed_.size(); i++)
    {
        auto nfi = nonFixed_[i];
        //        myq[i] += qfixed/nonFixed_.size();
        //if (nonFixed_.size() < static_cast<size_t>(natom_))
        //{
        //   myq[i] -= q_[myShell_.find(nfi)->second];
        //}
        q_[nfi] = myq[i] + qfixed/nonFixed_.size();
    }

    double qtot    = 0;
    for (int i = 0; i < natom_; i++)
    {
        qtot += q_[i];
    }
    if (fp)
    {
        fprintf(fp, "q: ");
        for (auto &qq : q_)
        {
            fprintf(fp, " %9g", qq);
        }
        fprintf(fp, "\n");
        if (fabs(qtot - qtotal_) > 1e-2)
        {
            fprintf(fp, "qtot = %g, it should be %d\n", qtot, qtotal_);
        }
    }
}

eQgen QgenAcm::generateCharges(FILE                         *fp,
                               const std::string            &molname,
                               const ForceField             *pd,
                               std::vector<ActAtom>         *atoms,
                               const std::vector<gmx::RVec> &x,
                               const std::vector<Bond>      &bonds)
{
    if (nonFixed_.empty())
    {
        return eQgen::OK;
    }
    if (fp)
    {
        fprintf(fp, "Generating charges for %s using %s algorithm\n",
                molname.c_str(), 
                chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str());
    }
    checkSupport(pd);
    if (eQgen::OK == eQGEN_)
    {
        updateParameters(pd, *atoms);
        updatePositions(x);
        calcJcc(epsilonr_, pd->yang(), pd->rappe());
        if (pd->interactionPresent(InteractionType::BONDCORRECTIONS) &&
            pd->chargeGenerationAlgorithm() == ChargeGenerationAlgorithm::SQE)
        {
            solveSQE(fp, pd, bonds);
        }
        else
        {
            calcRhs(*atoms, epsilonr_);
            solveEEM(fp);
        }
        if (eQgen::OK == eQGEN_)
        {
            copyChargesToAtoms(atoms);
            dump(fp, atoms);
        }
    }
    return eQGEN_;
}

} // namespace
