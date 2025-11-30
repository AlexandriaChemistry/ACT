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

#include "qgen_acm.h"

#include <cctype>
#include <numeric>

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

QgenAcm::QgenAcm(ForceField                 *pd,
                 const std::vector<ActAtom> &atoms,
                 const std::vector<Bond>    &bonds,
                 int                         qtotal)
{
    auto entype = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
    if (!pd->interactionPresent(entype))
    {
        return;
    }
    auto qt     = pd->findForces(InteractionType::ELECTROSTATICS);
    ChargeDistributionType_ = potentialToChargeDistributionType(qt->potential());
    bHaveShell_ = pd->polarizable();
    eQGEN_      = eQgen::OK;
    natom_      = atoms.size();
    qtotal_     = qtotal;
    if (!ffOption(*pd, InteractionType::ELECTROSTATICS, 
                  "epsilonr", &epsilonr_))
    {
        epsilonr_ = 1;
    }
    auto eem = pd->findForces(entype);
    auto bctype = InteractionType::BONDCORRECTIONS;
    bool haveBCC = pd->interactionPresent(bctype);
    std::vector<std::string> acmtypes;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        auto atype = pd->findParticleType(atoms[i].ffType());
        atomicNumber_.push_back(atype->atomnumber());
        auto qparm = atype->parameterConst("charge");
        if (atype->hasInteractionType(entype))
        {
            const auto &acmtype = atype->interactionTypeToIdentifier(entype);
            if (acmtype.id().empty())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("No ACM information for %s",
                                                               atype->id().id().c_str()).c_str()));
            }
            acmtypes.push_back(acmtype.id());
        }
        else
        {
            acmtypes.push_back("");
        }
        // For EEM we need only the Mutability::ACM particles, but if we use SQE
        // we need all atoms connected by bonds. As a result when using SQE the
        // non-ACM particles will be used in the SQE algorithm as well.
        if (atype->hasInteractionType(entype) && 
            ((qparm.mutability() == Mutability::ACM) ||
             (haveBCC && atoms[i].pType() == ActParticle::Atom)))
        {
            eta_.push_back(eem->findParameterTypeConst(acmtypes.back(), "eta").value());
            auto myrow = std::min(atype->row(), SLATER_MAX);
            row_.push_back(myrow);
            chi0_.push_back(eem->findParameterTypeConst(acmtypes.back(), "chi").value());
            nonFixed_.push_back(i);
            q_.push_back(0.0);
        }
        else
        {
            fixed_.push_back(i);
            q_.push_back(atoms[i].charge());
            charge_.push_back(atype->parameter("charge"));
            eta_.push_back(0.0);
            chi0_.push_back(0.0);
            row_.push_back(0);
        }
        if (qt->potential() != Potential::COULOMB_POINT)
        {
            auto qtype = atype->interactionTypeToIdentifier(InteractionType::ELECTROSTATICS);
            auto qdist = &(qt->findParameters(qtype)->find("zeta")->second);
            zeta_.push_back(qdist->value());
            qdist_id_.push_back(qdist);
        }
        else
        {
            zeta_.push_back(0.0);
        }
        if (atype->hasInteractionType(entype))
        {
            //! \todo this code does not make sense, if there is no BCC, the EEM
            // values will not be stored and updated.
            auto acmtype = atype->interactionTypeToIdentifier(entype);
            auto acm     = eem->findParameters(acmtype);
            acm_id_.push_back(acm);
        }
        else
        {
            //! \todo check whether this actually can occur
            acm_id_.push_back(nullptr);
        }
    }
    if (pd->interactionPresent(bctype))
    {
        auto fs = pd->findForces(bctype);
        for(auto &b : bonds)
        {
            if (!acmtypes.empty())
            {
                if (nonFixed_.empty() ||
                    static_cast<size_t>(b.aI()) >= nonFixed_.size() ||
                    static_cast<size_t>(b.aJ()) >= nonFixed_.size())
                {
                    GMX_THROW(gmx::InvalidInputError("Inconsistent EEM parameters in force field"));
                }
                auto ai = acmtypes[nonFixed_[b.aI()]];
                auto aj = acmtypes[nonFixed_[b.aJ()]];
                Identifier bccId( { ai, aj }, { b.bondOrder() }, fs->canSwap());
                double dcf = 1;
                if (!fs->parameterExists(bccId))
                {
                    if (CanSwap::Yes == fs->canSwap())
                    {
                        eQGEN_ = eQgen::NOSUPPORT;
                        return;
                    }
                    else
                    {
                        bccId = Identifier( { aj, ai }, { b.bondOrder() }, fs->canSwap());
                        dcf = -1;
                    }
                }
                if (fs->parameterExists(bccId))
                {
                    bcc_.push_back(fs->findParameters(bccId));
                    dchi_factor_.push_back(dcf);
                }
            }
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

void QgenAcm::updateParameters()
{
    if (!zeta_.empty())
    {
        // Zeta must exist for atoms and shells in this case
        for (size_t i = 0; i < qdist_id_.size(); i++)
        {
            zeta_[i] = qdist_id_[i]->value();
        }
    }
    // ACM typically for atoms only
    for (size_t i = 0; i < acm_id_.size(); i++)
    {
        if (acm_id_[i])
        {
            chi0_[i] = acm_id_[i]->find("chi")->second.value();
            eta_[i]  = acm_id_[i]->find("eta")->second.value();
        }
    }
    // Update fixed charges
    for (size_t i = 0 ; i < fixed_.size(); i++)
    {
        //! \todo if fixed charges are in the atoms struct we do not need to look them up
        q_[fixed_[i]] = charge_[i]->value();
    }
}

double QgenAcm::getQ(size_t atom)
{
    if (atom < natom_)
    {
        return q_[atom];
    }
    return 0;
}

int QgenAcm::getRow(size_t atom)
{
    if (atom < natom_)
    {
        return row_[atom];
    }
    return 0;

}

double QgenAcm::getZeta(size_t atom)
{
    if (atom < natom_)
    {
        return zeta_[atom];
    }
    return 0;
}

//! \return Coulomb for a point charge without units
static double Coulomb_PP(double r)
{
    return 1/r;
}

void QgenAcm::dump(MsgHandler                 *msg_handler,
                   const std::vector<ActAtom> *atoms) const
{
    if (msg_handler && msg_handler->debug() && eQGEN_ == eQgen::OK)
    {
        auto tw = msg_handler->twDebug();
        if (!tw)
        {
            return;
        }
        rvec  mu = { 0, 0, 0 };
        std::vector<std::string> msgs;
        msgs.push_back(gmx::formatString("Jcc_ matrix:"));
        for (size_t i = 0; i < nonFixed_.size(); i++)
        {
            std::string msg;
            for (size_t j = 0; j <= i; j++)
            {
                msg += gmx::formatString("  %6.2f", Jcc_[i][j]);
            }
            msgs.push_back(msg);
        }
        msgs.push_back("\n");
        msgs.push_back("RHS_ matrix:");
        for (size_t i = 0; i < nonFixed_.size(); i++)
        {
            msgs.push_back(gmx::formatString("%2zu RHS_ = %10g\n", i, rhs_[i]));
        }
        msgs.push_back("\n");

        msgs.push_back("                                          Core                 Shell\n");
        msgs.push_back("Res  Atom   Nr       J0    chi0  row   q   zeta        row    q     zeta\n");
        for (size_t i = 0; i < atoms->size(); i++)
        {
            for (auto m = 0; m < DIM; m++)
            {
                mu[m] += (*atoms)[i].charge() * x_[i][m] * ENM2DEBYE;
            }
            auto msg = gmx::formatString("%4s %4s%5lu %8g %8g",
                                             "",
                                         (*atoms)[i].name().c_str(), i+1, eta_[i], chi0_[i]);
            msg += gmx::formatString(" %3d %8.5f %8.4f\n", row_[i], q_[i], zeta_[i]);
            msgs.push_back(msg);
        }
        msgs.push_back("\n");
        msgs.push_back(gmx::formatString("|mu| = %8.3f ( %8.3f  %8.3f  %8.3f )\n\n",
                                         norm(mu), mu[XX], mu[YY], mu[ZZ]));
        for (const auto &msg: msgs)
        {
            tw->writeLine(msg);
        }
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
    
    if (r == 0 && ChargeDistributionType_ == ChargeDistributionType::Point)
    {
        gmx_fatal(FARGS, "Zero distance between atoms!\n");
    }
    switch (ChargeDistributionType_)
    {
    case ChargeDistributionType::Point:
        {
            eTot = Coulomb_PP(r);
            break;
        }
    case ChargeDistributionType::Gaussian:
        {
            eTot = Coulomb_GG(r, zetaI, zetaJ);
            break;
        }
    case ChargeDistributionType::Slater:
        {
            eTot = Coulomb_SS(r, rowI, rowJ, zetaI, zetaJ);
            break;
        }
    }
    return (ONE_4PI_EPS0*eTot)/(epsilonr*ELECTRONVOLT);
}

void QgenAcm::calcJcc(double epsilonr)
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
                Jcc_[n][m] = Jcc_[m][n] = (0.5 * Jcc);
            }
            else
            {
                auto j0 = eta_[j];
                Jcc_[n][n] = (j0 > 0) ? j0 : 0;
            }
        }
    }
}

double QgenAcm::calcJcs(MsgHandler                 *msg_handler,
                        int                         core_ndx,
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
            if (msg_handler && msg_handler->debug())
            {
                msg_handler->writeDebug(gmx::formatString("calcJcs: core_ndx: %d fixed_ndx: %d shell_charge: %0.1f", 
                                                          core_ndx, fixed_ndx, q_[fixed_ndx]));
            }
        }
    }
    if (msg_handler && msg_handler->debug())
    {
        msg_handler->writeDebug(gmx::formatString("Jcs:%0.3f\n", Jcs));
    }
    return Jcs;
}

void QgenAcm::calcRhs(MsgHandler                 *msg_handler,
                      const std::vector<ActAtom> &atoms,
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
        //! \todo Check this stuff: Delta_Eta * q0
        // rhs_[i] += eta_[nfi]*q0_[i];
        if (bHaveShell_)
        {
            // Core-Shell interaction
            rhs_[i] -= 0.5*calcJcs(msg_handler, nfi, atoms, epsilonr);
            // Hardness of atom multiplied by charge of shell(s)
            for(auto shellId : atoms[nfi].shells())
            {
                rhs_[i] -= eta_[nfi]*q_[shellId];
            }
        }
    }
    if (msg_handler && msg_handler->debug())
    {
        msg_handler->writeDebug(gmx::formatString("QgenACM: qtotal_ = %d qfixed = %g", qtotal_, qfixed));
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

void QgenAcm::solveEEM(MsgHandler *msg_handler)
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
        for (size_t i = 0; i < nelem; i++)
        {
            q_[nonFixed_[i]] = q[i];
        }
        double qtot = std::accumulate(q_.begin(), q_.end(), 0.0);
        
        if (msg_handler && (fabs(qtot - qtotal_) > 1e-2))
        {
            msg_handler->msg(ACTStatus::Warning,
                             gmx::formatString("qtot = %g, it should be %d rhs[%zu] = %g\n",
                                               qtot, qtotal_, nelem, rhs_[nelem]));
        }
    }
    else
    {
        eQGEN_ = eQgen::MATRIXSOLVER;
    }    
}

void QgenAcm::solveSQE(MsgHandler              *msg_handler,
                       const std::vector<Bond> &bonds)
{
    std::vector<double> myq(nonFixed_.size(), 0.0);
    int    info   = 0;
    size_t nbonds = bonds.size();
    if (nbonds > 0)
    {
        MatrixWrapper lhs(nbonds, nbonds);
        std::vector<double> rhs(nbonds, 0.0);

        std::vector<double> delta_chis;
        // First fill the matrix
        for (size_t bij = 0; bij < nbonds; bij++)
        {
            auto ai      = bonds[bij].aI();
            auto aj      = bonds[bij].aJ();
            
            double delta_chi = 0;
            double delta_eta = 0;
            if (bcc_.size() == nbonds)
            {
                delta_chi = bcc_[bij]->find("delta_chi")->second.value() * dchi_factor_[bij];
                delta_eta = bcc_[bij]->find("delta_eta")->second.value();
            }
            if (eQgen::OK != eQGEN_)
            {
                return;
            }
            if (msg_handler && msg_handler->debug())
            {
                msg_handler->writeDebug(gmx::formatString("Delta_chi %10g Delta_Eta %10g\n",
                                                          delta_chi, delta_eta));
            }
            // The bonds use the original numbering, to get to the ACM data
            // we have to map to numbers including shells.
            delta_chis.push_back(delta_chi);
            
            for (size_t bkl = 0; bkl < nbonds; bkl++)
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
        if (msg_handler && msg_handler->debug())
        {
            std::string msg = "Jsqe\n";
            for(size_t i = 0; i < nbonds; i++)
            {
                for(size_t j = 0; j <= i; j++)
                {
                    msg += gmx::formatString(" %7.2f", lhs.get(i, j));
                }
                msg += "\n";
            }
            msg_handler->writeDebug(msg);
        }
        // Now fill the right hand side
        std::vector<double> chi_corr;
        chi_corr.resize(nonFixed_.size(), 0.0);
        for(size_t i = 0; i < nonFixed_.size(); i++)
        {
            int  ai  = static_cast<int>(i);
            auto nfi = nonFixed_[ai];
            chi_corr[i] += chi0_[nfi];
            for (size_t bij = 0; bij < nbonds; bij++)
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
        for (size_t bij = 0; bij < nbonds; bij++)
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
        if (msg_handler && msg_handler->debug())
        {
            msg_handler->writeDebug("rhs: ");
            std::string msg;
            for(auto &p: rhs)
            {
                msg += gmx::formatString(" %8g", p);
            }
            msg_handler->writeDebug(msg);

            msg_handler->writeDebug("pij: ");
            msg.clear();
            for(auto &p: pij)
            {
                msg += gmx::formatString(" %8g", p);
            }
            msg_handler->writeDebug(msg);
        }
        for (size_t bij = 0; bij < nbonds; bij++)
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
    for (size_t i = 0; i < natom_; i++)
    {
        qtot += q_[i];
    }
    if (msg_handler && msg_handler->debug())
    {
        msg_handler->writeDebug("q: ");
        std::string msg;
        for (auto &qq : q_)
        {
            msg += gmx::formatString(" %9g", qq);
        }
        msg_handler->writeDebug(msg);
        if (fabs(qtot - qtotal_) > 1e-2)
        {
            msg_handler->msg(ACTStatus::Warning,
                             gmx::formatString("qtot = %g, it should be %d\n",
                                               qtot, qtotal_));
        }
    }
}

eQgen QgenAcm::generateCharges(MsgHandler                   *msg_handler,
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
    if (msg_handler)
    {
        msg_handler->msg(ACTStatus::Verbose,
                         gmx::formatString("Generating charges for %s using %s algorithm\n",
                                           molname.c_str(), 
                                           chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str()));
    }
    if (eQgen::OK == eQGEN_)
    {
        updateParameters();
        updatePositions(x);
        calcJcc(epsilonr_);
        if (pd->interactionPresent(InteractionType::BONDCORRECTIONS) &&
            pd->chargeGenerationAlgorithm() == ChargeGenerationAlgorithm::SQE)
        {
            solveSQE(msg_handler, bonds);
        }
        else
        {
            calcRhs(msg_handler, *atoms, epsilonr_);
            solveEEM(msg_handler);
        }
        if (eQgen::OK == eQGEN_)
        {
            copyChargesToAtoms(atoms);
            dump(msg_handler, atoms);
        }
    }
    return eQGEN_;
}

} // namespace
