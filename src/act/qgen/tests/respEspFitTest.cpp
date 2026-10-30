/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2026
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
 * Tests of fitting charges to the electrostatic potential with QgenResp,
 * using the quantum chemical potentials of the molecules in
 * src/act/molprop/tests/molprop.xml.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forces/forcecomputer.h"
#include "act/molprop/molprop_xml.h"
#include "act/qgen/qgen_resp.h"
#include "act/topology/actmol.h"
#include "act/topology/symmetrize_charges.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/*! \brief Tolerance for charges that should be the same
 * The total charge and the symmetry of the charges are imposed with a
 * large weight in the least squares fit, not exactly, and the total charge
 * is corrected afterwards, see QgenResp::optimizeCharges.
 */
const double symmetryTolerance = 1e-5;

//! A molecule with the charges fitted to its electrostatic potential
struct EspFit
{
    //! The molecule
    ACTMol            actmol;
    //! The force field
    ForceField       *pd       = nullptr;
    //! The RESP data that the charges were fitted with
    QgenResp         *qresp    = nullptr;
    //! Atoms with the same charge, each refers to the first of them
    std::vector<int>  symmetric;
    //! Relative dielectric constant
    double            epsilonr = 1;
};

/*! \brief Class to test fitting charges to the electrostatic potential
 */
class RespEspFitTest : public gmx::test::CommandLineTestBase
{
protected:
    //! Checking data structure
    gmx::test::TestReferenceChecker checker_;

    //! Set the tolerance
    RespEspFitTest() : checker_(this->rootChecker())
    {
        checker_.setDefaultTolerance(gmx::test::absoluteTolerance(1e-5));
    }

    /*! \brief Read a molecule and fit its charges to its electrostatic potential
     * \param[in]  ffName  Force field file in this directory, without .xml
     * \param[in]  molname Name of the molecule in the molprop file
     * \param[out] fit     The molecule and its fit
     * \return whether all went well
     */
    bool fitCharges(const std::string &ffName,
                    const std::string &molname,
                    EspFit            *fit)
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        fit->pd = getForceField(ffName);
        std::vector<MolProp> mps;
        auto mpFile = fileManager().getInputFilePath("../../molprop/tests/molprop.xml");
        MolPropRead(&msghandler, mpFile.c_str(), &mps);
        bool found = false;
        for(auto &mp : mps)
        {
            if (mp.getMolname() == molname)
            {
                fit->actmol.Merge(&mp);
                found = true;
            }
        }
        EXPECT_TRUE(found) << molname;
        if (!found)
        {
            return false;
        }
        auto &actmol = fit->actmol;
        actmol.GenerateTopology(&msghandler, fit->pd, missingParameters::Error);
        EXPECT_TRUE(msghandler.ok());
        if (!msghandler.ok())
        {
            return false;
        }
        std::map<MolPropObservable, iqmType> iqm = {
            { MolPropObservable::POTENTIAL, iqmType::QM }
        };
        actmol.getExpProps(&msghandler, fit->pd, iqm, 0, 100);
        EXPECT_TRUE(msghandler.ok());
        if (!ffOption(*fit->pd, InteractionType::ELECTROSTATICS, "epsilonr", &fit->epsilonr))
        {
            fit->epsilonr = 1;
        }
        ForceComputer          forceComp(1e-8, 100, fit->epsilonr);
        std::vector<gmx::RVec> coords = actmol.xOriginal();
        std::vector<gmx::RVec> forces(actmol.atomsConst().size());
        actmol.generateCharges(&msghandler, fit->pd, &forceComp, ChargeGenerationAlgorithm::ESP,
                               &coords, &forces, true);
        EXPECT_TRUE(msghandler.ok());
        // The charges are fitted to the first potential
        for(auto &qp : *actmol.qProps())
        {
            if (qp.qgenResp()->nEsp() > 0 && nullptr == fit->qresp)
            {
                fit->qresp = qp.qgenResp();
            }
        }
        EXPECT_NE(nullptr, fit->qresp);
        get_symmetrized_charges(actmol.topology(), fit->pd, nullptr, &fit->symmetric);
        return msghandler.ok() && nullptr != fit->qresp;
    }

    /*! \brief Compute the root mean square deviation from the ESP with other charges
     * \param[in] fit     The fit
     * \param[in] charges The charges of the atoms, those of the shells are fixed
     * \return the RMS deviation in Hartree/e
     */
    static double rmsWithCharges(const EspFit              &fit,
                                 const std::vector<double> &charges)
    {
        MsgHandler           msghandler;
        QgenResp             trial(*fit.qresp);
        std::vector<ActAtom> atoms = fit.actmol.atomsConst();
        for(size_t i = 0; i < atoms.size(); i++)
        {
            atoms[i].setCharge(charges[i]);
        }
        trial.updateAtomCharges(&msghandler, atoms);
        trial.calcPot(&msghandler, fit.epsilonr);
        real rrms, cosangle, mae, mse;
        return trial.getStatistics(&msghandler, &rrms, &cosangle, &mae, &mse);
    }

    /*! \brief Fit the charges and check them and the quality of the fit
     * \param[in] ffName  Force field file in this directory, without .xml
     * \param[in] molname Name of the molecule in the molprop file
     */
    void testFit(const std::string &ffName,
                 const std::string &molname)
    {
        EspFit fit;
        if (!fitCharges(ffName, molname, &fit))
        {
            return;
        }
        const auto &atoms = fit.actmol.atomsConst();
        ASSERT_EQ(atoms.size(), fit.qresp->natoms());
        ASSERT_EQ(atoms.size(), fit.symmetric.size());
        // The charges of the molecule are the fitted ones
        std::vector<double> charges;
        double              qtotal = 0;
        for(size_t i = 0; i < atoms.size(); i++)
        {
            EXPECT_EQ(fit.qresp->getCharge(i), atoms[i].charge()) << "atom " << i+1;
            charges.push_back(atoms[i].charge());
            qtotal += atoms[i].charge();
        }
        checker_.checkSequence(charges.begin(), charges.end(), "Charges");
        // The molecules are neutral, including the charges of the shells
        EXPECT_NEAR(0, qtotal, 1e-10);
        // Symmetric atoms have the same charge
        for(size_t i = 0; i < atoms.size(); i++)
        {
            EXPECT_NEAR(charges[fit.symmetric[i]], charges[i], symmetryTolerance)
                << "atom " << i+1 << " symmetric to atom " << fit.symmetric[i]+1;
        }
        // The potential of the fitted charges is close to the quantum chemical one
        MsgHandler msghandler;
        real       rrms, cosangle, mae, mse;
        double     rms = fit.qresp->getStatistics(&msghandler, &rrms, &cosangle, &mae, &mse);
        EXPECT_GT(cosangle, 0.95);
        EXPECT_LT(rrms, 0.3);
        checker_.checkReal(rms, "RMS");
        checker_.checkReal(rrms, "RRMS");
        checker_.checkReal(cosangle, "CosAngle");
        EXPECT_NEAR(rms, rmsWithCharges(fit, charges), 1e-12);
        // The fitted charges have the smallest deviation of all charges
        // with the same total charge and symmetry. Move charge from
        // each group of symmetric atoms to each other group.
        std::map<int, int> groupSize;
        for(size_t i = 0; i < atoms.size(); i++)
        {
            if (ActParticle::Atom == atoms[i].pType())
            {
                groupSize[fit.symmetric[i]] += 1;
            }
        }
        EXPECT_GT(groupSize.size(), 1u);
        double dq = 0.01;
        for(const auto &from : groupSize)
        {
            for(const auto &to : groupSize)
            {
                if (from.first == to.first)
                {
                    continue;
                }
                auto trial = charges;
                for(size_t i = 0; i < atoms.size(); i++)
                {
                    if (ActParticle::Atom != atoms[i].pType())
                    {
                        continue;
                    }
                    if (fit.symmetric[i] == from.first)
                    {
                        trial[i] -= dq/from.second;
                    }
                    else if (fit.symmetric[i] == to.first)
                    {
                        trial[i] += dq/to.second;
                    }
                }
                EXPECT_GT(rmsWithCharges(fit, trial), rms)
                    << "moving charge from atom " << from.first+1 << " to atom " << to.first+1;
            }
        }
    }

    /*! \brief Fit charges to the potential of the fitted charges
     * The potential at the points of the quantum chemical one is computed
     * from the fitted charges, then charges are fitted to it again, with
     * the given total charge. For the total charge of the molecule the
     * fitted charges are found again.
     * \param[in] ffName  Force field file in this directory, without .xml
     * \param[in] molname Name of the molecule in the molprop file
     * \param[in] qtotal  The total charge for the new fit
     */
    void testRefit(const std::string &ffName,
                   const std::string &molname,
                   int                qtotal)
    {
        EspFit fit;
        if (!fitCharges(ffName, molname, &fit))
        {
            return;
        }
        MsgHandler  msghandler;
        const auto &atoms = fit.actmol.atomsConst();
        QgenResp    refit;
        refit.setChargeDistributionType(fit.qresp->chargeType());
        refit.setAtomInfo(&msghandler, atoms, fit.pd, qtotal);
        refit.updateAtomCoords(fit.qresp->coords());
        for(const auto &ep : fit.qresp->espPoints())
        {
            refit.addEspPoint(ep.esp()[XX], ep.esp()[YY], ep.esp()[ZZ], ep.vCalc());
        }
        refit.setAtomSymmetry(fit.symmetric);
        refit.optimizeCharges(&msghandler, fit.epsilonr);
        refit.calcPot(&msghandler, fit.epsilonr);
        EXPECT_TRUE(msghandler.ok());
        double              qsum = 0;
        std::vector<double> charges;
        for(size_t i = 0; i < atoms.size(); i++)
        {
            charges.push_back(refit.getCharge(i));
            qsum += refit.getCharge(i);
            if (fit.actmol.totalCharge() == qtotal)
            {
                EXPECT_NEAR(fit.qresp->getCharge(i), refit.getCharge(i), symmetryTolerance) << "atom " << i+1;
            }
            EXPECT_NEAR(refit.getCharge(fit.symmetric[i]), refit.getCharge(i), symmetryTolerance)
                << "atom " << i+1 << " symmetric to atom " << fit.symmetric[i]+1;
        }
        EXPECT_NEAR(qtotal, qsum, 1e-10);
        checker_.checkSequence(charges.begin(), charges.end(), "Charges");
        real   rrms, cosangle, mae, mse;
        double rms = refit.getStatistics(&msghandler, &rrms, &cosangle, &mae, &mse);
        if (fit.actmol.totalCharge() == qtotal)
        {
            // Nothing is left of the deviation, apart from what
            // the small differences between symmetric charges give
            EXPECT_LT(rms, 1e-6);
            EXPECT_NEAR(1, cosangle, 1e-8);
        }
        else
        {
            // An extra charge cannot be spread such that the potential stays the same
            EXPECT_GT(rms, 1e-3);
        }
    }
};

TEST_F (RespEspFitTest, PointButanol)
{
    testFit("ESP-p", "1-butanol");
}

TEST_F (RespEspFitTest, PointAcetone)
{
    testFit("ESP-p", "acetone");
}

TEST_F (RespEspFitTest, PointWater)
{
    testFit("ESP-p", "water");
}

TEST_F (RespEspFitTest, PolarizableGaussianButanol)
{
    testFit("ESP-pg", "1-butanol");
}

TEST_F (RespEspFitTest, PolarizableGaussianAcetone)
{
    testFit("ESP-pg", "acetone");
}

TEST_F (RespEspFitTest, PolarizableGaussianWater)
{
    testFit("ESP-pg", "water");
}

TEST_F (RespEspFitTest, RefitPointButanol)
{
    testRefit("ESP-p", "1-butanol", 0);
}

TEST_F (RespEspFitTest, RefitPolarizableGaussianButanol)
{
    testRefit("ESP-pg", "1-butanol", 0);
}

TEST_F (RespEspFitTest, RefitPolarizableGaussianAcetone)
{
    testRefit("ESP-pg", "acetone", 0);
}

TEST_F (RespEspFitTest, RefitChargedPointWater)
{
    testRefit("ESP-p", "water", 1);
}

} // namespace

} // namespace alexandria
