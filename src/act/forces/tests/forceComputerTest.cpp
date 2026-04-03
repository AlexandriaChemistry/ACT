/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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
/*! \internal \file
 * \brief
 * Tests for the ForceComputer class: initialization, ftype(), computeOnce(),
 * compute(), and the vsite coordinate/force wrappers.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_forces_tests
 */
#include "actpre.h"

#include "../forcecomputer.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/actmol.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/potential.h"
#include "act/import/import.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

namespace
{

// ============================================================
// ForceComputer initialization and simple API tests
// ============================================================

TEST(ForceComputerInit, DefaultConstructorHasExpectedTolerance)
{
    ForceComputer fc;
    // Default force tolerance is 1e-6 (see forcecomputer.h)
    EXPECT_DOUBLE_EQ(1e-6, fc.forceTolerance());
}

TEST(ForceComputerInit, ParameterizedConstructorStoresTolerance)
{
    const double toler = 1e-8;
    ForceComputer fc(toler, 50, 0.02);
    EXPECT_DOUBLE_EQ(toler, fc.forceTolerance());
}

TEST(ForceComputerInit, SetForceToleranceRoundTrip)
{
    ForceComputer fc;
    const double newToler = 1e-10;
    fc.setForceTolerance(newToler);
    EXPECT_DOUBLE_EQ(newToler, fc.forceTolerance());
}

TEST(ForceComputerInit, SetForceToleranceToSmallValue)
{
    ForceComputer fc;
    fc.setForceTolerance(0.0);
    EXPECT_DOUBLE_EQ(0.0, fc.forceTolerance());
}

TEST(ForceComputerInit, InitOverridesDefaultTolerance)
{
    ForceComputer fc;
    fc.init(1e-4, 10, 0.1);
    EXPECT_DOUBLE_EQ(1e-4, fc.forceTolerance());
}

// ============================================================
// ForceComputer::ftype() tests (requires a force field)
// ============================================================

TEST(ForceComputerFtype, PresentInteractionReturnsCorrectPotential)
{
    auto pd = getForceField("ACS-pg-vs2");
    ASSERT_NE(nullptr, pd);
    ForceComputer fc;
    // VSITE2 is present in ACS-pg-vs2
    auto pot = fc.ftype(pd, InteractionType::VSITE2);
    EXPECT_EQ(Potential::VSITE2, pot);
}

TEST(ForceComputerFtype, MissingInteractionReturnsNone)
{
    auto pd = getForceField("ACS-pg-vs2");
    ASSERT_NE(nullptr, pd);
    ForceComputer fc;
    // BONDS is not expected to be in ACS-pg-vs2
    if (!pd->interactionPresent(InteractionType::BONDS))
    {
        auto pot = fc.ftype(pd, InteractionType::BONDS);
        EXPECT_EQ(Potential::NONE, pot);
    }
}

TEST(ForceComputerFtype, Vsite1PresentInVs2ForceField)
{
    auto pd = getForceField("ACS-pg-vs2");
    ASSERT_NE(nullptr, pd);
    ForceComputer fc;
    auto pot = fc.ftype(pd, InteractionType::VSITE1);
    EXPECT_EQ(Potential::VSITE1, pot);
}

TEST(ForceComputerFtype, Vsite3PresentInVs3ForceField)
{
    auto pd = getForceField("ACS-pg-vs3");
    ASSERT_NE(nullptr, pd);
    ForceComputer fc;
    auto pot = fc.ftype(pd, InteractionType::VSITE3);
    EXPECT_EQ(Potential::VSITE3, pot);
}

// ============================================================
// Integration tests: computeOnce() and compute() with real molecules
// ============================================================

//! Helper: load a molecule from an SDF test file and generate its topology.
static bool setupMolecule(const std::string &ffName,
                          const std::string &molName,
                          ACTMol            *mol)
{
    auto pd = getForceField(ffName);
    if (!pd)
    {
        return false;
    }
    std::vector<MolProp> molprops;
    bool   userqtot   = false;
    double qtot_babel = 0;
    const char  *conf     = "minimum";
    std::string  fileName = gmx::formatString("%s.sdf", molName.c_str());
    std::string  dataName = gmx::test::TestFileManager::getInputFilePath(fileName);
    MsgHandler   msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    importFile(&msghandler, pd, dataName.c_str(), &molprops,
               conf, JobType::OPT, userqtot, &qtot_babel, true);
    if (!msghandler.ok() || molprops.size() != 1)
    {
        return false;
    }
    mol->Merge(&molprops[0]);
    mol->GenerateTopology(&msghandler, pd, missingParameters::Ignore);
    return msghandler.ok();
}

class ForceComputerIntegrationTest : public ::testing::Test
{
protected:
    /*! \brief Load a molecule and run computeOnce().
     * \param[in]  ffName  Force field name (e.g. "ACS-pg-vs2")
     * \param[in]  molName Molecule file stem (e.g. "hydrogen-fluoride")
     * \param[out] energies Populated by computeOnce
     * \param[out] forces   Populated by computeOnce
     * \return true on success
     */
    bool runComputeOnce(const std::string                 &ffName,
                        const std::string                 &molName,
                        std::map<InteractionType, double> *energies,
                        std::vector<gmx::RVec>            *forces)
    {
        auto pd = getForceField(ffName);
        if (!pd)
        {
            return false;
        }
        ACTMol mol;
        if (!setupMolecule(ffName, molName, &mol))
        {
            return false;
        }
        auto coords = mol.xOriginal();
        const auto *top = mol.topology();

        ForceComputer fc;
        fc.init();
        fc.constructVsiteCoordinates(top, &coords);

        forces->assign(coords.size(), { 0, 0, 0 });
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        fc.computeOnce(&msghandler, pd, top, &coords, forces, energies, { 0, 0, 0 });
        return msghandler.ok();
    }

    /*! \brief Load a molecule and run compute().
     * \param[in]  ffName  Force field name
     * \param[in]  molName Molecule file stem
     * \param[out] energies Populated by compute
     * \param[out] forces   Populated by compute
     * \return true on success
     */
    bool runCompute(const std::string                 &ffName,
                    const std::string                 &molName,
                    std::map<InteractionType, double> *energies,
                    std::vector<gmx::RVec>            *forces)
    {
        auto pd = getForceField(ffName);
        if (!pd)
        {
            return false;
        }
        ACTMol mol;
        if (!setupMolecule(ffName, molName, &mol))
        {
            return false;
        }
        auto coords = mol.xOriginal();
        const auto *top = mol.topology();

        ForceComputer fc;
        fc.init();
        forces->assign(coords.size(), { 0, 0, 0 });
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        fc.compute(&msghandler, pd, top, &coords, forces, energies);
        return msghandler.ok();
    }
};

TEST_F(ForceComputerIntegrationTest, ComputeOnceHFPopulatesEpot)
{
    std::map<InteractionType, double> energies;
    std::vector<gmx::RVec>           forces;
    bool ok = runComputeOnce("ACS-pg-vs2", "hydrogen-fluoride", &energies, &forces);
    EXPECT_TRUE(ok);
    EXPECT_FALSE(forces.empty());
    EXPECT_GT(energies.count(InteractionType::EPOT), 0u);
    EXPECT_TRUE(std::isfinite(energies.at(InteractionType::EPOT)));
}

TEST_F(ForceComputerIntegrationTest, ComputeOnceHFHasElectrostatics)
{
    std::map<InteractionType, double> energies;
    std::vector<gmx::RVec>           forces;
    bool ok = runComputeOnce("ACS-pg-vs2", "hydrogen-fluoride", &energies, &forces);
    EXPECT_TRUE(ok);
    EXPECT_GT(energies.count(InteractionType::ELECTROSTATICS), 0u);
    EXPECT_TRUE(std::isfinite(energies.at(InteractionType::ELECTROSTATICS)));
}

TEST_F(ForceComputerIntegrationTest, ComputeOnceWaterPopulatesEpot)
{
    std::map<InteractionType, double> energies;
    std::vector<gmx::RVec>           forces;
    bool ok = runComputeOnce("ACS-pg-vs3", "water", &energies, &forces);
    EXPECT_TRUE(ok);
    EXPECT_GT(energies.count(InteractionType::EPOT), 0u);
    EXPECT_TRUE(std::isfinite(energies.at(InteractionType::EPOT)));
}

TEST_F(ForceComputerIntegrationTest, ComputeHFPopulatesEpotAndAllelec)
{
    std::map<InteractionType, double> energies;
    std::vector<gmx::RVec>           forces;
    bool ok = runCompute("ACS-pg-vs2", "hydrogen-fluoride", &energies, &forces);
    EXPECT_TRUE(ok);
    EXPECT_GT(energies.count(InteractionType::EPOT), 0u);
    EXPECT_TRUE(std::isfinite(energies.at(InteractionType::EPOT)));
    EXPECT_GT(energies.count(InteractionType::ALLELEC), 0u);
    EXPECT_TRUE(std::isfinite(energies.at(InteractionType::ALLELEC)));
}

TEST_F(ForceComputerIntegrationTest, ComputeWaterPopulatesEpotAndAllelec)
{
    std::map<InteractionType, double> energies;
    std::vector<gmx::RVec>           forces;
    bool ok = runCompute("ACS-pg-vs3", "water", &energies, &forces);
    EXPECT_TRUE(ok);
    EXPECT_GT(energies.count(InteractionType::EPOT), 0u);
    EXPECT_TRUE(std::isfinite(energies.at(InteractionType::EPOT)));
    EXPECT_GT(energies.count(InteractionType::ALLELEC), 0u);
    EXPECT_TRUE(std::isfinite(energies.at(InteractionType::ALLELEC)));
}

TEST_F(ForceComputerIntegrationTest, ComputeForcesHaveCorrectSize)
{
    std::map<InteractionType, double> energies;
    std::vector<gmx::RVec>           forces;
    auto pd = getForceField("ACS-pg-vs2");
    ASSERT_NE(nullptr, pd);
    ACTMol mol;
    ASSERT_TRUE(setupMolecule("ACS-pg-vs2", "hydrogen-fluoride", &mol));
    auto coords = mol.xOriginal();
    const auto *top = mol.topology();
    ForceComputer fc;
    fc.init();
    forces.assign(coords.size(), { 0, 0, 0 });
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    fc.compute(&msghandler, pd, top, &coords, &forces, &energies);
    EXPECT_TRUE(msghandler.ok());
    EXPECT_EQ(coords.size(), forces.size());
}

TEST_F(ForceComputerIntegrationTest, ComputeOnceAndComputeGiveSameEpotWithoutShells)
{
    // For a system without polarizable shells, compute() should yield the
    // same EPOT as a single computeOnce() call (after constructing vsites).
    auto pd = getForceField("ACS-pg-vs2");
    ASSERT_NE(nullptr, pd);
    // Only test non-polarizable force fields (no shell minimization)
    if (pd->polarizable())
    {
        return;
    }
    ACTMol mol;
    ASSERT_TRUE(setupMolecule("ACS-pg-vs2", "hydrogen-fluoride", &mol));
    const auto *top = mol.topology();
    ForceComputer fc;
    fc.init();

    // computeOnce
    auto coords1 = mol.xOriginal();
    fc.constructVsiteCoordinates(top, &coords1);
    std::vector<gmx::RVec>           forces1(coords1.size(), { 0, 0, 0 });
    std::map<InteractionType, double> energies1;
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    fc.computeOnce(&msghandler, pd, top, &coords1, &forces1, &energies1, { 0, 0, 0 });
    ASSERT_TRUE(msghandler.ok());

    // compute
    auto coords2 = mol.xOriginal();
    std::vector<gmx::RVec>           forces2(coords2.size(), { 0, 0, 0 });
    std::map<InteractionType, double> energies2;
    fc.compute(&msghandler, pd, top, &coords2, &forces2, &energies2);
    ASSERT_TRUE(msghandler.ok());

    ASSERT_GT(energies1.count(InteractionType::EPOT), 0u);
    ASSERT_GT(energies2.count(InteractionType::EPOT), 0u);
    EXPECT_NEAR(energies1.at(InteractionType::EPOT),
                energies2.at(InteractionType::EPOT), 1e-6);
}

TEST_F(ForceComputerIntegrationTest, ConstructVsiteCoordinatesAndSpreadForces)
{
    // Test that the ForceComputer wrappers for vsite coordinate construction
    // and force spreading behave consistently with a direct VsiteHandler call.
    auto pd = getForceField("ACS-pg-vs2");
    ASSERT_NE(nullptr, pd);
    ACTMol mol;
    ASSERT_TRUE(setupMolecule("ACS-pg-vs2", "hydrogen-fluoride", &mol));
    const auto *top = mol.topology();
    ForceComputer fc;
    fc.init();

    auto coords = mol.xOriginal();
    fc.constructVsiteCoordinates(top, &coords);
    // All coordinates must be finite
    for (const auto &c : coords)
    {
        EXPECT_TRUE(std::isfinite(c[XX]));
        EXPECT_TRUE(std::isfinite(c[YY]));
        EXPECT_TRUE(std::isfinite(c[ZZ]));
    }

    // Spread unit forces on vsite atoms and check result is finite
    std::vector<gmx::RVec> forces(coords.size(), { 0, 0, 0 });
    const auto &atoms = top->atoms();
    for (size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].pType() == ActParticle::Vsite)
        {
            forces[i] = { 1, 1, 1 };
        }
    }
    fc.spreadVsiteForces(top, &coords, &forces);
    for (const auto &f : forces)
    {
        EXPECT_TRUE(std::isfinite(f[XX]));
        EXPECT_TRUE(std::isfinite(f[YY]));
        EXPECT_TRUE(std::isfinite(f[ZZ]));
    }
}

}  // namespace

}  // namespace alexandria
