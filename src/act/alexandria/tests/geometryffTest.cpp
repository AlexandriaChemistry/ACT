/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2023,2026
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
 * Tests of the geometry_ff module.
 */

#include "actpre.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "act/alexandria/alex_modules.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class GeometryFFTest : public gmx::test::CommandLineTestBase
{
};

//! \return the parameters of a bond, given as text
static const ForceFieldParameterMap *bond(const ForceField &pd, const std::string &id)
{
    for (const auto &p : pd.findForcesConst(InteractionType::BONDS).parametersConst())
    {
        if (p.first.id() == id)
        {
            return &p.second;
        }
    }
    return nullptr;
}

// Test: The bond lengths of the compounds, with their standard deviation and
// number, are the new parameters, other bonds are removed, and the energies
// of harmonic bonds are negative
TEST_F(GeometryFFTest, AlcoholAndAcetone)
{
    const auto output = fileManager().getTemporaryFilePath("geometry.xml");
    gmx::test::CommandLine cmd;
    cmd.append("geometry_ff");
    cmd.addOption("-mp", fileManager().getInputFilePath("molprop/molprop.xml"));
    cmd.addOption("-ff", fileManager().getInputFilePath("forcefields/ACS-g1.xml"));
    cmd.addOption("-o", output);
    cmd.addOption("-g", fileManager().getTemporaryFilePath("geometry.log"));
    ASSERT_EQ(0, geometry_ff(cmd.argc(), cmd.argv()));
    ForceField pd;
    readForceField(output, &pd);
    auto ch = bond(pd, "c3_b~h_b");
    ASSERT_NE(nullptr, ch);
    // Nine C-H bonds in each of two structures of butanol, six in acetone
    auto &length = ch->find("bondlength")->second;
    EXPECT_EQ(24u, length.ntrain());
    EXPECT_NEAR(109.4, length.value(), 0.05);
    EXPECT_NEAR(0.5, length.uncertainty(), 0.05);
    EXPECT_LE(length.minimum(), length.value());
    EXPECT_GE(length.maximum(), length.value());
    auto &energy = ch->find("bondenergy")->second;
    EXPECT_LT(energy.value(), 0);
    EXPECT_FALSE(energy.nonNegative());
    EXPECT_LE(energy.minimum(), energy.value());
    EXPECT_GE(energy.maximum(), energy.value());
    // The bonds are made from the compounds only
    EXPECT_EQ(nullptr, bond(pd, "c3_b~cl_b"));
}

} // namespace

} // namespace alexandria
