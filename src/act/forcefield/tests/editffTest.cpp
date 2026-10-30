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
 * Tests of the edit_ff module.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "act/alexandria/alex_modules.h"
#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class EditFFTest : public gmx::test::CommandLineTestBase
{
protected:
    //! \return the mutability of chi for an ACM type
    static Mutability chiMutability(const ForceField &pd, const std::string &acmtype)
    {
        auto fs = pd.findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
        return fs.findParameterTypeConst(Identifier(acmtype), "chi").mutability();
    }
};

// Test: The mutability of a parameter of an interaction can be changed for
// several identifiers at once
TEST_F(EditFFTest, MutabilityOfSeveralIdentifiers)
{
    auto output = fileManager().getTemporaryFilePath("edited.xml");
    gmx::test::CommandLine cmd;
    cmd.append("edit_ff");
    cmd.addOption("-ff", fileManager().getInputFilePath("ACM-g.xml"));
    cmd.addOption("-o", output);
    cmd.addOption("-g", fileManager().getTemporaryFilePath("edit_ff.log"));
    cmd.addOption("-p", "chi");
    cmd.addOption("-a", "h_z f_z");
    cmd.addOption("-mut", "Fixed");
    ASSERT_EQ(0, edit_ff(cmd.argc(), cmd.argv()));
    // Without changes no file is written
    ASSERT_TRUE(std::filesystem::exists(output));
    ForceField pd;
    readForceField(output, &pd);
    EXPECT_EQ(Mutability::Fixed, chiMutability(pd, "h_z"));
    EXPECT_EQ(Mutability::Fixed, chiMutability(pd, "f_z"));
    EXPECT_EQ(Mutability::Bounded, chiMutability(pd, "c3_z"));
}

} // namespace

} // namespace alexandria
