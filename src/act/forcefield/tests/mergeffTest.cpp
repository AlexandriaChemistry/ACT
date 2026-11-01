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
 * Tests of the merge_ff module.
 */

#include "actpre.h"

#include <filesystem>
#include <string>
#include <vector>

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

class MergeFFTest : public gmx::test::CommandLineTestBase
{
protected:
    //! The force field the tests start from
    std::string input_;

    MergeFFTest() : input_(fileManager().getInputFilePath("ACM-g.xml")) {}

    //! \return a parameter of an ACM type, such as chi or eta
    static const ForceFieldParameter &acm(const ForceField &pd, const std::string &acmtype,
                                          const std::string &type)
    {
        auto &fs = pd.findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
        return fs.findParameterTypeConst(Identifier(acmtype), type);
    }

    /*! \brief A copy of the input with chi of h_z 4 and of f_z 6
     * \param[in] name   File name of the copy
     * \param[in] remove Identifiers whose EEM parameters are left out
     * \return the path of the copy
     */
    std::string modified(const std::string &name, const std::vector<std::string> &remove = {})
    {
        ForceField pd;
        readForceField(input_, &pd);
        auto eem = pd.findForces(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
        eem->findParameterType(Identifier("h_z"), "chi")->setValue(4);
        eem->findParameterType(Identifier("f_z"), "chi")->setValue(6);
        for (const auto &id : remove)
        {
            eem->parameters()->erase(Identifier(id));
        }
        auto fileName = fileManager().getTemporaryFilePath(name);
        writeForceField(fileName, &pd, false);
        return fileName;
    }

    //! \return the merged force field of the files, with the other options
    ForceField merge(const std::vector<std::string> &files, const std::vector<std::pair<std::string, std::string>> &options)
    {
        auto output = fileManager().getTemporaryFilePath("merged.xml");
        gmx::test::CommandLine cmd;
        cmd.append("merge_ff");
        cmd.addOption("-ff");
        for (const auto &f : files)
        {
            cmd.append(f);
        }
        cmd.addOption("-o", output);
        for (const auto &o : options)
        {
            cmd.addOption(o.first.c_str(), o.second);
        }
        EXPECT_EQ(0, merge_ff(cmd.argc(), cmd.argv()));
        ForceField pd;
        EXPECT_TRUE(std::filesystem::exists(output));
        readForceField(output, &pd);
        return pd;
    }
};

// Test: The values of a parameter are averaged, with their standard
// deviation as the uncertainty, and the others are those of the first file
TEST_F(MergeFFTest, Average)
{
    auto pd  = merge({ input_, modified("second.xml") }, { { "-merge", "chi" } });
    auto ref = ForceField();
    readForceField(input_, &ref);
    EXPECT_DOUBLE_EQ(3.75, acm(pd, "h_z", "chi").value());
    EXPECT_DOUBLE_EQ(0.25, acm(pd, "h_z", "chi").uncertainty());
    EXPECT_DOUBLE_EQ((acm(ref, "f_z", "chi").value() + 6)/2, acm(pd, "f_z", "chi").value());
    EXPECT_DOUBLE_EQ(acm(ref, "c3_z", "chi").value(), acm(pd, "c3_z", "chi").value());
    EXPECT_DOUBLE_EQ(0, acm(pd, "c3_z", "chi").uncertainty());
    EXPECT_DOUBLE_EQ(acm(ref, "h_z", "eta").value(), acm(pd, "h_z", "eta").value());
}

// Test: New bounds of the average plus and minus a number of standard
// deviations
TEST_F(MergeFFTest, Limits)
{
    auto pd = merge({ input_, modified("second.xml") }, { { "-merge", "chi" }, { "-limits", "2" } });
    auto p = acm(pd, "h_z", "chi");
    EXPECT_DOUBLE_EQ(3.25, p.minimum());
    EXPECT_DOUBLE_EQ(4.25, p.maximum());
}

// Test: The parameters are matched by their identifiers, such that a force
// field that lacks some still gives the averages of the others
TEST_F(MergeFFTest, MissingIdentifier)
{
    auto ref = ForceField();
    readForceField(input_, &ref);
    // An identifier that comes before h_z and f_z in the parameters
    auto &eem = ref.findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION).parametersConst();
    ASSERT_TRUE(eem.find(Identifier("br_z")) != eem.end());
    ASSERT_TRUE(Identifier("br_z") < Identifier("f_z"));
    auto pd = merge({ input_, modified("second.xml", { "br_z" }) }, { { "-merge", "chi" } });
    EXPECT_DOUBLE_EQ(3.75, acm(pd, "h_z", "chi").value());
    EXPECT_DOUBLE_EQ((acm(ref, "f_z", "chi").value() + 6)/2, acm(pd, "f_z", "chi").value());
    // Only in the first file
    EXPECT_DOUBLE_EQ(acm(ref, "br_z", "chi").value(), acm(pd, "br_z", "chi").value());
}

// Test: Without a list all parameters are merged, and a table is written
TEST_F(MergeFFTest, AllParametersAndTable)
{
    auto latex = fileManager().getTemporaryFilePath("table.tex");
    auto pd    = merge({ input_, modified("second.xml") }, { { "-latex", latex } });
    EXPECT_DOUBLE_EQ(3.75, acm(pd, "h_z", "chi").value());
    EXPECT_GT(std::filesystem::file_size(latex), 0u);
}

} // namespace

} // namespace alexandria
