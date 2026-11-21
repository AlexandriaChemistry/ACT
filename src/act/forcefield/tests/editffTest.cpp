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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <iterator>
#include <map>
#include <nlohmann/json.hpp>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/forcefield/combinationkernels.h"
#include "act/forcefield/edit_ff.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_validation.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/particletype.h"
#include "act/forcefield/potential.h"
#include "act/testutils/output_checks.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! The options of edit_ff besides the input force field
using Options = std::vector<std::pair<std::string, std::string>>;

class EditFFTest : public gmx::test::CommandLineTestBase
{
protected:
    //! The force field the tests start from
    std::string input_;

    EditFFTest() : input_(fileManager().getInputFilePath("forcefields/ACM-g1.xml")) {}

    //! \return the mutability of chi for an ACM type
    static Mutability chiMutability(const ForceField &pd, const std::string &acmtype)
    {
        const auto& fs = pd.findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
        return fs.findParameterTypeConst(Identifier(acmtype), "chi").mutability();
    }

    //! \return a parameter of an ACM type, such as chi or eta
    static const ForceFieldParameter &acm(const ForceField &pd, const std::string &acmtype,
                                          const std::string &type)
    {
        auto &fs = pd.findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
        return fs.findParameterTypeConst(Identifier(acmtype), type);
    }

    /*! \brief Run edit_ff
     * \param[in] ff      The force field to edit
     * \param[in] options The other options, an empty value for a flag
     * \param[in] output  Name of the output file, written with -o if not empty
     * \return the path of the output file
     */
    std::string run(const std::string &ff, const Options &options, const std::string &output = "edited.xml")
    {
        gmx::test::CommandLine cmd;
        cmd.append("edit_ff");
        cmd.addOption("-ff", ff);
        std::string out;
        if (!output.empty())
        {
            out = fileManager().getTemporaryFilePath(output);
            cmd.addOption("-o", out);
        }
        cmd.addOption("-g", fileManager().getTemporaryFilePath("edit_ff.log"));
        for (const auto &o : options)
        {
            if (o.second.empty())
            {
                cmd.append(o.first);
            }
            else
            {
                cmd.addOption(o.first.c_str(), o.second);
            }
        }
        EXPECT_EQ(0, edit_ff(cmd.argc(), cmd.argv()));
        return out;
    }

    //! \return the force field in a file, which has to exist
    static ForceField read(const std::string &fileName)
    {
        ForceField pd;
        EXPECT_TRUE(std::filesystem::exists(fileName)) << fileName;
        if (std::filesystem::exists(fileName))
        {
            readForceField(fileName, &pd);
        }
        return pd;
    }

    //! Expect a parameter value, which is stored as real
    static void expectRealEq(double expected, double actual)
    {
        if (GMX_DOUBLE)
        {
            EXPECT_DOUBLE_EQ(expected, actual);
        }
        else
        {
            EXPECT_FLOAT_EQ(expected, actual);
        }
    }

    //! Expect the minimum, value and maximum of a parameter
    static void expectRange(const ForceFieldParameter &p, double minimum, double value, double maximum)
    {
        expectRealEq(minimum, p.minimum());
        expectRealEq(value, p.value());
        expectRealEq(maximum, p.maximum());
    }
};

// Test: The mutability of a parameter of an interaction can be changed for
// several identifiers at once
/*! \brief Run edit_ff and compare what it changed with the reference data
 *
 * These tests fix the results of the tool with the given options, also
 * that nothing else changes, such that a change in how the options reach
 * the work shows.
 */
class EditFFResultsTest : public EditFFTest
{
protected:
    /*! \brief Edit a force field and check the changes
     * \param[in] options The options of the edit, an empty value for a flag
     * \param[in] input   The force field, the one of the other tests if empty
     */
    void checkChanges(const Options &options, const std::string &input = "")
    {
        const auto ff = input.empty() ? input_ : input;
        // The force field as edit_ff writes it, such that the numbers of
        // what did not change are written in the same way
        const auto before = read(run(ff, { { "-write", "" } }, "before.xml"));
        const auto after  = read(run(ff, options, "after.xml"));
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        checker.setDefaultTolerance(buildTolerance(1e-6, 1e-8, 1e-6));
        checkJson(&checker, forceFieldChanges(before, after), "changes");
    }
};

TEST_F(EditFFResultsTest, ValueAndBounds)
{
    checkChanges({ { "-p", "chi" }, { "-a", "h_z f_z" }, { "-min", "3" }, { "-max", "8" }, { "-val", "4.5" } });
}

TEST_F(EditFFResultsTest, Mutability)
{
    checkChanges({ { "-p", "chi" }, { "-a", "h_z f_z" }, { "-mut", "Fixed" } });
}

TEST_F(EditFFResultsTest, ScaleAll)
{
    // A parameter of all types
    checkChanges({ { "-p", "eta" }, { "-scale", "1.1" } });
}

TEST_F(EditFFResultsTest, Limits)
{
    checkChanges({ { "-p", "chi" }, { "-a", "h_z" }, { "-limits", "2" } });
}

TEST_F(EditFFResultsTest, ForcedCharge)
{
    checkChanges({ { "-p", "charge" }, { "-a", "Be2+" }, { "-val", "1.5" }, { "-force", "" } });
}

TEST_F(EditFFResultsTest, NonNegative)
{
    checkChanges({ { "-p", "eta" }, { "-nonneg", "no" } });
}

TEST_F(EditFFResultsTest, Remove)
{
    checkChanges({ { "-a", "f" }, { "-del", "" } });
}

TEST_F(EditFFResultsTest, CombinationRule)
{
    // A small force field, since all pairs of types change
    checkChanges({ { "-cr", "VANDERWAALS:sigma:Arithmetic" } },
                 fileManager().getInputFilePath("forcefields/TIP5P.xml"));
}

// The editing without a command line, as a program that uses ACT does it

TEST_F(EditFFTest, WithoutCommandLineSameAsTheProgram)
{
    // A change of bounds and value of two types, and a combination rule
    const auto tool = read(run(input_, { { "-p", "chi" }, { "-a", "h_z f_z" }, { "-min", "3" }, { "-max", "8" },
                                         { "-val", "4.5" }, { "-cr", "VANDERWAALS:sigma:Arithmetic" } },
                               "tool.xml"));
    EditFFOptions options;
    options.forceFieldFile   = input_;
    options.outputFile       = fileManager().getTemporaryFilePath("direct.xml");
    options.parameter        = "chi";
    options.particles        = "h_z f_z";
    options.minimum          = 3;
    options.maximum          = 8;
    options.value            = 4.5;
    options.combinationRules = "VANDERWAALS:sigma:Arithmetic";
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    runEditFF(&msghandler, options);
    const auto direct = read(options.outputFile);
    // As text: gtest cannot print a JSON value
    EXPECT_EQ("{}", forceFieldChanges(tool, direct).dump());
    // Something did change
    EXPECT_NE("{}", forceFieldChanges(read(input_), direct).dump());
}

TEST_F(EditFFTest, WithoutCommandLineWrongOptions)
{
    const auto valid = [this]()
    {
        EditFFOptions options;
        options.forceFieldFile = input_;
        options.outputFile     = fileManager().getTemporaryFilePath("wrong.xml");
        options.parameter      = "chi";
        return options;
    };
    std::vector<std::pair<EditFFOptions, std::string>> wrong;
    {
        auto o = valid();
        o.forceFieldFile.clear();
        wrong.emplace_back(o, "-ff");
    }
    {
        auto o = valid();
        o.mutability = "Changeable";
        wrong.emplace_back(o, "-mut");
    }
    {
        auto o = valid();
        o.scale = std::nan("");
        wrong.emplace_back(o, "-scale");
    }
    {
        auto o = valid();
        o.copyMissing = "EEM";
        wrong.emplace_back(o, "-ff2");
    }
    for (const auto &w : wrong)
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        try
        {
            runEditFF(&msghandler, w.first);
            ADD_FAILURE() << "Accepted options with a wrong " << w.second;
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find(w.second)) << ex.what();
        }
        EXPECT_FALSE(std::filesystem::exists(w.first.outputFile)) << w.second;
    }
}

TEST_F(EditFFTest, MutabilityOfSeveralIdentifiers)
{
    const auto output = fileManager().getTemporaryFilePath("edited.xml");
    gmx::test::CommandLine cmd;
    cmd.append("edit_ff");
    cmd.addOption("-ff", fileManager().getInputFilePath("forcefields/ACM-g1.xml"));
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

// Test: A combination rule from the command line replaces that of the
// force field, the others do not change
TEST_F(EditFFTest, CombinationRule)
{
    const auto pd = read(run(input_, { { "-cr", "VANDERWAALS:sigma:Arithmetic" } }));
    const auto &vdw = pd.findForcesConst(InteractionType::VDW);
    EXPECT_EQ(CombRule::Arithmetic, vdw.combinationRuleConst("sigma").rule());
    EXPECT_EQ(CombRule::Geometric, vdw.combinationRuleConst("epsilon").rule());
}

//! \return the charge parameter of a particle type
static const ForceFieldParameter &particleCharge(const ForceField &pd, const std::string &id)
{
    return pd.findParticleType(id)->parameterConst("charge");
}

// Test: The fixed charge of a particle type is only changed with -force,
// the other particle types keep theirs
TEST_F(EditFFTest, ParticleChargeForced)
{
    const auto ref = read(input_);
    const auto pd  = read(run(input_, { { "-p", "charge" }, { "-a", "Be2+" }, { "-val", "1.5" }, { "-force", "" } }));
    EXPECT_DOUBLE_EQ(1.5, particleCharge(pd, "Be2+").value());
    EXPECT_DOUBLE_EQ(particleCharge(ref, "Mg2+").value(), particleCharge(pd, "Mg2+").value());
}

// Test: The bounds of a fixed charge can be changed without -force, its
// value does not change
TEST_F(EditFFTest, ParticleChargeBounds)
{
    const auto ref = read(input_);
    const auto pd  = read(run(input_, { { "-p", "charge" }, { "-a", "Be2+" }, { "-min", "1" }, { "-max", "3" } }));
    expectRange(particleCharge(pd, "Be2+"), 1, particleCharge(ref, "Be2+").value(), 3);
}

// Test: A new value of one parameter of one identifier, the others do not change
TEST_F(EditFFTest, Value)
{
    const auto pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-val", "3.9" } }));
    expectRange(acm(pd, "h_z", "chi"), 2.8, 3.9, 4.375);
    EXPECT_DOUBLE_EQ(14, acm(pd, "h_z", "eta").value());
    const auto ref = read(input_);
    EXPECT_DOUBLE_EQ(acm(ref, "f_z", "chi").value(), acm(pd, "f_z", "chi").value());
}

// Test: A value outside the bounds widens them
TEST_F(EditFFTest, ValueOutsideBounds)
{
    auto pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-val", "5" } }));
    expectRange(acm(pd, "h_z", "chi"), 2.8, 5, 5);
    pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-val", "2" } }));
    expectRange(acm(pd, "h_z", "chi"), 2, 2, 4.375);
}

// Test: Scaling a value, of several identifiers at once
TEST_F(EditFFTest, Scale)
{
    const auto pd  = read(run(input_, { { "-p", "eta" }, { "-a", "h_z f_z" }, { "-scale", "1.1" } }));
    const auto ref = read(input_);
    expectRealEq(1.1*14, acm(pd, "h_z", "eta").value());
    expectRealEq(1.1*acm(ref, "f_z", "eta").value(), acm(pd, "f_z", "eta").value());
}

// Test: A new minimum moves the value into the bounds, also beyond the maximum
TEST_F(EditFFTest, Minimum)
{
    auto pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-min", "3.8" } }));
    expectRange(acm(pd, "h_z", "chi"), 3.8, 3.8, 4.375);
    pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-min", "5" } }));
    expectRange(acm(pd, "h_z", "chi"), 5, 5, 5);
}

// Test: A new maximum moves the value into the bounds, also below the minimum
TEST_F(EditFFTest, Maximum)
{
    auto pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-max", "3" } }));
    expectRange(acm(pd, "h_z", "chi"), 2.8, 3, 3);
    pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-max", "2" } }));
    expectRange(acm(pd, "h_z", "chi"), 2, 2, 2);
}

// Test: Limits around the value, given as a factor below or above one
TEST_F(EditFFTest, Limits)
{
    for (const char *factor : { "0.5", "2" })
    {
        SCOPED_TRACE(factor);
        const auto pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-limits", factor } }));
        expectRange(acm(pd, "h_z", "chi"), 1.75, 3.5, 7);
    }
}

// Test: Stretching the bounds of a parameter at its minimum or maximum
TEST_F(EditFFTest, Stretch)
{
    const auto atMinimum = run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-min", "3.5" } }, "minimum.xml");
    auto pd        = read(run(atMinimum, { { "-p", "chi" }, { "-a", "h_z" }, { "-stretch", "" } }));
    // The range of 3.5 to 4.375 is extended by a fifth below
    expectRange(acm(pd, "h_z", "chi"), 3.5 - 0.2*0.875, 3.5, 4.375);
    const auto atMaximum = run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-val", "7" } }, "maximum.xml");
    pd             = read(run(atMaximum, { { "-p", "chi" }, { "-a", "h_z" }, { "-stretch", "" } }, "stretched.xml"));
    // The range of 2.8 to 7 is extended by a fifth above
    expectRange(acm(pd, "h_z", "chi"), 2.8, 7, 7 + 0.2*4.2);
}

// Test: A fixed parameter only changes with -force, without changes no file
// is written
TEST_F(EditFFTest, FixedParameter)
{
    // The width of the charge distribution is a parameter of the particle type
    const auto unchanged = run(input_, { { "-p", "zeta" }, { "-a", "h" }, { "-val", "25" } });
    EXPECT_FALSE(std::filesystem::exists(unchanged));
    auto pd   = read(run(input_, { { "-p", "zeta" }, { "-a", "h" }, { "-val", "25" }, { "-force", "" } }));
    const auto zeta = pd.findParticleType("h")->parameterConst(zetaName);
    EXPECT_DOUBLE_EQ(25, zeta.value());
    EXPECT_EQ(Mutability::Fixed, zeta.mutability());
    // A file is written with -write even without changes
    EXPECT_TRUE(std::filesystem::exists(run(input_, { { "-write", "" } }, "copy.xml")));
}

// Test: Removing a particle type removes it and the interactions named
// after it, not those of its subtypes, which other types may share
TEST_F(EditFFTest, Remove)
{
    const auto pd  = read(run(input_, { { "-a", "f" }, { "-del", "" } }));
    const auto ref = read(input_);
    EXPECT_TRUE(ref.hasParticleType("f"));
    EXPECT_FALSE(pd.hasParticleType("f"));
    EXPECT_TRUE(pd.hasParticleType("f_s"));
    EXPECT_TRUE(ref.findForcesConst(InteractionType::VDW).parameterExists(Identifier("f")));
    EXPECT_FALSE(pd.findForcesConst(InteractionType::VDW).parameterExists(Identifier("f")));
    EXPECT_TRUE(pd.findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION).parameterExists(Identifier("f_z")));
}

// Test: The parameters of some identifiers are written to a file
TEST_F(EditFFTest, Dump)
{
    auto dump = fileManager().getTemporaryFilePath("params.dat");
    run(input_, { { "-ana", "EEM" }, { "-a", "h_z h" }, { "-dump", dump } }, "");
    std::ifstream            in(dump);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line); )
    {
        lines.push_back(line);
    }
    std::sort(lines.begin(), lines.end());
    // The EEM group includes the widths of the charges, which are those of
    // the particle types
    const std::vector<std::string> expected = { "chi|h_z|3.5", "eta|h_z|14", "zeta|h|30" };
    EXPECT_EQ(expected, lines);
}

/*! \brief A second force field, with other values and without some parameters
 * Chi of h_z is 4 and that of f_z 6, within their bounds, the EEM
 * parameters of i_z are missing unless withIodine, then its chi is 5.
 * \param[in] input      The force field to start from
 * \param[in] output     Where to write the result
 * \param[in] withIodine Whether to keep the EEM parameters of i_z
 */
static void writeSecondForceField(const std::string &input, const std::string &output,
                                  bool withIodine = false)
{
    ForceField pd;
    readForceField(input, &pd);
    auto eem = pd.findForces(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
    eem->findParameterType(Identifier("h_z"), "chi")->setValue(4);
    eem->findParameterType(Identifier("f_z"), "chi")->setValue(6);
    if (withIodine)
    {
        eem->findParameterType(Identifier("i_z"), "chi")->setValue(5);
    }
    else
    {
        eem->parameters()->erase(Identifier("i_z"));
    }
    writeForceField(output, &pd, false);
}

/*! \brief Write a force field with problems for -validate
 * In the force field of the tests: chi of h_z above its maximum, eta of
 * h_z twice, the charge of particle type h twice, the mass of h above a
 * new maximum, chi of i_z negative while it is marked as non-negative, and
 * eta of i_z, which is fixed, above a new maximum.
 * \param[in] input  The force field to start from
 * \param[in] output Where to write the result
 */
static void writeBrokenForceField(const std::string &input, const std::string &output)
{
    std::ifstream in(input);
    std::string   text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const auto replace = [&text](const std::string &from, const std::string &to)
          {
              auto pos = text.find(from);
              ASSERT_NE(std::string::npos, pos) << from;
              ASSERT_EQ(std::string::npos, text.find(from, pos + 1)) << "not unique: " << from;
              text.replace(pos, from.size(), to);
          };
    const std::string chi = "      <parameter type=\"chi\" unit=\"eV\" value=\"3.5\" uncertainty=\"0\" minimum=\"2.8\" maximum=\"4.375\" ntrain=\"1\" mutability=\"Bounded\" nonnegative=\"yes\"/>\n";
    const std::string eta = "      <parameter type=\"eta\" unit=\"eV/e\" value=\"14\" uncertainty=\"0\" minimum=\"11.2\" maximum=\"17.5\" ntrain=\"1\" mutability=\"Bounded\" nonnegative=\"yes\"/>\n";
    std::string eta13 = eta;
    eta13.replace(eta13.find("value=\"14\""), 10, "value=\"13\"");
    std::string chi5 = chi;
    chi5.replace(chi5.find("value=\"3.5\""), 11, "value=\"5\"");
    replace("<parameterlist identifier=\"h_z\">\n" + chi + eta,
            "<parameterlist identifier=\"h_z\">\n" + chi5 + eta + eta13);
    const std::string charge = "      <parameter type=\"charge\" unit=\"e\" value=\"0\" uncertainty=\"0\" minimum=\"-1\" maximum=\"1\" ntrain=\"0\" mutability=\"ACM\" nonnegative=\"no\"/>\n";
    std::string charge05 = charge;
    charge05.replace(charge05.find("value=\"0\""), 9, "value=\"0.5\"");
    replace("<option key=\"vdwtype\" value=\"h\"/>\n" + charge,
            "<option key=\"vdwtype\" value=\"h\"/>\n" + charge + charge05);
    replace("value=\"1.0079\" uncertainty=\"0\" minimum=\"0.80632\" maximum=\"1.20948\"",
            "value=\"1.0079\" uncertainty=\"0\" minimum=\"0.80632\" maximum=\"1\"");
    replace("<parameterlist identifier=\"i_z\">\n"
            "      <parameter type=\"chi\" unit=\"eV\" value=\"4\" uncertainty=\"0\" minimum=\"3.2\" maximum=\"5\" ntrain=\"1\" mutability=\"Bounded\" nonnegative=\"yes\"/>\n"
            "      <parameter type=\"eta\" unit=\"eV/e\" value=\"10\" uncertainty=\"0\" minimum=\"8\" maximum=\"12.5\" ntrain=\"1\" mutability=\"Bounded\" nonnegative=\"yes\"/>\n",
            "<parameterlist identifier=\"i_z\">\n"
            "      <parameter type=\"chi\" unit=\"eV\" value=\"-1\" uncertainty=\"0\" minimum=\"-2\" maximum=\"5\" ntrain=\"1\" mutability=\"Bounded\" nonnegative=\"yes\"/>\n"
            "      <parameter type=\"eta\" unit=\"eV/e\" value=\"10\" uncertainty=\"0\" minimum=\"8\" maximum=\"9\" ntrain=\"1\" mutability=\"Fixed\" nonnegative=\"yes\"/>\n");
    std::ofstream out(output);
    out << text;
}

//! \return the number of problems of a kind
static size_t countIssues(const std::vector<ForceFieldIssue> &issues, ForceFieldIssue::Kind kind)
{
    return std::count_if(issues.begin(), issues.end(),
                         [kind](const ForceFieldIssue &i) { return i.kind == kind; });
}

// Test: The force field of the tests has no problems
TEST_F(EditFFTest, ValidateClean)
{
    auto issues = validateForceField(read(input_));
    EXPECT_TRUE(issues.empty()) << issues.size() << " problems, the first: "
                                << (issues.empty() ? "" : issues[0].identifier + " " + issues[0].text);
    EXPECT_TRUE(repairedIssuesSummary(issues).empty());
}

// Test: Each kind of problem is found, with what reading did about it
TEST_F(EditFFTest, ValidateFindsProblems)
{
    const auto broken = fileManager().getTemporaryFilePath("broken.xml");
    writeBrokenForceField(input_, broken);
    const auto pd     = read(broken);
    const auto issues = validateForceField(pd);
    // The values of chi of h_z and of the mass of h were set to their maximum
    EXPECT_EQ(2u, countIssues(issues, ForceFieldIssue::Kind::ResetToBounds));
    EXPECT_EQ(2u, countIssues(issues, ForceFieldIssue::Kind::Duplicate));
    EXPECT_EQ(1u, countIssues(issues, ForceFieldIssue::Kind::OutOfBounds));
    EXPECT_EQ(1u, countIssues(issues, ForceFieldIssue::Kind::Negative));
    EXPECT_EQ(0u, countIssues(issues, ForceFieldIssue::Kind::MinimumAboveMaximum));
    EXPECT_EQ(6u, issues.size());
    for (const auto &issue : issues)
    {
        SCOPED_TRACE(issue.identifier + " " + issue.parameter + ": " + issue.text);
        bool byReading = (issue.kind == ForceFieldIssue::Kind::ResetToBounds ||
                          issue.kind == ForceFieldIssue::Kind::Duplicate);
        EXPECT_EQ(byReading, issue.repaired);
    }
    // Reading set chi to its maximum and used the first eta and charge
    expectRealEq(4.375, acm(pd, "h_z", "chi").value());
    EXPECT_TRUE(acm(pd, "h_z", "chi").wasReset());
    expectRealEq(5, acm(pd, "h_z", "chi").valueBeforeReset());
    expectRealEq(14, acm(pd, "h_z", "eta").value());
    expectRealEq(0, particleCharge(pd, "h").value());
    const auto summary = repairedIssuesSummary(issues);
    EXPECT_NE(std::string::npos, summary.find("2 parameter(s) were outside their bounds")) << summary;
    EXPECT_NE(std::string::npos, summary.find("2 were defined more than once")) << summary;
    EXPECT_NE(std::string::npos, summary.find("edit_ff -validate")) << summary;
}

// Test: edit_ff -validate writes a force field without the problems that
// reading repairs, the others are left
TEST_F(EditFFTest, ValidateWrites)
{
    const auto broken = fileManager().getTemporaryFilePath("broken.xml");
    writeBrokenForceField(input_, broken);
    auto pd     = read(run(broken, { { "-validate", "" } }, "validated.xml"));
    const auto issues = validateForceField(pd);
    EXPECT_EQ(2u, issues.size());
    EXPECT_EQ(1u, countIssues(issues, ForceFieldIssue::Kind::OutOfBounds));
    EXPECT_EQ(1u, countIssues(issues, ForceFieldIssue::Kind::Negative));
    expectRealEq(4.375, acm(pd, "h_z", "chi").value());
    EXPECT_FALSE(acm(pd, "h_z", "chi").wasReset());
    expectRealEq(14, acm(pd, "h_z", "eta").value());
    expectRealEq(1, pd.findParticleType("h")->parameterConst("mass").value());
    // The rest of the force field is as before
    const auto ref = read(input_);
    expectRealEq(acm(ref, "i_z", "eta").value(), acm(pd, "i_z", "eta").value());
    expectRealEq(acm(ref, "h_z", "eta").value(), acm(pd, "h_z", "eta").value());
}

//! \return the content of a JSON file, which has to exist
static nlohmann::json readJson(const std::string &fileName)
{
    std::ifstream in(fileName);
    EXPECT_TRUE(in.good()) << fileName;
    return in.good() ? nlohmann::json::parse(in) : nlohmann::json();
}

/*! \brief Compare a parameter in the JSON output with that of a force field
 * \param[in] j     The parameter in JSON
 * \param[in] type  The type of the parameter
 * \param[in] param The parameter
 */
static void expectParameter(const nlohmann::json &j, const std::string &type,
                            const ForceFieldParameter &param)
{
    SCOPED_TRACE("parameter " + type);
    EXPECT_EQ(type, j.at("type").get<std::string>());
    EXPECT_EQ(param.unit(), j.at("unit").get<std::string>());
    // The force field file has thirteen digits, the JSON output all
    const auto near = [](double expected, double actual)
          {
              EXPECT_NEAR(expected, actual, 1e-10*std::max(1.0, std::abs(expected)));
          };
    near(param.value(), j.at("value").get<double>());
    near(param.minimum(), j.at("minimum").get<double>());
    near(param.maximum(), j.at("maximum").get<double>());
    near(param.uncertainty(), j.at("uncertainty").get<double>());
    EXPECT_EQ(param.ntrain(), j.at("ntrain").get<unsigned int>());
    EXPECT_EQ(mutabilityName(param.mutability()), j.at("mutability").get<std::string>());
    EXPECT_EQ(param.nonNegative(), j.at("nonnegative").get<bool>());
}

/*! \brief Compare the JSON description of a force field with the force field
 * \param[in] ff The description
 * \param[in] pd The force field
 */
static void expectForceField(const nlohmann::json &ff, const ForceField &pd)
{
    EXPECT_EQ(pd.polarizable(), ff.at("polarizable").get<bool>());
    // The particle types, in the order of the force field
    const auto &ptypes = ff.at("particletypes").at("types");
    ASSERT_EQ(pd.particleTypesConst().size(), ptypes.size());
    size_t ipt = 0;
    for (const auto &pt : pd.particleTypesConst())
    {
        const auto &jpt = ptypes[ipt++];
        SCOPED_TRACE("particle type " + pt.second.id().id());
        EXPECT_EQ(pt.second.id().id(), jpt.at("identifier").get<std::string>());
        EXPECT_EQ(actParticleToString(pt.second.apType()), jpt.at("type").get<std::string>());
        EXPECT_EQ(pt.second.description(), jpt.at("description").get<std::string>());
        for (const auto &opt : pt.second.optionsConst())
        {
            EXPECT_EQ(opt.second, jpt.at("options").at(opt.first).get<std::string>());
        }
        ASSERT_EQ(pt.second.parametersConst().size(), jpt.at("parameters").size());
        size_t ip = 0;
        for (const auto &param : pt.second.parametersConst())
        {
            expectParameter(jpt.at("parameters")[ip++], param.first, param.second);
        }
    }
    // The interactions
    const auto &inters = ff.at("interactions");
    ASSERT_EQ(pd.forcesConst().size(), inters.size());
    size_t ifs     = 0;
    size_t nparams = 0;
    for (const auto &fs : pd.forcesConst())
    {
        const auto &jfs = inters[ifs++];
        SCOPED_TRACE("interaction " + interactionTypeToString(fs.first));
        EXPECT_EQ(interactionTypeToString(fs.first), jfs.at("type").get<std::string>());
        EXPECT_EQ(potentialToString(fs.second.potential()), jfs.at("function").get<std::string>());
        for (const auto &opt : fs.second.option())
        {
            EXPECT_EQ(opt.second, jfs.at("options").at(opt.first).get<std::string>());
        }
        // Objects, also when there is nothing in them
        EXPECT_TRUE(jfs.at("options").is_object());
        EXPECT_TRUE(jfs.at("combination_rules").is_object());
        EXPECT_EQ(fs.second.combinationRules().size(), jfs.at("combination_rules").size());
        for (const auto &cr : fs.second.combinationRules())
        {
            const auto &jcr = jfs.at("combination_rules").at(cr.first);
            EXPECT_EQ(combinationRuleName(cr.second.rule()), jcr.at("rule").get<std::string>());
            expectParameter(jcr.at("exponent"), "exponent", cr.second.ffplConst());
        }
        // As in the file, without the lists that are derived from others
        size_t ilist = 0;
        for (const auto &params : fs.second.parametersConst())
        {
            const bool dependent = std::all_of(params.second.begin(), params.second.end(),
                                               [](const auto &p) { return Mutability::Dependent == p.second.mutability(); });
            if (dependent)
            {
                continue;
            }
            ASSERT_LT(ilist, jfs.at("parameterlists").size());
            const auto &jlist = jfs.at("parameterlists")[ilist++];
            SCOPED_TRACE("identifier " + params.first.id());
            EXPECT_EQ(params.first.id(), jlist.at("identifier").get<std::string>());
            ASSERT_EQ(params.second.size(), jlist.at("parameters").size());
            size_t ip = 0;
            for (const auto &param : params.second)
            {
                expectParameter(jlist.at("parameters")[ip++], param.first, param.second);
                nparams++;
            }
        }
        EXPECT_EQ(ilist, jfs.at("parameterlists").size());
    }
    // There was something to compare
    EXPECT_GT(nparams, 50u);
}

// Test: -json describes the whole force field: the particle types and the
// interactions with their options, combination rules and parameters
TEST_F(EditFFTest, Json)
{
    for (const auto *ff : { "forcefields/ACM-g1.xml", "forcefields/ACS-pg1.xml" })
    {
        SCOPED_TRACE(ff);
        const auto input    = fileManager().getInputFilePath(ff);
        const auto jsonFile = fileManager().getTemporaryFilePath("ff.json");
        std::filesystem::remove(jsonFile);
        run(input, { { "-json", jsonFile } }, "");
        const auto j = readJson(jsonFile);
        ASSERT_TRUE(j.contains("edit_ff"));
        // Nothing was asked for, so nothing changed and nothing was written
        EXPECT_FALSE(j.at("edit_ff").at("changed").get<bool>());
        EXPECT_FALSE(j.at("edit_ff").at("written").get<bool>());
        const auto pd = read(input);
        expectForceField(j.at("edit_ff").at("forcefield"), pd);
        EXPECT_EQ(pd.checkSum(), j.at("edit_ff").at("forcefield").at("checksum").get<std::string>());
        // No problems in the force fields of the tests
        const auto &validation = j.at("edit_ff").at("validation");
        EXPECT_TRUE(validation.at("ok").get<bool>());
        EXPECT_TRUE(validation.at("issues").is_array());
        EXPECT_TRUE(validation.at("issues").empty());
    }
}

// Test: The JSON output has the force field after the changes, which is
// the one that is written with -o
TEST_F(EditFFTest, JsonAfterChange)
{
    const auto jsonFile = fileManager().getTemporaryFilePath("changed.json");
    const auto out      = run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-val", "3.6" },
                                        { "-json", jsonFile } }, "changed.xml");
    const auto pd = read(out);
    expectRealEq(3.6, acm(pd, "h_z", "chi").value());
    const auto root = readJson(jsonFile).at("edit_ff");
    EXPECT_TRUE(root.at("changed").get<bool>());
    EXPECT_TRUE(root.at("written").get<bool>());
    const auto &ff = root.at("forcefield");
    expectForceField(ff, pd);
    bool found = false;
    for (const auto &jfs : ff.at("interactions"))
    {
        if (jfs.at("type").get<std::string>() !=
            interactionTypeToString(InteractionType::ELECTRONEGATIVITYEQUALIZATION))
        {
            continue;
        }
        for (const auto &jlist : jfs.at("parameterlists"))
        {
            if (jlist.at("identifier").get<std::string>() != "h_z")
            {
                continue;
            }
            for (const auto &jp : jlist.at("parameters"))
            {
                if (jp.at("type").get<std::string>() == "chi")
                {
                    expectRealEq(3.6, jp.at("value").get<double>());
                    found = true;
                }
            }
        }
    }
    EXPECT_TRUE(found);
}

// Test: A change of a type that is not in the force field changes nothing:
// no file is written, and the JSON output says so
TEST_F(EditFFTest, JsonNothingChanged)
{
    const auto jsonFile = fileManager().getTemporaryFilePath("same.json");
    const auto output   = fileManager().getTemporaryFilePath("same.xml");
    std::filesystem::remove(output);
    gmx::test::CommandLine cmd;
    cmd.append("edit_ff");
    cmd.addOption("-ff", input_);
    cmd.addOption("-o", output);
    cmd.addOption("-p", "chi");
    cmd.addOption("-a", "nosuchtype_z");
    cmd.addOption("-val", "3.6");
    cmd.addOption("-g", fileManager().getTemporaryFilePath("edit_ff.log"));
    cmd.addOption("-json", jsonFile);
    ASSERT_EQ(0, edit_ff(cmd.argc(), cmd.argv()));
    EXPECT_FALSE(std::filesystem::exists(output));
    const auto root = readJson(jsonFile).at("edit_ff");
    EXPECT_FALSE(root.at("changed").get<bool>());
    EXPECT_FALSE(root.at("written").get<bool>());
    expectForceField(root.at("forcefield"), read(input_));
}

// Test: The problems that -validate prints are in the JSON output, with
// their kind and whether writing the force field removes them
TEST_F(EditFFTest, JsonValidation)
{
    const auto broken = fileManager().getTemporaryFilePath("broken.xml");
    writeBrokenForceField(input_, broken);
    const auto jsonFile = fileManager().getTemporaryFilePath("broken.json");
    run(broken, { { "-validate", "" }, { "-json", jsonFile } }, "");
    const auto validation = readJson(jsonFile).at("edit_ff").at("validation");
    EXPECT_FALSE(validation.at("ok").get<bool>());
    const auto &issues = validation.at("issues");
    ASSERT_EQ(6u, issues.size());
    std::map<std::string, int> kinds;
    int                        nrepaired = 0;
    for (const auto &issue : issues)
    {
        const auto kind = issue.at("kind").get<std::string>();
        kinds[kind] += 1;
        EXPECT_FALSE(issue.at("where").get<std::string>().empty());
        EXPECT_FALSE(issue.at("identifier").get<std::string>().empty());
        EXPECT_FALSE(issue.at("parameter").get<std::string>().empty());
        EXPECT_FALSE(issue.at("text").get<std::string>().empty());
        const bool repaired = issue.at("repaired").get<bool>();
        EXPECT_EQ(kind == "ResetToBounds" || kind == "Duplicate", repaired);
        nrepaired += repaired ? 1 : 0;
    }
    EXPECT_EQ(2, kinds["ResetToBounds"]);
    EXPECT_EQ(2, kinds["Duplicate"]);
    EXPECT_EQ(1, kinds["OutOfBounds"]);
    EXPECT_EQ(1, kinds["Negative"]);
    EXPECT_EQ(4, nrepaired);
    EXPECT_EQ(4, validation.at("repaired").get<int>());
    // chi of h_z is one of them
    const bool chi = std::any_of(issues.begin(), issues.end(), [](const nlohmann::json &i)
                                 { return i.at("identifier") == "h_z" && i.at("parameter") == "chi"; });
    EXPECT_TRUE(chi);
}

// Test: Without the option no JSON file is written
TEST_F(EditFFTest, NoJsonByDefault)
{
    std::remove("edit_ff.json");
    run(input_, {}, "");
    EXPECT_FALSE(std::ifstream("edit_ff.json").good());
}

// Test: -nonneg no lets parameters of an interaction be negative, only
// those that are selected, which removes the problem that -validate reports
TEST_F(EditFFTest, NonNegativeCleared)
{
    const auto broken = fileManager().getTemporaryFilePath("broken.xml");
    writeBrokenForceField(input_, broken);
    const auto pd = read(run(broken, { { "-p", "chi" }, { "-a", "i_z" }, { "-nonneg", "no" } }));
    EXPECT_FALSE(acm(pd, "i_z", "chi").nonNegative());
    EXPECT_TRUE(acm(pd, "h_z", "chi").nonNegative());
    EXPECT_TRUE(acm(pd, "i_z", "eta").nonNegative());
    expectRealEq(-1, acm(pd, "i_z", "chi").value());
    EXPECT_EQ(0u, countIssues(validateForceField(pd), ForceFieldIssue::Kind::Negative));
}

// Test: -nonneg yes makes all parameters of a type non-negative
TEST_F(EditFFTest, NonNegativeSet)
{
    const auto cleared = read(run(input_, { { "-p", "eta" }, { "-nonneg", "no" } }, "cleared.xml"));
    EXPECT_FALSE(acm(cleared, "h_z", "eta").nonNegative());
    EXPECT_FALSE(acm(cleared, "i_z", "eta").nonNegative());
    const auto set = read(run(fileManager().getTemporaryFilePath("cleared.xml"),
                              { { "-p", "eta" }, { "-nonneg", "yes" } }, "set.xml"));
    EXPECT_TRUE(acm(set, "h_z", "eta").nonNegative());
    EXPECT_TRUE(acm(set, "i_z", "eta").nonNegative());
}

// Test: Without problems that reading repairs, -validate writes no file
TEST_F(EditFFTest, ValidateNothingToWrite)
{
    const auto out = run(input_, { { "-validate", "" } }, "unchanged.xml");
    EXPECT_FALSE(std::filesystem::exists(out));
}

// Test: Parameters of a second force field that are missing in the first
// are copied, the others are not changed
TEST_F(EditFFTest, CopyMissing)
{
    const auto second = fileManager().getTemporaryFilePath("second.xml");
    writeSecondForceField(input_, second);
    const auto pd = read(run(second, { { "-ff2", input_ }, { "-copy_missing", "EEM" } }));
    const auto ref = read(input_);
    EXPECT_DOUBLE_EQ(acm(ref, "i_z", "chi").value(), acm(pd, "i_z", "chi").value());
    EXPECT_DOUBLE_EQ(4, acm(pd, "h_z", "chi").value());
}

// Test: All parameters of a group are replaced by those of a second force field
TEST_F(EditFFTest, Replace)
{
    auto second = fileManager().getTemporaryFilePath("second.xml");
    writeSecondForceField(input_, second, true);
    const auto pd = read(run(input_, { { "-ff2", second }, { "-replace", "EEM" } }));
    EXPECT_DOUBLE_EQ(4, acm(pd, "h_z", "chi").value());
    EXPECT_DOUBLE_EQ(6, acm(pd, "f_z", "chi").value());
    EXPECT_DOUBLE_EQ(5, acm(pd, "i_z", "chi").value());
}

// Test: Only the given parameter of the given identifier is implanted from a
// second force field
TEST_F(EditFFTest, Implant)
{
    auto second = fileManager().getTemporaryFilePath("second.xml");
    writeSecondForceField(input_, second);
    const auto pd  = read(run(input_, { { "-ff2", second }, { "-implant", "EEM" }, { "-a", "h_z" }, { "-p", "chi" } }));
    const auto ref = read(input_);
    EXPECT_DOUBLE_EQ(4, acm(pd, "h_z", "chi").value());
    EXPECT_DOUBLE_EQ(14, acm(pd, "h_z", "eta").value());
    EXPECT_DOUBLE_EQ(acm(ref, "f_z", "chi").value(), acm(pd, "f_z", "chi").value());
}

// Test: Analyzing and comparing force fields does not change them
TEST_F(EditFFTest, AnalyzeAndCompare)
{
    auto second = fileManager().getTemporaryFilePath("second.xml");
    writeSecondForceField(input_, second);
    EXPECT_FALSE(std::filesystem::exists(run(input_, { { "-ff2", second }, { "-p", "chi" } })));
    EXPECT_FALSE(std::filesystem::exists(run(input_, { { "-bounds", "" } }, "bounds.xml")));
}

/*! \brief Plot the interactions of a force field
 * The plots are written in the working directory, hence in a new one.
 * \param[in]  forceField The force field file
 * \param[in]  dir        Name of the directory for the plots
 * \param[out] status     The exit status of edit_ff
 * \return the directory with the plots
 */
static std::filesystem::path plotInteractions(const std::string &forceField,
                                              const std::string &dir,
                                              int               *status)
{
    std::filesystem::path plots(dir);
    std::filesystem::create_directories(plots);
    const auto cwd = std::filesystem::current_path();
    std::filesystem::current_path(plots);
    gmx::test::CommandLine cmd;
    cmd.append("edit_ff");
    cmd.addOption("-ff", forceField);
    cmd.addOption("-ana", "BONDS");
    cmd.append("-plot");
    *status = edit_ff(cmd.argc(), cmd.argv());
    std::filesystem::current_path(cwd);
    return plots;
}

//! \return the point with the lowest energy of a plot
static std::pair<double, double> minimum(const std::filesystem::path &xvg)
{
    std::ifstream in(xvg);
    std::pair<double, double> lowest = { 0, 0 };
    bool                      first  = true;
    for (double x, v; in >> x >> v; )
    {
        if (first || v < lowest.second)
        {
            lowest = { x, v };
            first  = false;
        }
    }
    EXPECT_FALSE(first) << xvg;
    return lowest;
}

// Test: The bond and van der Waals curves of SPC/E have their minimum
// at the bond length and at the minimum of the Lennard-Jones potential
TEST_F(EditFFTest, Plot)
{
    int  status;
    const auto plots = plotInteractions(fileManager().getInputFilePath("forcefields/SPCE.xml"),
                                        fileManager().getTemporaryFilePath("plots"), &status);
    ASSERT_EQ(0, status);
    const auto bond = minimum(plots / "BONDS_h_b~o3_b.xvg");
    EXPECT_NEAR(0.1, bond.first, 1e-3);
    EXPECT_NEAR(0, bond.second, 1e-6);
    const auto vdw = minimum(plots / "VANDERWAALS_o3~o3.xvg");
    double sigma = 0.3165558, epsilon = 0.65017032;
    EXPECT_NEAR(std::pow(2.0, 1.0/6.0)*sigma, vdw.first, 1e-3);
    EXPECT_NEAR(-epsilon, vdw.second, 1e-3);
    EXPECT_TRUE(std::filesystem::exists(plots / "COULOMB_ow~ow.xvg"));
    std::filesystem::remove_all(plots);
}

// Test: Impropers of which the bonds have no parameters are skipped, and
// all other interactions are plotted
TEST_F(EditFFTest, PlotSkipsImpropersWithoutBonds)
{
    int  status;
    const auto plots = plotInteractions(fileManager().getInputFilePath("forcefields/ACS-g1.xml"),
                                        fileManager().getTemporaryFilePath("plots"), &status);
    ASSERT_EQ(0, status);
    std::map<std::string, int> count;
    for (const auto &entry : std::filesystem::directory_iterator(plots))
    {
        const auto name = entry.path().filename().string();
        if (entry.path().extension() == ".xvg")
        {
            EXPECT_GT(std::filesystem::file_size(entry.path()), 0u) << name;
            count[name.substr(0, name.find('_'))] += 1;
        }
    }
    for (const char *itype : { "BONDS", "ANGLES", "LINEAR", "PROPER", "IMPROPER", "VANDERWAALS", "COULOMB" })
    {
        EXPECT_GT(count[itype], 0) << itype;
    }
    std::filesystem::remove_all(plots);
}

} // namespace

} // namespace alexandria
