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
#include <filesystem>
#include <fstream>
#include <map>
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

//! The options of edit_ff besides the input force field
using Options = std::vector<std::pair<std::string, std::string>>;

class EditFFTest : public gmx::test::CommandLineTestBase
{
protected:
    //! The force field the tests start from
    std::string input_;

    EditFFTest() : input_(fileManager().getInputFilePath("ACM-g.xml")) {}

    //! \return the mutability of chi for an ACM type
    static Mutability chiMutability(const ForceField &pd, const std::string &acmtype)
    {
        auto fs = pd.findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
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

// Test: A combination rule from the command line replaces that of the
// force field, the others do not change
TEST_F(EditFFTest, CombinationRule)
{
    auto pd = read(run(input_, { { "-cr", "VANDERWAALS:sigma:Arithmetic" } }));
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
    auto ref = read(input_);
    auto pd  = read(run(input_, { { "-p", "charge" }, { "-a", "Be2+" }, { "-val", "1.5" }, { "-force", "" } }));
    EXPECT_DOUBLE_EQ(1.5, particleCharge(pd, "Be2+").value());
    EXPECT_DOUBLE_EQ(particleCharge(ref, "Mg2+").value(), particleCharge(pd, "Mg2+").value());
}

// Test: The bounds of a fixed charge can be changed without -force, its
// value does not change
TEST_F(EditFFTest, ParticleChargeBounds)
{
    auto ref = read(input_);
    auto pd  = read(run(input_, { { "-p", "charge" }, { "-a", "Be2+" }, { "-min", "1" }, { "-max", "3" } }));
    expectRange(particleCharge(pd, "Be2+"), 1, particleCharge(ref, "Be2+").value(), 3);
}

// Test: A new value of one parameter of one identifier, the others do not change
TEST_F(EditFFTest, Value)
{
    auto pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-val", "3.9" } }));
    expectRange(acm(pd, "h_z", "chi"), 2.8, 3.9, 4.375);
    EXPECT_DOUBLE_EQ(14, acm(pd, "h_z", "eta").value());
    auto ref = read(input_);
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
    auto pd  = read(run(input_, { { "-p", "eta" }, { "-a", "h_z f_z" }, { "-scale", "1.1" } }));
    auto ref = read(input_);
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
        auto pd = read(run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-limits", factor } }));
        expectRange(acm(pd, "h_z", "chi"), 1.75, 3.5, 7);
    }
}

// Test: Stretching the bounds of a parameter at its minimum or maximum
TEST_F(EditFFTest, Stretch)
{
    auto atMinimum = run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-min", "3.5" } }, "minimum.xml");
    auto pd        = read(run(atMinimum, { { "-p", "chi" }, { "-a", "h_z" }, { "-stretch", "" } }));
    // The range of 3.5 to 4.375 is extended by a fifth below
    expectRange(acm(pd, "h_z", "chi"), 3.5 - 0.2*0.875, 3.5, 4.375);
    auto atMaximum = run(input_, { { "-p", "chi" }, { "-a", "h_z" }, { "-val", "7" } }, "maximum.xml");
    pd             = read(run(atMaximum, { { "-p", "chi" }, { "-a", "h_z" }, { "-stretch", "" } }, "stretched.xml"));
    // The range of 2.8 to 7 is extended by a fifth above
    expectRange(acm(pd, "h_z", "chi"), 2.8, 7, 7 + 0.2*4.2);
}

// Test: A fixed parameter only changes with -force, without changes no file
// is written
TEST_F(EditFFTest, FixedParameter)
{
    auto unchanged = run(input_, { { "-p", "zeta" }, { "-a", "h_z" }, { "-val", "25" } });
    EXPECT_FALSE(std::filesystem::exists(unchanged));
    auto pd   = read(run(input_, { { "-p", "zeta" }, { "-a", "h_z" }, { "-val", "25" }, { "-force", "" } }));
    auto zeta = pd.findForcesConst(InteractionType::ELECTROSTATICS).findParameterTypeConst(Identifier("h_z"), "zeta");
    EXPECT_DOUBLE_EQ(25, zeta.value());
    EXPECT_EQ(Mutability::Fixed, zeta.mutability());
    // A file is written with -write even without changes
    EXPECT_TRUE(std::filesystem::exists(run(input_, { { "-write", "" } }, "copy.xml")));
}

// Test: Removing a particle type removes it and the interactions named
// after it, not those of its subtypes, which other types may share
TEST_F(EditFFTest, Remove)
{
    auto pd  = read(run(input_, { { "-a", "f" }, { "-del", "" } }));
    auto ref = read(input_);
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
    run(input_, { { "-ana", "EEM" }, { "-a", "h_z" }, { "-dump", dump } }, "");
    std::ifstream            in(dump);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line); )
    {
        lines.push_back(line);
    }
    std::sort(lines.begin(), lines.end());
    // The EEM group includes the widths of the charges
    std::vector<std::string> expected = { "chi|h_z|3.5", "eta|h_z|14", "zeta|h_z|30" };
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

// Test: Parameters of a second force field that are missing in the first
// are copied, the others are not changed
TEST_F(EditFFTest, CopyMissing)
{
    auto second = fileManager().getTemporaryFilePath("second.xml");
    writeSecondForceField(input_, second);
    auto pd = read(run(second, { { "-ff2", input_ }, { "-copy_missing", "EEM" } }));
    auto ref = read(input_);
    EXPECT_DOUBLE_EQ(acm(ref, "i_z", "chi").value(), acm(pd, "i_z", "chi").value());
    EXPECT_DOUBLE_EQ(4, acm(pd, "h_z", "chi").value());
}

// Test: All parameters of a group are replaced by those of a second force field
TEST_F(EditFFTest, Replace)
{
    auto second = fileManager().getTemporaryFilePath("second.xml");
    writeSecondForceField(input_, second, true);
    auto pd = read(run(input_, { { "-ff2", second }, { "-replace", "EEM" } }));
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
    auto pd  = read(run(input_, { { "-ff2", second }, { "-implant", "EEM" }, { "-a", "h_z" }, { "-p", "chi" } }));
    auto ref = read(input_);
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
    auto cwd = std::filesystem::current_path();
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
    auto plots = plotInteractions(fileManager().getInputFilePath("../../properties/tests/SPCE.xml"),
                                  fileManager().getTemporaryFilePath("plots"), &status);
    ASSERT_EQ(0, status);
    auto bond = minimum(plots / "BONDS_h_b~o3_b.xvg");
    EXPECT_NEAR(0.1, bond.first, 1e-3);
    EXPECT_NEAR(0, bond.second, 1e-6);
    auto vdw = minimum(plots / "VANDERWAALS_o3~o3.xvg");
    double sigma = 0.3165558, epsilon = 0.65017032;
    EXPECT_NEAR(std::pow(2.0, 1.0/6.0)*sigma, vdw.first, 1e-3);
    EXPECT_NEAR(-epsilon, vdw.second, 1e-3);
    EXPECT_TRUE(std::filesystem::exists(plots / "COULOMB_o3_z~o3_z.xvg"));
    std::filesystem::remove_all(plots);
}

// Test: Impropers of which the bonds have no parameters are skipped, and
// all other interactions are plotted
TEST_F(EditFFTest, PlotSkipsImpropersWithoutBonds)
{
    int  status;
    auto plots = plotInteractions(fileManager().getInputFilePath("../../alexandria/tests/ACS-g.xml"),
                                  fileManager().getTemporaryFilePath("plots"), &status);
    ASSERT_EQ(0, status);
    std::map<std::string, int> count;
    for (const auto &entry : std::filesystem::directory_iterator(plots))
    {
        auto name = entry.path().filename().string();
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
