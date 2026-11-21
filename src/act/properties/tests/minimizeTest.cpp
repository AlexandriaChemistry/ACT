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
 * Tests of the alexandria minimize tool, in particular of the JSON output
 * that the Python interface and the graphical user interface read.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <map>
#include <nlohmann/json.hpp>
#include <set>
#include <string>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/testutils/output_checks.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! Change the working directory for the lifetime of the object
class ScopedDirectory
{
public:
    /*! \brief Change to a directory
     * \param[in] dir The directory to go to
     */
    explicit ScopedDirectory(const std::string &dir)
        : old_(std::filesystem::current_path())
    {
        std::filesystem::current_path(dir);
    }
    //! Go back to the original directory
    ~ScopedDirectory()
    {
        std::filesystem::current_path(old_);
    }
private:
    //! The original directory
    std::filesystem::path old_;
};

class MinimizeTest : public gmx::test::CommandLineTestBase
{
protected:
    //! \return path to an input file in src/act/testdata
    std::string alexandriaFile(const std::string &name)
    {
        return fileManager().getInputFilePath(name);
    }

    /*! \brief Run the minimize tool in the temporary directory of the tests
     * \param[in]  molFile  Molecule file in src/act/testdata/molecules, may be empty
     * \param[in]  ffFile   Force field file in src/act/testdata/forcefields
     * \param[in]  label    Unique label for the output files
     * \param[in]  extra    Extra command line options
     * \param[in]  jsonFile Name for the JSON file. If empty, the option
     *                      -json is not given.
     * \return The exit status of the tool
     */
    int runMinimize(const std::string              &molFile,
                    const std::string              &ffFile,
                    const std::string              &label,
                    const std::vector<std::string> &extra,
                    const std::string              &jsonFile)
    {
        gmx::test::CommandLine cmd;
        cmd.append("minimize");
        cmd.addOption("-ff", alexandriaFile("forcefields/" + ffFile));
        // No file for a molecule that is given in another way
        if (!molFile.empty())
        {
            cmd.addOption("-f", alexandriaFile("molecules/" + molFile));
        }
        cmd.addOption("-c", fileManager().getTemporaryFilePath(label + "-confout.pdb"));
        cmd.addOption("-g", fileManager().getTemporaryFilePath(label + ".log"));
        if (!jsonFile.empty())
        {
            cmd.addOption("-json", jsonFile);
        }
        for(const auto &e : extra)
        {
            cmd.append(e);
        }
        const ScopedDirectory cwd(fileManager().getOutputTempDirectory());
        return minimize(cmd.argc(), cmd.argv());
    }

    /*! \brief Run the minimize tool and read its JSON output
     * \param[in]  molFile Molecule file in src/act/testdata/molecules
     * \param[in]  ffFile  Force field file in src/act/testdata/forcefields
     * \param[in]  label   Unique label for the output files
     * \param[out] json    The content below the key minimize
     */
    void minimizeToJson(const std::string &molFile,
                        const std::string &ffFile,
                        const std::string &label,
                        nlohmann::json    *json)
    {
        const auto jsonFile = fileManager().getTemporaryFilePath(label + ".json");
        ASSERT_EQ(0, runMinimize(molFile, ffFile, label, {}, jsonFile));
        std::ifstream in(jsonFile);
        ASSERT_TRUE(in.good()) << "Cannot open " << jsonFile;
        const auto tree = nlohmann::json::parse(in, nullptr, false);
        ASSERT_TRUE(tree.is_object()) << jsonFile << " is not valid JSON";
        ASSERT_TRUE(tree.contains("minimize"));
        *json = tree["minimize"];
    }

    /*! \brief Run the minimize tool and compare what it finds and writes with the reference data
     *
     * These tests fix the results of the tool with the given options, such
     * that a change in how the options reach the minimization shows.
     * \param[in] molFile Molecule file in src/act/testdata/molecules
     * \param[in] ffFile  Force field file in src/act/testdata/forcefields
     * \param[in] label   Unique label for the output files
     * \param[in] extra   Extra command line options
     * \param[in] exitStatus What the tool returns, 1 if the minimization does not converge
     */
    void checkMinimization(const std::string              &molFile,
                           const std::string              &ffFile,
                           const std::string              &label,
                           const std::vector<std::string> &extra,
                           int                             exitStatus = 0)
    {
        const auto jsonFile = fileManager().getTemporaryFilePath(label + ".json");
        ASSERT_EQ(exitStatus, runMinimize(molFile, ffFile, label, extra, jsonFile));
        const auto tree = readJsonOutput(jsonFile);
        ASSERT_TRUE(tree.is_object() && tree.contains("minimize")) << jsonFile;
        const auto &min = tree["minimize"];
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        checkJson(&checker, min["status"], "status");
        checkJson(&checker, min["converged"], "converged");
        // The structures in single precision differ by about 1e-5 nm
        checker.setDefaultTolerance(buildTolerance(2e-4, 2e-6, 1e-4));
        checkJson(&checker, min["rmsd"], "rmsd");
        // The energies are written with four decimals. In single precision
        // they differ by the precision of the largest term.
        double largest = 0;
        for (const auto &when : { "Energies before", "Energies after" })
        {
            for (const auto &term : min[when].items())
            {
                largest = std::max(largest, std::abs(term.value()["value"].get<double>()));
            }
        }
        checker.setDefaultTolerance(buildTolerance(1e-3 + 2e-6*largest, 2e-4 + 1e-8*largest, 1e-5));
        checkJson(&checker, min["Energies before"], "Energies before");
        checkJson(&checker, min["Energies after"], "Energies after");
        // The charges depend on the structure
        checker.setDefaultTolerance(buildTolerance(2e-4, 2e-6, 1e-5));
        checkJson(&checker, min["molecule"], "molecule");
        checker.setDefaultTolerance(buildTolerance(2e-4, 1e-5, 1e-5));
        checkJson(&checker, min["coordinates"], "coordinates");
        // The structure file, in Angstrom with three decimals, without
        // the superposition of the coordinates in the JSON output, so with
        // the small motion of the whole molecule in single precision
        checker.setDefaultTolerance(buildTolerance(2e-2, 2e-3, 0));
        checkTextFile(&checker, fileManager().getTemporaryFilePath(label + "-confout.pdb"), {}, "confout");
    }
};

//! \return the root mean square difference of two sets of coordinates, without fitting
double plainRmsd(const nlohmann::json &a, const nlohmann::json &b)
{
    double sum = 0;
    for(size_t i = 0; i < a.size(); i++)
    {
        for(int m = 0; m < 3; m++)
        {
            const double d = a[i][m].get<double>() - b[i][m].get<double>();
            sum += d*d;
        }
    }
    return a.empty() ? 0 : std::sqrt(sum/a.size());
}

TEST_F(MinimizeTest, JsonHasStatusRmsdAndEnergies)
{
    nlohmann::json js;
    minimizeToJson("acetone.sdf", "ACS-g1.xml", "status", &js);
    ASSERT_TRUE(js.is_object());

    ASSERT_TRUE(js.contains("status"));
    EXPECT_EQ("OK", js["status"].get<std::string>());
    ASSERT_TRUE(js.contains("converged"));
    ASSERT_TRUE(js["converged"].is_boolean());
    EXPECT_TRUE(js["converged"].get<bool>());

    ASSERT_TRUE(js.contains("rmsd"));
    ASSERT_TRUE(js["rmsd"]["value"].is_number());
    EXPECT_EQ("nm", js["rmsd"]["unit"].get<std::string>());
    // The structure changes, but an optimized input structure not by much
    EXPECT_GT(js["rmsd"]["value"].get<double>(), 0);
    EXPECT_LT(js["rmsd"]["value"].get<double>(), 0.1);

    ASSERT_TRUE(js.contains("Energies before"));
    ASSERT_TRUE(js.contains("Energies after"));
    const auto &before = js["Energies before"]["EPOT"];
    const auto &after  = js["Energies after"]["EPOT"];
    ASSERT_TRUE(before["value"].is_number());
    ASSERT_TRUE(after["value"].is_number());
    EXPECT_EQ("kJ/mol", after["unit"].get<std::string>());
    EXPECT_LT(after["value"].get<double>(), before["value"].get<double>());
}

TEST_F(MinimizeTest, JsonDescribesTheMolecule)
{
    nlohmann::json js;
    minimizeToJson("acetone.sdf", "ACS-g1.xml", "molecule", &js);
    ASSERT_TRUE(js.is_object());
    ASSERT_TRUE(js.contains("molecule"));
    const auto &mol = js["molecule"];
    // The name may be empty, since the file need not have one
    EXPECT_TRUE(mol["name"].is_string());
    EXPECT_EQ("C3H6O", mol.value("formula", ""));
    ASSERT_TRUE(mol["charge"].is_number_integer());
    EXPECT_EQ(0, mol["charge"].get<int>());
    ASSERT_TRUE(mol["mass"]["value"].is_number());
    EXPECT_NEAR(58.08, mol["mass"]["value"].get<double>(), 0.01);

    // Acetone has ten atoms and this force field has no shells
    ASSERT_TRUE(mol["particles"].is_array());
    ASSERT_EQ(10U, mol["particles"].size());
    double                     qtot = 0;
    std::map<std::string, int> elements;
    for (const auto &p : mol["particles"])
    {
        EXPECT_FALSE(p.value("name", "").empty());
        EXPECT_FALSE(p.value("type", "").empty());
        EXPECT_EQ("Atom", p.value("particletype", ""));
        ASSERT_TRUE(p["atomicnumber"].is_number_integer());
        ASSERT_TRUE(p["mass"].is_number());
        EXPECT_GT(p["mass"].get<double>(), 1.0);
        ASSERT_TRUE(p["charge"].is_number());
        qtot += p["charge"].get<double>();
        elements[p.value("element", "")] += 1;
    }
    EXPECT_NEAR(0.0, qtot, 1e-4);
    EXPECT_EQ(3, elements["C"]);
    EXPECT_EQ(6, elements["H"]);
    EXPECT_EQ(1, elements["O"]);

    // Nine bonds, of which one is a double bond between carbon and oxygen
    ASSERT_TRUE(mol["bonds"].is_array());
    ASSERT_EQ(9U, mol["bonds"].size());
    int ndouble = 0;
    for (const auto &b : mol["bonds"])
    {
        ASSERT_TRUE(b["atoms"].is_array());
        ASSERT_EQ(2U, b["atoms"].size());
        const int ai = b["atoms"][0].get<int>();
        const int aj = b["atoms"][1].get<int>();
        ASSERT_TRUE(ai >= 0 && ai < 10 && aj >= 0 && aj < 10 && ai != aj);
        ASSERT_TRUE(b["order"].is_number());
        if (b["order"].get<double>() == 2)
        {
            ndouble += 1;
            const std::set<std::string> co = { mol["particles"][ai].value("element", ""),
                                               mol["particles"][aj].value("element", "") };
            EXPECT_EQ(std::set<std::string>({ "C", "O" }), co);
        }
    }
    EXPECT_EQ(1, ndouble);
}

TEST_F(MinimizeTest, JsonHasCoordinatesBeforeAndAfter)
{
    nlohmann::json js;
    minimizeToJson("acetone.sdf", "ACS-g1.xml", "coordinates", &js);
    ASSERT_TRUE(js.is_object());
    ASSERT_TRUE(js.contains("coordinates"));
    const auto &x = js["coordinates"];
    EXPECT_EQ("nm", x.value("unit", ""));
    ASSERT_TRUE(x["before"].is_array());
    ASSERT_TRUE(x["after"].is_array());
    ASSERT_EQ(js["molecule"]["particles"].size(), x["before"].size());
    ASSERT_EQ(x["before"].size(), x["after"].size());
    for (const auto *xx : { &x["before"], &x["after"] })
    {
        for (const auto &xi : *xx)
        {
            ASSERT_TRUE(xi.is_array());
            ASSERT_EQ(3U, xi.size());
            ASSERT_TRUE(xi[0].is_number() && xi[1].is_number() && xi[2].is_number());
        }
    }
    // The first bond has a sensible length in both structures
    const auto &bond = js["molecule"]["bonds"][0]["atoms"];
    for (const auto *xx : { &x["before"], &x["after"] })
    {
        double d2 = 0;
        for(int m = 0; m < 3; m++)
        {
            const double d = (*xx)[bond[0].get<int>()][m].get<double>() - (*xx)[bond[1].get<int>()][m].get<double>();
            d2 += d*d;
        }
        EXPECT_GT(std::sqrt(d2), 0.09);
        EXPECT_LT(std::sqrt(d2), 0.16);
    }
    // The final structure is superimposed on the initial one, such that the
    // RMSD in the file is that of the coordinates as they are written.
    const double rmsd = js["rmsd"]["value"].get<double>();
    EXPECT_GT(rmsd, 0);
    EXPECT_NEAR(rmsd, plainRmsd(x["before"], x["after"]), 1e-5);
}

TEST_F(MinimizeTest, JsonOfPolarizableModelHasShells)
{
    nlohmann::json js;
    minimizeToJson("acetone.sdf", "ACS-pg1.xml", "shells", &js);
    ASSERT_TRUE(js.is_object());
    const auto &particles = js["molecule"]["particles"];
    ASSERT_TRUE(particles.is_array());
    int natom = 0, nshell = 0;
    for (size_t i = 0; i < particles.size(); i++)
    {
        const auto &p = particles[i];
        if (p.value("particletype", "") == "Atom")
        {
            natom += 1;
        }
        else if (p.value("particletype", "") == "Shell")
        {
            nshell += 1;
            // A shell knows the atom it belongs to
            ASSERT_TRUE(p["parents"].is_array());
            ASSERT_EQ(1U, p["parents"].size());
            const int parent = p["parents"][0].get<int>();
            ASSERT_TRUE(parent >= 0 && parent < static_cast<int>(particles.size()));
            EXPECT_EQ("Atom", particles[parent].value("particletype", ""));
        }
    }
    EXPECT_EQ(10, natom);
    EXPECT_EQ(10, nshell);
    EXPECT_EQ(particles.size(), js["coordinates"]["after"].size());
}

TEST_F(MinimizeTest, JsonIsNotWrittenWithoutTheOption)
{
    {
        const ScopedDirectory cwd(fileManager().getOutputTempDirectory());
        std::filesystem::remove("minimize.json");
    }
    ASSERT_EQ(0, runMinimize("acetone.sdf", "ACS-g1.xml", "nojson", {}, ""));
    const ScopedDirectory cwd(fileManager().getOutputTempDirectory());
    EXPECT_FALSE(std::filesystem::exists("minimize.json"));
}

TEST_F(MinimizeTest, MoleculeFromSmiles)
{
    // Acetone from its SMILES string instead of from a file
    const auto jsonFile = fileManager().getTemporaryFilePath("smiles.json");
    ASSERT_EQ(0, runMinimize("", "ACS-g1.xml", "smiles", { "-smiles", "CC(=O)C" }, jsonFile));
    std::ifstream in(jsonFile);
    ASSERT_TRUE(in.good());
    const auto js = nlohmann::json::parse(in, nullptr, false);
    ASSERT_TRUE(js.is_object());
    const auto &min = js["minimize"];
    EXPECT_EQ("C3H6O", min["molecule"].value("formula", ""));
    EXPECT_EQ(10U, min["molecule"]["particles"].size());
    EXPECT_EQ(9U, min["molecule"]["bonds"].size());
    EXPECT_TRUE(min["converged"].get<bool>());
    // The generated structure is close to a minimum of the force field
    EXPECT_LT(min["rmsd"]["value"].get<double>(), 0.1);

    // About the same energy as from the file. The minimum can be another
    // one, since the methyl groups can be rotated, which makes a
    // difference of about 1 kJ/mol.
    nlohmann::json fromFile;
    minimizeToJson("acetone.sdf", "ACS-g1.xml", "smilesref", &fromFile);
    ASSERT_TRUE(fromFile.is_object());
    EXPECT_NEAR(fromFile["Energies after"]["EPOT"]["value"].get<double>(),
                min["Energies after"]["EPOT"]["value"].get<double>(), 5.0);
}

// The results of the tool for several sets of options, as reference data

TEST_F(MinimizeTest, ResultsDefaults)
{
    checkMinimization("acetone.sdf", "ACS-g1.xml", "results-defaults", {});
}

TEST_F(MinimizeTest, ResultsPolarizable)
{
    checkMinimization("acetone.sdf", "ACS-pg1.xml", "results-polarizable", {});
}

TEST_F(MinimizeTest, ResultsNewton)
{
    checkMinimization("acetone.sdf", "ACS-g1.xml", "results-newton",
                      { "-minalg", "Newton", "-toler", "1e-4" });
}

TEST_F(MinimizeTest, ResultsSteepestDescentLimited)
{
    // Does not converge in 5 steps, which the tool tells with its exit
    // status, and writes the structure with the lowest energy. After more
    // steps the path in single precision differs too much from that in
    // double precision.
    checkMinimization("acetone.sdf", "ACS-g1.xml", "results-steep",
                      { "-minalg", "Steep", "-maxiter", "5" }, 1);
}

TEST_F(MinimizeTest, ResultsConstrained)
{
    checkMinimization("acetone.sdf", "ACS-g1.xml", "results-constrained",
                      { "-constrain", "bonds" });
}

TEST_F(MinimizeTest, ResultsChargesOfTypes)
{
    checkMinimization("acetone.sdf", "ACS-g1.xml", "results-qnone",
                      { "-qalg", "None" });
}

TEST_F(MinimizeTest, WrongSmilesFails)
{
    const auto jsonFile = fileManager().getTemporaryFilePath("wrongsmiles.json");
    EXPECT_EQ(1, runMinimize("", "ACS-g1.xml", "wrongsmiles", { "-smiles", "C(C" }, jsonFile));
}

TEST_F(MinimizeTest, SmilesAndFileTogetherFail)
{
    const auto jsonFile = fileManager().getTemporaryFilePath("both.json");
    EXPECT_EQ(1, runMinimize("acetone.sdf", "ACS-g1.xml", "both", { "-smiles", "CC(=O)C" }, jsonFile));
}

TEST_F(MinimizeTest, MissingMoleculeFileFails)
{
    const auto jsonFile = fileManager().getTemporaryFilePath("missing.json");
    EXPECT_EQ(1, runMinimize("does-not-exist.sdf", "ACS-g1.xml", "missing", {}, jsonFile));
}

} // namespace

} // namespace alexandria
