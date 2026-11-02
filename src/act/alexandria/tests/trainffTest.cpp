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

/*! \internal \brief
 * Tests of the train_ff module with the MCMC optimizer on one core, using
 * input derived from the benchmarks in ACTdata/Performance: the
 * intermolecular and intramolecular training of alcohols and water, and
 * the hydrogen fluoride dimer with virtual sites.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "trainffTestUtil.h"

#include "gromacs/utility/exceptions.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

using trainfftest::TrainOutput;
using Options = std::vector<std::pair<std::string, std::string>>;

//! Options for MCMC on a single core, as in the benchmarks
Options mcmcOptions(int maxIter, const std::string &seed = "219463979")
{
    return {
        { "-seed", seed },
        { "-v", "3" },
        { "-optimizer", "MCMC" },
        { "-pop_size", "1" },
        { "-max_iter", std::to_string(maxIter) },
        { "-temp", "0.025" },
        { "-anneal_begin", "0.1" },
        { "-anneal_end", "0.9" },
        { "-sensitivity", "" },
        { "-printSP", "" },
        { "-fc_bound", "10000" },
        { "-ener_boltz_temp", "3000" }
    };
}

//! Intermolecular terms of alcohol and water dimers
Options interOptions()
{
    return {
        { "-fit", "sigma epsilon" },
        { "-fc_allelec", "1" },
        { "-fc_exch", "1" },
        { "-fc_disp", "1" }
    };
}

//! Intramolecular terms of alcohols
Options intraOptions()
{
    return {
        { "-fit", "De b c bondenergy kb bondlength angle kt c1 c2 c3 c4 c5" },
        { "-fc_epot", "1" },
        { "-fc_force", "0.1" }
    };
}

//! Electrostatics, induction and van der Waals of hydrogen fluoride
Options hfOptions()
{
    return {
        { "-fit", "charge eta chi zeta vs2a sigma epsilon gamma delta" },
        { "-fc_allelec", "1" },
        { "-fc_exch", "1" },
        { "-fc_disp", "1" },
        { "-fc_mu", "1" },
        { "-fc_quad", "1" }
    };
}

//! \return the options combined
Options combine(Options a, const Options &b)
{
    a.insert(a.end(), b.begin(), b.end());
    return a;
}

class TrainFFTest : public gmx::test::CommandLineTestBase
{
protected:
    /*! \brief Run train_ff and check the output
     * \param[in] ffName  Force field file
     * \param[in] molName Base name of the molprop and selection file
     * \param[in] options The options
     */
    void runAndCheck(const std::string &ffName,
                     const std::string &molName,
                     const Options     &options)
    {
        auto output = trainfftest::runTrainFF(&fileManager(), "train", ffName, molName, options);
        ASSERT_EQ(0, output.status);
        // MCMC never accepts a worse best parameter vector
        ASSERT_GT(output.initialChi2, 0);
        ASSERT_GT(output.finalChi2, 0);
        EXPECT_LE(output.finalChi2, output.initialChi2);
        // The trained force field is written
        EXPECT_NE(std::string::npos, output.forceField.find("</alexandria_chemistry>"));
        gmx::test::TestReferenceChecker checker(rootChecker());
        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-4));
        trainfftest::checkOutput(&checker, output);
    }
};

// Test: Training Lennard-Jones parameters on alcohol and water dimers
TEST_F(TrainFFTest, MCMCInterBooze)
{
    runAndCheck("lj-pc.xml", "booze-inter", combine(mcmcOptions(10), interOptions()));
}

// Test: Training bonded parameters on the energies and forces of alcohols
TEST_F(TrainFFTest, MCMCIntraBooze)
{
    runAndCheck("halgren-PG.xml", "booze-mono", combine(mcmcOptions(3), intraOptions()));
}

// Test: Training the hydrogen fluoride dimer, with virtual sites, on the
// interaction energy components and the multipoles of the monomer
TEST_F(TrainFFTest, MCMCHydrogenFluoride)
{
    runAndCheck("4S-lj14-7-pg-hvdw.xml", "hf-hf", combine(mcmcOptions(3), hfOptions()));
}

// Test: The files with the correlation of the energy components get the
// prefix given with -compound_file
TEST_F(TrainFFTest, CompoundFilePrefix)
{
    auto prefix  = fileManager().getTemporaryFilePath("hf-");
    auto options = combine(mcmcOptions(1), hfOptions());
    options.push_back({ "-compound_file", prefix });
    auto output = trainfftest::runTrainFF(&fileManager(), "prefix", "4S-lj14-7-pg-hvdw.xml", "hf-hf", options);
    ASSERT_EQ(0, output.status);
    for (const std::string term : { "ALLELEC", "DISPERSION", "EXCHANGE" })
    {
        EXPECT_TRUE(std::filesystem::exists(prefix + term + ".xvg")) << prefix + term + ".xvg";
    }
}

// Test: In verbose mode MCMC writes the convergence of each class of
// parameters and of the chi-squared, one line per iteration
TEST_F(TrainFFTest, ConvergenceFiles)
{
    auto options = combine(mcmcOptions(3), intraOptions());
    for (auto &o : options)
    {
        if (o.first == "-v")
        {
            o.second = "4";
        }
    }
    auto output = trainfftest::runTrainFF(&fileManager(), "conv", "halgren-PG.xml", "booze-mono", options);
    ASSERT_EQ(0, output.status);
    // The files are in the directory of the individual, named after the
    // file names given
    auto dir = std::filesystem::path(fileManager().getTemporaryFilePath("conv.log")).parent_path() / "inds";
    auto convName = std::filesystem::path(fileManager().getTemporaryFilePath("conv-conv.xvg")).filename().string();
    auto chi2Name = std::filesystem::path(fileManager().getTemporaryFilePath("conv-chi2.xvg")).filename().string();
    // The parameters of halgren-PG.xml that are trained
    std::vector<std::string> files = { "ind0-" + chi2Name };
    for (const char *pclass : { "De", "b", "c", "bondlength", "angle", "kt", "c1", "c2", "c3", "c4", "c5" })
    {
        files.push_back("ind0-" + std::string(pclass) + "-" + convName);
    }
    int nlines = -1;
    for (const auto &f : files)
    {
        std::ifstream in(dir / f);
        ASSERT_TRUE(in.good()) << f;
        int n = 0;
        for (std::string line; std::getline(in, line); )
        {
            if (!line.empty() && line[0] != '#' && line[0] != '@')
            {
                n++;
            }
        }
        // The same steps in each file
        EXPECT_GT(n, 0) << f;
        if (nlines >= 0)
        {
            EXPECT_EQ(nlines, n) << f;
        }
        nlines = n;
    }
    std::filesystem::remove_all(dir);
}

// Test: With -checkpoint MCMC writes the force field of each new minimum in
// the directory of the individual, named after the output file
TEST_F(TrainFFTest, Checkpoint)
{
    auto options = combine(mcmcOptions(3), intraOptions());
    options.push_back({ "-checkpoint", "" });
    auto output = trainfftest::runTrainFF(&fileManager(), "checkpoint", "halgren-PG.xml", "booze-mono", options);
    ASSERT_EQ(0, output.status);
    auto dir  = std::filesystem::path(fileManager().getTemporaryFilePath("checkpoint.log")).parent_path() / "inds";
    auto name = std::filesystem::path(fileManager().getTemporaryFilePath("checkpoint-ff.xml")).filename().string();
    auto file = dir / ("ind0-" + name);
    std::ifstream in(file);
    ASSERT_TRUE(in.good()) << file;
    std::stringstream ff;
    ff << in.rdbuf();
    EXPECT_NE(std::string::npos, ff.str().find("</alexandria_chemistry>"));
    std::filesystem::remove_all(dir);
}

// Test: The outliers are listed with their structure, from the largest
// to the smallest deviation
TEST_F(TrainFFTest, OutlierOverview)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "outliers", "lj-pc.xml", "booze-inter",
                                          combine(mcmcOptions(1), interOptions()));
    ASSERT_EQ(0, output.status);
    int  noutlier = 0;
    bool inTable  = false;
    double previous = 0;
    for (const auto &line : output.log)
    {
        if (line.find("Name") == 0 && line.find("Structure") != std::string::npos)
        {
            inTable  = true;
            previous = std::numeric_limits<double>::max();
            continue;
        }
        if (inTable)
        {
            std::istringstream iss(line);
            std::string        name, structure;
            double             ref, act, diff;
            if (!(iss >> name >> structure >> ref >> act >> diff))
            {
                inTable = false;
                continue;
            }
            EXPECT_NEAR(act - ref, diff, 1e-3*std::abs(diff) + 1e-3) << line;
            EXPECT_LE(std::abs(diff), previous) << line;
            previous  = std::abs(diff);
            noutlier += 1;
        }
    }
    EXPECT_GT(noutlier, 1);
}

// Test: The same seed gives the same parameters, a different one does not
TEST_F(TrainFFTest, Reproducible)
{
    std::vector<nlohmann::json> params;
    int                         run = 0;
    for (const std::string seed : { "1993", "1993", "1997" })
    {
        auto options = combine(mcmcOptions(3, seed), hfOptions());
        auto output = trainfftest::runTrainFF(&fileManager(), "seed" + std::to_string(run++),
                                              "4S-lj14-7-pg-hvdw.xml", "hf-hf", options);
        ASSERT_EQ(0, output.status);
        params.push_back(output.json["train_ff"]["parameters"]);
    }
    EXPECT_EQ(params[0], params[1]);
    EXPECT_NE(params[0], params[2]);
}

//! \return the number of points and the legends of the sets in an xvg file
static std::pair<int, std::vector<std::string>> xvgSets(const std::string &fileName)
{
    std::pair<int, std::vector<std::string>> sets = { 0, {} };
    for (const auto &line : trainfftest::readLines(fileName))
    {
        auto legend = line.find("legend \"");
        if (legend != std::string::npos)
        {
            auto begin = legend + std::string("legend \"").size();
            sets.second.push_back(line.substr(begin, line.rfind('"') - begin));
        }
        else if (!line.empty() && line[0] != '@' && line[0] != '#' && line[0] != '&')
        {
            sets.first += 1;
        }
    }
    return sets;
}

// Test: The frequencies of the optimized methanol in the training set are
// compared with the reference, and the forces of both sets, without
// optimizing. The minimized structures are written in a directory that
// is made for optimizing.
TEST_F(TrainFFTest, FrequenciesAndForces)
{
    auto freq   = fileManager().getTemporaryFilePath("freq.xvg");
    auto force  = fileManager().getTemporaryFilePath("force.xvg");
    auto output = trainfftest::runTrainFF(&fileManager(), "freq", "lj-pc-trained.xml", "booze-mono",
                                          { { "-nooptimize", "" }, { "-calc_frequencies", "" },
                                            { "-fit", "kb" }, { "-fc_epot", "1" }, { "-fc_force", "1" },
                                            { "-freqcorr", freq }, { "-forcecorr", force } });
    ASSERT_EQ(0, output.status);
    // The 12 vibrations of methanol
    auto freqSets = xvgSets(freq);
    EXPECT_EQ(12, freqSets.first);
    EXPECT_EQ(std::vector<std::string>({ "Train" }), freqSets.second);
    auto forceSets = xvgSets(force);
    EXPECT_GT(forceSets.first, 0);
    EXPECT_EQ(std::vector<std::string>({ "Train", "Test" }), forceSets.second);
    // Each molecule has a frequency analysis
    const auto &molecules = output.json["train_ff"]["molecules"];
    EXPECT_EQ(2u, molecules.size());
    for (const auto &mol : molecules)
    {
        EXPECT_TRUE(mol.contains("FrequencyAnalysis")) << mol["name"];
    }
}

// Test: Without optimizing, structures that change much in the minimization
// are written as well. With the crude starting force field methanol falls
// apart.
TEST_F(TrainFFTest, MinimizedStructures)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "minimized", "halgren-PG.xml", "booze-mono",
                                          { { "-nooptimize", "" }, { "-calc_frequencies", "" },
                                            { "-fit", "kb" }, { "-fc_epot", "1" } });
    ASSERT_EQ(0, output.status);
    // The structures are written in the working directory of train_ff
    auto dir = std::filesystem::path(fileManager().getTemporaryFilePath("minimized.log")).parent_path() / "inds";
    auto pdb = dir / "methanol-original-minimized.pdb";
    EXPECT_TRUE(std::filesystem::exists(pdb)) << pdb;
    std::filesystem::remove_all(dir);
}

// Test: Compounds with parameters that have too little training data are
// removed. Without methanol#water the water parameters have no data in the
// training set, after which the support for methanol#methanol is checked
// again.
TEST_F(TrainFFTest, MinimumData)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "mindata", "lj-pc.xml", "booze-inter",
                                          combine(combine(mcmcOptions(1), interOptions()),
                                                  { { "-mindata", "2" } }));
    ASSERT_EQ(0, output.status);
    auto hasLine = [&output](const std::string &text)
    {
        return std::any_of(output.log.begin(), output.log.end(),
                           [&text](const std::string &line) { return line.find(text) != std::string::npos; });
    };
    EXPECT_TRUE(hasLine("Found 2 molecules without sufficient support"));
    EXPECT_TRUE(hasLine("There are 1 molecules left to optimize"));
}

// Test: Without support for any compound the training set is empty
TEST_F(TrainFFTest, MinimumDataTooLarge)
{
    EXPECT_THROW(trainfftest::runTrainFF(&fileManager(), "mindata", "lj-pc.xml", "booze-inter",
                                         combine(combine(mcmcOptions(1), interOptions()),
                                                 { { "-mindata", "3" } })),
                 gmx::InvalidInputError);
}

// Test: For a diatomic compound the reference energy minus the
// electrostatic energy is written as a function of the distance, one
// point for each structure
TEST_F(TrainFFTest, Diatomics)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "diatomics", "4S-lj14-7-pg-hvdw.xml", "hf-hf",
                                          { { "-nooptimize", "" }, { "-diatomics", "" },
                                            { "-fit", "sigma epsilon" }, { "-fc_epot", "1" } });
    ASSERT_EQ(0, output.status);
    // The file is named after the compound, in the working directory
    auto fileName = (std::filesystem::path(fileManager().getTemporaryFilePath("x")).parent_path() /
                     "hydrogen-fluoride.xvg").string();
    auto lines    = trainfftest::readLines(fileName);
    int  npoints  = 0;
    for (const auto &line : lines)
    {
        if (line.empty() || line[0] == '#' || line[0] == '@')
        {
            continue;
        }
        std::istringstream iss(line);
        double             r, e;
        ASSERT_TRUE(static_cast<bool>(iss >> r >> e)) << line;
        // The bond lengths of the structures, in nm
        EXPECT_GT(r, 0.07);
        EXPECT_LT(r, 0.11);
        EXPECT_TRUE(std::isfinite(e));
        npoints += 1;
    }
    EXPECT_EQ(10, npoints);
    std::filesystem::remove(fileName);
}

// Test: The electrostatic potential, multipoles and polarizabilities of
// molecules with ESP data are analysed, the molecules with the energies
// written to a molprop file, and the structures of all energies dumped
// even though their data files include a directory
TEST_F(TrainFFTest, ChargesEspPolarizability)
{
    auto selFile = fileManager().getTemporaryFilePath("acm.dat");
    {
        std::ofstream sel(selFile);
        sel << "1-butanol|Train\nacetone|Train\nwater|Test\n";
    }
    std::map<std::string, std::string> files = {
        { "-alphacorr", "alpha.xvg" }, { "-isopol", "isopol.xvg" }, { "-anisopol", "anisopol.xvg" },
        { "-espcorr", "esp.xvg" }, { "-qhisto", "qhisto.xvg" }, { "-mpout", "mpout.xml" }
    };
    Options options = {
        { "-nooptimize", "" }, { "-fit", "chi eta" }, { "-fc_esp", "1" }, { "-fc_polar", "1" },
        { "-fc_mu", "1" }, { "-fc_epot", "1" }, { "-dump_outliers", "0" }
    };
    for (auto &f : files)
    {
        f.second = fileManager().getTemporaryFilePath(f.second);
        options.push_back(f);
    }
    auto output = trainfftest::runTrainFF(&fileManager(), "acm", "../ACS-pg.xml", "acm", options,
                                          fileManager().getInputFilePath("../../molprop/tests/molprop.xml"),
                                          selFile);
    ASSERT_EQ(0, output.status);
    for (const auto &f : files)
    {
        ASSERT_TRUE(std::filesystem::exists(f.second)) << f.first;
        EXPECT_GT(std::filesystem::file_size(f.second), 0u) << f.first;
    }
    int nesp = 0, nalpha = 0;
    for (const auto &line : output.log)
    {
        nesp   += line.find("ESP rms:") == 0;
        nalpha += line.find("Isotropic polarizability:") == 0;
    }
    // One potential for each compound, the polarizability for each of the
    // two levels of theory for butanol and water, one for acetone
    EXPECT_EQ(3, nesp);
    EXPECT_EQ(5, nalpha);
    auto mpout = trainfftest::readLines(files["-mpout"]);
    EXPECT_EQ(3, std::count_if(mpout.begin(), mpout.end(),
                               [](const std::string &l) { return l.find("<molecule ") != std::string::npos; }));
    // One structure for each energy, with the atoms
    auto dir    = std::filesystem::path(files["-mpout"]).parent_path();
    int  nxyz   = 0;
    for (const auto &entry : std::filesystem::directory_iterator(dir))
    {
        auto name = entry.path().filename().string();
        if (entry.path().extension() == ".xyz" && name.find("-EPOT-") != std::string::npos)
        {
            auto xyz = trainfftest::readLines(entry.path().string());
            ASSERT_FALSE(xyz.empty());
            EXPECT_EQ(xyz.size(), 2 + std::stoul(xyz[0])) << name;
            std::filesystem::remove(entry.path());
            nxyz += 1;
        }
    }
    EXPECT_EQ(5, nxyz);
}

} // namespace

} // namespace alexandria
