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
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "trainffTestUtil.h"

#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/particletype.h"
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
        { "-sensitivity", "1D" },
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

/*! \brief Set the value of an option
 * \param[in] options The options
 * \param[in] name    The option
 * \param[in] value   Its value
 * \return the options with the value replaced, or added if not there
 */
Options setOption(Options options, const std::string &name, const std::string &value)
{
    for (auto &opt : options)
    {
        if (opt.first == name)
        {
            opt.second = value;
            return options;
        }
    }
    options.emplace_back(name, value);
    return options;
}

/*! \brief Options to analyse a force field without training it
 * The Lennard-Jones parameters of the trained force field are analysed on
 * the alcohol and water dimers, see analysedForceField. These are sigma
 * and epsilon of c3, o3 and ow, the ones of the hydrogens are fixed.
 * \param[in] mode The sensitivity analysis
 * \return the options
 */
Options analyseOnly(const std::string &mode)
{
    auto options = setOption(combine(mcmcOptions(1), interOptions()), "-sensitivity", mode);
    options.emplace_back("-nooptimize", "");
    options.emplace_back("-norandom_init", "");
    return options;
}

//! The trained force field that is analysed with analyseOnly
const char *analysedForceField = "lj-pc-trained.xml";

//! The compounds it is analysed on
const char *analysedCompounds  = "booze-inter";

//! \return whether a line of the log file contains a text
bool logContains(const std::vector<std::string> &log, const std::string &text)
{
    for (const auto &line : log)
    {
        if (line.find(text) != std::string::npos)
        {
            return true;
        }
    }
    return false;
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
        const auto output = trainfftest::runTrainFF(&fileManager(), "train", ffName, molName, options);
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

//! \return the parameters listed as trained in the log, as identifier and type
static std::set<std::pair<std::string, std::string>> trainedParameters(const TrainOutput &output)
{
    std::set<std::pair<std::string, std::string>> trained;
    const auto &log = output.log;
    for (size_t i = 0; i < log.size(); i++)
    {
        if (log[i].find("parameters to train.") == std::string::npos)
        {
            continue;
        }
        // After the header, a line per parameter with its bounds
        for (size_t j = i + 2; j < log.size(); j++)
        {
            std::istringstream iss(log[j]);
            std::string        id, type;
            double             minimum, maximum;
            if (!(iss >> id >> type >> minimum >> maximum))
            {
                break;
            }
            trained.insert({ id, type });
        }
        break;
    }
    return trained;
}

/*! \brief A parameter of a particle type or of an interaction
 * \param[in] pd   The force field
 * \param[in] id   The identifier
 * \param[in] type The parameter type
 * \return the parameter, nullptr if there is none
 */
static const ForceFieldParameter *findParameter(const ForceField  &pd,
                                                const std::string &id,
                                                const std::string &type)
{
    if (pd.hasParticleType(id) && pd.findParticleType(id)->hasParameter(type))
    {
        return &pd.findParticleType(id)->parameterConst(type);
    }
    for (const auto &fs : pd.forcesConst())
    {
        for (const auto &p : fs.second.parametersConst())
        {
            const auto pp = p.second.find(type);
            if (p.first.id() == id && p.second.end() != pp)
            {
                return &pp->second;
            }
        }
    }
    return nullptr;
}

/*! \brief Check that the trained parameters have data, stay within their
 * bounds and that some changed, and that the chi squared decreased
 * \param[in] fileManager For temporary files
 * \param[in] output      The output of train_ff
 * \param[in] ffFile      The force field that was trained
 * \param[in] expected    The parameters that are trained
 */
static void expectTrained(gmx::test::TestFileManager                          *fileManager,
                          const TrainOutput                                   &output,
                          const std::string                                   &ffFile,
                          const std::set<std::pair<std::string, std::string>> &expected)
{
    ASSERT_EQ(0, output.status);
    ASSERT_FALSE(output.forceField.empty()) << "No force field written";
    EXPECT_EQ(expected, trainedParameters(output));
    ForceField input, trained;
    readForceField(ffFile, &input);
    const auto trainedFile = fileManager->getTemporaryFilePath("trained.xml");
    {
        std::ofstream out(trainedFile);
        out << output.forceField;
    }
    readForceField(trainedFile, &trained);
    int nchanged = 0;
    for (const auto &e : expected)
    {
        auto p0 = findParameter(input, e.first, e.second);
        auto p1 = findParameter(trained, e.first, e.second);
        ASSERT_TRUE(p0 && p1) << e.first << " " << e.second;
        EXPECT_EQ(Mutability::Bounded, p1->mutability()) << e.first << " " << e.second;
        EXPECT_LT(0u, p1->ntrain()) << e.first << " " << e.second;
        EXPECT_LE(p1->minimum(), p1->value()) << e.first << " " << e.second;
        EXPECT_GE(p1->maximum(), p1->value()) << e.first << " " << e.second;
        nchanged += (p0->value() != p1->value()) ? 1 : 0;
    }
    EXPECT_LT(0, nchanged);
    ASSERT_GT(output.initialChi2, 0);
    EXPECT_LT(output.finalChi2, output.initialChi2);
}

//! \return a selection file with these compounds in the training set
static std::string trainSelection(gmx::test::TestFileManager *fileManager, const std::string &label,
                                  const std::vector<std::string> &compounds)
{
    auto selFile = fileManager->getTemporaryFilePath(label + ".dat");
    std::ofstream sel(selFile);
    for (const auto &c : compounds)
    {
        sel << c << "|Train\n";
    }
    return selFile;
}

/*! \brief A chi squared component of the training set from the log
 * \param[in] output    The output of train_ff
 * \param[in] component The name of the component, as in the log
 * \param[in] best      Of the best parameters, otherwise of the initial ones
 * \return the chi squared, -1 if it is not in the log
 */
static double componentChi2(const TrainOutput &output, const std::string &component, bool best)
{
    const char *header = best ? "Chi2 components of the best parameter vector found (for Train)" :
        "initial parameter vector chi2 components";
    bool found = false;
    for (const auto &line : output.log)
    {
        if (line.find(header) != std::string::npos)
        {
            found = true;
        }
        else if (found && line.rfind(component + " ", 0) == 0)
        {
            std::istringstream iss(line);
            std::string        name;
            double             chi2;
            return (iss >> name >> chi2) ? chi2 : -1;
        }
        else if (found && line.rfind("TOT", 0) == 0)
        {
            break;
        }
    }
    return -1;
}

//! \return whether a chi squared component is in the log of the initial parameters
static bool hasInitialComponent(const TrainOutput &output, const std::string &component)
{
    return componentChi2(output, component, false) >= 0;
}

// Test: The widths of the charge distributions, which are parameters of the
// particle types, are trained, also those of shells, and reach the atoms
TEST_F(TrainFFTest, ParticleZeta)
{
    const auto ffFile = fileManager().getInputFilePath("forcefields/ACS-pol.xml");
    const auto output = trainfftest::runTrainFF(&fileManager(), "zeta", "", "booze-inter",
                                                combine(mcmcOptions(10),
                                                        { { "-fit", "zeta" }, { "-fc_allelec", "1" },
                                                    { "-fc_exch", "1" }, { "-fc_disp", "1" } }),
                                                "", "", ffFile);
    expectTrained(&fileManager(), output, ffFile,
                  { { "c3_s", "zeta" }, { "h1", "zeta" }, { "ho", "zeta" },
                    { "hw", "zeta" }, { "o3_s", "zeta" }, { "ow_s", "zeta" } });
}

// Test: Polarizabilities and the parameters of the induction correction are
// trained on the induction and its correction
TEST_F(TrainFFTest, PolarizabilityAndInductionCorrection)
{
    const auto ffFile = fileManager().getInputFilePath("forcefields/ACS-pol.xml");
    const auto output = trainfftest::runTrainFF(&fileManager(), "induction", "", "",
                                                combine(mcmcOptions(10),
                                                        { { "-fit", "alpha a1dexp a2dexp bdexp" },
                                                    { "-fc_induc", "1" }, { "-fc_deltahf", "1" } }),
                                                fileManager().getInputFilePath("molprop/methanol-ethane-deltaHF.xml"),
                                                trainSelection(&fileManager(), "induction", { "methanol#ethane" }),
                                                ffFile);
    std::set<std::pair<std::string, std::string>> expected = { { "c3_s", "alpha" }, { "o3_s", "alpha" } };
    for (const char *atype : { "c3", "h1", "hc", "ho", "o3" })
    {
        for (const char *ptype : { "a1dexp", "a2dexp", "bdexp" })
        {
            expected.insert({ atype, ptype });
        }
    }
    expectTrained(&fileManager(), output, ffFile, expected);
    EXPECT_TRUE(hasInitialComponent(output, "Induction"));
    EXPECT_TRUE(hasInitialComponent(output, "DeltaHF"));
}

// Test: The components of the interaction energy, electrostatics, exchange
// and induction together, and the total are training targets
TEST_F(TrainFFTest, InteractionComponents)
{
    const auto ffFile = fileManager().getInputFilePath("forcefields/ACS-pol.xml");
    const auto output = trainfftest::runTrainFF(&fileManager(), "components", "", "",
                                                combine(mcmcOptions(5),
                                                        { { "-fit", "zeta" }, { "-fc_elec", "1" },
                                                    { "-fc_exchind", "1" }, { "-fc_inter", "1" } }),
                                                fileManager().getInputFilePath("molprop/mpdimer.xml"),
                                                trainSelection(&fileManager(), "components",
                                                               { "1-propanol#ethane", "water#1-butanol", "guanidinium#water" }),
                                                ffFile);
    ASSERT_EQ(0, output.status);
    for (const char *component : { "Electrostatics", "ExchInd", "Interaction" })
    {
        EXPECT_TRUE(hasInitialComponent(output, component)) << component;
    }
    EXPECT_LT(0u, trainedParameters(output).size());
    ASSERT_GT(output.initialChi2, 0);
    EXPECT_LT(output.finalChi2, output.initialChi2);
}

// Test: With another loss function the chi squared of the same parameters
// differs, and training works with each of them
TEST_F(TrainFFTest, LossFunctions)
{
    std::map<std::string, double> initial;
    for (const char *lf : { "MSE", "MAE", "Huber", "Asinh" })
    {
        SCOPED_TRACE(lf);
        const auto output = trainfftest::runTrainFF(&fileManager(), std::string("loss") + lf, "lj-pc.xml", "booze-inter",
                                                    combine(combine(mcmcOptions(5), interOptions()),
                                                            { { "-loss_function", lf } }));
        ASSERT_EQ(0, output.status);
        ASSERT_GT(output.initialChi2, 0);
        EXPECT_LE(output.finalChi2, output.initialChi2);
        initial[lf] = output.initialChi2;
    }
    for (const char *lf : { "MAE", "Huber", "Asinh" })
    {
        EXPECT_NE(initial["MSE"], initial[lf]) << lf;
    }
}

// Test: The Morse bonds of the hydrogen halides are trained on their
// harmonic frequencies and intensities
TEST_F(TrainFFTest, FrequenciesAndIntensities)
{
    const auto output = trainfftest::runTrainFF(&fileManager(), "hhfreq", "hh-intra.xml", "hh-freq",
                                                combine(mcmcOptions(10),
                                                        { { "-fit", "De bondlength b c" },
                                                    { "-fc_freq", "1" }, { "-fc_inten", "1" } }));
    std::set<std::pair<std::string, std::string>> expected;
    for (const char *bond : { "br_b~hbr_b", "cl_b~hcl_b", "f_b~hf_b", "hi_b~i_b" })
    {
        for (const char *ptype : { "De", "bondlength", "b", "c" })
        {
            expected.insert({ bond, ptype });
        }
    }
    expectTrained(&fileManager(), output, fileManager().getInputFilePath("forcefields/hh-intra.xml"), expected);
    EXPECT_TRUE(hasInitialComponent(output, "FREQUENCY"));
    EXPECT_TRUE(hasInitialComponent(output, "INTENSITY"));
}

// Test: The Morse bonds, the ACM parameters, the charges and widths of the
// virtual sites and their positions of the hydrogen halides are trained on
// the dipoles, quadrupoles, harmonic frequencies and intensities, which all
// decrease
TEST_F(TrainFFTest, ElectrostaticsAndFrequencies)
{
    const auto output = trainfftest::runTrainFF(&fileManager(), "hhelec", "hh-intra.xml", "hh-freq",
                                                combine(mcmcOptions(10),
                                                        { { "-fit", "De bondlength b c eta chi chargeVS zeta zetaVS charge vs2a" },
                                                    { "-fc_mu", "1" }, { "-fc_quad", "5" },
                                                    { "-fc_freq", "1" }, { "-fc_inten", "1" } }));
    std::set<std::pair<std::string, std::string>> expected;
    for (const char *bond : { "br_b~hbr_b", "cl_b~hcl_b", "f_b~hf_b", "hi_b~i_b" })
    {
        for (const char *ptype : { "De", "bondlength", "b", "c" })
        {
            expected.insert({ bond, ptype });
        }
    }
    for (const char *vsite : { "br_b~hbr_b!v2br1", "cl_b~hcl_b!v2cl1", "f_b~hf_b!v2f1", "i_b~hi_b!v2i1" })
    {
        expected.insert({ vsite, "vs2a" });
    }
    for (const char *atom : { "br", "cl", "f", "i", "hbr", "hcl", "hf", "hi" })
    {
        expected.insert({ std::string(atom) + "_z", "chi" });
        expected.insert({ std::string(atom) + "_z", "eta" });
        // The virtual sites on top of the atoms
        expected.insert({ std::string("v1") + atom, "charge" });
        expected.insert({ std::string("v1") + atom, "zeta" });
    }
    for (const char *vsite : { "v2br1", "v2cl1", "v2f1", "v2i1" })
    {
        expected.insert({ vsite, "charge" });
    }
    expectTrained(&fileManager(), output, fileManager().getInputFilePath("forcefields/hh-intra.xml"), expected);
    for (const char *component : { "MU", "QUAD", "FREQUENCY", "INTENSITY" })
    {
        double initial = componentChi2(output, component, false);
        ASSERT_GT(initial, 0) << component;
        EXPECT_LT(componentChi2(output, component, true), initial) << component;
    }
}

// Test: The second charges of particles, which replace virtual sites on top
// of atoms when reading the force field, and the widths of their
// distributions are trained, see ForceField::convertVsite1 and
// https://github.com/AlexandriaChemistry/ACT/issues/870
TEST_F(TrainFFTest, SecondCharges)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "chargevs", "4S-lj14-7-pg-hvdw.xml", "hf-hf",
                                          combine(mcmcOptions(20),
                                                  { { "-fit", "chargeVS zetaVS" },
                                                    { "-fc_allelec", "1" }, { "-fc_mu", "1" } }));
    ASSERT_EQ(0, output.status);
    ASSERT_FALSE(output.forceField.empty()) << "No force field written";
    // The input force field has virtual sites, which become second charges
    ForceField input, trained;
    readForceField(fileManager().getInputFilePath("forcefields/4S-lj14-7-pg-hvdw.xml"), &input);
    const auto trainedFile = fileManager().getTemporaryFilePath("chargevs-trained.xml");
    {
        std::ofstream out(trainedFile);
        out << output.forceField;
    }
    readForceField(trainedFile, &trained);
    int nchanged = 0;
    for(const char *atype : { "f", "hf" })
    {
        ASSERT_TRUE(input.findParticleType(atype)->hasChargeVS()) << atype;
        for(const char *pname : { chargeVSName, zetaVSName })
        {
            // The parameter is listed with those that are trained
            EXPECT_TRUE(std::any_of(output.log.begin(), output.log.end(),
                                    [atype, pname](const std::string &line)
                                    {
                                        std::istringstream iss(line);
                                        std::string        id, type;
                                        return (iss >> id >> type) && id == atype && type == pname;
                                    })) << atype << " " << pname;
            const auto &p0 = input.findParticleType(atype)->parameterConst(pname);
            const auto &p1 = trained.findParticleType(atype)->parameterConst(pname);
            EXPECT_EQ(Mutability::Bounded, p1.mutability()) << atype << " " << pname;
            EXPECT_LT(0, p1.ntrain()) << atype << " " << pname;
            EXPECT_LE(p1.minimum(), p1.value()) << atype << " " << pname;
            EXPECT_GE(p1.maximum(), p1.value()) << atype << " " << pname;
            nchanged += (p0.value() != p1.value()) ? 1 : 0;
        }
    }
    EXPECT_LT(0, nchanged);
}

// Test: The second charges of particles and the widths of their
// distributions are trained with the charges and widths, see
// MolGen::fillIopt and https://github.com/AlexandriaChemistry/ACT/issues/1491
TEST_F(TrainFFTest, SecondChargesWithCharges)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "chargeall", "4S-lj14-7-pg-hvdw.xml", "hf-hf",
                                          combine(mcmcOptions(1),
                                                  { { "-fit", "charge zeta" },
                                                    { "-fc_allelec", "1" }, { "-fc_mu", "1" } }));
    ASSERT_EQ(0, output.status);
    for(const char *atype : { "f", "hf" })
    {
        for(const char *pname : { chargeVSName, zetaVSName })
        {
            // The parameter is listed with those that are trained
            EXPECT_TRUE(std::any_of(output.log.begin(), output.log.end(),
                                    [atype, pname](const std::string &line)
                                    {
                                        std::istringstream iss(line);
                                        std::string        id, type;
                                        return (iss >> id >> type) && id == atype && type == pname;
                                    })) << atype << " " << pname;
        }
    }
}

// Test: The files with the correlation of the energy components get the
// prefix given with -compound_file
TEST_F(TrainFFTest, CompoundFilePrefix)
{
    const auto prefix  = fileManager().getTemporaryFilePath("hf-");
    auto options = combine(mcmcOptions(1), hfOptions());
    options.emplace_back("-compound_file", prefix);
    const auto output = trainfftest::runTrainFF(&fileManager(), "prefix", "4S-lj14-7-pg-hvdw.xml", "hf-hf", options);
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
    const auto output = trainfftest::runTrainFF(&fileManager(), "conv", "halgren-PG.xml", "booze-mono", options);
    ASSERT_EQ(0, output.status);
    // The files are in the directory of the individual, named after the
    // file names given
    const auto dir = std::filesystem::path(fileManager().getTemporaryFilePath("conv.log")).parent_path() / "inds";
    const auto convName = std::filesystem::path(fileManager().getTemporaryFilePath("conv-conv.xvg")).filename().string();
    const auto chi2Name = std::filesystem::path(fileManager().getTemporaryFilePath("conv-chi2.xvg")).filename().string();
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
    options.emplace_back("-checkpoint", "");
    const auto output = trainfftest::runTrainFF(&fileManager(), "checkpoint", "halgren-PG.xml", "booze-mono", options);
    ASSERT_EQ(0, output.status);
    const auto dir  = std::filesystem::path(fileManager().getTemporaryFilePath("checkpoint.log")).parent_path() / "inds";
    const auto name = std::filesystem::path(fileManager().getTemporaryFilePath("checkpoint-ff.xml")).filename().string();
    const auto file = dir / ("ind0-" + name);
    const std::ifstream in(file);
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
    const auto output = trainfftest::runTrainFF(&fileManager(), "outliers", "lj-pc.xml", "booze-inter",
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
        const auto options = combine(mcmcOptions(3, seed), hfOptions());
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
    const auto freqSets = xvgSets(freq);
    EXPECT_EQ(12, freqSets.first);
    EXPECT_EQ(std::vector<std::string>({ "Train" }), freqSets.second);
    const auto forceSets = xvgSets(force);
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
    const auto output = trainfftest::runTrainFF(&fileManager(), "minimized", "halgren-PG.xml", "booze-mono",
                                                { { "-nooptimize", "" }, { "-calc_frequencies", "" },
                                            { "-fit", "kb" }, { "-fc_epot", "1" } });
    ASSERT_EQ(0, output.status);
    // The structures are written in the working directory of train_ff
    const auto dir = std::filesystem::path(fileManager().getTemporaryFilePath("minimized.log")).parent_path() / "inds";
    const auto pdb = dir / "methanol-original-minimized.pdb";
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
    const auto hasLine = [&output](const std::string &text)
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
    const auto output = trainfftest::runTrainFF(&fileManager(), "diatomics", "4S-lj14-7-pg-hvdw.xml", "hf-hf",
                                                { { "-nooptimize", "" }, { "-diatomics", "" },
                                            { "-fit", "sigma epsilon" }, { "-fc_epot", "1" } });
    ASSERT_EQ(0, output.status);
    // The file is named after the compound, in the working directory
    const auto fileName = (std::filesystem::path(fileManager().getTemporaryFilePath("x")).parent_path() /
                           "hydrogen-fluoride.xvg").string();
    const auto lines    = trainfftest::readLines(fileName);
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
    const auto selFile = fileManager().getTemporaryFilePath("acm.dat");
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
        options.emplace_back(f);
    }
    const auto output = trainfftest::runTrainFF(&fileManager(), "acm", "ACS-pg1.xml", "acm", options,
                                                fileManager().getInputFilePath("molprop/molprop.xml"),
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
    const auto dir    = std::filesystem::path(files["-mpout"]).parent_path();
    int  nxyz   = 0;
    for (const auto &entry : std::filesystem::directory_iterator(dir))
    {
        const auto name = entry.path().filename().string();
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

// Test: Without a sensitivity analysis there is no file and no summary
TEST_F(TrainFFTest, SensitivityNo)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "sensno", analysedForceField,
                                          analysedCompounds, analyseOnly("No"));
    ASSERT_EQ(0, output.status);
    EXPECT_TRUE(output.sensitivity.is_null());
    EXPECT_FALSE(logContains(output.log, "Sensitivity analysis ("));
    EXPECT_FALSE(logContains(output.log, "Details of the sensitivity analysis"));
    EXPECT_FALSE(output.json["train_ff"].contains("sensitivity"));
}

// Test: The one-dimensional analysis is written to its own file, the
// log file and train_ff.json get a summary
TEST_F(TrainFFTest, Sensitivity1DJson)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "sens1d", analysedForceField,
                                          analysedCompounds, analyseOnly("1D"));
    ASSERT_EQ(0, output.status);
    ASSERT_TRUE(output.sensitivity.contains("sensitivity"));
    const auto &sens = output.sensitivity["sensitivity"];
    EXPECT_EQ("1D", sens["mode"].get<std::string>());
    EXPECT_EQ("Train", sens["dataset"].get<std::string>());
    auto np = sens["n_parameters"].get<size_t>();
    EXPECT_GT(np, 0u);
    EXPECT_EQ(2*np, sens["n_evaluations"].get<size_t>());
    EXPECT_FALSE(sens.contains("pairs"));
    ASSERT_EQ(np, sens["parameters"].size());
    size_t nconv = 0;
    for (const auto &p : sens["parameters"])
    {
        const auto name = p["name"].get<std::string>();
        EXPECT_EQ(3u, p["points"].size()) << name;
        EXPECT_GE(p["p_orig"].get<double>(), p["lower"].get<double>()) << name;
        EXPECT_LE(p["p_orig"].get<double>(), p["upper"].get<double>()) << name;
        EXPECT_GT(p["chi2_orig"].get<double>(), 0) << name;
        // The force field is trained, all parameters change the chi2
        EXPECT_NE(0, p.value("a", 0.0)) << name;
        if (p["converged"].get<bool>())
        {
            nconv += 1;
        }
    }
    EXPECT_EQ(nconv, sens["n_parameters_converged"].get<size_t>());
    EXPECT_EQ(nconv == np, sens["minimum"].get<bool>());
    // The log file has the summary, not the parabola of each parameter
    EXPECT_TRUE(logContains(output.log, "Sensitivity analysis (1D) of " + std::to_string(np) + " parameters"));
    EXPECT_TRUE(logContains(output.log, "Converged: " + std::to_string(nconv) + "/" + std::to_string(np)));
    EXPECT_FALSE(logContains(output.log, "Fit to parabola"));
    // The summary in train_ff.json
    const auto &summary = output.json["train_ff"]["sensitivity"];
    EXPECT_EQ("1D", summary["mode"].get<std::string>());
    EXPECT_NE(std::string::npos, summary["file"].get<std::string>().find("sens1d-sensitivity.json"));
    EXPECT_EQ(np, trainfftest::jsonNumber(summary["n_parameters"]));
    EXPECT_EQ(nconv, trainfftest::jsonNumber(summary["n_parameters_converged"]));
    EXPECT_EQ(nconv == np ? "yes" : "no", summary["minimum"].get<std::string>());
}

// Test: The two-dimensional analysis has every pair of parameters, the log
// file and train_ff.json the numbers of pairs. The reference data are the
// same as those of TrainFFMpiTest.Sensitivity2DMiddlemen.
TEST_F(TrainFFTest, Sensitivity2DJson)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "sens2d", analysedForceField,
                                          analysedCompounds, analyseOnly("2D"));
    ASSERT_EQ(0, output.status);
    gmx::test::TestReferenceChecker checker(rootChecker());
    trainfftest::checkSensitivity2D(&checker, output);
    const auto &sens   = output.sensitivity["sensitivity"];
    auto        np     = sens["n_parameters"].get<size_t>();
    auto        npairs = sens["n_pairs"].get<size_t>();
    EXPECT_TRUE(logContains(output.log, "Sensitivity analysis (2D) of " + std::to_string(np) +
                            " parameters and " + std::to_string(npairs) + " pairs"));
    EXPECT_TRUE(logContains(output.log, "The Hessian of all parameters is"));
    const auto &summary = output.json["train_ff"]["sensitivity"];
    EXPECT_EQ("2D", summary["mode"].get<std::string>());
    EXPECT_EQ(npairs, trainfftest::jsonNumber(summary["n_pairs"]));
    EXPECT_EQ(sens["n_pairs_converged"].get<size_t>(), trainfftest::jsonNumber(summary["n_pairs_converged"]));
}

// Test: In the force field of hydrogen fluoride that training starts from,
// the Lennard-Jones parameters are zero and the virtual sites have no
// charge, such that the chi2 does not depend on them. The analysis reports
// them as flat, their pairs as degenerate, and says so in the log file.
TEST_F(TrainFFTest, SensitivityFlatParameters)
{
    auto options = setOption(combine(mcmcOptions(1), hfOptions()), "-sensitivity", "2D");
    options.emplace_back("-nooptimize", "");
    options.emplace_back("-norandom_init", "");
    auto output = trainfftest::runTrainFF(&fileManager(), "sensflat", "4S-lj14-7-pg-hvdw.xml",
                                          "hf-hf", options);
    ASSERT_EQ(0, output.status);
    ASSERT_TRUE(output.sensitivity.contains("sensitivity"));
    const auto &sens = output.sensitivity["sensitivity"];
    const std::set<std::string> expected = {
        "f delta", "f epsilon", "f gamma", "f sigma",
        "hf delta", "hf epsilon", "hf gamma", "hf sigma",
        "f_b~hf_b!v2f1 vs2a", "f_b~hf_b!v2f2 vs2a"
    };
    std::set<std::string> flat;
    for (const auto &p : sens["parameters"])
    {
        if (p["flat"].get<bool>())
        {
            flat.insert(p["name"].get<std::string>());
            EXPECT_FALSE(p["converged"].get<bool>()) << p["name"].get<std::string>();
        }
    }
    EXPECT_EQ(expected, flat);
    auto np = sens["n_parameters"].get<size_t>();
    EXPECT_EQ(expected.size(), sens["n_parameters_flat"].get<size_t>());
    EXPECT_FALSE(sens["minimum"].get<bool>());
    // Every pair with a flat parameter is degenerate
    for (const auto &pair : sens["pairs"])
    {
        if (flat.count(pair["names"][0].get<std::string>()) > 0 ||
            flat.count(pair["names"][1].get<std::string>()) > 0)
        {
            EXPECT_EQ("degenerate", pair["shape"].get<std::string>()) << pair["names"].dump();
            EXPECT_FALSE(pair["converged"].get<bool>()) << pair["names"].dump();
        }
    }
    EXPECT_TRUE(logContains(output.log, std::to_string(expected.size()) + " of the " + std::to_string(np) +
                            " parameters do not change the chi2 (flat)"));
    EXPECT_TRUE(logContains(output.log, "Their pairs are degenerate."));
    EXPECT_EQ(expected.size(), trainfftest::jsonNumber(output.json["train_ff"]["sensitivity"]["n_parameters_flat"]));
}

} // namespace

} // namespace alexandria
