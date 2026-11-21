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
 * Tests of the edit_mp module.
 */

#include "actpre.h"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <map>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/basics/msg_handler.h"
#include "act/molprop/composition.h"
#include "act/molprop/edit_mp.h"
#include "act/molprop/experiment.h"
#include "act/molprop/fragment.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "act/molprop/molpropobservable.h"
#include "act/testutils/output_checks.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class EditMPTest : public gmx::test::CommandLineTestBase
{
protected:
    //! The log file of the last run
    std::string log_;

    /*! \brief Run edit_mp
     * \param[in] inputs  The molprop files
     * \param[in] options The other options
     * \return the molecules in the output file
     */
    std::vector<MolProp> run(const std::vector<std::string>                        &inputs,
                             const std::vector<std::pair<std::string, std::string>> &options = {})
    {
        const auto output = fileManager().getTemporaryFilePath("out.xml");
        log_        = fileManager().getTemporaryFilePath("edit_mp.log");
        gmx::test::CommandLine cmd;
        cmd.append("edit_mp");
        cmd.addOption("-mp");
        for (const auto &in : inputs)
        {
            cmd.append(in);
        }
        cmd.addOption("-o", output);
        cmd.addOption("-g", log_);
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
        {
            // Some output is written in the working directory
            const auto old = std::filesystem::current_path();
            std::filesystem::current_path(std::filesystem::path(output).parent_path());
            EXPECT_EQ(0, edit_mp(cmd.argc(), cmd.argv()));
            std::filesystem::current_path(old);
        }
        std::vector<MolProp> mps;
        EXPECT_TRUE(std::filesystem::exists(output));
        MolPropRead(nullptr, output.c_str(), &mps);
        return mps;
    }

    //! \return the molecules in a file
    static std::vector<MolProp> read(const std::string &fileName)
    {
        std::vector<MolProp> mps;
        MolPropRead(nullptr, fileName.c_str(), &mps);
        return mps;
    }

    //! \return the number of experiments of each molecule
    static std::map<std::string, size_t> experiments(const std::vector<MolProp> &mps)
    {
        std::map<std::string, size_t> n;
        for (const auto &mp : mps)
        {
            n[mp.getMolname()] += mp.experimentConst().size();
        }
        return n;
    }
};

// The results of the tool for several sets of options, as reference data,
// such that a change in how the options reach the work shows

TEST_F(EditMPTest, ResultsRewrite)
{
    // The tool adds what the file of the test data lacks
    run({ fileManager().getInputFilePath("molprop/hf-hf.xml") });
    gmx::test::TestReferenceChecker checker(this->rootChecker());
    checker.setDefaultTolerance(buildTolerance(1e-6, 1e-8, 1e-6));
    checkTextFile(&checker, fileManager().getTemporaryFilePath("out.xml"), {}, "molprop");
}

TEST_F(EditMPTest, ResultsMerge)
{
    const auto json = fileManager().getTemporaryFilePath("merge.json");
    run({ fileManager().getInputFilePath("molprop/booze-mono.xml"),
          fileManager().getInputFilePath("molprop/water-butanol.xml") },
        { { "-json", json } });
    const auto tree = readJsonOutput(json);
    ASSERT_TRUE(tree.is_object() && tree.contains("edit_mp")) << json;
    gmx::test::TestReferenceChecker checker(this->rootChecker());
    checker.setDefaultTolerance(buildTolerance(1e-6, 1e-8, 1e-6));
    checkJson(&checker, tree["edit_mp"], "edit_mp");
}

TEST_F(EditMPTest, ResultsHistogram)
{
    // Other limits for the warnings and another temperature
    run({ fileManager().getInputFilePath("molprop/hf-hf.xml") },
        { { "-ehisto", "" }, { "-ewarnLow", "-15" }, { "-ewarnHi", "5" }, { "-temperature", "300" } });
    const auto histo = std::filesystem::path(fileManager().getTemporaryFilePath("out.xml")).parent_path() /
        "hydrogen-fluoride#hydrogen-fluoride-InteractionEnergy.xvg";
    gmx::test::TestReferenceChecker checker(this->rootChecker());
    checker.setDefaultTolerance(buildTolerance(1e-6, 1e-8, 1e-6));
    checkTextFile(&checker, histo.string(), { "#" }, "histogram");
    checkTextFile(&checker, fileManager().getTemporaryFilePath("out.xml"), {}, "molprop");
}

//! \return the contents of a file
static std::string fileContents(const std::string &fileName)
{
    const std::ifstream in(fileName);
    std::stringstream text;
    text << in.rdbuf();
    return text.str();
}

// The molprop files without a command line, as a program that uses ACT edits them

TEST_F(EditMPTest, WithoutCommandLineSameAsTheProgram)
{
    // A merge with validation and the description in JSON
    const auto json = fileManager().getTemporaryFilePath("tool.json");
    const std::vector<std::string> inputs = { fileManager().getInputFilePath("molprop/booze-mono.xml"),
                                              fileManager().getInputFilePath("molprop/water-butanol.xml") };
    run(inputs, { { "-json", json }, { "-validate", "" }, { "-nomolinfo", "" } });
    EditMPOptions options;
    options.inputFiles = inputs;
    options.outputFile = fileManager().getTemporaryFilePath("direct.xml");
    options.jsonFile   = fileManager().getTemporaryFilePath("direct.json");
    options.validate   = true;
    options.molinfo    = false;
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    runEditMP(&msghandler, options);
    EXPECT_EQ(fileContents(fileManager().getTemporaryFilePath("out.xml")), fileContents(options.outputFile));
    EXPECT_EQ(fileContents(json), fileContents(options.jsonFile));
}

TEST_F(EditMPTest, WithoutCommandLineWrongOptions)
{
    const auto valid = [this]()
    {
        EditMPOptions options;
        options.inputFiles = { fileManager().getInputFilePath("molprop/hf-hf.xml") };
        options.outputFile = fileManager().getTemporaryFilePath("wrong.xml");
        return options;
    };
    std::vector<std::pair<EditMPOptions, std::string>> wrong;
    {
        auto o = valid();
        o.inputFiles.clear();
        wrong.emplace_back(o, "-mp");
    }
    {
        auto o = valid();
        o.outputFile.clear();
        wrong.emplace_back(o, "-o");
    }
    {
        auto o = valid();
        o.temperature = 0;
        wrong.emplace_back(o, "-temperature");
    }
    {
        auto o = valid();
        o.energyWarnHigh = std::nan("");
        wrong.emplace_back(o, "-ewarnHi");
    }
    for (const auto &w : wrong)
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        try
        {
            runEditMP(&msghandler, w.first);
            ADD_FAILURE() << "Accepted options with a wrong " << w.second;
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find(w.second)) << ex.what();
        }
        EXPECT_FALSE(!w.first.outputFile.empty() && std::filesystem::exists(w.first.outputFile)) << w.second;
    }
}

// Test: Writing the molecules of a file gives the same molecules
TEST_F(EditMPTest, RoundTrip)
{
    const auto input = fileManager().getInputFilePath("molprop/molprop.xml");
    auto ref   = read(input);
    auto mps   = run({ input });
    ASSERT_EQ(ref.size(), mps.size());
    EXPECT_EQ(experiments(ref), experiments(mps));
    for (size_t i = 0; i < ref.size(); i++)
    {
        EXPECT_EQ(ref[i].getMolname(), mps[i].getMolname());
        EXPECT_EQ(ref[i].NAtom(), mps[i].NAtom());
        EXPECT_EQ(ref[i].formula(), mps[i].formula());
    }
}

// Test: Molecules of different files are combined, sorted by name
TEST_F(EditMPTest, MergeFiles)
{
    const auto first  = fileManager().getInputFilePath("molprop/molprop.xml");
    const auto second = fileManager().getInputFilePath("molprop/booze-mono.xml");
    const auto mps    = run({ first, second });
    std::vector<std::string> names;
    for (const auto &mp : mps)
    {
        names.push_back(mp.getMolname());
    }
    const std::vector<std::string> expected = { "1-butanol", "acetone", "ethanol", "methanol", "water" };
    EXPECT_EQ(expected, names);
}

// Test: The same molecule in two files is one molecule with the
// experiments of both
TEST_F(EditMPTest, MergeDuplicates)
{
    const auto input = fileManager().getInputFilePath("molprop/molprop.xml");
    const auto ref   = experiments(read(input));
    const auto mps   = run({ input, input });
    auto n     = experiments(mps);
    ASSERT_EQ(ref.size(), mps.size());
    for (const auto &r : ref)
    {
        EXPECT_EQ(2*r.second, n[r.first]) << r.first;
    }
}

// Test: With a force field the compounds are checked, the results are in
// the log file, and the output has all their data
TEST_F(EditMPTest, CheckWithForceField)
{
    const auto input = fileManager().getInputFilePath("molprop/molprop.xml");
    const auto mps   = run({ input }, { { "-ff", fileManager().getInputFilePath("forcefields/ACS-g1.xml") } });
    // Checking does not change the compounds
    EXPECT_EQ(read(input).size(), mps.size());
    EXPECT_EQ(experiments(read(input)), experiments(mps));
    const std::ifstream in(log_);
    std::stringstream text;
    text << in.rdbuf();
    EXPECT_NE(std::string::npos, text.str().find("Force field file")) << text.str();
}

// Test: A monomer and its dimer, which has no dipole or electrostatic potential of an
// optimized structure, are checked with a force field without bond charge corrections
TEST_F(EditMPTest, CheckDimersWithForceField)
{
    const auto input = fileManager().getInputFilePath("molprop/hf-hf.xml");
    const auto mps   = run({ input }, { { "-ff", fileManager().getInputFilePath("forcefields/HF-LJ14-7.xml") } });
    EXPECT_EQ(read(input).size(), mps.size());
    EXPECT_EQ(experiments(read(input)), experiments(mps));
    const std::ifstream in(log_);
    std::stringstream text;
    text << in.rdbuf();
    EXPECT_NE(std::string::npos, text.str().find("Succeed making 2 topologies, failed for 0 compounds")) << text.str();
}

//! \return the text of a file
static std::string fileText(const std::string &fileName)
{
    const std::ifstream in(fileName);
    std::stringstream   text;
    text << in.rdbuf();
    return text.str();
}

// Test: -validate finds no problems in a file with consistent dimers, and
// writes all molecules
TEST_F(EditMPTest, ValidateConsistent)
{
    const auto input = fileManager().getInputFilePath("molprop/ion-water.xml");
    const auto mps   = run({ input }, { { "-validate", "" } });
    EXPECT_EQ(read(input).size(), mps.size());
    EXPECT_EQ(experiments(read(input)), experiments(mps));
    const auto text = fileText(log_);
    EXPECT_EQ(std::string::npos, text.find("occurs as")) << text;
    EXPECT_EQ(std::string::npos, text.find("has its atoms in the order")) << text;
}

// Test: -validate accepts a file with a dimer as water#chloride and as
// chloride#water, if each is consistent
TEST_F(EditMPTest, ValidateDimerInBothOrders)
{
    auto mps   = read(fileManager().getInputFilePath("molprop/ion-water.xml"));
    auto again = read(fileManager().getInputFilePath("molprop/ion-water.xml"));
    for (auto &mp : again)
    {
        if (mp.getMolname() == "water#chloride")
        {
            mp.SetMolname("chloride#water");
            mps.push_back(std::move(mp));
            break;
        }
    }
    const auto input = fileManager().getTemporaryFilePath("both-orders.xml");
    MolPropWrite(input.c_str(), mps, false);
    run({ input }, { { "-validate", "" } });
    const auto text = fileText(log_);
    EXPECT_EQ(std::string::npos, text.find("Error")) << text;
}

//! \return the content of a JSON file, which has to exist
static nlohmann::json readJson(const std::string &fileName)
{
    std::ifstream in(fileName);
    EXPECT_TRUE(in.good()) << fileName;
    return in.good() ? nlohmann::json::parse(in) : nlohmann::json();
}

// Test: -json describes the molecules that are written: their fragments,
// the number of experiments and the properties, and the totals
TEST_F(EditMPTest, Json)
{
    const auto input    = fileManager().getInputFilePath("molprop/ion-water.xml");
    const auto jsonFile = fileManager().getTemporaryFilePath("mp.json");
    std::filesystem::remove(jsonFile);
    const auto mps = run({ input }, { { "-json", jsonFile } });
    const auto j   = readJson(jsonFile);
    ASSERT_TRUE(j.contains("edit_mp"));
    const auto &molecules = j.at("edit_mp").at("molecules");
    ASSERT_EQ(mps.size(), molecules.size());
    size_t nexp = 0, ndimer = 0, nprop = 0;
    for (size_t i = 0; i < mps.size(); i++)
    {
        const auto &mp = mps[i];
        const auto &jm = molecules[i];
        SCOPED_TRACE(mp.getMolname());
        EXPECT_EQ(mp.getMolname(), jm.at("name").get<std::string>());
        EXPECT_EQ(mp.formula(), jm.at("formula").get<std::string>());
        EXPECT_EQ(mp.experimentConst().size(), jm.at("experiments").get<size_t>());
        ASSERT_EQ(mp.fragments().size(), jm.at("fragments").size());
        int    charge = 0;
        size_t natom  = 0;
        for (size_t f = 0; f < mp.fragments().size(); f++)
        {
            const auto &frag = mp.fragments()[f];
            const auto &jf   = jm.at("fragments")[f];
            EXPECT_EQ(frag.inchi(), jf.at("inchi").get<std::string>());
            EXPECT_EQ(frag.formula(), jf.at("formula").get<std::string>());
            EXPECT_EQ(frag.charge(), jf.at("charge").get<int>());
            EXPECT_EQ(frag.multiplicity(), jf.at("multiplicity").get<int>());
            EXPECT_EQ(frag.atoms().size(), jf.at("natom").get<size_t>());
            charge += frag.charge();
            natom  += frag.atoms().size();
        }
        EXPECT_EQ(charge, jm.at("charge").get<int>());
        EXPECT_EQ(natom, jm.at("natom").get<size_t>());
        EXPECT_GT(natom, 0u);
        // The number of values of each property, over all experiments
        std::map<std::string, size_t> props;
        for (const auto &exp : mp.experimentConst())
        {
            for (const auto &prop : exp.propertiesConst())
            {
                props[mpo_name(prop.first)] += prop.second.size();
            }
        }
        EXPECT_EQ(props.size(), jm.at("properties").size());
        for (const auto &prop : props)
        {
            EXPECT_EQ(prop.second, jm.at("properties").at(prop.first).get<size_t>()) << prop.first;
            nprop += prop.second;
        }
        nexp   += mp.experimentConst().size();
        ndimer += mp.fragments().size() == 2 ? 1 : 0;
    }
    const auto &summary = j.at("edit_mp").at("summary");
    EXPECT_EQ(mps.size(), summary.at("molecules").get<size_t>());
    EXPECT_EQ(ndimer, summary.at("dimers").get<size_t>());
    EXPECT_EQ(nexp, summary.at("experiments").get<size_t>());
    size_t total = 0;
    for (const auto &prop : summary.at("properties").items())
    {
        total += prop.value().get<size_t>();
    }
    EXPECT_EQ(nprop, total);
    EXPECT_GT(nprop, 0u);
    EXPECT_GT(ndimer, 0u);
    // Not asked for
    EXPECT_FALSE(j.at("edit_mp").at("validation").at("checked").get<bool>());
}

// Test: The result of -validate is in the JSON output: no problems for a
// consistent file
TEST_F(EditMPTest, JsonValidationClean)
{
    const auto input    = fileManager().getInputFilePath("molprop/ion-water.xml");
    const auto jsonFile = fileManager().getTemporaryFilePath("clean.json");
    run({ input }, { { "-validate", "" }, { "-json", jsonFile } });
    const auto validation = readJson(jsonFile).at("edit_mp").at("validation");
    EXPECT_TRUE(validation.at("checked").get<bool>());
    EXPECT_TRUE(validation.at("ok").get<bool>());
    EXPECT_EQ(0, validation.at("problems").get<int>());
    EXPECT_TRUE(validation.at("invalid").is_array());
    EXPECT_TRUE(validation.at("invalid").empty());
    EXPECT_TRUE(validation.at("messages").is_array());
    EXPECT_TRUE(validation.at("messages").empty());
}

/*! \brief A molecule with one experiment, made of fragments
 * \param[in] name  The name
 * \param[in] frags For each fragment its InChI, its name, its formula and
 *                  the elements of its atoms
 * \return the molecule
 */
static MolProp makeMolecule(const std::string &name,
                            const std::vector<std::tuple<std::string, std::string, std::string, std::vector<std::string>>> &frags)
{
    MolProp    mp;
    mp.SetMolname(name);
    Experiment exp("ref", "min");
    int        natom = 0;
    for (const auto &f : frags)
    {
        std::vector<int> atoms;
        for (const auto &a : std::get<3>(f))
        {
            atoms.push_back(natom);
            natom += 1;
            CalcAtom ca(a, a, natom);
            ca.setCoords(0.1*natom, 0, 0);
            exp.AddAtom(std::move(ca));
        }
        Fragment frag(std::get<0>(f), 1.0, 0, 1, 1, std::get<2>(f), atoms);
        frag.setIupac(std::get<1>(f));
        mp.addFragment(frag);
    }
    mp.AddExperiment(std::move(exp));
    return mp;
}

// Test: A compound with its atoms in another order in one molecule is a
// problem, which the JSON output names with the molecule and the message
TEST_F(EditMPTest, JsonValidationProblem)
{
    const std::string              inchiWater = "InChI=1S/H2O/h1H2";
    const std::string              inchiMeoh  = "InChI=1S/CH4O/c1-2/h2H,1H3";
    const std::vector<std::string> ohh        = { "O", "H", "H" };
    const std::vector<std::string> hoh        = { "H", "O", "H" };
    const std::vector<std::string> meoh       = { "C", "O", "H", "H", "H", "H" };
    std::vector<MolProp>           mps;
    mps.push_back(makeMolecule("water#methanol", { { inchiWater, "water", "H2O", hoh }, { inchiMeoh, "methanol", "CH4O", meoh } }));
    mps.push_back(makeMolecule("water", { { inchiWater, "water", "H2O", ohh } }));
    mps.push_back(makeMolecule("water#water", { { inchiWater, "water", "H2O", ohh }, { inchiWater, "water", "H2O", ohh } }));
    const auto input = fileManager().getTemporaryFilePath("order.xml");
    MolPropWrite(input.c_str(), mps, false);
    const auto jsonFile = fileManager().getTemporaryFilePath("order.json");
    run({ input }, { { "-validate", "" }, { "-nomolinfo", "" }, { "-json", jsonFile } });
    const auto validation = readJson(jsonFile).at("edit_mp").at("validation");
    EXPECT_TRUE(validation.at("checked").get<bool>());
    EXPECT_FALSE(validation.at("ok").get<bool>());
    EXPECT_EQ(1, validation.at("problems").get<int>());
    ASSERT_EQ(1u, validation.at("invalid").size());
    EXPECT_EQ("water#methanol", validation.at("invalid")[0].get<std::string>());
    ASSERT_FALSE(validation.at("messages").empty());
    bool named = false;
    for (const auto &message : validation.at("messages"))
    {
        EXPECT_FALSE(message.at("level").get<std::string>().empty());
        EXPECT_TRUE(message.contains("kind"));
        named = named || message.at("text").get<std::string>().find("water#methanol") != std::string::npos;
    }
    EXPECT_TRUE(named);
}

// Test: Without the option no JSON file is written
TEST_F(EditMPTest, NoJsonByDefault)
{
    const auto input = fileManager().getInputFilePath("molprop/ion-water.xml");
    run({ input });
    // run() works in the directory of the output files
    const auto dir = std::filesystem::path(fileManager().getTemporaryFilePath("out.xml")).parent_path();
    EXPECT_FALSE(std::filesystem::exists(dir / "edit_mp.json"));
}

// Test: The histograms of the interaction energies are written for each
// compound that has them, the range is reported in the log file
TEST_F(EditMPTest, EnergyHistogram)
{
    const auto input = fileManager().getInputFilePath("molprop/hf-hf.xml");
    run({ input }, { { "-ehisto", "" } });
    const auto histo = std::filesystem::path(fileManager().getTemporaryFilePath("out.xml")).parent_path() /
        "hydrogen-fluoride#hydrogen-fluoride-InteractionEnergy.xvg";
    std::ifstream in(histo);
    ASSERT_TRUE(in.good()) << histo;
    int    npoints = 0;
    double total   = 0;
    for (std::string line; std::getline(in, line); )
    {
        if (line.empty() || line[0] == '#' || line[0] == '@')
        {
            continue;
        }
        std::istringstream iss(line);
        double x, y;
        ASSERT_TRUE(static_cast<bool>(iss >> x >> y)) << line;
        npoints++;
        total += y;
    }
    // Each energy is counted in one bin
    EXPECT_GT(npoints, 0);
    EXPECT_GT(total, 0);
    const std::ifstream logIn(log_);
    std::stringstream text;
    text << logIn.rdbuf();
    EXPECT_NE(std::string::npos, text.str().find("Range of energies for hydrogen-fluoride#hydrogen-fluoride")) << text.str();
    std::filesystem::remove(histo);
}

} // namespace

} // namespace alexandria
