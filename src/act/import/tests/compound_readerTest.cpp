/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2025,2026
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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <cstdlib>
#include <gtest/gtest.h>
#include <string>
#include <vector>

#include "act/basics/chargemodel.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/import/compound_reader.h"
#include "act/testutils/forcefield_utils.h"
#include "act/topology/actmol.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

class CompoundReaderTest : public gmx::test::CommandLineTestBase
{
protected:
    MsgHandler                       msghandler_;
    CompoundReader                   compR_;
    ForceComputer                    forceComp_;
    ForceField                      *pd;
    std::string                      qopt_;
    std::string                      fopt_;
    std::string                      charges_;
    gmx::test::TestReferenceChecker  checker_;

    //! Constructor that does initialization
    CompoundReaderTest() : checker_(this->rootChecker())
    {
        const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
        checker_.setDefaultTolerance(tolerance);
        pd = getForceField("ff-hh");
        qopt_.assign("-charges");
        fopt_.assign("-f");
        charges_ = gmx::test::TestFileManager::getInputFilePath("molprop/mp2-hh.xml");
    }

    //! Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    void hack_fnm(std::vector<t_filenm> *filenm)
    {
        for(auto &fn : *filenm)
        {
            if (qopt_.compare(fn.opt) == 0)
            {
                fn.fn    = charges_.c_str();
                fn.filenames.push_back(charges_);
                fn.flag &= ffSET;
            }
        }
    }
    void hack_pargs(std::vector<t_pargs> *pargs, std::string &infile)
    {
        for(auto &pa : *pargs)
        {
            if (fopt_.compare(pa.option) == 0 && !infile.empty())
            {
                *pa.u.c = infile.c_str();
                pa.bSet = true;
            }
        }
    }
    //! Do the actual testing
    void testCompoundReader ()
    {
        std::vector<ACTMol>       mols;
        std::vector<t_filenm>     filenm;
        std::vector<t_pargs>      pargs;
        std::vector<const char *> desc;
        CompoundReader            compR_;

        compR_.addOptions(&pargs, &filenm, &desc);
        hack_fnm(&filenm);
        compR_.optionsFinished(&msghandler_, filenm);
        checker_.checkString(compR_.qread(), "Charge type to read");
        checker_.checkString(chargeGenerationAlgorithmName(compR_.algorithm()).c_str(), "Algorithm");
        compR_.read(&msghandler_, *pd, &forceComp_, &mols);
        auto &qMap = compR_.chargeMapConst();
        checker_.checkInt64(qMap.size(), "Number of molecules");
        int i = 1;
        for(const auto &mol : qMap)
        {
            for(const auto &atom : mol.second)
            {
                const std::string aname = gmx::formatString("%s-%d",
                                                            atom.first.id().c_str(), i);
                checker_.checkDouble(atom.second, aname.c_str());
                i += 1;
            }
        }
    }

    //! Do the actual testing
    void testCompoundReaderFile (const std::string &input, bool hackFNM)
    {
        std::vector<ACTMol>       mols;
        std::vector<t_filenm>     filenm;
        std::vector<t_pargs>      pargs;
        std::vector<const char *> desc;
        CompoundReader            compR_;

        std::string dirinput("molecules/");
        dirinput += input;
        std::string infile = gmx::test::TestFileManager::getInputFilePath(dirinput);
        compR_.addOptions(&pargs, &filenm, &desc);
        hack_pargs(&pargs, infile);
        if (hackFNM)
        {
            hack_fnm(&filenm);
        }
        compR_.optionsFinished(&msghandler_, filenm);
        checker_.checkString(compR_.qread(), "Charge type to read");
        checker_.checkString(chargeGenerationAlgorithmName(compR_.algorithm()).c_str(), "Algorithm");
        compR_.read(&msghandler_, *pd, &forceComp_, &mols);
        checker_.checkInt64(mols.size(), "Number of molecules");
        for(const auto &mol : mols)
        {
            const auto &atoms = mol.atomsConst();
            int i = 1;
            for(const auto &atom : atoms)
            {
                const std::string aname = gmx::formatString("%s-%s-%d", mol.formula().c_str(), atom.name().c_str(), i);
                checker_.checkDouble(atom.charge(), aname.c_str());
                i += 1;
            }
        }
    }

    //! Clean the test data.
    static void TearDownTestCase()
    {
    }

};

TEST_F (CompoundReaderTest, Count){
    testCompoundReader();
}

TEST_F (CompoundReaderTest, HCl){
    testCompoundReaderFile("hydrogen-chloride.sdf", false);
}

TEST_F (CompoundReaderTest, Dimer){
    testCompoundReaderFile("hydrogen-bromide#hydrogen-fluoride.pdb", false);
}

TEST_F (CompoundReaderTest, Fluorane){
    testCompoundReaderFile("fluorane.sdf", false);
}

TEST_F (CompoundReaderTest, HCl_Charges){
    testCompoundReaderFile("hydrogen-chloride.sdf", true);
}

TEST_F (CompoundReaderTest, Dimer_Charges){
    testCompoundReaderFile("hydrogen-bromide#hydrogen-fluoride.pdb", true);
}

TEST_F (CompoundReaderTest, Fluorane_Charges){
    testCompoundReaderFile("fluorane.sdf", true);
}

// ---- The options of the compound reader from the command line ----
//
// What the command lines of gentop, minimize, nma and simulate give, such
// that a change in how the options are turned into settings shows.

/*! \brief Parse a command line into a compound reader as the tools do
 * \param[in] compR      The reader
 * \param[in] args       The arguments
 * \param[in] msghandler Gets what is wrong with the options
 * \param[in] cmd        The command line, which the reader refers to and
 *                       which must exist as long as the reader is used
 */
static void parseCompoundOptions(CompoundReader                  *compR,
                                 const std::vector<std::string>  &args,
                                 MsgHandler                      *msghandler,
                                 gmx::test::CommandLine          *cmd)
{
    std::vector<t_pargs>      pargs;
    std::vector<t_filenm>     fnms;
    std::vector<const char *> desc;
    compR->addOptions(&pargs, &fnms, &desc);
    cmd->append("test");
    for (const auto &a : args)
    {
        cmd->append(a);
    }
    gmx_output_env_t *oenv = nullptr;
    ASSERT_TRUE(parse_common_args(&cmd->argc(), cmd->argv(), 0, fnms.size(), fnms.data(),
                                  pargs.size(), pargs.data(), 0, nullptr, 0, nullptr, &oenv));
    output_env_done(oenv);
    compR->optionsFinished(msghandler, fnms);
}

TEST(CompoundReaderOptionsTest, Defaults)
{
    CompoundReader         compR;
    MsgHandler             msghandler;
    gmx::test::CommandLine cmd;
    parseCompoundOptions(&compR, {}, &msghandler, &cmd);
    EXPECT_TRUE(msghandler.ok());
    EXPECT_EQ(ChargeGenerationAlgorithm::SQE, compR.algorithm());
    EXPECT_STREQ("", compR.qread());
    EXPECT_STREQ("", compR.customCharges());
    EXPECT_FALSE(compR.userQtot());
    EXPECT_STREQ("", compR.fileName());
    EXPECT_STREQ("", compR.smiles());
    EXPECT_STREQ("", compR.dbName());
    EXPECT_TRUE(compR.chargeFile().empty());
    EXPECT_FALSE(compR.oneH());
    EXPECT_EQ(0, compR.totalCharge());
    EXPECT_EQ(0U, compR.molselect().nMol());
}

TEST(CompoundReaderOptionsTest, Input)
{
    const auto molprop = gmx::test::TestFileManager::getInputFilePath("molprop/molprop.xml");
    {
        CompoundReader         compR;
        MsgHandler             msghandler;
        gmx::test::CommandLine cmd;
        parseCompoundOptions(&compR, { "-f", "acetone.sdf", "-qtot", "-1", "-oneH" }, &msghandler, &cmd);
        EXPECT_TRUE(msghandler.ok());
        EXPECT_STREQ("acetone.sdf", compR.fileName());
        EXPECT_DOUBLE_EQ(-1, compR.totalCharge());
        EXPECT_TRUE(compR.oneH());
    }
    {
        CompoundReader         compR;
        MsgHandler             msghandler;
        gmx::test::CommandLine cmd;
        parseCompoundOptions(&compR, { "-smiles", "CC(=O)C" }, &msghandler, &cmd);
        EXPECT_TRUE(msghandler.ok());
        EXPECT_STREQ("CC(=O)C", compR.smiles());
        EXPECT_STREQ("", compR.fileName());
    }
    {
        CompoundReader         compR;
        MsgHandler             msghandler;
        gmx::test::CommandLine cmd;
        parseCompoundOptions(&compR, { "-db", "water", "-charges", molprop }, &msghandler, &cmd);
        EXPECT_TRUE(msghandler.ok());
        EXPECT_STREQ("water", compR.dbName());
        EXPECT_EQ(molprop, compR.chargeFile());
    }
    {
        // A selection of compounds from the charge file
        CompoundReader         compR;
        MsgHandler             msghandler;
        gmx::test::CommandLine cmd;
        parseCompoundOptions(&compR, { "-sel", gmx::test::TestFileManager::getInputFilePath("molprop/booze-mono.dat"),
                                       "-charges", molprop }, &msghandler, &cmd);
        EXPECT_TRUE(msghandler.ok());
        EXPECT_EQ(2U, compR.molselect().nMol());
    }
}

TEST(CompoundReaderOptionsTest, ChargeAlgorithm)
{
    for (auto alg : { ChargeGenerationAlgorithm::NONE, ChargeGenerationAlgorithm::EEM, ChargeGenerationAlgorithm::SQE,
                      ChargeGenerationAlgorithm::ESP, ChargeGenerationAlgorithm::Custom,
                      ChargeGenerationAlgorithm::Read })
    {
        const auto             name = chargeGenerationAlgorithmName(alg);
        CompoundReader         compR;
        MsgHandler             msghandler;
        gmx::test::CommandLine cmd;
        parseCompoundOptions(&compR, { "-qalg", name }, &msghandler, &cmd);
        EXPECT_TRUE(msghandler.ok()) << name;
        EXPECT_EQ(alg, compR.algorithm()) << name;
    }
}

TEST(CompoundReaderOptionsTest, ChargesWithoutAlgorithm)
{
    // Charges to read or custom charges choose the algorithm
    {
        CompoundReader         compR;
        MsgHandler             msghandler;
        gmx::test::CommandLine cmd;
        parseCompoundOptions(&compR, { "-qqm", "qMulliken" }, &msghandler, &cmd);
        EXPECT_TRUE(msghandler.ok());
        EXPECT_EQ(ChargeGenerationAlgorithm::Read, compR.algorithm());
        EXPECT_STREQ("qMulliken", compR.qread());
    }
    {
        CompoundReader         compR;
        MsgHandler             msghandler;
        gmx::test::CommandLine cmd;
        // A leading space, such that the value is not taken for an option
        parseCompoundOptions(&compR, { "-qcustom", " -0.8 0.4 0.4" }, &msghandler, &cmd);
        EXPECT_TRUE(msghandler.ok());
        EXPECT_EQ(ChargeGenerationAlgorithm::Custom, compR.algorithm());
        EXPECT_STREQ(" -0.8 0.4 0.4", compR.customCharges());
        EXPECT_TRUE(compR.userQtot());
    }
    {
        // An algorithm that is given goes first, with a warning
        CompoundReader         compR;
        MsgHandler             msghandler;
        gmx::test::CommandLine cmd;
        parseCompoundOptions(&compR, { "-qalg", "EEM", "-qqm", "qMulliken" }, &msghandler, &cmd);
        EXPECT_TRUE(msghandler.ok());
        EXPECT_EQ(ChargeGenerationAlgorithm::EEM, compR.algorithm());
    }
}

TEST(CompoundReaderOptionsTest, ConflictingOptions)
{
    const auto molprop   = gmx::test::TestFileManager::getInputFilePath("molprop/molprop.xml");
    const auto selection = gmx::test::TestFileManager::getInputFilePath("molprop/booze-mono.dat");
    const std::vector<std::vector<std::string>> wrong = {
        { "-f", "acetone.sdf", "-smiles", "CC(=O)C" },
        { "-f", "acetone.sdf", "-db", "water", "-charges", molprop },
        { "-smiles", "O", "-db", "water", "-charges", molprop },
        { "-sel", selection, "-f", "acetone.sdf" },
        { "-sel", selection, "-smiles", "O" },
        { "-db", "water" },
        { "-qcustom", " -0.8 0.4 0.4", "-qqm", "qMulliken" }
    };
    for (const auto &args : wrong)
    {
        CompoundReader         compR;
        MsgHandler             msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        gmx::test::CommandLine cmd;
        parseCompoundOptions(&compR, args, &msghandler, &cmd);
        std::string all;
        for (const auto &a : args)
        {
            all += a + " ";
        }
        EXPECT_FALSE(msghandler.ok()) << all;
    }
}

TEST(CompoundReaderOptionsTest, CopyOutlivesTheCommandLine)
{
    // The settings of a reader belong to it, not to the command line or to
    // the reader it was copied from
    CompoundReader copy;
    {
        CompoundReader         compR;
        MsgHandler             msghandler;
        gmx::test::CommandLine cmd;
        parseCompoundOptions(&compR, { "-f", "acetone.sdf", "-qqm", "qMulliken" }, &msghandler, &cmd);
        copy = compR;
    }
    EXPECT_STREQ("acetone.sdf", copy.fileName());
    EXPECT_STREQ("qMulliken", copy.qread());
    CompoundReader second;
    {
        CompoundReader set;
        set.setSmiles("CC(=O)C");
        second = set;
        set.setFileName("other.sdf");
    }
    EXPECT_STREQ("CC(=O)C", second.smiles());
    EXPECT_STREQ("", second.fileName());
}

TEST(CompoundReaderOptionsTest, UnknownChargeAlgorithm)
{
    CompoundReader         compR;
    MsgHandler             msghandler;
    gmx::test::CommandLine cmd;
    EXPECT_THROW(parseCompoundOptions(&compR, { "-qalg", "Magic" }, &msghandler, &cmd), gmx::InvalidInputError);
}

} // namespace alexandria
