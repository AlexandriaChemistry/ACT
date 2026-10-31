/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2026
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
#include <cmath>

#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <string>

#include <gtest/gtest.h>

#include "act/alexandria/openmm_xml.h"
#include "act/alexandria/pdbwriter.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/topology/actmol_util.h"

#include "testutils/cmdlinetest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace alexandria
{

namespace
{

class OpenMMXmlTest : public gmx::test::CommandLineTestBase
{
    protected:
    void testOpenMM(const std::string &forceField,
                    const std::string &fileName,
                    bool               numberAtypes,
                    double             mDrude = 0.1)
    {
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        auto ifm = gmx::test::TextFileMatch(gmx::test::ExactTextMatch()).createFileMatcher();
        gmx::test::TestFileManager tfm;
        MsgHandler                 msghandler; 
        msghandler.setPrintLevel(ACTStatus::Warning);
        auto ff         = getForceField(forceField);
        double rmsToler = 0.0000001;
        ForceComputer fcomp(rmsToler, 25, 1.0);
        std::vector<ACTMol> mps;
        std::string molfn("mols/");
        molfn += fileName;
        auto dataName = gmx::test::TestFileManager::getInputFilePath(molfn);
        initACTMol(&msghandler, dataName, ff, &fcomp, &mps);

        auto tmpFile  = tfm.getTemporaryFilePath("xml");
        int  ntrain   = 1;
        writeOpenMM(&msghandler, tmpFile, ff, mps, mDrude, ntrain, numberAtypes);
        EXPECT_TRUE(msghandler.ok());
        ifm->checkFile(tmpFile, &checker);
        // OpenMM finds the parent of a shell from its type, hence each
        // shell type can have only one parent type, see
        // https://github.com/AlexandriaChemistry/ACT/issues/747
        std::ifstream         xml(tmpFile);
        std::string           line;
        std::set<std::string> shellTypes;
        while (std::getline(xml, line))
        {
            auto pos = line.find("<Particle type1=\"");
            if (pos != std::string::npos)
            {
                auto begin = pos + std::strlen("<Particle type1=\"");
                auto type1 = line.substr(begin, line.find('"', begin) - begin);
                EXPECT_TRUE(shellTypes.insert(type1).second) << "Shell type " << type1 << " has more than one parent type";
            }
        }
    }
};

TEST_F(OpenMMXmlTest, Sulfate)
{
    testOpenMM("ACS-g", "sulfate.sdf", true);
}

TEST_F(OpenMMXmlTest, Furan)
{
    testOpenMM("ACS-g", "furan.sdf", true);
}

TEST_F(OpenMMXmlTest, Acetone)
{
    testOpenMM("ACS-g", "acetone.sdf", true);
}

TEST_F(OpenMMXmlTest, AcetonePolMDrude)
{
    testOpenMM("ACS-pg", "acetone.sdf", true, 0.2);
}

TEST_F(OpenMMXmlTest, SulfatePol)
{
    testOpenMM("ACS-pg", "sulfate.sdf", true);
}

TEST_F(OpenMMXmlTest, FuranPol)
{
    testOpenMM("ACS-pg", "furan.sdf", true);
}

TEST_F(OpenMMXmlTest, FuranPolNoNum)
{
    testOpenMM("ACS-pg", "furan.sdf", false);
}

TEST_F(OpenMMXmlTest, AcetonePol)
{
    testOpenMM("ACS-pg", "acetone.sdf", true);
}

TEST_F(OpenMMXmlTest, Merged)
{
    testOpenMM("ACS-g", "merged.pdb", true);
}

// Test: Water and methanol have the same shell type on their oxygens,
// which needs different OpenMM shell types for the two parents
TEST_F(OpenMMXmlTest, MethanolWaterPol)
{
    testOpenMM("ACS-pg", "methanol#water.pdb", true);
}

// Test: Without numbering the atom types this cannot be written
TEST(OpenMMXmlNoNumberTest, MethanolWaterPol)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    auto          ff = getForceField("ACS-pg");
    ForceComputer fcomp(1e-7, 25, 1.0);
    std::vector<ACTMol> mps;
    auto dataName = gmx::test::TestFileManager::getInputFilePath("mols/methanol#water.pdb");
    initACTMol(&msghandler, dataName, ff, &fcomp, &mps);
    ASSERT_TRUE(msghandler.ok());
    gmx::test::TestFileManager tfm;
    writeOpenMM(&msghandler, tfm.getTemporaryFilePath("xml"), ff, mps, 0.1, 1, false);
    EXPECT_FALSE(msghandler.ok());
}

TEST_F(OpenMMXmlTest, WaterHyperPol)
{
    testOpenMM("water_hyperpol", "water_dimer.pdb", true);
}

//! \return the residue names for methanol and water, with userNames
static std::map<std::string, std::string> methanolWaterNames(MsgHandler        *msghandler,
                                                             const std::string &userNames)
{
    auto          ff = getForceField("ACS-g");
    ForceComputer fcomp(1e-7, 25, 1.0);
    std::vector<ACTMol> mps;
    auto dataName = gmx::test::TestFileManager::getInputFilePath("mols/methanol#water.pdb");
    initACTMol(msghandler, dataName, ff, &fcomp, &mps);
    EXPECT_TRUE(msghandler->ok());
    auto names = openMMResidueNames(msghandler, mps, userNames);
    // The order of the compounds in the input
    std::map<std::string, std::string> byMolecule;
    const auto &ids = mps[0].fragmentHandler()->ids();
    if (ids.size() == 2 && names.size() == 2)
    {
        byMolecule.insert({ "methanol", names[ids[0]] });
        byMolecule.insert({ "water", names[ids[1]] });
    }
    return byMolecule;
}

TEST(OpenMMResidueNamesTest, Generated)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    auto names = methanolWaterNames(&msghandler, "");
    EXPECT_TRUE(msghandler.ok());
    // MET is methionine in OpenMM
    std::map<std::string, std::string> expected = { { "methanol", "ME1" }, { "water", "OXI" } };
    EXPECT_EQ(expected, names);
}

TEST(OpenMMResidueNamesTest, User)
{
    for (const auto &userNames : { "water=WTR methanol=MOH", "MOH,WTR", "water=WTR MOH" })
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        auto names = methanolWaterNames(&msghandler, userNames);
        EXPECT_TRUE(msghandler.ok()) << userNames;
        std::map<std::string, std::string> expected = { { "methanol", "MOH" }, { "water", "WTR" } };
        EXPECT_EQ(expected, names) << userNames;
    }
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    auto names = methanolWaterNames(&msghandler, "MOH");
    std::map<std::string, std::string> expected = { { "methanol", "MOH" }, { "water", "OXI" } };
    EXPECT_EQ(expected, names);
}

TEST(OpenMMResidueNamesTest, Errors)
{
    // SPC is water in OpenMM 8.6
    for (const auto &userNames : { "HOH", "wat", "SPC", "MOLECULE", "M-1", "MOH,MOH", "MOH WTR TIP",
                                   "methane=MOH", "water=WTR water=TIP" })
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        methanolWaterNames(&msghandler, userNames);
        EXPECT_FALSE(msghandler.ok()) << userNames;
    }
}

TEST(ParticleNamesTest, NumberedInResidue)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    auto          ff = getForceField("ACS-pg");
    ForceComputer fcomp(1e-7, 25, 1.0);
    std::vector<ACTMol> mps;
    auto dataName = gmx::test::TestFileManager::getInputFilePath("mols/methanol#water.pdb");
    initACTMol(&msghandler, dataName, ff, &fcomp, &mps);
    ASSERT_TRUE(msghandler.ok());
    // The shells get the number of their atom, the water starts at one again
    std::vector<std::string> expected = { "O1", "Os1", "H2", "Hs2", "C3", "Cs3", "H4", "Hs4",
                                          "H5", "Hs5", "H6", "Hs6",
                                          "O1", "Os1", "H2", "Hs2", "H3", "Hs3" };
    EXPECT_EQ(expected, particleNames(mps[0].atomsConst()));
}

TEST(IonChargesTest, SodiumWater)
{
    // The charge of the polarizable sodium is split over the core, the shell
    // and a virtual site, but only the atom gets the formal charge
    for (const auto &ffName : { "ACS-g", "PG-pol-ions" })
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        auto          ff = getForceField(ffName);
        ForceComputer fcomp(1e-7, 25, 1.0);
        std::vector<ACTMol> mps;
        auto dataName = gmx::test::TestFileManager::getInputFilePath("mols/sodium#water.pdb");
        initACTMol(&msghandler, dataName, ff, &fcomp, &mps);
        ASSERT_TRUE(msghandler.ok()) << ffName;
        const auto &atoms   = mps[0].atomsConst();
        auto        charges = ionCharges(atoms);
        ASSERT_EQ(atoms.size(), charges.size()) << ffName;
        for (size_t i = 0; i < atoms.size(); i++)
        {
            bool sodium = (atoms[i].pType() == ActParticle::Atom && atoms[i].element() == "Na");
            EXPECT_EQ(sodium ? "1+" : "", charges[i]) << ffName << " particle " << i;
        }
    }
}



}

}
