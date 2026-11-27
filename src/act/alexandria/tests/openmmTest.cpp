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

#include <cstring>
#include <fstream>
#include <gtest/gtest.h>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "act/alexandria/openmm_xml.h"
#include "act/alexandria/pdbwriter.h"
#include "act/basics/act_particle.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/testutils/actmol_util.h"
#include "act/testutils/forcefield_utils.h"
#include "gromacs/utility/stringutil.h"
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
        // In single precision the charges, which are generated, differ in
        // the last digits written, the numbers are compared within a
        // tolerance then
        auto ifm = GMX_DOUBLE ?
            gmx::test::TextFileMatch(gmx::test::ExactTextMatch()).createFileMatcher() :
            gmx::test::TextFileMatch(gmx::test::NumericTextMatch()).createFileMatcher();
        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-5));
        gmx::test::TestFileManager tfm;
        MsgHandler                 msghandler; 
        msghandler.setPrintLevel(ACTStatus::Warning);
        auto ff         = getForceField(forceField);
        double rmsToler = 0.0000001;
        ForceComputer fcomp(rmsToler, 25, 1.0);
        std::vector<ACTMol> mps;
        std::string molfn("molecules/");
        molfn += fileName;
        const auto dataName = gmx::test::TestFileManager::getInputFilePath(molfn);
        initACTMol(&msghandler, dataName, ff, &fcomp, &mps);

        const auto tmpFile  = tfm.getTemporaryFilePath("xml");
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
                const auto type1 = line.substr(begin, line.find('"', begin) - begin);
                EXPECT_TRUE(shellTypes.insert(type1).second) << "Shell type " << type1 << " has more than one parent type";
            }
        }
    }
};

TEST_F(OpenMMXmlTest, Sulfate)
{
    testOpenMM("ACS-g1", "sulfate.sdf", true);
}

TEST_F(OpenMMXmlTest, Furan)
{
    testOpenMM("ACS-g1", "furan.sdf", true);
}

TEST_F(OpenMMXmlTest, Acetone)
{
    testOpenMM("ACS-g1", "acetone.sdf", true);
}

TEST_F(OpenMMXmlTest, AcetonePolMDrude)
{
    testOpenMM("ACS-pg1", "acetone.sdf", true, 0.2);
}

TEST_F(OpenMMXmlTest, SulfatePol)
{
    testOpenMM("ACS-pg1", "sulfate.sdf", true);
}

TEST_F(OpenMMXmlTest, FuranPol)
{
    testOpenMM("ACS-pg1", "furan.sdf", true);
}

TEST_F(OpenMMXmlTest, FuranPolNoNum)
{
    testOpenMM("ACS-pg1", "furan.sdf", false);
}

TEST_F(OpenMMXmlTest, AcetonePol)
{
    testOpenMM("ACS-pg1", "acetone.sdf", true);
}

TEST_F(OpenMMXmlTest, Merged)
{
    testOpenMM("ACS-g1", "merged.pdb", true);
}

// Test: Water and methanol have the same shell type on their oxygens,
// which needs different OpenMM shell types for the two parents
TEST_F(OpenMMXmlTest, MethanolWaterPol)
{
    testOpenMM("ACS-pg1", "methanol#water.pdb", true);
}

// Test: Without numbering the atom types this cannot be written
TEST(OpenMMXmlNoNumberTest, MethanolWaterPol)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    auto          ff = getForceField("ACS-pg1");
    ForceComputer fcomp(1e-7, 25, 1.0);
    std::vector<ACTMol> mps;
    const auto dataName = gmx::test::TestFileManager::getInputFilePath("molecules/methanol#water.pdb");
    initACTMol(&msghandler, dataName, ff, &fcomp, &mps);
    ASSERT_TRUE(msghandler.ok());
    gmx::test::TestFileManager tfm;
    writeOpenMM(&msghandler, tfm.getTemporaryFilePath("xml"), ff, mps, 0.1, 1, false);
    EXPECT_FALSE(msghandler.ok());
}

TEST_F(OpenMMXmlTest, WaterHyperPol)
{
    testOpenMM("water_hyperpol1", "water_dimer.pdb", true);
}

// Test: Lennard-Jones and point charges, the NonbondedForce of OpenMM only
TEST_F(OpenMMXmlTest, WaterSPCE)
{
    testOpenMM("SPCE", "water_dimer.pdb", true);
}

// Test: The lone pairs of TIP5P, VSITE3OUT
TEST_F(OpenMMXmlTest, WaterTIP5P)
{
    testOpenMM("TIP5P", "water_dimer.pdb", true);
}

// Test: The Buckingham potential and VSITE3S of BK3
TEST_F(OpenMMXmlTest, WaterBK3)
{
    testOpenMM("BK3", "water_dimer.pdb", true);
}

// Test: The generalized Buckingham potential
TEST_F(OpenMMXmlTest, AcetoneGBHam)
{
    testOpenMM("ACS-pg-gbham", "acetone.sdf", true);
}

// Test: The buffered 14-7 potential
TEST_F(OpenMMXmlTest, HydrogenFluorideLJ147)
{
    testOpenMM("HF-LJ14-7", "hydrogen-fluoride#hydrogen-fluoride.pdb", true);
}

// Test: Urey-Bradley angles
TEST_F(OpenMMXmlTest, AcetoneUreyBradley)
{
    testOpenMM("ACS-g-dih", "acetone.sdf", true);
}

// Test: A linear angle, of which the parameter a depends on the direction
TEST_F(OpenMMXmlTest, AcetonitrileLinearAngle)
{
    testOpenMM("ACS-g1", "acetonitrile.sdf", true);
}

// Test: Flexible TIP4P water with the cross terms between the bonds and
// the angle
TEST_F(OpenMMXmlTest, WaterCrossBondAngles)
{
    testOpenMM("TIP4P-Flex", "water.sdf", true);
}

/*! \brief Write the OpenMM force field file of a compound
 * \param[in] forceField The name of the force field
 * \param[in] fileName   The file with the compound in the test data
 * \return the text of the file, empty if it could not be written
 */
static std::string openMMText(const std::string &forceField,
                              const std::string &fileName)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    auto          ff = getForceField(forceField);
    ForceComputer fcomp(1e-7, 25, 1.0);
    std::vector<ACTMol> mps;
    const auto dataName = gmx::test::TestFileManager::getInputFilePath("molecules/" + fileName);
    initACTMol(&msghandler, dataName, ff, &fcomp, &mps);
    EXPECT_TRUE(msghandler.ok());
    gmx::test::TestFileManager tfm;
    const auto tmpFile = tfm.getTemporaryFilePath("xml");
    writeOpenMM(&msghandler, tmpFile, ff, mps, 0.1, 1, true);
    EXPECT_TRUE(msghandler.ok());
    const std::ifstream xml(tmpFile);
    std::stringstream text;
    text << xml.rdbuf();
    return text.str();
}

//! \return how often a text is found in another text
static size_t count(const std::string &text, const std::string &what)
{
    size_t n = 0;
    for (auto pos = text.find(what); pos != std::string::npos; pos = text.find(what, pos + what.size()))
    {
        n++;
    }
    return n;
}

// Test: The linear angle has its energy expression and its parameters, and
// a script turns the force into one that OpenMM can evaluate, with the
// atoms in the order of the parameter a.
TEST(OpenMMAngleForceTest, LinearAngle)
{
    const auto text = openMMText("ACS-g1", "acetonitrile.sdf");
    ASSERT_FALSE(text.empty());
    EXPECT_EQ(1u, count(text, "<CustomAngleForce energy=\"0.5*klin*(a*d12^2+(1-a)*d23^2-a*(1-a)*d13^2); "
                        "d12=distance(p1,p2); d23=distance(p2,p3); d13=distance(p1,p3)\">"));
    EXPECT_EQ(1u, count(text, "<GlobalParameter name=\"pot-LINEAR_ANGLES\" defaultValue=\"1\"/>"));
    EXPECT_EQ(1u, count(text, "<PerAngleParameter name=\"a\"/>"));
    EXPECT_EQ(1u, count(text, "<PerAngleParameter name=\"klin\"/>"));
    EXPECT_EQ(1u, count(text, "<Angle class1=\"c3_b\" class2=\"c1_b\" class3=\"n1_b\" a=\""));
    EXPECT_EQ(1u, count(text, "<Script>"));
    EXPECT_EQ(1u, count(text, "openmm.CustomCompoundBondForce(3, energy)"));
    EXPECT_EQ(1u, count(text, "    \"LINEAR_ANGLES\": [\n        ( \"c3_b\", \"c1_b\", \"n1_b\" ),\n    ],\n"));
    for (const auto &type : { "c3_1", "c1_1", "n1_1" })
    {
        EXPECT_EQ(1u, count(text, gmx::formatString("    \"%s\": \"", type))) << type;
    }
    // The hydrogen atoms are not in an angle with a direction
    EXPECT_EQ(0u, count(text, "    \"h_1\": \""));
    // No force with an energy of zero
    EXPECT_EQ(0u, count(text, "CustomAngleForce energy=\"0\""));
}

// Test: The Urey-Bradley angle has its energy expression and its
// parameters. It is the same in both directions.
TEST(OpenMMAngleForceTest, UreyBradley)
{
    const auto text = openMMText("ACS-g-dih", "acetone.sdf");
    ASSERT_FALSE(text.empty());
    EXPECT_EQ(1u, count(text, "<CustomAngleForce energy=\"0.5*kt*(angle(p1,p2,p3)-ub_angle)^2"
                        "+0.5*kub*(distance(p1,p3)-r13)^2\">"));
    EXPECT_EQ(1u, count(text, "<GlobalParameter name=\"pot-UREY_BRADLEY_ANGLES\" defaultValue=\"1\"/>"));
    for (const auto &param : { "kt", "ub_angle", "r13", "kub" })
    {
        EXPECT_EQ(1u, count(text, gmx::formatString("<PerAngleParameter name=\"%s\"/>", param))) << param;
    }
    EXPECT_EQ(4u, count(text, " ub_angle=\""));
    EXPECT_EQ(1u, count(text, "<Script>"));
    EXPECT_EQ(1u, count(text, "ordered = {\n}"));
    EXPECT_EQ(1u, count(text, "atom_class = {\n}"));
}

// Test: The cross terms between the bonds and the angle of flexible TIP4P
// water have their energy expression and the six parameters of the angle.
// The angle of water is the same in both directions.
TEST(OpenMMAngleForceTest, CrossBondAngles)
{
    const auto text = openMMText("TIP4P-Flex", "water.sdf");
    ASSERT_FALSE(text.empty());
    EXPECT_EQ(1u, count(text, "<CustomAngleForce energy=\"krt*(d13-r3e)*(d12-r1e+d23-r2e)"
                        "+krr*(d12-r1e)*(d23-r2e)+0.5*kub*(d13-r3e)^2; "
                        "d12=distance(p1,p2); d23=distance(p2,p3); d13=distance(p1,p3)\">"));
    EXPECT_EQ(1u, count(text, "<GlobalParameter name=\"pot-CROSS_BOND_ANGLES\" defaultValue=\"1\"/>"));
    for (const auto &param : { "r1e", "r2e", "r3e", "krt", "krr", "kub" })
    {
        EXPECT_EQ(1u, count(text, gmx::formatString("<PerAngleParameter name=\"%s\"/>", param))) << param;
    }
    // Lengths in nm and force constants in kJ/mol nm^2, as OpenMM has them
    EXPECT_EQ(1u, count(text, "r1e=\"0.09347\" r2e=\"0.09347\" r3e=\"0.1526\" krt=\"-81901.5\" "
                        "krr=\"52995.1\" kub=\"103098.98\"/>"));
    EXPECT_EQ(1u, count(text, "<Script>"));
    EXPECT_EQ(1u, count(text, "openmm.CustomCompoundBondForce(3, energy)"));
    EXPECT_EQ(1u, count(text, "ordered = {\n}"));
}

// Test: A force without angles is not written, and neither is the script
// if no force needs it.
TEST(OpenMMAngleForceTest, NoAngles)
{
    const auto text = openMMText("ACS-g1", "acetone.sdf");
    ASSERT_FALSE(text.empty());
    EXPECT_EQ(1u, count(text, "<HarmonicAngleForce"));
    EXPECT_EQ(0u, count(text, "CustomAngleForce"));
    EXPECT_EQ(0u, count(text, "<Script>"));
}

// Test: A proper dihedral has the periodicity, the phase and the force
// constant under the names that the PeriodicTorsionForce of OpenMM reads
TEST(OpenMMDihedralTest, ProperDihedrals)
{
    const auto text = openMMText("ACS-g1", "acetone.sdf");
    ASSERT_FALSE(text.empty());
    const auto nproper = count(text, "<Proper ");
    EXPECT_GT(nproper, 0u);
    EXPECT_EQ(nproper, count(text, "\" periodicity1=\"3\" phase1=\""));
    EXPECT_EQ(nproper, count(text, "\" k1=\"5\"/>"));
    for (const auto &actName : { " phi0=\"", " kp=\"", " mult=\"" })
    {
        EXPECT_EQ(0u, count(text, actName)) << actName;
    }
}

// Test: The harmonic improper dihedrals are a CustomTorsionForce with
// their energy expression, and the central atom is the first one, also
// if the identifier in the force field has it as the last one
TEST(OpenMMDihedralTest, ImproperDihedrals)
{
    const auto text = openMMText("ACS-g1", "furan.sdf");
    ASSERT_FALSE(text.empty());
    EXPECT_EQ(1u, count(text, "<CustomTorsionForce energy=\"0.5*kimp*theta^2\" ordering=\"charmm\">"));
    EXPECT_EQ(1u, count(text, "<PerTorsionParameter name=\"kimp\"/>"));
    EXPECT_EQ(2u, count(text, "<Improper "));
    // The carbon atoms next to the oxygen atom
    EXPECT_EQ(1u, count(text, "<Improper class1=\"c2_b\" class2=\"h_b\" class3=\"c2_b\" class4=\"o2_b\" kimp=\""));
    // The other two, of which the identifier is c2_b:c2_b~h_b:c2_b
    EXPECT_EQ(1u, count(text, "<Improper class1=\"c2_b\" class2=\"h_b\" class3=\"c2_b\" class4=\"c2_b\" kimp=\""));
}

// Test: No force for improper dihedrals if the compound has none
TEST(OpenMMDihedralTest, NoImproperDihedrals)
{
    const auto text = openMMText("ACS-g1", "acetonitrile.sdf");
    ASSERT_FALSE(text.empty());
    EXPECT_EQ(0u, count(text, "CustomTorsionForce"));
    EXPECT_EQ(0u, count(text, "<Improper"));
}

// Test: The induction correction
TEST_F(OpenMMXmlTest, WaterInductionCorrection)
{
    testOpenMM("ACS-pol", "water_dimer.pdb", true);
}

// Test: VSITE1 on the polarizable sodium
TEST_F(OpenMMXmlTest, SodiumWaterVsite1)
{
    testOpenMM("PG-pol-ions", "sodium#water.pdb", true);
}

// Test: Hua bonds
TEST_F(OpenMMXmlTest, AcetateHua)
{
    testOpenMM("PC+GV-elec", "acetate.sdf", true);
}

// Test: VSITE3 on water
TEST_F(OpenMMXmlTest, WaterVsite3)
{
    testOpenMM("ACS-pg-vs3", "water_dimer.pdb", true);
}

// Test: OpenMM has no four particle virtual sites, hence the amine of
// methylamine cannot be written
TEST(OpenMMXmlVsite4Test, Methylamine)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    auto          ff = getForceField("ACS-pg-vs3");
    ForceComputer fcomp(1e-7, 25, 1.0);
    std::vector<ACTMol> mps;
    const auto dataName = gmx::test::TestFileManager::getInputFilePath("molecules/methylamine.sdf");
    initACTMol(&msghandler, dataName, ff, &fcomp, &mps);
    ASSERT_TRUE(msghandler.ok());
    gmx::test::TestFileManager tfm;
    writeOpenMM(&msghandler, tfm.getTemporaryFilePath("xml"), ff, mps, 0.1, 1, true);
    EXPECT_FALSE(msghandler.ok());
}

//! \return the residue names for methanol and water, with userNames
static std::map<std::string, std::string> methanolWaterNames(MsgHandler        *msghandler,
                                                             const std::string &userNames)
{
    auto          ff = getForceField("ACS-g1");
    ForceComputer fcomp(1e-7, 25, 1.0);
    std::vector<ACTMol> mps;
    const auto dataName = gmx::test::TestFileManager::getInputFilePath("molecules/methanol#water.pdb");
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
    const auto names = methanolWaterNames(&msghandler, "");
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
        const auto names = methanolWaterNames(&msghandler, userNames);
        EXPECT_TRUE(msghandler.ok()) << userNames;
        std::map<std::string, std::string> expected = { { "methanol", "MOH" }, { "water", "WTR" } };
        EXPECT_EQ(expected, names) << userNames;
    }
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    const auto names = methanolWaterNames(&msghandler, "MOH");
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
    auto          ff = getForceField("ACS-pg1");
    ForceComputer fcomp(1e-7, 25, 1.0);
    std::vector<ACTMol> mps;
    const auto dataName = gmx::test::TestFileManager::getInputFilePath("molecules/methanol#water.pdb");
    initACTMol(&msghandler, dataName, ff, &fcomp, &mps);
    ASSERT_TRUE(msghandler.ok());
    // The shells get the number of their atom, the water starts at one again
    const std::vector<std::string> expected = { "O1", "Os1", "H2", "Hs2", "C3", "Cs3", "H4", "Hs4",
                                          "H5", "Hs5", "H6", "Hs6",
                                          "O1", "Os1", "H2", "Hs2", "H3", "Hs3" };
    EXPECT_EQ(expected, particleNames(mps[0].atomsConst()));
}

TEST(IonChargesTest, SodiumWater)
{
    // The charge of the polarizable sodium is split over the core, the shell
    // and a virtual site, but only the atom gets the formal charge
    for (const auto &ffName : { "ACS-g1", "PG-pol-ions" })
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        auto          ff = getForceField(ffName);
        ForceComputer fcomp(1e-7, 25, 1.0);
        std::vector<ACTMol> mps;
        const auto dataName = gmx::test::TestFileManager::getInputFilePath("molecules/sodium#water.pdb");
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
