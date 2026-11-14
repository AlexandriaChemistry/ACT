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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include <fstream>

#include "../allmols.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class AllmolsTest : public gmx::test::CommandLineTestBase
{
    protected:
    gmx::test::TestReferenceChecker checker_;
    AllmolsTest () : checker_(this->rootChecker())
    {
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    static void TearDownTestCase()
    {
    }
    
    void checkMol(const AlexandriaMol *amol)
    {
        if (amol)
        {
            checker_.checkString(amol->iupac, "IUPAC");
            checker_.checkString(amol->formula, "Formula");
            checker_.checkInteger(amol->charge, "Charge");
            checker_.checkDouble(amol->mass, "Mass");
            checker_.checkInteger(amol->mult, "Multiplicity");
            checker_.checkString(amol->cas, "CAS");
            checker_.checkInteger(amol->csid, "ChemSpiderID");
            checker_.checkInteger(amol->pubid, "PubChemID");
            checker_.checkString(amol->inchikey, "InChiKey");
            checker_.checkString(amol->pointGroup, "PointGroup");
            checker_.checkInteger(amol->symmetryNumber, "SymmetryNumber");
        }
    }

    void runTest(const std::string &inchi)
    {
        const AlexandriaMols amols;

        auto *amol = amols.findInChi(inchi);
        EXPECT_FALSE(nullptr == amol);
        checkMol(amol);
    }

    void runMol(const std::string &mol)
    {
        const AlexandriaMols amols;

        auto *amol = amols.findMol(mol);
        EXPECT_FALSE(nullptr == amol);
        checkMol(amol);
    }

    void runUnknown(const std::string &mol)
    {
        const AlexandriaMols amols;

        auto *amol = amols.findMol(mol);
        EXPECT_TRUE(nullptr == amol);
    }
};

TEST_F(AllmolsTest, Water) 
{
    runTest("InChI=1S/H2O/h1H2");
}

TEST_F(AllmolsTest, Ammonium) 
{
    runTest("InChI=1S/H3N/h1H3/p+1");
}

TEST_F(AllmolsTest, HydrogenSulfate) 
{
    runTest("InChI=1S/HO4S/c1-5(2,3)4/h(H,1,2,3)");
}

TEST_F(AllmolsTest, WaterMol) 
{
    runMol("water");
}

TEST_F(AllmolsTest, OxidaneMol) 
{
    runMol("oxidane");
}

TEST_F(AllmolsTest, AmmoniumMol) 
{
    runMol("ammonium");
}

TEST_F(AllmolsTest, HydrogenSulfateMol) 
{
    runMol("hydrogen sulfate");
}

TEST_F(AllmolsTest, Yoghurt) 
{
    runUnknown("yoghurt");
}

TEST_F(AllmolsTest, Marsipulami) 
{
    runUnknown("Marsipulami");
}

TEST(AllmolsJsonTest, ReadJson)
{
    // Molecules of alexandria.json, see src/act/python/mol_json.py, with a
    // plus-minus sign written as UTF-8 and as an escape, and no weight
    const std::string json =
        "[\n"
        "{\"alexandria\": 52, \"parentid\": 0, \"iupac\": \"1,2-dichloropropane\", \"formula\": \"C3H6Cl2\","
        " \"charge\": 0, \"multiplicity\": 1, \"weight\": 112.985, \"cas\": \"78-87-5\", \"csid\": \"6316\","
        " \"pubchem\": \"6564\", \"inchi\": \"InChI=1S/C3H6Cl2/c1-3(5)2-4/h3H,2H2,1H3\","
        " \"inchikey\": \"KNKRKFALVUDBJE-UHFFFAOYSA-N\", \"classification\": [\"primary_carbon\", \"halogenated compound\"],"
        " \"synonyms\": [\"propylene chloride\", \"1,2-Dichloropropane, (\xc2\xb1)-\"], \"filename\": \"12-dichloropropane.sdf\","
        " \"point_group\": \"C1\", \"symmetry_number\": 1, \"rotatable_bonds\": 1, \"html_name\": \"\", \"latex_name\": \"\","
        " \"author\": \"maaren\", \"timestamp\": \"2018-03-31 15:19:31\"},\n"
        "{\"alexandria\": 150, \"parentid\": 0, \"iupac\": \"2-bromobutane\", \"formula\": \"C4H9Br\","
        " \"charge\": -1, \"multiplicity\": 2, \"cas\": \"\", \"csid\": \"\","
        " \"pubchem\": \"6554\", \"inchi\": \"InChI=1S/C4H9Br/c1-3-4(2)5/h4H,3H2,1-2H3\","
        " \"inchikey\": \"UPSXAPQYNGXVBF-UHFFFAOYSA-N\", \"classification\": [],"
        " \"synonyms\": [\"2-Bromobutane, (\\u00b1)-\"], \"filename\": \"2-bromobutane.sdf\","
        " \"point_group\": \"Cs\", \"symmetry_number\": 2, \"rotatable_bonds\": 1, \"html_name\": \"\", \"latex_name\": \"\","
        " \"author\": \"maaren\", \"timestamp\": \"2018-03-31 15:25:56\"}\n"
        "]\n";
    gmx::test::TestFileManager fileManager;
    const auto fileName = fileManager.getTemporaryFilePath("alexandria.json");
    std::ofstream(fileName) << json;
    const AlexandriaMols amols(fileName);

    auto dcp = amols.findInChi("InChI=1S/C3H6Cl2/c1-3(5)2-4/h3H,2H2,1H3");
    ASSERT_TRUE(dcp != nullptr);
    EXPECT_EQ("1,2-dichloropropane", dcp->iupac);
    EXPECT_EQ("C3H6Cl2", dcp->formula);
    EXPECT_EQ(0, dcp->charge);
    EXPECT_EQ(1, dcp->mult);
    EXPECT_EQ(112.985, dcp->mass);
    EXPECT_EQ("78-87-5", dcp->cas);
    EXPECT_EQ(6316, dcp->csid);
    EXPECT_EQ(6564, dcp->pubid);
    EXPECT_EQ("KNKRKFALVUDBJE-UHFFFAOYSA-N", dcp->inchikey);
    EXPECT_EQ("C1", dcp->pointGroup);
    EXPECT_EQ(1, dcp->symmetryNumber);
    // The classes and synonyms with spaces also with hyphens, the IUPAC name and the file name
    EXPECT_EQ((std::set<std::string>{ "primary_carbon", "halogenated compound", "halogenated-compound" }), dcp->classid);
    EXPECT_EQ((std::set<std::string>{ "1,2-dichloropropane", "propylene chloride", "propylene-chloride",
                                      "1,2-Dichloropropane, (\xc2\xb1)-", "1,2-Dichloropropane,-(\xc2\xb1)-",
                                      "12-dichloropropane" }), dcp->synonyms);

    // Without weight and identifiers, as an empty column was read before
    auto bb = amols.findMol("2-Bromobutane, (\xc2\xb1)-");
    ASSERT_TRUE(bb != nullptr);
    EXPECT_EQ("2-bromobutane", bb->iupac);
    EXPECT_EQ(-1, bb->charge);
    EXPECT_EQ(2, bb->mult);
    EXPECT_EQ(0, bb->mass);
    EXPECT_EQ("", bb->cas);
    EXPECT_EQ(0, bb->csid);
    EXPECT_EQ(2, bb->symmetryNumber);
    EXPECT_TRUE(bb->classid.empty());
    EXPECT_TRUE(amols.findMol("2-bromobutane") == bb);
}

} // namespace

} // namespace
