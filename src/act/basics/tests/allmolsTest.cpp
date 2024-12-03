/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2023
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

#include "../allmols.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

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
        }
    }

    void runTest(const std::string &inchi)
    {
        AlexandriaMols amols;

        auto *amol = amols.findInChi(inchi);
        EXPECT_FALSE(nullptr == amol);
        checkMol(amol);
    }

    void runMol(const std::string &mol)
    {
        AlexandriaMols amols;

        auto *amol = amols.findMol(mol);
        EXPECT_FALSE(nullptr == amol);
        checkMol(amol);
    }

    void runUnknown(const std::string &mol)
    {
        AlexandriaMols amols;

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

} // namespace

} // namespace
