/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2024,2025
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
#include <cstdlib>

#include "../atomization_energy.h"

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace
{

//! Class to test a rotation algorithm
class AtomizationEnergyTest : public gmx::test::CommandLineTestBase
{
protected:
    //! Checking data structure
    gmx::test::TestReferenceChecker checker_;
    alexandria::AtomizationEnergy   atomenergy_;
    
    //! Init set tolecrance
    AtomizationEnergyTest () : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        atomenergy_.read();
    }
    void dump()
    {
        atomenergy_.dump("test.csv");
    }
    /*! \brief test the atomization energies
     * This is done by rotating tensor p on tensor q and computing the
     * RMSD after rotating.
     * \param[in] p First tensor
     * \param[in] q Second tensor
     */
    void test(const std::string &elem, double T, int charge = 0)
    {
        std::map<std::string, std::vector<std::string>> props = 
            {
                { "exp", { "H(0)-H(T)", "S0(T)", "DHf(T)" } },
                { "B3LYP/aug-cc-pvtz", { "E0" } },
                { "G4", { "G4(0K)" } },
                { "W1U", { "W1U(0K)" } }
            };
        std::string unit, ref;
        for (const auto &source : props)
        {
            auto schk = checker_.checkCompound("source", source.first.c_str());
            for (const auto &prop : source.second)
            {
                auto pchk = schk.checkCompound("property", prop.c_str());
                double temp = T;
                if ("exp" != source.first)
                {
                    temp = 0;
                }
                auto term = atomenergy_.term(elem, charge, source.first, prop, temp, &unit, &ref);
                if (term != 0)
                {
                    pchk.checkString(elem, "element");
                    pchk.checkInt64(charge, "charge");
                    pchk.checkDouble(T, "temperature");
                    pchk.checkDouble(term, "value");
                    pchk.checkString(unit, "unit");
                    pchk.checkString(ref, "reference");
                }
            }
        }
    }
};

//TEST_F (AtomizationEnergyTest, Dump)
//{
//    dump();
//}

TEST_F (AtomizationEnergyTest, OxygenRT)
{
    test("O", 298.15);
}

TEST_F (AtomizationEnergyTest, OxyideRT)
{
    test("O", 298.15, -2);
}

TEST_F (AtomizationEnergyTest, Oxygen500)
{
    test("O", 500);
}

TEST_F (AtomizationEnergyTest, NitrogenRT)
{
    test("N", 298.15);
}

TEST_F (AtomizationEnergyTest, Nitrogen700)
{
    test("N", 700);
}

}
