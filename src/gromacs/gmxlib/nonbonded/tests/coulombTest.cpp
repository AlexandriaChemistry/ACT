/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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

#include "../nb_generic.h"

#include <vector>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class CoulombTest : public gmx::test::CommandLineTestBase
{
protected:
    gmx::test::TestReferenceChecker checker_;
    std::vector<real>               r_ = { 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 };
    
    CoulombTest () : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    static void TearDownTestCase()
    {
    }
    
    void runTest(real iq, real jq, real izeta, real jzeta)
    {
        std::vector<real> energy;
        std::vector<real> force;
        
        for (const auto &r : r_)
        {
            real velec, felec;
            coulomb_gaussian(iq*jq, izeta, jzeta, r, &velec, &felec); 
            energy.push_back(velec);
            force.push_back(felec);
        }
        checker_.checkReal(iq, "iq");
        checker_.checkReal(jq, "jq");
        checker_.checkReal(izeta, "izeta");
        checker_.checkReal(jzeta, "jzeta");
        checker_.checkSequence(r_.begin(), r_.end(), "Distance");
        checker_.checkSequence(energy.begin(), energy.end(), "Energy");
        checker_.checkSequence(force.begin(), force.end(), "Force");
    }
};

TEST_F(CoulombTest, EnergyForceNaCl) 
{
    runTest(1, -1, 0.0, 0.0);
}

TEST_F(CoulombTest, EnergyForceNaClZeta100) 
{
    runTest(1, -1, 100.0, 100.0);
}

TEST_F(CoulombTest, EnergyForce0) 
{
    runTest(0.3, 0.5, 0.0, 0.0);
}

TEST_F(CoulombTest, EnergyForce1) 
{
    runTest(0.3, 0.5, 8.0, 12.0);
}

TEST_F(CoulombTest, EnergyForce2) 
{
    runTest(0.2, -0.4, 12.0, 14.0);
}

TEST_F(CoulombTest, EnergyForce3) 
{
    runTest(-2.0, -0.4, 4.0, 10.0);
}

}

} // namespace
