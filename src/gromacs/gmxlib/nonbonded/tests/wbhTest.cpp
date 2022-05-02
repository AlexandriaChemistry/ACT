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

class WBHTest : public gmx::test::CommandLineTestBase
{
protected:
    gmx::test::TestReferenceChecker checker_;
    std::vector<real>               r_;
    
    WBHTest () : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        for(int i = 50; i<=250; i++)
        {
            r_.push_back(0.002*i);
        }
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    static void TearDownTestCase()
    {
    }
    
    void runTest(real sigma, real epsilon, real gamma)
    {
        std::vector<real> energy;
        std::vector<real> force;
        
        for (const auto &r : r_)
        {
            real rsq  = r*r;
            real rinv = 1.0/r;
            real vvdw, fvdw;
            wang_buckingham(sigma, epsilon, gamma, rsq, rinv, &vvdw, &fvdw);
            energy.push_back(vvdw);
            force.push_back(fvdw);
        }
        checker_.checkReal(sigma, "Sigma");
        checker_.checkReal(epsilon, "Epsilon");
        checker_.checkReal(gamma, "Gamma");
        checker_.checkSequence(r_.begin(), r_.end(), "Distance");
        checker_.checkSequence(energy.begin(), energy.end(), "Energy");
        checker_.checkSequence(force.begin(), force.end(), "Force");
        std::vector<real> mindEdr;
        for (size_t i = 1; i < energy.size()-1; i++)
        {
            mindEdr.push_back(force[i] + (energy[i+1]-energy[i-1])/(r_[i+1]-r_[i-1]));
        }
        checker_.checkSequence(mindEdr.begin(), mindEdr.end(), "Force+dE/dr");
    }
};

TEST_F(WBHTest, EnergyForce1) 
{
    runTest(0.15, 0.5, 8.0);
}

TEST_F(WBHTest, EnergyForce2) 
{
    runTest(0.2, 0.4, 12.0);
}

}

} // namespace
