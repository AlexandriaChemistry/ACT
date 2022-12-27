/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022
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

#include <math.h>

#include <gtest/gtest.h>

#include "alexandria/secondvirial.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class SphereIntegratorTest : public gmx::test::CommandLineTestBase
{
protected:
    void test(const std::vector<std::pair<double, double>> &data,
              const double expectedValue)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        // Test integration routine.
        std::vector<double> dB;
        double              B = 0;
        for(size_t i = 0; i < data.size()-1; i++)
        {
            double ddB = sphereIntegrator(data[i].first, data[i+1].first, data[i].second, data[i+1].second);
            dB.push_back(ddB);
            B += ddB;
        }
        checker_.checkReal(B, "Numerical Integral");
        checker_.checkReal(expectedValue, "Analytical Integral");
        checker_.checkSequence(dB.begin(), dB.end(), "Components");
    }
};

TEST_F (SphereIntegratorTest, Constant)
{
    std::vector<std::pair<double, double>> data;
    for(size_t i = 0; i <= 100; i++)
    {
        data.push_back({ i*0.01, 1.0 });
    }
    // Expected analytical integral 4/3 Pi r^3
    double expected = 4*M_PI/3.0;
    test(data, expected);
}

TEST_F (SphereIntegratorTest, Exponential)
{
    std::vector<std::pair<double, double>> data;
    for(size_t i = 0; i <= 100; i++)
    {
        double x = i*0.01;
        data.push_back({ x, std::exp(-x) });
    }
    // Expected analytical integral (from Mathematica)
    // 4 (2 - 5/E) * Pi
    double expected = 4*M_PI*(2-5/std::exp(1));
    test(data, expected);
}

// #endif

} // namespace

} // namespace alexandria
