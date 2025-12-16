/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022,2023
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

#include "actpre.h"

#include <vector>

#include <gtest/gtest.h>

#include "act/forces/combinationrules.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class CombinationRuleTest : public gmx::test::CommandLineTestBase
{
protected:
    void testLow(CombRule crule, double x1, double x2)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        auto value = combineTwo(crule, x1, x2);
        checker_.checkReal(x1, "x1");
        checker_.checkReal(x2, "x2");
        checker_.checkReal(value, combinationRuleName(crule).c_str());
    }
    
    void testWE(double s1, double s2, double e1, double e2)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        auto value = combineWaldmanEpsilon(e1, e2, s1, s2);
        checker_.checkReal(e1, "epsilon1");
        checker_.checkReal(e2, "epsilon2");
        checker_.checkReal(s1, "sigma1");
        checker_.checkReal(s2, "sigma2");
        checker_.checkReal(value, combinationRuleName(CombRule::WaldmanEpsilon).c_str());
    }
    
    void testHS(double e1, double e2, double g1, double g2, double s1, double s2)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        auto value = combineHogervorstSigma(e1, e2, g1, g2, s1, s2);
        checker_.checkReal(e1, "epsilon1");
        checker_.checkReal(e2, "epsilon2");
        checker_.checkReal(g1, "gamma1");
        checker_.checkReal(g2, "gamma2");
        checker_.checkReal(s1, "sigma1");
        checker_.checkReal(s2, "sigma2");
        checker_.checkReal(value, combinationRuleName(CombRule::HogervorstSigma).c_str());
    }

    void testGM(double x1, double x2, double exponent)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        auto value = combineGeneralizedMean(x1, x2, exponent);
        checker_.checkReal(x1, "x1");
        checker_.checkReal(x2, "x2");
        checker_.checkReal(exponent, "exponent");
        checker_.checkReal(value, combinationRuleName(CombRule::GeneralizedMean).c_str());
    }
};

TEST_F (CombinationRuleTest, Geometric_1_2)
{
    testLow(CombRule::Geometric, 1.0, 2.0);
}

TEST_F (CombinationRuleTest, Geometric_0_2)
{
    testLow(CombRule::Geometric, 0.0, 2.0);
}

TEST_F (CombinationRuleTest, Arithmetic_1_2)
{
    testLow(CombRule::Arithmetic, 1.0, 2.0);
}

TEST_F (CombinationRuleTest, Arithmetic_0_2)
{
    testLow(CombRule::Arithmetic, 0.0, 2.0);
}

TEST_F (CombinationRuleTest, Volumetric_1_2)
{
    testLow(CombRule::Volumetric, 1.0, 2.0);
}

TEST_F (CombinationRuleTest, Volumetric_0_2)
{
    testLow(CombRule::Volumetric, 0.0, 2.0);
}

TEST_F (CombinationRuleTest, InverseSquare_1_2)
{
    testLow(CombRule::InverseSquare, 1.0, 2.0);
}

TEST_F (CombinationRuleTest, InverseSquare_0_2)
{
    testLow(CombRule::InverseSquare, 0.0, 2.0);
}

TEST_F (CombinationRuleTest, Yang_1_2)
{
    testLow(CombRule::Yang, 1.0, 2.0);
}

TEST_F (CombinationRuleTest, Yang_0_2)
{
    testLow(CombRule::Yang, 0.0, 2.0);
}

TEST_F (CombinationRuleTest, HogervorstEpsilon_1_2)
{
    testLow(CombRule::HogervorstEpsilon, 1.0, 2.0);
}

TEST_F (CombinationRuleTest, HogervorstEpsilon_0_2)
{
    testLow(CombRule::HogervorstEpsilon, 0.0, 2.0);
}

TEST_F (CombinationRuleTest, WaldmanSigma_1_2)
{
    testLow(CombRule::WaldmanSigma, 1.0, 2.0);
}

TEST_F (CombinationRuleTest, WaldmanSigma_0_2)
{
    testLow(CombRule::WaldmanSigma, 0.0, 2.0);
}

TEST_F (CombinationRuleTest, QiSigma_1_2)
{
    testLow(CombRule::QiSigma, 1.0, 2.0);
}

TEST_F (CombinationRuleTest, QiSigma_0_2)
{
    testLow(CombRule::QiSigma, 0.0, 2.0);
}

TEST_F (CombinationRuleTest, HalgrenEpsilon_1_2)
{
    testLow(CombRule::HalgrenEpsilon, 1.0, 2.0);
}

TEST_F (CombinationRuleTest, HalgrenEpsilon_0_2)
{
    testLow(CombRule::HalgrenEpsilon, 0.0, 2.0);
}

TEST_F (CombinationRuleTest, WaldmanEpsilon_1_2_3_4)
{
    testWE(1,2,3,4);
}

TEST_F (CombinationRuleTest, WaldmanEpsilon_1_2_3_0)
{
    testWE(1,2,3,0);
}

TEST_F (CombinationRuleTest, WaldmanEpsilon_1_0_3_4)
{
    testWE(1,0,3,4);
}

TEST_F (CombinationRuleTest, HogervorstSigma_1_2_15_16_3_4)
{
    testHS(1,2,15,16,3,4);
}

TEST_F (CombinationRuleTest, HogervorstSigma_1_2_15_16_3_0)
{
    testHS(1,2,15,16,3,0);
}

TEST_F (CombinationRuleTest, HogervorstSigma_1_0_15_16_3_4)
{
    testHS(1,0,15,16,3,4);
}

TEST_F (CombinationRuleTest, GeneralizedMean_1_2_1)
{
    testGM(1.0, 2.0, 1.0);
}

TEST_F (CombinationRuleTest, GeneralizedMean_1_4_0)
{
    testGM(1.0, 4.0, 0.0);
}

TEST_F (CombinationRuleTest, GeneralizedMean_1_4_2)
{
    testGM(1.0, 4.0, 2.0);
}

TEST_F (CombinationRuleTest, GeneralizedMean_3_4_2)
{
    testGM(3.0, 4.0, 2.0);
}

TEST_F (CombinationRuleTest, GeneralizedMean_3_4__2)
{
    testGM(3.0, 4.0, -2.0);
}

} // namespace

} // namespace alexandria
