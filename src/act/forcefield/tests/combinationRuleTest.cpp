/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022,2023,2026
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

#include <cmath>
#include <map>
#include <vector>

#include <gtest/gtest.h>

#include "act/basics/msg_handler.h"
#include "act/forcefield/combinationrules.h"
#include "act/forcefield/combruleutil.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/generate_dependent.h"
#include "act/forcefield/potential.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/exceptions.h"
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
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
        checker_.setDefaultTolerance(tolerance);
        auto value = combineTwo(crule, x1, x2);
        checker_.checkReal(x1, "x1");
        checker_.checkReal(x2, "x2");
        checker_.checkReal(value, combinationRuleName(crule).c_str());
    }
    
    void testWE(double s1, double s2, double e1, double e2)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
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
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
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
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
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

TEST_F (CombinationRuleTest, GeneralizedMean_1_4_small)
{
    testGM(1.0, 4.0, 1e-5);
}

TEST_F (CombinationRuleTest, GeneralizedMean_1_4_very_small)
{
    testGM(1.0, 4.0, -1e-7);
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

TEST_F (CombinationRuleTest, GeneralizedMean_0_4_2)
{
    testGM(0.0, 4.0, 2.0);
}

TEST_F (CombinationRuleTest, GeneralizedMean_0_4_Small)
{
    testGM(0.0, 4.0, 0.01);
}

TEST_F (CombinationRuleTest, GeneralizedMean_0_4_QuiteSmall)
{
    testGM(0.0, 4.0, 0.11);
}

TEST_F (CombinationRuleTest, GeneralizedMean_1_4_VerySmall)
{
    testGM(1.0, 4.0, 0.01);
}

TEST_F (CombinationRuleTest, GeneralizedMean_1_4_QuiteSmall)
{
    testGM(1.0, 4.0, 0.11);
}

// Test: The number of combination rules from the command line that change
// those of the force field
/*! \brief Combine the van der Waals parameters of the hydrogen fluoride force field with a rule for epsilon
 * \param[in]  rule    The combination rule of epsilon
 * \param[out] mixed   The epsilon of the pair f and hf
 * \param[out] same    The epsilon of the pair f and f
 */
void combineEpsilon(const std::string &rule, double *mixed, double *same)
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/HF-LJ14-7.xml"), &pd);
    auto *vdw = pd.findForces(InteractionType::VDW);
    ASSERT_EQ(Potential::LJ14_7, vdw->potential());
    // Different sigma and epsilon for the two particle types
    const std::map<std::string, std::pair<double, double>> sigmaEpsilon = {
        { "f", { 0.32, 0.6 } }, { "hf", { 0.21, 0.15 } }
    };
    for (const auto &se : sigmaEpsilon)
    {
        auto *params = vdw->findParameters(Identifier(se.first));
        ASSERT_TRUE(params != nullptr) << se.first;
        params->find("sigma")->second.forceSetValue(se.second.first);
        params->find("epsilon")->second.forceSetValue(se.second.second);
    }
    vdw->addCombinationRule("epsilon", rule);
    MsgHandler msghandler;
    generateDependentParameter(&msghandler, &pd, true);
    EXPECT_TRUE(msghandler.ok());
    auto *pairFHf = vdw->findParameters(Identifier({ "f", "hf" }, { 1 }, CanSwap::Yes));
    auto *pairFF  = vdw->findParameters(Identifier({ "f", "f" }, { 1 }, CanSwap::Yes));
    ASSERT_TRUE(pairFHf != nullptr && pairFF != nullptr);
    *mixed = pairFHf->find("epsilon")->second.value();
    *same  = pairFF->find("epsilon")->second.value();
}

// Test: QiEpsilon, the epsilon of Qi et al. (Eqn. cr_whepsilon of the
// manual), from the sigma and epsilon of both particle types, as
// WaldmanEpsilon, when the pairs of a force field are made
TEST(CombinationRuleForceFieldTest, QiEpsilon)
{
    double s1 = 0.32, e1 = 0.6, s2 = 0.21, e2 = 0.15;
    double s13 = s1*s1*s1, s23 = s2*s2*s2;
    double expected = std::sqrt(e1*e2)*2*s13*s23/(s13*s13 + s23*s23);
    double mixed = 0, same = 0;
    combineEpsilon("QiEpsilon", &mixed, &same);
    EXPECT_NEAR(expected, mixed, 1e-12*expected);
    // The same particle types keep their epsilon
    EXPECT_NEAR(e1, same, 1e-12*e1);
    double waldman = 0, waldmanSame = 0;
    combineEpsilon("WaldmanEpsilon", &waldman, &waldmanSame);
    EXPECT_EQ(waldman, mixed);
    // Clearly not the geometric mean
    EXPECT_GT(std::abs(mixed - std::sqrt(e1*e2)), 0.01*mixed);
}

// Test: A potential without sigma and epsilon cannot use QiEpsilon
TEST(CombinationRuleForceFieldTest, QiEpsilonNeedsSigma)
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/HF-LJ14-7.xml"), &pd);
    ForceFieldParameterList buckingham("BUCKINGHAM", CanSwap::Yes);
    for (const std::string type : { "f", "hf" })
    {
        buckingham.addParameter(Identifier(type), "Abh",
                                ForceFieldParameter("kJ/mol", 5e4, 0, 1, 1e4, 1e5, Mutability::Bounded, true, true));
        buckingham.addParameter(Identifier(type), "bbh",
                                ForceFieldParameter("1/nm", 35, 0, 1, 30, 45, Mutability::Bounded, true, true));
        buckingham.addParameter(Identifier(type), "c6bh",
                                ForceFieldParameter("kJ/mol nm6", 2e-3, 0, 1, 1e-3, 4e-3, Mutability::Bounded, true, true));
    }
    buckingham.addCombinationRule("Abh", "QiEpsilon");
    buckingham.addCombinationRule("bbh", "Arithmetic");
    buckingham.addCombinationRule("c6bh", "Geometric");
    (*pd.forces())[InteractionType::VDW] = buckingham;
    MsgHandler msghandler;
    EXPECT_THROW(generateDependentParameter(&msghandler, &pd, true), gmx::InvalidInputError);
}

TEST(CombRuleUtilTest, Extract)
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACM-g1.xml"), &pd);
    auto count = [&pd](const char *rules)
    {
        CombRuleUtil         crule;
        std::vector<t_pargs> pa;
        crule.addPargs(&pa);
        *pa[0].u.c = rules;
        return crule.extract(&pd);
    };
    // Sigma is Geometric in ACM-g
    EXPECT_EQ(0, count("VANDERWAALS:sigma:Geometric"));
    EXPECT_EQ(1, count("VANDERWAALS:sigma:Arithmetic"));
    EXPECT_EQ(CombRule::Arithmetic, pd.findForcesConst(InteractionType::VDW).combinationRuleConst("sigma").rule());
    EXPECT_EQ(0, count("VANDERWAALS:sigma:Arithmetic"));
    EXPECT_EQ(2, count("VANDERWAALS:epsilon:Arithmetic VANDERWAALS:gamma:Arithmetic"));
}

} // namespace

} // namespace alexandria
