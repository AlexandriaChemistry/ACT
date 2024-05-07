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

static void myCheckReal(gmx::test::TestReferenceChecker *checker,
                        double                           value,
                        const char                      *label,
                        int                              index)
{
    std::string mylabel = gmx::formatString("%s-%d", label, index);
    checker->checkReal(value, mylabel.c_str());
}

class CombinationRuleTest : public gmx::test::CommandLineTestBase
{
protected:
    void test(int crule, Potential ftype)
    {
        auto CombinationRule = oldCombinationRule(ecomb_names[crule], ftype);

        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        
        for(const auto &cr : CombinationRule)
        {
            checker_.checkString(combinationRuleName(cr.second).c_str(), cr.first.c_str());
        }
        std::vector<double> sigma   = { 0.2, 0.4 };
        std::vector<double> epsilon = { 0.5, 1.0 };
        std::vector<double> gamma   = { 10, 13   };
        std::vector<double> delta   = { 0, 0.2   };
        int index = 0;
        std::string csigma("sigma");
        if (CombinationRule.end() != CombinationRule.find("rmin"))
        {
            csigma.assign("rmin");
        }
        bool haveGamma = (CombinationRule.end() != CombinationRule.find("gamma"));
        bool haveDelta = (CombinationRule.end() != CombinationRule.find("delta"));

        for(auto &si : sigma)
        {
            for(auto &ei : epsilon)
            {
                bool sigEpsDone = false;
                for(auto &gi : gamma)
                {
                    bool sigEpsGamDone = false;
                    for(auto &di : delta)
                    {
                        ForceFieldParameterMap ivdw;
                        ivdw.insert({ csigma, ForceFieldParameter("", si, 0, 1, si, si, Mutability::Bounded, false, false) });
                        ivdw.insert({ "epsilon", ForceFieldParameter("", ei, 0, 1, ei, ei, Mutability::Bounded, false, false) });
                        if (haveGamma)
                        {
                            ivdw.insert({ "gamma", ForceFieldParameter("", gi, 0, 1, gi, gi, Mutability::Bounded, false, false) });
                        }
                        if (haveDelta)
                        {
                            ivdw.insert({ "delta", ForceFieldParameter("", di, 0, 1, di, di, Mutability::Bounded, false, false) });
                        }
                        for(auto &sj: sigma)
                        {
                            for(auto &ej : epsilon)
                            {
                                for(auto &gj: gamma)
                                {
                                    for(auto &dj : delta)
                                    {
                                        ForceFieldParameterMap jvdw;
                                        jvdw.insert({ csigma, ForceFieldParameter("", sj, 0, 1, sj, sj, Mutability::Bounded, false, false) });
                                        jvdw.insert({ "epsilon", ForceFieldParameter("", ej, 0, 1, ej, ej, Mutability::Bounded, false, false) });
                                        if (haveGamma)
                                        {
                                            jvdw.insert({ "gamma", ForceFieldParameter("", gj, 0, 1, gj, gj, Mutability::Bounded, false, false) });
                                        
                                            if (haveDelta)
                                            {
                                                jvdw.insert({ "delta", ForceFieldParameter("", dj, 0, 1, dj, dj, Mutability::Bounded, false, false) });
                                            }
                                        }
                                        if ((haveGamma && haveDelta) ||
                                            (haveGamma && !sigEpsGamDone) ||
                                            (!sigEpsDone))
                                        {
                                            auto pmap = evalCombinationRule(ftype, CombinationRule, ivdw, jvdw);
                                            for(const auto &pm : pmap)
                                            {
                                                myCheckReal(&checker_, pm.second.value(), pm.first.c_str(), index);
                                            }
                                            index++;
                                            sigEpsGamDone = true;
                                            sigEpsDone = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
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
};

TEST_F (CombinationRuleTest, LjGeometric)
{
    test(eCOMB_GEOMETRIC, Potential::LJ12_6);
}

TEST_F (CombinationRuleTest, LjArithmetic)
{
    test(eCOMB_ARITHMETIC, Potential::LJ12_6);
}

TEST_F (CombinationRuleTest, LjLorentzBerthelot)
{
    test(eCOMB_LORENTZ_BERTHELOT, Potential::LJ12_6);
}

TEST_F (CombinationRuleTest, BhamGeometric)
{
    test(eCOMB_GEOMETRIC, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, BhamArithmetic)
{
    test(eCOMB_ARITHMETIC, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, BhamKongMason)
{
    test(eCOMB_KONG_MASON, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, BhamHogerVorst)
{
    test(eCOMB_HOGERVORST, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, BhamYang)
{
    test(eCOMB_YANG, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, BhamQi)
{
    test(eCOMB_QI, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, BhamQi2)
{
    test(eCOMB_QI_2, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, BhamQyQy)
{
    test(eCOMB_QYQY, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, BhamQKmQG)
{
    test(eCOMB_QKmQG, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, BhamWaldmanHagler)
{
    test(eCOMB_WALDMAN_HAGLER, Potential::WANG_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamGeometric)
{
    test(eCOMB_GEOMETRIC, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamArithmetic)
{
    test(eCOMB_ARITHMETIC, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamKongMason)
{
    test(eCOMB_KONG_MASON, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamHogerVorst)
{
    test(eCOMB_HOGERVORST, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamYang)
{
    test(eCOMB_YANG, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamQi)
{
    test(eCOMB_QI, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamQi2)
{
    test(eCOMB_QI_2, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamQyQy)
{
    test(eCOMB_QYQY, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamQKmQG)
{
    test(eCOMB_QKmQG, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, GbhamWaldmanHagler)
{
    test(eCOMB_WALDMAN_HAGLER, Potential::GENERALIZED_BUCKINGHAM);
}

TEST_F (CombinationRuleTest, LJ14_7Geometric)
{
    test(eCOMB_GEOMETRIC, Potential::LJ14_7);
}

TEST_F (CombinationRuleTest, LJ14_7Arithmetic)
{
    test(eCOMB_ARITHMETIC, Potential::LJ14_7);
}

TEST_F (CombinationRuleTest, LJ14_7KongMason)
{
    test(eCOMB_KONG_MASON, Potential::LJ14_7);
}

TEST_F (CombinationRuleTest, LJ14_7HogerVorst)
{
    test(eCOMB_HOGERVORST, Potential::LJ14_7);
}

TEST_F (CombinationRuleTest, LJ14_7Yang)
{
    test(eCOMB_YANG, Potential::LJ14_7);
}

TEST_F (CombinationRuleTest, LJ14_7Qi)
{
    test(eCOMB_QI, Potential::LJ14_7);
}

TEST_F (CombinationRuleTest, LJ14_7Qi2)
{
    test(eCOMB_QI_2, Potential::LJ14_7);
}

TEST_F (CombinationRuleTest, LJ14_7QyQy)
{
    test(eCOMB_QYQY, Potential::LJ14_7);
}

TEST_F (CombinationRuleTest, LJ14_7QKmQG)
{
    test(eCOMB_QKmQG, Potential::LJ14_7);
}

TEST_F (CombinationRuleTest, LJ14_7WaldmanHagler)
{
    test(eCOMB_WALDMAN_HAGLER, Potential::LJ14_7);
}

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

TEST_F (CombinationRuleTest, HogervorstSigma_1_0_15_6_3_4)
{
    EXPECT_THROW(testHS(1,0,15,6,3,4), gmx::InvalidInputError);
}

} // namespace

} // namespace alexandria
