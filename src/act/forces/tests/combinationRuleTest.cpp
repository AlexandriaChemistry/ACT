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

#include "actpre.h"

#include <vector>

#include <gtest/gtest.h>

#include "act/forces/combinationrules.h"
#include "gromacs/mdtypes/md_enums.h"
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
    void testLJ(int CombinationRule)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        
        checker_.checkString(ecomb_names[CombinationRule], "CombinationRule");
        std::vector<double> sigma   = { 0.2, 0.4 };
        std::vector<double> epsilon = { 0.5, 1.0 };
        int index = 0;
        for(auto &si : sigma)
        {
            for(auto &sj: sigma)
            {
                for(auto &ei : epsilon)
                {
                    for(auto &ej : epsilon)
                    {
                        double cc6, cc12;
                        CombineLJ(CombinationRule,
                                  si, sj, ei, ej, &cc6, &cc12);
                        myCheckReal(&checker_, si, "sigmaI", index);
                        myCheckReal(&checker_, sj, "sigmaJ", index);
                        myCheckReal(&checker_, ei, "epsilonI", index);
                        myCheckReal(&checker_, ej, "epsilonJ", index);
                        myCheckReal(&checker_, cc6,  "c6", index);
                        myCheckReal(&checker_, cc12, "c12", index);
                        index++;
                    }
                }
            }
        }
    }
    
    void testBham(int CombinationRule)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        
        checker_.checkString(ecomb_names[CombinationRule], "CombinationRule");
        std::vector<double> sigma   = { 0.2, 0.4 };
        std::vector<double> epsilon = { 0.5, 1.0 };
        std::vector<double> gamma   = { 10, 15 };
        int index = 0;
        for(auto &si : sigma)
        {
            for(auto &sj: sigma)
            {
                for(auto &ei : epsilon)
                {
                    for(auto &ej : epsilon)
                    {
                        for(auto &gi : gamma)
                        {
                            for(auto &gj : gamma)
                            {
                                double sij, eij, gij;
                                CombineBham(CombinationRule,
                                            si, sj, ei, ej, gi, gj,
                                            &sij, &eij, &gij);
                                myCheckReal(&checker_, si,  "sigmaI", index);
                                myCheckReal(&checker_, sj,  "sigmaJ", index);
                                myCheckReal(&checker_, ei,  "epsilonI", index);
                                myCheckReal(&checker_, ej,  "epsilonJ", index);
                                myCheckReal(&checker_, gi,  "gammaI", index);
                                myCheckReal(&checker_, gj,  "gammaJ", index);
                                myCheckReal(&checker_, sij, "sigmaIJ", index);
                                myCheckReal(&checker_, eij, "epsilonIJ", index);
                                myCheckReal(&checker_, gij, "gammaIJ", index);
                                index++;
                            }
                        }
                    }
                }
            }
        }
    }
};

TEST_F (CombinationRuleTest, LjGeometric)
{
    testLJ(eCOMB_GEOMETRIC);
}

TEST_F (CombinationRuleTest, LjArithmetic)
{
    testLJ(eCOMB_ARITHMETIC);
}

TEST_F (CombinationRuleTest, LjLorentzBerthelot)
{
    testLJ(eCOMB_LORENTZ_BERTHELOT);
}

TEST_F (CombinationRuleTest, BhamGeometric)
{
    testBham(eCOMB_GEOMETRIC);
}

TEST_F (CombinationRuleTest, BhamArithmetic)
{
    testBham(eCOMB_ARITHMETIC);
}

TEST_F (CombinationRuleTest, BhamKongMason)
{
    testBham(eCOMB_KONG_MASON);
}

TEST_F (CombinationRuleTest, BhamHogerVorst)
{
    testBham(eCOMB_HOGERVORST);
}

TEST_F (CombinationRuleTest, BhamYang)
{
    testBham(eCOMB_YANG);
}

TEST_F (CombinationRuleTest, BhamQi)
{
    testBham(eCOMB_QI);
}

TEST_F (CombinationRuleTest, BhamWaldmanHagler)
{
    testBham(eCOMB_WALDMAN_HAGLER);
}

} // namespace

} // namespace alexandria
