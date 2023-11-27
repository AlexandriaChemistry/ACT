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
    void test(int crule, int ftype)
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
                                            auto pmap = evalCombinationRule(CombinationRule, ivdw, jvdw);
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
};

TEST_F (CombinationRuleTest, LjGeometric)
{
    test(eCOMB_GEOMETRIC, F_LJ);
}

TEST_F (CombinationRuleTest, LjArithmetic)
{
    test(eCOMB_ARITHMETIC, F_LJ);
}

TEST_F (CombinationRuleTest, LjLorentzBerthelot)
{
    test(eCOMB_LORENTZ_BERTHELOT, F_LJ);
}

TEST_F (CombinationRuleTest, BhamGeometric)
{
    test(eCOMB_GEOMETRIC, F_BHAM);
}

TEST_F (CombinationRuleTest, BhamArithmetic)
{
    test(eCOMB_ARITHMETIC, F_BHAM);
}

TEST_F (CombinationRuleTest, BhamKongMason)
{
    test(eCOMB_KONG_MASON, F_BHAM);
}

TEST_F (CombinationRuleTest, BhamHogerVorst)
{
    test(eCOMB_HOGERVORST, F_BHAM);
}

TEST_F (CombinationRuleTest, BhamYang)
{
    test(eCOMB_YANG, F_BHAM);
}

TEST_F (CombinationRuleTest, BhamQi)
{
    test(eCOMB_QI, F_BHAM);
}

TEST_F (CombinationRuleTest, BhamWaldmanHagler)
{
    test(eCOMB_WALDMAN_HAGLER, F_BHAM);
}

TEST_F (CombinationRuleTest, GbhamGeometric)
{
    test(eCOMB_GEOMETRIC, F_GBHAM);
}

TEST_F (CombinationRuleTest, GbhamArithmetic)
{
    test(eCOMB_ARITHMETIC, F_GBHAM);
}

TEST_F (CombinationRuleTest, GbhamKongMason)
{
    test(eCOMB_KONG_MASON, F_GBHAM);
}

TEST_F (CombinationRuleTest, GbhamHogerVorst)
{
    test(eCOMB_HOGERVORST, F_GBHAM);
}

TEST_F (CombinationRuleTest, GbhamYang)
{
    test(eCOMB_YANG, F_GBHAM);
}

TEST_F (CombinationRuleTest, GbhamQi)
{
    test(eCOMB_QI, F_GBHAM);
}

TEST_F (CombinationRuleTest, GbhamWaldmanHagler)
{
    test(eCOMB_WALDMAN_HAGLER, F_GBHAM);
}

} // namespace

} // namespace alexandria
