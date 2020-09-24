/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2020
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

#include "programs/alexandria/forcefieldparameterlist.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class ForceFieldParameterListTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;

        ForceFieldParameterListTest () : checker_(this->rootChecker())
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

        void checkOption(const ForceFieldParameterList &ff, const std::string &option)
        {
            checker_.checkString(ff.function(), "function");
            checker_.checkString(ff.optionValue(option), "optionValue");
        }
        void runTest(ForceFieldParameterList *ff, bool modify)
        {
            std::vector<std::string> atoms = { "c2", "h3", "c3" };
            
            checker_.checkString(ff->function(), "function");
            for(const auto &atom : atoms)
            {
                if (ff->parameterExists(atom))
                {
                    if (modify)
                    {
                        auto params = ff->findParameter(atom);
                        for (auto &p : params)
                        {
                            switch (p.mutability()) 
                            {
                            case Mutability::Free:
                                {
                                    p.setValue(2*p.value());
                                    break;
                                }
                            case Mutability::Bounded:
                                {
                                    p.setValue(0.5*(p.minimum()+p.maximum()));
                                    checker_.checkDouble(p.value(), "IntermediateValue");
                                    if (p.strict())
                                    {
                                        EXPECT_THROW(p.setValue(2*p.maximum()), gmx::InvalidInputError);
                                        EXPECT_THROW(p.setValue(0.5*p.minimum()), gmx::InvalidInputError);
                                    }
                                    else
                                    {
                                        p.setValue(2*p.maximum());
                                        checker_.checkDouble(p.value(), "MaxValue");
                                    }
                                    break;
                                }
                            case Mutability::Fixed:
                                {
                                    if (p.strict())
                                    {
                                        EXPECT_THROW(p.setValue(2*p.maximum()), gmx::InvalidInputError);
                                        EXPECT_THROW(p.setValue(0.5*p.minimum()), gmx::InvalidInputError);
                                    }
                                    else
                                    {
                                        p.setValue(2*p.maximum());
                                        checker_.checkDouble(p.value(), "UnchangedValue");
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    else
                    {
                        auto p = ff->searchParameter(atom);
                        checker_.checkInteger(p.size(), atom.c_str());
                    }
                }
                else
                {
                    EXPECT_THROW(auto p = ff->searchParameter(atom), gmx::InvalidInputError);
                }
            }
        }

};

TEST_F (ForceFieldParameterListTest, AddOptionWorks) {
    ForceFieldParameterList ff("BONDS");
    std::string option("Order");
    ff.addOption(option, "No");
    EXPECT_TRUE(ff.optionExists(option));
    checkOption(ff, option);
}

TEST (ForceFieldParameterListSimpleTest, OptionMissing) {
    ForceFieldParameterList ff("ANGLES");
    std::string option("Order");
    EXPECT_FALSE(ff.optionExists(option));
}

TEST (ForceFieldParameterListSimpleTest, OptionMissingThrows) {
    ForceFieldParameterList ff("ANGLES");
    std::string option("Order");
    EXPECT_THROW(auto val = ff.optionValue(option), gmx::InvalidInputError);
}

TEST_F (ForceFieldParameterListTest, AddParameter) {
    ForceFieldParameterList ff("Torsion");
    ff.addParameter("h3", ForceFieldParameter("sigma", "nm", 12.0, 0.3, 10.0, 18.0, Mutability::Free, false));
    ff.addParameter("c2", ForceFieldParameter("gamma", "", 11.0, 0.25, 8.0, 15.0, Mutability::Fixed, true));
    ff.addParameter("h3", ForceFieldParameter("epsilon", "kJ/mol", 0.2, 0.3, 10.0, 18.0, Mutability::Fixed, false));
    runTest(&ff, false);
}

TEST_F (ForceFieldParameterListTest, ModifyParameter) {
    ForceFieldParameterList ff("Improper");
    ff.addParameter("h3", ForceFieldParameter("sigma", "nm", 12.0, 0.3, 10.0, 18.0, Mutability::Free, false));
    ff.addParameter("c2", ForceFieldParameter("gamma", "", 11.0, 0.25, 8.0, 15.0, Mutability::Fixed, true));
    ff.addParameter("h3", ForceFieldParameter("epsilon", "kJ/mol", 0.2, 0.3, 10.0, 18.0, Mutability::Bounded, false));
    runTest(&ff, true);
}

}

}
