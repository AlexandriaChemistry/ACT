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

#include "alexandria/act/forcefieldparameter.h"

#include <map>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class ForceFieldParameterTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;

        ForceFieldParameterTest () : checker_(this->rootChecker())
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

        void runTest(const ForceFieldParameterMap &fff)
        {
            for(auto &ff : fff)
            {
                checker_.checkString(ff.first, "type");
                checker_.checkString(ff.second.unit(), "unit");
                checker_.checkDouble(ff.second.value(), "value");
                checker_.checkDouble(ff.second.originalValue(), "originalValue");
                checker_.checkDouble(ff.second.minimum(), "minimum");
                checker_.checkDouble(ff.second.maximum(), "maximum");
                checker_.checkDouble(ff.second.uncertainty(), "uncertainty");
                checker_.checkDouble(ff.second.originalUncertainty(), "originalUncertainty");
                checker_.checkUInt64(ff.second.ntrain(), "ntrain");
                checker_.checkString(mutabilityName(ff.second.mutability()), "mutability");
                checker_.checkBoolean(ff.second.strict(), "strict");
            }
        }

};

TEST_F (ForceFieldParameterTest, FreeNotStrict) {
    ForceFieldParameter    fp("nm", 12.0, 0.3, 13, 10.0, 18.0, Mutability::Free, false);
    ForceFieldParameterMap ff = { { "sigma", fp } } ;
    runTest(ff);
}

TEST_F (ForceFieldParameterTest, BoundedStrict) {
    ForceFieldParameter    fp("kJ/mol", 11.0, 0.25, 24, 8.0, 15.0, Mutability::Bounded, true);
    ForceFieldParameterMap ff = { { "epsilon", fp } };
    runTest(ff);
    EXPECT_THROW(fp.setValue(17.0), gmx::InvalidInputError);
}

TEST_F (ForceFieldParameterTest, FixedStrict) {
    ForceFieldParameter    fp("", 11.0, 0.25, 57, 8.0, 15.0, Mutability::Fixed, true);
    ForceFieldParameterMap ff = { { "gamma", fp } };
    runTest(ff);
    EXPECT_THROW(fp.setValue(13.0), gmx::InvalidInputError);
}

TEST_F (ForceFieldParameterTest, BoundedNotStrict) {
    ForceFieldParameter    fp("1/nm", 11.0, 0.25, 913, 8.0, 15.0, Mutability::Bounded, false);
    ForceFieldParameterMap ff = { { "beta", fp } };
    runTest(ff);
    fp.setValue(17.0);
}

TEST_F (ForceFieldParameterTest, BoundedNotStrictOriginal) {
    ForceFieldParameter    fp("1/nm", 11.0, 0.25, 35, 8.0, 15.0, Mutability::Bounded, false);
    ForceFieldParameterMap ff = { {"beta", fp } };
    fp.setValue(17.0);
    runTest(ff);
}

TEST_F (ForceFieldParameterTest, FixedNotStrict) {
    ForceFieldParameter    fp("nm2", 11.0, 0.25, 3, 8.0, 15.0, Mutability::Fixed, false);
    ForceFieldParameterMap ff = { { "kappa", fp } };
    runTest(ff);
    fp.setValue(13.0);
}

TEST_F (ForceFieldParameterTest, FixedNotStrictOriginal) {
    ForceFieldParameter    fp("nm2", 11.0, 0.25, 1, 8.0, 15.0, Mutability::Fixed, false);
    ForceFieldParameterMap ff = { { "kappa", fp } };
    fp.setValue(13.0);
    runTest(ff);
}

TEST_F (ForceFieldParameterTest, UncertaintyNotStrict) {
    ForceFieldParameter    fp("ps", 11.0, 0.25, 17, 8.0, 15.0, Mutability::Fixed, false);
    ForceFieldParameterMap ff = { { "ypsilon", fp } };
    fp.setUncertainty(1.0);
}

TEST_F (ForceFieldParameterTest, UncertaintyNotStrictOriginal) {
    ForceFieldParameter    fp("fs", 11.0, 0.25, 123957, 8.0, 15.0, Mutability::Free, false);
    ForceFieldParameterMap ff = { { "ypsilon", fp } };
    fp.setUncertainty(1.0);
    runTest(ff);
}

TEST(ForceFieldParameterSimpleTest, UncertaintyStrict) {
    ForceFieldParameter    fp("kJ/mol", 11.0, 0.25, 45, 8.0, 15.0, Mutability::Fixed, true);
    ForceFieldParameterMap ff = { { "zeta", fp } };
    EXPECT_THROW(fp.setUncertainty(3.0), gmx::InternalError);
}

TEST_F (ForceFieldParameterTest, NtrainNotStrict) {
    ForceFieldParameter    fp("ps", 11.0, 0.25, 14, 8.0, 15.0, Mutability::Fixed, false);
    ForceFieldParameterMap ff = { { "ypsilon", fp } };
    fp.setNtrain(173);
    runTest(ff);
}

TEST(ForceFieldParameterSimpleTest, NtrainStrict) {
    ForceFieldParameter    fp("kJ/mol", 11.0, 0.25, 45, 8.0, 15.0, Mutability::Fixed, true);
    ForceFieldParameterMap ff = { { "zeta", fp } };
    EXPECT_THROW(fp.setNtrain(30), gmx::InternalError);
}

TEST(ForceFieldParameterSimpleTest, NameToMutability) {
    Mutability m;
    EXPECT_TRUE(nameToMutability("Free", &m));
    EXPECT_TRUE(nameToMutability("Fixed", &m));
    EXPECT_TRUE(nameToMutability("Bounded", &m));
    EXPECT_FALSE(nameToMutability("BounDed", &m));
    EXPECT_FALSE(nameToMutability("Bread", &m));
    EXPECT_FALSE(nameToMutability("spples", &m));
}

}

}
