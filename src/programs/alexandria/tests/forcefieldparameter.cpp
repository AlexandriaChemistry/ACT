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

#include "programs/alexandria/forcefieldparameter.h"

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

        void runTest(const ForceFieldParameter &ff)
        {
            checker_.checkString(ff.identifier(), "identifier");
            checker_.checkString(ff.type(), "type");
            checker_.checkDouble(ff.value(), "value");
            checker_.checkDouble(ff.originalValue(), "originalValue");
            checker_.checkDouble(ff.minimum(), "minimum");
            checker_.checkDouble(ff.maximum(), "maximum");
            checker_.checkDouble(ff.uncertainty(), "uncertainty");
            checker_.checkDouble(ff.originalUncertainty(), "originalUncertainty");
            checker_.checkString(mutabilityName(ff.mutability()), "mutability");
            checker_.checkBoolean(ff.strict(), "strict");
        }

};

TEST_F (ForceFieldParameterTest, FreeNotStrict) {
    ForceFieldParameter ff("h3", "sigma", 12.0, 0.3, 10.0, 18.0, Mutability::Free, false);
    runTest(ff);
}

TEST_F (ForceFieldParameterTest, BoundedStrict) {
    ForceFieldParameter ff("c2", "epsilon", 11.0, 0.25, 8.0, 15.0, Mutability::Bounded, true);
    runTest(ff);
    EXPECT_THROW(ff.setValue(17.0), gmx::InvalidInputError);
}

TEST_F (ForceFieldParameterTest, FixedStrict) {
    ForceFieldParameter ff("c2", "gamma", 11.0, 0.25, 8.0, 15.0, Mutability::Fixed, true);
    runTest(ff);
    EXPECT_THROW(ff.setValue(13.0), gmx::InvalidInputError);
}

TEST_F (ForceFieldParameterTest, BoundedNotStrict) {
    ForceFieldParameter ff("c2", "beta", 11.0, 0.25, 8.0, 15.0, Mutability::Bounded, false);
    runTest(ff);
    ff.setValue(17.0);
}

TEST_F (ForceFieldParameterTest, BoundedNotStrictOriginal) {
    ForceFieldParameter ff("c2", "beta", 11.0, 0.25, 8.0, 15.0, Mutability::Bounded, false);
    ff.setValue(17.0);
    runTest(ff);
}

TEST_F (ForceFieldParameterTest, FixedNotStrict) {
    ForceFieldParameter ff("c2", "kappa", 11.0, 0.25, 8.0, 15.0, Mutability::Fixed, false);
    runTest(ff);
    ff.setValue(13.0);
}

TEST_F (ForceFieldParameterTest, FixedNotStrictOriginal) {
    ForceFieldParameter ff("c2", "kappa", 11.0, 0.25, 8.0, 15.0, Mutability::Fixed, false);
    ff.setValue(13.0);
    runTest(ff);
}

TEST_F (ForceFieldParameterTest, UncertaintyNotStrict) {
    ForceFieldParameter ff("c2", "ypsilon", 11.0, 0.25, 8.0, 15.0, Mutability::Fixed, false);
    ff.setUncertainty(1.0);
}

TEST_F (ForceFieldParameterTest, UncertaintyNotStrictOriginl) {
    ForceFieldParameter ff("c2", "ypsilon", 11.0, 0.25, 8.0, 15.0, Mutability::Fixed, false);
    ff.setUncertainty(1.0);
    runTest(ff);
}

TEST_F (ForceFieldParameterTest, UncertaintyStrict) {
    ForceFieldParameter ff("c2", "zeta", 11.0, 0.25, 8.0, 15.0, Mutability::Fixed, true);
    EXPECT_THROW(ff.setUncertainty(3.0), gmx::InternalError);
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
