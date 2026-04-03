/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
 * Tests for chargemodel.h/cpp (ChargeDistributionType and
 * ChargeGenerationAlgorithm enums and conversion functions).
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include "../chargemodel.h"

#include "testutils/testasserts.h"

#include "gromacs/utility/exceptions.h"

namespace alexandria
{

namespace
{

TEST(ChargeDistributionTypeTest, qdnamesHasThreeEntries)
{
    auto names = qdnames();
    EXPECT_EQ(3u, names.size());
}

TEST(ChargeDistributionTypeTest, qdnamesContainsExpectedNames)
{
    auto names = qdnames();
    EXPECT_EQ("Point",    names[0]);
    EXPECT_EQ("Gaussian", names[1]);
    EXPECT_EQ("Slater",   names[2]);
}

TEST(ChargeDistributionTypeTest, chargeDistributionTypeName)
{
    EXPECT_EQ("Point",    chargeDistributionTypeName(ChargeDistributionType::Point));
    EXPECT_EQ("Gaussian", chargeDistributionTypeName(ChargeDistributionType::Gaussian));
    EXPECT_EQ("Slater",   chargeDistributionTypeName(ChargeDistributionType::Slater));
}

TEST(ChargeDistributionTypeTest, name2ChargeDistributionTypeValid)
{
    EXPECT_EQ(ChargeDistributionType::Point,    name2ChargeDistributionType("Point"));
    EXPECT_EQ(ChargeDistributionType::Gaussian, name2ChargeDistributionType("Gaussian"));
    EXPECT_EQ(ChargeDistributionType::Slater,   name2ChargeDistributionType("Slater"));
}

TEST(ChargeDistributionTypeTest, name2ChargeDistributionTypeInvalidThrows)
{
    EXPECT_THROW(name2ChargeDistributionType("InvalidType"), gmx::InvalidInputError);
    // Input is case-sensitive
    EXPECT_THROW(name2ChargeDistributionType("point"),       gmx::InvalidInputError);
    EXPECT_THROW(name2ChargeDistributionType("gaussian"),    gmx::InvalidInputError);
    EXPECT_THROW(name2ChargeDistributionType(""),            gmx::InvalidInputError);
}

TEST(ChargeDistributionTypeTest, roundTrip)
{
    for (auto ct : { ChargeDistributionType::Point,
                     ChargeDistributionType::Gaussian,
                     ChargeDistributionType::Slater })
    {
        EXPECT_EQ(ct, name2ChargeDistributionType(chargeDistributionTypeName(ct)));
    }
}

TEST(ChargeGenerationAlgorithmTest, chargeGenerationAlgorithmName)
{
    EXPECT_EQ("None",   chargeGenerationAlgorithmName(ChargeGenerationAlgorithm::NONE));
    EXPECT_EQ("EEM",    chargeGenerationAlgorithmName(ChargeGenerationAlgorithm::EEM));
    EXPECT_EQ("SQE",    chargeGenerationAlgorithmName(ChargeGenerationAlgorithm::SQE));
    EXPECT_EQ("ESP",    chargeGenerationAlgorithmName(ChargeGenerationAlgorithm::ESP));
    EXPECT_EQ("Custom", chargeGenerationAlgorithmName(ChargeGenerationAlgorithm::Custom));
    EXPECT_EQ("Read",   chargeGenerationAlgorithmName(ChargeGenerationAlgorithm::Read));
}

TEST(ChargeGenerationAlgorithmTest, nameToChargeGenerationAlgorithmValid)
{
    EXPECT_EQ(ChargeGenerationAlgorithm::NONE,   nameToChargeGenerationAlgorithm("None"));
    EXPECT_EQ(ChargeGenerationAlgorithm::EEM,    nameToChargeGenerationAlgorithm("EEM"));
    EXPECT_EQ(ChargeGenerationAlgorithm::SQE,    nameToChargeGenerationAlgorithm("SQE"));
    EXPECT_EQ(ChargeGenerationAlgorithm::ESP,    nameToChargeGenerationAlgorithm("ESP"));
    EXPECT_EQ(ChargeGenerationAlgorithm::Custom, nameToChargeGenerationAlgorithm("Custom"));
    EXPECT_EQ(ChargeGenerationAlgorithm::Read,   nameToChargeGenerationAlgorithm("Read"));
}

TEST(ChargeGenerationAlgorithmTest, nameToChargeGenerationAlgorithmInvalidThrows)
{
    EXPECT_THROW(nameToChargeGenerationAlgorithm("InvalidAlgorithm"), gmx::InvalidInputError);
    // Input is case-sensitive
    EXPECT_THROW(nameToChargeGenerationAlgorithm("eem"),              gmx::InvalidInputError);
    EXPECT_THROW(nameToChargeGenerationAlgorithm("none"),             gmx::InvalidInputError);
    EXPECT_THROW(nameToChargeGenerationAlgorithm(""),                 gmx::InvalidInputError);
}

TEST(ChargeGenerationAlgorithmTest, roundTrip)
{
    for (auto cg : { ChargeGenerationAlgorithm::NONE,
                     ChargeGenerationAlgorithm::EEM,
                     ChargeGenerationAlgorithm::SQE,
                     ChargeGenerationAlgorithm::ESP,
                     ChargeGenerationAlgorithm::Custom,
                     ChargeGenerationAlgorithm::Read })
    {
        EXPECT_EQ(cg, nameToChargeGenerationAlgorithm(chargeGenerationAlgorithmName(cg)));
    }
}

} // namespace

} // namespace alexandria
