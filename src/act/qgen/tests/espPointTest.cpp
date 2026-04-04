/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2025
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
 * Tests electrostatic potential point representation and RESP charge
 * generation initialization.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <gtest/gtest.h>

#include "act/qgen/qgen_resp.h"

#include "gromacs/math/vectypes.h"

#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

TEST(EspPointTest, ConstructorAndGetters)
{
    gmx::RVec coords = { 1.0, 2.0, 3.0 };
    double     potential = 0.5;
    EspPoint   ep(coords, potential);

    EXPECT_DOUBLE_EQ(1.0, ep.esp()[XX]);
    EXPECT_DOUBLE_EQ(2.0, ep.esp()[YY]);
    EXPECT_DOUBLE_EQ(3.0, ep.esp()[ZZ]);
    EXPECT_DOUBLE_EQ(0.5, ep.v());
}

TEST(EspPointTest, DefaultVCalcIsZero)
{
    gmx::RVec coords = { 0.0, 0.0, 0.0 };
    EspPoint   ep(coords, 1.0);
    EXPECT_DOUBLE_EQ(0.0, ep.vCalc());
}

TEST(EspPointTest, DefaultRhoIsZero)
{
    gmx::RVec coords = { 0.0, 0.0, 0.0 };
    EspPoint   ep(coords, 1.0);
    EXPECT_DOUBLE_EQ(0.0, ep.rho());
}

TEST(EspPointTest, SetV)
{
    gmx::RVec coords = { 0.0, 0.0, 0.0 };
    EspPoint   ep(coords, 1.0);
    EXPECT_DOUBLE_EQ(1.0, ep.v());
    ep.setV(42.0);
    EXPECT_DOUBLE_EQ(42.0, ep.v());
}

TEST(EspPointTest, SetVCalc)
{
    gmx::RVec coords = { 0.0, 0.0, 0.0 };
    EspPoint   ep(coords, 1.0);
    ep.setVCalc(3.14);
    EXPECT_DOUBLE_EQ(3.14, ep.vCalc());
}

TEST(EspPointTest, SetRho)
{
    gmx::RVec coords = { 0.0, 0.0, 0.0 };
    EspPoint   ep(coords, 0.0);
    ep.setRho(2.71828);
    EXPECT_DOUBLE_EQ(2.71828, ep.rho());
}

TEST(EspPointTest, NegativeValues)
{
    gmx::RVec coords = { -1.5, -2.5, -3.5 };
    EspPoint   ep(coords, -0.123);

    EXPECT_DOUBLE_EQ(-1.5, ep.esp()[XX]);
    EXPECT_DOUBLE_EQ(-2.5, ep.esp()[YY]);
    EXPECT_DOUBLE_EQ(-3.5, ep.esp()[ZZ]);
    EXPECT_DOUBLE_EQ(-0.123, ep.v());

    ep.setVCalc(-0.456);
    EXPECT_DOUBLE_EQ(-0.456, ep.vCalc());

    ep.setRho(-0.789);
    EXPECT_DOUBLE_EQ(-0.789, ep.rho());
}

TEST(EspPointTest, MultipleSetCalls)
{
    gmx::RVec coords = { 0.0, 0.0, 0.0 };
    EspPoint   ep(coords, 0.0);

    // Set and override multiple times
    ep.setV(1.0);
    ep.setV(2.0);
    EXPECT_DOUBLE_EQ(2.0, ep.v());

    ep.setVCalc(10.0);
    ep.setVCalc(20.0);
    EXPECT_DOUBLE_EQ(20.0, ep.vCalc());

    ep.setRho(100.0);
    ep.setRho(200.0);
    EXPECT_DOUBLE_EQ(200.0, ep.rho());
}

// ============================================================
// Tests for QgenResp basic construction and properties
// ============================================================

TEST(QgenRespTest, DefaultConstruction)
{
    QgenResp resp;
    EXPECT_EQ(ChargeDistributionType::Point, resp.chargeType());
    EXPECT_DOUBLE_EQ(0.0, resp.getMolecularCharge());
    EXPECT_EQ(0u, resp.nEsp());
    EXPECT_EQ(0u, resp.natoms());
    EXPECT_TRUE(resp.espPoints().empty());
    EXPECT_TRUE(resp.getStoichiometry().empty());
}

TEST(QgenRespTest, SetChargeDistributionType)
{
    QgenResp resp;
    resp.setChargeDistributionType(ChargeDistributionType::Gaussian);
    EXPECT_EQ(ChargeDistributionType::Gaussian, resp.chargeType());

    resp.setChargeDistributionType(ChargeDistributionType::Slater);
    EXPECT_EQ(ChargeDistributionType::Slater, resp.chargeType());

    resp.setChargeDistributionType(ChargeDistributionType::Point);
    EXPECT_EQ(ChargeDistributionType::Point, resp.chargeType());
}

TEST(QgenRespTest, AddEspPoints)
{
    QgenResp resp;
    EXPECT_EQ(0u, resp.nEsp());

    resp.addEspPoint(1.0, 2.0, 3.0, 0.5);
    EXPECT_EQ(1u, resp.nEsp());

    resp.addEspPoint(4.0, 5.0, 6.0, -0.3);
    EXPECT_EQ(2u, resp.nEsp());

    // Verify first point
    const auto &ep0 = resp.espPoint(0);
    EXPECT_DOUBLE_EQ(1.0, ep0.esp()[XX]);
    EXPECT_DOUBLE_EQ(2.0, ep0.esp()[YY]);
    EXPECT_DOUBLE_EQ(3.0, ep0.esp()[ZZ]);
    EXPECT_DOUBLE_EQ(0.5, ep0.v());

    // Verify second point
    const auto &ep1 = resp.espPoint(1);
    EXPECT_DOUBLE_EQ(4.0, ep1.esp()[XX]);
    EXPECT_DOUBLE_EQ(5.0, ep1.esp()[YY]);
    EXPECT_DOUBLE_EQ(6.0, ep1.esp()[ZZ]);
    EXPECT_DOUBLE_EQ(-0.3, ep1.v());
}

TEST(QgenRespTest, EspPointsReference)
{
    QgenResp resp;
    resp.addEspPoint(1.0, 2.0, 3.0, 10.0);
    resp.addEspPoint(4.0, 5.0, 6.0, 20.0);

    const auto &points = resp.espPoints();
    EXPECT_EQ(2u, points.size());
    EXPECT_DOUBLE_EQ(10.0, points[0].v());
    EXPECT_DOUBLE_EQ(20.0, points[1].v());
}

}

}
