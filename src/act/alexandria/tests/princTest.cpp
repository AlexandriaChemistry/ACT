/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024-2026
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
 * Tests for principal component and coordinate geometry functions.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <cmath>

#include <gtest/gtest.h>

#include "act/alexandria/princ.h"
#include "act/alexandria/topology.h"
#include "gromacs/math/vec.h"

#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

class PrincTest : public ::testing::Test
{
protected:
    std::vector<gmx::RVec> coords_;
    std::vector<ActAtom>   atoms_;
    std::vector<int>       index_;

    PrincTest()
    {
        // Set up a simple water-like molecule:
        // O at (0,0,0), H at (0.1,0,0), H at (-0.03,0.095,0)
        coords_ = {
            { 0.0, 0.0, 0.0 },
            { 0.1, 0.0, 0.0 },
            { -0.03, 0.095, 0.0 }
        };
        atoms_.emplace_back("O", "O", "OW", ActParticle::Atom, 8, 16.0, -0.82, 2);
        atoms_.emplace_back("H1", "H", "HW", ActParticle::Atom, 1, 1.008, 0.41, 1);
        atoms_.emplace_back("H2", "H", "HW", ActParticle::Atom, 1, 1.008, 0.41, 1);
        index_ = { 0, 1, 2 };
    }
};

TEST_F(PrincTest, CalcXcmMassWeighted)
{
    gmx::RVec xcm;
    real tm = calc_xcm(coords_, index_, atoms_, &xcm, false);
    // Total mass = 16.0 + 1.008 + 1.008 = 18.016
    EXPECT_NEAR(tm, 18.016, 1e-6);
    // Center of mass should be close to the oxygen atom (heaviest)
    // xcm_x = (16*0 + 1.008*0.1 + 1.008*(-0.03))/18.016
    double expected_x = (16.0*0.0 + 1.008*0.1 + 1.008*(-0.03)) / 18.016;
    double expected_y = (16.0*0.0 + 1.008*0.0 + 1.008*0.095) / 18.016;
    EXPECT_NEAR(xcm[XX], expected_x, 1e-6);
    EXPECT_NEAR(xcm[YY], expected_y, 1e-6);
    EXPECT_NEAR(xcm[ZZ], 0.0, 1e-6);
}

TEST_F(PrincTest, CalcXcmChargeWeighted)
{
    gmx::RVec xcm;
    real tm = calc_xcm(coords_, index_, atoms_, &xcm, true);
    // Total absolute charge = |−0.82| + |0.41| + |0.41| = 1.64
    EXPECT_NEAR(tm, 1.64, 1e-6);
    // Charge-weighted center
    double expected_x = (0.82*0.0 + 0.41*0.1 + 0.41*(-0.03)) / 1.64;
    double expected_y = (0.82*0.0 + 0.41*0.0 + 0.41*0.095) / 1.64;
    EXPECT_NEAR(xcm[XX], expected_x, 1e-6);
    EXPECT_NEAR(xcm[YY], expected_y, 1e-6);
    EXPECT_NEAR(xcm[ZZ], 0.0, 1e-6);
}

TEST_F(PrincTest, SubXcmAndAddXcmRoundTrip)
{
    // Save original coordinates
    std::vector<gmx::RVec> orig = coords_;

    // Subtract center of mass
    gmx::RVec xcm;
    sub_xcm(&coords_, index_, atoms_, &xcm, false);

    // After subtracting, center of mass should be at origin
    gmx::RVec xcm_check;
    calc_xcm(coords_, index_, atoms_, &xcm_check, false);
    EXPECT_NEAR(xcm_check[XX], 0.0, 1e-10);
    EXPECT_NEAR(xcm_check[YY], 0.0, 1e-10);
    EXPECT_NEAR(xcm_check[ZZ], 0.0, 1e-10);

    // Add it back
    add_xcm(&coords_, index_, xcm);

    // Should match original coordinates
    for (size_t i = 0; i < coords_.size(); i++)
    {
        EXPECT_NEAR(coords_[i][XX], orig[i][XX], 1e-10);
        EXPECT_NEAR(coords_[i][YY], orig[i][YY], 1e-10);
        EXPECT_NEAR(coords_[i][ZZ], orig[i][ZZ], 1e-10);
    }
}

TEST_F(PrincTest, RotateAtomsIdentity)
{
    // Identity matrix should not change coordinates
    matrix identity;
    identity[XX][XX] = 1; identity[XX][YY] = 0; identity[XX][ZZ] = 0;
    identity[YY][XX] = 0; identity[YY][YY] = 1; identity[YY][ZZ] = 0;
    identity[ZZ][XX] = 0; identity[ZZ][YY] = 0; identity[ZZ][ZZ] = 1;

    std::vector<gmx::RVec> orig = coords_;
    rotate_atoms(index_, &coords_, identity);

    for (size_t i = 0; i < coords_.size(); i++)
    {
        EXPECT_NEAR(coords_[i][XX], orig[i][XX], 1e-10);
        EXPECT_NEAR(coords_[i][YY], orig[i][YY], 1e-10);
        EXPECT_NEAR(coords_[i][ZZ], orig[i][ZZ], 1e-10);
    }
}

TEST_F(PrincTest, RotateAtoms90DegreesAroundZ)
{
    // 90 degree rotation around Z axis: x' = -y, y' = x
    matrix rot;
    rot[XX][XX] =  0; rot[XX][YY] = -1; rot[XX][ZZ] = 0;
    rot[YY][XX] =  1; rot[YY][YY] =  0; rot[YY][ZZ] = 0;
    rot[ZZ][XX] =  0; rot[ZZ][YY] =  0; rot[ZZ][ZZ] = 1;

    std::vector<gmx::RVec> orig = coords_;
    rotate_atoms(index_, &coords_, rot);

    for (size_t i = 0; i < index_.size(); i++)
    {
        int ii = index_[i];
        EXPECT_NEAR(coords_[ii][XX], -orig[ii][YY], 1e-10);
        EXPECT_NEAR(coords_[ii][YY],  orig[ii][XX], 1e-10);
        EXPECT_NEAR(coords_[ii][ZZ],  orig[ii][ZZ], 1e-10);
    }
}

TEST_F(PrincTest, PrincipalCompInertia)
{
    // Center atoms at origin first
    gmx::RVec xcm;
    sub_xcm(&coords_, index_, atoms_, &xcm, false);

    std::vector<real> masses;
    for (const auto &a : atoms_)
    {
        masses.push_back(a.mass());
    }

    matrix    trans;
    gmx::RVec inertia;
    principal_comp(index_, masses, coords_, &trans, &inertia);

    // For a planar molecule, the largest moment of inertia
    // should be around the axis perpendicular to the plane (Z)
    // Moments should be sorted in ascending order
    EXPECT_TRUE(inertia[0] <= inertia[1]);
    EXPECT_TRUE(inertia[1] <= inertia[2]);

    // All moments should be non-negative
    EXPECT_GE(inertia[0], 0.0);
    EXPECT_GE(inertia[1], 0.0);
    EXPECT_GE(inertia[2], 0.0);
}

TEST_F(PrincTest, CalcXcmEmptyIndex)
{
    // With empty index, all atoms are used
    std::vector<int> emptyIndex;
    gmx::RVec xcm1, xcm2;

    real tm1 = calc_xcm(coords_, emptyIndex, atoms_, &xcm1, false);
    real tm2 = calc_xcm(coords_, index_,      atoms_, &xcm2, false);

    // Results should be the same since index_ contains all atoms
    EXPECT_NEAR(tm1, tm2, 1e-10);
    EXPECT_NEAR(xcm1[XX], xcm2[XX], 1e-10);
    EXPECT_NEAR(xcm1[YY], xcm2[YY], 1e-10);
    EXPECT_NEAR(xcm1[ZZ], xcm2[ZZ], 1e-10);
}

} // namespace

} // namespace alexandria
