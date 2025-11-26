/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
/*! \internal \file
 * \brief
 * Implements test of bonded force routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_forces_tests
 */
#include "actpre.h"

#include "../forcecomputerutils.h"

#include <cmath>

#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
//#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! Number of atoms used in these tests.
#define NATOMS 4

class ForceComputerUtilsTest : public ::testing::Test
{
protected:
    gmx::test::TestReferenceData    refData_;
    gmx::test::TestReferenceChecker checker_;
    ForceComputerUtilsTest( ) :
        checker_(refData_.rootChecker())
    {
        gmx::test::FloatingPointTolerance tolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-6));
        checker_.setDefaultTolerance(tolerance);
    }

    void testBondAngle(const std::vector<gmx::RVec> &x)
    {
        rvec  r_ij, r_kj;
        real  cosine_angle, angle;
        
        angle = bond_angle(x[0], x[1], x[2],
                           r_ij, r_kj, &cosine_angle);
        checker_.checkReal(angle, "angle");
        checker_.checkReal(cosine_angle, "cosine_angle");
    }
    
    void testDihedralAngle(const std::vector<gmx::RVec> &x)
    {
        rvec  r_ij, r_kj, r_kl, m, n;
        real  angle;

        angle = dih_angle(x[0], x[1], x[2], x[3],
                          r_ij, r_kj, r_kl, m, n);
        
        checker_.checkReal(angle, "angle");
    }
};
    
TEST_F (ForceComputerUtilsTest, NinetyDegreesAngle)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0 },
        { 0, 0, 1 },
        { 0, 1, 1 }
    };
    testBondAngle(x);
}

TEST_F (ForceComputerUtilsTest, LinearAngle)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0 },
        { 0, 0, 1 },
        { 0, 0, 2 }
    };
    testBondAngle(x);
}

TEST_F (ForceComputerUtilsTest, ZeroAngle)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0 },
        { 0, 0, 1 },
        { 0, 0, 0 }
    };
    testBondAngle(x);
}

TEST_F (ForceComputerUtilsTest, DihedralAnglePbcNone)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0 },
        { 0, 0, 1 },
        { 0, 1, 1 },
        { 1, 1, 1 }
    };
    testDihedralAngle(x);
}

}  // namespace

}  // namespace alexandria
