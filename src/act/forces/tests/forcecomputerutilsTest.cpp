/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements test of bonded force routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_listed-forces
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
