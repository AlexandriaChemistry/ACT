/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2023
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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include "act/forcefield/act_checksum.h"
#include "act/forcefield/forcefield_utils.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "gromacs/topology/atoms.h"

namespace alexandria
{

namespace
{

class ActCheckSumTest : public gmx::test::CommandLineTestBase
{
protected:
    gmx::test::TestReferenceChecker checker_;
    
    ActCheckSumTest () : checker_(this->rootChecker())
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
    
    void TestOne(const std::string &fileName)
    {
        std::string dataName = gmx::test::TestFileManager::getInputFilePath(fileName+".xml");
        auto pd      = getForceField(fileName);
        auto orig    = pd->checkSum();
        auto sumpd   = forcefieldCheckSum(pd);
        EXPECT_TRUE(orig == sumpd);
        checker_.checkString(sumpd, "ForceField checksum");
    }
};

TEST_F (ActCheckSumTest, ACM_g) {
    TestOne("ACM-g");
}

TEST_F (ActCheckSumTest, ACM_pg) {
    TestOne("ACM-pg");
}

}

} // namespace
