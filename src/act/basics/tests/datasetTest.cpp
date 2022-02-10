/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022
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

#include "../dataset.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

class DataSetTest : public gmx::test::CommandLineTestBase
{
    protected:
    gmx::test::TestReferenceChecker checker_;
    DataSetTest () : checker_(this->rootChecker())
    {
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    static void TearDownTestCase()
    {
    }
    
    void runTest()
    {
        std::vector<std::string> names;
        for(auto &imsn : iMolSelectNames())
        {
            names.push_back(iMolSelectName(imsn.first));
        }
        checker_.checkSequence(names.begin(), names.end(), "iMolSelectNames");
    }
};

TEST(DataSetSimpleTest, WrongName) 
{
    iMolSelect ims;
    EXPECT_FALSE(name2molselect("koko", &ims));
}

TEST(DataSetSimpleTest, CorrectNames) 
{
    for(auto &imsn : iMolSelectNames())
    {
        iMolSelect ims;
        EXPECT_TRUE(name2molselect(imsn.second, &ims));
    }
}

TEST_F(DataSetTest, CheckNames) 
{
    runTest();
}

}

} // namespace
