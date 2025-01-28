/*
 * This source file is part of the Alexandria program.
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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include "../msg_handler.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

class MsgHandlerTest : public gmx::test::CommandLineTestBase
{
    protected:
    gmx::test::TestReferenceChecker checker_;
    MsgHandler mh;
    MsgHandlerTest () : checker_(this->rootChecker())
    {
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    static void TearDownTestCase()
    {
    }

    void listMessages()
    {
        int i = 0;
        for(const auto &actm : ACTMessages)
        {
            checker_.checkString(actMessage(actm.first), std::to_string(i).c_str());
            i += 1;
        }
    }
        
    void runTest(bool verbose)
    {
        mh.setFilePointer(stdout);
        mh.setVerbosity(verbose);
        int i = 0;
        for(const auto &actm : ACTMessages)
        {
            mh.warning(actm.first, std::to_string(i));
            i += 1;
        }
        for(const auto &actm : ACTMessages)
        {
            checker_.checkInt64(mh.warningCount(actm.first), actMessage(actm.first));
        }        
    }
};

TEST_F(MsgHandlerTest, List)
{
    listMessages();
}

TEST_F(MsgHandlerTest, Verbose) 
{
    runTest(true);
}

TEST_F(MsgHandlerTest, NonVerbose) 
{
    runTest(false);
}

} // namespace

} // namespace
