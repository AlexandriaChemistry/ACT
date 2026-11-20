/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2025,2026
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
#include <string>

#include "../msg_handler.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"

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
            const auto aaa = actMessage(actm.first);
            if (!aaa.empty())
            {
                checker_.checkString(aaa, std::to_string(i).c_str());
            }
            i += 1;
        }
    }
        
    void runTest(ACTStatus level)
    {
        mh.setPrintLevel(level);
        int i = 0;
        for(const auto &actm : ACTMessages)
        {
            mh.msg(level, actm.first, std::to_string(i));
            i += 1;
        }
        for(const auto &actm : ACTMessages)
        {
            checker_.checkInt64(mh.warningCount(actm.first), actMessage(actm.first).c_str());
        }        
    }
};

TEST_F(MsgHandlerTest, List)
{
    listMessages();
}

TEST_F(MsgHandlerTest, Verbose) 
{
    runTest(ACTStatus::Warning);
}

TEST_F(MsgHandlerTest, NonVerbose) 
{
    runTest(ACTStatus::Error);
}

// Test: Messages can be kept, for programs that report them in another
// way than in the log file, from a level of seriousness on
TEST(MsgHandlerRecordTest, Record)
{
    MsgHandler                      mh;
    std::vector<ACTRecordedMessage> kept;
    // Nothing is printed
    mh.setPrintLevel(ACTStatus::Fatal);
    mh.record(&kept, ACTStatus::Warning);
    mh.msg(ACTStatus::Error, ACTMessage::NoInChI, "first");
    mh.msg(ACTStatus::Warning, "second");
    mh.msg(ACTStatus::Info, "not serious enough");
    ASSERT_EQ(2u, kept.size());
    EXPECT_EQ(ACTStatus::Error, kept[0].level);
    EXPECT_EQ(ACTMessage::NoInChI, kept[0].actm);
    EXPECT_EQ("first", kept[0].text);
    EXPECT_EQ(ACTStatus::Warning, kept[1].level);
    EXPECT_EQ(ACTMessage::Silent, kept[1].actm);
    EXPECT_EQ("second", kept[1].text);
    // No longer
    mh.record(nullptr);
    mh.msg(ACTStatus::Error, "after");
    EXPECT_EQ(2u, kept.size());
    EXPECT_STREQ("Error", actStatusName(ACTStatus::Error));
    EXPECT_STREQ("Warning", actStatusName(ACTStatus::Warning));
}

TEST_F(MsgHandlerTest, Fatal) 
{
    EXPECT_THROW(runTest(ACTStatus::Fatal), gmx::InvalidInputError);
}

} // namespace

} // namespace
