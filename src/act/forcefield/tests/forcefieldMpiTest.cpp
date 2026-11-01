/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2026
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
 * Tests of sending a force field to other MPI ranks.
 */

#include "actpre.h"

#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "act/forcefield/act_checksum.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class ForceFieldMpiTest : public ::testing::Test
{
protected:
    //! For the files
    gmx::test::TestFileManager fileManager_;
    //! The communication between the ranks
    CommunicationRecord        cr_;
    //! The force field as read by each rank
    ForceField                 ref_;

    ForceFieldMpiTest()
    {
        cr_.init(cr_.size());
        readForceField(fileManager_.getInputFilePath("ACM-pg.xml"), &ref_);
    }

    //! \return the contents of the force field written to a file of this rank
    std::string written(ForceField *pd, const std::string &name)
    {
        auto fileName = fileManager_.getTemporaryFilePath(gmx::formatString("%s%d.xml", name.c_str(), cr_.rank()));
        writeForceField(fileName, pd, false);
        std::ifstream     in(fileName);
        std::stringstream ss;
        ss << in.rdbuf();
        return ss.str();
    }

    //! Check that a force field is the same as the one read
    void check(ForceField *pd)
    {
        EXPECT_EQ(forcefieldCheckSum(&ref_), forcefieldCheckSum(pd));
        EXPECT_EQ(ref_.polarizable(), pd->polarizable());
        EXPECT_EQ(written(&ref_, "ref"), written(pd, "received"));
    }
};

// Test: The force field the master sends to each rank is the same as the
// one read from the file
TEST_F(ForceFieldMpiTest, SendReceive)
{
    // With one rank there is nothing to send
    if (cr_.isMaster())
    {
        for (int dest = 1; dest < cr_.size(); dest++)
        {
            EXPECT_EQ(CommunicationStatus::OK, ref_.Send(&cr_, dest));
        }
    }
    else
    {
        ForceField pd;
        EXPECT_EQ(CommunicationStatus::OK, pd.Receive(&cr_, 0));
        check(&pd);
    }
}

// Test: The force field the master broadcasts is the same on all ranks
TEST_F(ForceFieldMpiTest, BroadCast)
{
    ForceField pd;
    if (cr_.isMaster())
    {
        pd = ref_;
    }
    EXPECT_EQ(CommunicationStatus::OK, pd.BroadCast(&cr_, 0, cr_.comm_world()));
    check(&pd);
}

} // namespace

} // namespace alexandria
