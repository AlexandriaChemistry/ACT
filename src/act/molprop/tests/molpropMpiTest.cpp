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
 * Tests of sending molecular properties to other MPI ranks.
 */

#include "actpre.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/basics/msg_handler.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class MolpropMpiTest : public ::testing::Test
{
protected:
    //! For the files
    gmx::test::TestFileManager fileManager_;
    //! The communication between the ranks
    CommunicationRecord        cr_;
    //! The molecules as read by each rank
    std::vector<MolProp>       ref_;

    MolpropMpiTest()
    {
        cr_.init(cr_.size());
        MsgHandler msghandler;
        const auto fileName = fileManager_.getInputFilePath("molprop/molprop.xml");
        MolPropRead(&msghandler, fileName.c_str(), &ref_);
        // Alcohols with frequencies and intensities
        std::vector<MolProp> booze;
        const auto boozeName = fileManager_.getInputFilePath("molprop/booze-mono.xml");
        MolPropRead(&msghandler, boozeName.c_str(), &booze);
        for (auto &mp : booze)
        {
            ref_.push_back(std::move(mp));
        }
        // The molecules in the file have no categories
        for (auto &mp : ref_)
        {
            mp.AddCategory("Test");
            mp.AddCategory(mp.getMolname() + " category");
        }
    }

    //! \return the contents of the molecules written to a file of this rank
    std::string written(const std::vector<MolProp> &mps, const std::string &name)
    {
        const auto fileName = fileManager_.getTemporaryFilePath(gmx::formatString("%s%d.xml", name.c_str(), cr_.rank()));
        MolPropWrite(fileName.c_str(), mps, false);
        const std::ifstream in(fileName);
        std::stringstream ss;
        ss << in.rdbuf();
        return ss.str();
    }
};

// Test: The molecules the master broadcasts are the same on all ranks
TEST_F(MolpropMpiTest, BroadCast)
{
    // The master broadcasts the molecules it read, the others overwrite theirs
    std::vector<MolProp> mps(ref_.size());
    auto                &received = cr_.isMaster() ? ref_ : mps;
    for (auto &mp : received)
    {
        EXPECT_EQ(CommunicationStatus::OK, mp.BroadCast(&cr_, 0, cr_.comm_world()));
    }
    if (!cr_.isMaster())
    {
        EXPECT_EQ(written(ref_, "ref"), written(mps, "received"));
    }
}

} // namespace

} // namespace alexandria
