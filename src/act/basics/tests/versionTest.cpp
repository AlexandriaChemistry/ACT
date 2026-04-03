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
 * Tests for version.h/cpp (act_welcome and act_goodbye functions).
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include "../version.h"

#include "testutils/testasserts.h"

namespace
{

TEST(VersionTest, actWelcomeNotEmpty)
{
    EXPECT_FALSE(act_welcome().empty());
}

TEST(VersionTest, actGoodbyeNotEmpty)
{
    EXPECT_FALSE(act_goodbye().empty());
}

TEST(VersionTest, actWelcomeContainsAlexandria)
{
    EXPECT_NE(std::string::npos, act_welcome().find("Alexandria"));
}

TEST(VersionTest, actGoodbyeContainsAlexandria)
{
    EXPECT_NE(std::string::npos, act_goodbye().find("Alexandria"));
}

} // namespace
