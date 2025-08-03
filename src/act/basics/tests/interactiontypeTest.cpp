/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2020-2025
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

#include "../interactiontype.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "gromacs/utility/exceptions.h"

namespace alexandria
{

namespace
{

TEST(InteractionTypeTest, stringToInteractionType) {
    InteractionType iType;
    EXPECT_TRUE(stringToInteractionType("BONDCORRECTIONS", &iType));
    EXPECT_TRUE(iType == InteractionType::BONDCORRECTIONS);
    EXPECT_TRUE(stringToInteractionType("BONDS", &iType));
    EXPECT_FALSE(iType == InteractionType::ANGLES);
    EXPECT_FALSE(stringToInteractionType("FOO", &iType));
}

}

}
