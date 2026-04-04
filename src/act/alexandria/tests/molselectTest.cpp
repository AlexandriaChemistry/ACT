/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024-2026
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
 * Tests for the MolSelect and IMolSelect classes.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <string>

#include <gtest/gtest.h>

#include "act/alexandria/molselect.h"
#include "act/basics/dataset.h"

#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

// ---- IMolSelect tests ----

TEST(MolSelectTest, IMolSelectConstructorAndGetters)
{
    IMolSelect ims("water", iMolSelect::Train, 5);
    EXPECT_EQ(ims.iupac(), "water");
    EXPECT_EQ(ims.status(), iMolSelect::Train);
    EXPECT_EQ(ims.index(), 5);
}

TEST(MolSelectTest, IMolSelectIsNotDimer)
{
    IMolSelect ims("ethanol", iMolSelect::Test, 0);
    EXPECT_FALSE(ims.isDimer());
}

TEST(MolSelectTest, IMolSelectIsDimer)
{
    // Dimers use '#' separator in IUPAC name
    IMolSelect ims("water#methanol", iMolSelect::Train, 1);
    EXPECT_TRUE(ims.isDimer());
}

// ---- MolSelect tests ----

TEST(MolSelectTest, MolSelectEmpty)
{
    MolSelect ms;
    EXPECT_EQ(ms.nMol(), 0u);
    EXPECT_EQ(ms.countDimers(), 0u);
}

TEST(MolSelectTest, MolSelectAddAndRetrieve)
{
    MolSelect ms;
    ms.addOne("water", 0, iMolSelect::Train);
    ms.addOne("ethanol", 1, iMolSelect::Test);
    ms.addOne("methane", 2, iMolSelect::Ignore);

    EXPECT_EQ(ms.nMol(), 3u);
}

TEST(MolSelectTest, MolSelectStatusLookup)
{
    MolSelect ms;
    ms.addOne("water", 0, iMolSelect::Train);
    ms.addOne("ethanol", 1, iMolSelect::Test);

    iMolSelect status;
    EXPECT_TRUE(ms.status("water", &status));
    EXPECT_EQ(status, iMolSelect::Train);

    EXPECT_TRUE(ms.status("ethanol", &status));
    EXPECT_EQ(status, iMolSelect::Test);

    EXPECT_FALSE(ms.status("unknown", &status));
}

TEST(MolSelectTest, MolSelectIndexLookup)
{
    MolSelect ms;
    ms.addOne("water", 42, iMolSelect::Train);
    ms.addOne("methane", 7, iMolSelect::Test);

    int index;
    EXPECT_TRUE(ms.index("water", &index));
    EXPECT_EQ(index, 42);

    EXPECT_TRUE(ms.index("methane", &index));
    EXPECT_EQ(index, 7);

    EXPECT_FALSE(ms.index("unknown_mol", &index));
}

TEST(MolSelectTest, MolSelectCount)
{
    MolSelect ms;
    ms.addOne("water", 0, iMolSelect::Train);
    ms.addOne("ethanol", 1, iMolSelect::Train);
    ms.addOne("methane", 2, iMolSelect::Test);
    ms.addOne("benzene", 3, iMolSelect::Ignore);

    EXPECT_EQ(ms.count(iMolSelect::Train), 2);
    EXPECT_EQ(ms.count(iMolSelect::Test), 1);
    EXPECT_EQ(ms.count(iMolSelect::Ignore), 1);
}

TEST(MolSelectTest, MolSelectCountDimers)
{
    MolSelect ms;
    ms.addOne("water", 0, iMolSelect::Train);
    ms.addOne("water#methanol", 1, iMolSelect::Train);
    ms.addOne("ethanol#acetone", 2, iMolSelect::Test);
    ms.addOne("methane", 3, iMolSelect::Test);

    EXPECT_EQ(ms.countDimers(), 2u);
}

TEST(MolSelectTest, MolSelectDimerReverseLookup)
{
    MolSelect ms;
    ms.addOne("water#methanol", 0, iMolSelect::Train);

    // findIupac should find the dimer with reversed order too
    iMolSelect status;
    EXPECT_TRUE(ms.status("water#methanol", &status));
    EXPECT_EQ(status, iMolSelect::Train);

    // Try reversed order
    EXPECT_TRUE(ms.status("methanol#water", &status));
    EXPECT_EQ(status, iMolSelect::Train);
}

TEST(MolSelectTest, MolSelectImolSelectVector)
{
    MolSelect ms;
    ms.addOne("water", 0, iMolSelect::Train);
    ms.addOne("ethanol", 1, iMolSelect::Test);

    const auto &vec = ms.imolSelect();
    ASSERT_EQ(vec.size(), 2u);
    EXPECT_EQ(vec[0].iupac(), "water");
    EXPECT_EQ(vec[1].iupac(), "ethanol");
}

} // namespace

} // namespace alexandria
