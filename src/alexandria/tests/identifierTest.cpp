/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2020
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

#include "alexandria/identifier.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

TEST(IdentifierSimpleTest, ConstructorOne) {
    EXPECT_THROW(Identifier id({"C", "H"}, { 1 }, CanSwap::Yes), gmx::InvalidInputError);
}

TEST(IdentifierSimpleTest, BSmallerThanA) {
    Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    Identifier b({"C", "H"}, { 1.0 }, CanSwap::Yes);
    bool compare = b < a;
    EXPECT_TRUE((compare));
}

TEST(IdentifierSimpleTest, CompareDifferentLength) {
    Identifier a({"C", "H"}, { 1.0 }, CanSwap::Yes);
    Identifier b({"C", "O", "H"}, { 1.0, 1.0 }, CanSwap::Yes);
    bool compare = b < a;
    EXPECT_FALSE(compare);
}

TEST(IdentifierSimpleTest, BSmallerThanASecondAtom) {
    Identifier a({"H", "P"}, { 1.0 }, CanSwap::No);
    Identifier b({"H", "C"}, { 1.0 }, CanSwap::No);
    bool compare = a < b;
    EXPECT_FALSE(compare);
}

TEST(IdentifierSimpleTest, BSmallerThanABondOrder) {
    Identifier a({"H", "C"}, { 2 }, CanSwap::No);
    Identifier b({"H", "C"}, { 1 }, CanSwap::No);
    bool compare = b < a;
    EXPECT_TRUE(compare);
}

TEST(IdentifierSimpleTest, NotASmallerThanB) {
    Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    Identifier b({"C", "H"}, { 1.0 }, CanSwap::Yes);
    bool compare = a < b;
    EXPECT_FALSE(compare);
}

TEST(IdentifierSimpleTest, AEqualToB) {
    Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    Identifier b({"P", "H"}, { 1.0 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE(equal);
}

TEST(IdentifierSimpleTest, SwappedAEqualToB) {
    Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    Identifier b({"H", "P"}, { 1.0 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE(equal);
}

TEST(IdentifierSimpleTest, ANotEqualToB) {
    Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    Identifier b({"C", "H"}, { 1.0 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_FALSE(equal);
}

TEST(IdentifierSimpleTest, BSmallerThanABondOrder1) {
    Identifier a({"P", "H"}, { 1 }, CanSwap::Yes);
    Identifier b({"C", "H"}, { 1 }, CanSwap::Yes);
    bool compare = b < a;
    EXPECT_TRUE((compare));
}

TEST(IdentifierSimpleTest, NotASmallerThanBBondOrder2) {
    Identifier a({"P", "H"}, { 2 }, CanSwap::Yes);
    Identifier b({"C", "H"}, { 2 }, CanSwap::Yes);
    bool compare = a < b;
    EXPECT_FALSE((compare));
}

TEST(IdentifierSimpleTest, AEqualToBBondOrder1) {
    Identifier a({"P", "H"}, { 1 }, CanSwap::Yes);
    Identifier b({"P", "H"}, { 1 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, AEqualToBBondOrderDiff) {
    Identifier a({"P", "H"}, { 1 }, CanSwap::Yes);
    Identifier b({"P", "H"}, { 3 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_FALSE((equal));
}

TEST(IdentifierSimpleTest, SwappedAEqualToBBondOrder2) {
    Identifier a({"P", "H"}, { 2 }, CanSwap::Yes);
    Identifier b({"H", "P"}, { 2 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, EqualNotSwappableAB) {
    Identifier a({"P", "H"}, { 1.0 }, CanSwap::No);
    Identifier b({"H", "P"}, { 1.0 }, CanSwap::No);
    bool equal = !(a < b) && !(b < a);
    EXPECT_FALSE((equal));
}

TEST(IdentifierSimpleTest, EqualNotSwappableA) {
    Identifier a({"P", "H"}, { 1.0 }, CanSwap::No);
    Identifier b({"H", "P"}, { 1.0 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_FALSE((equal));
}

TEST(IdentifierSimpleTest, EqualNotSwappableB) {
    Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    Identifier b({"H", "P"}, { 1.0 }, CanSwap::No);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, ANotEqualToBBondOrder2) {
    Identifier a({"P", "H"}, { 2 }, CanSwap::Yes);
    Identifier b({"C", "H"}, { 2 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_FALSE((equal));
}

}

}
