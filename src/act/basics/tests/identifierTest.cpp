/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2020-2024
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

#include "../identifier.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

TEST(IdentifierSimpleTest, ConstructorOne) {
    EXPECT_THROW(Identifier id({"C", "H"}, { 1, 2 }, CanSwap::Yes), gmx::InvalidInputError);
}

TEST(IdentifierSimpleTest, BSmallerThanA) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    const Identifier b({"C", "H"}, { 1.0 }, CanSwap::Yes);
    bool compare = b < a;
    EXPECT_TRUE((compare));
}

TEST(IdentifierSimpleTest, SortsCanSwapLinear) {
    const Identifier a({"H", "C", "N"}, { 1.0, 3.0 }, CanSwap::Linear);
    const Identifier b({"N", "C", "H"}, { 3.0, 1.0 }, CanSwap::Linear);
    EXPECT_TRUE(a.id() == "H~C#N");
    EXPECT_TRUE(a.swapped() == "N#C~H");
    EXPECT_TRUE(a == b);
}

TEST(IdentifierSimpleTest, SortsCanSwapVsite2) {
    const Identifier a({"P", "H", "V"}, { 1.0, 9 }, CanSwap::Vsite2);
    EXPECT_TRUE(a.id() == "P~H!V");
    EXPECT_TRUE(a.swapped() == "H~P!V");
}

TEST(IdentifierSimpleTest, SortsCanSwapYes) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    EXPECT_TRUE(a.id() == "H~P");
}

TEST(IdentifierSimpleTest, SortsCanSwapNo) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::No);
    EXPECT_TRUE(a.id() == "P~H");
}

TEST(IdentifierSimpleTest, CompareDifferentLength) {
    const Identifier a({"C", "H"}, { 1.0 }, CanSwap::Yes);
    const Identifier b({"C", "O", "H"}, { 1.0, 1.0 }, CanSwap::Yes);
    bool compare = b < a;
    EXPECT_FALSE(compare);
}

TEST(IdentifierSimpleTest, BSmallerThanASecondAtom) {
    const Identifier a({"H", "P"}, { 1.0 }, CanSwap::No);
    const Identifier b({"H", "C"}, { 1.0 }, CanSwap::No);
    bool compare = a < b;
    EXPECT_FALSE(compare);
}

TEST(IdentifierSimpleTest, BSmallerThanABondOrder) {
    const Identifier a({"H", "C"}, { 2 }, CanSwap::No);
    const Identifier b({"H", "C"}, { 1 }, CanSwap::No);
    // Since bond order is ignored this should yield false
    bool compare = b < a;
    EXPECT_FALSE(compare);
}

TEST(IdentifierSimpleTest, NotASmallerThanB) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    const Identifier b({"C", "H"}, { 1.0 }, CanSwap::Yes);
    bool compare = a < b;
    EXPECT_FALSE(compare);
}

TEST(IdentifierSimpleTest, AEqualToB) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    const Identifier b({"P", "H"}, { 1.0 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE(equal);
}

TEST(IdentifierSimpleTest, SwappedAEqualToB) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    const Identifier b({"H", "P"}, { 1.0 }, CanSwap::Yes);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE(equal);
}

TEST(IdentifierSimpleTest, ANotEqualToB) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    const Identifier b({"C", "H"}, { 1.0 }, CanSwap::Yes);
    EXPECT_FALSE(a < b);
    EXPECT_TRUE(b < a);
    EXPECT_FALSE(a == b);
}

TEST(IdentifierSimpleTest, BSmallerThanABondOrder1) {
    const Identifier a({"P", "H"}, { 1 }, CanSwap::Yes);
    const Identifier b({"C", "H"}, { 1 }, CanSwap::Yes);
    EXPECT_TRUE(b < a);
    EXPECT_FALSE(a < b);
    EXPECT_FALSE(a == b);
}

TEST(IdentifierSimpleTest, NotASmallerThanBBondOrder2) {
    const Identifier a({"P", "H"}, { 2 }, CanSwap::Yes);
    const Identifier b({"C", "H"}, { 2 }, CanSwap::Yes);
    EXPECT_FALSE(a < b);
}

TEST(IdentifierSimpleTest, AEqualToBBondOrder1) {
    const Identifier a({"P", "H"}, { 1 }, CanSwap::Yes);
    const Identifier b({"P", "H"}, { 1 }, CanSwap::Yes);
    EXPECT_TRUE(a == b);
}

TEST(IdentifierSimpleTest, AEqualToBBondOrderDiff) {
    const Identifier a({"P", "H"}, { 1 }, CanSwap::Yes);
    const Identifier b({"P", "H"}, { 3 }, CanSwap::Yes);
    EXPECT_FALSE(a == b);
}

TEST(IdentifierSimpleTest, SwappedAEqualToBBondOrder2) {
    const Identifier a({"P", "H"}, { 2 }, CanSwap::Yes);
    const Identifier b({"H", "P"}, { 2 }, CanSwap::Yes);
    EXPECT_FALSE(a < b);
    EXPECT_FALSE(b < a);
    EXPECT_TRUE(a == b);
}

TEST(IdentifierSimpleTest, EqualNotSwappableAB) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::No);
    const Identifier b({"H", "P"}, { 1.0 }, CanSwap::No);
    EXPECT_FALSE(a < b);
    EXPECT_TRUE(b < a);
}

TEST(IdentifierSimpleTest, EqualNotSwappableA) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::No);
    const Identifier b({"H", "P"}, { 1.0 }, CanSwap::Yes);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(b < a);
    EXPECT_FALSE(a == b);
}

TEST(IdentifierSimpleTest, EqualNotSwappableB) {
    const Identifier a({"P", "H"}, { 1.0 }, CanSwap::Yes);
    const Identifier b({"H", "P"}, { 1.0 }, CanSwap::No);
    EXPECT_FALSE(a < b);
    EXPECT_TRUE(b < a);
    EXPECT_FALSE(a == b);
}

TEST(IdentifierSimpleTest, ANotEqualToBBondOrder2) {
    const Identifier a({"P", "H"}, { 2 }, CanSwap::Yes);
    const Identifier b({"C", "H"}, { 2 }, CanSwap::Yes);
    EXPECT_FALSE(a < b);
    EXPECT_TRUE(b < a);
    EXPECT_FALSE(a == b);
}

TEST(IdentifierSimpleTest, DiffLengthCanSwap) {
    const Identifier a({"P" }, { }, CanSwap::No);
    const Identifier b({"C", "H"}, { 2 }, CanSwap::Yes);
    EXPECT_FALSE(a == b);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(b < a);
}

TEST(IdentifierSimpleTest, DiffLength) {
    const Identifier a({"P" }, { }, CanSwap::Yes);
    const Identifier b({"C", "H"}, { 2 }, CanSwap::Yes);
    EXPECT_FALSE(a == b);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(b < a);
}

TEST(IdentifierSimpleTest, LinearNotEqual) {
    const Identifier a({"P", "H", "X"}, { 1, 9 }, CanSwap::Linear);
    const Identifier b({"C", "H", "X"}, { 2, 9 }, CanSwap::Linear);
    bool equal = !(a < b) && !(b < a);
    EXPECT_FALSE((equal));
}

TEST(IdentifierSimpleTest, LinearEqual) {
    const Identifier a({"P", "H", "X"}, { 1, 9 }, CanSwap::Linear);
    const Identifier b({"P", "H", "X"}, { 1, 9 }, CanSwap::Linear);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, LinearSwappedEqual) {
    const Identifier a({"H", "P", "O"}, { 1, 2 }, CanSwap::Linear);
    const Identifier b({"O", "P", "H"}, { 2, 1 }, CanSwap::Linear);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, Vsite2NotEqual) {
    const Identifier a({"P", "H", "X"}, { 1, 9 }, CanSwap::Vsite2);
    const Identifier b({"C", "H", "X"}, { 2, 9 }, CanSwap::Vsite2);
    bool equal = !(a < b) && !(b < a);
    EXPECT_FALSE((equal));
}

TEST(IdentifierSimpleTest, Vsite2Equal) {
    const Identifier a({"P", "H", "X"}, { 1, 9 }, CanSwap::Vsite2);
    const Identifier b({"P", "H", "X"}, { 1, 9 }, CanSwap::Vsite2);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, Vsite2SwappedEqual) {
    const Identifier a({"H", "P", "X"}, { 1, 9 }, CanSwap::Vsite2);
    const Identifier b({"H", "P", "X"}, { 1, 9 }, CanSwap::Vsite2);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, Idih1a) {
    const std::vector<double> bOa = { 1, 1.5, 1.5 };
    const std::vector<double> bOb = { 1.5, 1, 1.5 };
    const std::vector<double> bOc = { 1.5, 1.5, 1 };
    const std::vector<std::string> aa = {"C", "H", "C", "C"};
    const std::vector<std::string> ab = {"C", "C", "H", "C"};
    const std::vector<std::string> ac = {"C", "C", "C", "H"};
    const Identifier a(aa, bOa, CanSwap::Idih);
    const Identifier b(ab, bOb, CanSwap::Idih);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(bOb == b.bondOrders());
    EXPECT_FALSE(bOa == a.bondOrders());
    EXPECT_TRUE(bOc == a.bondOrders());
    EXPECT_TRUE(bOc == b.bondOrders());
    EXPECT_FALSE(aa == a.atoms());
    EXPECT_FALSE(ab == b.atoms());
    EXPECT_TRUE(ac == a.atoms());
    EXPECT_TRUE(ac == b.atoms());
}

TEST(IdentifierSimpleTest, Idih1b) {
    const Identifier a({"C", "H", "C", "C"}, { 1, 1.5, 1.5 }, CanSwap::Idih);
    const Identifier b({"C", "C", "C", "H"}, { 1.5, 1.5, 1 }, CanSwap::Idih);
    bool equal = a == b;
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, Idih1c) {
    const Identifier a({"C", "C", "H", "C"}, { 1.5, 1, 1.5 }, CanSwap::Idih);
    const Identifier b({"C", "C", "C", "H"}, { 1.5, 1.5, 1 }, CanSwap::Idih);
    bool equal = a == b;
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, Idih1d) {
    const Identifier a({"C", "C", "H", "C"}, { 1.5, 1, 1.5 }, CanSwap::Idih);
    const Identifier b({"C", "C", "C", "H"}, { 1.5, 1.5, 1 }, CanSwap::Idih);
    bool equal = !(a < b) && !(b < a);
    EXPECT_TRUE((equal));
}

TEST(IdentifierSimpleTest, Idih2) {
    const Identifier a({"C", "H", "C", "C"}, { 1, 1.5, 1.5 }, CanSwap::Idih);
    const Identifier b({"C", "C", "F", "C"}, { 1.5, 1, 1.5 }, CanSwap::Idih);
    bool equal = a == b;
    EXPECT_FALSE((equal));
}

TEST(IdentifierSimpleTest, Idih3) {
    const Identifier a({"C", "H", "C", "C"}, { 1, 1.5, 1.5 }, CanSwap::Idih);
    const Identifier b({"C", "C", "F", "C"}, { 1.5, 1, 1.5 }, CanSwap::Idih);
    bool smaller = b < a;
    EXPECT_TRUE((smaller));
}

TEST(IdentifierSimpleTest, Idih4) {
    const Identifier a({"C", "H", "C", "C"}, { 1, 1.5, 1.5 }, CanSwap::Idih);
    const Identifier b({"C", "C", "F", "C"}, { 1.5, 1, 1.5 }, CanSwap::Idih);
    bool smaller = a < b;
    EXPECT_FALSE((smaller));
}

TEST(IdentifierSimpleTest, Idih5) {
    const Identifier a({"C", "H", "C", "C"}, { 1, 1.5, 1.5 }, CanSwap::Idih);
    const Identifier b({"C", "C", "F", "C"}, { 1.5, 1, 1.5 }, CanSwap::Idih);
    bool equal = a == b;
    EXPECT_FALSE((equal));
}

}

}
