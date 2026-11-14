/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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
 * Tests of the JSON reader, see jsonreader.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include "act/utility/jsonreader.h"
#include "gromacs/utility/exceptions.h"

namespace alexandria
{

namespace
{

TEST(JsonReaderTest, Scalars)
{
    EXPECT_TRUE(parseJson("null").isNull());
    EXPECT_TRUE(parseJson(" true ").asBool());
    EXPECT_FALSE(parseJson("false").asBool());
    EXPECT_EQ(42, parseJson("42").asInt());
    EXPECT_EQ(-7, parseJson("-7").asInt());
    EXPECT_EQ(112.985, parseJson("112.985").asNumber());
    EXPECT_EQ(1.5e-3, parseJson("1.5e-3").asNumber());
    EXPECT_EQ(2e10, parseJson("2E+10").asNumber());
    EXPECT_EQ("abc", parseJson("\"abc\"").asString());
}

TEST(JsonReaderTest, Escapes)
{
    EXPECT_EQ("a\"b\\c/d\n\t", parseJson("\"a\\\"b\\\\c\\/d\\n\\t\"").asString());
    // The plus-minus sign as an escape and as UTF-8 are the same
    EXPECT_EQ("(\xc2\xb1)", parseJson("\"(\\u00b1)\"").asString());
    EXPECT_EQ("(\xc2\xb1)", parseJson("\"(\xc2\xb1)\"").asString());
    // Three bytes, and a pair of surrogates for four bytes
    EXPECT_EQ("\xe2\x82\xac", parseJson("\"\\u20AC\"").asString());
    EXPECT_EQ("\xf0\x9f\x98\x80", parseJson("\"\\ud83d\\ude00\"").asString());
}

TEST(JsonReaderTest, Structures)
{
    const auto doc = parseJson("{ \"name\": \"water\", \"atoms\": [ 8, 1, 1 ], \"empty\": {}, \"none\": [] }");
    ASSERT_EQ(JsonValue::Type::Object, doc.type());
    EXPECT_EQ(4u, doc.asObject().size());
    EXPECT_EQ("water", doc.find("name")->asString());
    const auto &atoms = doc.find("atoms")->asArray();
    ASSERT_EQ(3u, atoms.size());
    EXPECT_EQ(8, atoms[0].asInt());
    EXPECT_EQ(0u, doc.find("empty")->asObject().size());
    EXPECT_EQ(0u, doc.find("none")->asArray().size());
    EXPECT_TRUE(doc.find("missing") == nullptr);
    // The members in their order
    EXPECT_EQ("atoms", doc.asObject()[1].first);
}

TEST(JsonReaderTest, Errors)
{
    for (const char *bad : { "", "{", "[1,]", "{\"a\" 1}", "{\"a\":1,}", "tru", "01", "1.", "-",
                             "\"abc", "\"\\x\"", "\"\\ud83d\"", "[1] 2", "\"a\tb\"" })
    {
        EXPECT_THROW(parseJson(bad), gmx::InvalidInputError) << bad;
    }
    // The wrong type
    EXPECT_THROW(parseJson("\"1\"").asNumber(), gmx::InvalidInputError);
    EXPECT_THROW(parseJson("1.5").asInt(), gmx::InvalidInputError);
    EXPECT_THROW(parseJson("[]").find("a"), gmx::InvalidInputError);
}

} // namespace

} // namespace alexandria
