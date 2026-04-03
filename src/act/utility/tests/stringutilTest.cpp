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
 * Tests for the string utility functions in stringutil.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "../stringutil.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

// ---- split() tests ----

TEST(StringUtilTest, SplitByCommaBasic)
{
    auto result = split("a,b,c", ',');
    ASSERT_EQ(3u, result.size());
    EXPECT_EQ("a", result[0]);
    EXPECT_EQ("b", result[1]);
    EXPECT_EQ("c", result[2]);
}

TEST(StringUtilTest, SplitByCommaPreservesEmptyTokens)
{
    auto result = split("a,,b", ',');
    ASSERT_EQ(3u, result.size());
    EXPECT_EQ("a", result[0]);
    EXPECT_EQ("", result[1]);
    EXPECT_EQ("b", result[2]);
}

TEST(StringUtilTest, SplitBySpaceSquashesMultipleSpaces)
{
    // Space is whitespace so empty tokens are skipped
    auto result = split("a  b  c", ' ');
    ASSERT_EQ(3u, result.size());
    EXPECT_EQ("a", result[0]);
    EXPECT_EQ("b", result[1]);
    EXPECT_EQ("c", result[2]);
}

TEST(StringUtilTest, SplitBySpaceSingleWords)
{
    auto result = split("hello world", ' ');
    ASSERT_EQ(2u, result.size());
    EXPECT_EQ("hello", result[0]);
    EXPECT_EQ("world", result[1]);
}

TEST(StringUtilTest, SplitEmptyStringByCommaGivesNoTokens)
{
    auto result = split("", ',');
    EXPECT_TRUE(result.empty());
}

TEST(StringUtilTest, SplitSingleTokenNoDelimiter)
{
    auto result = split("hello", ',');
    ASSERT_EQ(1u, result.size());
    EXPECT_EQ("hello", result[0]);
}

TEST(StringUtilTest, SplitWithElemsAccumulates)
{
    std::vector<std::string> elems;
    elems.push_back("existing");
    split("a,b", ',', elems);
    ASSERT_EQ(3u, elems.size());
    EXPECT_EQ("existing", elems[0]);
    EXPECT_EQ("a", elems[1]);
    EXPECT_EQ("b", elems[2]);
}

// ---- gmx_ftoa() tests ----

TEST(StringUtilTest, GmxFtoaInRangeUsesFFormat)
{
    // fabs(f) > 1 && fabs(f) < 100 -> uses %.3f
    EXPECT_EQ("5.000", gmx_ftoa(5.0));
    EXPECT_EQ("50.500", gmx_ftoa(50.5));
    EXPECT_EQ("-5.000", gmx_ftoa(-5.0));
    EXPECT_EQ("1.500", gmx_ftoa(1.5));
}

TEST(StringUtilTest, GmxFtoaOutOfRangeUsesGFormat)
{
    // fabs(f) == 1.0 is NOT > 1, so uses %g
    EXPECT_EQ("1", gmx_ftoa(1.0));
    // fabs(f) == 100.0 is NOT < 100, so uses %g
    EXPECT_EQ("100", gmx_ftoa(100.0));
    // fabs(f) < 1
    EXPECT_EQ("0.5", gmx_ftoa(0.5));
    // Zero
    EXPECT_EQ("0", gmx_ftoa(0.0));
}

// ---- gmx_itoa() tests ----

TEST(StringUtilTest, GmxItoaPositive)
{
    EXPECT_EQ("42", gmx_itoa(42));
    EXPECT_EQ("0", gmx_itoa(0));
    EXPECT_EQ("1000", gmx_itoa(1000));
}

TEST(StringUtilTest, GmxItoaNegative)
{
    EXPECT_EQ("-7", gmx_itoa(-7));
}

// ---- my_atof() tests ----

TEST(StringUtilTest, MyAtofParsesDouble)
{
    EXPECT_DOUBLE_EQ(3.14, my_atof("3.14", "pi"));
    EXPECT_DOUBLE_EQ(-2.5, my_atof("-2.5", "neg"));
    EXPECT_DOUBLE_EQ(0.0, my_atof("0.0", "zero"));
}

TEST(StringUtilTest, MyAtofParsesIntegerString)
{
    EXPECT_DOUBLE_EQ(42.0, my_atof("42", "value"));
}

TEST(StringUtilTest, MyAtofEmptyStringReturnsZero)
{
    // std::atof on empty string returns 0.0
    EXPECT_DOUBLE_EQ(0.0, my_atof("", "empty"));
}

// ---- my_atoi() tests ----

TEST(StringUtilTest, MyAtoiParsesInt)
{
    EXPECT_EQ(7, my_atoi("7", "count"));
    EXPECT_EQ(-3, my_atoi("-3", "neg"));
    EXPECT_EQ(0, my_atoi("0", "zero"));
}

TEST(StringUtilTest, MyAtoiLargeNumber)
{
    EXPECT_EQ(12345, my_atoi("12345", "large"));
}

// ---- stringEqual() tests ----

TEST(StringUtilTest, StringEqualSameCase)
{
    EXPECT_TRUE(stringEqual("hello", "hello"));
    EXPECT_TRUE(stringEqual("ABC", "ABC"));
    EXPECT_TRUE(stringEqual("", ""));
}

TEST(StringUtilTest, StringEqualDifferentCase)
{
    EXPECT_TRUE(stringEqual("Hello", "hello"));
    EXPECT_TRUE(stringEqual("WORLD", "world"));
    EXPECT_TRUE(stringEqual("MiXeD", "mixed"));
}

TEST(StringUtilTest, StringEqualDifferentStrings)
{
    EXPECT_FALSE(stringEqual("hello", "world"));
    EXPECT_FALSE(stringEqual("abc", "abcd"));
    EXPECT_FALSE(stringEqual("abcd", "abc"));
}

TEST(StringUtilTest, StringEqualEmptyVsNonEmpty)
{
    EXPECT_FALSE(stringEqual("", "a"));
    EXPECT_FALSE(stringEqual("a", ""));
}

// ---- get_option() tests ----

TEST(StringUtilTest, GetOptionNullptrReturnsMinusOne)
{
    EXPECT_EQ(-1, get_option(nullptr));
}

TEST(StringUtilTest, GetOptionFoundAtFirstOption)
{
    const char *opts[] = {"A", "A", "B", "C", nullptr};
    EXPECT_EQ(0, get_option(opts));
}

TEST(StringUtilTest, GetOptionFoundAtSecondOption)
{
    const char *opts[] = {"B", "A", "B", "C", nullptr};
    EXPECT_EQ(1, get_option(opts));
}

TEST(StringUtilTest, GetOptionFoundAtThirdOption)
{
    const char *opts[] = {"C", "A", "B", "C", nullptr};
    EXPECT_EQ(2, get_option(opts));
}

TEST(StringUtilTest, GetOptionNotFoundReturnsZero)
{
    const char *opts[] = {"X", "A", "B", "C", nullptr};
    EXPECT_EQ(0, get_option(opts));
}

TEST(StringUtilTest, GetOptionCaseInsensitive)
{
    const char *opts[] = {"hello", "HELLO", "WORLD", nullptr};
    EXPECT_EQ(0, get_option(opts));
}
