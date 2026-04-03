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
 * Tests for the memory_check utility functions.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "../memory_check.h"

#include <cstdio>
#include <string>

#include <gtest/gtest.h>

// ---- memory_usage_low() tests ----

TEST(MemoryCheckTest, MemoryUsageLowReturnsValidContent)
{
    std::string result = memory_usage_low("testfile.cpp", 10);
    // The result is either empty (RSS unavailable) or contains "VMEM"
    EXPECT_TRUE(result.empty() || result.find("VMEM") != std::string::npos);
}

TEST(MemoryCheckTest, MemoryUsageLowContainsFileAndLine)
{
    std::string result = memory_usage_low("myfile.cpp", 99);
    if (!result.empty())
    {
        EXPECT_NE(std::string::npos, result.find("myfile.cpp"));
        EXPECT_NE(std::string::npos, result.find("99"));
    }
}

TEST(MemoryCheckTest, MemoryUsageMacroWorks)
{
    // Just verify the macro expands and doesn't crash
    std::string result = memory_usage();
    EXPECT_TRUE(result.empty() || result.find("VMEM") != std::string::npos);
}

// ---- print_memory_usage_low() tests ----

TEST(MemoryCheckTest, PrintMemoryUsageLowNullFpDoesNotCrash)
{
    EXPECT_NO_THROW(print_memory_usage_low(nullptr, "test.cpp", 1));
}

TEST(MemoryCheckTest, PrintMemoryUsageLowWritesToFile)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    print_memory_usage_low(fp, "source.cpp", 42);
    fflush(fp);
    rewind(fp);
    char buf[512] = {};
    fgets(buf, sizeof(buf), fp);
    fclose(fp);
    std::string content(buf);
    // Either nothing was written (RSS unavailable) or the line contains VMEM
    EXPECT_TRUE(content.empty() || content.find("VMEM") != std::string::npos);
    if (!content.empty())
    {
        EXPECT_NE(std::string::npos, content.find("source.cpp"));
        EXPECT_NE(std::string::npos, content.find("42"));
    }
}

TEST(MemoryCheckTest, PrintMemoryUsageMacroWorks)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    // print_memory_usage expands to print_memory_usage_low(fp, __FILE__, __LINE__)
    print_memory_usage(fp);
    fflush(fp);
    fclose(fp);
    // Just verify no crash occurs
}
