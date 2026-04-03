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
 * Tests for the LongTable LaTeX utility class.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "../latex_util.h"

#include <cstdio>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"

namespace alexandria
{

namespace
{

//! Read all content from a rewound FILE* into a string.
static std::string readFile(FILE *fp)
{
    rewind(fp);
    std::string content;
    char buf[256];
    while (fgets(buf, sizeof(buf), fp) != nullptr)
    {
        content += buf;
    }
    return content;
}

// ---- setColumns() tests ----

TEST(LongTableTest, SetColumnsIntOneColumn)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setColumns(1);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.printHeader();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("{l}"));
}

TEST(LongTableTest, SetColumnsIntThreeColumns)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setColumns(3);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.printHeader();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("{lcc}"));
}

TEST(LongTableTest, SetColumnsString)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setColumns("lcrc");
    table.setCaption("cap");
    table.setLabel("lbl");
    table.printHeader();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("{lcrc}"));
}

// ---- printLine() escape logic tests ----

TEST(LongTableTest, PrintLineNoSpecialChars)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("hello world");
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    // No escaping: output contains the literal text + \\
    EXPECT_NE(std::string::npos, content.find("hello world"));
    EXPECT_NE(std::string::npos, content.find("\\\\"));
}

TEST(LongTableTest, PrintLineEscapesUnderscore)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("a_b");
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    // Underscore should be replaced by \_{...}
    EXPECT_NE(std::string::npos, content.find("\\_{"));
    // Literal underscore must not appear
    EXPECT_EQ(std::string::npos, content.find("a_b"));
}

TEST(LongTableTest, PrintLineEscapesHash)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("a#b");
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    // Hash is treated the same as underscore: \_{...}
    EXPECT_NE(std::string::npos, content.find("\\_{"));
    EXPECT_EQ(std::string::npos, content.find("a#b"));
}

TEST(LongTableTest, PrintLineUnderscoreAtEnd)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("abc_");
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\_{"));
    // The group must be closed
    EXPECT_NE(std::string::npos, content.find("}"));
}

TEST(LongTableTest, PrintLineMultipleUnderscores)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("a_b_c");
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\_{"));
    EXPECT_EQ(std::string::npos, content.find("a_b_c"));
}

// ---- printColumns() tests ----

TEST(LongTableTest, PrintColumnsTwoColumns)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    std::vector<std::string> cols = {"Col1", "Col2"};
    table.printColumns(cols);
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("Col1 & Col2"));
}

TEST(LongTableTest, PrintColumnsThreeColumns)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    std::vector<std::string> cols = {"Name", "Value", "Unit"};
    table.printColumns(cols);
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("Name & Value & Unit"));
}

TEST(LongTableTest, PrintColumnsSingleColumn)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    std::vector<std::string> cols = {"OnlyColumn"};
    table.printColumns(cols);
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("OnlyColumn"));
    // No & separator for a single column
    EXPECT_EQ(std::string::npos, content.find("&"));
}

// ---- printHeader() and printFooter() tests ----

TEST(LongTableTest, PrintHeaderContainsLongtable)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setCaption("My caption");
    table.setLabel("tab:my");
    table.setColumns(2);
    table.printHeader();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\begin{longtable}"));
    EXPECT_NE(std::string::npos, content.find("My caption"));
    EXPECT_NE(std::string::npos, content.find("tab:my"));
    EXPECT_NE(std::string::npos, content.find("\\endfirsthead"));
    EXPECT_NE(std::string::npos, content.find("\\endhead"));
}

TEST(LongTableTest, PrintHeaderLandscape)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, true, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\begin{landscape}"));
}

TEST(LongTableTest, PrintHeaderNoLandscape)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_EQ(std::string::npos, content.find("\\begin{landscape}"));
}

TEST(LongTableTest, PrintHeaderWithFont)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, "footnotesize");
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(2);
    table.printHeader();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\begin{footnotesize}"));
}

TEST(LongTableTest, PrintHeaderWithSpacing)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr, true);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\begin{spacing}"));
}

TEST(LongTableTest, PrintFooterContainsEndLongtable)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    table.printFooter();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\end{longtable}"));
}

TEST(LongTableTest, PrintFooterLandscape)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, true, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    table.printFooter();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\end{landscape}"));
}

TEST(LongTableTest, PrintFooterWithFont)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, "footnotesize");
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    table.printFooter();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\end{footnotesize}"));
}

TEST(LongTableTest, PrintHLine)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printHLine();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("\\hline"));
}

// ---- addHeadLine() tests ----

TEST(LongTableTest, AddHeadLineAppearsInHeader)
{
    FILE *fp = tmpfile();
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(2);
    table.addHeadLine("Compound & Energy");
    table.printHeader();
    fflush(fp);
    std::string content = readFile(fp);
    fclose(fp);
    EXPECT_NE(std::string::npos, content.find("Compound & Energy"));
}

// ---- NullFpConstructorThrows ----

TEST(LongTableTest, NullFpConstructorThrows)
{
    EXPECT_THROW((LongTable(nullptr, false, nullptr)), gmx::FileIOError);
}

} // namespace

} // namespace alexandria
