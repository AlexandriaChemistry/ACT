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

#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! Read all content from a named file into a string.
static std::string readFile(const std::string &path)
{
    std::ifstream ifs(path);
    std::string   content;
    std::string   line;
    while (std::getline(ifs, line))
    {
        content += line + "\n";
    }
    return content;
}

// ---- setColumns() tests ----

TEST(LongTableTest, SetColumnsIntOneColumn)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setColumns(1);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.printHeader();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("{l}"));
}

TEST(LongTableTest, SetColumnsIntThreeColumns)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setColumns(3);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.printHeader();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("{lcc}"));
}

TEST(LongTableTest, SetColumnsString)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setColumns("lcrc");
    table.setCaption("cap");
    table.setLabel("lbl");
    table.printHeader();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("{lcrc}"));
}

// ---- printLine() escape logic tests ----

TEST(LongTableTest, PrintLineNoSpecialChars)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("hello world");
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    // No escaping: output contains the literal text + \\
    EXPECT_NE(std::string::npos, content.find("hello world"));
    EXPECT_NE(std::string::npos, content.find("\\\\"));
}

TEST(LongTableTest, PrintLineEscapesUnderscore)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("a_b");
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    // Underscore should be replaced by \_{...}
    EXPECT_NE(std::string::npos, content.find("\\_{"));
    // Literal underscore must not appear
    EXPECT_EQ(std::string::npos, content.find("a_b"));
}

TEST(LongTableTest, PrintLineEscapesHash)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("a#b");
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    // Hash is treated the same as underscore: \_{...}
    EXPECT_NE(std::string::npos, content.find("\\_{"));
    EXPECT_EQ(std::string::npos, content.find("a#b"));
}

TEST(LongTableTest, PrintLineUnderscoreAtEnd)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("abc_");
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\_{"));
    // The group must be closed
    EXPECT_NE(std::string::npos, content.find("}"));
}

TEST(LongTableTest, PrintLineMultipleUnderscores)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printLine("a_b_c");
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\_{"));
    EXPECT_EQ(std::string::npos, content.find("a_b_c"));
}

// ---- printColumns() tests ----

TEST(LongTableTest, PrintColumnsTwoColumns)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    std::vector<std::string> cols = {"Col1", "Col2"};
    table.printColumns(cols);
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("Col1 & Col2"));
}

TEST(LongTableTest, PrintColumnsThreeColumns)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    std::vector<std::string> cols = {"Name", "Value", "Unit"};
    table.printColumns(cols);
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("Name & Value & Unit"));
}

TEST(LongTableTest, PrintColumnsSingleColumn)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    std::vector<std::string> cols = {"OnlyColumn"};
    table.printColumns(cols);
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("OnlyColumn"));
    // No & separator for a single column
    EXPECT_EQ(std::string::npos, content.find("&"));
}

// ---- printHeader() and printFooter() tests ----

TEST(LongTableTest, PrintHeaderContainsLongtable)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setCaption("My caption");
    table.setLabel("tab:my");
    table.setColumns(2);
    table.printHeader();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\begin{longtable}"));
    EXPECT_NE(std::string::npos, content.find("My caption"));
    EXPECT_NE(std::string::npos, content.find("tab:my"));
    EXPECT_NE(std::string::npos, content.find("\\endfirsthead"));
    EXPECT_NE(std::string::npos, content.find("\\endhead"));
}

TEST(LongTableTest, PrintHeaderLandscape)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, true, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\begin{landscape}"));
}

TEST(LongTableTest, PrintHeaderNoLandscape)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_EQ(std::string::npos, content.find("\\begin{landscape}"));
}

TEST(LongTableTest, PrintHeaderWithFont)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, "footnotesize");
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(2);
    table.printHeader();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\begin{footnotesize}"));
}

TEST(LongTableTest, PrintHeaderWithSpacing)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr, true);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\begin{spacing}"));
}

TEST(LongTableTest, PrintFooterContainsEndLongtable)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    table.printFooter();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\end{longtable}"));
}

TEST(LongTableTest, PrintFooterLandscape)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, true, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    table.printFooter();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\end{landscape}"));
}

TEST(LongTableTest, PrintFooterWithFont)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, "footnotesize");
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(1);
    table.printHeader();
    table.printFooter();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\end{footnotesize}"));
}

TEST(LongTableTest, PrintHLine)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.printHLine();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("\\hline"));
}

// ---- addHeadLine() tests ----

TEST(LongTableTest, AddHeadLineAppearsInHeader)
{
    gmx::test::TestFileManager files;
    std::string path = files.getTemporaryFilePath(".tex");
    FILE *fp = fopen(path.c_str(), "w");
    ASSERT_NE(nullptr, fp);
    LongTable table(fp, false, nullptr);
    table.setCaption("cap");
    table.setLabel("lbl");
    table.setColumns(2);
    table.addHeadLine("Compound & Energy");
    table.printHeader();
    fflush(fp);
    fclose(fp);
    std::string content = readFile(path);
    EXPECT_NE(std::string::npos, content.find("Compound & Energy"));
}

// ---- NullFpConstructorThrows ----

TEST(LongTableTest, NullFpConstructorThrows)
{
    EXPECT_THROW((LongTable(nullptr, false, nullptr)), gmx::FileIOError);
}

} // namespace

} // namespace alexandria
