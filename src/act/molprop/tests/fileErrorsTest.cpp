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
 * Tests that a molprop file that cannot be read or written gives an
 * exception that tells what is wrong, and leaves the process alive: a
 * program that uses ACT as a library, such as the Python module, catches it.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <fstream>
#include <gtest/gtest.h>
#include <iterator>
#include <libxml/xmlversion.h>
#include <string>
#include <vector>

#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/*! \brief Read a molprop file that cannot be read
 * \param[in] fileName The file
 * \return The message of the exception for a file that cannot be read,
 * empty if there is no such exception
 */
std::string fileErrorOfReading(const std::string &fileName)
{
    std::vector<MolProp> mps;
    try
    {
        MolPropRead(nullptr, fileName.c_str(), &mps);
    }
    catch (const gmx::FileIOError &ex)
    {
        return ex.what();
    }
    return "";
}

//! \return The compounds of a molprop file of the test data
std::vector<MolProp> testCompounds()
{
    std::vector<MolProp> mps;
    MolPropRead(nullptr, gmx::test::TestFileManager::getInputFilePath("molprop/molprop.xml").c_str(), &mps);
    return mps;
}

/*! \brief Tell whether a file is compressed
 * \param[in] fileName The file
 * \return Whether the file starts as a gzip file does
 */
bool isCompressed(const std::string &fileName)
{
    std::ifstream in(fileName, std::ios::binary);
    const int     first  = in.get();
    const int     second = in.get();
    return first == 0x1f && second == 0x8b;
}

// Test: A molprop file that does not exist, in a program and in a library
TEST(MolpropFileErrorsTest, MissingFile)
{
    const std::string name = "no-such-molprop-7f3a.xml";
    auto text = fileErrorOfReading(name);
    EXPECT_NE(std::string::npos, text.find(name)) << text;
    {
        const gmx::ThrowOnFatalError library;
        text = fileErrorOfReading(name);
        EXPECT_NE(std::string::npos, text.find(name)) << text;
    }
    // The next file is read as if nothing happened
    EXPECT_EQ(3U, testCompounds().size());
}

// Test: A file that is not XML, and one that ends halfway
TEST(MolpropFileErrorsTest, BrokenFile)
{
    gmx::test::TestFileManager fileManager;
    const std::string noXml = fileManager.getTemporaryFilePath("not-xml.xml");
    {
        std::ofstream out(noXml);
        out << "water|Train\n";
    }
    auto text = fileErrorOfReading(noXml);
    EXPECT_NE(std::string::npos, text.find("not-xml.xml")) << text;

    const std::string dataName = gmx::test::TestFileManager::getInputFilePath("molprop/molprop.xml");
    std::ifstream     in(dataName);
    const std::string xml((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    ASSERT_GT(xml.size(), 1000U);
    const std::string half = fileManager.getTemporaryFilePath("half.xml");
    {
        std::ofstream out(half);
        out << xml.substr(0, xml.size()/2);
    }
    {
        const gmx::ThrowOnFatalError library;
        text = fileErrorOfReading(half);
        EXPECT_NE(std::string::npos, text.find("half.xml")) << text;
    }
    EXPECT_EQ(3U, testCompounds().size());
}

// Test: A molprop file cannot be written where there is no directory
TEST(MolpropFileErrorsTest, CannotWrite)
{
    const auto        mps     = testCompounds();
    const std::string nowhere = "/no-such-directory-7f3a/sub/mp.xml";
    EXPECT_THROW(MolPropWrite(nowhere.c_str(), mps, false), gmx::FileIOError);
    {
        const gmx::ThrowOnFatalError library;
        try
        {
            MolPropWrite(nowhere.c_str(), mps, false);
            ADD_FAILURE() << "Writing to a directory that does not exist did not throw";
        }
        catch (const gmx::FileIOError &ex)
        {
            const std::string text = ex.what();
            EXPECT_NE(std::string::npos, text.find(nowhere)) << text;
        }
    }
    gmx::test::TestFileManager fileManager;
    const std::string          good = fileManager.getTemporaryFilePath("written.xml");
    EXPECT_NO_THROW(MolPropWrite(good.c_str(), mps, false));
    std::vector<MolProp> again;
    MolPropRead(nullptr, good.c_str(), &again);
    EXPECT_EQ(mps.size(), again.size());
}

// Test: A compressed molprop file is read as one that is not, with any version of libxml2
TEST(MolpropFileErrorsTest, CompressedFile)
{
    const auto                 mps = testCompounds();
    gmx::test::TestFileManager fileManager;
    const std::string          packed = fileManager.getTemporaryFilePath("packed.xml");
    MolPropWrite(packed.c_str(), mps, true);
#ifdef LIBXML_ZLIB_ENABLED
    EXPECT_TRUE(isCompressed(packed));
#endif
    const std::string plain = fileManager.getTemporaryFilePath("plain.xml");
    MolPropWrite(plain.c_str(), mps, false);
    EXPECT_FALSE(isCompressed(plain));
    std::vector<MolProp> again;
    MolPropRead(nullptr, packed.c_str(), &again);
    EXPECT_EQ(mps.size(), again.size());
    EXPECT_GT(again.size(), 0U);
}

//! \return The contents of a file
std::string contents(const std::string &fileName)
{
    std::ifstream in(fileName);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

// Test: A copy of compounds is a copy of all they have, including the
// properties of their experiments, and does not share them
TEST(MolpropCopyTest, CopiesAreDeep)
{
    const auto           mps  = testCompounds();
    std::vector<MolProp> copy = mps;
    gmx::test::TestFileManager fileManager;
    const std::string          original = fileManager.getTemporaryFilePath("original.xml");
    const std::string          copied   = fileManager.getTemporaryFilePath("copied.xml");
    MolPropWrite(original.c_str(), mps, false);
    MolPropWrite(copied.c_str(), copy, false);
    EXPECT_EQ(contents(original), contents(copied));
    // Without the properties of the copy the original keeps its own
    size_t before = 0;
    for (const auto &experiment : mps[0].experimentConst())
    {
        before += experiment.propertiesConst().size();
    }
    ASSERT_GT(before, 0U);
    for (auto &experiment : *copy[0].experiment())
    {
        experiment.properties()->clear();
    }
    size_t after = 0;
    for (const auto &experiment : mps[0].experimentConst())
    {
        after += experiment.propertiesConst().size();
    }
    EXPECT_EQ(before, after);
    // And a copy that is assigned
    std::vector<MolProp> assigned;
    assigned = mps;
    MolPropWrite(copied.c_str(), assigned, false);
    EXPECT_EQ(contents(original), contents(copied));
}

// Test: A file is written and read one piece at a time. What is written
// is what is read, also when a file is not laid out as ACT writes it
TEST(MolpropStreamTest, WrittenAsRead)
{
    const auto                 mps = testCompounds();
    gmx::test::TestFileManager fileManager;
    const std::string          first  = fileManager.getTemporaryFilePath("first.xml");
    const std::string          second = fileManager.getTemporaryFilePath("second.xml");
    MolPropWrite(first.c_str(), mps, false);
    const std::string text = contents(first);
    const std::string head = "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
                             "<!DOCTYPE molprops.dtd PUBLIC \"molprops.dtd\" \"molprops.dtd\">\n"
                             "<molecules>\n  <molecule molname=\"";
    const std::string tail = "\n  </molecule>\n</molecules>\n";
    ASSERT_GT(text.size(), head.size() + tail.size());
    EXPECT_EQ(head, text.substr(0, head.size()));
    EXPECT_EQ(tail, text.substr(text.size() - tail.size()));
    {
        std::vector<MolProp> again;
        MolPropRead(nullptr, first.c_str(), &again);
        ASSERT_EQ(mps.size(), again.size());
        MolPropWrite(second.c_str(), again, false);
        EXPECT_EQ(text, contents(second));
    }
    // Comments around and between the compounds and the experiments, and
    // no line breaks between them
    std::string odd = text;
    for (const std::string &element : { "<molecules>", "<molecule ", "<experiment ", "</molecule>", "</molecules>" })
    {
        size_t at = 0;
        while ((at = odd.find(element, at)) != std::string::npos)
        {
            const std::string note = "<!-- a note -->";
            odd.insert(at, note);
            at += note.size() + element.size();
        }
    }
    size_t at = 0;
    while ((at = odd.find(">\n", at)) != std::string::npos)
    {
        odd.erase(at + 1, 1);
    }
    {
        std::ofstream out(second);
        out << odd;
    }
    {
        std::vector<MolProp> again;
        MolPropRead(nullptr, second.c_str(), &again);
        ASSERT_EQ(mps.size(), again.size());
        MolPropWrite(second.c_str(), again, false);
        EXPECT_EQ(text, contents(second));
    }
}

// Test: A file without compounds, and names with characters that XML and
// the encoding of the file cannot have as they are
TEST(MolpropStreamTest, EmptyAndSpecial)
{
    gmx::test::TestFileManager fileManager;
    const std::string          none = fileManager.getTemporaryFilePath("none.xml");
    MolPropWrite(none.c_str(), {}, false);
    EXPECT_NE(std::string::npos, contents(none).find("<molecules/>")) << contents(none);
    std::vector<MolProp> again;
    EXPECT_NO_THROW(MolPropRead(nullptr, none.c_str(), &again));
    EXPECT_TRUE(again.empty());

    auto              mps  = testCompounds();
    const std::string name = "caf\xc3\xa9 \xce\xb2 & <x> \"q\" 's";
    mps[0].SetMolname(name);
    const std::string special = fileManager.getTemporaryFilePath("special.xml");
    MolPropWrite(special.c_str(), mps, false);
    MolPropRead(nullptr, special.c_str(), &again);
    ASSERT_EQ(mps.size(), again.size());
    EXPECT_EQ(name, again[0].getMolname());
    EXPECT_EQ(mps[1].getMolname(), again[1].getMolname());
}

// Test: Of a file that ends halfway nothing is kept, and what was there stays
TEST(MolpropStreamTest, HalfAFileAddsNothing)
{
    gmx::test::TestFileManager fileManager;
    auto                       mps    = testCompounds();
    const size_t               before = mps.size();
    const std::string          whole  = fileManager.getTemporaryFilePath("whole.xml");
    MolPropWrite(whole.c_str(), mps, false);
    const std::string text = contents(whole);
    // The end is in the last compound, after others were read
    const std::string half = fileManager.getTemporaryFilePath("most.xml");
    {
        std::ofstream out(half);
        out << text.substr(0, text.size() - 200);
    }
    {
        const gmx::ThrowOnFatalError library;
        EXPECT_THROW(MolPropRead(nullptr, half.c_str(), &mps), gmx::FileIOError);
    }
    EXPECT_EQ(before, mps.size());
    // Reading adds to what is there
    MolPropRead(nullptr, whole.c_str(), &mps);
    EXPECT_EQ(2*before, mps.size());
}

} // namespace

} // namespace alexandria
