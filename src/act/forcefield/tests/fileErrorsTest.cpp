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
 * Tests that a force field file that cannot be read or written gives an
 * exception that tells what is wrong, and leaves the process alive: a
 * program that uses ACT as a library, such as the Python module, catches it.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <fstream>
#include <gtest/gtest.h>
#include <string>

#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/*! \brief Read a force field file that cannot be read
 * \param[in] fileName The file
 * \return The message of the exception for a file that cannot be read,
 * empty if there is no such exception
 */
std::string fileErrorOfReading(const std::string &fileName)
{
    ForceField pd;
    try
    {
        readForceField(fileName, &pd);
    }
    catch (const gmx::FileIOError &ex)
    {
        return ex.what();
    }
    return "";
}

//! \return The number of particle types of a force field that is read from the test data
size_t particleTypesOfTestForceField()
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACM-g1.xml"), &pd);
    return pd.particleTypesConst().size();
}

// Test: A force field file that does not exist, in a program and in a library
TEST(ForceFieldFileErrorsTest, MissingFile)
{
    const std::string name = "no-such-force-field-7f3a.xml";
    auto text = fileErrorOfReading(name);
    EXPECT_NE(std::string::npos, text.find(name)) << text;
    {
        const gmx::ThrowOnFatalError library;
        text = fileErrorOfReading(name);
        EXPECT_NE(std::string::npos, text.find(name)) << text;
    }
    // The next force field is read as if nothing happened
    EXPECT_GT(particleTypesOfTestForceField(), 0U);
}

// Test: A file that is not XML, and one that ends halfway
TEST(ForceFieldFileErrorsTest, BrokenFile)
{
    gmx::test::TestFileManager fileManager;
    const std::string noXml = fileManager.getTemporaryFilePath("not-xml.xml");
    {
        std::ofstream out(noXml);
        out << "sigma 0.3\nepsilon 0.5\n";
    }
    auto text = fileErrorOfReading(noXml);
    EXPECT_NE(std::string::npos, text.find("not-xml.xml")) << text;

    // The first half of a real force field
    const std::string dataName = gmx::test::TestFileManager::getInputFilePath("forcefields/ACM-g1.xml");
    std::ifstream     in(dataName);
    const std::string xml((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    ASSERT_GT(xml.size(), 1000U);
    const std::string half = fileManager.getTemporaryFilePath("half.xml");
    {
        std::ofstream out(half);
        out << xml.substr(0, xml.size()/2);
    }
    text = fileErrorOfReading(half);
    EXPECT_NE(std::string::npos, text.find("half.xml")) << text;
    EXPECT_GT(particleTypesOfTestForceField(), 0U);
}

// Test: A force field cannot be written where there is no directory
TEST(ForceFieldFileErrorsTest, CannotWrite)
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACM-g1.xml"), &pd);
    const std::string nowhere = "/no-such-directory-7f3a/sub/ff.xml";
    EXPECT_THROW(writeForceField(nowhere, &pd, false), gmx::FileIOError);
    {
        const gmx::ThrowOnFatalError library;
        try
        {
            writeForceField(nowhere, &pd, false);
            ADD_FAILURE() << "Writing to a directory that does not exist did not throw";
        }
        catch (const gmx::FileIOError &ex)
        {
            const std::string text = ex.what();
            EXPECT_NE(std::string::npos, text.find(nowhere)) << text;
        }
    }
    // And where it can be written, it is
    gmx::test::TestFileManager fileManager;
    const std::string          good = fileManager.getTemporaryFilePath("written.xml");
    EXPECT_NO_THROW(writeForceField(good, &pd, false));
    ForceField again;
    readForceField(good, &again);
    EXPECT_EQ(pd.particleTypesConst().size(), again.particleTypesConst().size());
}

} // namespace

} // namespace alexandria
