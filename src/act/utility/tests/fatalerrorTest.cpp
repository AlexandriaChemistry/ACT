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
 * Tests that fatal errors can be caught when ACT is used as a library,
 * as the Python module does: the process must stay alive.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <gtest/gtest.h>
#include <string>

#include "act/basics/msg_handler.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

namespace
{

//! \return The message of an exception that a statement throws, empty if it throws none
template<typename Function>
std::string messageOf(Function function)
{
    try
    {
        function();
    }
    catch (const std::exception &ex)
    {
        return ex.what();
    }
    return "";
}

// Test: A program ends on a fatal error, a library throws when it is told to
TEST(FatalErrorTest, SwitchIsOffUnlessAsked)
{
    EXPECT_FALSE(gmx_fatal_error_throws());
    {
        const gmx::ThrowOnFatalError library;
        EXPECT_TRUE(gmx_fatal_error_throws());
        {
            // Nested, as when a function of the library calls another one
            const gmx::ThrowOnFatalError inner;
            EXPECT_TRUE(gmx_fatal_error_throws());
        }
        EXPECT_TRUE(gmx_fatal_error_throws());
    }
    EXPECT_FALSE(gmx_fatal_error_throws());
    gmx_set_fatal_error_throws(true);
    EXPECT_TRUE(gmx_fatal_error_throws());
    gmx_set_fatal_error_throws(false);
    EXPECT_FALSE(gmx_fatal_error_throws());
}

// Test: gmx_fatal throws an exception with its message
TEST(FatalErrorTest, FatalThrows)
{
    const gmx::ThrowOnFatalError library;
    EXPECT_THROW(gmx_fatal(FARGS, "Cannot find atomtype %s in %d molecules", "ow", 3), gmx::FatalError);
    // All exceptions of the library have one base class and are standard exceptions
    EXPECT_THROW(gmx_fatal(FARGS, "No"), gmx::GromacsException);
    EXPECT_THROW(gmx_fatal(FARGS, "No"), std::exception);
    EXPECT_EQ("Cannot find atomtype ow in 3 molecules",
              messageOf([]() { gmx_fatal(FARGS, "Cannot find atomtype %s in %d molecules", "ow", 3); }));
    // The process is alive, and the next error is one of its own
    EXPECT_EQ("Second", messageOf([]() { gmx_fatal(FARGS, "Second"); }));
}

// Test: The other fatal error routines throw the exception that fits their kind
TEST(FatalErrorTest, ErrorsByKindThrow)
{
    const gmx::ThrowOnFatalError library;
    EXPECT_THROW(gmx_file("ff.xml"), gmx::FileIOError);
    EXPECT_THROW(gmx_open("ff.xml"), gmx::FileIOError);
    EXPECT_THROW(gmx_input("No such option"), gmx::InvalidInputError);
    EXPECT_THROW(gmx_incons("Sizes differ"), gmx::InternalError);
    EXPECT_THROW(gmx_call("Not here"), gmx::InternalError);
    EXPECT_THROW(gmx_impl("Not yet"), gmx::NotImplementedError);
    EXPECT_THROW(gmx_mem("No more"), gmx::InternalError);
    EXPECT_THROW(gmx_comm("No answer"), gmx::InternalError);
    EXPECT_EQ("Sizes differ", messageOf([]() { gmx_incons("Sizes differ"); }));
    int       value = 7;
    EXPECT_THROW(range_check(value, 0, 5), gmx::RangeError);
    EXPECT_NO_THROW(range_check(value, 0, 8));
    EXPECT_NE(std::string::npos, messageOf([value]() { range_check(value, 0, 5); }).find("has value 7"));
}

// Test: A failed assertion throws, with the condition and the message
TEST(FatalErrorTest, AssertionThrows)
{
    const gmx::ThrowOnFatalError library;
    int                          n = 3;
    EXPECT_NO_THROW(GMX_RELEASE_ASSERT(n == 3, "Three it is"));
    EXPECT_THROW(GMX_RELEASE_ASSERT(n == 4, "The number of atoms must be four"), gmx::InternalError);
    const auto text = messageOf([n]() { GMX_RELEASE_ASSERT(n == 4, "The number of atoms must be four"); });
    EXPECT_NE(std::string::npos, text.find("The number of atoms must be four")) << text;
    EXPECT_NE(std::string::npos, text.find("n == 4")) << text;
}

// Test: An exception that a program would end on is passed on as it is
TEST(FatalErrorTest, CaughtExceptionIsPassedOn)
{
    const gmx::ThrowOnFatalError library;
    const auto endsProgram = []()
    {
        try
        {
            GMX_THROW(gmx::FileIOError("Cannot read ff.xml"));
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    };
    EXPECT_THROW(endsProgram(), gmx::FileIOError);
    EXPECT_EQ("Cannot read ff.xml", messageOf(endsProgram));
}

// Test: A fatal message of ACT throws with the message itself, such that
// the caller knows what went wrong without the output of the program
TEST(FatalErrorTest, MessageHandlerThrowsItsMessage)
{
    const gmx::ThrowOnFatalError library;
    MsgHandler                   msghandler;
    EXPECT_THROW(msghandler.fatal("The force field has no particle types"), gmx::InvalidInputError);
    MsgHandler again;
    const auto text = messageOf([&again]() { again.fatal("The force field has no particle types"); });
    EXPECT_NE(std::string::npos, text.find("The force field has no particle types")) << text;
    EXPECT_FALSE(again.ok());
}

} // namespace

} // namespace alexandria
