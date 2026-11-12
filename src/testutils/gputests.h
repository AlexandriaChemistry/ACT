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
/*! \libinternal \file
 * \brief
 * Skipping the tests that need a GPU, such that CTest reports them as skipped.
 *
 * The tests that need a GPU have names that end in Gpu. CTest runs them
 * separately from the other tests of a test program, see
 * act_add_gpu_unit_test in TestMacros.cmake, and reports that entry as
 * skipped if the output contains gmx::test::c_gpuTestsSkipped. Since CTest
 * then ignores whether tests failed, the text may only be printed for a
 * reason that holds for all GPU tests of the test program alike, such as
 * the absence of a GPU. If a GPU test cannot run for a reason of its own,
 * for instance since the GPU cannot compute its molecules, it fails.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef ACT_TESTUTILS_GPUTESTS_H
#define ACT_TESTUTILS_GPUTESTS_H

#include <cstdio>

#include <gtest/gtest.h>

namespace gmx
{

namespace test
{

//! The text by which CTest knows that the GPU tests were skipped, see act_add_gpu_unit_test
static const char *const c_gpuTestsSkipped = "ACT GPU tests skipped";

/*! \brief Whether a test that needs a GPU can run, if not it says so
 *
 * Without googletest's GTEST_SKIP (googletest 1.8), the test then returns
 * without testing anything, and CTest reports its entry as skipped.
 * \param[in] canRun Whether the test can run
 * \param[in] why    Why not, for all GPU tests of the test program alike
 * \return canRun
 */
static inline bool gpuTestRuns(bool canRun, const char *why)
{
    if (!canRun)
    {
        printf("%s: %s: %s\n", c_gpuTestsSkipped,
               ::testing::UnitTest::GetInstance()->current_test_info()->name(), why);
        fflush(stdout);
    }
    return canRun;
}

} // namespace test

} // namespace gmx

#endif
