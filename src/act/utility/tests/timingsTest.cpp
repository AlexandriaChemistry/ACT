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
 * Tests of the timings, the counting of floating point operations and
 * the report of the floating point performance.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <chrono>
#include <cstdlib>
#include <string>

#include <gtest/gtest.h>

#include "act/forces/simd_information.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/jsontree.h"
#include "act/utility/timings.h"

namespace alexandria
{

namespace
{

/*! \brief Find an object in a JsonTree
 * \param[in] tree The tree
 * \param[in] key  The key of the object
 * \return the object or nullptr if not found
 */
const JsonTree *findObject(const JsonTree &tree, const std::string &key)
{
    for (const auto &obj : tree.objects())
    {
        if (obj.key() == key)
        {
            return &obj;
        }
    }
    return nullptr;
}

/*! \brief The number in a JsonTree, or in its value if it has a unit
 * \param[in] tree The tree
 * \param[in] key  The key of the number
 * \return the number
 */
double findNumber(const JsonTree &tree, const std::string &key)
{
    auto obj = findObject(tree, key);
    EXPECT_NE(nullptr, obj) << key;
    if (nullptr == obj)
    {
        return 0;
    }
    auto value = findObject(*obj, "value");
    return std::stod(value ? value->value() : obj->value());
}

//! Spend some time, such that a region has a measurable duration
void spin()
{
    auto start = std::chrono::steady_clock::now();
    while (std::chrono::steady_clock::now() - start < std::chrono::milliseconds(5))
    {
    }
}

// Test: Operations count for the innermost active region, or outside of regions
TEST(TimingsTest, AddFlopsToInnermostRegion)
{
    auto &timings = Timings::instance();
    timings.reset();
    auto outer = timings.index("Test outer");
    auto inner = timings.index("Test inner");
    timings.addFlops(2);
    timings.start(outer);
    timings.addFlops(10);
    timings.start(inner);
    timings.addFlops(5);
    timings.stop(inner);
    timings.addFlops(1);
    timings.stop(outer);
    EXPECT_EQ(11, timings.flops()[outer]);
    EXPECT_EQ(5, timings.flops()[inner]);
    EXPECT_EQ(2, timings.untimedFlops());
    timings.reset();
    EXPECT_EQ(0, timings.flops()[outer]);
    EXPECT_EQ(0, timings.flops()[inner]);
    EXPECT_EQ(0, timings.untimedFlops());
}

// Test: The performance is reported per core on average and in the
// regions with operations, and as a fraction of the peak
TEST(TimingsTest, PerformanceInJson)
{
    auto &timings = Timings::instance();
    timings.reset();
    auto kernel = timings.index("Test kernel");
    timings.start(kernel);
    timings.addFlops(3e9);
    spin();
    timings.stop(kernel);
    double seconds = timings.seconds()[kernel];
    ASSERT_GT(seconds, 0);

    CommunicationRecord cr;
    cr.init(cr.size());
    JsonTree        jtree("Test");
    PeakPerformance peak = { 100, "test processor" };
    double          wall = 2;
    printTimings(&cr, nullptr, wall, &jtree, &peak);
    auto jtimings = findObject(jtree, "Timings");
    ASSERT_NE(nullptr, jtimings);
    auto jperf = findObject(*jtimings, "Performance");
    ASSERT_NE(nullptr, jperf);
    EXPECT_NEAR(3e9, findNumber(*jperf, "Floating point operations"), 1e3);
    // 3e9 operations in 2 s
    EXPECT_NEAR(1.5, findNumber(*jperf, "Average per core"), 1e-4);
    double inKernel = 3/seconds;
    EXPECT_NEAR(inKernel, findNumber(*jperf, "Per core in the counted parts"), 1e-3*inKernel);
    EXPECT_NEAR(100, findNumber(*jperf, "Theoretical peak per core"), 1e-4);
    EXPECT_EQ("test processor", findObject(*jperf, "Theoretical peak basis")->value());
    // The peak is for double precision, as are all calculations in a double build
    EXPECT_EQ("double", findObject(*jperf, "Theoretical peak precision")->value());
#if GMX_DOUBLE
    EXPECT_EQ("double", findObject(*jperf, "Precision")->value());
#endif
    EXPECT_NEAR(1.5, findNumber(*jperf, "Fraction of the peak on average"), 1e-3);
    EXPECT_NEAR(inKernel, findNumber(*jperf, "Fraction of the peak in the counted parts"), 1e-3*inKernel);
    // The rate of each part
    auto jparts = findObject(*jtimings, "Parts");
    ASSERT_NE(nullptr, jparts);
    auto jkernel = findObject(*jparts, "Test kernel");
    ASSERT_NE(nullptr, jkernel);
    EXPECT_NEAR(inKernel, findNumber(*jkernel, "GFlop/s"), 1e-3*inKernel);
    timings.reset();
}

// Test: Without counted operations there is no performance report, and
// without a peak there is no fraction of it
TEST(TimingsTest, PerformanceOnlyIfCounted)
{
    auto &timings = Timings::instance();
    CommunicationRecord cr;
    cr.init(cr.size());
    auto region = timings.index("Test region");
    timings.reset();
    timings.start(region);
    timings.stop(region);
    JsonTree none("Test");
    printTimings(&cr, nullptr, 1, &none);
    EXPECT_EQ(nullptr, findObject(*findObject(none, "Timings"), "Performance"));

    timings.reset();
    timings.start(region);
    timings.addFlops(1e6);
    timings.stop(region);
    JsonTree nopeak("Test");
    printTimings(&cr, nullptr, 1, &nopeak);
    auto jperf = findObject(*findObject(nopeak, "Timings"), "Performance");
    ASSERT_NE(nullptr, jperf);
    EXPECT_NE(nullptr, findObject(*jperf, "Average per core"));
    EXPECT_EQ(nullptr, findObject(*jperf, "Theoretical peak per core"));
    timings.reset();
}

// Test: The peak can be given in the environment, otherwise it is that of
// the processor if it is known
TEST(PeakPerformanceTest, Environment)
{
    const char *name  = "ACT_PEAK_GFLOPS";
    const char *saved = std::getenv(name);
    std::string old   = saved ? saved : "";
    setenv(name, "12.5", 1);
    auto peak = theoreticalPeakPerCore();
    EXPECT_EQ(12.5, peak.gflops);
    EXPECT_EQ("given by ACT_PEAK_GFLOPS", peak.basis);
    unsetenv(name);
    auto detected = theoreticalPeakPerCore();
    EXPECT_FALSE(detected.basis.empty());
    if (detected.gflops > 0)
    {
        EXPECT_NE(std::string::npos, detected.basis.find("flop/cycle")) << detected.basis;
    }
    if (saved)
    {
        setenv(name, old.c_str(), 1);
    }
}

} // namespace

} // namespace alexandria
