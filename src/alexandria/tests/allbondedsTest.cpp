/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2021
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <math.h>

#include <map>

#include <gtest/gtest.h>

#include "alexandria/allbondeds.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "act/poldata/poldata_utils.h"

namespace alexandria
{

namespace
{

class AllBondedsTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
 
        AllBondedsTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
            //Poldata *mypd     = getPoldata("ACM-g");
        }

        void runOBstats(const std::vector<double> &xx)
        {
            Identifier atom("h");
            OneBonded ob(atom);
            for(const auto &x : xx)
            {
                ob.addPoint(x);
            }
            real   average, sigma;
            size_t N;
            eStats ok = ob.getAverageSigmaN(&average, &sigma, &N);
            if (eStats::OK == ok)
            {
                checker_.checkReal(average, "average");
                checker_.checkReal(sigma, "sigma");
                checker_.checkInteger(static_cast<int>(N), "N");
            }
            else
            {
                checker_.checkString(gmx_stats_message(ok), "Error");
            }
            
        }
        static void TearDownTestCase()
        {
        }

};

TEST_F (AllBondedsTest, OneBondedStats){
    std::vector<double> xx = { 103, 104, 105, 106.5, 107, 109, 106, 104, 107 };
    runOBstats(xx);
}

TEST_F (AllBondedsTest, OneBondedStatsEmpty){
    std::vector<double> xx;
    runOBstats(xx);
}

}

} // namespace
