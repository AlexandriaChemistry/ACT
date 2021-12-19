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

#include "actpre.h"

#include <math.h>

#include <cstdlib>

#include <gtest/gtest.h>

#include "gromacs/statistics/statistics.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace 
{

class StatisticsTest : public gmx::test::CommandLineTestBase
{
    std::vector<double> x, y, dx, dy;
protected:
    StatisticsTest()
    {
        x  = {  3 , 5,  8, 11, 13, 17, 21, 4,   6,  9, 12 };
        y  = { 13, 51, 18, 21, 16, 27, 26, 40, 36, 19, 42 };
        dx = {  1,  2,  3,  1,  1,1.3,  1,  2,  2,  1,  2 };
        dy = {  1,  2,1.3,  1,  2,  3,  1,  2,  3,  1,  2 };
    }

    //! Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    void analyze(gmx_stats *gs)
    {
        gmx::test::TestReferenceChecker myCheck(this->rootChecker());
        real a, b, da, db, chi2, Rfit;
        gs->get_ab(elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &Rfit);
        myCheck.checkDouble(a, "a");
        myCheck.checkDouble(b, "b");
        myCheck.checkDouble(da, "da");
        myCheck.checkDouble(db, "db");
        myCheck.checkDouble(chi2, "chi2");
        myCheck.checkDouble(Rfit, "Rfit");
        real R, rmsd;
        gs->get_corr_coeff(&R);
        gs->get_rmsd(&rmsd);
        myCheck.checkDouble(R, "R");
        myCheck.checkDouble(rmsd, "rmsd");
        real mse, mae;
        gs->get_mse_mae(&mse, &mae);
        myCheck.checkDouble(mae, "mean absolute error");
        myCheck.checkDouble(mse, "mean signed error");
    }
    
    //! Do the actual testing
    void testXYDXDY ()
    {
        gmx_stats gs;
        for(size_t i = 0; i < x.size(); i++)
        {
            gs.add_point(x[i], y[i], dx[i], dy[i]);
        }
        analyze(&gs);
    }
    void testYDY ()
    {
        gmx_stats gs;
        for(size_t i = 0; i < y.size(); i++)
        {
            gs.add_point_ydy(y[i], dy[i]);
        }
        analyze(&gs);
    }
    void testHisto(eHisto ehisto, bool normalized)
    {
        gmx::test::TestReferenceChecker myCheck(this->rootChecker());
        gmx_stats gs;
        gs.add_points(x.size(), x.data(), y.data(), dx.data(), dy.data());
        int N = gs.get_npoints();
        myCheck.checkInteger(N, "Number of points");
        std::vector<double> xx, yy;
        int nbins=0;
        eStats ok = gs.make_histogram(1.0, &nbins, ehisto, normalized, &xx, &yy);
        myCheck.checkString(gmx_stats_message(ok), "Histogram");
        myCheck.checkInteger(nbins, "nbins");
        myCheck.checkSequenceArray(xx.size(), xx.data(), "x-histo");
        myCheck.checkSequenceArray(yy.size(), yy.data(), "y-histo");
    }
    //! Clean the test data.
    static void TearDownTestCase()
    {
    }
    
};

TEST_F (StatisticsTest, XYDXDY){
    testXYDXDY();
}

TEST_F (StatisticsTest, YDY){
    testYDY();
}

TEST_F (StatisticsTest, HistoX){
    testHisto(eHisto::X, false);
}

TEST_F (StatisticsTest, HistoXNormalized){
    testHisto(eHisto::X, true);
}

TEST_F (StatisticsTest, HistoY){
    testHisto(eHisto::Y, false);
}

TEST_F (StatisticsTest, HistoYNormalized){
    testHisto(eHisto::Y, true);
}

} // anonymous
