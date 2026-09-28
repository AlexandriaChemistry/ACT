/*
 * This source file is part of the Alexandria program.
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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#include "../simpson.h"

#include <cmath>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

TEST(SimpsonTest, OddVectorExp)
{
    std::vector<double> x = { 0, 1, 2, 3, 4, 5, 6 };
    std::vector<double> y;
    for(size_t i = 0; i < x.size(); i++)
    {
        y.push_back(std::exp(-x[i]));
    }
    alexandria::MsgHandler mh;
    auto integral = simpsonIntegrate(&mh, x, y);
    EXPECT_TRUE(mh.ok());
    EXPECT_NEAR(integral, 1+std::exp(-6.0), 1e-4);
}

TEST(SimpsonTest, OddVectorExpOffset)
{
    std::vector<double> x = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
    std::vector<double> y;
    double offset = std::exp(-x[x.size()-1]);
    for(size_t i = 0; i < x.size(); i++)
    {
        y.push_back(std::exp(-x[i])-offset);
    }
    alexandria::MsgHandler mh;
    auto integral = simpsonIntegrate(&mh, x, y);
    EXPECT_TRUE(mh.ok());
    double x0     = x[x.size()-1];
    double result = 1 - std::exp(-x0) - x0*std::exp(-x0);
    EXPECT_NEAR(integral, result, 1e-4);
}

TEST(SimpsonTest, OddVectorInvSquare)
{
    std::vector<double> x;
    for(size_t i = 0; i < 81; i++)
    {
        x.push_back(1+0.1*i);
    }
    std::vector<double> y;
    for(size_t i = 0; i < x.size(); i++)
    {
        y.push_back(1.0/(x[i]));
    }
    alexandria::MsgHandler mh;
    auto integral = simpsonIntegrate(&mh, x, y);
    EXPECT_TRUE(mh.ok());
    EXPECT_NEAR(integral, std::log(9.0), 1e-4);
}


TEST(SimpsonTest, OddVectorCosine)
{
    std::vector<double> x, y;
    for(size_t i = 0; i < 91; i++)
    {
        double xx = i/90.0;
        x.push_back(xx);
        y.push_back(1.0+std::cos(M_PI*xx));
    }
    alexandria::MsgHandler mh;
    auto integral = simpsonIntegrate(&mh, x, y);
    EXPECT_TRUE(mh.ok());
    EXPECT_NEAR(integral, 1, 1e-4);
}

TEST(SimpsonTest, EvenVectorCosine)
{
    std::vector<double> x, y;
    for(size_t i = 0; i < 90; i++)
    {
        double xx = i/90.0;
        x.push_back(xx);
        y.push_back(1.0+std::cos(M_PI*xx));
    }
    alexandria::MsgHandler mh;
    auto integral = simpsonIntegrate(&mh, x, y, false);
    EXPECT_FALSE(mh.ok());
    mh.resetStatus();
    integral = simpsonIntegrate(&mh, x, y, true);
    EXPECT_TRUE(mh.ok());
    EXPECT_NEAR(integral, 1, 1e-4);
}


