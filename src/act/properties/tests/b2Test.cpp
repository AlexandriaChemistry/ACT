/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022,2023
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
#include <cmath>

#include <map>

#include <gtest/gtest.h>

#include "act/properties/b2data.h"
#include "act/properties/secondvirial.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{


namespace
{

class SecondVirialTest : public gmx::test::CommandLineTestBase
{
protected:
    gmx::test::TestReferenceChecker checker_;
    
    SecondVirialTest () : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
        //ForceField *mypd     = getForceField("ACM-g");
    }
    
    void testB2(bool LJ,
                double sig, double eps, double mass, 
                const std::vector<double> &Temperature)
    {
        gmx_output_env_t       *oenv    = nullptr;
        std::vector<gmx::RVec>  inertia = { { 0, 0, 0 }, { 0, 0, 0 } };
        
        output_env_init_default(&oenv);
 
        // First data point at this distance
        double x0    = 0.1;
        double bw    = 0.0001;
        double rmax  = 10;
        // We will add data points until rmax+x0
        int    irmax = (x0+rmax)/bw;
        B2Data b2data(irmax, bw, Temperature);
        // We start computing from x0, however we need the points from zero
        // for integration later. The integer index times the binwidth is
        // used as the distance in the integration algorithm.
        for(int i = int(x0/bw); i < irmax; i++)
        {
            if (LJ)
            {
                double x = bw*i;
                double y = 4*eps*(std::pow(sig/x, 12) - std::pow(sig/x, 6));
                double f = 4*(eps/sig)*(12*std::pow(sig/x, 13) - 6*std::pow(sig/x, 7));
                gmx::RVec ttt = { 0,0,0 };
                for(size_t iTemp = 0; iTemp < Temperature.size(); iTemp++)
                {
                    double beta  = 1.0/(BOLTZ*Temperature[iTemp]);
                    double g0_12 = std::exp(-beta*y);
                    b2data.addData(iTemp, i, g0_12-1,
                                   g0_12*f*f, g0_12*f*f, ttt, ttt);
                }
            }
        }
        if (LJ)
        {
            std::vector<double> m2 = { mass, mass };
            std::map<b2Type, std::vector<double>> b2t = {
                { b2Type::Classical, {} },
                { b2Type::Force, {} },
                { b2Type::Torque1, {} },
                { b2Type::Torque2, {} },
                { b2Type::Total, {} }
            };
            // Conversion to regular units cm^3/mol.
            double fac  = AVOGADRO*1e-21;
            for(size_t iTemp = 0; iTemp < Temperature.size(); iTemp++)
            {
                double Bclass, BqmForce, BqmTorque1, BqmTorque2;
                double beta = 1.0/(BOLTZ*Temperature[iTemp]);
                b2data.fillToXmin(iTemp, x0, bw);
                b2data.integrate(iTemp, bw, beta, m2, inertia,
                                 &Bclass, &BqmForce,
                                 &BqmTorque1, &BqmTorque2);
                b2t[b2Type::Classical].push_back(Bclass*fac);
                b2t[b2Type::Force].push_back(BqmForce*fac);
                b2t[b2Type::Torque1].push_back(BqmTorque1*fac);
                b2t[b2Type::Torque2].push_back(BqmTorque2*fac);
                b2t[b2Type::Total].push_back((Bclass+BqmForce+0.5*(BqmTorque1+BqmTorque2))*fac);
            };
            
            for(const auto &b2 : b2Type2str)
            {
                checker_.checkSequence(b2t[b2.first].begin(),
                                       b2t[b2.first].end(),
                                       b2.second.c_str());
            }
        }
    }
    
    static void TearDownTestCase()
    {
    }
};

TEST_F (SecondVirialTest, LennardJonesArgon) 
{
    double                        Targon  = 119.8;
    const std::vector<double>    &Temperature = { 80, 90, 100, 110, Targon };
    double                        mass = 39.948; // Argon
    // Numbers from Hirschfelder, Curtiss & Bird
    double sig = 0.3405;
    double eps = Targon*BOLTZ;
    testB2(true, sig, eps, mass, Temperature);
}

TEST_F (SecondVirialTest, LennardJonesXenon) 
{
    double                        Txenon  = 221;
    const std::vector<double>    &Temperature = { 298.2, 348.2, 373.2, 498.2, 548.3, 573.3 };
    double                        mass = 131.293; // Xenon
    // Numbers from Hirschfelder, Curtiss & Bird
    double sig = 0.41;
    double eps = Txenon*BOLTZ;
    // Results tabulated are for each of the temperatures (p. 167)
    // -128.4, -94,5, -81.5, -38.9, -28.0, -23.4
    // Those numbers were computed manually but our values are very close indeed.
    // In addition, the tabulated values did not include quantum corrections.
    testB2(true, sig, eps, mass, Temperature);
}

}

} // namespace
