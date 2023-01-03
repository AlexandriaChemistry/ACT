/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022
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
#include <math.h>

#include <map>

#include <gtest/gtest.h>

#include "alexandria/secondvirial.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "gromacs/math/units.h"
#include "gromacs/utility/stringutil.h"

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
        //Poldata *mypd     = getPoldata("ACM-g");
    }
    
    void testB2(bool LJ,
                double sig, double eps, double mass, 
                const std::vector<double> &Temperature)
    {
        FILE                                *logFile = nullptr;
        gmx_stats                            edist;
        gmx_output_env_t                    *oenv    = nullptr;
        gmx::RVec                            inertia[2] = { { 0, 0, 0 }, { 0, 0, 0 } };
        std::vector<gmx::RVec>               force1;
        std::vector<std::vector<gmx::RVec>>  torque1;
        
        output_env_init_default(&oenv);
        
        double x0 = 0.15;
        for(int i = 0; i < 2581; i++)
        {
            if (LJ)
            {
                double x = x0 + 0.002*i;
                double y = 4*eps*(std::pow(sig/x, 12) - std::pow(sig/x, 6));
                double f = 4*eps*(12*std::pow(sig/x, 13) - 6*std::pow(sig/x, 7));
                edist.add_point(x, y, 0, 0);
                force1.push_back({ 0, 0, f });
                std::vector<gmx::RVec> ttt;
                // Torque is needed for two compounds
                ttt.push_back({ 0, 0, 0 });
                torque1.push_back(ttt);
                torque1.push_back(ttt);
            }
        }
        if (LJ)
        {
            ReRunner rerun;
            std::vector<t_filenm> fnm;
            rerun.setTemperatures(Temperature);
            rerun.computeB2(logFile, edist,
                            mass, inertia, force1, torque1, fnm);
            auto b2t = rerun.b2Temp();
            checker_.checkSequence(b2t.begin(), b2t.end(), "B2(T)");
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
