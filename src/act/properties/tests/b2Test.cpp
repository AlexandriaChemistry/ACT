/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022-2026
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
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/basics/msg_handler.h"
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
        B2Data b2data(irmax, Temperature);
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
                    b2data.addData(iTemp, i, x, g0_12-1,
                                   g0_12*f*f, g0_12*f*f, ttt, ttt);
                }
            }
        }
        MsgHandler msghandler;
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
                b2data.integrate(&msghandler,
                                 iTemp, beta, m2, inertia,
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

//! Results of the B2Data integration
struct B2Result
{
    //! Classical second virial (nm^3)
    double classical = 0;
    //! Translational quantum correction (nm^3)
    double force     = 0;
    //! Rotational quantum correction for molecule 1 (nm^3)
    double torque1   = 0;
    //! Rotational quantum correction for molecule 2 (nm^3)
    double torque2   = 0;
};

/*! \brief Integrate synthetic data that is constant up to rmax
 * \param[in] mayer   Mayer function value for r < sigma, zero beyond
 * \param[in] sigma   Range of the Mayer function (nm)
 * \param[in] F2      Weighted square force on each molecule
 * \param[in] tau2    Weighted square torque on each molecule
 * \param[in] mass    Mass of each molecule
 * \param[in] inertia Moments of inertia of each molecule
 * \param[in] beta    1/kT
 * \return The integrated components
 */
B2Result integrateSynthetic(double                        mayer,
                            double                        sigma,
                            const double                  F2[2],
                            const gmx::RVec               tau2[2],
                            const std::vector<double>    &mass,
                            const std::vector<gmx::RVec> &inertia,
                            double                        beta)
{
    // Bins of 1 pm up to 1 nm
    int    nbins = 1001;
    double bw    = 0.001;
    B2Data b2data(nbins, { 1.0/(BOLTZ*beta) });
    for(int i = 0; i < nbins; i++)
    {
        double r = i*bw;
        b2data.addData(0, i, r, r < sigma ? mayer : 0, F2[0], F2[1], tau2[0], tau2[1]);
    }
    MsgHandler msghandler;
    B2Result   res;
    b2data.integrate(&msghandler, 0, beta, mass, inertia,
                     &res.classical, &res.force, &res.torque1, &res.torque2);
    EXPECT_TRUE(msghandler.ok());
    return res;
}

//! \return the prefactor of the quantum corrections, beta^3 hbar^2/24
double hbarFactor(double beta)
{
    return beta*gmx::square(beta*PLANCK/(2*M_PI))/24;
}

//! \return 4 pi R^3/3 for R = 1 nm, the integral of a constant over the grid
double sphereVolume()
{
    return 4*M_PI/3;
}

//! Relative tolerance for the numerical integration of a constant
const double integrationTolerance = 1e-5;

// Test: The classical second virial of hard spheres is 2 pi sigma^3/3
TEST(B2DataTest, HardSpheres)
{
    double    F2[2]      = { 0, 0 };
    gmx::RVec tau2[2]    = { { 0, 0, 0 }, { 0, 0, 0 } };
    double    sigma      = 0.3;
    auto      res        = integrateSynthetic(-1, sigma, F2, tau2, { 18, 18 },
                                              { { 1, 1, 1 }, { 1, 1, 1 } }, 1.0/(BOLTZ*300));
    EXPECT_NEAR(2*M_PI*std::pow(sigma, 3)/3, res.classical, 1e-2*2*M_PI*std::pow(sigma, 3)/3);
    EXPECT_EQ(0, res.force);
    EXPECT_EQ(0, res.torque1);
    EXPECT_EQ(0, res.torque2);
}

// Test: The translational quantum correction, Schenter Eqn. 9
TEST(B2DataTest, ForceCorrection)
{
    double    beta    = 1.0/(BOLTZ*300);
    double    F2[2]   = { 100, 300 };
    gmx::RVec tau2[2] = { { 0, 0, 0 }, { 0, 0, 0 } };
    std::vector<double> mass = { 18, 36 };
    auto res      = integrateSynthetic(0, 0, F2, tau2, mass,
                                       { { 1, 1, 1 }, { 1, 1, 1 } }, beta);
    double expect = 0.5*hbarFactor(beta)*sphereVolume()*(F2[0]/mass[0] + F2[1]/mass[1]);
    EXPECT_NEAR(expect, res.force, integrationTolerance*expect);
}

// Test: The rotational quantum corrections, Gray and Gubbins Eqn. 3.282
TEST(B2DataTest, TorqueCorrection)
{
    double    beta    = 1.0/(BOLTZ*300);
    double    F2[2]   = { 0, 0 };
    gmx::RVec tau2    = { 50, 70, 90 };
    gmx::RVec tau[2]  = { tau2, tau2 };
    gmx::RVec I       = { 0.006, 0.012, 0.018 };
    // Identical molecules give identical torque corrections
    auto res     = integrateSynthetic(0, 0, F2, tau, { 18, 18 }, { I, I }, beta);
    double expect = hbarFactor(beta)*sphereVolume()*(tau2[XX]/I[XX] + tau2[YY]/I[YY] + tau2[ZZ]/I[ZZ]);
    EXPECT_NEAR(expect, res.torque1, integrationTolerance*expect);
    EXPECT_DOUBLE_EQ(res.torque1, res.torque2);
    // Doubling the moments of inertia of molecule 2 halves its correction
    gmx::RVec I2  = { 2*I[XX], 2*I[YY], 2*I[ZZ] };
    auto res2     = integrateSynthetic(0, 0, F2, tau, { 18, 18 }, { I, I2 }, beta);
    EXPECT_DOUBLE_EQ(res.torque1, res2.torque1);
    EXPECT_NEAR(0.5*res.torque2, res2.torque2, 1e-12*res.torque2);
}

// Test: Axes with zero moment of inertia (linear molecules) do not contribute
TEST(B2DataTest, LinearMolecule)
{
    double    beta    = 1.0/(BOLTZ*300);
    double    F2[2]   = { 0, 0 };
    gmx::RVec tau2    = { 50, 70, 90 };
    gmx::RVec tau[2]  = { tau2, tau2 };
    gmx::RVec I       = { 0, 0.012, 0.012 };
    auto res      = integrateSynthetic(0, 0, F2, tau, { 44, 44 }, { I, I }, beta);
    double expect = hbarFactor(beta)*sphereVolume()*(tau2[YY] + tau2[ZZ])/I[YY];
    EXPECT_TRUE(std::isfinite(res.torque1));
    EXPECT_NEAR(expect, res.torque1, integrationTolerance*expect);
}

}

} // namespace
