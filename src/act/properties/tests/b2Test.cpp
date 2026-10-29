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
#include "act/properties/dimergenerator.h"
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
        int    i0    = std::round(x0/bw);
        // The bins start at x0, as for the dimers of b2, below the Mayer
        // function is -1.
        B2Data b2data(irmax - i0, Temperature);
        for(int i = i0; i < irmax; i++)
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
                    b2data.addData(iTemp, i - i0, x, g0_12-1,
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


// Test: Data points count with their weight, as for a grid of orientations
TEST(B2DataTest, WeightedAverage)
{
    double    beta    = 1.0/(BOLTZ*300);
    gmx::RVec tau2    = { 50, 70, 90 };
    gmx::RVec I       = { 0.006, 0.012, 0.018 };
    std::vector<double>    mass    = { 18, 18 };
    std::vector<gmx::RVec> inertia = { I, I };
    // A Mayer function of -1 with weight 3 and of 0 with weight 1 below
    // 0.3 nm, the same as -0.75 with weight one.
    int    nbins = 1001;
    double bw    = 0.001;
    B2Data weighted(nbins, { 300 }), single(nbins, { 300 });
    gmx::RVec zero = { 0, 0, 0 };
    for(int i = 0; i < nbins; i++)
    {
        double r = i*bw;
        weighted.addData(0, i, r, r < 0.3 ? -1 : 0, 10, 10, tau2, tau2, 3);
        weighted.addData(0, i, r, 0, 2, 2, zero, zero, 1);
        single.addData(0, i, r, r < 0.3 ? -0.75 : 0, 8, 8, { 37.5, 52.5, 67.5 }, { 37.5, 52.5, 67.5 });
    }
    MsgHandler msghandler;
    B2Result   rw, rs;
    weighted.integrate(&msghandler, 0, beta, mass, inertia, &rw.classical, &rw.force, &rw.torque1, &rw.torque2);
    single.integrate(&msghandler, 0, beta, mass, inertia, &rs.classical, &rs.force, &rs.torque1, &rs.torque2);
    EXPECT_TRUE(msghandler.ok());
    EXPECT_NEAR(rs.classical, rw.classical, 1e-12*std::abs(rs.classical));
    EXPECT_NEAR(rs.force, rw.force, 1e-12*rs.force);
    EXPECT_NEAR(rs.torque1, rw.torque1, 1e-12*rs.torque1);
    EXPECT_NEAR(rs.torque2, rw.torque2, 1e-12*rs.torque2);
}

// Test: The grid of orientations as in the program Orient. Gauss-Legendre
// integration over cos(beta) with n points is exact for polynomials up
// to order 2n-1, and the equally spaced angles average cosines to zero.
TEST(GaussGridTest, Orientations)
{
    int nalpha = 3, nbeta = 4, ngamma = 2;
    DimerGenerator gendimers;
    gendimers.setGaussGrid(nalpha, nbeta, ngamma);
    EXPECT_EQ(RandAlg::Gauss, gendimers.randAlg());
    int nori = nalpha*nbeta*nbeta*ngamma*ngamma;
    EXPECT_EQ(nori, gendimers.numberOfOrientations(1024));
    gendimers.generateRandomNumbers(1024);
    const auto &q = gendimers.allRandom();
    ASSERT_EQ(static_cast<size_t>(nori), q.size());
    double wsum = 0, c2 = 0, c4 = 0, c22 = 0, cg1 = 0, ca2 = 0;
    for(const auto &qq : q)
    {
        ASSERT_EQ(7u, qq.size());
        for(int j = 0; j < 6; j++)
        {
            EXPECT_GE(qq[j], 0);
            EXPECT_LT(qq[j], 1);
        }
        // Rotating both molecules about the line between them does not
        // change the dimer, alpha of the first molecule is zero
        EXPECT_EQ(0, qq[0]);
        double w  = qq[6];
        double b1 = 2*qq[1] - 1;
        double b2 = 2*qq[4] - 1;
        wsum += w;
        c2   += w*b1*b1;
        c4   += w*b2*b2*b2*b2;
        c22  += w*b1*b1*b2*b2;
        cg1  += w*std::cos(2*M_PI*qq[2]);
        ca2  += w*std::cos(2*M_PI*qq[3]);
    }
    // The Gauss-Legendre weights add up to two for each molecule
    EXPECT_NEAR(4.0*nalpha*ngamma*ngamma, wsum, 1e-12*wsum);
    EXPECT_NEAR(1.0/3, c2/wsum, 1e-12);
    EXPECT_NEAR(1.0/5, c4/wsum, 1e-12);
    EXPECT_NEAR(1.0/9, c22/wsum, 1e-12);
    EXPECT_NEAR(0, cg1/wsum, 1e-12);
    EXPECT_NEAR(0, ca2/wsum, 1e-12);
}

}

} // namespace
