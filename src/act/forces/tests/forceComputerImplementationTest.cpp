/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2025
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
/*! \internal \file
 * \brief
 * Implements test of bonded force routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_listed-forces
 */
#include "actpre.h"

#include "../forcecomputerimpl.h"

#include <cmath>

#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "act/forcefield/forcefield_parametername.h"
#include "act/molprop/topologyentry.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

//! Number of atoms used in these tests.
#define NATOMS 4

class ForceComputerImplementationTest : public ::testing::Test
{
protected:
    gmx::test::TestReferenceData    refData_;
    gmx::test::TestReferenceChecker checker_;
    std::vector<ActAtom>            atoms_;

    ForceComputerImplementationTest( ) :
        checker_(refData_.rootChecker())
    {
        gmx::test::FloatingPointTolerance tolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-6));
        checker_.setDefaultTolerance(tolerance);
        std::vector<double> q = { 1, 1, -1, -2 };
        for (int i = 0; i < 4; i++)
        {
            // Note the row value must be <= 2 in order to work for Slaters
            atoms_.push_back({ "C", "C", "C", ActParticle::Atom, 6, 12, q[i], 2 });
        }
    }

    void testPot(Potential                     p,
                 const TopologyEntryVector    &top,
                 const std::vector<gmx::RVec> *coordinates)
    {
        // Generate force array with zeroes
        gmx::RVec                         fzero = { 0, 0, 0 };
        std::vector<gmx::RVec>            forces(coordinates->size(), fzero);
        // Energy map
        std::map<InteractionType, double> energies;

        // The correct force computer
        auto bfc = getBondForceComputer(p);

        auto epot = bfc(top, atoms_, coordinates, &forces, &energies);
        checker_.checkReal(epot, "Epot");
        checker_.checkSequence(forces.begin(), forces.end(), "Forces");
        for (const auto &e : energies)
        {
            checker_.checkReal(e.second, interactionTypeToString(e.first).c_str());
        }
        // Now check whether force is negative derivative of the potential
        // only for potentials with two atoms
        if (coordinates->size() != 2)
        {
            return;
        }
        double dx = 1e-6;
        double ener[2];
        for (int k = 0; k < 2; k++)
        {
            // Additional force and energy structures
            std::vector<gmx::RVec>            forces(coordinates->size(), fzero);
            std::map<InteractionType, double> energies;
            // Modify the coordinates
            auto newx = *coordinates;
            newx[1][ZZ] += (2*k-1)*dx;
            ener[k] = bfc(top, atoms_, &newx, &forces, &energies);
        }
        // Compute numerical derivative
        double fz = -(ener[1]-ener[0])/(2*dx);
        double toler = 5e-6;
        // It should be identical to the forces computed at the central point
        checker_.checkReal(fz, "fz");
        double diff = abs(fz - forces[1][ZZ]);
        if (diff >= toler)
        {
            fprintf(stderr, "Numeric force %g, analytical force %g, diff %g\n", fz, forces[1][ZZ], diff);
        }
        EXPECT_TRUE(diff < toler);
    }    
};

TEST_F (ForceComputerImplementationTest, LJ8_6)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[lj8_6SIGMA]   = 0.5;
    params[lj8_6EPSILON] = 0.25;
    top[0]->setParams(params);

    testPot(Potential::LJ8_6, top, &x);
}

TEST_F (ForceComputerImplementationTest, LJ12_6)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[lj12_6SIGMA]   = 0.5;
    params[lj12_6EPSILON] = 0.25;
    top[0]->setParams(params);

    testPot(Potential::LJ12_6, top, &x);
}

TEST_F (ForceComputerImplementationTest, LJ12_6_4)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(3);
    params[lj12_6_4SIGMA]   = 0.5;
    params[lj12_6_4EPSILON] = 0.25;
    params[lj12_6_4GAMMA]   = 0.5;
    top[0]->setParams(params);

    testPot(Potential::LJ12_6_4, top, &x);
}

TEST_F (ForceComputerImplementationTest, LJ14_7)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(4);
    params[lj14_7SIGMA]   = 0.5;
    params[lj14_7EPSILON] = 1;
    params[lj14_7GAMMA]   = 0.5;
    params[lj14_7DELTA]   = 0.1;
    top[0]->setParams(params);

    testPot(Potential::LJ14_7, top, &x);
}

TEST_F (ForceComputerImplementationTest, BUCKINGHAM)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(3);
    params[bhA]  = 5000;
    params[bhB]  = 20;
    params[bhC6] = 0.001;
    top[0]->setParams(params);

    testPot(Potential::BUCKINGHAM, top, &x);
}

TEST_F (ForceComputerImplementationTest, WANG_BUCKINGHAM)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(3);
    params[wbhSIGMA]   = 0.5;
    params[wbhEPSILON] = 0.25;
    params[wbhGAMMA]   = 10.0;
    top[0]->setParams(params);

    testPot(Potential::WANG_BUCKINGHAM, top, &x);
}

TEST_F (ForceComputerImplementationTest, GENERALIZED_BUCKINGHAM)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(4);
    params[gbhRMIN]    = 0.5;
    params[gbhEPSILON] = 0.25;
    params[gbhGAMMA]   = 0.5;
    params[gbhDELTA]   = 0.1;
    top[0]->setParams(params);

    testPot(Potential::GENERALIZED_BUCKINGHAM, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_POINT)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_POINT, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_GAUSSIAN)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_GAUSSIAN, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_GAUSSIAN_NEG)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    atoms_[1].setCharge(-1);
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_GAUSSIAN, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_GAUSSIAN_ZERO)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.0 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_GAUSSIAN, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_SLATER)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_SLATER, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_SLATER_ZETA0)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[coulZETA]  = 0;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_SLATER, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_SLATER_ZETA0_NEG)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    atoms_[1].setCharge(-1);
    std::vector<double> params(2);
    params[coulZETA]  = 0;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_SLATER, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_SLATER_CLOSE)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.2 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_SLATER, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_SLATER_CLOSE_NEG)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.2 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    atoms_[1].setCharge(-1);
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_SLATER, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_SLATER_ZERO)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0 },
        { 0, 0, 0 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_SLATER, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_SLATER_ZERO_NEG)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0 },
        { 0, 0, 0 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    atoms_[1].setCharge(-1);
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_SLATER, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_SLATER_NEG)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0    },
        { 0, 0, 0.5  }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    atoms_[1].setCharge(-1);
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_SLATER, top, &x);
}

TEST_F (ForceComputerImplementationTest, COULOMB_SLATER_NEG2)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0    },
        { 0, 0, 0.5  }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    atoms_[1].setCharge(-2);
    std::vector<double> params(2);
    params[coulZETA]  = 10;
    params[coulZETA2] = 6;
    top[0]->setParams(params);

    testPot(Potential::COULOMB_SLATER, top, &x);
}

TEST_F (ForceComputerImplementationTest, TT2bExch)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(6, 0.0);
    params[tt2bA]     = 1000;
    params[tt2bBexch] = 10;
    params[tt2bBdisp] = 10;
    top[0]->setParams(params);

    testPot(Potential::TT2b, top, &x);
}

TEST_F (ForceComputerImplementationTest, TT2bDispC6)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(6, 0.0);
    params[tt2bA]     = 0;
    params[tt2bBexch] = 10;
    params[tt2bBdisp] = 10;
    params[tt2bC6]    = 0.001;
    top[0]->setParams(params);

    testPot(Potential::TT2b, top, &x);
}

TEST_F (ForceComputerImplementationTest, TT2bDispC8)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(6, 0.0);
    params[tt2bA]     = 0;
    params[tt2bBexch] = 10;
    params[tt2bBdisp] = 10;
    params[tt2bC8]    = 0.0001;
    top[0]->setParams(params);

    testPot(Potential::TT2b, top, &x);
}

TEST_F (ForceComputerImplementationTest, TT2bDispC10)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(6, 0.0);
    params[tt2bA]     = 0;
    params[tt2bBexch] = 10;
    params[tt2bBdisp] = 10;
    params[tt2bC10]   = 0.001;
    top[0]->setParams(params);

    testPot(Potential::TT2b, top, &x);
}

TEST_F (ForceComputerImplementationTest, TT2bAll)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(6, 0.0);
    params[tt2bA]     = 1000;
    params[tt2bBexch] = 10;
    params[tt2bBdisp] = 20;
    params[tt2bC6]    = 0.001;
    params[tt2bC8]    = 0.001;
    params[tt2bC10]   = 0.001;
    top[0]->setParams(params);

    testPot(Potential::TT2b, top, &x);
}

TEST_F (ForceComputerImplementationTest, SLATER_ISA_TT)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(6, 0.0);
    params[tt2bA]     = 1000;
    params[tt2bBexch] = 10;
    params[tt2bBdisp] = 20;
    params[tt2bC6]    = 0.001;
    params[tt2bC8]    = 0.001;
    params[tt2bC10]   = 0.001;
    top[0]->setParams(params);

    testPot(Potential::SLATER_ISA_TT, top, &x);
}

TEST_F (ForceComputerImplementationTest, TANG_TOENNIES)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(5, 0.0);
    params[ttA]   = 1000;
    params[ttB]   = 10;
    params[ttC6]  = 0.001;
    params[ttC8]  = 0.001;
    params[ttC10] = 0.001;
    top[0]->setParams(params);

    testPot(Potential::TANG_TOENNIES, top, &x);
}

TEST_F (ForceComputerImplementationTest, BORN_MAYER)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[expA] = 10000;
    params[expB] = 8;
    top[0]->setParams(params);

    testPot(Potential::BORN_MAYER, top, &x);
}

TEST_F (ForceComputerImplementationTest, SLATER_ISA)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(2);
    params[expA] = 10000;
    params[expB] = 8;
    top[0]->setParams(params);

    testPot(Potential::SLATER_ISA, top, &x);
}

TEST_F (ForceComputerImplementationTest, MACDANIEL_SCHMIDT)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(3);
    params[dexpA1] = 10000;
    params[dexpA1] = -3000;
    params[dexpB]  = 8;
    top[0]->setParams(params);

    testPot(Potential::MACDANIEL_SCHMIDT, top, &x);
}

TEST_F (ForceComputerImplementationTest, POLARIZATION)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(3);
    params[polALPHA]   = 0.1;
    params[polRHYPER]  = 0.2;
    params[polFCHYPER] = 1e4;
    top[0]->setParams(params);

    testPot(Potential::POLARIZATION, top, &x);
}

TEST_F (ForceComputerImplementationTest, HARMONIC_BONDS)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(3);
    params[bondKB]     = 100000;
    params[bondLENGTH] = 0.4;
    params[bondENERGY] = 10;
    top[0]->setParams(params);

    testPot(Potential::HARMONIC_BONDS, top, &x);
}

TEST_F (ForceComputerImplementationTest, CUBIC_BONDS)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(4);
    params[cubicDE]     = 100;
    params[cubicLENGTH] = 0.4;
    params[cubicRMAX]   = 0.6;
    params[cubicKB]     = 60000;
    top[0]->setParams(params);

    testPot(Potential::CUBIC_BONDS, top, &x);
}

TEST_F (ForceComputerImplementationTest, HUA_BONDS)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(4);
    params[huaDE]     = 100;
    params[huaLENGTH] = 0.4;
    params[huaB]      = 20;
    params[huaC]      = 0.1;
    top[0]->setParams(params);

    testPot(Potential::HUA_BONDS, top, &x);
}

TEST_F (ForceComputerImplementationTest, MORSE_BONDS)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.5 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    top.push_back(AtomPair(0, 1));
    std::vector<double> params(4);
    params[morseBETA]   = 12;
    params[morseDE]     = 100;
    params[morseD0]     = 40;
    params[morseLENGTH] = 0.6;
    top[0]->setParams(params);

    testPot(Potential::MORSE_BONDS, top, &x);
}

TEST_F (ForceComputerImplementationTest, HARMONIC_ANGLES)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.2 },
        { 0, 0.1, 0.3 },
    };
    // Generate topology info
    TopologyEntryVector top{};
    Bond b1(0, 1, 1);
    Bond b2(1, 2, 1);
    top.push_back(Angle(b1, b2));
    std::vector<double> params(2);
    params[angleKT]    = 100;
    params[angleANGLE] = 100;
    top[0]->setParams(params);

    testPot(Potential::HARMONIC_ANGLES, top, &x);
}

TEST_F (ForceComputerImplementationTest, LINEAR_ANGLES)
{
    std::vector<gmx::RVec> x = {
        { 0, 0,    0    },
        { 0, 0.01, 0.15 },
        { 0, 0,    0.29  },
    };
    // Generate topology info
    TopologyEntryVector top{};
    Bond b1(0, 1, 1);
    Bond b2(1, 2, 1);
    top.push_back(Angle(b1, b2));
    std::vector<double> params(2);
    params[linangA]    = 0.5;
    params[linangKLIN] = 10000;
    top[0]->setParams(params);

    testPot(Potential::LINEAR_ANGLES, top, &x);
}

TEST_F (ForceComputerImplementationTest, UREY_BRADLEY_ANGLES)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.2 },
        { 0, 0.1, 0.3 },
    };
    // Generate topology info
    TopologyEntryVector top{};
    Bond b1(0, 1, 1);
    Bond b2(1, 2, 1);
    top.push_back(Angle(b1, b2));
    std::vector<double> params(4);
    params[ubKT]    = 100;
    params[ubANGLE] = 100;
    params[ubR13]   = 0.33;
    params[ubKUB]   = 20;
    top[0]->setParams(params);

    testPot(Potential::UREY_BRADLEY_ANGLES, top, &x);
}

TEST_F (ForceComputerImplementationTest, HARMONIC_DIHEDRALS)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.2 },
        { 0, 0.1, 0.3 },
        { 0.1, 0.15, 0.36 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    Bond b1(0, 1, 1);
    Bond b2(0, 2, 1);
    Bond b3(0, 3, 1);
    top.push_back(Improper(b1, b2, b3));
    std::vector<double> params(1);
    params[idihKPHI] = 10;
    top[0]->setParams(params);

    testPot(Potential::HARMONIC_DIHEDRALS, top, &x);
}

TEST_F (ForceComputerImplementationTest, PROPER_DIHEDRALS)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.2 },
        { 0, 0.1, 0.3 },
        { 0.1, 0.15, 0.36 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    Bond b1(0, 1, 1);
    Bond b2(1, 2, 1);
    Bond b3(2, 3, 1);
    top.push_back(Proper(b1, b2, b3));
    std::vector<double> params(3);
    params[pdihANGLE] = 30;
    params[pdihKP]    = 10;
    params[pdihMULT]  = 4;
    top[0]->setParams(params);

    testPot(Potential::PROPER_DIHEDRALS, top, &x);
}

TEST_F (ForceComputerImplementationTest, FOURIER_DIHEDRALS)
{
    std::vector<gmx::RVec> x = {
        { 0, 0, 0   },
        { 0, 0, 0.2 },
        { 0, 0.1, 0.3 },
        { 0.1, 0.15, 0.36 }
    };
    // Generate topology info
    TopologyEntryVector top{};
    Bond b1(0, 1, 1);
    Bond b2(1, 2, 1);
    Bond b3(2, 3, 1);
    top.push_back(Proper(b1, b2, b3));
    std::vector<double> params = { 1, 2, -3, 4, -5, 6 };
    top[0]->setParams(params);

    testPot(Potential::FOURIER_DIHEDRALS, top, &x);
}

}  // namespace

}  // namespace alexandria
