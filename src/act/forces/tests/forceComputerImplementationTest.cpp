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
 * \ingroup group_forces_tests
 */
#include "actpre.h"

#include "../forcecomputerimpl.h"

#include <cmath>

#include <functional>
#include <memory>
#include <string>
#include <tuple>
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

/*! \brief Atomic coordinates and charge configuration for a force computer test case */
struct ForceComputerCoordParams
{
    //! Name suffix used in test identification (empty = use potential name only)
    std::string                             name;
    //! Factory function that builds the atomic coordinates for this test
    std::function<std::vector<gmx::RVec>()> coordsBuilder;
    //! Charge for atoms_[1] (default 1.0 matches base fixture)
    double                                  charge1 = 1.0;
};

/*! \brief Potential type and topology parameters for a force computer test case */
struct ForceComputerPotParams
{
    //! Base name used in test identification (maps to refdata file)
    std::string                          name;
    //! The potential to test
    Potential                            potential;
    //! Factory function that builds the topology/parameters for this test
    std::function<TopologyEntryVector()> topBuilder;
};

//! Combined test parameter: (potential params, coordinate config)
using FCTestParam = std::tuple<ForceComputerPotParams, ForceComputerCoordParams>;

class ForceComputerImplementationTest : public ::testing::TestWithParam<FCTestParam>
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

        auto epot = bfc(nullptr, top, atoms_, coordinates, &forces, &energies);
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
            ener[k] = bfc(nullptr, top, atoms_, &newx, &forces, &energies);
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

TEST_P(ForceComputerImplementationTest, RunTest)
{
    const auto &pot   = std::get<0>(GetParam());
    const auto &coord = std::get<1>(GetParam());
    atoms_[1].setCharge(coord.charge1);
    auto x   = coord.coordsBuilder();
    auto top = pot.topBuilder();
    testPot(pot.potential, top, &x);
}

// ============================================================
// Coordinate parameter sets
// ============================================================

//! Standard 2-atom pair coord: atoms at distance 0.5 nm, charge1 = +1
static const ForceComputerCoordParams c_stdPairCoord {
    "", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.5 } }; }, 1.0
};

//! Coulomb Gaussian non-standard coord variants (NEG charge, ZERO distance)
static const std::vector<ForceComputerCoordParams> c_coulGaussCoords = {
    { "NEG",  []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.5 } }; }, -1.0 },
    { "ZERO", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.0 } }; },  1.0 },
};

//! Coulomb Slater (zeta=10) non-standard coord variants
static const std::vector<ForceComputerCoordParams> c_coulSlaterCoords = {
    { "NEG",       []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.5 } }; }, -1.0 },
    { "NEG2",      []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.5 } }; }, -2.0 },
    { "CLOSE",     []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.2 } }; },  1.0 },
    { "CLOSE_NEG", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.2 } }; }, -1.0 },
    { "ZERO",      []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0   } }; },  1.0 },
    { "ZERO_NEG",  []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0   } }; }, -1.0 },
};

//! Coulomb Slater (zeta=0) non-standard coord variants
static const std::vector<ForceComputerCoordParams> c_coulSlaterZ0Coords = {
    { "NEG", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.5 } }; }, -1.0 },
};

//! Standard 3-atom angle coord (used by HARMONIC_ANGLES and UREY_BRADLEY_ANGLES)
static const ForceComputerCoordParams c_angleCoord {
    "", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.2 }, { 0, 0.1, 0.3 } }; }, 1.0
};

//! Near-linear 3-atom angle coord (used by LINEAR_ANGLES)
static const ForceComputerCoordParams c_linAngleCoord {
    "", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0.01, 0.15 }, { 0, 0, 0.29 } }; }, 1.0
};

//! Standard 4-atom dihedral coord (used by all dihedral potentials)
static const ForceComputerCoordParams c_dihedralCoord {
    "", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.2 }, { 0, 0.1, 0.3 }, { 0.1, 0.15, 0.36 } }; }, 1.0
};

//! Position restraint coord variants (z = 0.5, 2, -2 nm)
static const std::vector<ForceComputerCoordParams> c_posResCoords = {
    { "0p5", []() { return std::vector<gmx::RVec>{ { 0, 0,  0.5 } }; }, 1.0 },
    { "2",   []() { return std::vector<gmx::RVec>{ { 0, 0,  2   } }; }, 1.0 },
    { "n2",  []() { return std::vector<gmx::RVec>{ { 0, 0, -2   } }; }, 1.0 },
};

// ============================================================
// Potential parameter sets
// ============================================================

//! Standard pair potentials: all use c_stdPairCoord (2 atoms, 0.5 nm, charge=+1)
static const std::vector<ForceComputerPotParams> c_stdPotentials = {
    // LJ8_6
    { "LJ8_6", Potential::LJ8_6, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[lj8_6SIGMA]   = 0.5;
        params[lj8_6EPSILON] = 0.25;
        top[0]->setParams(params);
        return top;
    } },
    // LJ12_6
    { "LJ12_6", Potential::LJ12_6, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[lj12_6SIGMA]   = 0.5;
        params[lj12_6EPSILON] = 0.25;
        top[0]->setParams(params);
        return top;
    } },
    // LJ12_6_4
    { "LJ12_6_4", Potential::LJ12_6_4, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(3);
        params[lj12_6_4SIGMA]   = 0.5;
        params[lj12_6_4EPSILON] = 0.25;
        params[lj12_6_4GAMMA]   = 0.5;
        top[0]->setParams(params);
        return top;
    } },
    // LJ14_7
    { "LJ14_7", Potential::LJ14_7, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(4);
        params[lj14_7SIGMA]   = 0.5;
        params[lj14_7EPSILON] = 1;
        params[lj14_7GAMMA]   = 0.5;
        params[lj14_7DELTA]   = 0.1;
        top[0]->setParams(params);
        return top;
    } },
    // BUCKINGHAM
    { "BUCKINGHAM", Potential::BUCKINGHAM, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(3);
        params[bhA]  = 5000;
        params[bhB]  = 20;
        params[bhC6] = 0.001;
        top[0]->setParams(params);
        return top;
    } },
    // WANG_BUCKINGHAM
    { "WANG_BUCKINGHAM", Potential::WANG_BUCKINGHAM, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(3);
        params[wbhSIGMA]   = 0.5;
        params[wbhEPSILON] = 0.25;
        params[wbhGAMMA]   = 10.0;
        top[0]->setParams(params);
        return top;
    } },
    // GENERALIZED_BUCKINGHAM
    { "GENERALIZED_BUCKINGHAM", Potential::GENERALIZED_BUCKINGHAM, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(4);
        params[gbhRMIN]    = 0.5;
        params[gbhEPSILON] = 0.25;
        params[gbhGAMMA]   = 0.5;
        params[gbhDELTA]   = 0.1;
        top[0]->setParams(params);
        return top;
    } },
    // COULOMB_POINT
    { "COULOMB_POINT", Potential::COULOMB_POINT, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[coulZETA]  = 10;
        params[coulZETA2] = 6;
        top[0]->setParams(params);
        return top;
    } },
    // COULOMB_GAUSSIAN (default: distance 0.5, charge +1)
    { "COULOMB_GAUSSIAN", Potential::COULOMB_GAUSSIAN, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[coulZETA]  = 10;
        params[coulZETA2] = 6;
        top[0]->setParams(params);
        return top;
    } },
    // COULOMB_SLATER zeta=10 (default: distance 0.5, charge +1)
    { "COULOMB_SLATER", Potential::COULOMB_SLATER, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[coulZETA]  = 10;
        params[coulZETA2] = 6;
        top[0]->setParams(params);
        return top;
    } },
    // COULOMB_SLATER zeta=0 (default: distance 0.5, charge +1)
    { "COULOMB_SLATER_ZETA0", Potential::COULOMB_SLATER, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[coulZETA]  = 0;
        params[coulZETA2] = 6;
        top[0]->setParams(params);
        return top;
    } },
    // TT2bExch
    { "TT2bExch", Potential::TT2b, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(6, 0.0);
        params[tt2bA]     = 1000;
        params[tt2bBexch] = 10;
        params[tt2bBdisp] = 10;
        top[0]->setParams(params);
        return top;
    } },
    // TT2bDispC6
    { "TT2bDispC6", Potential::TT2b, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(6, 0.0);
        params[tt2bA]     = 0;
        params[tt2bBexch] = 10;
        params[tt2bBdisp] = 10;
        params[tt2bC6]    = 0.001;
        top[0]->setParams(params);
        return top;
    } },
    // TT2bDispC8
    { "TT2bDispC8", Potential::TT2b, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(6, 0.0);
        params[tt2bA]     = 0;
        params[tt2bBexch] = 10;
        params[tt2bBdisp] = 10;
        params[tt2bC8]    = 0.0001;
        top[0]->setParams(params);
        return top;
    } },
    // TT2bDispC10
    { "TT2bDispC10", Potential::TT2b, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(6, 0.0);
        params[tt2bA]     = 0;
        params[tt2bBexch] = 10;
        params[tt2bBdisp] = 10;
        params[tt2bC10]   = 0.001;
        top[0]->setParams(params);
        return top;
    } },
    // TT2bAll
    { "TT2bAll", Potential::TT2b, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(6, 0.0);
        params[tt2bA]     = 1000;
        params[tt2bBexch] = 10;
        params[tt2bBdisp] = 20;
        params[tt2bC6]    = 0.001;
        params[tt2bC8]    = 0.001;
        params[tt2bC10]   = 0.001;
        top[0]->setParams(params);
        return top;
    } },
    // SLATER_ISA_TT
    { "SLATER_ISA_TT", Potential::SLATER_ISA_TT, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(6, 0.0);
        params[tt2bA]     = 1000;
        params[tt2bBexch] = 10;
        params[tt2bBdisp] = 20;
        params[tt2bC6]    = 0.001;
        params[tt2bC8]    = 0.001;
        params[tt2bC10]   = 0.001;
        top[0]->setParams(params);
        return top;
    } },
    // TANG_TOENNIES
    { "TANG_TOENNIES", Potential::TANG_TOENNIES, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(5, 0.0);
        params[ttA]   = 1000;
        params[ttB]   = 10;
        params[ttC6]  = 0.001;
        params[ttC8]  = 0.001;
        params[ttC10] = 0.001;
        top[0]->setParams(params);
        return top;
    } },
    // BORN_MAYER
    { "BORN_MAYER", Potential::BORN_MAYER, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[expA] = 10000;
        params[expB] = 8;
        top[0]->setParams(params);
        return top;
    } },
    // SLATER_ISA
    { "SLATER_ISA", Potential::SLATER_ISA, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[expA] = 10000;
        params[expB] = 8;
        top[0]->setParams(params);
        return top;
    } },
    // MACDANIEL_SCHMIDT
    { "MACDANIEL_SCHMIDT", Potential::MACDANIEL_SCHMIDT, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(3);
        params[dexpA1] = 10000;
        params[dexpA1] = -3000; // intentional overwrite, preserved from original test to keep refdata unchanged
        params[dexpB]  = 8;
        top[0]->setParams(params);
        return top;
    } },
    // POLARIZATION
    { "POLARIZATION", Potential::POLARIZATION, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(3);
        params[polALPHA]   = 0.1;
        params[polRHYPER]  = 0.2;
        params[polFCHYPER] = 1e4;
        top[0]->setParams(params);
        return top;
    } },
    // HARMONIC_BONDS
    { "HARMONIC_BONDS", Potential::HARMONIC_BONDS, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(3);
        params[bondKB]     = 100000;
        params[bondLENGTH] = 0.4;
        params[bondENERGY] = 10;
        top[0]->setParams(params);
        return top;
    } },
    // CUBIC_BONDS
    { "CUBIC_BONDS", Potential::CUBIC_BONDS, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(4);
        params[cubicDE]     = 100;
        params[cubicLENGTH] = 0.4;
        params[cubicRMAX]   = 0.6;
        params[cubicKB]     = 60000;
        top[0]->setParams(params);
        return top;
    } },
    // HUA_BONDS
    { "HUA_BONDS", Potential::HUA_BONDS, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(4);
        params[huaDE]     = 100;
        params[huaLENGTH] = 0.4;
        params[huaB]      = 20;
        params[huaC]      = 0.1;
        top[0]->setParams(params);
        return top;
    } },
    // MORSE_BONDS
    { "MORSE_BONDS", Potential::MORSE_BONDS, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(4);
        params[morseBETA]   = 12;
        params[morseDE]     = 100;
        params[morseD0]     = 40;
        params[morseLENGTH] = 0.6;
        top[0]->setParams(params);
        return top;
    } },
};

//! Coulomb Gaussian potential (shared topology for default + non-standard coord variants)
static const ForceComputerPotParams c_coulGaussPot {
    "COULOMB_GAUSSIAN", Potential::COULOMB_GAUSSIAN, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[coulZETA]  = 10;
        params[coulZETA2] = 6;
        top[0]->setParams(params);
        return top;
    }
};

//! Coulomb Slater potential with zeta=10 (shared topology for non-standard coord variants)
static const ForceComputerPotParams c_coulSlaterPot {
    "COULOMB_SLATER", Potential::COULOMB_SLATER, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[coulZETA]  = 10;
        params[coulZETA2] = 6;
        top[0]->setParams(params);
        return top;
    }
};

//! Coulomb Slater potential with zeta=0 (shared topology for non-standard coord variants)
static const ForceComputerPotParams c_coulSlaterZ0Pot {
    "COULOMB_SLATER_ZETA0", Potential::COULOMB_SLATER, []() {
        TopologyEntryVector top;
        top.push_back(AtomPair(0, 1));
        std::vector<double> params(2);
        params[coulZETA]  = 0;
        params[coulZETA2] = 6;
        top[0]->setParams(params);
        return top;
    }
};

//! Angle potentials sharing the same 3-atom coordinate set
static const std::vector<ForceComputerPotParams> c_anglePotentials = {
    // HARMONIC_ANGLES
    { "HARMONIC_ANGLES", Potential::HARMONIC_ANGLES, []() {
        TopologyEntryVector top;
        Bond b1(0, 1, 1);
        Bond b2(1, 2, 1);
        top.push_back(Angle(b1, b2));
        std::vector<double> params(2);
        params[angleKT]    = 100;
        params[angleANGLE] = 100;
        top[0]->setParams(params);
        return top;
    } },
    // UREY_BRADLEY_ANGLES
    { "UREY_BRADLEY_ANGLES", Potential::UREY_BRADLEY_ANGLES, []() {
        TopologyEntryVector top;
        Bond b1(0, 1, 1);
        Bond b2(1, 2, 1);
        top.push_back(Angle(b1, b2));
        std::vector<double> params(4);
        params[ubKT]    = 100;
        params[ubANGLE] = 100;
        params[ubR13]   = 0.33;
        params[ubKUB]   = 20;
        top[0]->setParams(params);
        return top;
    } },
};

//! Linear angle potential (uses a near-linear 3-atom coordinate set)
static const ForceComputerPotParams c_linAnglePot {
    "LINEAR_ANGLES", Potential::LINEAR_ANGLES, []() {
        TopologyEntryVector top;
        Bond b1(0, 1, 1);
        Bond b2(1, 2, 1);
        top.push_back(Angle(b1, b2));
        std::vector<double> params(2);
        params[linangA]    = 0.5;
        params[linangKLIN] = 10000;
        top[0]->setParams(params);
        return top;
    }
};

//! Dihedral potentials sharing the same 4-atom coordinate set
static const std::vector<ForceComputerPotParams> c_dihedralPotentials = {
    // HARMONIC_DIHEDRALS
    { "HARMONIC_DIHEDRALS", Potential::HARMONIC_DIHEDRALS, []() {
        TopologyEntryVector top;
        Bond b1(0, 1, 1);
        Bond b2(0, 2, 1);
        Bond b3(0, 3, 1);
        top.push_back(Improper(b1, b2, b3));
        std::vector<double> params(1);
        params[idihKPHI] = 10;
        top[0]->setParams(params);
        return top;
    } },
    // PROPER_DIHEDRALS
    { "PROPER_DIHEDRALS", Potential::PROPER_DIHEDRALS, []() {
        TopologyEntryVector top;
        Bond b1(0, 1, 1);
        Bond b2(1, 2, 1);
        Bond b3(2, 3, 1);
        top.push_back(Proper(b1, b2, b3));
        std::vector<double> params(3);
        params[pdihANGLE] = 30;
        params[pdihKP]    = 10;
        params[pdihMULT]  = 4;
        top[0]->setParams(params);
        return top;
    } },
    // FOURIER_DIHEDRALS
    { "FOURIER_DIHEDRALS", Potential::FOURIER_DIHEDRALS, []() {
        TopologyEntryVector top;
        Bond b1(0, 1, 1);
        Bond b2(1, 2, 1);
        Bond b3(2, 3, 1);
        top.push_back(Proper(b1, b2, b3));
        std::vector<double> params = { 1, 2, -3, 4, -5, 6 };
        top[0]->setParams(params);
        return top;
    } },
};

//! Position restraint potential (flat-bottomed, K=100, R0=1)
static const ForceComputerPotParams c_posResPot {
    "FBPOSRES", Potential::POSITION_RESTRAINT, []() {
        TopologyEntryVector top;
        top.push_back(SingleAtom(0));
        std::vector<double> params(2);
        params[fbprK]  = 100;
        params[fbprR0] = 1;
        top[0]->setParams(params);
        return top;
    }
};

// ============================================================
// Name printer: combines pot name and coord name suffix
// ============================================================
static std::string fcTestName(const ::testing::TestParamInfo<FCTestParam> &info)
{
    const auto &pot   = std::get<0>(info.param);
    const auto &coord = std::get<1>(info.param);
    return coord.name.empty() ? pot.name : pot.name + "_" + coord.name;
}

// ============================================================
// Instantiations
// ============================================================

//! Standard 2-atom pair potentials (2 atoms, 0.5 nm separation, charge +1)
INSTANTIATE_TEST_CASE_P(All, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::ValuesIn(c_stdPotentials),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

//! Coulomb Gaussian with non-standard charge or zero-distance coord
INSTANTIATE_TEST_CASE_P(CoulGauss, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulGaussPot),
                                            ::testing::ValuesIn(c_coulGaussCoords)),
                         fcTestName);

//! Coulomb Slater (zeta=10) with non-standard charge/distance coord variants
INSTANTIATE_TEST_CASE_P(CoulSlater, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulSlaterPot),
                                            ::testing::ValuesIn(c_coulSlaterCoords)),
                         fcTestName);

//! Coulomb Slater (zeta=0) with negative charge
INSTANTIATE_TEST_CASE_P(CoulSlaterZ0, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulSlaterZ0Pot),
                                            ::testing::ValuesIn(c_coulSlaterZ0Coords)),
                         fcTestName);

//! Angle potentials (HARMONIC and UREY-BRADLEY share the same 3-atom coord set)
INSTANTIATE_TEST_CASE_P(Angle, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::ValuesIn(c_anglePotentials),
                                            ::testing::Values(c_angleCoord)),
                         fcTestName);

//! Linear angle potential with near-linear 3-atom coord set
INSTANTIATE_TEST_CASE_P(LinearAngle, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_linAnglePot),
                                            ::testing::Values(c_linAngleCoord)),
                         fcTestName);

//! Dihedral potentials (3 types sharing the same 4-atom coord set)
INSTANTIATE_TEST_CASE_P(Dihedral, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::ValuesIn(c_dihedralPotentials),
                                            ::testing::Values(c_dihedralCoord)),
                         fcTestName);

//! Flat-bottomed position restraint with 3 different single-atom positions
INSTANTIATE_TEST_CASE_P(PosRes, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_posResPot),
                                            ::testing::ValuesIn(c_posResCoords)),
                         fcTestName);

}  // namespace
}  // namespace alexandria
