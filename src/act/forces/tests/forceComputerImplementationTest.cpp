/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2025,2026
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

/*! \brief Force field parameter values for a force computer test case.
 *  Kept separate from ForceComputerPotParams so that a single potential can be
 *  tested with multiple parameter sets via ::testing::Combine. */
struct ForceComputerParamParams
{
    //! Name suffix appended to the potential name in the test identifier (empty for default)
    std::string         nameSuffix;
    //! Force field parameter values passed to TopologyEntry::setParams()
    std::vector<double> params;
};

/*! \brief Potential type and topology structure for a force computer test case.
 *  Does NOT set force field parameters; those come from ForceComputerParamParams. */
struct ForceComputerPotParams
{
    //! Base name used in test identification (maps to refdata file)
    std::string                          name;
    //! The potential to test
    Potential                            potential;
    //! Factory function that builds the topology entries (without setting params)
    std::function<TopologyEntryVector()> topBuilder;
};

//! Combined test parameter: (potential topology, force field params, coordinate config)
using FCTestParam = std::tuple<ForceComputerPotParams, ForceComputerParamParams, ForceComputerCoordParams>;

/*! \brief Returns a refdata filename stem without the INSTANTIATE prefix.
 *
 * For a parameterised test `INST/FixtureName` + `Method/Param`, Google Test
 * would normally produce `INST_FixtureName_Method_Param`.  This helper strips
 * the instantiation prefix, yielding `FixtureName_Method_Param`, so that the
 * XML files are named `ForceComputerImplementationTest_All_<pot>.xml`.
 */
static std::string fcStrippedRefDataFilename()
{
    const ::testing::TestInfo *info =
        ::testing::UnitTest::GetInstance()->current_test_info();
    std::string caseName = info->test_case_name();
    auto        slash    = caseName.find('/');
    if (slash != std::string::npos)
    {
        caseName = caseName.substr(slash + 1);
    }
    std::string stem = caseName + "_" + info->name();
    std::replace(stem.begin(), stem.end(), '/', '_');
    return stem;
}

class ForceComputerImplementationTest : public ::testing::TestWithParam<FCTestParam>
{
protected:
    gmx::test::TestReferenceData    refData_;
    gmx::test::TestReferenceChecker checker_;
    std::vector<ActAtom>            atoms_;

    ForceComputerImplementationTest( ) :
        refData_(fcStrippedRefDataFilename()),
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
        // Add check that energy is the same if we change the sign of coordinates 
        std::vector<gmx::RVec> forces2(coordinates->size(), fzero);
        std::vector<gmx::RVec> coords2(coordinates->size());
        for (size_t i = 0; i < coordinates->size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                coords2[i][m] = -(*coordinates)[i][m];
            }
        }
        // Energy map
        std::map<InteractionType, double> energies2;
        (void) bfc(nullptr, top, atoms_, &coords2, &forces2, &energies2);
        double ediff = 1e-8;
        for(auto &ee2 : energies2)
        { 
            EXPECT_TRUE(std::abs(ee2.second-energies[ee2.first]) < ediff);
        }
        double fdiff = 1e-4;
        for(size_t i = 0; i < forces2.size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                EXPECT_TRUE(std::abs(forces2[i][m] + forces[i][m]) < fdiff);
            }
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

TEST_P(ForceComputerImplementationTest, All)
{
    const auto &pot   = std::get<0>(GetParam());
    const auto &par   = std::get<1>(GetParam());
    const auto &coord = std::get<2>(GetParam());
    atoms_[1].setCharge(coord.charge1);
    auto x   = coord.coordsBuilder();
    auto top = pot.topBuilder();
    if (!top.empty())
    {
        top[0]->setParams(par.params);
    }
    testPot(pot.potential, top, &x);
}

// ============================================================
// Coordinate parameter sets
// ============================================================

//! Standard 2-atom pair coord: atoms at distance 0.5 nm, charge1 = +1
static const ForceComputerCoordParams c_stdPairCoord {
    "", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.5 } }; }, 1.0
};

//! Standard 3-atom angle coord
static const ForceComputerCoordParams c_angleCoord {
    "", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.2 }, { 0, 0.1, 0.3 } }; }, 1.0
};

//! Near-linear 3-atom angle coord (used by LINEAR_ANGLES)
static const ForceComputerCoordParams c_linAngleCoord {
    "", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0.01, 0.15 }, { 0, 0, 0.29 } }; }, 1.0
};

//! Standard 4-atom dihedral coord
static const ForceComputerCoordParams c_dihedralCoord {
    "", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.2 }, { 0, 0.1, 0.3 }, { 0.1, 0.15, 0.36 } }; }, 1.0
};

//! Coulomb non-standard coord variants: negative charge on atom 1, distance 0.5 nm
static const ForceComputerCoordParams c_coulNegCoord {
    "NEG", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.5 } }; }, -1.0
};

//! Coulomb SLATER non-standard charge: charge = -2
static const ForceComputerCoordParams c_coulNeg2Coord {
    "NEG2", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.5 } }; }, -2.0
};

//! Coulomb non-standard coord: zero interatomic distance
static const ForceComputerCoordParams c_coulZeroCoord {
    "ZERO", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.0 } }; }, 1.0
};

//! Coulomb SLATER close-distance coord: 0.2 nm
static const ForceComputerCoordParams c_coulCloseCoord {
    "CLOSE", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.2 } }; }, 1.0
};

//! Coulomb SLATER close-distance + negative charge
static const ForceComputerCoordParams c_coulCloseNegCoord {
    "CLOSE_NEG", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0.2 } }; }, -1.0
};

//! Coulomb SLATER zero-distance + negative charge
static const ForceComputerCoordParams c_coulZeroNegCoord {
    "ZERO_NEG", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0, 0, 0 } }; }, -1.0
};

//! Position restraint coord variants (z = 0.5, 2, -2 nm)
static const std::vector<ForceComputerCoordParams> c_posResCoords = {
    { "0p5", []() { return std::vector<gmx::RVec>{ { 0, 0,  0.5 } }; }, 1.0 },
    { "2",   []() { return std::vector<gmx::RVec>{ { 0, 0,  2   } }; }, 1.0 },
    { "n2",  []() { return std::vector<gmx::RVec>{ { 0, 0, -2   } }; }, 1.0 },
};

// ============================================================
// Topology builders (structure only; params applied separately)
// ============================================================

//! Atom-pair topology for standard 2-body potentials
static TopologyEntryVector makePairTop()
{
    TopologyEntryVector top;
    top.push_back(AtomPair(0, 1));
    return top;
}

//! Angle topology (bond 0-1, bond 1-2)
static TopologyEntryVector makeAngleTop()
{
    Bond b1(0, 1, 1);
    Bond b2(1, 2, 1);
    TopologyEntryVector top;
    top.push_back(Angle(b1, b2));
    return top;
}

//! Improper dihedral topology (bonds 0-1, 0-2, 0-3)
static TopologyEntryVector makeImproperTop()
{
    Bond b1(0, 1, 1);
    Bond b2(0, 2, 1);
    Bond b3(0, 3, 1);
    TopologyEntryVector top;
    top.push_back(Improper(b1, b2, b3));
    return top;
}

//! Proper dihedral topology (bonds 0-1, 1-2, 2-3)
static TopologyEntryVector makeProperTop()
{
    Bond b1(0, 1, 1);
    Bond b2(1, 2, 1);
    Bond b3(2, 3, 1);
    TopologyEntryVector top;
    top.push_back(Proper(b1, b2, b3));
    return top;
}

//! Single-atom topology (for position restraints)
static TopologyEntryVector makeSingleAtomTop()
{
    TopologyEntryVector top;
    top.push_back(SingleAtom(0));
    return top;
}

// ============================================================
// Potential descriptors (topology structure only, no params)
// ============================================================

static const ForceComputerPotParams c_lj8_6Pot        { "LJ8_6",                Potential::LJ8_6,                makePairTop };
static const ForceComputerPotParams c_lj12_6Pot       { "LJ12_6",               Potential::LJ12_6,               makePairTop };
static const ForceComputerPotParams c_lj12_6_4Pot     { "LJ12_6_4",             Potential::LJ12_6_4,             makePairTop };
static const ForceComputerPotParams c_lj14_7Pot       { "LJ14_7",               Potential::LJ14_7,               makePairTop };
static const ForceComputerPotParams c_buckinghamPot   { "BUCKINGHAM",           Potential::BUCKINGHAM,           makePairTop };
static const ForceComputerPotParams c_wangBuckPot     { "WANG_BUCKINGHAM",      Potential::WANG_BUCKINGHAM,      makePairTop };
static const ForceComputerPotParams c_genBuckPot      { "GENERALIZED_BUCKINGHAM", Potential::GENERALIZED_BUCKINGHAM, makePairTop };
static const ForceComputerPotParams c_coulPointPot    { "COULOMB_POINT",        Potential::COULOMB_POINT,        makePairTop };
static const ForceComputerPotParams c_coulGaussPot    { "COULOMB_GAUSSIAN",     Potential::COULOMB_GAUSSIAN,     makePairTop };
static const ForceComputerPotParams c_coulSlaterPot   { "COULOMB_SLATER",       Potential::COULOMB_SLATER,       makePairTop };
//! COULOMB_SLATER with zeta=0 — distinct name to keep separate refdata from the zeta=10 variant
static const ForceComputerPotParams c_coulSlaterZ0Pot { "COULOMB_SLATER_ZETA0", Potential::COULOMB_SLATER,       makePairTop };
static const ForceComputerPotParams c_tt2bPot         { "TT2b",                 Potential::TT2b,                 makePairTop };
static const ForceComputerPotParams c_slaterIsaTtPot  { "SLATER_ISA_TT",        Potential::SLATER_ISA_TT,        makePairTop };
static const ForceComputerPotParams c_tangToennPot    { "TANG_TOENNIES",        Potential::TANG_TOENNIES,        makePairTop };
static const ForceComputerPotParams c_bornMayerPot    { "BORN_MAYER",           Potential::BORN_MAYER,           makePairTop };
static const ForceComputerPotParams c_slaterIsaPot    { "SLATER_ISA",           Potential::SLATER_ISA,           makePairTop };
static const ForceComputerPotParams c_macdanielPot    { "MACDANIEL_SCHMIDT",    Potential::MACDANIEL_SCHMIDT,    makePairTop };
static const ForceComputerPotParams c_polPot          { "POLARIZATION",         Potential::POLARIZATION,         makePairTop };
static const ForceComputerPotParams c_harmBondsPot    { "HARMONIC_BONDS",       Potential::HARMONIC_BONDS,       makePairTop };
static const ForceComputerPotParams c_cubicBondsPot   { "CUBIC_BONDS",          Potential::CUBIC_BONDS,          makePairTop };
static const ForceComputerPotParams c_huaBondsPot     { "HUA_BONDS",            Potential::HUA_BONDS,            makePairTop };
static const ForceComputerPotParams c_morseBondsPot   { "MORSE_BONDS",          Potential::MORSE_BONDS,          makePairTop };
static const ForceComputerPotParams c_harmAnglePot    { "HARMONIC_ANGLES",      Potential::HARMONIC_ANGLES,      makeAngleTop };
static const ForceComputerPotParams c_ubAnglePot      { "UREY_BRADLEY_ANGLES",  Potential::UREY_BRADLEY_ANGLES,  makeAngleTop };
static const ForceComputerPotParams c_linAnglePot     { "LINEAR_ANGLES",        Potential::LINEAR_ANGLES,        makeAngleTop };
static const ForceComputerPotParams c_harmDihPot      { "HARMONIC_DIHEDRALS",   Potential::HARMONIC_DIHEDRALS,   makeImproperTop };
static const ForceComputerPotParams c_properDihPot    { "PROPER_DIHEDRALS",     Potential::PROPER_DIHEDRALS,     makeProperTop };
static const ForceComputerPotParams c_fourierDihPot   { "FOURIER_DIHEDRALS",    Potential::FOURIER_DIHEDRALS,    makeProperTop };
static const ForceComputerPotParams c_posResPot       { "FBPOSRES",             Potential::POSITION_RESTRAINT,   makeSingleAtomTop };

// ============================================================
// Force field parameter sets
// ============================================================

static const ForceComputerParamParams c_lj8_6Params     { "", { 0.5, 0.25 } };   // sigma, epsilon
static const ForceComputerParamParams c_lj12_6Params    { "", { 0.5, 0.25 } };   // sigma, epsilon
static const ForceComputerParamParams c_lj12_6_4Params  { "", { 0.5, 0.25, 0.5 } };  // sigma, epsilon, gamma
static const ForceComputerParamParams c_lj14_7Params    { "", { 0.5, 1, 0.5, 0.1 } };  // sigma, epsilon, gamma, delta
static const ForceComputerParamParams c_buckinghamParams { "", { 5000, 20, 0.001 } };   // A, B, C6
static const ForceComputerParamParams c_wangBuckParams  { "", { 0.5, 0.25, 10.0 } };   // sigma, epsilon, gamma
static const ForceComputerParamParams c_genBuckParams   { "", { 0.5, 0.25, 0.5, 0.1 } };  // Rmin, epsilon, gamma, delta
//! Coulomb parameters: zeta=10, zeta2=6 (shared by COULOMB_POINT, COULOMB_GAUSSIAN, COULOMB_SLATER with zeta=10)
static const ForceComputerParamParams c_coulZeta10Params { "", { 10, 6 } };
//! Coulomb parameters: zeta=0, zeta2=6 (for COULOMB_SLATER with zeta=0)
static const ForceComputerParamParams c_coulZeta0Params  { "", { 0, 6 } };
//! TT2b: exchange repulsion only
static const ForceComputerParamParams c_tt2bExchParams   { "Exch",    { 1000, 10, 10, 0, 0, 0 } };
//! TT2b: dispersion C6 only
static const ForceComputerParamParams c_tt2bDispC6Params { "DispC6",  { 0, 10, 10, 0.001, 0, 0 } };
//! TT2b: dispersion C8 only
static const ForceComputerParamParams c_tt2bDispC8Params { "DispC8",  { 0, 10, 10, 0, 0.0001, 0 } };
//! TT2b: dispersion C10 only
static const ForceComputerParamParams c_tt2bDispC10Params{ "DispC10", { 0, 10, 10, 0, 0, 0.001 } };
//! TT2b: all terms
static const ForceComputerParamParams c_tt2bAllParams    { "All",     { 1000, 10, 20, 0.001, 0.001, 0.001 } };
// Note: TT2b param order is: tt2bA, tt2bBexch, tt2bBdisp, tt2bC6, tt2bC8, tt2bC10
static const std::vector<ForceComputerParamParams> c_tt2bParamSets = {
    c_tt2bExchParams, c_tt2bDispC6Params, c_tt2bDispC8Params, c_tt2bDispC10Params, c_tt2bAllParams
};
//! SLATER_ISA_TT uses the same TT2b parameter structure with all terms active
static const ForceComputerParamParams c_slaterIsaTtParams { "", { 1000, 10, 20, 0.001, 0.001, 0.001 } };
//! TANG_TOENNIES: A, B, C6, C8, C10
static const ForceComputerParamParams c_tangToennParams { "", { 1000, 10, 0.001, 0.001, 0.001 } };
//! Born-Mayer / SLATER_ISA exponential repulsion: A, B
static const ForceComputerParamParams c_expParams       { "", { 10000, 8 } };
//! MACDANIEL_SCHMIDT: dexpA1, dexpA2, dexpB
static const ForceComputerParamParams c_macdanielParams { "",
    []() {
        std::vector<double> p(3);
        p[dexpA1] = 10000;
        p[dexpA2] = 13000;
        p[dexpB]  = 8;
        return p;
    }()
};
static const ForceComputerParamParams c_polParams       { "", { 0.1, 0.2, 1e4 } };   // alpha, rHyper, fcHyper
static const ForceComputerParamParams c_harmBondsParams { "", { 100000, 0.4, 10 } };  // kB, length, energy
static const ForceComputerParamParams c_cubicBondsParams{ "", { 0.4, 0.6, 60000, 100 } };  // length, Rmax, kB, DE (cubicLENGTH, cubicRMAX, cubicKB, cubicDE)
static const ForceComputerParamParams c_huaBondsParams  { "", { 0.4, 100, 20, 0.1 } };    // length, DE, b, c (huaLENGTH, huaDE, huaB, huaC)
static const ForceComputerParamParams c_morseBondsParams{ "", { 12, 100, 40, 0.6 } };     // beta, DE, D0, length
static const ForceComputerParamParams c_harmAngleParams { "", { 100, 100 } };     // kT, angle
static const ForceComputerParamParams c_ubAngleParams   { "", { 100, 100, 0.33, 20 } };   // kT, angle, r13, kUB
static const ForceComputerParamParams c_linAngleParams  { "", { 0.5, 10000 } };   // a, kLin
static const ForceComputerParamParams c_harmDihParams   { "", { 10 } };           // kPhi
static const ForceComputerParamParams c_properDihParams { "", { 30, 10, 4 } };    // angle, kP, mult
static const ForceComputerParamParams c_fourierDihParams{ "", { 1, 2, -3, 4, -5, 6 } };
static const ForceComputerParamParams c_posResParams    { "", { 100, 1 } };        // k, R0

// ============================================================
// Name printer: pot.name + par.nameSuffix [+ "_" + coord.name]
// ============================================================
static std::string fcTestName(const ::testing::TestParamInfo<FCTestParam> &info)
{
    const auto &pot   = std::get<0>(info.param);
    const auto &par   = std::get<1>(info.param);
    const auto &coord = std::get<2>(info.param);
    std::string name  = pot.name + par.nameSuffix;
    if (!coord.name.empty())
    {
        name += "_" + coord.name;
    }
    return name;
}

// ============================================================
// Instantiations — one per potential family, using ::testing::Combine
// ============================================================

// --- LJ potentials ---
INSTANTIATE_TEST_CASE_P(LJ8_6, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_lj8_6Pot),
                                            ::testing::Values(c_lj8_6Params),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(LJ12_6, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_lj12_6Pot),
                                            ::testing::Values(c_lj12_6Params),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(LJ12_6_4, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_lj12_6_4Pot),
                                            ::testing::Values(c_lj12_6_4Params),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(LJ14_7, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_lj14_7Pot),
                                            ::testing::Values(c_lj14_7Params),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

// --- Buckingham potentials ---
INSTANTIATE_TEST_CASE_P(BUCKINGHAM, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_buckinghamPot),
                                            ::testing::Values(c_buckinghamParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(WANG_BUCKINGHAM, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_wangBuckPot),
                                            ::testing::Values(c_wangBuckParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(GENERALIZED_BUCKINGHAM, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_genBuckPot),
                                            ::testing::Values(c_genBuckParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

// --- Coulomb potentials (each also tested with multiple coord/charge variants) ---

//! COULOMB_POINT: only standard coord
INSTANTIATE_TEST_CASE_P(COULOMB_POINT, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulPointPot),
                                            ::testing::Values(c_coulZeta10Params),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

//! COULOMB_GAUSSIAN: standard coord + negative-charge + zero-distance variants
INSTANTIATE_TEST_CASE_P(CoulGauss, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulGaussPot),
                                            ::testing::Values(c_coulZeta10Params),
                                            ::testing::Values(c_stdPairCoord,
                                                              c_coulNegCoord,
                                                              c_coulZeroCoord)),
                         fcTestName);

//! COULOMB_SLATER zeta=10: standard coord + charge/distance variants
INSTANTIATE_TEST_CASE_P(CoulSlater, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulSlaterPot),
                                            ::testing::Values(c_coulZeta10Params),
                                            ::testing::Values(c_stdPairCoord,
                                                              c_coulNegCoord,
                                                              c_coulNeg2Coord,
                                                              c_coulCloseCoord,
                                                              c_coulCloseNegCoord,
                                                              c_coulZeroCoord,
                                                              c_coulZeroNegCoord)),
                         fcTestName);

//! COULOMB_SLATER zeta=0: standard coord + negative-charge variant
INSTANTIATE_TEST_CASE_P(CoulSlaterZ0, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulSlaterZ0Pot),
                                            ::testing::Values(c_coulZeta0Params),
                                            ::testing::Values(c_stdPairCoord,
                                                              c_coulNegCoord)),
                         fcTestName);

// --- TT2b family: ONE potential tested with FIVE different parameter sets ---
//! Demonstrates the Combine benefit: adding a new param set here tests TT2b under new conditions
INSTANTIATE_TEST_CASE_P(TT2b, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_tt2bPot),
                                            ::testing::ValuesIn(c_tt2bParamSets),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

// --- Other pair potentials ---
INSTANTIATE_TEST_CASE_P(SLATER_ISA_TT, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_slaterIsaTtPot),
                                            ::testing::Values(c_slaterIsaTtParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(TANG_TOENNIES, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_tangToennPot),
                                            ::testing::Values(c_tangToennParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(BORN_MAYER, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_bornMayerPot),
                                            ::testing::Values(c_expParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(SLATER_ISA, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_slaterIsaPot),
                                            ::testing::Values(c_expParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(MACDANIEL_SCHMIDT, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_macdanielPot),
                                            ::testing::Values(c_macdanielParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(POLARIZATION, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_polPot),
                                            ::testing::Values(c_polParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

// --- Bond potentials ---
INSTANTIATE_TEST_CASE_P(HARMONIC_BONDS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_harmBondsPot),
                                            ::testing::Values(c_harmBondsParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(CUBIC_BONDS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_cubicBondsPot),
                                            ::testing::Values(c_cubicBondsParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(HUA_BONDS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_huaBondsPot),
                                            ::testing::Values(c_huaBondsParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(MORSE_BONDS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_morseBondsPot),
                                            ::testing::Values(c_morseBondsParams),
                                            ::testing::Values(c_stdPairCoord)),
                         fcTestName);

// --- Angle potentials ---
INSTANTIATE_TEST_CASE_P(HARMONIC_ANGLES, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_harmAnglePot),
                                            ::testing::Values(c_harmAngleParams),
                                            ::testing::Values(c_angleCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(UREY_BRADLEY_ANGLES, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_ubAnglePot),
                                            ::testing::Values(c_ubAngleParams),
                                            ::testing::Values(c_angleCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(LINEAR_ANGLES, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_linAnglePot),
                                            ::testing::Values(c_linAngleParams),
                                            ::testing::Values(c_linAngleCoord)),
                         fcTestName);

// --- Dihedral potentials ---
INSTANTIATE_TEST_CASE_P(HARMONIC_DIHEDRALS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_harmDihPot),
                                            ::testing::Values(c_harmDihParams),
                                            ::testing::Values(c_dihedralCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(PROPER_DIHEDRALS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_properDihPot),
                                            ::testing::Values(c_properDihParams),
                                            ::testing::Values(c_dihedralCoord)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(FOURIER_DIHEDRALS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_fourierDihPot),
                                            ::testing::Values(c_fourierDihParams),
                                            ::testing::Values(c_dihedralCoord)),
                         fcTestName);

// --- Position restraint ---
INSTANTIATE_TEST_CASE_P(PosRes, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_posResPot),
                                            ::testing::Values(c_posResParams),
                                            ::testing::ValuesIn(c_posResCoords)),
                         fcTestName);

}  // namespace
}  // namespace alexandria
