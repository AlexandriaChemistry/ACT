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
 * Implements test of the force routines of each potential, see
 * computePotential in forces_simd.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_forces_tests
 */
#include "actpre.h"

#include "../forces_simd.h"
#include "../forces_explicit_simd.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/molprop/topologyentry.h"
#include "act/topology/coords_energy.h"
#include "act/topology/topology.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/stringutil.h"

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

/*! \brief Compute the energy and forces of one potential for one structure
 * \param[in]  p           The potential
 * \param[in]  top         The interactions
 * \param[in]  atoms       The particles
 * \param[in]  coordinates The coordinates
 * \param[out] forces      The forces
 * \param[out] energies    The energy terms, without the total
 * \param[in]  epsilonr    Relative dielectric constant
 * \return The total energy
 */
static double computeOne(Potential                          p,
                         const TopologyEntryVector         &top,
                         const std::vector<ActAtom>        &atoms,
                         const std::vector<gmx::RVec>      *coordinates,
                         std::vector<gmx::RVec>            *forces,
                         std::map<InteractionType, double> *energies,
                         double                             epsilonr)
{
    ConformationData cdata;
    cdata.confs.resize(1);
    cdata.setX(0, *coordinates);
    computePotential(nullptr, p, top, atoms, &cdata, epsilonr, false);
    *forces = cdata.f(0).toRVec();
    energies->clear();
    double epot = 0;
    for(const auto &e : cdata.energies.energies(0))
    {
        if (InteractionType::EPOT == e.first)
        {
            epot = e.second;
        }
        else
        {
            energies->insert({ e.first, e.second });
        }
    }
    return epot;
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
        const gmx::test::FloatingPointTolerance tolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 1e-6));
        checker_.setDefaultTolerance(tolerance);
        std::vector<double> q = { 1, 1, -1, -2 };
        for (int i = 0; i < 4; i++)
        {
            // Note the row value must be <= 2 in order to work for Slaters
            atoms_.emplace_back("C", "C", "C", ActParticle::Atom, 6, 12, q[i], 2);
        }
    }

    void testPot(Potential                     p,
                 const TopologyEntryVector    &top,
                 const std::vector<gmx::RVec> *coordinates,
                 const ForceField             *pd)
    {
        // Generate force array with zeroes
        gmx::RVec                         fzero = { 0, 0, 0 };
        std::vector<gmx::RVec>            forces(coordinates->size(), fzero);
        // Energy map
        std::map<InteractionType, double> energies;

        double epsilonr = pd ? pd->epsilonr() : 1;
        auto epot = computeOne(p, top, atoms_, coordinates, &forces, &energies, epsilonr);
        // In single precision the coordinates are rounded, by up to
        // GMX_FLOAT_EPS times their magnitude. That changes the forces by
        // the stiffness times this, also where the reference force is zero,
        // and the energies by the forces times this. For two atoms the
        // stiffness along z is estimated from the forces at displaced
        // coordinates.
        double xround    = 0;
        double stiffness = 0;
        if (!GMX_DOUBLE)
        {
            for (const auto &xi : *coordinates)
            {
                for (int m = 0; m < DIM; m++)
                {
                    xround += GMX_FLOAT_EPS*std::abs(xi[m]);
                }
            }
            if (coordinates->size() == 2)
            {
                double h = 1.0/1024;
                double fz[2];
                for (int k = 0; k < 2; k++)
                {
                    std::vector<gmx::RVec>            fk(coordinates->size(), fzero);
                    std::map<InteractionType, double> ek;
                    auto                              newx = *coordinates;
                    newx[1][ZZ] += (2*k-1)*h;
                    (void) computeOne(p, top, atoms_, &newx, &fk, &ek, epsilonr);
                    fz[k] = fk[1][ZZ];
                }
                stiffness = std::abs(fz[1] - fz[0])/(2*h);
            }
        }
        double fmax = 0;
        for (const auto &f : forces)
        {
            for (int m = 0; m < DIM; m++)
            {
                fmax = std::max(fmax, static_cast<double>(std::abs(f[m])));
            }
        }
        // Tolerance with an absolute part for the rounded coordinates in
        // single precision, for values of both float and double type
        const auto roundedTolerance = [](double absolute)
        {
            return gmx::test::FloatingPointTolerance(absolute, absolute, 1e-5, 1e-5,
                                                     UINT64_MAX, UINT64_MAX, false);
        };
        if (!GMX_DOUBLE)
        {
            checker_.setDefaultTolerance(roundedTolerance(2*fmax*xround));
        }
        checker_.checkReal(epot, "Epot");
        if (!GMX_DOUBLE)
        {
            // Components that cancel have errors relative to the largest
            // force. Forces that cancel to zero, as beyond the turning point
            // of a cubic bond, have errors relative to the energy over the
            // distance between the particles, while the reference is the
            // rounding error of double precision.
            double extent = 0;
            for (const auto &xi : *coordinates)
            {
                for (const auto &xj : *coordinates)
                {
                    const gmx::RVec dx = xi - xj;
                    extent = std::max(extent, static_cast<double>(norm(dx)));
                }
            }
            double cancel = extent > 0 ? 4*GMX_FLOAT_EPS*std::abs(epot)/extent : 0;
            checker_.setDefaultTolerance(roundedTolerance(2*stiffness*xround + 1e-6*fmax + cancel));
        }
        checker_.checkSequence(forces.begin(), forces.end(), "Forces");
        if (!GMX_DOUBLE)
        {
            checker_.setDefaultTolerance(roundedTolerance(2*fmax*xround));
        }
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
            for (int m = 0; m < DIM; m++)
            {
                coords2[i][m] = -(*coordinates)[i][m];
            }
        }
        // Energy map
        std::map<InteractionType, double> energies2;
        (void) computeOne(p, top, atoms_, &coords2, &forces2, &energies2, epsilonr);
        double ediff = 1e-8;
        // First, check that both maps have the same number of entries
        ASSERT_EQ(energies.size(), energies2.size());
        // For each entry in energies2, check that the key exists in energies and compare values
        for (const auto &ee2 : energies2)
        {
            const auto it = energies.find(ee2.first);
            ASSERT_NE(it, energies.end());
            EXPECT_TRUE(std::abs(ee2.second - it->second) < ediff);
        }
        // Also check the other direction: no extra entries in energies
        for (const auto &ee : energies)
        {
            const auto it2 = energies2.find(ee.first);
            ASSERT_NE(it2, energies2.end());
            EXPECT_TRUE(std::abs(ee.second - it2->second) < ediff);
        }
        double fdiff = 1e-4;
        for (size_t i = 0; i < forces2.size(); i++)
        {
            for (int m = 0; m < DIM; m++)
            {
                EXPECT_TRUE(std::abs(forces2[i][m] + forces[i][m]) < fdiff);
            }
        }
#if GMX_DOUBLE
        // A power of two, such that the displaced coordinates are exact on
        // the grid of the integer coordinates, see coordinate.h
        double dx = 1.0/1048576;
        double ener[2];
        rvec   dxvector[2];
        for (int k = 0; k < 2; k++)
        {
            // Additional force and energy structures
            std::vector<gmx::RVec>            forces(coordinates->size(), fzero);
            std::map<InteractionType, double> energies;
            // Modify the coordinates
            auto newx = *coordinates;
            newx[1][ZZ] += (2*k-1)*dx;
            rvec_sub((*coordinates)[1], newx[1], dxvector[k]);
            ener[k] = computeOne(p, top, atoms_, &newx, &forces, &energies, epsilonr);
        }
        // Compute numerical derivative
        double fz = -(ener[1]-ener[0])/(norm(dxvector[0])+norm(dxvector[1]));
#else
        // Central differences of fourth order, with a larger step, since the
        // rounding error of the energies is divided by it, and a power of
        // two, such that the displaced coordinates are exact
        double h    = 1.0/1024;
        double emax = std::abs(epot);
        const auto shifted = [&](double dz)
        {
            std::vector<gmx::RVec>            forces(coordinates->size(), fzero);
            std::map<InteractionType, double> energies;
            auto                              newx = *coordinates;
            newx[1][ZZ] += dz;
            double e = computeOne(p, top, atoms_, &newx, &forces, &energies, epsilonr);
            emax     = std::max(emax, std::abs(e));
            return e;
        };
        double e1 = shifted(h) - shifted(-h);
        double e2 = shifted(2*h) - shifted(-2*h);
        double e4 = shifted(4*h) - shifted(-4*h);
        double fz = -(8*e1 - e2)/(12*h);
        // The truncation error, of order h^4, follows from the same
        // stencil with step 2h (Richardson).
        double fz2h  = -(8*e2 - e4)/(24*h);
        double trunc = std::abs(fz2h - fz)/15;
        // The rounding errors of the energies enter fz with the sum of the
        // absolute weights 18/(12h). They are up to about 8 GMX_FLOAT_EPS
        // each, from exp and pow and from terms that partly cancel. The
        // reference force has the error from the rounded coordinates.
        double fzround = trunc + 12*GMX_FLOAT_EPS*emax/h + 2*stiffness*xround;
#endif
        double toler = 3e-5;
        // It should be identical to the forces computed at the central point
        // A difference of two energies, in which the last bits of the
        // energies, that can differ between platforms, are amplified
#if GMX_DOUBLE
        checker_.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-7));
#else
        checker_.setDefaultTolerance(gmx::test::FloatingPointTolerance(fzround, fzround, 1e-4, 1e-4,
                                                                       UINT64_MAX, UINT64_MAX, false));
#endif
        checker_.checkReal(fz, "fz");
        checker_.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 1e-6));
        double diff = std::abs(fz - forces[1][ZZ]);
        double f2   = std::abs(fz + forces[1][ZZ]);
#if !GMX_DOUBLE
        // Differences within the rounding errors are accepted
        if (diff < fzround)
        {
            diff = 0;
        }
#endif
        if ((f2 > 0) && (diff >= toler))
        {
            // Make it relative
            diff = diff/f2;
        }
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
    const auto x   = coord.coordsBuilder();
    auto top = pot.topBuilder();
    if (Potential::COULOMB_POINT == pot.potential || Potential::COULOMB_GAUSSIAN == pot.potential ||
        Potential::COULOMB_SLATER == pot.potential)
    {
        // The widths of the charge distributions belong to the atoms
        for(size_t i = 0; i < par.params.size(); i++)
        {
            atoms_[i].setZeta(par.params[i]);
        }
    }
    else if (!top.empty())
    {
        top[0]->setParams(par.params);
    }
    testPot(pot.potential, top, &x, nullptr);
}

// ============================================================
// Coordinate parameter sets
// ============================================================

/*! \brief Five standard 2-atom pair coordinate sets with slightly different separations
 *         and off-axis displacements, all with charge1 = +1.
 *
 *  - "r0p5":  atoms along Z, separation 0.5 nm  (original geometry)
 *  - "r0p4":  atoms along Z, separation 0.4 nm  (shorter bond)
 *  - "r0p6d": atoms at separation ~0.6 nm with a small X offset (off-axis)
 *  - "r0p3":  atoms along Z, separation 0.3 nm  (even shorter distance)
 *  - "r0pxd": atoms offset in X and Y, separation approximately 0.5 nm
 */
static const std::vector<ForceComputerCoordParams> c_stdPairCoords = {
    { "r0p5",  []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0,    0,    0.5   } }; }, 1.0 },
    { "r0p4",  []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0,    0,    0.4   } }; }, 1.0 },
    { "r0p6d", []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0.05, 0.05, 0.59  } }; }, 1.0 },
    { "r0p3",  []() { return std::vector<gmx::RVec>{ { 0, 0, 0 }, { 0.0,  0.0,  0.3   } }; }, 1.0 },
    { "r0pxd", []() { return std::vector<gmx::RVec>{ { 1, 0, 0 }, { 0.95, 0.05, 0.495 } }; }, 1.0 }
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
    top.emplace_back(AtomPair(0, 1));
    return top;
}

//! Angle topology (bond 0-1, bond 1-2)
static TopologyEntryVector makeAngleTop()
{
    const Bond b1(0, 1, 1);
    const Bond b2(1, 2, 1);
    TopologyEntryVector top;
    top.emplace_back(Angle(b1, b2));
    return top;
}

//! Improper dihedral topology (bonds 0-1, 0-2, 0-3)
static TopologyEntryVector makeImproperTop()
{
    const Bond b1(0, 1, 1);
    const Bond b2(0, 2, 1);
    const Bond b3(0, 3, 1);
    TopologyEntryVector top;
    top.emplace_back(Improper(b1, b2, b3));
    return top;
}

//! Proper dihedral topology (bonds 0-1, 1-2, 2-3)
static TopologyEntryVector makeProperTop()
{
    const Bond b1(0, 1, 1);
    const Bond b2(1, 2, 1);
    const Bond b3(2, 3, 1);
    TopologyEntryVector top;
    top.emplace_back(Proper(b1, b2, b3));
    return top;
}

//! Single-atom topology (for position restraints)
static TopologyEntryVector makeSingleAtomTop()
{
    TopologyEntryVector top;
    top.emplace_back(SingleAtom(0));
    return top;
}

/*! \brief Tang-Toennies damping function f_n(x) = 1 - exp(-x) sum_k=0^n x^k/k!
 * Summed to convergence from the terms beyond n for small x, independent
 * of the implementation in the force computer.
 */
static double ttDamping(int n, double x)
{
    double term = 1;
    for (int k = 1; k <= n; k++)
    {
        term *= x/k;
    }
    if (x > 20)
    {
        double sum = 1, t = 1;
        for (int k = 1; k <= n; k++)
        {
            t   *= x/k;
            sum += t;
        }
        return 1 - std::exp(-x)*sum;
    }
    double tail = 0;
    for (int k = n+1; term > 1e-20*tail; k++)
    {
        term *= x/k;
        tail += term;
    }
    return std::exp(-x)*tail;
}

/* The damping functions for small br are differences of nearly equal
 * numbers, which lost up to 1e-5 relative in double precision for br = 0.5.
 * The energies and forces are compared to the damping functions summed to
 * convergence, for each of C6, C8 and C10 and for br on both sides of the
 * switch to the sum of the terms beyond 2m.
 */
TEST(TangToenniesDampingTest, ShortDistances)
{
    std::vector<ActAtom> atoms;
    for (int i = 0; i < 2; i++)
    {
        atoms.emplace_back("C", "C", "C", ActParticle::Atom, 6, 12, 0, 2);
    }
    // Exact in single precision
    double r = 0.25;
    const std::vector<gmx::RVec> x = { { 0, 0, 0 }, { 0, 0, static_cast<real>(r) } };
    // In single precision, the direct formula for br >= ttTailMax has an
    // error of a few GMX_FLOAT_EPS divided by the damping function there.
    double toler = gmx::test::precisionTolerance(1e-13, 3*GMX_FLOAT_EPS/ttDamping(10, ttTailMax));
    for (int m = 3; m <= 5; m++)
    {
        int n = 2*m;
        for (double br : { 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 5.9, 6.1, 8.0, 12.0 })
        {
            double b = br/r;
            // Parameter order: tt2bA, tt2bBexch, tt2bBdisp, tt2bC6, tt2bC8, tt2bC10
            std::vector<double> params = { 0, 10, b, 0, 0, 0 };
            double              cn     = 1e-3;
            params[m] = cn;
            auto top = makePairTop();
            top[0]->setParams(params);
            std::vector<gmx::RVec>            forces(2, { 0, 0, 0 });
            std::map<InteractionType, double> energies;
            double epot = computeOne(Potential::TT2b, top, atoms, &x, &forces, &energies, 1);
            double fn   = ttDamping(n, br);
            double eref = -cn*fn/std::pow(r, n);
            // dE/dr with dfn/dx = exp(-x) x^n/n!
            double dfn  = std::exp(-br)*std::pow(br, n)/std::tgamma(n+1);
            double fref = cn*(b*dfn - n*fn/r)/std::pow(r, n);
            // The two terms of the force cancel for small br
            double fmag = cn*(b*dfn + n*fn/r)/std::pow(r, n);
            SCOPED_TRACE(gmx::formatString("n %d br %g", n, br));
            EXPECT_NEAR(epot, eref, toler*std::abs(eref));
            EXPECT_NEAR(forces[1][ZZ], fref, toler*fmag);
            EXPECT_NEAR(forces[0][ZZ], -fref, toler*fmag);
        }
    }
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
//! COULOMB_SLATER with both widths zero, two point charges
static const ForceComputerPotParams c_coulSlaterPointPot { "COULOMB_SLATER_POINT", Potential::COULOMB_SLATER,    makePairTop };
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
//! Coulomb parameters: both widths zero, point charges
static const ForceComputerParamParams c_coulPointParams  { "", { 0, 0 } };
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
static const ForceComputerParamParams c_huaBondsParams  { "", { 100, 20, 0.4, 0.1 } };    // DE, b, bondlength, c (huaDE, huaB, huaLENGTH, huaC)
static const ForceComputerParamParams c_morseBondsParams{ "", { 12, 100, 40, 0.6 } };     // beta, DE, D0, length
static const ForceComputerParamParams c_harmAngleParams { "", { 100, 100 } };     // kT, angle
static const ForceComputerParamParams c_ubAngleParams   { "", { 100, 100, 0.33, 20 } };   // kT, angle, r13, kUB
static const ForceComputerParamParams c_linAngleParams  { "", { 0.5, 10000 } };   // a, kLin
static const ForceComputerParamParams c_harmDihParams   { "", { 10 } };           // kPhi
static const ForceComputerParamParams c_properDihParams { "", { 30, 10, 4 } };    // angle, kP, mult
static const ForceComputerParamParams c_fourierDihParams{ "", { 1, 2, -3, 4, -5, 6 } };
static const ForceComputerParamParams c_posResParams    { "", { 100, 1, 0, 0, 0 } }; // k, R0, center

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
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(LJ12_6, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_lj12_6Pot),
                                            ::testing::Values(c_lj12_6Params),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(LJ12_6_4, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_lj12_6_4Pot),
                                            ::testing::Values(c_lj12_6_4Params),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(LJ14_7, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_lj14_7Pot),
                                            ::testing::Values(c_lj14_7Params),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

// --- Buckingham potentials ---
INSTANTIATE_TEST_CASE_P(BUCKINGHAM, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_buckinghamPot),
                                            ::testing::Values(c_buckinghamParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(WANG_BUCKINGHAM, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_wangBuckPot),
                                            ::testing::Values(c_wangBuckParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(GENERALIZED_BUCKINGHAM, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_genBuckPot),
                                            ::testing::Values(c_genBuckParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

// --- Coulomb potentials (each also tested with multiple coord/charge variants) ---

//! COULOMB_POINT: standard coords only
INSTANTIATE_TEST_CASE_P(COULOMB_POINT, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulPointPot),
                                            ::testing::Values(c_coulZeta10Params),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

//! COULOMB_GAUSSIAN: standard coords + negative-charge + zero-distance variants
static const std::vector<ForceComputerCoordParams> c_coulGaussCoords = []() {
    std::vector<ForceComputerCoordParams> v(c_stdPairCoords.begin(), c_stdPairCoords.end());
    v.push_back(c_coulNegCoord);
    v.push_back(c_coulZeroCoord);
    return v;
}();

INSTANTIATE_TEST_CASE_P(CoulGauss, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulGaussPot),
                                            ::testing::Values(c_coulZeta10Params),
                                            ::testing::ValuesIn(c_coulGaussCoords)),
                         fcTestName);

//! COULOMB_SLATER zeta=10: standard coords + charge/distance variants
static const std::vector<ForceComputerCoordParams> c_coulSlaterCoords = []() {
    std::vector<ForceComputerCoordParams> v(c_stdPairCoords.begin(), c_stdPairCoords.end());
    v.push_back(c_coulNegCoord);
    v.push_back(c_coulNeg2Coord);
    v.push_back(c_coulCloseCoord);
    v.push_back(c_coulCloseNegCoord);
    v.push_back(c_coulZeroCoord);
    v.push_back(c_coulZeroNegCoord);
    return v;
}();

INSTANTIATE_TEST_CASE_P(CoulSlater, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulSlaterPot),
                                            ::testing::Values(c_coulZeta10Params),
                                            ::testing::ValuesIn(c_coulSlaterCoords)),
                         fcTestName);

//! COULOMB_SLATER zeta=0: standard coords + negative-charge variant
static const std::vector<ForceComputerCoordParams> c_coulSlaterZ0Coords = []() {
    std::vector<ForceComputerCoordParams> v(c_stdPairCoords.begin(), c_stdPairCoords.end());
    v.push_back(c_coulNegCoord);
    return v;
}();

INSTANTIATE_TEST_CASE_P(CoulSlaterZ0, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulSlaterZ0Pot),
                                            ::testing::Values(c_coulZeta0Params),
                                            ::testing::ValuesIn(c_coulSlaterZ0Coords)),
                         fcTestName);

//! COULOMB_SLATER with two point charges, which interact as with COULOMB_POINT
INSTANTIATE_TEST_CASE_P(CoulSlaterPoint, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_coulSlaterPointPot),
                                            ::testing::Values(c_coulPointParams),
                                            ::testing::ValuesIn(c_coulSlaterZ0Coords)),
                         fcTestName);

// --- TT2b family: ONE potential tested with FIVE different parameter sets ---
//! Demonstrates the Combine benefit: adding a new param set here tests TT2b under new conditions
INSTANTIATE_TEST_CASE_P(TT2b, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_tt2bPot),
                                            ::testing::ValuesIn(c_tt2bParamSets),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

// --- Other pair potentials ---
INSTANTIATE_TEST_CASE_P(SLATER_ISA_TT, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_slaterIsaTtPot),
                                            ::testing::Values(c_slaterIsaTtParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(TANG_TOENNIES, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_tangToennPot),
                                            ::testing::Values(c_tangToennParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(BORN_MAYER, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_bornMayerPot),
                                            ::testing::Values(c_expParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(SLATER_ISA, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_slaterIsaPot),
                                            ::testing::Values(c_expParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(MACDANIEL_SCHMIDT, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_macdanielPot),
                                            ::testing::Values(c_macdanielParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(POLARIZATION, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_polPot),
                                            ::testing::Values(c_polParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

// --- Bond potentials ---
INSTANTIATE_TEST_CASE_P(HARMONIC_BONDS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_harmBondsPot),
                                            ::testing::Values(c_harmBondsParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(CUBIC_BONDS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_cubicBondsPot),
                                            ::testing::Values(c_cubicBondsParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(HUA_BONDS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_huaBondsPot),
                                            ::testing::Values(c_huaBondsParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
                         fcTestName);

INSTANTIATE_TEST_CASE_P(MORSE_BONDS, ForceComputerImplementationTest,
                         ::testing::Combine(::testing::Values(c_morseBondsPot),
                                            ::testing::Values(c_morseBondsParams),
                                            ::testing::ValuesIn(c_stdPairCoords)),
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

TEST_F(ForceComputerImplementationTest, GaussianEpsilonR)
{
    // Get forcefield
    auto pd  = getForceField("ACS-pg-vs2");
    ForceField mypd = *pd;
    mypd.setEpsilonr(2.0);
    const auto &pot   = c_coulGaussPot;
    const auto &par   = c_coulZeta10Params;
    const auto &coord = c_stdPairCoords[0];
    atoms_[1].setCharge(coord.charge1);
    const auto x   = coord.coordsBuilder();
    const auto top = pot.topBuilder();
    // The widths of the charge distributions belong to the atoms
    for(size_t i = 0; i < par.params.size(); i++)
    {
        atoms_[i].setZeta(par.params[i]);
    }
    testPot(pot.potential, top, &x, &mypd);
}

}  // namespace
}  // namespace alexandria
