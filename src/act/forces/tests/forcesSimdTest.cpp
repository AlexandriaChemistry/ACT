/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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
/*! \internal \file
 * \brief
 * Implements test of the potential routines in forces_simd.cpp
 *
 * Each potential is computed for a number of conformations, both with
 * the conformations in the outer loop and in the inner, vectorized, loop.
 * The results should be the same, and the forces should be consistent
 * with the energies.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_forces_tests
 */
#include "actpre.h"

#include "../forces_simd.h"

#include <cmath>
#include <set>

#include <gtest/gtest.h>

#include "act/basics/act_particle.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/potential.h"
#include "act/forces/simd_information.h"
#include "act/topology/actatom.h"
#include "act/topology/coords_energy.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

//! Atom indices and parameters of one interaction
typedef struct
{
    std::vector<int>    atoms;
    std::vector<double> params;
} Interaction;

class ForcesSimdTest : public ::testing::Test
{
protected:
    gmx::test::TestReferenceData    refData_;
    gmx::test::TestReferenceChecker checker_;
    //! Number of conformations, odd to test the remainder of vectorized loops
    const size_t                    nconf_ = 7;
    //! The atoms
    std::vector<ActAtom>            atoms_;
    //! The conformations
    std::vector<std::vector<gmx::RVec>> x_;
    //! More conformations, to test the unrolled loops for all SIMD widths
    std::vector<std::vector<gmx::RVec>> xMany_;
    //! Number of conformations in xMany_, more than four times the widest SIMD plus a remainder
    const size_t                        nMany_ = 37;

    ForcesSimdTest() : checker_(refData_.rootChecker())
    {
        // Vector math libraries may differ in the last bits
#if GMX_DOUBLE
        checker_.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1, 1e-10));
#else
        checker_.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1, 1e-5));
#endif
        // A chain of four atoms, two atoms nearby for non-bonded
        // interactions and a shell close to atom 4.
        std::vector<gmx::RVec> x0 = {
            { 0.000,  0.000,  0.000 },
            { 0.150,  0.000,  0.000 },
            { 0.200,  0.140,  0.000 },
            { 0.320,  0.180,  0.090 },
            { 0.050,  0.300,  0.180 },
            { 0.400, -0.150,  0.250 },
            { 0.058,  0.304,  0.183 }
        };
        std::vector<double> q   = { 0.4, -0.3, 0.2, -0.3, 0.5, -0.5, -1.2 };
        std::vector<int>    row = { 1, 2, 2, 3, 1, 2, 2 };
        for (size_t i = 0; i < x0.size(); i++)
        {
            atoms_.push_back(ActAtom("A", "C", "c", ActParticle::Atom,
                                     6, 12.0, q[i], row[i]));
        }
        // Deterministic perturbations of up to 0.01 nm
        unsigned int seed = 1993;
        for (size_t c = 0; c < nconf_ + nMany_; c++)
        {
            auto x = x0;
            for (auto &xi : x)
            {
                for (int m = 0; m < DIM; m++)
                {
                    seed   = 1103515245*seed + 12345;
                    xi[m] += 0.01*(((seed >> 8) % 2001)/1000.0 - 1.0);
                }
            }
            if (c < nconf_)
            {
                x_.push_back(x);
            }
            else
            {
                xMany_.push_back(x);
            }
        }
    }

    //! Make topology entries from the interactions
    TopologyEntryVector makeEntries(const std::vector<Interaction> &inter)
    {
        TopologyEntryVector entries;
        for (const auto &i : inter)
        {
            TopologyEntry te;
            for (auto a : i.atoms)
            {
                te.addAtom(a);
            }
            te.setParams(i.params);
            entries.push_back(BaseContainer<TopologyEntry>(std::move(te)));
        }
        return entries;
    }

    //! Make the conformations from coordinates
    std::vector<CoordsEnergy> makeCoords(const std::vector<std::vector<gmx::RVec>> &x)
    {
        std::vector<CoordsEnergy> coords(x.size());
        for (size_t c = 0; c < x.size(); c++)
        {
            coords[c].x = x[c];
        }
        return coords;
    }

    //! \return the potential energy of conformation c
    double epot(const std::vector<CoordsEnergy> &coords, size_t c)
    {
        return coords[c].eMap.find(InteractionType::EPOT)->second.eact();
    }

    //! \return the potential energy of all conformations with the given atoms
    std::vector<double> energies(Potential                       mypot,
                                 const std::vector<Interaction> &inter,
                                 const std::vector<ActAtom>     &atoms)
    {
        auto coords = makeCoords(x_);
        computePotential(nullptr, mypot, makeEntries(inter), atoms, &coords, 1, true);
        std::vector<double> e;
        for (size_t c = 0; c < nconf_; c++)
        {
            e.push_back(epot(coords, c));
        }
        return e;
    }

    /*! \brief Test Coulomb with second charges on some of the atoms
     * Besides the usual tests the energy is compared to the sum of the
     * four combinations of the charges of both atoms in each pair, computed
     * with one charge per atom. The atoms are first or second in the pairs,
     * such that the charges can be chosen per role.
     * \param[in] mypot The potential
     * \param[in] inter Pairs with the widths of the first charges
     */
    void testSecondCharges(Potential mypot, const std::vector<Interaction> &inter)
    {
        const std::vector<double> qvs   = { -0.8, 0.6, 0.0, -0.4, -1.1, 0.9, 0.0 };
        const std::vector<double> zvs   = { 8.0, 11.0, 0.0, 7.0, 9.0, 10.0, 0.0 };
        const std::vector<int>    rowvs = { 1, 2, 1, 1, 2, 1, 1 };
        auto original = atoms_;
        for (size_t i = 0; i < atoms_.size(); i++)
        {
            atoms_[i].setChargeVS(qvs[i], zvs[i], rowvs[i]);
        }
        test(mypot, inter);
        auto eSecond = energies(mypot, inter, atoms_);
        atoms_       = original;
        // The atoms that are first and second in the pairs
        std::set<int> first, second;
        for (const auto &i : inter)
        {
            first.insert(i.atoms[0]);
            second.insert(i.atoms[1]);
        }
        for (auto a : first)
        {
            ASSERT_EQ(0u, second.count(a));
        }
        std::vector<double> eSum(nconf_, 0.0);
        for (int a = 0; a < 2; a++)
        {
            for (int c = 0; c < 2; c++)
            {
                // Atoms with only the first or only the second charge
                std::vector<ActAtom> atoms;
                for (size_t k = 0; k < original.size(); k++)
                {
                    const auto &o   = original[k];
                    bool        vs  = (first.count(k) && a == 1) || (second.count(k) && c == 1);
                    atoms.push_back(ActAtom(o.name(), o.element(), o.ffType(), o.pType(), o.atomicNumber(),
                                            o.mass(), vs ? qvs[k] : o.charge(), vs ? rowvs[k] : o.row()));
                }
                auto myinter = inter;
                for (auto &i : myinter)
                {
                    if (a == 1)
                    {
                        i.params[0] = zvs[i.atoms[0]];
                    }
                    if (c == 1)
                    {
                        i.params[1] = zvs[i.atoms[1]];
                    }
                }
                auto e = energies(mypot, myinter, atoms);
                for (size_t conf = 0; conf < nconf_; conf++)
                {
                    eSum[conf] += e[conf];
                }
            }
        }
        for (size_t conf = 0; conf < nconf_; conf++)
        {
            EXPECT_NEAR(eSum[conf], eSecond[conf], 1e-10*(1 + std::abs(eSum[conf]))) << "conformation " << conf;
        }
    }

    /*! \brief Compare the routines for many conformations
     * With explicit SIMD, if available, with the loops vectorized by the
     * compiler and one conformation at a time. The energies and forces
     * have to be the same, except for the last bits from vector math
     * libraries. Future optimizations of one of the routines should not
     * change the results.
     * \param[in] mypot    The potential
     * \param[in] entries  The interactions
     * \param[in] epsilonr The relative dielectric constant
     * \param[in] tolSame  The tolerance
     */
    void compareRoutines(Potential                              mypot,
                         const TopologyEntryVector             &entries,
                         double                                 epsilonr,
                         const gmx::test::FloatingPointTolerance &tolSame)
    {
        auto cOne      = makeCoords(xMany_);
        computePotential(nullptr, mypot, entries, atoms_, &cOne, epsilonr, false);
        bool explicitSimd = setUseExplicitSimd(false);
        auto cCompiler = makeCoords(xMany_);
        computePotential(nullptr, mypot, entries, atoms_, &cCompiler, epsilonr, true);
        setUseExplicitSimd(true);
        auto cExplicit = makeCoords(xMany_);
        computePotential(nullptr, mypot, entries, atoms_, &cExplicit, epsilonr, true);
        // Without forces, as in training, the energies are exactly the same,
        // with explicit SIMD and with the loops vectorized by the compiler
        for (bool useExplicit : { false, true })
        {
            setUseExplicitSimd(useExplicit);
            auto        cNoForces = makeCoords(xMany_);
            computePotential(nullptr, mypot, entries, atoms_, &cNoForces, epsilonr, true, false);
            const auto &withForces = useExplicit ? cExplicit : cCompiler;
            for (size_t c = 0; c < nMany_; c++)
            {
                ASSERT_EQ(withForces[c].eMap.size(), cNoForces[c].eMap.size());
                for (const auto &e : withForces[c].eMap)
                {
                    auto ei = cNoForces[c].eMap.find(e.first);
                    ASSERT_TRUE(ei != cNoForces[c].eMap.end());
                    EXPECT_EQ(e.second.eact(), ei->second.eact())
                        << interactionTypeToString(e.first) << " without forces, conformation " << c
                        << (useExplicit ? ", explicit SIMD" : ", vectorized by the compiler");
                }
            }
        }
        setUseExplicitSimd(explicitSimd);
        for (const auto *other : { &cCompiler, &cExplicit })
        {
            for (size_t c = 0; c < nMany_; c++)
            {
                SCOPED_TRACE(gmx::formatString("conformation %zu, %s", c,
                                               other == &cCompiler ? "vectorized by the compiler" : "explicit SIMD"));
                ASSERT_EQ(cOne[c].eMap.size(), (*other)[c].eMap.size());
                for (const auto &e : cOne[c].eMap)
                {
                    auto ei = (*other)[c].eMap.find(e.first);
                    ASSERT_TRUE(ei != (*other)[c].eMap.end());
                    EXPECT_REAL_EQ_TOL(e.second.eact(), ei->second.eact(), tolSame);
                }
                for (size_t a = 0; a < atoms_.size(); a++)
                {
                    for (int m = 0; m < DIM; m++)
                    {
                        EXPECT_REAL_EQ_TOL(cOne[c].fcalc[a][m], (*other)[c].fcalc[a][m], tolSame);
                    }
                }
            }
        }
    }
    void test(Potential mypot, const std::vector<Interaction> &inter)
    {
        auto   entries  = makeEntries(inter);
        double epsilonr = 1;
        // Conformations in the outer loop
        auto   cOuter   = makeCoords(x_);
        computePotential(nullptr, mypot, entries, atoms_, &cOuter, epsilonr, false);
        // Conformations in the inner loop
        auto   cInner   = makeCoords(x_);
        computePotential(nullptr, mypot, entries, atoms_, &cInner, epsilonr, true);
#if GMX_DOUBLE
        // Vector math libraries may differ from the scalar ones in the last bits
        auto   tolSame  = gmx::test::relativeToleranceAsFloatingPoint(1, 1e-12);
        double h        = 1e-6;
        double tolFD    = 1e-5;
#else
        auto   tolSame  = gmx::test::relativeToleranceAsFloatingPoint(1, 1e-5);
        double h        = 1e-3;
        double tolFD    = 1e-2;
#endif
        compareRoutines(mypot, entries, epsilonr, tolSame);
        for (size_t c = 0; c < nconf_; c++)
        {
            ASSERT_EQ(cOuter[c].eMap.size(), cInner[c].eMap.size());
            for (const auto &e : cOuter[c].eMap)
            {
                auto ei = cInner[c].eMap.find(e.first);
                ASSERT_TRUE(ei != cInner[c].eMap.end());
                EXPECT_REAL_EQ_TOL(e.second.eact(), ei->second.eact(), tolSame);
            }
            for (size_t a = 0; a < atoms_.size(); a++)
            {
                for (int m = 0; m < DIM; m++)
                {
                    EXPECT_REAL_EQ_TOL(cOuter[c].fcalc[a][m], cInner[c].fcalc[a][m], tolSame);
                }
            }
        }
        // Check the forces against finite differences of the energy,
        // using all conformations at once.
        double fmax = 0;
        for (const auto &ce : cInner)
        {
            for (const auto &f : ce.fcalc)
            {
                for (int m = 0; m < DIM; m++)
                {
                    fmax = std::max(fmax, std::abs(f[m]));
                }
            }
        }
        for (size_t a = 0; a < atoms_.size(); a++)
        {
            for (int m = 0; m < DIM; m++)
            {
                auto xplus  = x_;
                auto xminus = x_;
                for (size_t c = 0; c < nconf_; c++)
                {
                    xplus[c][a][m]  += h;
                    xminus[c][a][m] -= h;
                }
                auto cPlus  = makeCoords(xplus);
                auto cMinus = makeCoords(xminus);
                computePotential(nullptr, mypot, entries, atoms_, &cPlus, epsilonr, true);
                computePotential(nullptr, mypot, entries, atoms_, &cMinus, epsilonr, true);
                for (size_t c = 0; c < nconf_; c++)
                {
                    double fnum = -(epot(cPlus, c) - epot(cMinus, c))/(2*h);
                    EXPECT_NEAR(fnum, cInner[c].fcalc[a][m], tolFD*std::max(1.0, fmax))
                        << "conformation " << c << " atom " << a << " dim " << m;
                }
            }
        }
        // Store the energies
        for (size_t c = 0; c < nconf_; c++)
        {
            for (const auto &e : cInner[c].eMap)
            {
                std::string label = gmx::formatString("%s-%zu",
                                                      interactionTypeToString(e.first).c_str(), c);
                checker_.checkReal(e.second.eact(), label.c_str());
            }
        }
    }
};

//! Non-bonded pairs
static std::vector<Interaction> pairs(const std::vector<double> &params)
{
    return { { { 0, 4 }, params }, { { 1, 5 }, params },
             { { 2, 4 }, params }, { { 3, 5 }, params },
             { { 0, 5 }, params } };
}

//! Bonds
static std::vector<Interaction> bonds(const std::vector<double> &params)
{
    return { { { 0, 1 }, params }, { { 1, 2 }, params }, { { 2, 3 }, params } };
}

//! Angles
static std::vector<Interaction> angles(const std::vector<double> &params)
{
    return { { { 0, 1, 2 }, params }, { { 1, 2, 3 }, params } };
}

//! Dihedrals
static std::vector<Interaction> dihedrals(const std::vector<double> &params)
{
    return { { { 0, 1, 2, 3 }, params }, { { 3, 2, 1, 0 }, params } };
}

TEST_F(ForcesSimdTest, PositionRestraint)
{
    std::vector<double> params = { 5000, 0.2 };
    test(Potential::POSITION_RESTRAINT,
         { { { 1 }, params }, { { 2 }, params }, { { 3 }, params }, { { 4 }, params } });
}

TEST_F(ForcesSimdTest, LJ12_6)
{
    test(Potential::LJ12_6, pairs({ 0.3, 0.5 }));
}

TEST_F(ForcesSimdTest, LJ12_6_4)
{
    test(Potential::LJ12_6_4, pairs({ 0.3, 0.5, 0.0002 }));
}

TEST_F(ForcesSimdTest, LJ8_6)
{
    test(Potential::LJ8_6, pairs({ 0.3, 0.5 }));
}

TEST_F(ForcesSimdTest, LJ14_7)
{
    auto inter = pairs({ 0.3, 0.5, 0.12, 0.07 });
    // No parameters for one pair
    inter[1].params = { 0, 0, 0.12, 0.07 };
    test(Potential::LJ14_7, inter);
}

TEST_F(ForcesSimdTest, BornMayer)
{
    test(Potential::BORN_MAYER, pairs({ 20000, 30, 0 }));
}

TEST_F(ForcesSimdTest, SlaterISA)
{
    test(Potential::SLATER_ISA, pairs({ 20000, 30, 0 }));
}

TEST_F(ForcesSimdTest, DoubleExponential)
{
    test(Potential::MACDANIEL_SCHMIDT, pairs({ 30000, 10000, 30 }));
}

TEST_F(ForcesSimdTest, WangBuckingham)
{
    auto inter = pairs({ 0.3, 0.5, 10 });
    inter[1].params = { 0.3, 0, 10 };
    test(Potential::WANG_BUCKINGHAM, inter);
}

TEST_F(ForcesSimdTest, Buckingham)
{
    auto inter = pairs({ 20000, 30, 0.002 });
    inter[1].params = { 20000, 30, 0 };
    test(Potential::BUCKINGHAM, inter);
}

TEST_F(ForcesSimdTest, TangToennies)
{
    test(Potential::TANG_TOENNIES, pairs({ 20000, 30, 0.002, 0.0002, 0.00002 }));
}

TEST_F(ForcesSimdTest, TT2b)
{
    test(Potential::TT2b, pairs({ 20000, 30, 25, 0.002, 0.0002, 0.00002 }));
}

TEST_F(ForcesSimdTest, SlaterISATT)
{
    test(Potential::SLATER_ISA_TT, pairs({ 20000, 30, 25, 0.002, 0.0002, 0.00002 }));
}

TEST_F(ForcesSimdTest, GeneralizedBuckingham)
{
    auto inter = pairs({ 0.35, 0.5, 12, 8 });
    inter[1].params = { 0.35, 0.5, 0, 8 };
    test(Potential::GENERALIZED_BUCKINGHAM, inter);
}

TEST_F(ForcesSimdTest, CoulombPoint)
{
    test(Potential::COULOMB_POINT, pairs({ 0, 0 }));
}

TEST_F(ForcesSimdTest, CoulombGaussian)
{
    test(Potential::COULOMB_GAUSSIAN, pairs({ 12, 9 }));
}

TEST_F(ForcesSimdTest, CoulombSlater)
{
    test(Potential::COULOMB_SLATER, pairs({ 25, 20 }));
}

TEST_F(ForcesSimdTest, CoulombPointSecondCharges)
{
    testSecondCharges(Potential::COULOMB_POINT, pairs({ 0, 0 }));
}

TEST_F(ForcesSimdTest, CoulombGaussianSecondCharges)
{
    testSecondCharges(Potential::COULOMB_GAUSSIAN, pairs({ 12, 9 }));
}

TEST_F(ForcesSimdTest, CoulombSlaterSecondCharges)
{
    testSecondCharges(Potential::COULOMB_SLATER, pairs({ 25, 20 }));
}

TEST_F(ForcesSimdTest, Polarization)
{
    std::vector<double> params = { 0.001, 0.01, 1e8 };
    test(Potential::POLARIZATION, { { { 4, 6 }, params }, { { 0, 6 }, { 0, 0.01, 1e8 } } });
}

TEST_F(ForcesSimdTest, HarmonicBonds)
{
    test(Potential::HARMONIC_BONDS, bonds({ 200000, 0.15, -400 }));
}

TEST_F(ForcesSimdTest, CubicBonds)
{
    test(Potential::CUBIC_BONDS, bonds({ 0.15, 0.2, 400000, 400 }));
}

TEST_F(ForcesSimdTest, MorseBonds)
{
    test(Potential::MORSE_BONDS, bonds({ 20, 400, 10, 0.15 }));
}

TEST_F(ForcesSimdTest, HuaBonds)
{
    test(Potential::HUA_BONDS, bonds({ 400, 20, 0.15, 0.1 }));
}

TEST_F(ForcesSimdTest, HarmonicAngles)
{
    test(Potential::HARMONIC_ANGLES, angles({ 400, 1.9 }));
}

TEST_F(ForcesSimdTest, LinearAngles)
{
    test(Potential::LINEAR_ANGLES, angles({ 0.4, 300000 }));
}

TEST_F(ForcesSimdTest, UreyBradleyAngles)
{
    test(Potential::UREY_BRADLEY_ANGLES, angles({ 400, 1.9, 0.24, 20000 }));
}

TEST_F(ForcesSimdTest, ProperDihedrals)
{
    test(Potential::PROPER_DIHEDRALS, dihedrals({ 0.3, 5, 3 }));
}

TEST_F(ForcesSimdTest, ImproperDihedrals)
{
    test(Potential::HARMONIC_DIHEDRALS, dihedrals({ 50 }));
}

TEST_F(ForcesSimdTest, FourierDihedrals)
{
    test(Potential::FOURIER_DIHEDRALS, dihedrals({ 3, -2, 1, 4, -1, 0.5 }));
}

TEST_F(ForcesSimdTest, FourierDihedralsFewParameters)
{
    test(Potential::FOURIER_DIHEDRALS, dihedrals({ 3, -2, 1.5 }));
}

} // namespace

} // namespace alexandria
