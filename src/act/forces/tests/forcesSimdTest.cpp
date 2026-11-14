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
#include <cstdio>
#include <map>
#include <set>

#include <gtest/gtest.h>

#include "act/basics/act_particle.h"
#include "act/coulombintegrals/gaussian_integrals.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/potential.h"
#include "act/forces/forces_explicit_simd.h"
#include "act/forces/forces_gpu.h"
#include "act/forces/simd_information.h"
#include "act/topology/actatom.h"
#include "act/topology/coords_energy.h"
#include "act/gpu/gpu.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"

#include "testutils/gputests.h"
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
        const std::vector<gmx::RVec> x0 = {
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
            atoms_.emplace_back("A", "C", "c", ActParticle::Atom,
                                6, 12.0, q[i], row[i]);
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
            entries.emplace_back(std::move(te));
        }
        return entries;
    }

    /*! \brief Give the atoms the widths of their charge distributions
     * The Coulomb interactions list the widths of their first and second
     * atom, which have to be the same in all interactions of an atom.
     * \param[in] inter The interactions
     */
    void setWidths(const std::vector<Interaction> &inter)
    {
        std::map<int, double> zeta;
        for (const auto &i : inter)
        {
            ASSERT_EQ(2u, i.atoms.size());
            ASSERT_EQ(2u, i.params.size());
            for (int k = 0; k < 2; k++)
            {
                const auto z = zeta.find(i.atoms[k]);
                if (zeta.end() != z)
                {
                    ASSERT_EQ(z->second, i.params[k]) << "atom " << i.atoms[k];
                }
                zeta[i.atoms[k]] = i.params[k];
                atoms_[i.atoms[k]].setZeta(i.params[k]);
            }
        }
    }

    //! Make the conformations from coordinates
    ConformationData makeCoords(const std::vector<std::vector<gmx::RVec>> &x)
    {
        ConformationData coords;
        coords.confs.resize(x.size());
        for (size_t c = 0; c < x.size(); c++)
        {
            coords.setX(c, x[c]);
        }
        return coords;
    }

    //! \return the potential energy of conformation c
    double epot(const ConformationData &coords, size_t c)
    {
        return coords.energies.energies(c).at(InteractionType::EPOT);
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
        std::vector<double> eSum(nconf_, 0.0), eAbs(nconf_, 0.0);
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
                    atoms.emplace_back(o.name(), o.element(), o.ffType(), o.pType(), o.atomicNumber(),
                                       o.mass(), vs ? qvs[k] : o.charge(), vs ? rowvs[k] : o.row());
                    atoms.back().setZeta(vs ? zvs[k] : o.zeta());
                }
                auto e = energies(mypot, inter, atoms);
                for (size_t conf = 0; conf < nconf_; conf++)
                {
                    eSum[conf] += e[conf];
                    eAbs[conf] += std::abs(e[conf]);
                }
            }
        }
        // In single precision relative to the terms, which may cancel
        double toler = gmx::test::precisionTolerance(1e-10, 1e-6);
        for (size_t conf = 0; conf < nconf_; conf++)
        {
            EXPECT_NEAR(eSum[conf], eSecond[conf], toler*(1 + eAbs[conf])) << "conformation " << conf;
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
    /*! \brief The tolerance for comparing forces
     * In single precision a force is accurate relative to the largest
     * force, since the terms of the forces cancel. The explicit SIMD
     * computes them in single precision as well, with other roundings and
     * math functions, as in compareSinglePrecision.
     * \param[in] f       The forces
     * \param[in] tolSame The tolerance in double precision
     * \return The tolerance
     */
    gmx::test::FloatingPointTolerance forceTolerance(const std::vector<gmx::RVec>           &f,
                                                     const gmx::test::FloatingPointTolerance &tolSame)
    {
        if (GMX_DOUBLE)
        {
            return tolSame;
        }
        double fmax = 1;
        for (const auto &fi : f)
        {
            for (int m = 0; m < DIM; m++)
            {
                fmax = std::max<double>(fmax, std::abs(fi[m]));
            }
        }
        return gmx::test::relativeToleranceAsFloatingPoint(fmax, 1e-4);
    }

    /*! \brief The coordinates of xMany_ in the layout of explicit SIMD
     * \tparam    T      The floating point type
     * \param[out] stride The number of conformations, including the padding
     *                    to a multiple of the SIMD width
     * \return the coordinates in units of coordinateUnit, see computeExplicitSimd
     */
    template <typename T>
    std::vector<Coordinate> manyArrays(size_t *stride)
    {
        size_t width = explicitSimdWidth<T>();
        *stride      = ((nMany_ + width - 1)/width)*width;
        std::vector<Coordinate> x(DIM*atoms_.size()*(*stride), 0);
        for (size_t c = 0; c < *stride; c++)
        {
            // The padding repeats the last conformation
            const auto &xc = xMany_[std::min(c, nMany_-1)];
            for (size_t a = 0; a < atoms_.size(); a++)
            {
                for (int m = 0; m < DIM; m++)
                {
                    x[(DIM*a+m)*(*stride) + c] = toCoordinate(xc[a][m]);
                }
            }
        }
        return x;
    }

    /*! \brief Compute one potential for xMany_ with explicit SIMD
     * \tparam    T        The floating point type of the kernel
     * \param[in]  mypot    The potential
     * \param[in]  entries  The interactions
     * \param[in]  epsilonr Relative dielectric constant
     * \param[out] stride   The stride of the arrays, see computeExplicitSimd
     * \param[out] f        The forces
     * \param[out] energies The energy terms, in double precision for both types
     * \return false if the potential is not available
     */
    template <typename T>
    bool explicitSimd(Potential                  mypot,
                      const TopologyEntryVector &entries,
                      double                     epsilonr,
                      size_t                    *stride,
                      std::vector<T>            *f,
                      std::vector<double>       *energies)
    {
        auto x = manyArrays<T>(stride);
        f->assign(x.size(), 0);
        energies->assign(explicitSimdEnergyTerms*(*stride), 0);
        return computeExplicitSimd(mypot, entries, atoms_, epsilonr, *stride, *stride,
                                   x.data(), f->data(), energies->data(), true);
    }

    /*! \brief Compare explicit SIMD in single precision with double precision
     * \param[in] mypot    The potential
     * \param[in] entries  The interactions
     * \param[in] epsilonr Relative dielectric constant
     */
    void compareSinglePrecision(Potential                  mypot,
                                const TopologyEntryVector &entries,
                                double                     epsilonr)
    {
        if (!useExplicitSimd())
        {
            return;
        }
        size_t              sd, sf;
        std::vector<double> fd, ed, ef;
        std::vector<float>  ff;
        bool                haveDouble = explicitSimd(mypot, entries, epsilonr, &sd, &fd, &ed);
        ASSERT_EQ(haveDouble, explicitSimd(mypot, entries, epsilonr, &sf, &ff, &ef));
        if (!haveDouble)
        {
            return;
        }
        // Within the precision of floats, the forces relative to the largest
        for (size_t c = 0; c < nMany_; c++)
        {
            double fmax = 0;
            for (size_t i = 0; i < DIM*atoms_.size(); i++)
            {
                fmax = std::max(fmax, std::abs(fd[i*sd + c]));
            }
            for (size_t i = 0; i < DIM*atoms_.size(); i++)
            {
                EXPECT_NEAR(fd[i*sd + c], ff[i*sf + c], 1e-4*(1 + fmax))
                    << potentialToString(mypot) << " conformation " << c << " force " << i;
            }
            for (int k = 0; k < explicitSimdEnergyTerms; k++)
            {
                EXPECT_NEAR(ed[k*sd + c], ef[k*sf + c], 2e-5*(1 + std::abs(ed[k*sd + c])))
                    << potentialToString(mypot) << " conformation " << c << " energy term " << k;
            }
        }
        if (Potential::COULOMB_POINT == mypot || Potential::COULOMB_GAUSSIAN == mypot)
        {
            // The energies for charges of one, all four terms of each interaction
            std::vector<double>        kd;
            std::vector<float>         kf;
            const std::vector<unsigned char> terms(entries.size(), 0xf);
            coulombKernels(entries, terms, &sd, &kd);
            coulombKernels(entries, terms, &sf, &kf);
            for (size_t row = 0; row < 4*entries.size(); row++)
            {
                for (size_t c = 0; c < nMany_; c++)
                {
                    EXPECT_NEAR(kd[row*sd + c], kf[row*sf + c], 2e-5*(1 + std::abs(kd[row*sd + c])))
                        << "Coulomb kernel " << row << " conformation " << c;
                }
            }
        }
    }

    /*! \brief Compare the GPU with explicit SIMD in double precision for xMany_
     * Only call this if gpu::useGpu() and useExplicitSimd() return true.
     * The GPU computes in single precision and sums the energies in double
     * precision, such that the tolerances are those of compareSinglePrecision.
     * \param[in] mypot    The potential
     * \param[in] entries  The interactions
     * \param[in] epsilonr Relative dielectric constant
     */
    void compareGpu(Potential                  mypot,
                    const TopologyEntryVector &entries,
                    double                     epsilonr)
    {
        size_t              sd, sf;
        std::vector<double> fd, ed;
        ASSERT_TRUE(explicitSimd(mypot, entries, epsilonr, &sd, &fd, &ed));
        // The interactions of the tests are less work than the GPU gets otherwise
        size_t              minWork = setGpuMinimumWork(0);
        auto                x = manyArrays<float>(&sf);
        std::vector<float>  ff(x.size(), 0);
        std::vector<double> ef(explicitSimdEnergyTerms*sf, 0);
        ASSERT_TRUE(computeGpu(mypot, entries, atoms_, epsilonr, sf, sf, x.data(), ff.data(), ef.data(), true));
        // Without forces only the energies are computed
        std::vector<float>  fNone(x.size(), 0);
        std::vector<double> eNone(explicitSimdEnergyTerms*sf, 0);
        ASSERT_TRUE(computeGpu(mypot, entries, atoms_, epsilonr, sf, sf, x.data(), fNone.data(), eNone.data(), false));
        // The same results in every run, and for fewer conformations
        std::vector<float>  fAgain(x.size(), 0);
        std::vector<double> eAgain(explicitSimdEnergyTerms*sf, 0);
        size_t              nfew = sf/2;
        ASSERT_TRUE(computeGpu(mypot, entries, atoms_, epsilonr, sf, nfew, x.data(), fAgain.data(), eAgain.data(), true));
        setGpuMinimumWork(minWork);
        for (size_t c = 0; c < nfew; c++)
        {
            for (size_t i = 0; i < DIM*atoms_.size(); i++)
            {
                EXPECT_EQ(ff[i*sf + c], fAgain[i*sf + c]) << "conformation " << c << " force " << i;
            }
            EXPECT_EQ(ef[c], eAgain[c]) << "conformation " << c;
        }
        for (size_t c = 0; c < nMany_; c++)
        {
            double fmax = 0;
            for (size_t i = 0; i < DIM*atoms_.size(); i++)
            {
                fmax = std::max(fmax, std::abs(fd[i*sd + c]));
            }
            for (size_t i = 0; i < DIM*atoms_.size(); i++)
            {
                EXPECT_NEAR(fd[i*sd + c], ff[i*sf + c], 1e-4*(1 + fmax))
                    << potentialToString(mypot) << " conformation " << c << " force " << i;
                EXPECT_EQ(0, fNone[i*sf + c]) << "force " << i << " without forces";
            }
            for (int k = 0; k < explicitSimdEnergyTerms; k++)
            {
                EXPECT_NEAR(ed[k*sd + c], ef[k*sf + c], 2e-5*(1 + std::abs(ed[k*sd + c])))
                    << potentialToString(mypot) << " conformation " << c << " energy term " << k;
                // The energy of an interaction with several terms can differ in
                // the last bits without forces
                EXPECT_NEAR(ed[k*sd + c], eNone[k*sf + c], 2e-5*(1 + std::abs(ed[k*sd + c])))
                    << potentialToString(mypot) << " conformation " << c << " energy term " << k
                    << " without forces";
            }
        }
    }

    /*! \brief The Coulomb energies for charges of one of xMany_ with explicit SIMD
     * \tparam    T       The floating point type of the kernel
     * \param[in]  entries The interactions
     * \param[in]  terms   The terms to compute, see coulombKernelsExplicitSimd
     * \param[out] stride  The stride of the arrays
     * \param[out] kernel  The energies
     */
    template <typename T>
    void coulombKernels(const TopologyEntryVector        &entries,
                        const std::vector<unsigned char> &terms,
                        size_t                           *stride,
                        std::vector<T>                   *kernel)
    {
        auto x = manyArrays<T>(stride);
        kernel->assign(4*entries.size()*(*stride), 0);
        ASSERT_TRUE(coulombKernelsExplicitSimd(entries, atoms_, *stride, *stride, x.data(), terms.data(),
                                               kernel->data()));
    }

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
        // Without forces, as in training, the energies are the same, with
        // explicit SIMD and with the loops vectorized by the compiler. They
        // can differ in the last bits, since compilers may combine a
        // multiplication and an addition to a fused multiply-add
        // differently without the operations for the forces.
        for (bool useExplicit : { false, true })
        {
            setUseExplicitSimd(useExplicit);
            auto        cNoForces = makeCoords(xMany_);
            computePotential(nullptr, mypot, entries, atoms_, &cNoForces, epsilonr, true, false);
            const auto &withForces = useExplicit ? cExplicit : cCompiler;
            for (size_t c = 0; c < nMany_; c++)
            {
                const auto eWith = withForces.energies.energies(c);
                auto eNo   = cNoForces.energies.energies(c);
                ASSERT_EQ(eWith.size(), eNo.size());
                for (const auto &e : eWith)
                {
                    const auto ei = eNo.find(e.first);
                    ASSERT_TRUE(ei != eNo.end());
                    EXPECT_REAL_EQ_TOL(e.second, ei->second, tolSame)
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
                const auto eOne   = cOne.energies.energies(c);
                auto eOther = other->energies.energies(c);
                ASSERT_EQ(eOne.size(), eOther.size());
                for (const auto &e : eOne)
                {
                    const auto ei = eOther.find(e.first);
                    ASSERT_TRUE(ei != eOther.end());
                    EXPECT_REAL_EQ_TOL(e.second, ei->second, tolSame);
                }
                const auto tolForce = forceTolerance(cOne.f(c).toRVec(), tolSame);
                for (size_t a = 0; a < atoms_.size(); a++)
                {
                    for (int m = 0; m < DIM; m++)
                    {
                        EXPECT_REAL_EQ_TOL(cOne.f(c)[a][m], other->f(c)[a][m], tolForce);
                    }
                }
            }
        }
    }
    void test(Potential mypot, const std::vector<Interaction> &inter)
    {
        const auto entries  = makeEntries(inter);
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
        const auto tolSame  = gmx::test::relativeToleranceAsFloatingPoint(1, 1e-5);
        double h        = 1e-3;
        double tolFD    = 1e-2;
#endif
        compareRoutines(mypot, entries, epsilonr, tolSame);
        compareSinglePrecision(mypot, entries, epsilonr);
        for (size_t c = 0; c < nconf_; c++)
        {
            const auto eOuter = cOuter.energies.energies(c);
            auto eInner = cInner.energies.energies(c);
            ASSERT_EQ(eOuter.size(), eInner.size());
            for (const auto &e : eOuter)
            {
                const auto ei = eInner.find(e.first);
                ASSERT_TRUE(ei != eInner.end());
                EXPECT_REAL_EQ_TOL(e.second, ei->second, tolSame);
            }
            const auto tolForce = forceTolerance(cOuter.f(c).toRVec(), tolSame);
            for (size_t a = 0; a < atoms_.size(); a++)
            {
                for (int m = 0; m < DIM; m++)
                {
                    EXPECT_REAL_EQ_TOL(cOuter.f(c)[a][m], cInner.f(c)[a][m], tolForce);
                }
            }
        }
        // Check the forces against finite differences of the energy,
        // using all conformations at once.
        double fmax = 0;
        for (size_t c = 0; c < cInner.nconf(); c++)
        {
            for (const auto &f : cInner.f(c).toRVec())
            {
                for (int m = 0; m < DIM; m++)
                {
                    fmax = std::max<double>(fmax, std::abs(f[m]));
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
                    EXPECT_NEAR(fnum, cInner.f(c)[a][m], tolFD*std::max(1.0, fmax))
                        << "conformation " << c << " atom " << a << " dim " << m;
                }
            }
        }
        // Store the energies
        for (size_t c = 0; c < nconf_; c++)
        {
            for (const auto &e : cInner.energies.energies(c))
            {
                const std::string label = gmx::formatString("%s-%zu",
                                                            interactionTypeToString(e.first).c_str(), c);
                checker_.checkReal(e.second, label.c_str());
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
    // Force constant, radius and the center of the restraints
    std::vector<double> params = { 5000, 0.2, 0, 0, 0 };
    test(Potential::POSITION_RESTRAINT,
         { { { 1 }, params }, { { 2 }, params }, { { 3 }, params }, { { 4 }, params } });
}

// Test: Position restraints around a center that is not the origin
TEST_F(ForcesSimdTest, PositionRestraintCenter)
{
    std::vector<double> params = { 5000, 0.2, 0.1, -0.05, 0.2 };
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

TEST_F(ForcesSimdTest, TT2bShortDistances)
{
    // With this exponent of the damping, br is between 4.2 and 7.4, such
    // that the damping functions are computed both from the terms beyond
    // 2m (br < ttTailMax) and as 1 - exp(-br) sum_k=0^2m (br)^k/k!
    test(Potential::TT2b, pairs({ 20000, 30, 15, 0.002, 0.0002, 0.00002 }));
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
    setWidths(pairs({ 0, 0 }));
    test(Potential::COULOMB_POINT, pairs({ 0, 0 }));
}

TEST_F(ForcesSimdTest, CoulombGaussian)
{
    setWidths(pairs({ 12, 9 }));
    test(Potential::COULOMB_GAUSSIAN, pairs({ 12, 9 }));
}

// Test: The force between Gaussian charges at short distances, where the
// formula is a difference of nearly equal terms, for each routine. The
// expected value is from the series of the force divided by the distance,
// with all the terms needed in double precision.
TEST_F(ForcesSimdTest, CoulombGaussianShortDistance)
{
    // Combined width of 12 and 9 is 7.2/nm
    const double zeta = 12.0*9.0/std::sqrt(12.0*12.0 + 9.0*9.0);
    const double qq   = ONE_4PI_EPS0*atoms_[0].charge()*atoms_[1].charge();
    setWidths({ { { 0, 1 }, { 12, 9 } } });
    const auto   entries = makeEntries({ { { 0, 1 }, { 12, 9 } } });
    const std::vector<double> dist = { 1e-5, 1e-4, 1e-3, 1e-2, 0.03, 0.1 };
    std::vector<std::vector<gmx::RVec>> x;
    for (auto r : dist)
    {
        auto xc = x_[0];
        xc[0]   = { 0, 0, 0 };
        xc[1]   = { static_cast<real>(r), 0, 0 };
        x.push_back(xc);
    }
    bool explicitSimd = setUseExplicitSimd(false);
    for (int routine = 0; routine < 3; routine++)
    {
        setUseExplicitSimd(routine == 2);
        // In single precision the explicit SIMD computes the difference of
        // the terms in single precision as well
        const double toler = gmx::test::precisionTolerance(1e-10, routine == 2 ? 1e-5 : 1e-6);
        auto coords = makeCoords(x);
        computePotential(nullptr, Potential::COULOMB_GAUSSIAN, entries, atoms_, &coords, 1, routine > 0);
        for (size_t c = 0; c < dist.size(); c++)
        {
            // The distance as stored
            double r  = coords.x(c)[1][XX];
            double x2 = gmx::square(zeta*r);
            // (erf(x)/x - 2 exp(-x^2)/sqrt(pi))/x^2 as a series in x^2
            double series = 0, term = 1;
            for (int n = 0; n < 30; n++)
            {
                int k   = n + 1;
                series += ((k % 2) ? -1.0 : 1.0)*(1.0/(std::tgamma(k + 1.0)*(2*k + 1)) - 1.0/std::tgamma(k + 1.0))*term;
                term   *= x2;
            }
            double fOverR = qq*zeta*zeta*zeta*2*invsqrt_pi*series;
            // The force on the first atom, at the origin, points away from the second
            double expect = -fOverR*r;
            EXPECT_NEAR(expect, coords.f(c)[0][XX], toler*std::abs(expect)) << "routine " << routine << " r " << r;
            EXPECT_NEAR(-expect, coords.f(c)[1][XX], toler*std::abs(expect)) << "routine " << routine << " r " << r;
        }
    }
    setUseExplicitSimd(explicitSimd);
}

// Test: Coulomb on the GPU, if there is one, against explicit SIMD in
// double precision. Without a GPU these tests pass without testing.
TEST_F(ForcesSimdTest, CoulombPointGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && useExplicitSimd(), "no GPU or no explicit SIMD"))
    {
        return;
    }
    setWidths(pairs({ 0, 0 }));
    compareGpu(Potential::COULOMB_POINT, makeEntries(pairs({ 0, 0 })), 1);
}

TEST_F(ForcesSimdTest, CoulombGaussianGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && useExplicitSimd(), "no GPU or no explicit SIMD"))
    {
        return;
    }
    setWidths(pairs({ 12, 9 }));
    compareGpu(Potential::COULOMB_GAUSSIAN, makeEntries(pairs({ 12, 9 })), 1);
    compareGpu(Potential::COULOMB_GAUSSIAN, makeEntries(pairs({ 12, 9 })), 2.5);
}

TEST_F(ForcesSimdTest, CoulombGaussianSecondChargesGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && useExplicitSimd(), "no GPU or no explicit SIMD"))
    {
        return;
    }
    const std::vector<double> qvs   = { -0.8, 0.6, 0.0, -0.4, -1.1, 0.9, 0.0 };
    const std::vector<double> zvs   = { 8.0, 11.0, 0.0, 7.0, 9.0, 10.0, 0.0 };
    for (size_t i = 0; i < atoms_.size(); i++)
    {
        atoms_[i].setChargeVS(qvs[i], zvs[i], 1);
    }
    setWidths(pairs({ 12, 9 }));
    compareGpu(Potential::COULOMB_GAUSSIAN, makeEntries(pairs({ 12, 9 })), 1);
}

// Test: the other pair potentials on the GPU, as CoulombGaussianGpu
TEST_F(ForcesSimdTest, PairPotentialsGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && useExplicitSimd(), "no GPU or no explicit SIMD"))
    {
        return;
    }
    auto lj147 = pairs({ 0.3, 0.5, 0.12, 0.07 });
    lj147[1].params = { 0, 0, 0.12, 0.07 };
    auto wbh   = pairs({ 0.3, 0.5, 10 });
    wbh[1].params = { 0.3, 0, 10 };
    const std::vector<std::pair<Potential, std::vector<Interaction>>> tests = {
        { Potential::LJ12_6,            pairs({ 0.3, 0.5 }) },
        { Potential::LJ14_7,            lj147 },
        { Potential::WANG_BUCKINGHAM,   wbh },
        { Potential::MACDANIEL_SCHMIDT, pairs({ 30000, 10000, 30 }) },
        { Potential::POLARIZATION,      { { { 4, 6 }, { 0.001, 0.01, 1e8 } }, { { 0, 6 }, { 0, 0.01, 1e8 } } } },
        { Potential::HARMONIC_BONDS,    bonds({ 200000, 0.15, -400 }) },
        { Potential::CUBIC_BONDS,       bonds({ 0.15, 0.2, 400000, 400 }) },
        { Potential::MORSE_BONDS,       bonds({ 20, 400, 10, 0.15 }) },
        { Potential::HUA_BONDS,         bonds({ 400, 20, 0.15, 0.1 }) }
    };
    for (const auto &t : tests)
    {
        SCOPED_TRACE(potentialToString(t.first));
        compareGpu(t.first, makeEntries(t.second), 1);
    }
}

// Test: several molecules computed at once on the GPU, with different
// strides and numbers of slots, give the same results as each alone
TEST_F(ForcesSimdTest, CoulombGaussianBatchGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && useExplicitSimd(), "no GPU or no explicit SIMD"))
    {
        return;
    }
    setWidths(pairs({ 12, 9 }));
    auto   entries = makeEntries(pairs({ 12, 9 }));
    size_t stride;
    auto   x1      = manyArrays<float>(&stride);
    // The second molecule has a larger stride and its conformations in
    // the reverse order
    size_t stride2 = stride + 8;
    std::vector<Coordinate> x2(DIM*atoms_.size()*stride2, 0);
    for (size_t row = 0; row < DIM*atoms_.size(); row++)
    {
        for (size_t c = 0; c < stride; c++)
        {
            x2[row*stride2 + c] = x1[row*stride + stride - 1 - c];
        }
    }
    size_t minWork = setGpuMinimumWork(0);
    // Each alone
    std::vector<float>  f1(x1.size(), 0), f2(x2.size(), 0);
    std::vector<double> e1(explicitSimdEnergyTerms*stride, 0), e2(explicitSimdEnergyTerms*stride2, 0);
    ASSERT_TRUE(computeGpu(Potential::COULOMB_GAUSSIAN, entries, atoms_, 1, stride, stride,
                           x1.data(), f1.data(), e1.data(), true));
    ASSERT_TRUE(computeGpu(Potential::COULOMB_GAUSSIAN, entries, atoms_, 1, stride2, stride - 3,
                           x2.data(), f2.data(), e2.data(), true));
    // Together, the first molecule twice
    GpuBatch batch;
    batch.add(Potential::COULOMB_GAUSSIAN, entries, atoms_, 1, stride, stride, x1.data());
    batch.add(Potential::COULOMB_GAUSSIAN, entries, atoms_, 1, stride2, stride - 3, x2.data());
    batch.add(Potential::COULOMB_GAUSSIAN, entries, atoms_, 1, stride, stride, x1.data());
    ASSERT_TRUE(batch.compute(true));
    setGpuMinimumWork(minWork);
    ASSERT_EQ(3u, batch.size());
    for (size_t k = 0; k < batch.size(); k++)
    {
        const auto         &fk = (k == 1) ? f2 : f1;
        const auto         &ek = (k == 1) ? e2 : e1;
        size_t              sk = (k == 1) ? stride2 : stride;
        size_t              nk = (k == 1) ? stride - 3 : stride;
        std::vector<float>  f(fk.size(), 0);
        std::vector<double> e(ek.size(), 0);
        batch.addResults(k, f.data(), e.data());
        for (size_t c = 0; c < nk; c++)
        {
            EXPECT_EQ(ek[c], e[c]) << "molecule " << k << " conformation " << c;
            for (size_t row = 0; row < DIM*atoms_.size(); row++)
            {
                EXPECT_EQ(fk[row*sk + c], f[row*sk + c]) << "molecule " << k << " conformation " << c << " row " << row;
            }
        }
    }
}

// Test: Gaussian charges on top of each other, as shells on their cores,
// and nearly so, where the series of the force is used
TEST_F(ForcesSimdTest, CoulombGaussianOverlapGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && useExplicitSimd(), "no GPU or no explicit SIMD"))
    {
        return;
    }
    for (size_t c = 0; c < nMany_; c++)
    {
        // The same position, or at up to 0.03 nm
        auto &xc = xMany_[c];
        xc[4]    = xc[0];
        xc[4][XX] += 0.001*(c % 31);
    }
    setWidths(pairs({ 12, 9 }));
    compareGpu(Potential::COULOMB_GAUSSIAN, makeEntries(pairs({ 12, 9 })), 1);
}

TEST_F(ForcesSimdTest, CoulombSlater)
{
    setWidths(pairs({ 25, 20 }));
    test(Potential::COULOMB_SLATER, pairs({ 25, 20 }));
}

TEST_F(ForcesSimdTest, CoulombPointSecondCharges)
{
    setWidths(pairs({ 0, 0 }));
    testSecondCharges(Potential::COULOMB_POINT, pairs({ 0, 0 }));
}

TEST_F(ForcesSimdTest, CoulombGaussianSecondCharges)
{
    setWidths(pairs({ 12, 9 }));
    testSecondCharges(Potential::COULOMB_GAUSSIAN, pairs({ 12, 9 }));
}

TEST_F(ForcesSimdTest, CoulombSlaterSecondCharges)
{
    setWidths(pairs({ 25, 20 }));
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
