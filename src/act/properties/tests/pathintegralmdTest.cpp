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

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <random>
#include <vector>

#include "act/properties/fragmentwall.h"
#include "act/properties/pathintegralmd.h"
#include "act/properties/pathintegralthermo.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"

namespace alexandria
{

namespace
{

/*! \brief Compounds of a few particles at random positions
 * \param[in]  nfrag     The number of compounds
 * \param[in]  spread    How far the compounds are apart in nm
 * \param[out] atomStart The first particle of each compound and the number of particles
 * \param[out] mass      The masses, the last particle of each compound has none
 * \param[out] x         The coordinates
 */
void makeFragments(size_t                  nfrag,
                   double                  spread,
                   std::vector<size_t>    *atomStart,
                   std::vector<double>    *mass,
                   std::vector<gmx::DVec> *x)
{
    std::mt19937 gen(1234 + nfrag);
    const auto   uniform = [&gen]() { return static_cast<double>(gen() % 100000)/100000.0; };
    atomStart->assign(1, 0);
    mass->clear();
    x->clear();
    for(size_t f = 0; f < nfrag; f++)
    {
        const gmx::DVec center = { spread*uniform(), spread*uniform(), spread*uniform() };
        size_t          natom  = 3 + f;
        for(size_t i = 0; i < natom; i++)
        {
            // A particle without mass, as a shell or a virtual site
            mass->push_back(i + 1 == natom ? 0.0 : 1.0 + 15*uniform());
            x->emplace_back(center[XX] + 0.1*uniform(), center[YY] + 0.1*uniform(), center[ZZ] + 0.1*uniform());
        }
        atomStart->push_back(atomStart->back() + natom);
    }
}

class FragmentWallTest : public ::testing::TestWithParam<size_t>
{
};

TEST_P(FragmentWallTest, ForcesAreTheDerivativeOfTheEnergy)
{
    std::vector<size_t>    atomStart;
    std::vector<double>    mass;
    std::vector<gmx::DVec> x, f, fdum;
    makeFragments(GetParam(), 1.0, &atomStart, &mass, &x);
    const FragmentWall wall(atomStart, mass, 0.15, 800);
    ASSERT_TRUE(wall.active());
    ASSERT_EQ(GetParam(), wall.nfragments());
    double energy = wall.compute(x, &f);
    ASSERT_GT(energy, 0);
    ASSERT_EQ(x.size(), f.size());
    double    delta = 1e-5;
    gmx::DVec ftot  = { 0, 0, 0 };
    double    fmax  = 0;
    for(size_t i = 0; i < x.size(); i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            auto xp = x;
            auto xm = x;
            xp[i][m] += delta;
            xm[i][m] -= delta;
            double fnum = -(wall.compute(xp, &fdum) - wall.compute(xm, &fdum))/(2*delta);
            EXPECT_NEAR(fnum, f[i][m], 1e-6*(1 + std::abs(fnum))) << "particle " << i;
            ftot[m] += f[i][m];
            fmax     = std::max(fmax, std::abs(f[i][m]));
        }
    }
    EXPECT_GT(fmax, 1.0);
    for(int m = 0; m < DIM; m++)
    {
        // No force on all compounds together
        EXPECT_NEAR(0.0, ftot[m], 1e-10*fmax);
    }
    // The acceleration is the same for all particles of a compound, and a
    // particle without mass has no force
    for(size_t k = 0; k + 1 < atomStart.size(); k++)
    {
        size_t i0 = atomStart[k];
        for(size_t i = i0; i < atomStart[k + 1]; i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                if (mass[i] > 0)
                {
                    EXPECT_NEAR(f[i0][m]/mass[i0], f[i][m]/mass[i], 1e-10*fmax);
                }
                else
                {
                    EXPECT_EQ(0.0, f[i][m]);
                }
            }
        }
    }
}

TEST_P(FragmentWallTest, NoForceInsideTheWall)
{
    std::vector<size_t>    atomStart;
    std::vector<double>    mass;
    std::vector<gmx::DVec> x, f;
    makeFragments(GetParam(), 1.0, &atomStart, &mass, &x);
    const FragmentWall wall(atomStart, mass, 5.0, 800);
    ASSERT_TRUE(wall.active());
    EXPECT_GT(wall.distance(x), 0.0);
    EXPECT_LT(wall.distance(x), 5.0);
    EXPECT_EQ(0.0, wall.compute(x, &f));
    for(const auto &fi : f)
    {
        EXPECT_EQ(0.0, fi[XX]*fi[XX] + fi[YY]*fi[YY] + fi[ZZ]*fi[ZZ]);
    }
}

INSTANTIATE_TEST_CASE_P(Fragments, FragmentWallTest, ::testing::Values(2, 3, 5));

TEST(FragmentWallSimpleTest, TwoAtoms)
{
    // Two atoms 1 nm apart, a wall at 0.75 nm with force constant 100
    const std::vector<gmx::DVec> x = { { 0, 0, 0 }, { 0, 0.6, 0.8 } };
    std::vector<gmx::DVec>       f;
    const FragmentWall wall({ 0, 1, 2 }, { 4.0, 40.0 }, 0.75, 100);
    EXPECT_NEAR(1.0, wall.distance(x), 1e-14);
    EXPECT_NEAR(0.5*100*0.25*0.25, wall.compute(x, &f), 1e-12);
    // Towards each other with k times the distance beyond the wall
    EXPECT_NEAR(0.0, f[0][XX], 1e-12);
    EXPECT_NEAR(25*0.6, f[0][YY], 1e-12);
    EXPECT_NEAR(25*0.8, f[0][ZZ], 1e-12);
    EXPECT_NEAR(-25*0.6, f[1][YY], 1e-12);
    EXPECT_NEAR(-25*0.8, f[1][ZZ], 1e-12);
}

TEST(FragmentWallSimpleTest, DistanceOfMoreFragments)
{
    // Three atoms of equal mass on a line: the largest distance to the
    // center of mass of all, which is at 1
    const std::vector<gmx::DVec> x = { { 0, 0, 0 }, { 0.5, 0, 0 }, { 2.5, 0, 0 } };
    const FragmentWall wall({ 0, 1, 2, 3 }, { 2.0, 2.0, 2.0 }, 0, 0);
    EXPECT_NEAR(1.5, wall.distance(x), 1e-14);
    EXPECT_FALSE(wall.active());
}

TEST(FragmentWallSimpleTest, NoWall)
{
    const std::vector<gmx::DVec> x = { { 0, 0, 0 }, { 0, 0, 3 } };
    std::vector<gmx::DVec>       f;
    // No radius, no force constant, or one compound
    for(const auto &wall : { FragmentWall({ 0, 1, 2 }, { 1.0, 1.0 }, 0, 100),
                             FragmentWall({ 0, 1, 2 }, { 1.0, 1.0 }, 1, 0),
                             FragmentWall({ 0, 2 }, { 1.0, 1.0 }, 1, 100) })
    {
        EXPECT_FALSE(wall.active());
        EXPECT_EQ(0.0, wall.compute(x, &f));
        ASSERT_EQ(2U, f.size());
        EXPECT_EQ(0.0, f[1][ZZ]);
    }
    EXPECT_EQ(0.0, FragmentWall({ 0, 2 }, { 1.0, 1.0 }, 1, 100).distance(x));
}

TEST(FragmentWallSimpleTest, WrongInput)
{
    // A negative radius or force constant
    EXPECT_THROW(FragmentWall({ 0, 1, 2 }, { 1.0, 1.0 }, -1, 100), gmx::InvalidInputError);
    EXPECT_THROW(FragmentWall({ 0, 1, 2 }, { 1.0, 1.0 }, 1, -100), gmx::InvalidInputError);
    // The compounds do not have all particles
    EXPECT_THROW(FragmentWall({ 0, 1 }, { 1.0, 1.0 }, 1, 100), gmx::InvalidInputError);
    // A compound without mass, and one without particles
    EXPECT_THROW(FragmentWall({ 0, 1, 2 }, { 1.0, 0.0 }, 1, 100), gmx::InvalidInputError);
    EXPECT_THROW(FragmentWall({ 0, 2, 2 }, { 1.0, 1.0 }, 1, 100), gmx::InvalidInputError);
}

TEST(BlockAveragerTest, AverageAndError)
{
    // Eight blocks of two values each, with the block averages 1 to 8
    BlockAverager ba(16, 8);
    for(int b = 1; b <= 8; b++)
    {
        ba.add(b - 0.5);
        ba.add(b + 0.5);
    }
    const auto ae = ba.result();
    EXPECT_EQ(16U, ae.n);
    EXPECT_NEAR(4.5, ae.average, 1e-14);
    // The variance of 1 to 8 is 6, the standard error is its root over the eight blocks
    EXPECT_NEAR(std::sqrt(6.0/8.0), ae.error, 1e-14);
}

TEST(BlockAveragerTest, CorrelatedValuesHaveALargerError)
{
    // The same values in each block: no error. The values of a block
    // equal: the error of the block averages.
    BlockAverager same(800), blocks(800);
    for(int i = 0; i < 800; i++)
    {
        same.add(i % 100);
        blocks.add(i/100);
    }
    EXPECT_NEAR(0.0, same.result().error, 1e-12);
    EXPECT_NEAR(std::sqrt(6.0/8.0), blocks.result().error, 1e-12);
    EXPECT_NEAR(49.5, same.result().average, 1e-12);
}

TEST(BlockAveragerTest, IncompleteSeries)
{
    // Three of eight blocks are complete, the fourth is not used for the error
    BlockAverager ba(80, 8);
    for(int i = 0; i < 35; i++)
    {
        ba.add(i < 10 ? 1.0 : (i < 20 ? 2.0 : (i < 30 ? 3.0 : 100.0)));
    }
    const auto ae = ba.result();
    EXPECT_EQ(35U, ae.n);
    EXPECT_NEAR((10 + 20 + 30 + 500)/35.0, ae.average, 1e-12);
    EXPECT_NEAR(std::sqrt(1.0/3.0), ae.error, 1e-12);
    // Without two complete blocks there is no error, and without values no average
    BlockAverager few(80, 8);
    few.add(3.0);
    EXPECT_EQ(0.0, few.result().error);
    EXPECT_EQ(3.0, few.result().average);
    EXPECT_EQ(0U, BlockAverager(10).result().n);
    // More values than expected go to the last block
    BlockAverager more(2, 2);
    more.add(1.0);
    more.add(2.0);
    more.add(4.0);
    EXPECT_EQ(3U, more.result().n);
    EXPECT_NEAR(7.0/3.0, more.result().average, 1e-14);
    EXPECT_NEAR(1.0, more.result().error, 1e-14);
}

TEST(BindingEnergyTest, EnergyAndEnthalpy)
{
    PathIntegralResults cluster;
    cluster.temperature  = 100;
    cluster.etot.average = -30;
    cluster.etot.error   = 0.3;
    std::vector<PathIntegralResults> alone(2);
    alone[0].etot.average = 5;
    alone[0].etot.error   = 0.4;
    alone[1].etot.average = -2;
    alone[1].etot.error   = 1.2;
    AverageError deltaH;
    const auto   deltaU = bindingEnergy(cluster, alone, &deltaH);
    EXPECT_NEAR(-33.0, deltaU.average, 1e-12);
    EXPECT_NEAR(1.3, deltaU.error, 1e-12);
    // One compound binds: minus k_B T
    EXPECT_NEAR(-33.0 - 0.0083144626*100, deltaH.average, 1e-6);
    EXPECT_NEAR(1.3, deltaH.error, 1e-12);
    // A third compound takes another k_B T
    alone.push_back(alone[0]);
    alone[2].etot.error = 0;
    (void) bindingEnergy(cluster, alone, &deltaH);
    EXPECT_NEAR(-38.0 - 2*0.0083144626*100, deltaH.average, 1e-6);
}


TEST(CovarianceAveragerTest, CovarianceAndError)
{
    // b = 2a + noise that does not correlate with a
    CovarianceAverager cov(8000);
    std::mt19937       gen(77);
    std::normal_distribution<double> normal(0.0, 1.0);
    for(int i = 0; i < 8000; i++)
    {
        double a = normal(gen);
        cov.add(a + 10, 2*a + 3*normal(gen) - 5);
    }
    const auto res = cov.result();
    EXPECT_EQ(8000U, res.n);
    // The variance of a is one
    EXPECT_NEAR(2.0, res.average, 0.15);
    EXPECT_GT(res.error, 0.01);
    EXPECT_LT(res.error, 0.15);
    // Without values there is nothing, with one block no error
    EXPECT_EQ(0U, CovarianceAverager(10).result().n);
    CovarianceAverager one(4, 1);
    for(int i = 0; i < 4; i++)
    {
        one.add(i, i);
    }
    EXPECT_NEAR(1.25, one.result().average, 1e-12);
    EXPECT_EQ(0.0, one.result().error);
}

TEST(GaussLegendreTest, PointsAndWeights)
{
    std::vector<double> x, w;
    gaussLegendre(1, &x, &w);
    ASSERT_EQ(1U, x.size());
    EXPECT_NEAR(0.5, x[0], 1e-14);
    EXPECT_NEAR(1.0, w[0], 1e-14);
    gaussLegendre(2, &x, &w);
    ASSERT_EQ(2U, x.size());
    EXPECT_NEAR(0.5 - 0.5/std::sqrt(3.0), x[0], 1e-14);
    EXPECT_NEAR(0.5 + 0.5/std::sqrt(3.0), x[1], 1e-14);
    EXPECT_NEAR(0.5, w[0], 1e-14);
    EXPECT_NEAR(0.5, w[1], 1e-14);
    for(int n = 3; n <= 16; n++)
    {
        gaussLegendre(n, &x, &w);
        ASSERT_EQ(static_cast<size_t>(n), x.size());
        // Polynomials up to the degree 2n-1 are integrated exactly
        for(int degree : { 0, 1, n, 2*n - 1 })
        {
            double sum = 0;
            for(int i = 0; i < n; i++)
            {
                EXPECT_GT(x[i], 0.0);
                EXPECT_LT(x[i], 1.0);
                if (i > 0)
                {
                    EXPECT_GT(x[i], x[i - 1]);
                }
                sum += w[i]*std::pow(x[i], degree);
            }
            EXPECT_NEAR(1.0/(degree + 1), sum, 1e-13) << n << " points, degree " << degree;
        }
    }
    EXPECT_THROW(gaussLegendre(0, &x, &w), gmx::InvalidInputError);
}

//! Frequencies in 1/ps of a bend and two stretches as those of water, and one that is not a vibration
const std::vector<double> waterLike = { 48.0, 110.0, 113.0 };

TEST(HarmonicBeadEnergyTest, ClassicalAndQuantumLimits)
{
    const double T  = 298.15;
    const double kT = BOLTZ*T;
    double       cv = 0;
    // One bead is the classical oscillator
    EXPECT_NEAR(3*kT, harmonicBeadEnergy(waterLike, T, 1, 1.0, &cv), 1e-12);
    EXPECT_NEAR(3*BOLTZ, cv, 1e-14);
    // Frequencies that are not positive are left out
    EXPECT_NEAR(3*kT, harmonicBeadEnergy({ 48.0, 0.0, 110.0, -5.0, 113.0 }, T, 1, 1.0), 1e-12);
    // Infinitely many beads: the quantum oscillator
    double energy = 0, cvq = 0;
    for(auto nu : waterLike)
    {
        const double x = 0.5*PLANCK*nu/kT;
        energy += 0.5*PLANCK*nu/std::tanh(x);
        cvq    += BOLTZ*(x/std::sinh(x))*(x/std::sinh(x));
    }
    EXPECT_NEAR(energy, harmonicBeadEnergy(waterLike, T, 0, 1.0, &cv), 1e-10);
    EXPECT_NEAR(cvq, cv, 1e-12);
    // which many beads approach from below, the heat capacity from above
    double previous = 3*kT, cvPrevious = 0;
    for(int nbeads : { 8, 32, 128, 2048 })
    {
        const double e = harmonicBeadEnergy(waterLike, T, nbeads, 1.0, &cv);
        EXPECT_GT(e, previous);
        EXPECT_LT(e, energy);
        if (nbeads > 8)
        {
            EXPECT_LT(cv, cvPrevious);
        }
        EXPECT_GT(cv, cvq);
        previous   = e;
        cvPrevious = cv;
    }
    EXPECT_NEAR(energy, previous, 1e-3);
    // With 32 beads water lacks close to 2 kJ/mol and has 11 J/mol K too much
    EXPECT_NEAR(1.76, energy - harmonicBeadEnergy(waterLike, T, 32, 1.0, &cv), 0.1);
    EXPECT_NEAR(-0.0112, cvq - cv, 0.001);
    // A small constant of Planck gives the classical oscillator, also the series for it
    EXPECT_NEAR(3*kT, harmonicBeadEnergy(waterLike, T, 0, 1e-6, &cv), 1e-9);
    EXPECT_NEAR(3*BOLTZ, cv, 1e-12);
    EXPECT_NEAR(3*kT, harmonicBeadEnergy(waterLike, T, 16, 1e-6), 1e-9);
    // Scaling the constant of Planck is the same as scaling the frequencies
    EXPECT_NEAR(harmonicBeadEnergy({ 24.0, 55.0, 56.5 }, T, 16, 1.0),
                harmonicBeadEnergy(waterLike, T, 16, 0.5), 1e-12);
    // The heat capacity is the derivative of the energy to the temperature
    for(int nbeads : { 0, 1, 16 })
    {
        const double dT = 0.01;
        const double numeric = (harmonicBeadEnergy(waterLike, T + dT, nbeads, 1.0) -
                                harmonicBeadEnergy(waterLike, T - dT, nbeads, 1.0))/(2*dT);
        harmonicBeadEnergy(waterLike, T, nbeads, 1.0, &cv);
        EXPECT_NEAR(numeric, cv, 1e-8) << nbeads;
    }
}

TEST(HarmonicBeadEnergyTest, FreeEnergyIsTheIntegralOfTheKineticEnergy)
{
    // F_q - F_cl = int_0^1 2 (K(g) - K_cl)/g dg, here with many points
    // since nothing is taken off the integrand
    for(double T : { 100.0, 298.15, 1000.0 })
    {
        const double        kT = BOLTZ*T;
        std::vector<double> x, w;
        gaussLegendre(40, &x, &w);
        double integral = 0;
        for(size_t i = 0; i < x.size(); i++)
        {
            integral += w[i]*2*(0.5*harmonicBeadEnergy(waterLike, T, 0, x[i]) - 1.5*kT)/x[i];
        }
        double expected = 0;
        for(auto nu : waterLike)
        {
            const double xx = 0.5*PLANCK*nu/kT;
            expected += kT*std::log(std::sinh(xx)/xx);
        }
        EXPECT_NEAR(expected, harmonicQuantumFreeEnergy(waterLike, T), 1e-10) << T;
        EXPECT_NEAR(expected, integral, 1e-6*expected) << T;
    }
    // Low frequencies, for which the series is used
    const double kT = BOLTZ*300;
    const double xx = 0.5*PLANCK*0.01/kT;
    EXPECT_NEAR(kT*xx*xx*(1.0/6 - xx*xx/180), harmonicQuantumFreeEnergy({ 0.01 }, 300), 1e-18);
    EXPECT_NEAR(kT*std::log(std::sinh(xx)/xx), harmonicQuantumFreeEnergy({ 0.01 }, 300), 1e-14);
    EXPECT_EQ(0.0, harmonicQuantumFreeEnergy({ 0.0, -3.0 }, 300));
}

/*! \brief A compound of three atoms that is harmonic, for computeThermochemistry
 * \param[in] T The temperature
 * \return the reference with the quantum values of harmonic oscillators
 */
HarmonicReference harmonicCompound(double T)
{
    HarmonicReference ref;
    ref.epot               = -123.4;
    ref.frequencies        = waterLike;
    ref.natom              = 3;
    ref.linear             = false;
    ref.symmetryNumber     = 2;
    ref.entropyTranslation = 144.8;
    ref.entropyRotation    = 43.7;
    const double kT        = BOLTZ*T;
    double       cv        = 0;
    ref.internalEnergy     = harmonicBeadEnergy(waterLike, T, 0, 1.0, &cv) + 3*kT;
    ref.heatCapacity       = KILO*(cv + 3*BOLTZ);
    ref.entropy            = ref.entropyTranslation + ref.entropyRotation;
    for(auto nu : waterLike)
    {
        const double x       = 0.5*PLANCK*nu/kT;
        ref.zeroPointEnergy += 0.5*PLANCK*nu;
        ref.entropy         += RGAS*(x/std::tanh(x) - std::log(2*std::sinh(x)));
    }
    return ref;
}

/*! \brief What simulations of a harmonic compound give, without noise
 * \param[in]  ref    The compound
 * \param[in]  T      The temperature
 * \param[in]  nbeads The number of beads
 * \param[in]  nquad  The number of points of the integrals
 * \param[out] ti     The simulations for the entropy
 * \return the simulation at the temperature
 */
PathIntegralResults harmonicSimulations(const HarmonicReference &ref,
                                        double                   T,
                                        int                      nbeads,
                                        int                      nquad,
                                        ThermoIntegration       *ti)
{
    const double        kT = BOLTZ*T;
    PathIntegralResults main;
    main.nbeads      = nbeads;
    main.temperature = T;
    double cv        = 0;
    // The vibrations with the beads, and the classical rotation and translation
    main.etot.average         = ref.epot + harmonicBeadEnergy(ref.frequencies, T, nbeads, 1.0, &cv) + 3*kT;
    main.etot.n               = 1000;
    main.heatCapacity.average = cv + 3*BOLTZ;
    gaussLegendre(nquad, &ti->nodes, &ti->weights);
    for(int i = 0; i < nquad; i++)
    {
        const double g = ti->nodes[i];
        AverageError epot, ekin;
        epot.average = ref.epot + 1.5*kT*g;
        ti->classicalEpot.push_back(epot);
        const int nb = std::max(1, static_cast<int>(g*nbeads));
        ti->quantumBeads.push_back(nb);
        ekin.average = 0.5*harmonicBeadEnergy(ref.frequencies, T, nb, g) + 3*kT;
        ti->quantumEkin.push_back(ekin);
    }
    return main;
}

TEST(ComputeThermochemistryTest, HarmonicCompound)
{
    // For a harmonic compound the simulations give the values of the
    // normal mode analysis back, whatever the number of beads
    for(double T : { 150.0, 298.15 })
    {
        for(int nbeads : { 4, 32 })
        {
            const auto        ref = harmonicCompound(T);
            ThermoIntegration ti;
            const auto        main   = harmonicSimulations(ref, T, nbeads, 4, &ti);
            const auto        thermo = computeThermochemistry(ref, main, &ti);
            EXPECT_EQ(nbeads, thermo.nbeads);
            EXPECT_NEAR(T, thermo.temperature, 1e-12);
            EXPECT_GT(thermo.beadEnergy, 0);
            EXPECT_LT(thermo.beadHeatCapacity, 0);
            EXPECT_NEAR(ref.internalEnergy, thermo.internalEnergy.average, 1e-9) << T << " " << nbeads;
            EXPECT_NEAR(ref.internalEnergy + BOLTZ*T, thermo.enthalpy.average, 1e-9);
            EXPECT_NEAR(ref.heatCapacity, thermo.heatCapacity.average, 1e-9);
            ASSERT_TRUE(thermo.haveEntropy);
            EXPECT_NEAR(0.0, thermo.anharmonic.average, 1e-10);
            EXPECT_NEAR(thermo.quantumHarmonic, thermo.quantum.average, 1e-9);
            EXPECT_NEAR(ref.entropy, thermo.entropy.average, 1e-7) << T << " " << nbeads;
            EXPECT_NEAR(ref.internalEnergy + BOLTZ*T - T*ref.entropy/KILO, thermo.gibbsEnergy.average, 1e-8);
            EXPECT_EQ(0.0, thermo.entropy.error);
        }
    }
}

TEST(ComputeThermochemistryTest, AnharmonicityAndErrors)
{
    const double      T   = 298.15;
    const double      kT  = BOLTZ*T;
    const auto        ref = harmonicCompound(T);
    ThermoIntegration ti;
    auto              main = harmonicSimulations(ref, T, 32, 4, &ti);
    const auto        harmonic = computeThermochemistry(ref, main, &ti);
    // A classical energy that is a (kT)^2 higher gives -a (kT)^2 in the
    // free energy and 2 a k^2 T in the entropy
    const double a = 0.002;
    for(size_t i = 0; i < ti.nodes.size(); i++)
    {
        ti.classicalEpot[i].average += a*kT*kT*ti.nodes[i]*ti.nodes[i];
        ti.classicalEpot[i].error    = 0.01;
        ti.quantumEkin[i].error      = 0.02;
    }
    main.etot.average      += a*kT*kT;
    main.etot.error         = 0.05;
    main.heatCapacity.error = 0.001;
    const auto thermo = computeThermochemistry(ref, main, &ti);
    EXPECT_NEAR(-a*kT*kT, thermo.anharmonic.average, 1e-10);
    EXPECT_NEAR(harmonic.entropy.average + KILO*2*a*BOLTZ*kT, thermo.entropy.average, 1e-7);
    EXPECT_NEAR(harmonic.internalEnergy.average + a*kT*kT, thermo.internalEnergy.average, 1e-10);
    EXPECT_NEAR(0.05, thermo.internalEnergy.error, 1e-12);
    EXPECT_NEAR(1.0, thermo.heatCapacity.error, 1e-12);
    // The errors of the integrals, with their weights
    double va = 0, vq = 0;
    for(size_t i = 0; i < ti.nodes.size(); i++)
    {
        va += std::pow(ti.weights[i]*0.01/(ti.nodes[i]*ti.nodes[i]), 2);
        vq += std::pow(2*ti.weights[i]*0.02/ti.nodes[i], 2);
    }
    EXPECT_NEAR(std::sqrt(va), thermo.anharmonic.error, 1e-12);
    EXPECT_NEAR(std::sqrt(vq), thermo.quantum.error, 1e-12);
    EXPECT_NEAR(KILO*std::sqrt(0.05*0.05 + va + vq)/T, thermo.entropy.error, 1e-9);
    EXPECT_NEAR(std::sqrt(va + vq), thermo.gibbsEnergy.error, 1e-12);
}

TEST(ComputeThermochemistryTest, OneBeadIsClassical)
{
    const double      T   = 298.15;
    const double      kT  = BOLTZ*T;
    const auto        ref = harmonicCompound(T);
    ThermoIntegration ti;
    const auto        main = harmonicSimulations(ref, T, 1, 3, &ti);
    // No simulations with a scaled constant of Planck are needed
    ti.quantumEkin.clear();
    ti.quantumBeads.clear();
    const auto thermo = computeThermochemistry(ref, main, &ti);
    EXPECT_EQ(0.0, thermo.beadEnergy);
    EXPECT_EQ(0.0, thermo.quantum.average);
    EXPECT_NEAR(6*kT, thermo.internalEnergy.average, 1e-10);
    EXPECT_NEAR(ref.classicalEnergy(T), thermo.internalEnergy.average, 1e-10);
    EXPECT_NEAR(KILO*6*BOLTZ, thermo.heatCapacity.average, 1e-9);
    ASSERT_TRUE(thermo.haveEntropy);
    // The entropy of classical oscillators, which is negative for these frequencies
    EXPECT_NEAR(ref.classicalEntropy(T), thermo.entropy.average, 1e-9);
    double vib = 0;
    for(auto nu : waterLike)
    {
        vib += RGAS*(1 - std::log(PLANCK*nu/kT));
    }
    EXPECT_NEAR(ref.entropyTranslation + ref.entropyRotation + vib, thermo.entropy.average, 1e-9);
}

TEST(ComputeThermochemistryTest, WithoutEntropyAndWrongInput)
{
    const double      T   = 298.15;
    const auto        ref = harmonicCompound(T);
    ThermoIntegration ti;
    const auto        main = harmonicSimulations(ref, T, 16, 4, &ti);
    const auto        thermo = computeThermochemistry(ref, main, nullptr);
    EXPECT_FALSE(thermo.haveEntropy);
    EXPECT_NEAR(ref.internalEnergy, thermo.internalEnergy.average, 1e-9);
    // Simulations that are missing
    auto fewer = ti;
    fewer.quantumEkin.pop_back();
    EXPECT_THROW(computeThermochemistry(ref, main, &fewer), gmx::InvalidInputError);
    fewer = ti;
    fewer.classicalEpot.clear();
    EXPECT_THROW(computeThermochemistry(ref, main, &fewer), gmx::InvalidInputError);
    fewer = ti;
    fewer.quantumBeads.clear();
    EXPECT_THROW(computeThermochemistry(ref, main, &fewer), gmx::InvalidInputError);
}

TEST(HarmonicReferenceTest, DegreesOfFreedom)
{
    HarmonicReference ref;
    ref.natom = 1;
    EXPECT_EQ(0, ref.nrotation());
    EXPECT_EQ(0, ref.nvibration());
    EXPECT_NEAR(1.5*BOLTZ*300, ref.classicalEnergy(300), 1e-12);
    ref.natom  = 2;
    ref.linear = true;
    EXPECT_EQ(2, ref.nrotation());
    EXPECT_EQ(1, ref.nvibration());
    EXPECT_NEAR(3.5*BOLTZ*300, ref.classicalEnergy(300), 1e-12);
    ref.natom  = 3;
    EXPECT_EQ(4, ref.nvibration());
    ref.linear = false;
    EXPECT_EQ(3, ref.nrotation());
    EXPECT_EQ(3, ref.nvibration());
    EXPECT_NEAR(6*BOLTZ*300, ref.classicalEnergy(300), 1e-12);
}

} // namespace

} // namespace alexandria
