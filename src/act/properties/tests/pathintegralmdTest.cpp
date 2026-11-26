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

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <random>
#include <vector>

#include "act/properties/fragmentwall.h"
#include "act/properties/pathintegralmd.h"
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

} // namespace

} // namespace alexandria
