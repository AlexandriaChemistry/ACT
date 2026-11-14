/*
 * This file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

/*! \file
 * \brief
 * Tests for MCMCMutator class.
 *
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 * \ingroup module_alexandria
 */

#include "actpre.h"

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/mcmcmutator.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/acmindividual.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/ga/genome.h"
#include "act/ga/ga_test_helper.h"
#include "gromacs/utility/stringutil.h"
#include "act/basics/msg_handler.h"

#include "gromacs/utility/exceptions.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

// Test class that follows the project's testing patterns
class MCMCMutatorTest :  public gmx::test::CommandLineTestBase
{
protected:
    gmx::test::TestReferenceChecker *checker_;
    gmx_output_env_t     *oenv_;
    // Methods
    void SetUp() override
    {
        // Setup code for tests if needed
        checker_ = new gmx::test::TestReferenceChecker(this->rootChecker());
        const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_->setDefaultTolerance(tolerance);
    }

    void TearDown() override
    {
        // Cleanup code for tests if needed
    }

    void doMutate(ga::Genome                     &before,
                  int                             maxiter,
                  const std::vector<std::string> &fitstrings,
                  std::vector<alexandria::eRMS>  &eRms)
    {
        int                       maxgenerations = 10;
        ChargeGenerationAlgorithm alg            = ChargeGenerationAlgorithm::EEM;
        int nmiddlemen = 1;
        // The first seed is checked at the root level, the others in a compound
        for(int seed : { 1993, 1997 })
        {
            gmx::test::TestReferenceChecker checker(*checker_);
            if (seed != 1993)
            {
                checker = checker_->checkCompound("Seed", gmx::formatString("Seed%d", seed).c_str());
            }
            // Start from the same genome for each seed
            ga::Genome         mybefore = before;
            BayesConfigHandler bch;
            auto gth   = new GaTestHelper(nmiddlemen, fitstrings, eRms);
            auto mymut = new MCMCMutator(seed, &bch, gth->fitnessComputer,
                                         gth->sii, false, alg, maxgenerations);
            bch.setMaxIter(maxiter);
            checker.checkInt64(bch.maxIter(), "MaxIter");
            checker.checkReal(mybefore.fitness(iMolSelect::Train), "Fitness-before");
            checker.checkSequence(mybefore.bases().begin(),
                                  mybefore.bases().end(), "Genome-before");
            checker.checkSequence(gth->sii->paramNames().begin(),
                                  gth->sii->paramNames().end(), "Parameters");
            ga::Genome after = mybefore;
            mymut->mutate(gth->msghandler, &mybefore, &after, 0.0);
            checker.checkReal(after.fitness(iMolSelect::Train), "Fitness-after");
            checker.checkSequence(after.bases().begin(),
                                  after.bases().end(), "Genome-after");
            checker.checkInt64(mymut->numberObjectiveFunctionCalls(), "numberObjectiveFunctionCalls");
            checker.checkInt64(gth->forceComputer->numEvaluations(), "numEvaluations ForceComputer");
            checker.checkInt64(gth->forceComputer->numShellIterations(), "numShellIterations ForceComputer");
            checker.checkInt64(gth->forceComputer->numShellConvergenceFailed(), "numShellConvergenceFailed ForceComputer");
            auto allRand = mymut->allRand();
            checker.checkSequence(allRand.begin(), allRand.end(), "Random numbers");
            auto randInd = mymut->randIndices();
            checker.checkSequence(randInd.begin(), randInd.end(), "Random parameter indices");
            delete mymut;
            delete gth;
        }
   }
};

TEST_F(MCMCMutatorTest, MutateCharges10)
{
    // This is the core functionality of the MCMCMutator
    const std::vector<std::string> fitstrings = { "eta", "chi" };
    // Train on the dipole
    std::vector<alexandria::eRMS> eRms = { eRMS::MU };
    // Bases correspond to chi eta chi eta chi eta
    const std::vector<double> bases = {4.5, 12, 3.3, 10, 5, 13 };
    const FitnessMap    fm    = { {iMolSelect::Train, 10.0} };
    ga::Genome before(bases, fm);
    int maxiter = 2;
    doMutate(before, maxiter, fitstrings, eRms);
}

TEST_F(MCMCMutatorTest, MutateChargesESP50)
{
    // This is the core functionality of the MCMCMutator
    const std::vector<std::string> fitstrings = { "eta", "chi" };
    // Train on the electrostatic potential
    std::vector<alexandria::eRMS> eRms = { eRMS::ESP };
    // Bases correspond to chi eta chi eta chi eta
    const std::vector<double> bases = {4.5, 12, 3.3, 10, 5, 13 };
    const FitnessMap    fm    = { {iMolSelect::Train, 10000.0} };
    ga::Genome before(bases, fm);
    int maxiter = 5;
    doMutate(before, maxiter, fitstrings, eRms);
}

TEST_F(MCMCMutatorTest, MutateBonds10)
{
    // This is the core functionality of the MCMCMutator
    const std::vector<std::string> fitstrings = { "kb", "bondlength", "bondenergy" };
    // Train on internal energy
    std::vector<alexandria::eRMS> eRms = { eRMS::EPOT };
    // Bases correspond to bondenergy bondlength kb for each bond
    const std::vector<double> bases = { -250, 130, 3e5,  -200, 100, 3.2e5,  -300, 120, 3e5, -400, 100, 3.7e5 };
    const FitnessMap    fm    = { {iMolSelect::Train, 1e8} };
    ga::Genome before(bases, fm);
    int maxiter = 2;
    doMutate(before, maxiter, fitstrings, eRms);
}

TEST_F(MCMCMutatorTest, MutateBonds100)
{
    // This is the core functionality of the MCMCMutator
    const std::vector<std::string> fitstrings = { "kb", "bondlength", "bondenergy" };
    // Train on internal energy
    std::vector<alexandria::eRMS> eRms = { eRMS::EPOT };
    // Bases correspond to bondenergy bondlength kb for each bond
    const std::vector<double> bases = { -250, 130, 3e5,  -200, 100, 3.2e5,  -300, 120, 3e5, -400, 100, 3.7e5 };
    const FitnessMap    fm    = { {iMolSelect::Train, 1e8} };
    ga::Genome before(bases, fm);
    int maxiter = 10;
    doMutate(before, maxiter, fitstrings, eRms);
}

TEST_F(MCMCMutatorTest, MutateAngles)
{
    const std::vector<std::string> fitstrings = { "kt", "angle" };
    // Train on internal energy
    std::vector<alexandria::eRMS> eRms = { eRMS::EPOT };
    // Bases correspond to angle kt for each angle
    const std::vector<double> bases = { 105, 320, 111, 290, 103, 300, 120, 280, 107, 300, 108, 340 };
    const FitnessMap    fm    = { {iMolSelect::Train, 1e8} };
    ga::Genome before(bases, fm);
    int maxiter = 10;
    doMutate(before, maxiter, fitstrings, eRms);
}

TEST_F(MCMCMutatorTest, MutateDihedrals)
{
    const std::vector<std::string> fitstrings = { "kp", "phi0" };
    // Train on internal energy
    std::vector<alexandria::eRMS> eRms = { eRMS::EPOT };
    // Bases correspond to kp phi0 for each dihedral
    const std::vector<double> bases = { 6.5, 180, 5, 320, 6, 290, 7, 300, 10, 280, 8, 300, 4, 340 };
    const FitnessMap    fm    = { {iMolSelect::Train, 1e8} };
    ga::Genome before(bases, fm);
    int maxiter = 10;
    doMutate(before, maxiter, fitstrings, eRms);
}

TEST_F(MCMCMutatorTest, MutateLennardJones)
{
    // This is the core functionality of the MCMCMutator
    const std::vector<std::string> fitstrings = { "sigma", "epsilon" };
    // Train on internal energy
    std::vector<alexandria::eRMS> eRms = { eRMS::EPOT };
    // Bases correspond to epsilon sigma for each atom type
    const std::vector<double> bases = { 0.5, 0.45, 0.09, 0.4, 0.52, 0.3 };
    const FitnessMap    fm    = { {iMolSelect::Train, 1e8} };
    ga::Genome before(bases, fm);
    int maxiter = 10;
    doMutate(before, maxiter, fitstrings, eRms);
}

// Test: Each parameter is chosen with the same probability. The random
// numbers of a fine grid in [0, 1), away from the edges of the intervals of
// the indices, give each index equally often, and the ends give the first
// and the last index. The last index was chosen 1.5 times as often before.
TEST(MCMCMutatorIndexTest, IndexOfRandomUniform)
{
    for (size_t n : { 1, 2, 3, 7, 60 })
    {
        SCOPED_TRACE(gmx::formatString("%zu parameters", n));
        const size_t        perIndex = 1000;
        std::vector<size_t> count(n, 0);
        for (size_t i = 0; i < n*perIndex; i++)
        {
            double r = (i + 0.5)/(n*perIndex);
            size_t k = MCMCMutator::indexOfRandom(r, n);
            ASSERT_LT(k, n);
            count[k] += 1;
        }
        for (size_t k = 0; k < n; k++)
        {
            EXPECT_EQ(perIndex, count[k]) << "index " << k;
        }
        EXPECT_EQ(0u, MCMCMutator::indexOfRandom(0.0, n));
        EXPECT_EQ(n - 1, MCMCMutator::indexOfRandom(std::nextafter(1.0, 0.0), n));
        // Also if a random number generator gives 1
        EXPECT_EQ(n - 1, MCMCMutator::indexOfRandom(1.0, n));
    }
}

// Test: The random numbers of the MCMC steps only depend on the seed and the
// number of steps, three for each step, not on the objective function, such
// that they can be drawn before the steps are computed, for instance on a GPU
TEST_F(MCMCMutatorTest, RandomNumbersIndependentOfFitness)
{
    const std::vector<std::string> fitstrings = { "eta", "chi" };
    const std::vector<double>      bases      = { 4.5, 12, 3.3, 10, 5, 13 };
    const int                      maxiter    = 3;
    std::vector<std::vector<double>> numbers;
    std::vector<std::vector<int>>    indices;
    std::vector<double>              fitness;
    // The dipole and the electrostatic potential give other values of the
    // objective function, and therefore other steps that are accepted
    for (auto erms : { eRMS::MU, eRMS::ESP })
    {
        const std::vector<alexandria::eRMS> eRms = { erms };
        BayesConfigHandler bch;
        bch.setMaxIter(maxiter);
        auto gth   = new GaTestHelper(1, fitstrings, eRms);
        auto mymut = new MCMCMutator(1993, &bch, gth->fitnessComputer, gth->sii, false,
                                     ChargeGenerationAlgorithm::EEM, 10);
        ga::Genome genome(bases, { { iMolSelect::Train, 10.0 } });
        ga::Genome best = genome;
        mymut->mutate(gth->msghandler, &genome, &best, 0.0);
        size_t nsteps = maxiter*gth->sii->nParam();
        EXPECT_EQ(3*nsteps, mymut->allRand().size());
        EXPECT_EQ(nsteps, mymut->randIndices().size());
        numbers.push_back(mymut->allRand());
        indices.push_back(mymut->randIndices());
        fitness.push_back(best.fitness(iMolSelect::Train));
        delete mymut;
        delete gth;
    }
    EXPECT_NE(fitness[0], fitness[1]);
    EXPECT_EQ(numbers[0], numbers[1]);
    EXPECT_EQ(indices[0], indices[1]);
}

} // namespace alexandria
