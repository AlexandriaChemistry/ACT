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

#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/mcmcmutator.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/acmindividual.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/ga/genome.h"
#include "act/ga/ga_test_helper.h"
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
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
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
        BayesConfigHandler        bch;
        int                       seed           = 1993;
        int                       maxgenerations = 10;
        ChargeGenerationAlgorithm alg            = ChargeGenerationAlgorithm::EEM;
        int nmiddlemen = 1;
        auto gth   = new GaTestHelper(nmiddlemen, fitstrings, eRms);
        auto mymut = new MCMCMutator(seed, &bch, gth->fitnessComputer,
                                     gth->sii, false, alg, maxgenerations);
        bch.setMaxIter(maxiter);
        checker_->checkInt64(bch.maxIter(), "MaxIter");
        checker_->checkReal(before.fitness(iMolSelect::Train), "Fitness-before");
        checker_->checkSequence(before.bases().begin(),
                                before.bases().end(), "Genome-before");
        checker_->checkSequence(gth->sii->paramNames().begin(),
                                gth->sii->paramNames().end(), "Parameters");
        ga::Genome after = before;
        mymut->mutate(gth->msghandler, &before, &after, 0.0);
        checker_->checkReal(after.fitness(iMolSelect::Train), "Fitness-after");
        checker_->checkSequence(after.bases().begin(),
                                after.bases().end(), "Genome-after");
        checker_->checkInt64(mymut->numberObjectiveFunctionCalls(), "numberObjectiveFunctionCalls");
        checker_->checkInt64(gth->forceComputer->numEvaluations(), "numEvaluations ForceComputer");
        checker_->checkInt64(gth->forceComputer->numShellIterations(), "numShellIterations ForceComputer");
        checker_->checkInt64(gth->forceComputer->numShellConvergenceFailed(), "numShellConvergenceFailed ForceComputer");
        auto allRand = mymut->allRand();
        checker_->checkSequence(allRand.begin(), allRand.end(), "Random numbers");
        delete mymut;
        delete gth;
   }
};

TEST_F(MCMCMutatorTest, MutateCharges10)
{
    // This is the core functionality of the MCMCMutator
    const std::vector<std::string> fitstrings = { "eta", "chi" };
    // Train on the dipole
    std::vector<alexandria::eRMS> eRms = { eRMS::MU };
    // Bases correspond to chi eta chi eta chi eta
    std::vector<double> bases = {4.5, 12, 3.3, 10, 5, 13 };
    FitnessMap          fm    = { {iMolSelect::Train, 10.0} };
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
    std::vector<double> bases = {4.5, 12, 3.3, 10, 5, 13 };
    FitnessMap          fm    = { {iMolSelect::Train, 10000.0} };
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
    std::vector<double> bases = { -250, 130, 3e5,  -200, 100, 3.2e5,  -300, 120, 3e5, -400, 100, 3.7e5 };
    FitnessMap          fm    = { {iMolSelect::Train, 1e8} };
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
    std::vector<double> bases = { -250, 130, 3e5,  -200, 100, 3.2e5,  -300, 120, 3e5, -400, 100, 3.7e5 };
    FitnessMap          fm    = { {iMolSelect::Train, 1e8} };
    ga::Genome before(bases, fm);
    int maxiter = 10;
    doMutate(before, maxiter, fitstrings, eRms);
}

} // namespace alexandria
