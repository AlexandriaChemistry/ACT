/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
#include <math.h>

#include <map>

#include "mpi.h"

#include <gtest/gtest.h>

#include "../GeneticAlgorithm.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "alexandria/acm_ga.h"
#include "alexandria/bayes.h"
#include "alexandria/acmfitnesscomputer.h"
#include "alexandria/acminitializer.h"
#include "alexandria/acthelper.h"
#include "alexandria/actmiddleman.h"
#include "alexandria/confighandler.h"
#include "alexandria/mcmcmutator.h"
#include "alexandria/molgen.h"
#include "alexandria/molselect.h"
#include "alexandria/npointcrossover.h"
#include "alexandria/percentmutator.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"
#include "../FitnessComputer.h"
#include "../Sorter.h"
#include "../Terminator.h"

namespace ga
{

namespace
{

class GeneticAlgorithmTest : public gmx::test::CommandLineTestBase
{
    protected:
 
        GeneticAlgorithmTest ()
        {
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }
 
        static void TearDownTestCase()
        {
        }
    
        void testIt(alexandria::OptimizerAlg alg,
                    int nElites, int popSize,
                    bool verbose, int ncrossovers, 
                    const std::vector<std::string> &fitstrings,
                    int seed, std::vector<alexandria::eRMS> erms)
        {
            // GA stuff
            alexandria::GAConfigHandler      gach;
            gach.setPopSize(popSize);
            gach.setCrossovers(ncrossovers);
            gach.setOptimizerAlg(alg);
            int nmiddlemen = gach.popSize();
            alexandria::CommunicationRecord  cr;
            cr.init(nmiddlemen);
            gmx_output_env_t    *oenv;
            output_env_init_default(&oenv);
            
            // Create static individual
            alexandria::StaticIndividualInfo sii(&cr);
            std::string ffName("ACS-g.xml");
            std::string ffDataName = gmx::test::TestFileManager::getInputFilePath(ffName);
            sii.fillPoldata(nullptr, ffDataName.c_str());
            std::string selName("testAlcohol.dat");
            std::string selDataName = gmx::test::TestFileManager::getInputFilePath(selName);
            alexandria::MolSelect gms;
            gms.read(selDataName.c_str());
            
            // Now read the molprop file
            alexandria::MolGen               molgen(&cr);
            std::string mpName("testAlcohol.xml");
            std::string mpDataName = gmx::test::TestFileManager::getInputFilePath(mpName);
            for(const auto &fs : fitstrings)
            {
                molgen.addFitOption(fs);
            }
            // TODO check return value
            (void) molgen.Read(nullptr, mpDataName.c_str(), sii.poldata(),
                               gms, false);
            // Continue filling the shared individual
            sii.generateOptimizationIndex(nullptr, &molgen);
            sii.fillVectors(molgen.mindata());
            std::string xvgconv("param_conv.xvg"), xvgepot("param_epot.xvg");
            std::vector<std::string> paramClass;
            for(const auto &fm : molgen.typesToFit())
            {
                paramClass.push_back(fm.first);
            }
            sii.setOutputFiles(xvgconv.c_str(), paramClass, xvgepot.c_str());
            sii.assignParamClassIndex();
            for(const auto &er : erms)
            {
                sii.target(iMolSelect::Train, er)->setWeight(1.0);
            }
            sii.computeWeightedTemperature(true);
            sii.propagateWeightFittingTargets();
            // One more config handler
            alexandria::BayesConfigHandler bch;
            bch.setMaxIter(20);
            bch.setAnneal(0.5);
            bch.setTemperature(10);
            bch.setSeed(seed);
            bch.setStep(0.1);
            if (cr.isMaster())
            {
                // Checker, Master only
                gmx::test::TestReferenceChecker  checker_(this->rootChecker());
                auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-4);
                checker_.setDefaultTolerance(tolerance);
            
                checker_.checkInt64(nElites, "nElites");
                checker_.checkBoolean(verbose, "verbose");
                checker_.checkInt64(gach.popSize(), "popSize");
                checker_.checkInt64(gach.nCrossovers(), "ncrossovers");
                checker_.checkString(optimizerAlgToString(gach.optimizer()), "algorithm");
                for(auto &ttf : molgen.typesToFit())
                {
                    checker_.checkBoolean(ttf.second, ttf.first.c_str());
                }
                for(auto &io : molgen.iopt())
                {
                    checker_.checkBoolean(io.second, interactionTypeToString(io.first).c_str());
                }
                for(auto &ft : sii.fittingTargetsConst(iMolSelect::Train))
                {
                    double w = ft.second.weight();
                    if (w > 0)
                    {
                        std::string name = gmx::formatString("Weight_%s", rmsName(ft.first));
                        checker_.checkReal(w, name.c_str());
                    }
                }
                checker_.checkSequence(paramClass.begin(), paramClass.end(), "paramClass");
                std::vector<std::string> molnames;
                for(const auto &mm: molgen.mymols())
                {
                    molnames.push_back(mm.getMolname());
                }
                checker_.checkSequence(molnames.begin(), molnames.end(), "MoleculeName");
                checker_.checkSequence(sii.weightedTemperature().begin(),
                                   sii.weightedTemperature().end(), "weightedTemperatures");
                checker_.checkSequence(sii.paramNames().begin(), sii.paramNames().end(), "paramNames");
                // Now the rest of the classes
                std::string outputFile("GeneticAlgorithmTest.dat");
                // bool randInit     = false;
                
                // Adjust the seed that gets passed around to components of the optimizer
                // Create random number generator and feed it the global seed
                std::random_device rd;  // Will be used to obtain a seed for the random number engine
                std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
                std::uniform_int_distribution<int> dis(0); // Default constructor to cover all available (positive) range
                gen.seed(seed);
                seed = dis(gen);

                // Initializer
                auto *initializer = new alexandria::ACMInitializer(&sii, gach.randomInit(), seed);
                auto *fitComp     = new alexandria::ACMFitnessComputer(nullptr, &sii, &molgen, false, false);

                auto probComputer = new RankProbabilityComputer(gach.popSize());
                // Selector
                auto selector     = new ga::RouletteSelector(seed);
                auto crossover    = new ga::NPointCrossover(sii.nParam(),
                                                            gach.nCrossovers(),
                                                            seed);
            
                // Mutator
                Mutator *mut;
                if (alg == alexandria::OptimizerAlg::GA)
                {
                    mut = new alexandria::PercentMutator(&sii, seed, gach.percent());
                }
                else
                {
                    // mut = new alexandria::MCMCMutator(stdout, false, seed, &bch, fitComp, &sii, bch.evaluateTestset());
                    mut = new alexandria::MCMCMutator(nullptr, false, seed, &bch, fitComp, &sii, bch.evaluateTestset());
                }

                // Terminator
                std::vector<Terminator *> terminators;
                terminators.push_back(new ga::GenerationTerminator(gach.maxGenerations(), nullptr));
                GeneticAlgorithm *ga;
                if (alexandria::OptimizerAlg::MCMC == alg)
                {
                    checker_.checkInt64(bch.maxIter(), "bch.maxIter");
                    checker_.checkInt64(bch.seed(), "bch.seed");
                    checker_.checkReal(bch.temperature(), "bch.temperature");

                    ga = new ga::MCMC(stdout, initializer, fitComp, mut, &sii, &gach);
                }
                else
                {
                    if (alg == alexandria::OptimizerAlg::GA)
                    {
                        checker_.checkReal(gach.percent(), "gach.percent");
                    }
                    else if (alg == alexandria::OptimizerAlg::HYBRID)
                    {
                        checker_.checkReal(bch.step(), "bch.step");
                    }
                    ga = new ga::HybridGAMC(stdout, initializer, fitComp,
                                            probComputer, selector, crossover,
                                            mut, &terminators,
                                            &sii, &gach, bch.seed());
                }
                checker_.checkInt64(gach.maxGenerations(), "Maximum Number of Generations");
                checker_.checkReal(gach.prCross(), "Probability for Crossover");
                checker_.checkReal(gach.prMut(), "Probability for Mutation");
                Genome best;
                ga->evolve(&best);
                checker_.checkReal(best.fitness(iMolSelect::Train), "Training fitness");
                if (cr.nmiddlemen() > 1)  // If we have more middlemen than the master, stop them
                {
                    for(auto &dest : cr.middlemen())
                    {
                        if (cr.rank() != dest)
                        {
                            cr.send_done(dest);
                        }
                    }
                }

                auto imstr = iMolSelect::Train;
                checker_.checkReal(best.fitness(imstr), "Training fitness after evolve");
                
                auto bestDuplicate = best;
                fitComp->compute(&bestDuplicate, imstr);
                EXPECT_EQ(best.fitness(imstr), bestDuplicate.fitness(imstr));
                if (best.nBase() > 0)
                {
                    checker_.checkSequence(best.bases().begin(),
                                           best.bases().end(), "bestParam");
                }

                // Stop MASTER's helpers
                std::vector<double> dummy;
                fitComp->calcDeviation(&dummy, alexandria::CalcDev::Final, imstr);
            }
            else if (cr.isMiddleMan())
            {
                // Run middleman-like code.
                alexandria::ACTMiddleMan middleman(&molgen, 
                                                   &sii, &gach, &bch,
                                                   false, oenv);
                middleman.run();
            }
            else if (cr.isHelper())
            {
                // Run helper-like code
                alexandria::ACTHelper helper(&sii, &molgen);
                helper.run();
            }
        }
    
};


TEST_F (GeneticAlgorithmTest, PopSix)  // HYBRID
{
    GMX_MPI_TEST(6);
    testIt(alexandria::OptimizerAlg::HYBRID, 0, 6, true, 1,
           { "chi", "zeta" }, 1993, { alexandria::eRMS::QUAD });
}

TEST_F (GeneticAlgorithmTest, PopTwo)  // GA
{
    GMX_MPI_TEST(6);
    testIt(alexandria::OptimizerAlg::GA, 0, 2, true, 1,
           { "chi", "zeta" }, 1993, { alexandria::eRMS::QUAD });
}

TEST_F (GeneticAlgorithmTest, PopOneMCMC)  // MCMC
{
    GMX_MPI_TEST(6);
    testIt(alexandria::OptimizerAlg::MCMC, 0, 1, false, 0,
           { "chi", "jaa" }, 1993, { alexandria::eRMS::MU });
}

TEST_F (GeneticAlgorithmTest, PopTwoEspGA)  // GA
{
    GMX_MPI_TEST(6);
    testIt(alexandria::OptimizerAlg::GA, 0, 2, true, 1,
           { "chi", "jaa" }, 1993, { alexandria::eRMS::ESP });
}

TEST_F (GeneticAlgorithmTest, MCMCThreeVariables)  // MCMC
{
    GMX_MPI_TEST(6);
    testIt(alexandria::OptimizerAlg::MCMC, 0, 1, false, 0,
           { "chi", "jaa", "zeta" }, 1997, { alexandria::eRMS::MU });
}

TEST_F (GeneticAlgorithmTest, MCMCTwoVariablesQuad)  // MCMC
{
    GMX_MPI_TEST(6);
    testIt(alexandria::OptimizerAlg::MCMC, 0, 1, false, 0,
           { "jaa", "zeta" }, 1997, { alexandria::eRMS::QUAD });
}

TEST_F (GeneticAlgorithmTest, MCMCEpot)  // MCMC
{
    GMX_MPI_TEST(6);
    testIt(alexandria::OptimizerAlg::MCMC, 0, 1, false, 0,
           { "Dm", "sigma" }, 1991, { alexandria::eRMS::EPOT });
}

TEST_F (GeneticAlgorithmTest, GAEpot)  // GA
{
    GMX_MPI_TEST(6);
    testIt(alexandria::OptimizerAlg::GA, 0, 6, false, 0,
           { "epsilon", "kt" }, 1991, { alexandria::eRMS::EPOT } );
}

TEST_F (GeneticAlgorithmTest, GAEpotESP)  // GA
{
    GMX_MPI_TEST(6);
    testIt(alexandria::OptimizerAlg::GA, 0, 6, false, 0,
           { "epsilon", "kt", "jaa", "chi" }, 1991, { alexandria::eRMS::EPOT, alexandria::eRMS::ESP } );
}


}

} 
