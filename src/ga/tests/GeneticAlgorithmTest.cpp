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

#include <gtest/gtest.h>

#include "ga/GeneticAlgorithm.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "alexandria/acm_ga.h"
#include "alexandria/acmfitnesscomputer.h"
#include "alexandria/acminitializer.h"
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
#include "ga/FitnessComputer.h"
#include "ga/Sorter.h"
#include "ga/Terminator.h"

namespace ga
{

namespace
{

class GeneticAlgorithmTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
 
        GeneticAlgorithmTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }
 
        static void TearDownTestCase()
        {
        }

        void testIt(int nElites, int popSize, double tolerance,
                    bool verbose, int nrep, int ncrossovers, bool hybrid_gamc,
                    const std::vector<std::string> &fitstrings,
                    int seed)
        {
            checker_.checkInt64(nElites, "nElites");
            checker_.checkReal(tolerance, "tolerance");
            checker_.checkBoolean(verbose, "verbose");
            checker_.checkInt64(nrep, "nrep");
            alexandria::GAConfigHandler      gach;
            gach.setPopSize(popSize);
            gach.setCrossovers(ncrossovers);
            checker_.checkInt64(gach.popSize(), "popSize");
            checker_.checkInt64(gach.nCrossovers(), "ncrossovers");
            
            alexandria::CommunicationRecord  cr;
            cr.init(popSize);
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
            for(auto &ttf : molgen.typesToFit())
            {
                checker_.checkBoolean(ttf.second, ttf.first.c_str());
            }
            (void) molgen.Read(nullptr, mpDataName.c_str(), sii.poldata(), false, gms, nullptr, false);
            for(auto &io : molgen.iopt())
            {
                checker_.checkBoolean(io.second, interactionTypeToString(io.first).c_str());
            }
            std::vector<std::string> molnames;
            for(const auto &mm: molgen.mymols())
            {
                molnames.push_back(mm.getMolname());
            }
            checker_.checkSequence(molnames.begin(), molnames.end(), "MoleculeName");
            
            // Continue filling the shared individual
            sii.generateOptimizationIndex(nullptr, &molgen);
            sii.fillVectors(molgen.mindata());
            std::string xvgconv("param_conv.xvg"), xvgepot("param_epot.xvg");
            std::vector<std::string> paramClass;
            for(const auto &fm : molgen.typesToFit())
            {
                paramClass.push_back(fm.first);
            }
            checker_.checkSequence(paramClass.begin(), paramClass.end(), "paramClass");
            sii.setOutputFiles(xvgconv.c_str(), paramClass, xvgepot.c_str());
            sii.assignParamClassIndex();
            sii.target(iMolSelect::Train, alexandria::eRMS::MU)->setWeight(1.0);
            sii.computeWeightedTemperature(true);
            sii.propagateWeightFittingTargets();
            checker_.checkSequence(sii.weightedTemperature().begin(),
                                   sii.weightedTemperature().end(), "weightedTemperatures");
            checker_.checkSequence(sii.paramNames().begin(), sii.paramNames().end(), "paramNames");
	        
            alexandria::BayesConfigHandler bch;
            bch.setMaxIter(5);
            bch.setSeed(seed);
            // Now the rest of the classes
            std::string outputFile("GeneticAlgorithmTest.dat");
            bool randInit     = false;
            auto init         = new alexandria::ACMInitializer(&sii, randInit, bch.seed());
            auto fit          = new alexandria::ACMFitnessComputer(nullptr, &sii, &molgen, false, verbose, false);
            auto probComputer = new RankProbabilityComputer(gach.popSize());
            // Selector
            auto selector     = new ga::RouletteSelector();
            auto crossover    = new alexandria::NPointCrossover(gach.popSize(),
                                                                gach.nCrossovers());
            
            // Terminator
            auto terminator   = new ga::GenerationTerminator(gach.maxGenerations());
            GeneticAlgorithm *ga;
            if (hybrid_gamc)
            {
                auto mutator = new alexandria::PercentMutator(&sii, gach.percent());
                checker_.checkInt64(gach.percent(), "gach.percent");
                ga = new ga::HybridGAMC(nullptr, init, 
                                        fit, probComputer,
                                        selector, crossover, mutator, terminator,
                                        &sii, &gach);
            }
            else
            {
                checker_.checkInt64(bch.maxIter(), "bch.maxIter");
                checker_.checkInt64(bch.seed(), "bch.seed");
                checker_.checkReal(bch.temperature(), "bch.temperature");
                auto mutator      = new alexandria::MCMCMutator(nullptr, false, &bch, fit, &sii);

                ga = new ga::MCMC(nullptr, init, 
                                  fit, probComputer,
                                  selector, crossover, mutator, terminator,
                                  &sii, &gach, false);
            }
            checker_.checkInt64(gach.maxGenerations(), "Maximum Number of Generations");
            checker_.checkReal(gach.prCross(), "Probability for Crossover");
            checker_.checkReal(gach.prMut(), "Probability for Mutation");
            if (cr.isMaster())
            {
                Genome best;
                ga->evolve(&best);
                if (best.nBase() > 0)
                {
                    checker_.checkSequence(best.bases().begin(),
                                           best.bases().end(), "bestParam");
                }
            }
            else
            {
                // Run middleman like code.
            }
        }
    
};

TEST_F (GeneticAlgorithmTest, EmptyPop)
{
    GMX_MPI_TEST(3);
    testIt(0, 0, 0.1, false, 1, 1, true, { "sigma", "alpha" }, 1993);
}

TEST_F (GeneticAlgorithmTest, PopFour)
{
    GMX_MPI_TEST(3);
    testIt(0, 4, 0.1, false, 1, 1, true, { "epsilon", "gamma" }, 1993);
}

TEST_F (GeneticAlgorithmTest, PopOneMCMC)
{
    GMX_MPI_TEST(3);
    testIt(0, 1, 0.1, false, 1, 1, false, { "epsilon", "gamma" }, 1993);
}

}

} 
