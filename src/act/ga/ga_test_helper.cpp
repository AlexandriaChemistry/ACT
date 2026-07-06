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
 * Implements testing of the genetic algorithm
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "config.h"

#include "ga_test_helper.h"

#include "act/alexandria/acm_ga.h"
#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/acminitializer.h"
#include "act/alexandria/mcmcmutator.h"
#include "act/alexandria/molselect.h"
#include "act/alexandria/percentmutator.h"
#include "act/basics/msg_handler.h"
#include "act/forces/forcecomputer.h"
#include "act/ga/fitness_computer.h"
#include "act/ga/genetic_algorithm.h"
#include "act/ga/sorter.h"
#include "act/ga/terminator.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

GaTestHelper::GaTestHelper(int                                  nmiddlemen,
                           const std::vector<std::string>      &fitstrings,
                           const std::vector<alexandria::eRMS> &erms)
{
    cr.init(nmiddlemen);
    msghandler = new MsgHandler;
    msghandler->setPrintLevel(alexandria::ACTStatus::Warning);

    sii = new StaticIndividualInfo(&cr);
    std::string ffName("ACS-g.xml");
    std::string ffDataName = gmx::test::TestFileManager::getInputFilePath(ffName);
    sii->fillForceField(msghandler, ffDataName.c_str());
    molgen = new MolGen(&cr);
    std::string mpName("../../ga/tests/testAlcohol.xml");
    std::string mpDataName = gmx::test::TestFileManager::getInputFilePath(mpName);
    for(const auto &fs : fitstrings)
    {
        molgen->addFitOption(fs);
    }
    // Read molprops
    std::vector<const char *>  desc;
    std::vector<t_pargs>       pa;
    std::vector<t_filenm>      fnms;
    compR.addOptions(&pa, &fnms, &desc);
    // We have an old molprop file
    compR.setOneH(true);
    // Add molgen options
    molgen->addFilenames(&fnms);
    // Hack structure to add a file name
    fnms.back().filenames.push_back(mpDataName);
    // Selection file
    std::string selName("../../ga/tests/testAlcohol.dat");
    std::string selDataName = gmx::test::TestFileManager::getInputFilePath(selName);
    alexandria::MolSelect gms;
    // ForceComputer
    forceComputer   = new ForceComputer();
    gms.read(selDataName.c_str());
    //! \todo  check return value
    (void) molgen->Read(msghandler, fnms, sii->forcefield(),
                        gms, sii->fittingTargetsConst(iMolSelect::Train),
                        &compR, forceComputer);
    // Continue filling the shared individual
    sii->generateOptimizationIndex(msghandler, molgen, sii->commRec());
    sii->fillVectors(molgen->mindata());
    std::string xvgconv("param_conv.xvg"), xvgepot("param_epot.xvg");
    std::vector<std::string> paramClass;
    for(const auto &fm : molgen->typesToFit())
    {
        paramClass.push_back(fm.first);
    }
    sii->setOutputFiles(xvgconv.c_str(), paramClass, xvgepot.c_str());
    sii->assignParamClassIndex();
    for(const auto &er : erms)
    {
        sii->target(iMolSelect::Train, er)->setWeight(1.0);
    }
    sii->computeWeightedTemperature(true);
    sii->propagateWeightFittingTargets();
    // Random number stuff
    // Adjust the seed that gets passed around to components of the optimizer
    // Create random number generator and feed it the global seed
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> dis(0); // Default constructor to cover all available (positive) range
    int seed = 1993;
    
    // Initializer
    initializer     = new ACMInitializer(sii, false, seed);
    fitnessComputer = new ACMFitnessComputer();
    fitnessComputer->init(msghandler, sii, molgen, false, forceComputer,
                          sii->forcefield()->chargeGenerationAlgorithm());
}

} // namespace
