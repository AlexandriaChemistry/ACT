/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */
#ifndef ALEXANDRIA_ACM_GA_H
#define ALEXANDRIA_ACM_GA_H

#include <cstdio>

#include <vector>

#include "act/ga/crossover.h"
#include "act/ga/fitness_computer.h"
#include "act/ga/genetic_algorithm.h"
#include "act/ga/initializer.h"
#include "act/ga/probability_computer.h"
#include "act/ga/selector.h"
#include "act/ga/terminator.h"

#include "confighandler.h"
#include "staticindividualinfo.h"

struct gmx_output_env;

namespace ga
{

class MsgHandler;

class HybridGAMC : public GeneticAlgorithm
{
private:
    //! Who am I?
    alexandria::StaticIndividualInfo *sii_;
    //! GAConfigHandler pointer
    alexandria::GAConfigHandler      *gach_;
    //! Output filename for fitness files
    const char                       *fitnessFile_;
    //! Gene pool input (may be null pointer)
    const char                       *gpin_;
    //! Gene pool output
    const char                       *gpout_;
    //! seed for random numbers
    int                               seed_;
public:
    /*!
     * \brief Constructor for self-building
     */
    HybridGAMC(Initializer                         *initializer,
               FitnessComputer                     *fitnessComputer,
               ProbabilityComputer                 *probComputer,
               Selector                            *selector,
               Crossover                           *crossover,
               Mutator                             *mutator,
               std::vector<Terminator*>            *terminators,
               std::vector<Penalizer*>             *penalizers,
               alexandria::StaticIndividualInfo    *sii,
               alexandria::GAConfigHandler         *gach,
               const char                          *fitnessFileName,
               const char                          *genePoolIn,
               const char                          *genePoolOut,
               int                                  seed)
    : GeneticAlgorithm(initializer, fitnessComputer, probComputer, selector, crossover,
                       mutator, terminators, penalizers, gach->popSize()),
      sii_(sii), gach_(gach), fitnessFile_(fitnessFileName),
      gpin_(genePoolIn), gpout_(genePoolOut), seed_(seed)
    {}

    //! \copydoc ga::GeneticAlgorithm::evolve
    virtual bool evolve(alexandria::MsgHandler       *msghandler,
                        std::map<iMolSelect, Genome> *bestGenome);

};

class MCMC : public GeneticAlgorithm
{
private:
    //! Who am I?
    alexandria::StaticIndividualInfo *sii_;
    //! GAConfigHandler pointer
    alexandria::GAConfigHandler      *gach_;
public:
    /*!
     * \brief Constructor for self-building
     */
    MCMC(Initializer                         *initializer,
         FitnessComputer                     *fitnessComputer,
         Mutator                             *mutator,
         alexandria::StaticIndividualInfo    *sii,
         alexandria::GAConfigHandler         *gach)
    : GeneticAlgorithm(initializer, fitnessComputer, nullptr, nullptr, nullptr,
                       mutator, nullptr, nullptr, gach->popSize()),
      sii_(sii), gach_(gach)
    {}

    //! \copydoc ga::GeneticAlgorithm::evolve
    virtual bool evolve(alexandria::MsgHandler       *msghandler,
                        std::map<iMolSelect, Genome> *bestGenome);

};

} // namespace ga

#endif
