/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */
#ifndef ALEXANDRIA_MCMC_H
#define ALEXANDRIA_MCMC_H

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
