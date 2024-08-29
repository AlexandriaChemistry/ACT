/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2023
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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "genetic_algorithm.h"

#include <cstdio>

#include "crossover.h"
#include "fitness_computer.h"
#include "initializer.h"
#include "mutator.h"
#include "probability_computer.h"
#include "terminator.h"

#include "act/basics/dataset.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace ga
{

Terminator *GeneticAlgorithm::terminator(const int index)
{
    GMX_RELEASE_ASSERT(index < static_cast<int>(terminators_->size()), "Index of terminator is out of bounds!");
    return terminators_->at(index);
}

void GeneticAlgorithm::updateGenePool(const GenePool &gpin)
{
}

bool GeneticAlgorithm::terminate(const GenePool *pool,
                                 const int       generationNumber)
{
    GMX_RELEASE_ASSERT(
        terminators_ != nullptr,
        "GA does not have a pointer to a vector of Terminator."
    );
    // Even if a terminator says we must halt, we will ask them all just in case
    // there is more than one stopping message to be printed
    bool halt = false;
    for (auto *term : *terminators_)
    {
        if (term->terminate(pool, generationNumber))
        {
            halt = true;
        }
    }
    return halt;
}

bool GeneticAlgorithm::penalize(      GenePool *pool,
                                const int       generationNumber)
{
    GMX_RELEASE_ASSERT(
        penalizers_ != nullptr,
        "GA does not have a pointer to a vector of Penalizer."
    );
    bool penalized = false;
    for (auto pen : *penalizers_)
    {
        if (pen->penalize(pool, generationNumber))
        {
            penalized = true;
        }
    }
    return penalized;
}

void GeneticAlgorithm::openFitnessFiles(const std::string &filename,
                                        iMolSelect         ims)
{
    std::string defname("ga_fitness");
    if (!filename.empty())
    {
        auto rpos = filename.rfind(".");
        if (std::string::npos != rpos)
        {
            defname = filename.substr(0, rpos);
        }
    }
    auto imsNames = iMolSelectNames();
    std::string fn = defname + gmx::formatString("_%s.dat", imsNames[ims]);
    fileFitness_.insert({ims, gmx_fio_fopen(fn.c_str(), "w")});
    GMX_RELEASE_ASSERT(fileFitness_[ims] != NULL, "Could not open file");
}

void GeneticAlgorithm::closeFitnessFiles()
{
    for(const auto &ff : fileFitness_)
    {
        gmx_fio_fclose(ff.second);
    }
}

void GeneticAlgorithm::fprintFitness(const GenePool &pool)
{
    for (const auto &ff : fileFitness_)
    {
        for (const auto &genome : pool.genePool())
        {
            if (genome.hasFitness(ff.first))
            {
                fprintf(ff.second, "%f ", genome.fitness(ff.first));
            }
        }
        fprintf(ff.second, "\n");
    }
}

}  //namespace ga
