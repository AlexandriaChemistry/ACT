/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */


#include "acminitializer.h"

#include "act/ga/individual.h"
#include "acmindividual.h"


namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMInitializer                *
* * * * * * * * * * * * * * * * * * * */

ACMInitializer::ACMInitializer(StaticIndividualInfo *sii,
                               bool                  randInit,
                               int                   seed)
    : gen_(rd_()), dis_(std::uniform_real_distribution<double>(0.0, 1.0)),
      sii_(sii), randInit_(randInit)
{
    gen_.seed(seed);
}

ACMIndividual *ACMInitializer::initialize()
{
    int id   = sii_->commRec()->middleManOrdinal();
    auto ind = new ACMIndividual(id, sii_);
    if (randInit_)
    // Insert random value in range
    {
        for (size_t i = 0; i < sii_->nParam(); i++)
        {
            ind->addParam(dis_(gen_)*(sii_->upperBoundAtIndex(i) - sii_->lowerBoundAtIndex(i))
                          + sii_->lowerBoundAtIndex(i));
        }
    }
    else  // Insert default values in StaticIndividualInfo
    {
        for (const double val : sii_->defaultParam())
        {
            ind->addParam(val);
        }
    }
    ind->setBestGenome(ind->genome());
    return ind;
}

void ACMInitializer::randomizeGenome(ga::Genome *genome)
{
    for (size_t i = 0; i < sii_->nParam(); i++)
    {
        genome->setBase(
            i,
            dis_(gen_)*(sii_->upperBoundAtIndex(i) - sii_->lowerBoundAtIndex(i)) + sii_->lowerBoundAtIndex(i)
        );
    }
}

/* * * * * * * * * * * * * * * * * * * *
* END: ACMInitializer                  *
* * * * * * * * * * * * * * * * * * * */


} // namespace alexandria
