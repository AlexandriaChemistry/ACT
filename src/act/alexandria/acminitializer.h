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


#ifndef ALEXANDRIA_ACMINITIALIZER_H
#define ALEXANDRIA_ACMINITIALIZER_H

#include "act/ga//initializer.h"
#include "staticindividualinfo.h"

#include <time.h>
#include <random>

namespace alexandria
{

/*!
 * \brief Initializes Individual instances as ACMIndividual objects
 */
class ACMInitializer : public ga::Initializer
{

private:

    // Random number generation
    std::random_device                       rd_;
    std::mt19937                             gen_;
    std::uniform_real_distribution<double>   dis_;
    //! StaticIndividualInfo pointer
    StaticIndividualInfo                    *sii_;
    //! Whether we do random initialization or not.
    bool                                     randInit_;
    //! Seeds for random number generation
    // std::vector<int>                         seeds_;

public: 

    /*!
     * \brief Property constructor
     * \param[in] sii           pointer to StaticIndividualInfo instance
     * \param[in] randInit      whether we initialize the force field parameters randomly
     * \param[in] seed          Seed for random number initialization,
     *                          if zero it will be generated.
     */
    ACMInitializer(StaticIndividualInfo   *sii,
                   bool                    randInit,
                   int                     seed);

    ga::Individual *initialize() override;

    void randomizeGenome(ga::Genome *genome) override;

};


} // namespace alexandria


#endif // ALEXANDRIA_ACMINITIALIZER_H
