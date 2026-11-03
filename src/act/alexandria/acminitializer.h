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

#include <memory>
#include <random>

#include "act/alexandria/acmindividual.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/ga/initializer.h"

namespace alexandria
{

/*!
 * \brief Initializes Individual instances as ACMIndividual objects
 */
class ACMInitializer : public ga::Initializer
{

private:

    //! Random device used to seed the random number generator
    std::random_device                       rd_;
    //! Mersenne twister random number generator
    std::mt19937                             gen_;
    //! Uniform distribution on [0, 1) used to draw random parameters
    std::uniform_real_distribution<double>   dis_;
    //! StaticIndividualInfo pointer
    StaticIndividualInfo                    *sii_;
    //! Whether we do random initialization or not.
    bool                                     randInit_;

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

    /*! \brief Create a new ACMIndividual and fill its genome
     * The parameters are drawn uniformly between the lower and upper
     * bounds if random initialization is requested, otherwise the
     * default parameter values are used. The best genome is set
     * to the initial genome.
     * \return pointer to a new individual
     */
    ACMIndividual *initialize() override;

    /*! \brief Set all bases of a genome to random values drawn
     * uniformly between the lower and upper bounds of each parameter
     * \param[inout] genome The Genome to randomize
     */
    void randomizeGenome(ga::Genome *genome) override;

};


} // namespace alexandria


#endif // ALEXANDRIA_ACMINITIALIZER_H
