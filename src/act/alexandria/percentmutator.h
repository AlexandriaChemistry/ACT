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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */
#ifndef ALEXANDRIA_PERCENT_MUTATOR
#define ALEXANDRIA_PERCENT_MUTATOR

#include "act/ga/genome.h"
#include "act/ga/mutator.h"

#include "acmindividual.h"

namespace alexandria
{


/*!
 * Changes values of the genes by a maximum of \p percent * their allowed range
 */
class PercentMutator : public ga::Mutator
{

private:
    //! IndividualInfo  pointer
    StaticIndividualInfo *sii_;
    //! The maximum change allowed as percent/100 of the range
    double                percent_;
public:

    /*!
     * Constructor
     * \param[in] sii     Pointer to StaticindividualInfo instance
     * \param[in] seed    seed passed by the middle man creating this
     * \param[in] percent Maximum allowed change
     */
    PercentMutator(StaticIndividualInfo *sii,
                   int                   seed,
                   double                percent)
    : ga::Mutator(seed), sii_(sii), percent_(percent) {}                

    /*! \brief Do the actual mutation
     * \param[inout] genome     The genome to mutate
     * \param[out]   bestGenome The best genome found
     * \param[in]    prMut      Probability for mutation
     */
    virtual void mutate(ga::Genome *genome,
                        ga::Genome *bestGenome,
                        double      prMut);

    virtual void sensitivityAnalysis(ga::Genome *, iMolSelect) {}
    
    //! \return whether a minimum was found
    bool foundMinimum() { return false; }
};


} //namespace alexandria


#endif //ALEXANDRIA_PERCENT_MUTATOR
