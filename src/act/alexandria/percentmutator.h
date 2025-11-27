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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */
#ifndef ALEXANDRIA_PERCENT_MUTATOR
#define ALEXANDRIA_PERCENT_MUTATOR

#include "act/basics/msg_handler.h"
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
     * \param[in] sii       Pointer to StaticindividualInfo instance
     * \param[in] seed      seed passed by the middle man creating this
     * \param[in] algorithm The charge generation algorithm
     * \param[in] percent   Maximum allowed change
     */
    PercentMutator(StaticIndividualInfo                  *sii,
                   int                                    seed,
                   alexandria::ChargeGenerationAlgorithm  algorithm,
                   double                                 percent)
        : ga::Mutator(seed, algorithm), sii_(sii), percent_(percent) {}                

    //! \copydoc ga::Mutator::mutate
    virtual void mutate(MsgHandler *msghandler,
                        ga::Genome *genome,
                        ga::Genome *bestGenome,
                        double      prMut);

    //! \copydoc ga::Mutator::sensitivityAnalysis
    virtual void sensitivityAnalysis(gmx_unused MsgHandler *msghandler,
                                     gmx_unused ga::Genome *bestGenome,
                                     gmx_unused iMolSelect  ims) {}
    
    //! \return whether a minimum was found
    bool foundMinimum() { return false; }
};


} //namespace alexandria


#endif //ALEXANDRIA_PERCENT_MUTATOR
