/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2025
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
#ifndef GA_MUTATOR_H
#define GA_MUTATOR_H

#include <random>
#include <time.h>

#include "act/basics/chargemodel.h"
#include "act/basics/msg_handler.h"
#include "act/ga/genome.h"

namespace ga
{
/*!
 * \brief Abstract class for gene mutation of Genome objects
 */
class Mutator
{
private:
    //! Ranom number generation
    std::random_device                      rd_base;
    std::mt19937                            gen_base;
    std::uniform_real_distribution<double>  dis_base;
    //! Charge generation
    alexandria::ChargeGenerationAlgorithm   algorithm_;
protected:

    /*! \brief Constructor
     * \param[in] seed      seed for the random number generator
     * \param[in] algorithm The algorithm used to generate partial charges
     */
    Mutator(int                                    seed,
            alexandria::ChargeGenerationAlgorithm  algorithm)
    : gen_base(rd_base()), dis_base(std::uniform_real_distribution<double>(0.0, 1.0))
    {
        gen_base.seed(seed);
        algorithm_ = algorithm;
    }

    //! \return a random number in \f$ [0, 1] \f$
    double randNum() { return dis_base(gen_base); }

public:

    virtual ~Mutator() = default;

    /*!
     * \brief Mutate genes of a Genome (in place)
     * \param[in]  msghandler   The message and status handler
     * \param[inout] genome     Pointer to the genome to mutate
     * \param[out]   bestGenome Pointer to the best genome found
     * \param[in]    prMut      Probability of mutating a gene
     */
    virtual void mutate(alexandria::MsgHandler *msghandler,
                        Genome                 *genome,
                        Genome                 *bestGenome,
                        double                  prMut) = 0;

    /*!
     * \brief Perform a sensitivity analysis by systematically changing all parameters and
     * re-evaluating the \f$ \chi^2 \f$.
     * \param[in] msghandler The message and status handler
     * \param[in] bestGenome Pointer to genome
     * \param[in] ims        Dataset to perform sensitivity analysis on
     */
    virtual void sensitivityAnalysis(alexandria::MsgHandler *msghandler,
                                     Genome                 *bestGenome,
                                     iMolSelect              ims) = 0;

    //! \return whether a minimum was found
    virtual bool foundMinimum() = 0;

    //! \return the ChargeGenerationAlgorithm
    alexandria::ChargeGenerationAlgorithm algorithm() const { return algorithm_; }
};


} //namespace ga


#endif //GA_MUTATOR_H
