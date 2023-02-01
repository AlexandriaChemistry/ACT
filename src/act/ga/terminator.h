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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_TERMINATOR_H
#define GA_TERMINATOR_H


#include <vector>
#include <limits>

#include "gene_pool.h"


namespace ga
{


/*!
 * \brief Abstract class to check for evolution termination conditions
 */
class Terminator
{

protected:

    //! File to print stuff to (may be stdout, or nullptr)
    FILE *outfile_ = nullptr;

    /*!
     * \brief Constructor
     * \param[in] outfile the file to write stuff to
     */
    Terminator(FILE *outfile)
    : outfile_(outfile) {}

public:

    /*!
     * \brief Check whether the evolution should be terminated
     * \param[in] pool              The GenePool
     * \param[in] generationNumber  The generation number
     * \return true if we should terminate, false otherwise
     */
    virtual bool terminate(const GenePool *pool,
                           const int       generationNumber) = 0;

};

/*!
 * \brief Terminator which stops evolution after a given amount of generations.
 */
class GenerationTerminator : public Terminator
{

private:

    //! Maximum allowed amount of generations
    int maxGenerations_;

public:

    /*!
     * \brief Constructor
     * \param[in] maxGenerations the maximum allowed amount of generations
     * \param[in] outfile        file for printing stuff
     */
    GenerationTerminator(const int maxGenerations, FILE *outfile)
    : Terminator(outfile), maxGenerations_(maxGenerations) {}

    /*!
     * \brief Will return true when \p generationNumber \f$\geq\f$ \p maxGenerations, and false otherwise.
     * \param[in] pool             The gene pool
     * \param[in] generationNumber The generation number
     * \return true if we should terminate, false otherwise
     */
    virtual bool terminate(const GenePool *pool,
                           const int       generationNumber);

};

/*!
 * \brief Terminator which stops evolution after the best test fitness has not
 * improved in a given number of generations.
 */
class TestGenTerminator : public Terminator
{

private:

    //! Amount of generations we allow the test fitness to not improve
    int generations_;
    //! Remaining generations to beat the best test fitness
    int remaining_;
    //! Best test fitness so far
    double bestFitness_;

public:

    /*!
     * \brief Constructor
     * \param[in] generations max amount of generations we allow the test
     *                        fitness to not improve
     * \param[in] outfile     file for printing stuff
     */
    TestGenTerminator(const int generations, FILE *outfile)
    : Terminator(outfile), generations_(generations), remaining_(generations),
      bestFitness_(std::numeric_limits<double>::max()) {}

    /*!
     * \brief Will return true when the test fitness has not improved over
     * the given amount of last generations, and false otherwise.
     * \param[in] pool             The gene pool
     * \param[in] generationNumber The generation number
     * \return true if we should terminate, false otherwise
     */
    virtual bool terminate(const GenePool *pool,
                           const int       generationNumber);

};


} //namespace ga


#endif //GA_TERMINATOR_H
