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
#ifndef GA_PENALIZER_H
#define GA_PENALIZER_H

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */

#include "gene_pool.h"
#include "initializer.h"

#include "gromacs/fileio/xvgr.h"

namespace gmx
{
class TextWriter;
}

namespace ga
{

/*!
 * \brief Penalizes a population
 */
class Penalizer
{

protected:

    //! \brief Default constructor
    Penalizer() {}

public:

    //! \brief Default base destructor
    virtual ~Penalizer() {}

    /*!
     * \brief Penalize a population
     * \param[in]    tw         TextWriter
     * \param[inout] pool       the collection of genomes
     * \param[in]    generation the current generation number
     * \return true if the population has been penalized, false otherwise
     */
    virtual bool penalize(gmx::TextWriter *tw,
                          GenePool        *pool,
                          const int        generation) = 0;

}; // class Penalizer

/*!
 * \brief Penalizer that punishes the population if its volume is not large enough
 *
 * We create a parallelogram whose side lengths are defined by the maximum and
 * minimum value of each parameter in the population.
 * If the volume of the parallelogram divided by the total volume of the parameter
 * space is smaller than a given fraction, we reinitialize a given fraction of the
 * worst genomes.
 */
class VolumeFractionPenalizer : public Penalizer
{

private:

    //! Whether we compute the volume in log scale
    bool logVolume_;
    //! Total volume of the parameter space (may be log scale)
    double totalVolume_;
    //! Limit of fraction of volume
    double volFracLimit_;
    //! Fraction of the population to reinitialize (worst genomes)
    double popFrac_;
    //! Initializer to randomize genomes
    Initializer *initializer_;

    //! Output file
    FILE *outfile_ = nullptr;

    //! \return the volume of a population
    double getPoolVolume(const GenePool &pool) const;

public:

    //! \brief Destructor: closes output file
    ~VolumeFractionPenalizer();

    /*!
     * \brief Create a new VolumeFractionPenalizer
     * \param[in] oenv         gromacs output environment
     * \param[in] logVolume    true if volume in log scale, false otherwise
     * \param[in] logfile      file to print log stuff to (may be nullptr)
     * \param[in] totalVolume  total volume of the parameter space
     * \param[in] volFracLimit if the volume of the population divided by
     *                         \p totalVolume is smaller than \p volFracLimit then
     *                         we punish
     * \param[in] popFrac      fraction of the population to penalize (worst genomes)
     * \param[in] initializer  Initializer to randomize genomes
     */
    VolumeFractionPenalizer(      gmx_output_env_t  *oenv,
                            const bool               logVolume,
                            const double             totalVolume,
                            const double             volFracLimit,
                            const double             popFrac,
                                  Initializer       *initializer);

    bool penalize(gmx::TextWriter *tw,
                  GenePool        *pool,
                  const int       generation) override;

}; // class VolumeFactorPenalizer

/*!
 * \brief A catastrophe happens and a fraction of the population (at random) is
 * randomized
 * 
 * Each X generations a fraction Y (at random) of the population is randomized
 * TODO: if needed, make this penalizer print the indices of the individuals it is
 * randomizing to a file
 */
class CatastrophePenalizer : public Penalizer
{

private:

    // Random number stuff (for shuffling)
    std::random_device  rd;
    std::mt19937        gen;

    //! Interval of generations between the penalizer triggers
    int genInterval_;
    //! Fraction of the population to reinitialize
    double popFrac_;
    //! Initializer to randomize genomes
    Initializer *initializer_;
    //! List of all existing indices for genomes in the population (0 to popSize-1)
    std::vector<size_t> availableIndices_;

public:

    /*!
     * \brief Create a new CatastrophePenalizer
     * \param[in] seed        seed for the random number generation
     * \param[in] genInterval number of generations as interval between penalties
     * \param[in] popFrac     fraction of the population to penalize
     * \param[in] initializer Initializer to randomize genomes
     * \param[in] popSize     Size of the population
     */
    CatastrophePenalizer(const int          seed,
                         const int          genInterval,
                         const double       popFrac,
                               Initializer *initializer,
                         const size_t       popSize);

    bool penalize(gmx::TextWriter *tw,
                  GenePool        *pool,
                  const int        generation) override;

}; // class CatastrophePenalizer

} // namespace ga


#endif // GA_PENALIZER_H
