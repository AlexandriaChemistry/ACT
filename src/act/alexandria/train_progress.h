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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ALEXANDRIA_TRAIN_PROGRESS_H
#define ALEXANDRIA_TRAIN_PROGRESS_H

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <map>
#include <string>

#include "act/basics/dataset.h"
#include "act/ga/genome.h"

namespace alexandria
{

/*! \brief Writes the progress of training to a file, for a program that follows it
 *
 * The file gets lines with one JSON object each, and each line is on disk
 * as soon as it is written. The key "event" says what a line is:
 *
 * - "start", once, when the first population has its chi squared: the
 *   optimizer, the number of generations ("ngenerations"), of Monte Carlo
 *   iterations in each generation ("niterations", zero for the genetic
 *   algorithm on its own), of individuals ("population") and of parameters
 *   ("nparameters"), and the best chi squared of the population ("chi2").
 * - "step", for steps of the Monte Carlo chain of the master process, a
 *   few times per second at most: the generation, from one, the iteration
 *   in that generation with the fraction of the parameters that were
 *   changed, the chi squared of the chain ("chi2") and the lowest one of
 *   the chain so far ("best").
 * - "generation", at the end of each generation: the generation, from one,
 *   the number of generations, the lowest chi squared of all individuals
 *   so far ("best") and whether this generation lowered it ("improved").
 *   Monte Carlo on its own is one generation.
 *
 * A chi squared is an object with the value for the data set "train" and,
 * when the optimizer computes it, for "test". The lines for steps and generations have
 * the time since the start in "elapsed", with its unit.
 *
 * Only the master process has such an object. The file is made by start(),
 * such that there is none when nothing is trained.
 */
class TrainProgress
{
private:
    //! The clock
    using Clock = std::chrono::steady_clock;
    //! The name of the file
    std::string       fileName_;
    //! The name of the optimizer
    std::string       optimizer_;
    //! The number of generations
    int               maxGenerations_;
    //! The number of Monte Carlo iterations in a generation
    int               maxIter_;
    //! The number of individuals
    int               popSize_;
    //! The number of parameters
    size_t            nParam_;
    //! Whether the individuals have a chi squared for the test set that is up to date
    bool              testSet_;
    //! The file, nullptr before start()
    FILE             *fp_           = nullptr;
    //! The generation of the steps, from one
    int               generation_   = 1;
    //! The shortest time between two lines for steps (s)
    double            stepInterval_ = 0.2;
    //! Whether a line for a step was written
    bool              haveStep_     = false;
    //! When start() was called
    Clock::time_point start_;
    //! When the last line for a step was written
    Clock::time_point lastStep_;

    //! \return the time since start() (s)
    double elapsed() const;

public:
    /*! \brief Constructor, which does not make the file yet
     * \param[in] fileName       The name of the file
     * \param[in] optimizer      The name of the optimizer
     * \param[in] maxGenerations The number of generations
     * \param[in] maxIter        The number of Monte Carlo iterations in a generation
     * \param[in] popSize        The number of individuals
     * \param[in] nParam         The number of parameters
     * \param[in] testSet        Whether the optimizer computes the chi squared of the
     *                           test set for the individuals, without that only the one
     *                           of the training set is reported
     */
    TrainProgress(const std::string &fileName,
                  const std::string &optimizer,
                  int                maxGenerations,
                  int                maxIter,
                  int                popSize,
                  size_t             nParam,
                  bool               testSet);

    //! Destructor, closes the file
    ~TrainProgress();

    //! No copies, there is one file
    TrainProgress(const TrainProgress &) = delete;
    //! No copies, there is one file
    TrainProgress &operator=(const TrainProgress &) = delete;

    /*! \brief The chi squared of a genome, as the functions below take it
     * \param[in] genome The genome
     * \return the fitness for the data sets train and test, as far as the genome has it
     */
    static std::map<iMolSelect, double> chi2(const ga::Genome &genome);

    /*! \brief Make the file and write the line for the start
     * \param[in] chi2 The best chi squared of the first population, for the
     *                 data sets for which it is known
     */
    void start(const std::map<iMolSelect, double> &chi2);

    /*! \brief Set the generation of the steps that follow
     * \param[in] generation The generation, from one
     */
    void setGeneration(int generation) { generation_ = generation; }

    /*! \brief Set the shortest time between two lines for steps
     * \param[in] seconds The time (s), zero for a line for every step
     */
    void setStepInterval(double seconds) { stepInterval_ = seconds; }

    /*! \brief Write a line for a step of the Monte Carlo chain, unless the last one is too recent
     * Does nothing before start().
     * \param[in] iteration The iteration in the generation, with the fraction of the parameters done
     * \param[in] chi2      The chi squared of the chain after this step
     * \param[in] best      The lowest chi squared of the chain so far, of which the
     *                      one of the training set is reported
     */
    void step(double                              iteration,
              const std::map<iMolSelect, double> &chi2,
              const std::map<iMolSelect, double> &best);

    /*! \brief Write the line for the end of a generation
     * Does nothing before start().
     * \param[in] generation The generation, from one
     * \param[in] best       The lowest chi squared of all individuals so far
     * \param[in] improved   Whether this generation lowered it
     */
    void generation(int                                 generation,
                    const std::map<iMolSelect, double> &best,
                    bool                                improved);
};

} // namespace alexandria

#endif
