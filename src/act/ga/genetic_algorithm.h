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
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#ifndef GA_GENETICALGORITHM_H
#define GA_GENETICALGORITHM_H

#include <cstdlib>

#include <map>

#include "gene_pool.h"
#include "genome.h"
#include "initializer.h"
#include "fitness_computer.h"
#include "probability_computer.h"
#include "selector.h"
#include "crossover.h"
#include "mutator.h"
#include "terminator.h"
#include "penalizer.h"

struct gmx_output_env_t;
struct t_commrec;

namespace ga
{

/*!
 * \brief Class which encapsulates a genetic algorithm
 */
class GeneticAlgorithm
{

private:
    //! The population size
    int popSize_          = 0;

    //! Files for fitness output
    std::map<iMolSelect, FILE *>  fileFitness_;
    //! Initializes each individual in the population
    Initializer                  *initializer_  = nullptr;
    //! Computes fitness for each individual in the population
    FitnessComputer              *fitComputer_  = nullptr;
    //! Computes the probability of selection of each individual
    ProbabilityComputer          *probComputer_ = nullptr;
    //! Selects an individual from the population based on its probability
    Selector                     *selector_     = nullptr;
    //! Grabs 2 individuals and crosses their genes to generate 2 new individuals
    Crossover                    *crossover_    = nullptr;
    //! Grabs 1 individual and mutates its genes
    Mutator                      *mutator_      = nullptr;
    //! Checks if the evolution should continue or be terminated
    std::vector<Terminator*>     *terminators_  = nullptr;
    //! Penalizes the population
    std::vector<Penalizer*>      *penalizers_   = nullptr;

protected:
    //! Last population
    GenePool                      lastPop_;

public:

    //! \brief Default constructor
    GeneticAlgorithm() {}

    /*!
     * \brief Constructor for self-building
     */
    GeneticAlgorithm(Initializer                         *initializer,
                     FitnessComputer                     *fitnessComputer,
                     ProbabilityComputer                 *probComputer,
                     Selector                            *selector,
                     Crossover                           *crossover,
                     Mutator                             *mutator,
                     std::vector<Terminator*>            *terminators,
                     std::vector<Penalizer*>             *penalizers,
                     int                                  popSize) :
        popSize_(popSize), initializer_(initializer),
        fitComputer_(fitnessComputer), probComputer_(probComputer),
        selector_(selector), crossover_(crossover), mutator_(mutator),
        terminators_(terminators), penalizers_(penalizers)
    {
        if (terminators_)
        {
            GMX_RELEASE_ASSERT(!terminators_->empty(), "There are no terminators!");
        }
    }
    /*! \brief Update the genepool 
     * \param gpin New gene pool
     */
    void updateGenePool(const GenePool &gpin);
 
    /*! \brief Evolve the initial population
     * \param[out] bestGenome The best genome(s) found during the evolution (for different datasets, if applicable).
     *                        Comes in as an empty map, and subclasses of GeneticAlgorithm will fill it for the
     *                        datasets according to their configuration
     * \return whether a genome with better fitness (for training set) was found.
     */
    virtual bool evolve(std::map<iMolSelect, Genome> *bestGenome) = 0;

    /*! \brief Retrieve the last population (to be called after evolve,
     * otherwise undersired behavior will occur)
     * \return a constant reference to the last population
     */
    const GenePool &getLastPop() { return lastPop_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and Setters               *
    * * * * * * * * * * * * * * * * * * * * * */
    //! \return the population size
    int populationSize() const { return popSize_; }

    /*! Open fitness output
     * \param[in] filename Name for output files (data set will be inserted in file name)
     * \param[in] ims      Selection category
     */
    void openFitnessFiles(const std::string &filename,
                          iMolSelect         ims);

    //! And close them when the time is due
    void closeFitnessFiles();

    //! \return the fitness computer
    FitnessComputer *fitnessComputer() { return fitComputer_; }

    //! \return the probability computer
    ProbabilityComputer *probabilityComputer() { return probComputer_; }

    //! \return the initializer
    Initializer *initializer() { return initializer_; }

    //! \return the crossover
    Crossover *crossover() { return crossover_; }

    //! \return the mutator
    Mutator *mutator() { return mutator_; }

    //! \return the selector
    Selector *selector() { return selector_; }

    //! \return the terminators
    std::vector<Terminator*> *terminators() { return terminators_; }

    /*!
     * \param[in] index the index of the terminator
     * \return the terminator at the given index
     */
    Terminator *terminator(const int index);

    //! \return the penalizers
    std::vector<Penalizer*> *penalizers() { return penalizers_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Output routines                   *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! Return a FILE pointer for a data set
     * \param[in] ims The corresponding data set
     * \return fitness file for training
     */
    FILE *fitnessFile(iMolSelect ims) { return fileFitness_[ims]; }
    
    /*!
     * \brief Print the fitness of each genome in a pool
     * \param[in] pool the genome pool
     */
    void fprintFitness(const GenePool &pool);
    
    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Output routines                     *
    * * * * * * * * * * * * * * * * * * * * * */  

    /*!
     * \brief Check if we have to stop the evolution by asking each terminator.
     * Never call this method if no terminators were given to the GA.
     * \param[in] pool             the GenePool    
     * \param[in] generationNumber the current generation number
     * \return true if we stop the evolution, false otherwise
     */
    bool terminate(const GenePool *pool,
                   const int       generationNumber);

    /*!
     * \brief Penalize the population.
     * Never call this method if no terminators were given to the GA.
     * \param[in] pool             the GenePool    
     * \param[in] generationNumber the current generation number
     * \return true if the population has been penalized, false otherwise
     */
    bool penalize(      GenePool *pool,
                  const int       generationNumber);

};

} //namespace ga

#endif // GA_GENETICALGORITHM_H
