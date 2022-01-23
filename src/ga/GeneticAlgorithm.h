/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#ifndef GA_GENETICALGORITHM_H
#define GA_GENETICALGORITHM_H

#include <cstdlib>

//#include "GenePool.h"
#include "Genome.h"
#include "Initializer.h"
#include "FitnessComputer.h"
#include "Sorter.h"
#include "ProbabilityComputer.h"
#include "Selector.h"
#include "Crossover.h"
#include "Mutator.h"
#include "Terminator.h"

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

     //! File for fitness train output
    FILE *fileFitnessTrain_ = nullptr;
    //! File for fitness test output
    FILE *fileFitnessTest_  = nullptr;
    //! The best genome
    Genome                  bestGenome_;
    //! Initializes each individual in the population
    Initializer            *initializer_ = nullptr;
    //! Computes fitness for each individual in the population
    FitnessComputer        *fitComputer_ = nullptr;
    //! Sorts the individuals based on their fitness
    Sorter                 *sorter_ = nullptr;
    //! Computes the probability of selection of each individual
    ProbabilityComputer    *probComputer_ = nullptr;
    //! Selects an individual from the population based on its probability
    Selector               *selector_ = nullptr;
    //! Grabs 2 individuals and crosses their genes to generate 2 new individuals
    Crossover              *crossover_ = nullptr;
    //! Mutates the genes of the individuals
    Mutator                *mutator_ = nullptr;
    //! Checks if the evolution should continue or be terminated
    Terminator             *terminator_ = nullptr;

    // FIXME: Something could be done about generalizing the evolution.
    // We could make all Individuals
    // have their own parameter and fitness convergence files. 
    // IDK if this makes sense...
public:

    //! \brief Default constructor
    GeneticAlgorithm() {}

    /*!
     * \brief Constructor for self-building
     */
    GeneticAlgorithm(Initializer                         *initializer,
                     FitnessComputer                     *fitnessComputer,
                     Sorter                              *sorter,
                     ProbabilityComputer                 *probComputer,
                     Selector                            *selector,
                     Crossover                           *crossover,
                     Mutator                             *mutator,
                     Terminator                          *terminator,
                     int                                  popSize) :
        popSize_(popSize), initializer_(initializer),
        fitComputer_(fitnessComputer), sorter_(sorter), probComputer_(probComputer),
        selector_(selector), crossover_(crossover), mutator_(mutator), terminator_(terminator) {}

 
    /*! \brief Evolve the initial population
     * \param[out] bestGenome The best genome found during the evolution
     * \return whether a genome with better fitness was found.
     */
    virtual bool evolve(Genome *bestGenome) = 0;

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and Setters               *
    * * * * * * * * * * * * * * * * * * * * * */
    //! \return the population size
    int populationSize() const { return popSize_; }

    //! Open fitness output
    void openFitnessFiles();

    //! And close them when the time is due
    void closeFitnessFiles();

    //! \return pointer to bestGenome
    ga::Genome *bestGenomePtr() { return &bestGenome_; }
    
    //! \return constant best genome
    const ga::Genome &bestGenome() { return bestGenome_; }
    
    //! \return the mutator
    Mutator *mutator() { return mutator_; }

    //! \return the fitness computer
    FitnessComputer *fitnessComputer() { return fitComputer_; }

    //! \return the probability computer
    ProbabilityComputer *probabilityComputer() { return probComputer_; }

    //! \return the initializer
    Initializer *initializer() { return initializer_; }

    //! \return the crossover
    Crossover *crossover() { return crossover_; }

    //! \return the selector
    Selector *selector() { return selector_; }

    //! \return the sorter
    Sorter *sorter() { return sorter_; }

    //! \return the terminator
    Terminator *terminator() { return terminator_; }

    /* * * * * * * * * * * * * * * * * * * * * *
     * BEGIN: Output routines                  *
     * * * * * * * * * * * * * * * * * * * * * */

    //! \return fitness file for training
    FILE *fitnessTrain() { return fileFitnessTrain_; }
    
    //! \return fitness file for testing
    FILE *fitnessTest() { return fileFitnessTest_; }
    
    /* * * * * * * * * * * * * * * * * * * * * *
     * END: Output routines                  *
     * * * * * * * * * * * * * * * * * * * * * */  
};

} //namespace ga

#endif // GA_GENETICALGORITHM_H
