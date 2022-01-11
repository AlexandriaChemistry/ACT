/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#ifndef GA_GENETICALGORITHM_H
#define GA_GENETICALGORITHM_H

#include <cstdlib>

#include "Individual.h"
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
    //! Whether or not to evaluate the test set
    bool evaluateTestSet_ = false;

    //! Logfile for logging info
    FILE *logfile_          = nullptr;
     //! File for fitness train output
    FILE *fileFitnessTrain_ = nullptr;
    //! File for fitness test output
    FILE *fileFitnessTest_  = nullptr;

    //! Output environment (GROMACS)
    struct gmx_output_env_t *oenv_ = nullptr;

    //! Old population
    std::vector<Individual*>  oldPop_;
    //! New population, which emerges from the old population
    std::vector<Individual*>  newPop_;
    //! The best individual
    Individual               *bestInd_ = nullptr;

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
    GeneticAlgorithm(FILE                                *logFile,
                     struct gmx_output_env_t             *oenv,
                     Initializer                         *initializer,
                     FitnessComputer                     *fitnessComputer,
                     Sorter                              *sorter,
                     ProbabilityComputer                 *probComputer,
                     Selector                            *selector,
                     Crossover                           *crossover,
                     Mutator                             *mutator,
                     Terminator                          *terminator,
                     int                                  popSize,
                     bool                                 evaluateTestSet);
 
    //! \brief Evolve the initial population
    virtual void evolve() = 0;

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and Setters               *
    * * * * * * * * * * * * * * * * * * * * * */
    //! \return the population size
    int populationSize() const { return popSize_; }

    //! \return whether to evaluate the test set
    bool evaluateTestSet() const { return evaluateTestSet_; }

    //! \return the best individual
    Individual *bestInd() { return bestInd_; }

    /*! \brief Set the new best individual
     * \param[in] ind New best individual that will be moved to the storage.
     */
    void setBestIndividual(Individual *ind);

    //! \return the logfile
    FILE *logFile() { return logfile_; }

    //! Open fitness output
    void openFitnessFiles();

    //! And close them when the time is due
    void closeFitnessFiles();

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

    //! Return the output environment
    struct gmx_output_env_t *oenv() const { return oenv_; }

    //! \return a constant reference to \p oldPop_
    const std::vector<Individual*> &oldPop() const { return oldPop_; }

    //! \return a pointer to \p oldPop_
    std::vector<Individual*> *oldPopPtr() { return &oldPop_; }

    //! \return a pointer to \p newPop_
    std::vector<Individual*> *newPopPtr() { return &newPop_; }

    //! \brief Swap the old and new populations
    void swapOldNewPopulations();

    //! \brief Replace new by old populations
    void copyOldToNewPopulations();

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
     * BEGIN: Output routines                  *
     * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Print population to log file
    void fprintPop() const;

    //! \brief Print best individual to log file
    void fprintBestInd() const;

    //! \brief Print best individual (in current population) to log file
    void fprintBestIndInPop() const;

    //! \return the index of the Individual with the best fitness. FIXME: make this general. Now, the lower the fitness the better
    int findBestIndex() const;

    //! \brief Print the probability of each individual
    void fprintProbability() const;

    //! \brief Print the fitness of the population to the output files \p fileFitnessTrain_ and \p fileFitnessTest_
    void fprintFitness() const;

    /* * * * * * * * * * * * * * * * * * * * * *
     * END: Output routines                  *
     * * * * * * * * * * * * * * * * * * * * * */  
};

} //namespace ga

#endif // GA_GENETICALGORITHM_H
