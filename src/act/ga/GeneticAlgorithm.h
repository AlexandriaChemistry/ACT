/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#ifndef GA_GENETICALGORITHM_H
#define GA_GENETICALGORITHM_H

#include <cstdlib>

#include <map>

#include "GenePool.h"
#include "Genome.h"
#include "Initializer.h"
#include "FitnessComputer.h"
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

     //! Files for fitness output
    std::map<iMolSelect, FILE *>  fileFitness_;
    //! The best genome
    Genome                        bestGenome_;
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
    Terminator                   *terminator_   = nullptr;

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
                     ProbabilityComputer                 *probComputer,
                     Selector                            *selector,
                     Crossover                           *crossover,
                     Mutator                             *mutator,
                     Terminator                          *terminator,
                     int                                  popSize) :
        popSize_(popSize), initializer_(initializer),
        fitComputer_(fitnessComputer), probComputer_(probComputer),
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

    //! \return the terminator
    Terminator *terminator() { return terminator_; }

    /* * * * * * * * * * * * * * * * * * * * * *
     * BEGIN: Output routines                  *
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
     * END: Output routines                  *
     * * * * * * * * * * * * * * * * * * * * * */  
};

} //namespace ga

#endif // GA_GENETICALGORITHM_H
