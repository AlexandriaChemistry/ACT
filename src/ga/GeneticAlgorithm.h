/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef GA_GENETICALGORITHM_H
#define GA_GENETICALGORITHM_H


#include <cstdlib>

#include "Initializer.h"
#include "FitnessComputer.h"
#include "Sorter.h"
#include "ProbabilityComputer.h"
#include "Selector.h"
#include "Crossover.h"
#include "Mutator.h"
#include "Terminator.h"

#include "alexandria/confighandler.h"
#include "alexandria/sharedindividualinfo.h"
#include "alexandria/acmfitnesscomputer.h"
#include "alexandria/acminitializer.h"
#include "alexandria/mcmcmutator.h"
#include "alexandria/percentmutator.h"
#include "alexandria/npointcrossover.h"


namespace ga
{


/*!
 * \brief Class which encapsulates a genetic algorithm
 */
class GeneticAlgorithm
{

private:

    //! BayesConfigHandler pointer
    alexandria::BayesConfigHandler *bch_;
    //! GAConfigHandler pointer
    alexandria::GAConfigHandler *gach_;
    //! Communcations record
    t_commrec *cr_;

    //! Logfile for logging info
    FILE *logfile_;
    //! Name for fitness train output file
    const char *const filenameFitnessTrain_ = "ga_fitness_train.txt";
    //! File for fitness train output
    FILE *fileFitnessTrain_ = NULL;
    //! Name for fitness test output file
    const char *const filenameFitnessTest_ = "ga_fitness_test.txt";
    //! File for fitness test output
    FILE *fileFitnessTest_ = NULL;

    //! Output environment (GROMACS)
    gmx_output_env_t *oenv_;

    //! Old population
    std::vector<Individual*> oldPop_;
    //! New population, which emerges from the old population
    std::vector<Individual*> newPop_;
    //! Temporal storage to swap "oldPop" and "newPop" after each generation
    std::vector<Individual*> tmpPop_;
    //! The best individual
    Individual* bestInd_;


    //! Initializes each individual in the population
    Initializer            *initializer_;
    //! Computes fitness for each individual in the population
    FitnessComputer        *fitComputer_;
    //! Sorts the individuals based on their fitness
    Sorter                 *sorter_;
    //! Computes the probability of selection of each individual
    ProbabilityComputer    *probComputer_;
    //! Selects an individual from the population based on its probability
    Selector               *selector_;
    //! Grabs 2 individuals and crosses their genes to generate 2 new individuals
    Crossover              *crossover_;
    //! Mutates the genes of the individuals
    Mutator                *mutator_;
    //! Checks if the evolution should continue or be terminated
    Terminator             *terminator_;


    // FIXME: Something could be done about generalizing the evolution. We could make all Individuals
    // have their own parameter and fitness convergence files. IDK if this makes sense...

    //! \brief Pure MCMC evaluation
    void evolveMCMC();

    //! \brief Regular GA evolution (could be HYBRID too)
    void evolveGA();

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

public:

    //! \brief Default constructor
    GeneticAlgorithm() {}

    /*!
     * \brief Constructor for self-building
     */
    GeneticAlgorithm(const  bool                                 verbose,
                     const  bool                                 removeMol,
                     const  bool                                 fullQuadrupole,
                            t_commrec                           *cr,
                            MolGen                              *mg,
                            FILE                                *logFile,
                            gmx_output_env_t                    *oenv,
                            alexandria::BayesConfigHandler      *bch,
                            alexandria::SharedIndividualInfo    *sii,
                            alexandria::GAConfigHandler         *gach,
                     const  std::string                         &outputFile)
    : bch_(bch), gach_(gach), cr_(cr), logfile_(logFile), oenv_(oenv),
      oldPop_(gach->popSize()), newPop_(gach->popSize())
    {

        // Initializer
        initializer_ = new alexandria::ACMInitializer(mg->mindata(), sii, gach->randomInit(), outputFile);
        
        // FitnessComputer
        ACMFitnessComputer* tmpACMFitComp = new alexandria::ACMFitnessComputer(cr, logFile, sii, mg, removeMol, verbose, fullQuadrupole);
        fitComputer_ = tmpACMFitComp;
        
        // Mutator
        if (strcmp(gach->optimizer(), "GA") == 0)
        {
            mutator_ = new alexandria::PercentMutator(sii, gach->percent());
        }
        else
        {
            mutator_ = new alexandria::MCMCMutator(logFile, verbose, bch, tmpACMFitComp, sii, sii->nParam());
        }

        // If GA or HYBRID have been selected as optimizers, intialize the rest of the elements
        if (strcmp(gach->optimizer(), "MCMC") != 0)
        {

            // Sorter
            if (strcmp(gach->sorter(), "QUICK") == 0)
            {
                sorter_ = new QuickSorter(false);
            }
            else if (strcmp(gach->sorter(), "MERGE") == 0)
            {
                sorter_ = new MergeSorter(gach->popSize(), false);
            }
            else  // No sorting requested
            {
                sorter_ = new EmptySorter();
            }

            // ProbabilityComputer
            if (strcmp(gach->probComputer(), "RANK") == 0)
            {
                probComputer_ = new RankProbabilityComputer(gach->popSize());
            }
            else if (strcmp(gach->probComputer(), "FITNESS") == 0)
            {
                probComputer_ = new FitnessProbabilityComputer();
            }
            else  // BOLTZMANN
            {
                probComputer_ = new BoltzmannProbabilityComputer(gach->popSize(),
                                                                 gach->boltzTemp());
            }

            // Selector
            selector_ = new RouletteSelector();

            // Crossover
            GMX_RELEASE_ASSERT( ( (unsigned int) gach->nCrossovers() ) < sii->nParam(),
                               gmx::formatString("The order of the crossover operator should be smaller than the amount of parameters. You chose -nCrossovers %i, but there are %lu parameters. Please adjust -nCrossovers.", gach->nCrossovers(), sii->nParam()).c_str() );
            crossover_ = new NPointCrossover(sii->nParam(), gach->nCrossovers());

            // Terminator
            terminator_ = new GenerationTerminator(gach->maxGenerations());

        }

        // Create directories for each individual
        for (int i = 0; i <= gach->popSize(); i++)
        {
            const std::string command = "mkdir ind" + std::to_string(i);
            system(command.c_str());
            // std::filesystem::create_directory(dirName.c_str());
        }

    }

    //! \brief Evolve the initial population
    void evolve();

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and Setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return the best individual
    Individual *bestInd() { return bestInd_; }

    //! \return the mutator
    Mutator *mutator() { return mutator_; }

    //! \return the fitness computer
    FitnessComputer *fitComputer() { return fitComputer_; }

    //! \return a constant reference to \p oldPop_
    const std::vector<Individual*> &oldPop() const { return oldPop_; }

    //! \return a pointer to \p oldPop_
    std::vector<Individual*> *oldPopPtr() { return &oldPop_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace ga


#endif //GA_GENETICALGORITHM_H
