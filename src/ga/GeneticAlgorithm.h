#ifndef GA_GENETICALGORITHM_H
#define GA_GENETICALGORITHM_H


// #include "aliases.h"

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


namespace ga
{


/*!
 * Class which encapsulates a genetic algorithm
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

    //! Output environment (GROMACS)
    gmx_output_env_t *oenv_;

    //! Old population
    std::vector<Individual*> oldPop_;
    //! New population, which emerges from the old population
    std::vector<Individual*> newPop_;
    //! Temporal storage to swap "oldPop" and "newPop" after each generation
    std::vector<Individual*> tmpPop_;
    //! Probability of selection for each individual
    std::vector<double> probability_;

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

    //! Pure MCMC evaluation
    evolveMCMC();

public:

    //! Default constructor
    GeneticAlgorithm() {}

    /*!
     * Constructor for ACT
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
      probability_(gach->popSize()), oldPop_(gach->popSize()), newPop_(gach->popSize())
    {
        initializer_ = new alexandria::ACMInitializer(mg->mindata(), sii, gach->randomInit(), outputFile);
        fitComputer_ = new alexandria::ACMFitnessComputer(cr, logFile, sii, mg, removeMol, verbose, fullQuadrupole);
        mutator_ = new alexandria::MCMCMutator(logFile, verbose, bch, fitComputer, sii, sii.nParam());
    }

    //! \brief Evolve the initial population
    void evolve();

};


} //namespace ga


#endif //GA_GENETICALGORITHM_H
