/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef ALEXANDRIA_MCMCMUTATOR_H
#define ALEXANDRIA_MCMCMUTATOR_H

#include "ga/Genome.h"
#include "ga/Mutator.h"
#include "confighandler.h"
#include "acmfitnesscomputer.h"
#include "acmindividual.h"

namespace alexandria
{


/*!
 * \brief Mutator which applies MCMC optimization to an ACMIndividual.
 * It can also conduct sensitivity analysis if requested.
 * FIXME: shouldn't we move sensitivity analysis somewhere else?
 */
class MCMCMutator : public ga::Mutator
{

private:
    //! Did we find a minimum yet?
    bool                  bMinimum_ = false;
    //! Pointer to BayesConfigHandler
    BayesConfigHandler   *bch_     = nullptr;
    //! Pointer to ACMFitnessComputer
    ACMFitnessComputer   *fitComp_ = nullptr;
    //! Pointer to Individual
    StaticIndividualInfo *sii_     = nullptr;
    //! The initial genome is stored for printing
    ga::Genome            initialGenome_;
    //! Attempted changes for each parameter
    std::vector<int>      attemptedMoves_;
    //! Accepted changes for each parameter
    std::vector<int>      acceptedMoves_;
    //! Mean of each parameter
    std::vector<double>   pMean_;
    //! Standard deviation of each parameter
    std::vector<double>   pSigma_;
    //! Convergence file for each parameter type
    std::vector<FILE*>    fpc_;
    //! Convergence file for Chi2
    FILE                 *fpe_ = nullptr;
    //! Pointer to log file (may be nullptr)
    FILE                 *logfile_ = nullptr;
    /*! Flush output immediately rather than letting the OS buffer it.
     * Don't use for production simulations.
     */
    bool                  verbose_ = false;

    // Random number generation
    std::random_device                      rd;
    std::mt19937                            gen;
    std::uniform_int_distribution<size_t>   dis;

    /*!
     * \brief Change force field parameter at a given index for an individual
     * \param[in] genome   pointer to the genome
     * \param[in] j     index of the parameter to change
     */
    void changeParam(ga::Genome *genome,
                     size_t  j);

    /*!
     * \brief Print new minimum to log file and, if necessary, print params to debug file
     * \param[in] chi2              the new minimum
     * \param[in] bEvaluate_testset true if test set is evaluated, false otherwise
     * \param[in] xiter             fractional iteration. E.g., if we are halfway through iteration 3, it is 3.5
     */
    void printNewMinimum(const std::map<iMolSelect, double> &chi2,
                         bool                                bEvaluate_testset,
                         double                              xiter);

    /*!
     * \brief Print parameter values to their respective surveillance files
     * \param[in] genome Pointer to the genome
     * \param[in] xiter  Fractional iteration. E.g., if we are halfway through iteration 3, it is 3.5
     */
    void printParameterStep(ga::Genome          *genome,
                            double               xiter);                                          

    /*!
     * \brief Write \f$ \chi^2 \f$ value of a genome to its convergence file, if it exists
     * \param[in] chi2                  chi2 for train and test sets
     * \param[in] bEvaluate_testset     true if test set is evaluated, false otherwise
     * \param[in] xiter                 fractional iteration. E.g., if we are halfway through iteration 3, it is 3.5
     */
    void printChi2Step(const std::map<iMolSelect, double> &chi2,
                       bool                                bEvaluate_testset,
                       double                              xiter);

    //! \return a random index of the force field parameter vector
    size_t randIndex() { return dis(gen); }

    /*!
     * \brief Compute mean and standard deviation for each force field parameter
     * \param[in] sum           over "nsum" iterations, the sum of each parameter
     * \param[in] nsum          number of iterations to compute statistics over
     * \param[in] sum_of_sq     over "nsum" iterations, the sum of each parameter squared
     */
    void computeMeanSigma(const std::vector<double>    &sum,
                          const int                     nsum,
                                std::vector<double>    *sum_of_sq);

    /*!
     * \brief Take a step of MCMC by attempting to alter a parameter
     * \param[inout] genome         pointer to a genome
     * \param[out] bestGenome       pointer to the best genome, to be filled
     * \param[in] changed           a reference to a vector which has true for parameters that change and false otherwise
     * \param[in] prevEval          pointer to a map with the previous \f$ \chi^2 \f$
     * \param[in] minEval           pointer to the minimum \f$ \chi^2 \f$ found so far for the training set
     * \param[in] bEvaluate_testset true if evaluation should be done on test set, false otherwise
     * \param[in] pp                index of inner loop over number of parameters
     * \param[in] iter              current iteration number
     * \param[in] beta0             pointer to beta for annealing
     */
    void stepMCMC(ga::Genome                   *genome,
                  ga::Genome                   *bestGenome,
                  std::vector<bool>            *changed,
                  std::map<iMolSelect, double> *prevEval,
                  std::map<iMolSelect, double> *minEval,
                  bool                          evaluate_testset,
                  size_t                        pp,
                  int                           iter,
                  double                       *beta0);

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: File stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Create a working directory
    void makeWorkDir();


    /* * * * * * * * * * * * * * * * * * * * * *
    * END: File stuff                          *
    * * * * * * * * * * * * * * * * * * * * * */
public:

    /*!
     * \brief Constructor
     * \param[in] logfile   pointer to log file (may be nullptr)
     * \param[in] verbose   Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
     * \param[in] bch       pointer to BayesConfigHandler object
     * \param[in] fitComp   pointer to ACMFitnessComputer object
     * \param[in] sii       pointer to StaticIndividualInfo object
     * \param[in] nParam    size of the force field parameter vector
     */
    MCMCMutator(FILE               *logfile,
                bool                verbose,
                BayesConfigHandler *bch,
                ACMFitnessComputer *fitComp,
                StaticIndividualInfo *sii);
 
    /*!
     * \brief Run the Markov chain Monte carlo (MCMC) simulation
     * \param[in]  genome     pointer to the genome
     * \param[out] bestGenome pointer to the best genome
     * \param[in]  prMut      Probability for mutation. 
     *                   Abused as a boolean for evaluating the test set here, if > 0.
     */
    virtual void mutate(ga::Genome *genome,
                        ga::Genome *bestGenome,
                        double      prMut);

    //! \return the number of calls to the objective function MCMCMutator::MCMC() routine
    size_t numberObjectiveFunctionCalls() const
    {
        return 1 + bch_->maxIter() * sii_->nParam();
    }

    //! \return whether a minimum was found
    bool foundMinimum() { return bMinimum_; }

    /*!
     * \brief Print the MC statistics to a file.
     * \param[in] fp         File pointer to print to
     * \param[in] bestGenome The genome with highest fitness
     */
    void printMonteCarloStatistics(FILE             *fp,
                                   const ga::Genome *bestGenome);

    /*!
     * \brief Perform a sensitivity analysis by systematically changing all parameters and re-evaluating the \f$ \chi^2 \f$.
     * \param[in] genome Pointer to genome
     * \param[in] ims    Dataset to perform sensitivity analysis on
     */
    void sensitivityAnalysis(ga::Genome *genome,
                             iMolSelect  ims);

    /*!
     * \brief Open parameter convergence files
     * \param[in] oenv the GROMACS output environment
     */
    void openParamConvFiles(const gmx_output_env_t *oenv);

    /*!
     * \brief Open a \f$ \chi^2 \f$ convergence file
     * \param[in] oenv              the GROMACS output environment
     * \param[in] bEvaluate_testset whether the test set will be evaluated in MCMC
     */
    void openChi2ConvFile(const gmx_output_env_t    *oenv,
                          const bool                 bEvaluate_testset);
    
    void stopHelpers();

    //! Close \f$ \chi^2 \f$ and parameter convergence files
    void finalize();
};


} //namespace alexandria


#endif //ALEXANDRIA_MCMCMUTATOR_H
