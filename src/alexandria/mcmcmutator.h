/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef ALEXANDRIA_MCMCMUTATOR_H
#define ALEXANDRIA_MCMCMUTATOR_H


#include "ga/Mutator.h"
#include "confighandler.h"
#include "acmfitnesscomputer.h"
#include "acmindividual.h"


namespace alexandria
{


class MCMCMutator : public ga::Mutator
{

private:

    //! Pointer to BayesConfigHandler
    BayesConfigHandler *bch_;
    //! Pointer to ACMFitnessComputer
    ACMFitnessComputer *fitComp_;
    //! Pointer to SharedIndividualInfo
    SharedIndividualInfo *sii_;
    //! Pointer to log file (may be nullptr)
    FILE *logfile_;
    //! Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
    bool verbose_;

    std::random_device                      rd;
    std::mt19937                            gen;
    std::uniform_int_distribution<size_t>   dis;

    /*! \brief
     * Change parameter \p j in an individual based on a random number
     * obtained from a uniform distribution.
     * \param[in] ind   pointer to the individual
     * \param[in] j     index of the parameter to change
     */
    void changeParam(ACMIndividual *ind,
                     size_t         j);

    /*!
     * Print new minimum to log file and, if necessary, print params to debug file
     * @param ind                   pointer to individual
     * @param bEvaluate_testset     true if test set is evaluated, false otherwise
     * @param xiter                 fractional iteration. E.g, if we are halfway through iteration 3 it is 3.5
     * @param currEval              current chi2 in training set
     * @param currEval_testset      current chi2 in test set
     */
    void fprintNewMinimum(      ACMIndividual  *ind,
                          const bool            bEvaluate_testset,
                          const double          xiter,
                          const double          currEval,
                          const double          currEval_testset);     

    /*!
     * Print parameter values to their respective surveillance files
     * @param ind   pointer to the individual
     * @param xiter fractional iteration (e.g. 3.5, 2.89, ...)
     */
    void fprintParameterStep(      ACMIndividual   *ind,
                             const double           xiter);                                          

    /*!
     * Write chi2 value to surveillance file, if it exists
     * @param ind                   pointer to individual
     * @param bEvaluate_testset     true if test set is evaluated, false otherwise
     * @param xiter                 fractional iteration (3.6, 3.89, ...)
     * @param prevEval              chi2 fro training set
     * @param prevEval_testset      chi2 for test set
     */
    void fprintChi2Step(      ACMIndividual    *ind,
                        const bool              bEvaluate_testset,
                        const double            xiter,
                        const double            prevEval,
                        const double            prevEval_testset);

    //! \return a random index of the parameter vector
    size_t randIndex() { return dis(gen); }

    /*!
     * Compute mean (pmean_) and standard deviation (psigma_) for each parameter
     * @param pmean         pointer to \p pmean_ vector in the individual
     * @param psigma        pointer to \p psigma_ vector in the individual
     * @param nParam        number of parameters in the system
     * @param sum           over "nsum" iterations, the sum of each parameter
     * @param nsum          number of iterations to compute statistics over
     * @param sum_of_sq     over "nsum" iterations, the sum of each parameter squared
     */
    void computeMeanSigma(      std::vector<double>    *pmean,
                                std::vector<double>    *psigma,
                          const size_t                  nParam,
                          const std::vector<double>    &sum,
                          const int                     nsum,
                                std::vector<double>    *sum_of_sq);

    /*!
     * Take a step of MCMC by attempting to alter a parameter
     * @param ind               pointer to individual
     * @param param             pointer to parameter vector of the individual
     * @param changed           a reference to a vector which has true for parameters that change and false otherwise
     * @param prevEval          pointer to a double storage with the previous chi2 for training set
     * @param prevEval_testset  a pointer to a double storage with the previous chi2 for test set
     * @param bEvaluate_testset true if evaluation should be done on test set, false otherwise
     * @param pp                index of inner loop over number of parameters
     * @param iter              current iteration number
     * @param beta0             pointer to beta for annealing
     * @param nParam            number of parameters in the model
     * @param minEval           pointer to the minimum chi2 found so far for the training set
     * @param paramClassIndex   class (by index) of each parameter in the model
     */
    void stepMCMC(      ACMIndividual          *ind,
                        std::vector<double>    *param,
                        std::vector<bool>      *changed,
                        double                 *prevEval,
                        double                 *prevEval_testset,
                  const bool                    evaluate_testset,
                  const size_t                  pp,
                  const int                     iter,
                        double                 *beta0,
                  const size_t                  nParam,
                        double                 *minEval,
                  const std::vector<size_t>    &paramClassIndex);

    /*!
     * Perform a mutation step
     * @param ind               pointer to individual
     * @param param             pointer to parameter vector of the individual
     * @param changed           a reference to a vector which has true for parameters that change and false otherwise
     * @param prevEval          pointer to a double storage with the previous chi2 for training set
     * @param pp                index of inner loop over number of parameters
     * @param iter              current iteration number
     * @param beta0             pointer to beta for annealing
     * @param nParam            number of parameters in the model
     * @param paramClassIndex   class (by index) of each parameter in the model
     */
    void stepMutation(      ACMIndividual          *ind,
                            std::vector<double>    *param,
                            std::vector<bool>      *changed,
                            double                 *prevEval,
                      const size_t                  pp,
                      const int                     iter,
                            double                 *beta0,
                      const size_t                  nParam,
                      const std::vector<size_t>    &paramClassIndex);

public:

    /*!
     * Constructor of MCMCMutator
     * @param logfile   pointer to log file (may be nullptr)
     * @param verbose   Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
     * @param bch       pointer to BayesConfigHandler object
     * @param fitComp   pointer to ACMFitnessComputer object
     * @param sii       pointer to SharedIndividualInfo object
     * @param nParam    size of the parameter vector
     */
    MCMCMutator(      FILE                     *logfile,
                const bool                      verbose,
                      BayesConfigHandler       *bch,
                      ACMFitnessComputer       *fitComp,
                      SharedIndividualInfo     *sii,
                const size_t                    nParam)
    : Mutator(), gen(rd()), dis(std::uniform_int_distribution<size_t>(0, nParam-1))
    {
        gen.seed(::time(NULL));

        logfile_ = logfile;
        verbose_ = verbose;
        bch_     = bch;
        fitComp_ = fitComp;
        sii_     = sii;
    };

    /*!
     * Mutate an individual's genes (in place). Only used when combined with GA.
     * @param ind       pointer to the individual to mutate
     * @param prMut     probability of mutating a gene
     */
    virtual void mutate(      ga::Individual   *individual,
                        const double            prMut);

    /*! \brief
     * Run the Markov chain Monte carlo (MCMC) simulation
     * \param[in]  ind              pointer to the individual
     * \param[in]  evaluate_testset If true, evaluate the energy on
     *                              the test set.
     */
    bool MCMC(      ACMIndividual  *ind,
              const bool            evaluate_testset);

    /*! \brief
     * Return the number of calls to the objective function
     * that will be made by the MCMC routine
     */
    size_t numberObjectiveFunctionCalls() const
    {
        return 1 + bch_->maxIter() * sii_->nParam();
    }

    /*! \brief
     * Print the MC statistics to a file.
     * \param[in] ind   pointer to the indivifdual
     * \param[in] fp    File pointer to print to
     */
    void printMonteCarloStatistics(ACMIndividual   *ind,
                                   FILE            *fp);

    /*! \brief
     * Perform a sensitivity analysis by systematically changing
     * all parameters and re-evaluating the chi2.
     * \param[in] ind   pointer to individual
     * \param[in] ims   Data set to perform sensitivity analysis on
     */
    void sensitivityAnalysis(ACMIndividual  *ind,
                             iMolSelect      ims);

};


} //namespace alexandria


#endif //ALEXANDRIA_MCMCMUTATOR_H