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


/*!
 * \brief Mutator which applies MCMC optimization to an ACMIndividual.
 * It can also conduct sensitivity analysis if requested.
 * FIXME: shouldn't we move sensitivity analysis somewhere else?
 */
class MCMCMutator : public ga::Mutator
{

private:
    //! Did we find a minimum yet?
    bool                bMinimum_ = false;
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

    // Random number generation
    std::random_device                      rd;
    std::mt19937                            gen;
    std::uniform_int_distribution<size_t>   dis;

    /*!
     * \brief Change force field parameter at a given index for an individual
     * \param[in] ind   pointer to the individual
     * \param[in] j     index of the parameter to change
     */
    void changeParam(ACMIndividual *ind,
                     size_t         j);

    /*!
     * \brief Print new minimum to log file and, if necessary, print params to debug file
     * \param[in] ind                   pointer to individual
     * \param[in] bEvaluate_testset     true if test set is evaluated, false otherwise
     * \param[in] xiter                 fractional iteration. E.g., if we are halfway through iteration 3, it is 3.5
     * \param[in] currEval              current chi2 in training set
     * \param[in] currEval_testset      current chi2 in test set
     */
    void fprintNewMinimum(      ACMIndividual  *ind,
                          const bool            bEvaluate_testset,
                          const double          xiter,
                          const double          currEval,
                          const double          currEval_testset);     

    /*!
     * \brief Print parameter values to their respective surveillance files
     * \param[in] ind   pointer to the individual
     * \param[in] xiter fractional iteration. E.g., if we are halfway through iteration 3, it is 3.5
     */
    void fprintParameterStep(      ACMIndividual   *ind,
                             const double           xiter);                                          

    /*!
     * \brief Write \f$ \chi^2 \f$ value of an individual to its convergence file, if it exists
     * \param[in] ind                   pointer to individual
     * \param[in] bEvaluate_testset     true if test set is evaluated, false otherwise
     * \param[in] xiter                 fractional iteration. E.g., if we are halfway through iteration 3, it is 3.5
     * \param[in] prevEval              \f$ \chi^2 \f$ for training set
     * \param[in] prevEval_testset      \f$ \chi^2 \f$ for test set
     */
    void fprintChi2Step(      ACMIndividual    *ind,
                        const bool              bEvaluate_testset,
                        const double            xiter,
                        const double            prevEval,
                        const double            prevEval_testset);

    //! \return a random index of the force field parameter vector
    size_t randIndex() { return dis(gen); }

    /*!
     * \brief Compute mean and standard deviation for each force field parameter
     * \param[in] pmean         pointer to mean vector in the individual
     * \param[in] psigma        pointer to standard deviation vector in the individual
     * \param[in] nParam        number of parameters in the system
     * \param[in] sum           over "nsum" iterations, the sum of each parameter
     * \param[in] nsum          number of iterations to compute statistics over
     * \param[in] sum_of_sq     over "nsum" iterations, the sum of each parameter squared
     */
    void computeMeanSigma(      std::vector<double>    *pmean,
                                std::vector<double>    *psigma,
                          const size_t                  nParam,
                          const std::vector<double>    &sum,
                          const int                     nsum,
                                std::vector<double>    *sum_of_sq);

    /*!
     * \brief Take a step of MCMC by attempting to alter a parameter
     * \param[in] ind               pointer to individual
     * \param[in] param             pointer to parameter vector of the individual
     * \param[in] changed           a reference to a vector which has true for parameters that change and false otherwise
     * \param[in] prevEval          pointer to a double storage with the previous \f$ \chi^2 \f$ for training set
     * \param[in] prevEval_testset  a pointer to a double storage with the previous \f$ \chi^2 \f$ for test set
     * \param[in] bEvaluate_testset true if evaluation should be done on test set, false otherwise
     * \param[in] pp                index of inner loop over number of parameters
     * \param[in] iter              current iteration number
     * \param[in] beta0             pointer to beta for annealing
     * \param[in] nParam            number of parameters in the model
     * \param[in] minEval           pointer to the minimum \f$ \chi^2 \f$ found so far for the training set
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
                        double                 *minEval);

    /*!
     * \brief Perform a mutation step
     * \param[in] ind               pointer to individual
     * \param[in] param             pointer to parameter vector of the individual
     * \param[in] changed           a reference to a vector which has true for parameters that change and false otherwise
     * \param[in] prevEval          pointer to a double storage with the previous \f$ \chi^2 \f$ for training set
     * \param[in] pp                index of inner loop over number of parameters
     * \param[in] iter              current iteration number
     * \param[in] beta0             pointer to beta for annealing
     * \param[in] nParam            number of parameters in the model
     */
    void stepMutation(      ACMIndividual          *ind,
                            std::vector<double>    *param,
                            std::vector<bool>      *changed,
                            double                 *prevEval,
                      const size_t                  pp,
                      const int                     iter,
                            double                 *beta0,
                      const size_t                  nParam);
   void mutateOld(      ga::Individual   *individual,
                        const double            prMut);

public:

    /*!
     * \brief Constructor
     * \param[in] logfile   pointer to log file (may be nullptr)
     * \param[in] verbose   Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
     * \param[in] bch       pointer to BayesConfigHandler object
     * \param[in] fitComp   pointer to ACMFitnessComputer object
     * \param[in] sii       pointer to SharedIndividualInfo object
     * \param[in] nParam    size of the force field parameter vector
     */
    MCMCMutator(      FILE                     *logfile,
                const bool                      verbose,
                      BayesConfigHandler       *bch,
                      ACMFitnessComputer       *fitComp,
                      SharedIndividualInfo     *sii,
                const size_t                    nParam)
    : Mutator(), gen(rd()), dis(std::uniform_int_distribution<size_t>(0, nParam-1))
    {
        if (bch->seed() == 0)
        {
            gen.seed(::time(NULL));
        }
        else
        {
            gen.seed(bch->seed());
        }
        logfile_ = logfile;
        verbose_ = verbose;
        bch_     = bch;
        fitComp_ = fitComp;
        sii_     = sii;
    };

 
    /*!
     * \brief Run the Markov chain Monte carlo (MCMC) simulation
     * \param[in] ind               pointer to the individual
     * \param[in] evaluate_testset  If true, evaluate the energy on
     *                              the test set.
     */
    virtual void mutate(ga::Individual *ind,
                        double          prMut);

    //! \return the number of calls to the objective function MCMCMutator::MCMC() routine
    size_t numberObjectiveFunctionCalls() const
    {
        return 1 + bch_->maxIter() * sii_->nParam();
    }

    //! \return whether a minimum was found
    bool foundMinimum() const { return bMinimum_; }

    /*!
     * \brief Print the MC statistics to a file.
     * \param[in] ind   pointer to the indivifdual
     * \param[in] fp    File pointer to print to
     */
    void printMonteCarloStatistics(ACMIndividual   *ind,
                                   FILE            *fp);

    /*!
     * \brief Perform a sensitivity analysis by systematically changing all parameters and re-evaluating the \f$ \chi^2 \f$.
     * \param[in] ind   pointer to individual
     * \param[in] ims   Dataset to perform sensitivity analysis on
     */
    void sensitivityAnalysis(ACMIndividual  *ind,
                             iMolSelect      ims);

};


} //namespace alexandria


#endif //ALEXANDRIA_MCMCMUTATOR_H