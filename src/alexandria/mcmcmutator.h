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
    //! Whether we are in verbose mode
    bool verbose_;

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

public:

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

};


} //namespace alexandria


#endif //ALEXANDRIA_MCMCMUTATOR_H