/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julianramon.marradesfurquet.8049@student.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef ALEXANDRIA_ACMINDIVIDUAL_H
#define ALEXANDRIA_ACMINDIVIDUAL_H


#include <cstdio>
#include <map>
#include <vector>
#include <string>

#include "molgen.h"
#include "poldata.h"
#include "sharedindividualinfo.h"

#include "ga/Individual.h"


namespace alexandria
{


class ACMIndividual : public ga::Individual
{

private:

    //! ID of the individual
    int id_;
    //! Pointer to shared individual information
    SharedIndividualInfo *sii_;
    //! Fitting targets for each dataset and eRMS
    std::map<iMolSelect, std::map<eRMS, FittingTarget>> targets_;
    //! Force field data structure
    Poldata pd_;
    //! Initial parameter vector
    std::vector<double> initialParam_;
    //! Parameter vector
    std::vector<double> param_;
    //! Best parameter vector
    std::vector<double> bestParam_;
    //! Mean of each parameter
    std::vector<double> pmean_;
    //! Standard deviation of each parameter
    std::vector<double> psigma_;
    //! Attempted changes for each parameter
    std::vector<int> attemptedMoves_;
    //! Accepted changes for each parameter
    std::vector<int> acceptedMoves_;
    //! Convergence file for each parameter type
    std::vector<FILE*> fpc_;
    //! Convergence file for Chi2
    FILE *fpe_;
    //! Force field file name (to be used in saveState())
    std::string outputFile_;

public:

    ACMIndividual(const int                     id,
                        SharedIndividualInfo   *sii,
                  const std::string            &outputFile)
    : ga::Individual()
    {
        id_ = id;
        sii_ = sii;
        outputFile_ = "ind" + std::to_string(id_) + "-" + outputFile;

        // Initialize vectors for statistics and bestParam_.
        // initialParam_ and param_ will be initialized later
        size_t nParam = sii_->nParam();
        pmean_.resize(nParam, 0.0);
        psigma_.resize(nParam, 0.0);
        attemptedMoves_.resize(nParam, 0);
        acceptedMoves_.resize(nParam, 0);
        bestParam_.resize(nParam, 0.0);

        // Copy targets_ from sii_
        targets_ = sii_->targets();  // This should make a deep copy if
                                     // the copy constructors are well made

        // Copy poldata from sii_
        pd_ = sii_->poldataConst();
        
    }

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Adding parameters                 *
    * * * * * * * * * * * * * * * * * * * * * */

    void addParam(const real val) { initialParam_.push_back(val); param_.push_back(val); }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Adding parameters                   *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: File stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * Open parameter convergence surveillance files
     */
    void openParamConvFiles(const gmx_output_env_t *oenv);

    /*!
     * Open a chi2 surveillance file
     * @param bEvaluate_testset     whether the test set will be evaluated
     */
    void openChi2ConvFile(const gmx_output_env_t    *oenv,
                          const bool                 bEvaluate_testset);
    
    /*!
     * Close chi2 and parameter convergence files
     */
    void closeConvFiles();

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: File stuff                          *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Chi2 stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief Set the chiSquared to zero.
     * \param[in] ims The selection to reset
     */
    void resetChiSquared(iMolSelect ims)
    {
        auto fts = fittingTargets(ims);
        if (fts != nullptr)
        {
            for (auto &ft : *fts)
            {
                ft.second.reset();
            }
        }
    }

    /*! \brief 
     * Sum over the energies of the cores if desired.
     * Also multiplies the terms by the weighting factors.
     * \param[in] cr        Pointer to communication record
     * \param[in] parallel  Whether or not to sum in parallel
     * \param[in] ims       The selection to sum
     */
    void sumChiSquared(t_commrec *cr, bool parallel, iMolSelect ims);

    /*! \brief Print the chiSquared components.
     * \param[in] fp  File pointer to print to, may be nullptr
     * \param[in] ims The selection to print
     */  
    void printChiSquared(t_commrec *cr, FILE *fp, iMolSelect ims) const;

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Chi2 stuff                          *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Output stuff                      *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Save the current state of the Force field
    void saveState();

    /*! \brief
     * Print the parameters to a file
     * \param[in] fp File pointer to open file
     */
    void printParameters(FILE *fp) const;

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Output stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Poldata stuff                     *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief
     * Copy the optimization parameters to the poldata structure
     * \param[in] changed List over the parameters that have changed.
     */
    void toPoldata(const std::vector<bool> &changed);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Poldata stuff                       *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: FittingTarget queries             *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief Return the fitting targets for editing
     * \param[in] ims The selection to return
     * \return The map of fittingtargets or nullptr
     */
    std::map<eRMS, FittingTarget> *fittingTargets(iMolSelect ims)
    {
        auto tt = targets_.find(ims);
        if (targets_.end() == tt)
        {
            return nullptr;
        }
        else
        {
            return &tt->second;
        }
    }

    /*! \brief Return the fitting targets as constant
     * \param[in] ims The selection to return
     * \return The map of fittingtargets or nullptr
     */
    const std::map<eRMS, FittingTarget> &fittingTargetsConst(iMolSelect ims) const
    {
        auto tt = targets_.find(ims);
        GMX_RELEASE_ASSERT(targets_.end() != tt, gmx::formatString("Cannot find selection %s", iMolSelectName(ims)).c_str());
        return tt->second;
    }

    /*! \brief return appropriate fitting target
     * \param[in] ims The selection
     * \param[in] rms The contributor to the chi squared
     * \return FittingTarget class or nullptr if not found
     */
    FittingTarget *target(iMolSelect ims, eRMS rms)
    {
        auto itarget = targets_.find(ims);
        if (itarget == targets_.end())
        {
            return nullptr;
        }
        auto ift = itarget->second.find(rms);
        if (ift == itarget->second.end())
        {
            return nullptr;
        }
        return &ift->second;
    }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: FittingTarget queries               *
    * * * * * * * * * * * * * * * * * * * * * */
    
    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and Setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * Get the ID of the individual
     * @returns the id
     */
    int id() const { return id_; }

    //! \brief Return the poldata as pointe to const variable
    const Poldata *poldata() const { return &pd_; }

    //! \brief Return the poldata as const reference
    const Poldata &poldataConst() const { return pd_; }
    
    //! \brief Return pointer to the poldata
    Poldata *poldata() { return &pd_; }

    //! \brief Return the number of parameters
    size_t nParam() const { return param_.size(); }

    /*! \brief
     * Set parameter j to a new value
     * \param[j]   Index
     * \param[val] The new value
     */
    void setParam(size_t j, real val)
    {
        GMX_RELEASE_ASSERT(j < param_.size(), "Parameter out of range");
        param_[j] = val;
    }

    /*! \brief
     * Set all parameters to the array passed
     */
    void setParam(std::vector<double> param)
    {
        GMX_RELEASE_ASSERT(param.size() == param_.size() || param_.empty(),
                           "Incorrect size of input parameters");
        param_ = param;
    }

    /*! \brief
     * Returns the current vector of parameters.
     */
    const std::vector<double> &initialParam() const { return initialParam_; }

    /*! \brief
     * Returns the current vector of parameters.
     */
    const std::vector<double> &param() const { return param_; }

    //! \return pointer to \p param_
    std::vector<double> *paramPtr() { return &param_; }

    /*! \brief
     * Returns the vector of best found value for each parameter.
     */
    const std::vector<double> &bestParam() const { return bestParam_; }

    /*!
     * Set a new best parameter vector
     * \param[in] param the new best parameter vector
     */
    void setBestParam(const std::vector<double> &param) { bestParam_ = param; }

    /*! \brief
     * Returns the vector of mean value calculated for each parameter.
     */
    const std::vector<double> &pMean() const { return pmean_; }

    //! \return a pointer to \p pmean_
    std::vector<double> *pMeanPtr() { return &pmean_; }

    /*! \brief
     * Returns the vector of standard deviation calculated for each parameter.
     */
    const std::vector<double> &pSigma() const { return psigma_; }

    //! \return a pointer to \p psigma_
    std::vector<double> *pSigmaPtr() { return &psigma_; }

    /*! \brief
     * Return the vector of number of attempted moves for each parameter
     */
    const std::vector<int> &attemptedMoves() const { return attemptedMoves_; }

    //! \return a pointer to \p attemptedMoves_
    std::vector<int> *attemptedMovesPtr() { return &attemptedMoves_; }

    /*! \brief
     * Return the vector of number of accepted moves for each parameter
     */
    const std::vector<int> &acceptedMoves() const { return acceptedMoves_; }

    //! \return a pointer to \p acceptedMoves_
    std::vector<int> *acceptedMovesPtr() { return &acceptedMoves_; }

    //! \return a constant \p fpc_ reference
    const std::vector<FILE*> &fpc() const { return fpc_; }

    //! \return a \p fpe_ file for Chi2
    FILE *fpe() { return fpe_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace alexandria


#endif //ALEXANDRIA_ACMINDIVIDUAL_H