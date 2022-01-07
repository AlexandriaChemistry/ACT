/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
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

/*!
 * \brief An individual for the Alexandria Charge Model (ACM)
 * It has its own force field parameters and handles its output files.
 */
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
    FILE *fpe_ = nullptr;
    //! Force field file name (to be used in saveState())
    std::string outputFile_;

public:
    /*!
     * \brief Property constructor
     * \param[in] id            the ID of the individual
     * \param[in] sii           pointer to SharedIndividualInfo instance
     * \param[in] outputFile    the base name for Force Field output files
     */
    ACMIndividual(const int                     id,
                        SharedIndividualInfo   *sii,
                  const std::string            &outputFile)
    : ga::Individual()
    {
        id_ = id;
        sii_ = sii;
        outputFile_ = "ind" + std::to_string(id_) + "/ind" + std::to_string(id_) + "-" + outputFile;

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
    * BEGIN: Cloning                           *
    * * * * * * * * * * * * * * * * * * * * * */

    virtual void copyGenome(Individual *other) { param_ = static_cast<ACMIndividual*>(other)->param(); }

    virtual ga::Individual *clone() { return new ACMIndividual(*this); }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Cloning                             *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Adding parameters                 *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Add a force field parameter
     * \param[in] val the value of the parameter
     */
    void addParam(const real val) { initialParam_.push_back(val); param_.push_back(val); }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Adding parameters                   *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: File stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    virtual void fprintSelf(FILE *fp);

    /*!
     * \brief Print individual header to a file
     * For individual with ID X, the header is: \nIndividual X\n
     * \param[in] fp the file to print to
     */
    void printHeader(FILE *fp);

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
    
    //! Close \f$ \chi^2 \f$ and parameter convergence files
    void closeConvFiles();

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: File stuff                          *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Chi2 stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Set the \f$ \chi^2 \f$ to 0 for a given data set.
     * \param[in] ims The data set to reset
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

    /*!
     * \brief Add up the \f$ \chi^2 \f$ components
     * Sum over the energies of the cores if desired.
     * Also multiplies the terms by the weighting factors.
     * \param[in] cr        Pointer to communication record
     * \param[in] parallel  Whether or not to sum in parallel
     * \param[in] ims       The selection dataset to sum
     */
    void sumChiSquared(t_commrec *cr, bool parallel, iMolSelect ims);

    /*!
     * \brief Print the \f$ \chi^2 \f$ components.
     * \param[in] fp  File pointer to print to, may be nullptr
     * \param[in] ims The selection dataset to print
     */  
    void printChiSquared(t_commrec *cr, FILE *fp, iMolSelect ims) const;

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Chi2 stuff                          *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Output stuff                      *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Save the current state of the Force Field to the output file
    void saveState();

    /*!
     * \brief Print the Force Field parameters to a file
     * \param[in] fp File pointer to open file
     */
    void printParameters(FILE *fp) const;

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Output stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Poldata stuff                     *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Copy the Force Field parameters to the Poldata structure
     * \param[in] changed List over the parameters that have changed.
     */
    void toPoldata(const std::vector<bool> &changed);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Poldata stuff                       *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: FittingTarget queries             *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Return the fitting targets for editing for a given dataset
     * \param[in] ims The selection dataset to return
     * \return The map of eRMS to FittingTarget, or nullptr
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

    /*!
     * \brief Return the fitting targets for a given dataset as constant
     * \param[in] ims The selection dataset to return
     * \return The map of eRMS to FittingTarget or nullptr
     */
    const std::map<eRMS, FittingTarget> &fittingTargetsConst(iMolSelect ims) const
    {
        auto tt = targets_.find(ims);
        GMX_RELEASE_ASSERT(targets_.end() != tt, gmx::formatString("Cannot find selection %s", iMolSelectName(ims)).c_str());
        return tt->second;
    }

    /*!
     * \brief get a pointer to a FittingTarget given a dataset choice and \f$ \chi^2 \f$ component
     * \param[in] ims The selection dataset
     * \param[in] rms The \f$ \chi^2 \f$ component
     * \return FittingTarget pointer or nullptr if not found
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

    //! \return the ID
    int id() const { return id_; }

    /*!
     * \brief Set a new value for the ID
     * Also adjusts the Force Field output file.<br>
     * CAREFUL! This does not change convergence files! You will have to call the open file routines again to fix it!<br>
     * That is, ACMIndividual::openParamConvFiles() and ACMIndividual::openChi2ConvFile()
     * \param[in] id the new ID
     */
    void setId(const int id)
    {
        id_ = id;
        const size_t firstIndex = outputFile_.find("-");  // We need to discard the existing indX- part
        const size_t strLength = outputFile_.size();
        outputFile_ = "ind" + std::to_string(id_) + "/ind" + std::to_string(id_) + outputFile_.substr(firstIndex, strLength - firstIndex);
    }

    //! \return a pointer to the SharedIndividualInfo instance
    SharedIndividualInfo *sii() { return sii_; }

    //! \return a pointer to the SharedIndividualInfo instance (for const objects)
    SharedIndividualInfo *siiConst() const { return sii_; }

    //! \return a constant reference of the fitting targets
    const std::map<iMolSelect, std::map<eRMS, FittingTarget>> &targets() const { return targets_; }

    //! \return the Poldata structure as pointr to a const variable
    const Poldata *poldata() const { return &pd_; }

    //! \return the Poldata structure as const reference
    const Poldata &poldataConst() const { return pd_; }
    
    //! \return pointer to the Poldata structure
    Poldata *poldata() { return &pd_; }

    //! \return the number of Force Field parameters
    size_t nParam() const { return param_.size(); }

    /*!
     * \brief Set parameter at a given index to a new value
     * \param[in] j   the index
     * \param[in] val the new value
     */
    void setParam(const size_t j, const real val)
    {
        GMX_RELEASE_ASSERT(j < param_.size(), "Parameter out of range");
        param_[j] = val;
    }

    /*!
     * \brief Get the value of a parameter by index
     * \param[in] j the index
     * \return the value of the parameter at index \p j
     */
    double paramAtIndex(const size_t j) const { return param_[j]; }

    /*!
     * \brief Set all parameters to new values
     * \param[in] param the new values
     */
    void setParam(std::vector<double> param)
    {
        GMX_RELEASE_ASSERT(param.size() == param_.size() || param_.empty(),
                           "Incorrect size of input parameters");
        param_ = param;
    }

    //! \return the initial vector of parameters as a const reference
    const std::vector<double> &initialParam() const { return initialParam_; }

    //! \return the current vector of parameters as a const reference
    const std::vector<double> &param() const { return param_; }

    //! \return a pointer to the current vector of parameters
    std::vector<double> *paramPtr() { return &param_; }

    //! \return the vector of best parameters as a const reference
    const std::vector<double> &bestParam() const { return bestParam_; }

    /*!
     * \brief Set a new best parameter vector
     * \param[in] param the new best parameter vector
     */
    void setBestParam(const std::vector<double> &param) { bestParam_ = param; }

    //! \return the vector of mean value calculated for each parameter as const reference
    const std::vector<double> &pMean() const { return pmean_; }

    //! \return a pointer to the vector of mean value calculated for each parameter
    std::vector<double> *pMeanPtr() { return &pmean_; }

    //! the vector of standard deviation calculated for each parameter as const reference
    const std::vector<double> &pSigma() const { return psigma_; }

    //! \return a pointer to vector of standard deviation calculated for each parameter
    std::vector<double> *pSigmaPtr() { return &psigma_; }

    //! \return the vector of number of attempted moves for each parameter as const reference
    const std::vector<int> &attemptedMoves() const { return attemptedMoves_; }

    //! \return a pointer to the vector of number of attempted moves for each parameter
    std::vector<int> *attemptedMovesPtr() { return &attemptedMoves_; }

    //! \return the vector of number of accepted moves for each parameter as const reference
    const std::vector<int> &acceptedMoves() const { return acceptedMoves_; }

    //! \return a pointer to the vector of number of accepted moves for each parameter
    std::vector<int> *acceptedMovesPtr() { return &acceptedMoves_; }

    //! \return the vector of force field parameter convergence files as const reference
    const std::vector<FILE*> &fpc() const { return fpc_; }

    //! \return a pointer to the \f$ \chi^2 \f$ convergence file
    FILE *fpe() { return fpe_; }

    //! \return a pointer to the \f$ \chi^2 \f$ convergence file (for const objects)
    FILE *fpeConst() const { return fpe_; }

    //! \return the name of the Force Field output file as const reference
    const std::string &outputFile() const { return outputFile_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace alexandria


#endif //ALEXANDRIA_ACMINDIVIDUAL_H
