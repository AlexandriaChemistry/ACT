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
    SharedIndividualInfo sii_;
    //! Fitness for training dataset
    double fitnessTrain_;
    //! Fitness for test set
    double fitnessTest_;
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

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Output stuff                      *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Save the current state of the Force field
    void saveState();

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

    /*!
     * Get the fitness of the individual in training set
     * @returns the fitness
     */
    double fitnessTrain() const { return fitnessTrain_; }

    /*!
     * Get a pointer to \p fitnessTrain_
     * @returns the pointer
     */
    double *fitnessTrainPtr() { return &fitnessTrain_; }

    /*!
     * Set the fitness of the individual in training set
     * @param fitnessTrain the fitness
     */
    void setFitnessTrain(const double fitnessTrain) { fitnessTrain_ = fitnessTrain; }

    /*!
     * Get the fitness of the individual in test set
     * @returns the fitness
     */
    double fitnessTest() const { return fitnessTest_; }

    /*!
     * Get a pointer to \p fitnessTest_
     * @returns the pointer
     */
    double *fitnessTestPtr() { return &fitnessTest_; }

    /*!
     * Set the fitness of the individual in test set
     * @param fitnessTest the fitness
     */
    void setFitnessTest(const double fitnessTest) { fitnessTest_ = fitnessTest; }

    //! \brief Return the poldata as const variable
    const Poldata *poldata() const { return &pd_; }
    
    //! \brief Return the poldata
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

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace alexandria


#endif //ALEXANDRIA_ACMINDIVIDUAL_H