#ifndef ALEXANDRIA_SHAREDINDIVIDUALINFO_H
#define ALEXANDRIA_SHAREDINDIVIDUALINFO_H


#include <vector>
#include <string>
#include <map>

#include "mutability.h"
#include "molgen.h"


namespace alexandria
{


class SharedIndividualInfo
{

private:

    //! Communication record
    t_commrec *cr_;
    //! Base Poldata from which to make copies
    Poldata pd_;
    //! Base targets_ from which to make copies
    std::map<iMolSelect, std::map<eRMS, FittingTarget>> targets_;
    //! Training steps per parameter
    std::vector<int> ntrain_;
    //! Lower bound per parameter
    std::vector<double> lowerBound_;
    //! Upper bound per parameter
    std::vector<double> upperBound_;
    //! Mutability per parameter
    std::vector<Mutability> mutability_; 
    //! Weighted temperature
    std::vector<double> weightedTemperature_;
    //! Optimization index for each parameter
    std::vector<OptimizationIndex> optIndex_;
    //! Class-name for each parameter
    std::vector<std::string> paramNames_;
    //! Classes in the system
    std::vector<std::string> paramClass_;
    //! Class-index for each parameter
    std::vector<int> paramClassIndex_;  // FIXME: shouldn't this be a vector of size_t
    //! Base name for parameter convergence files
    std::string xvgconv_;
    //! Base name for Chi2 convergence file
    std::string xvgepot_;
    //! Base name for Force field output file
    std::string outputFile_;

    //! \brief Fills the \p targets_ map
    void fillFittingTargets();

public:

    /*!
     * Constructor
     * \param[in] cr The communications record
     */
    SharedIndividualInfo(t_commrec *cr)
    : cr_(cr)
    {
        fillFittingTargets();
    }

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Poldata stuff                     *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Fill the Poldata attribute from a file
     * \param[in] fp        File pointer for printing information
     * \param[in] pd_fn     name of the gentop file
     */
    void fillPoldata(      FILE *fp,
                     const char *pd_fn);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Poldata stuff                       *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: FittingTarget stuff               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Propage weights from training set to other sets
    void propagateWeightFittingTargets();

    // FIXME: How can we avoid having to decleare the query functions both in SharedIndividualInfo
    // and in Individual, we put them around the FittingTarget class?

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
    * END: FittingTarget stuff                 *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Weighted temperature stuff        *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * Compute weighted temperature for each parameter
     * If weighting is not required, it will just be a vector of 1s
     * @param tempWeight    whether weighting is required
     */
    void computeWeightedTemperature(const bool tempWeight);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Weighted temperature stuff          *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: OptimizationIndex stuff           *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief Generate \p optIndex_
     * \param[in] fp File to print logging information to. May be nullptr.
     * \param[in] mg MolGen pointer
     */
    void generateOptimizationIndex(FILE   *fp,
                                   MolGen *mg);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: OptimizationIndex stuff             *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Vector stuff                      *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! Fills \p ntrain_ \p lowerBound_ \p upperBound_ \p mutability_ and \p paramNames_
     * @param mindata   mininum number of datapoints to consider a parameter
     */
    void fillVectors(const int mindata);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Vector stuff                      *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: ParamClassIndex stuff             *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * Assign a class (by index) to each parameter
     */
    void assignParamClassIndex();

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: ParamClassIndex stuff               *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: File stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief Set the output file names.
    *
    * The parameter values are split over
    * a number of files in order to make it easier to visualize the
    * results. The parameter classes should therefore match the
    * parameter names. E.g. a class could be alpha, another zeta.
    *
    * \param[in] xvgconv    The parameter convergence base name
    * \param[in] paramClass The parameter classes (e.g. zeta, alpha)
    * \param[in] xvgepot    The base filename to print the chi2 value
    */
    void setOutputFiles(const char                     *xvgconv,
                        const std::vector<std::string> &paramClass,
                        const char                     *xvgepot);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: File stuff                          *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Get a constant \p optIndex_ reference
    const std::vector<OptimizationIndex> &optIndex() const { return optIndex_; }

    //! \brief Get a pointer to \p optIndex_
    std::vector<OptimizationIndex> *optIndexPtr() { return &optIndex_; }
    
    //! \brief Get a constant \p paramNames_ reference
    const std::vector<std::string> &paramNames() const { return paramNames_; }

    //! \brief Return the number of parameters
    size_t nParam() const { return paramNames_.size(); }

    //! \brief Get a constant \p paramClass_ reference
    const std::vector<std::string> &paramClass() const { return paramClass_; }

    //! \brief Get a constant \p paramClassIndex_ reference
    const std::vector<int> &paramClassIndex() const { return paramClassIndex_; }

    /*! \brief
     * Returns the current vector of lower bounds
     * @return the current vector of lower bounds
     */
    const std::vector<double> &lowerBound() const { return lowerBound_; }

    /*! \brief
     * Returns the current vector of upper bounds
     * @return the current vector of upper bounds
     */
    const std::vector<double> &upperBound() const { return upperBound_; }

    /*! \brief
     * Returns the number of training points per parameter
     */
    const std::vector<int> &nTrain() const { return ntrain_; }

    //! \brief Get a constant \p mutability_ reference
    const std::vector<Mutability> &mutability() const { return mutability_; }

    //! \brief Return xvg file for convergence information
    const std::string &xvgConv() const { return xvgconv_; }

    //! \brief Return xvg file for epot information
    const std::string &xvgEpot() const { return xvgepot_; }

    //! \brief Get a constant \p outputFile_ reference
    const std::string &outputFile() const { return outputFile_; }

    //! \brief Get a constant \p weightedTemperature_ reference
    const std::vector<double> &weightedTemperature() const { return weightedTemperature_; }

    //! \brief Get a constant \p targets_ reference
    const std::map<iMolSelect, std::map<eRMS, FittingTarget>> &targets() const { return targets_; }

    //! \return pointer to \p targets_
    std::map<iMolSelect, std::map<eRMS, FittingTarget>> *targetsPtr() { return &targets_; }
    
    //! \brief Return the poldata as pointer to const variable
    const Poldata *poldata() const { return &pd_; }

    //! \brief Return the poldata as const reference
    const Poldata &poldataConst() const { return pd_; }
    
    //! \brief Return pointer to the poldata
    Poldata *poldata() { return &pd_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace alexandria


#endif //ALEXANDRIA_SHAREDINDIVIDUALINFO_H