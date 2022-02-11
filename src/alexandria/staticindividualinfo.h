/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef ALEXANDRIA_STATICINDIVIDUALINFO_H
#define ALEXANDRIA_STATICINDIVIDUALINFO_H

#include <vector>
#include <string>
#include <map>
#include "act/ga//Genome.h"

#include "act/utility/communicationrecord.h"
#include "act/basics/mutability.h"
#include "molgen.h"

namespace alexandria
{
/*!
 * \brief Contains all information that is shared among ACMIndividual objects and other classes that manage them
 */
class StaticIndividualInfo
{

private:
    //! My identy
    int                                                  id_ = -1;
    //! Communication record
    const CommunicationRecord                           *cr_;
    //! Poldata for this individual
    Poldata                                              pd_;
    //! Base targets_ from which to make copies
    std::map<iMolSelect, std::map<eRMS, FittingTarget>>  targets_;
    //! Default parameter values as specified by input file
    std::vector<double>                                  defaultParam_;
    //! Training datapoints per parameter
    std::vector<int>                                     ntrain_;
    //! Lower bound per parameter
    std::vector<double>                                  lowerBound_;
    //! Upper bound per parameter
    std::vector<double>                                  upperBound_;
    //! Mutability per parameter
    std::vector<Mutability>                              mutability_; 
    //! Weighted temperature
    std::vector<double>                                  weightedTemperature_;
    //! Optimization index for each parameter
    std::vector<OptimizationIndex>                       optIndex_;
    //! Class-name for each parameter
    std::vector<std::string>                             paramNames_;
    //! Classes in the system
    std::vector<std::string>                             paramClass_;
    //! Class-index for each parameter
    std::vector<size_t>                                  paramClassIndex_;
    //! Base name for parameter convergence files
    std::string                                          xvgconv_;
    //! Base name for Chi2 convergence file
    std::string                                          xvgepot_;
    //! Base name for Force field output file
    std::string                                          outputFile_;
    //! Output file prefix
    std::string                                          prefix_;

    //! \brief Fills the fitting targets data structure
    void fillFittingTargets();

public:

    /*!
     * Constructor
     * \param[in] cr The communication record
     */
    StaticIndividualInfo(const CommunicationRecord *cr);

    /*! \brief Fill the \p id_ and \p prefix_ fields
     * (to be called after the communication record structure has been initialized)
     */
    void fillIdAndPrefix();
    
    
    /*! \brief Set the output file name
     * \param[in] outputFile The base force field file name for output
     */
    void setOutputFile(const std::string &outputFile);

    /*! \brief Set the output file name for best output
     * \param[in] outputFile The base force field file name for output
     */
    void setFinalOutputFile(const std::string &outputFile);
    
    /*!
     * \brief Set the \f$ \chi^2 \f$ to 0 for a given data set.
     * \param[in] ims The data set to reset
     */
    void resetChiSquared(iMolSelect ims);
    
    /*!
     * \brief Add up the \f$ \chi^2 \f$ components
     * Sum over the energies of the cores if desired.
     * Also multiplies the terms by the weighting factors.
     * \param[in] cr        Pointer to communication record
     * \param[in] parallel  Whether or not to sum in parallel
     * \param[in] ims       The selection dataset to sum
     */
    void sumChiSquared(bool             parallel,
                       iMolSelect       ims);

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Poldata stuff                     *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Fill the Poldata attribute by reading from a file
     * \param[in] fp        File pointer for printing information
     * \param[in] pd_fn     name of the gentop (Force Field) file
     */
    void fillPoldata(      FILE *fp,
                     const char *pd_fn);

    /*!
     * \brief Copy the Force Field parameters to the Poldata structure
     * \param[in] changed List over the parameters that have changed.
     * \param[in] genome  The parameters in the current genome
     */
    void updatePoldata(const std::vector<bool> &changed,
                       const ga::Genome        *genome);
    
    /*! \brief Save the current state of the Force Field to the output file
     * \param[in] updateCheckSum If true, the checksum is updated, typically
     *                           this should only be done at the end of a run
     *                           since it is expensive.
     */
    void saveState(bool updateCheckSum);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Poldata stuff                       *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return the communication record
    const CommunicationRecord *commRec() const { return cr_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: FittingTarget stuff               *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief In the  fitting targets structure, 
     * propagate the \f$ \chi^2 \f$ component weights from 
     * the training set to the other sets.
     */
    void propagateWeightFittingTargets();

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
    * END: FittingTarget stuff                 *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Weighted temperature stuff        *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Compute weighted temperature for each parameter (for MCMC optimization)
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

    /*!
     * \brief Generate the vector of OptimizationIndex instances
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

    /*!
     * \brief Fills some vector structures for parameter information.
     * Per parameter, we have: its default value, required amount of training steps, lower bound, upper bound, Mutability, and name<br>
     * Dev: Fills \p defaultParam_ \p ntrain_ \p lowerBound_ \p upperBound_ \p mutability_ and \p paramNames_
     * \param[in] mindata mininum number of existing datapoints to consider a parameter for optimization
     */
    void fillVectors(const int mindata);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Vector stuff                      *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: ParamClassIndex stuff             *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Assign a class index to each parameter
    void assignParamClassIndex();

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: ParamClassIndex stuff               *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: File stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
    * \brief Set the output file names.
    *
    * The parameter values are split over
    * a number of files in order to make it easier to visualize the
    * results. The parameter classes should therefore match the
    * parameter names. E.g. a class could be alpha, another zeta.
    *
    * \param[in] xvgconv    The parameter convergence base name
    * \param[in] paramClass The parameter classes (e.g. zeta, alpha)
    * \param[in] xvgepot    The base filename to print the \f$ \chi^2 \f$ value
    */
    void setOutputFiles(const char                     *xvgconv,
                        const std::vector<std::string> &paramClass,
                        const char                     *xvgepot);

    //! \return the output file prefix
    const std::string &prefix() const { return prefix_; }
    
    //! \return my id
    int id() const { return id_; }
    
    /* * * * * * * * * * * * * * * * * * * * * *
    * END: File stuff                          *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return the vector of OptimizationIndex instances as a const reference
    const std::vector<OptimizationIndex> &optIndex() const { return optIndex_; }

    //! \return a pointer to the vector of OptimizationIndex instances
    std::vector<OptimizationIndex> *optIndexPtr() { return &optIndex_; }
    
    //! \return the vector of default parameter values as const reference
    const std::vector<double> &defaultParam() const { return defaultParam_; }

    //! \return the vector of parameter names as a const reference
    const std::vector<std::string> &paramNames() const { return paramNames_; }

    //! \return the number of parameters to optimize
    size_t nParam() const { return paramNames_.size(); }

    //! \return the vector of existing parameter classes as const reference
    const std::vector<std::string> &paramClass() const { return paramClass_; }

    //! \return a vector which links parameters to a class via the class index
    const std::vector<size_t> &paramClassIndex() const { return paramClassIndex_; }

    //! \return the vector of lower bounds as a const reference
    const std::vector<double> &lowerBound() const { return lowerBound_; }

    /*!
     * Get the lower bound of a given parameter
     * \param[in] i the index of the parameter
     * \return the lower bound of the parameter at index \p i
     */
    double lowerBoundAtIndex(const size_t i) const { return lowerBound_[i]; }

    //! \return the vector of upper bounds as a const reference
    const std::vector<double> &upperBound() const { return upperBound_; }

    /*!
     * Get the upper bound of a given parameter
     * \param[in] i the index of the parameter
     * \return the upper bound of the parameter at index \p i
     */
    double upperBoundAtIndex(const size_t i) const { return upperBound_[i]; }

    //! \return the vector of training datapoints as a const reference
    const std::vector<int> &nTrain() const { return ntrain_; }

    //! \return the vector of Mutability instances as a const reference
    const std::vector<Mutability> &mutability() const { return mutability_; }

    //! \return the base name for parameter convergence files as a const reference
    const std::string &xvgConv() const { return xvgconv_; }

    //! \return the base name for the \f$ \chi_2 \f$ convergence file as a const reference
    const std::string &xvgEpot() const { return xvgepot_; }

    //! \return the base name for Force Field output files as a const reference
    const std::string &outputFile() const { return outputFile_; }

    //! \return the vector of weighted temperatures as a const reference
    const std::vector<double> &weightedTemperature() const { return weightedTemperature_; }

    //! \return the fitting targets data structure as a const reference
    const std::map<iMolSelect, std::map<eRMS, FittingTarget>> &targets() const { return targets_; }

    //! \return a pointer to the fitting targets data structure
    std::map<iMolSelect, std::map<eRMS, FittingTarget>> *targetsPtr() { return &targets_; }
    
    //! \return a pointer to the const Poldata structure (for const objects)
    const Poldata *poldata() const { return &pd_; }

    //! \return the Poldata structure as a const reference
    const Poldata &poldataConst() const { return pd_; }
    
    //! \return a pointer to the Poldata structure
    Poldata *poldata() { return &pd_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace alexandria


#endif //ALEXANDRIA_STATICINDIVIDUALINFO_H
