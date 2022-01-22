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
#include "staticindividualinfo.h"

#include "ga/Genome.h"
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
    int                   id_ = 0;
    //! Pointer to static individual information
    StaticIndividualInfo *sii_ = nullptr;
    //! Initial genome
    ga::Genome            initialGenome_;
    //! Parameter vector
    ga::Genome            genome_;
    //! Best parameter vector
    ga::Genome            bestGenome_;
    //! Force field file name (to be used in saveState())
    //std::string           outputFile_;

public:
    /*!
     * \brief Property constructor
     * \param[in] id            the ID of the individual
     * \param[in] sii           pointer to StaticIndividualInfo instance
     * \param[in] outputFile    the base name for Force Field output files
     */
    ACMIndividual(const int                     id,
                        StaticIndividualInfo   *sii,
                  const std::string            &outputFile);

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Cloning                           *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief Copy the genome to this individual
     * \param[in] genome Complete genome, only the bases are copied.
     */
    virtual void copyGenome(ga::Genome *genome);

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
    void addParam(const real val);

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Adding parameters                   *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Output stuff                      *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Print the Force Field parameters to a file
     * \param[in] fp File pointer to open file
     */
    void printParameters(FILE *fp) const;

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Output stuff                        *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: FittingTarget queries             *
    * * * * * * * * * * * * * * * * * * * * * */

#ifdef OLD
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
#endif
    
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
    //void setId(const int id)
    //{
    //  id_ = id;
    //  const size_t firstIndex = outputFile_.find("-");  // We need to discard the existing indX- part
    //  const size_t strLength = outputFile_.size();
    //  outputFile_ = "ind" + std::to_string(id_) + "/ind" + std::to_string(id_) + outputFile_.substr(firstIndex, strLength - firstIndex);
    //}

    //! \return a pointer to the StaticIndividualInfo instance
    StaticIndividualInfo *sii() { return sii_; }

    //! \return a pointer to the StaticIndividualInfo instance (for const objects)
    StaticIndividualInfo *siiConst() const { return sii_; }

    //! \return a constant reference of the fitting targets
    // const std::map<iMolSelect, std::map<eRMS, FittingTarget>> &targets() const { return targets_; }

    //! \return the Poldata structure as pointr to a const variable
    //const Poldata *poldata() const { return &sii_.pd_; }

    //! \return the Poldata structure as const reference
    //const Poldata &poldataConst() const { return pd_; }
    
    //! \return pointer to the Poldata structure
    //Poldata *poldata() { return &pd_; }

    //! \return the number of Force Field parameters
    //size_t nParam() const { return genome_.nBase(); }

    /*!
     * \brief Set parameter at a given index to a new value
     * \param[in] j   the index
     * \param[in] val the new value
     */
    //void setParam(const size_t j, const real val)
    //{
    //   GMX_RELEASE_ASSERT(j < genome_.nParam(), "Parameter out of range");
    //  genome_setParam(j, val);
    //}

    /*!
     * \brief Get the value of a parameter by index
     * \param[in] j the index
     * \return the value of the parameter at index \p j
     */
    //double paramAtIndex(const size_t j) const { return genome_[j]; }

    /*!
     * \brief Set all parameters to new values
     * \param[in] param the new values
     */
    //void setParam(std::vector<double> param)
    //{
    //  GMX_RELEASE_ASSERT(param.size() == genome_.size() || genome_.empty(),
    //                     "Incorrect size of input parameters");
    //  genome_ = param;
    //}

    //! \return the initial vector of parameters as a const reference
    const ga::Genome &initialGenome() const { return initialGenome_; }

    //! \return the current genome as a const reference
    const ga::Genome &genome() const { return genome_; }

    //! \return a pointer to the current vector of parameters
    ga::Genome *genomePtr() { return &genome_; }

    //! \return the vector of best parameters as a const reference
    const ga::Genome &bestGenome() const { return bestGenome_; }

    //! \return a pointer to the current vector of parameters
    //ga::Genome *bestGenomePtr() { return &bestGenome_; }

    /*!
     * \brief Set a new best parameter vector
     * \param[in] param the new best parameter vector
     */
    void setBestGenome(const ga::Genome &genome) { bestGenome_ = genome; }

    //! \return the name of the Force Field output file as const reference
    //const std::string &outputFile() const { return outputFile_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and Setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};


} //namespace alexandria


#endif //ALEXANDRIA_ACMINDIVIDUAL_H
