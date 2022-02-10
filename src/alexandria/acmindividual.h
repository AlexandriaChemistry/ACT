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
#include "act/poldata/poldata.h"
#include "staticindividualinfo.h"

#include "act/ga//Genome.h"
#include "act/ga//Individual.h"

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

public:
    /*!
     * \brief Property constructor
     * \param[in] id            the ID of the individual
     * \param[in] sii           pointer to StaticIndividualInfo instance
     * \param[in] outputFile    the base name for Force Field output files
     */
    ACMIndividual(const int             id,
                  StaticIndividualInfo *sii) : ga::Individual(), id_(id), sii_(sii) {}

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Cloning                           *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief Copy the genome to this individual
     * \param[in] genome Complete input genome
     */
    virtual void copyGenome(const ga::Genome &genome);

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
    * BEGIN: Getters and Setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return the ID
    int id() const { return id_; }

    //! \return a pointer to the StaticIndividualInfo instance
    StaticIndividualInfo *sii() { return sii_; }

    //! \return a pointer to the StaticIndividualInfo instance (for const objects)
    StaticIndividualInfo *siiConst() const { return sii_; }

    //! \return the initial vector of parameters as a const reference
    const ga::Genome &initialGenome() const { return initialGenome_; }

    //! \return the current genome as a const reference
    const ga::Genome &genome() const { return genome_; }

    //! \return a pointer to the current vector of parameters
    ga::Genome *genomePtr() { return &genome_; }

    /*!
     * \brief Set a new parameter vector
     * \param[in] param the new parameter vector
     */
    void setGenome(const ga::Genome &genome) { genome_ = genome; }

    //! \return the vector of best parameters as a const reference
    const ga::Genome &bestGenome() const { return bestGenome_; }

    //! \return a pointer to the current vector of parameters
    ga::Genome *bestGenomePtr() { return &bestGenome_; }

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
