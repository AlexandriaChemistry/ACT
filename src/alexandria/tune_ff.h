/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#ifndef ALEXANDRIA_TUNE_EEM_H
#define ALEXANDRIA_TUNE_EEM_H

#include <cstdio>

#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "molgen.h"
#include "bayes.h"
#include "staticindividualinfo.h"
#include "acmfitnesscomputer.h"
#include "acminitializer.h"
#include "act/utility/communicationrecord.h"
#include "act/ga//Mutator.h"
#include "act/ga//GeneticAlgorithm.h"


namespace alexandria
{

/*! \brief Wrapper for closing a file
 * Will print a warning if something is wrong when closing.
 * \param[in] fp The file pointer
 */
void my_fclose(FILE *fp);

/*! \brief FIXME: The class that does all the optimization work.
 * This class inherits the MolGen class that holds molecule data and
 * properties, and the Bayes class, that does the Monte Carlo steps.
 *
 * The class can run as a Master process or as a Helper process.
 */
class OptACM : public ConfigHandler
{

private:
    //! Whether or not to remove molecules that fail to converge in the shell minimization
    bool bRemoveMol_ = true;
    //! Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
    bool flush_ = false;
    //! Print extra stuff during the optimization
    bool verbose_ = false;
    //! Pointer to log file
    gmx::unique_cptr<FILE, my_fclose> fplog_ = nullptr;
    //! ACT Communication data structure
    CommunicationRecord   commRec_;
    //! GROMACS output environment
    gmx_output_env_t     *oenv_ = nullptr;
    //! MolGen instance
    MolGen                mg_;
    //! BayesConfigHandler instance
    BayesConfigHandler    bch_;
    //! GAConfigHandler instance
    GAConfigHandler       gach_;
    //! StaticIndividualInfo instance
    StaticIndividualInfo *sii_;
    //! Base name for FF output
    std::string           baseOutputFileName_;

    // This is for MASTER node
    //! GeneticAlgorithm instance
    ga::GeneticAlgorithm *ga_          = nullptr;

    //! Pointer to ACMFitnessComputer since it will be initialized later
    ACMFitnessComputer   *fitComp_     = nullptr;

    /*!
     * \brief Print to log file (if it exists), the estimated number of times
     * we will call calcDeviation per dataset
     */
    void printNumCalcDevEstimate();

    /*!
     * \brief Print a table about the parameters in a best genome and a population
     * to the logfile, if it exists
     * \param[in] genome the best Genome object
     * \param[in] pop    the population to define summary statistics
     */
    void printGenomeTable(const ga::Genome   &genome,
                          const ga::GenePool &pop);

public:

    //! Constructor
    OptACM() : mg_(&commRec_)
    {
        sii_ = new StaticIndividualInfo(&commRec_);
    }

    virtual void add_pargs(std::vector<t_pargs> *pargs);

    virtual void check_pargs();

    /*! \brief Routine to be called after processing options
     * \param[in] outputFile The force field target file
     */
    void optionsFinished(const std::string &outputFile);

    /*! \brief Routine that opens a log file
     * \param[in] logfileName The log file name to open
     */
    void openLogFile(const char *logfileName);

    //! \return a file pointer to the open logfile
    FILE *logFile();

    /*! \brief Initialize charge generation
     * \param[in] ims The data set to do the work for
     */
    void initChargeGeneration(iMolSelect ims);

    /*! \brief Do the actual optimization.
     * \param[in] optimize    If true an optimization will be done
     * \param[in] sensitivity If true, a sensitivity analysis will be done
     * \return true if better parameters were found.
     */
    bool runMaster(bool optimize,
                   bool sensitivity);

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Initializing stuff                *
    * * * * * * * * * * * * * * * * * * * * * */

    /*! \brief Initialize the main components of the
     * Genetic Algorithm, just on the master.
     */
    void initMaster();

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Initializing stuff                  *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return whether or not we remove problematic compounds
    bool removeMol() const { return bRemoveMol_; }

    //! \return whether or not we are in verbose mode
    bool verbose() const { return verbose_; }

    const CommunicationRecord *commRec() const { return &commRec_; }
    
    /*! \brief Set the output environment pointer \p oenv_
     * \param[in] oenv the reference pointer
     */
    void setOenv(gmx_output_env_t *oenv) { oenv_ = oenv; }

    //! \brief Get the output environment \p oenv_ pointer
    gmx_output_env_t *oenv() { return oenv_; }

    //! \brief Get the MolGen \p mg_ pointer
    MolGen *mg() { return &mg_; }

    //! \brief Get the BayesConfigHandler \p bch_ pointer
    BayesConfigHandler *bch() { return &bch_; }

    //! \brief Get the GAConfigHandler \p gach_ pointer
    GAConfigHandler *gach() { return &gach_; }

    //! \brief Get the StaticIndividualInfo \p sii_ pointer
    StaticIndividualInfo *sii() { return sii_; }

    //! \brief Get the GeneticAlgorithm \p ga_ pointer
    ga::GeneticAlgorithm *ga() { return ga_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};

}


#endif //ALEXANDRIA_TUNE_EEM_H
