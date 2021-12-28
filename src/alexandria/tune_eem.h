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
#include "sharedindividualinfo.h"
#include "acmfitnesscomputer.h"
#include "acminitializer.h"
#include "mcmcmutator.h"
#include "ga/GeneticAlgorithm.h"


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
    //! Whether or not to use off-diagonal elements of the quadrupole for fitting
    bool bFullQuadrupole_ = false;
    //! Whether or not to remove molecules that fail to converge in the shell minimization
    bool bRemoveMol_ = true;
    //! Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
    bool verbose_ = false;
    //! Pointer to log file
    gmx::unique_cptr<FILE, my_fclose> fplog_ = nullptr;
    //! File name for the output force field file
    std::string outputFile_;
    //! GROMACS communication data structure
    t_commrec *cr_ = nullptr;
    //! GROMACS output environment
    gmx_output_env_t *oenv_;
    //! MolGen instance
    MolGen mg_;
    //! BayesConfigHandler instance
    BayesConfigHandler bch_;
    //! GAConfigHandler instance
    GAConfigHandler gach_;
    //! SharedIndividualInfo instance
    SharedIndividualInfo sii_;

    // This is for MASTER node
    //! GeneticAlgorithm instance
    ga::GeneticAlgorithm *ga_;
    //! A pointer to a future best individual
    ACMIndividual *bestInd_;

    // This is for the HELPER node (and MASTER when sensitivity)
    //! Pointer to ACMFitnessComputer since it will be initialized later
    ACMFitnessComputer *fitComp_;
    //! Pointer to ACMInitializer since it will be initialized later
    ACMInitializer *initializer_;
    //! Pointer to MCMCMutator instance
    MCMCMutator *mutator_;
    //! Pointer to the Base class individual
    ga::Individual *baseInd_;
    //! Pointer to ACMIndividual which will be casted from \p baseInd_
    ACMIndividual *ind_;

public:

    //! Constructor
    OptACM()
    : cr_(init_commrec()), mg_(cr_), sii_(cr_) {}

    //! Destructor
    ~OptACM()
    {
        if (cr_) done_commrec(cr_);
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

    /*! \brief
     * Do the actual optimization.
     * \param[in] xvgconv     Output file monitoring parameters
     * \param[in] xvgepot     Output file monitoring penalty function
     * \param[in] optimize    If true an optimization will be done
     * \param[in] sensitivity If true, a sensitivity analysis will be done
     * \return true if better parameters were found.
     */
    bool runMaster(const char *xvgconv,
                   const char *xvgepot,
                   bool optimize,
                   bool sensitivity);

    /*! \brief
     * For the helper nodes.
     */
    void runHelper();

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Initializing stuff                *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \brief Initialize the ACMFitnessComputer
    void initFitComp()
    {
        fitComp_ = new ACMFitnessComputer(cr_, logFile(), &sii_, &mg_, bRemoveMol_, verbose_, bFullQuadrupole_);
    }

    //! \brief Initialize the ACMInitializer
    void initInitializer()
    {
        initializer_ = new ACMInitializer(mg_.mindata(), &sii_, gach_.randomInit(), outputFile_);
    }

    //! \brief initialize the Individual \p baseInd_ and ACMIndividual \p ind_
    void initIndividual()
    {
        initializer_->initialize(&baseInd_);
        ind_ = static_cast<ACMIndividual*>(baseInd_);
    }

    //! \brief Initialize the MCMCMutator
    void initMutator()
    {
        mutator_ = new MCMCMutator(logFile(), verbose_, &bch_, fitComp_, &sii_, sii_.nParam());
    }

    //! \brief Initialize the Genetic Algorithm
    void initGA()
    {
        ga_ = new ga::GeneticAlgorithm(verbose_, bRemoveMol_, bFullQuadrupole_, cr_, &mg_,
                                       logFile(), oenv_, &bch_, &sii_, &gach_, outputFile_);
    }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Initializing stuff                  *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return whether or not we use the full quadrupol
    bool fullQuadrupole() const { return bFullQuadrupole_; }

    //! \return whether or not we remove problematic compounds
    bool removeMol() const { return bRemoveMol_; }

    //! \return whether or not we are in verbose mode
    bool verbose() const { return verbose_; }

    //! \brief Get the constant communications record \p cr_ constant pointer
    const t_commrec *commrec() const { return cr_; }

    //! \brief Get the communications record \p cr_ pointer
    t_commrec *commrec() { return cr_; }

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

    //! \brief Get the SharedIndividualInfo \p sii_ pointer
    SharedIndividualInfo *sii() { return &sii_; }

    //! \brief Get the ACMFitnessComputer \p fitComp_ pointer
    ACMFitnessComputer *fitComp() { return fitComp_; }

    //! \brief Get the ACMIndividual \p ind_ pointer
    ACMIndividual *ind() { return ind_; }

    //! \brief Get the GeneticAlgorithm \p ga_ pointer
    ga::GeneticAlgorithm *ga() { return ga_; }

    //! \brief Get the \p bestInd_ pointer
    ACMIndividual *bestInd() { return bestInd_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};

}


#endif //ALEXANDRIA_TUNE_EEM_H
