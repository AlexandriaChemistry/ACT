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
#include "devcomputer.h"


namespace alexandria
{

/*! \brief Wrapper for closing a file
 * Will print a warning if something is wrong when closing.
 * \param[in] fp The file pointer
 */
void my_fclose(FILE *fp);

/*! \brief The class that does all the optimization work.
 * This class inherits the MolGen class that holds molecule data and
 * properties, and the Bayes class, that does the Monte Carlo steps.
 *
 * The class can run as a Master process or as a Helper process.
 */
class OptACM : public MolGen, public Bayes
{
    //! Abbreviation to make clear this is not a random vector
    using param_type = std::vector<double>;

private:
    //! Whether or not to use off-diagonal elements of the quadrupole for fitting
    bool bFullQuadrupole_ = false;
    //! Whether or not to remove molecules that fail to converge in the shell minimization
    bool bRemoveMol_ = true;
    //! How many times has calcDeviation been called
    int numberCalcDevCalled_ = 0;
    //! Pointer to log file
    gmx::unique_cptr<FILE, my_fclose> fplog_ = nullptr;
    //! File name for the output force field file
    std::string outputFile_;
    //! BoundsDevComputer, if required by the user
    BoundsDevComputer *bdc_ = nullptr;
    //! Vector of non-bound DevComputers for the different components of chi-squared
    std::vector<DevComputer*> devComputers_;

public:
    //! Constructor
    OptACM() {}

    //! \return whether or not we use the full quadrupol
    bool fullQuadrupole() const { return bFullQuadrupole_; }

    //! \return whether or not we remove problematic compounds
    bool removeMol() const { return bRemoveMol_; }

    //! \return whether or not we are in verbose mode
    bool verbose() { return configHandlerPtr()->verbose(); }

    /*! \brief Add our command line parameters to the array.
     * This also calls the addOptions routine of the child class Bayes.
     * \param[inout] pargs The vector of parameters
     */
    void add_pargs(std::vector<t_pargs> *pargs);

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
     * \param[in] oenv        Output environment for managing xvg files etc.
     * \param[in] xvgconv     Output file monitoring parameters
     * \param[in] xvgepot     Output file monitoring penalty function
     * \param[in] optimize    If true an optimization will be done
     * \param[in] sensitivity If true, a sensitivity analysis will be done
     * \param[in] bEvaluate_testset Whether or not to evaluate the test set
     * \return true if better parameters were found.
     */
    bool runMaster(const gmx_output_env_t *oenv,
                   const char *xvgconv,
                   const char *xvgepot,
                   bool optimize,
                   bool sensitivity,
                   bool bEvaluate_testset);

    /*! \brief
     * For the helper nodes.
     */
    void runHelper();
};

}


#endif //ALEXANDRIA_TUNE_EEM_H
