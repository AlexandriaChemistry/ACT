#ifndef ACT_TUNE_EEM_H
#define ACT_TUNE_EEM_H

#include "molgen.h"
#include "optparam.h"

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/unique_cptr.h"

#include "alex_modules.h"
#include "gentop_core.h"
#include "gmx_simple_comm.h"
#include "memory_check.h"
#include "molgen.h"
#include "molprop_util.h"
#include "mymol_low.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_tables.h"
#include "poldata_xml.h"
#include "tuning_utility.h"
#include "units.h"

namespace alexandria {

/*! \brief The class that does all the optimization work.
 * This class inherits the MolGen class that holds molecule data and
 * properties, and the Bayes class, that does the Monte Carlo steps.
 *
 * The class can run as a Master process or as a Helper process.
 */
class OptACM : public MolGen, Bayes {
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

public:
    //! Constructor
    OptACM() {}

    //! \return whether or not we use the full quadrupol
    bool fullQuadrupole() const { return bFullQuadrupole_; }

    //! \return whether or not we remove problematic compounds
    bool removeMol() const { return bRemoveMol_; }

    //! \return whether or not we are in verbose mode
    bool verbose() { return Bayes::verbose(); }

    //! \brief This function will store the current state of the force field
    void saveState();

    /*! \brief Add our command line parameters to the array.
     * This also calls the addOptions routine of the child class Bayes.
     * \param[inout] pargs The vector of parameters
     */
    void add_pargs(std::vector <t_pargs> *pargs) {
        t_pargs pa[] =
                {
                        {"-fullQuadrupole", FALSE, etBOOL, {&bFullQuadrupole_},
                                "Consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization"},
                        {"-removemol",      FALSE, etBOOL, {&bRemoveMol_},
                                "Remove a molecule from training set if shell minimzation does not converge."},
                };
        for (int i = 0; i < asize(pa); i++) {
            pargs->push_back(pa[i]);
        }
        addOptions(pargs, eTune::EEM);
        Bayes::add_pargs(pargs);
    }

    /*! \brief Routine to be called after processing options
     * \param[in] outputFile The force field target file
     */
    void optionsFinished(const std::string &outputFile) {
        MolGen::optionsFinished();
        outputFile_ = outputFile;
    }

    /*! \brief Routine that opens a log file
     * \param[in] logfileName The log file name to open
     */
    void openLogFile(const char *logfileName) {
        fplog_.reset(gmx_ffopen(logfileName, "w"));
    }

    //! \return a filepointer to the open logfile
    FILE *logFile() {
        if (fplog_) {
            return fplog_.get();
        } else {
            return nullptr;
        }
    }

    /*! \brief Initialize charge generation
     * \param[in] ims The data set to do the work for
     */
    void initChargeGeneration(iMolSelect ims);

    /*! \brief
     *
     * Fill parameter vector based on Poldata.
     * \param[in] bRandom Generate random initial values for parameters if true
     */
    void InitOpt(bool bRandom);

    /*! \brief
     * Copy the optimization parameters to the poldata structure
     * \param[in] changed List over the parameters that have changed.
     */
    virtual void toPoldata(const std::vector<bool> &changed);

    /*! \brief
     * Copy the optimization parameters to the poldata structure
     * @param param     vector of parameters
     * @param psigma    standard deviation of each parameter
     */
    virtual void toPoldata(const std::vector<double> &param,
                           const std::vector<double> &psigma);

    /*! \brief
     * Computes deviation from target
     * \param[in] verbose Whether or not to print a lot
     * \param[in] calcDev The type of calculation to do
     * \param[in] ims     The data set to do computations on
     * \return the square deviation
     */
    virtual double calcDeviation(bool verbose,
                                 CalcDev calcDev,
                                 iMolSelect ims);

    /*! \brief Compute penalty for variables that are out of bounds
     * \param[in] x       The actual value
     * \param[in] min     The minimum allowed value
     * \param[in] max     The maximum allowed value
     * \param[in] label   String to print if verbose
     * \param[in] verbose Whether or not to print deviations
     * \return 0 when in bounds, square deviation from bounds otherwise.
     */
    double l2_regularizer(double x,
                          double min,
                          double max,
                          const std::string &label,
                          bool verbose);

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


#endif //ACT_TUNE_EEM_H
