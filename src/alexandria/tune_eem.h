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

    /*!
     * \brief Handle out of bounds variables (calcDeviation) in MASTER node before calculating the rest of the deviation
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param verbose   whether we are in verbose mode
     */
    void handleBoundsCD(      std::map<eRMS, FittingTarget>    *targets,
                        const double                            verbose);

    /*!
     * \brief Handle (CM5) charge (calcDeviation)
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param mymol     the molecule
     */
    void handleChargeCM5CD(std::map<eRMS, FittingTarget> *targets,
                           MyMol                          mymol);

    /*!
     * \brief Handle electrostatic potential (calcDeviation)
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param mymol     the molecule
     */
    void handleEspCD(std::map<eRMS, FittingTarget>   *targets,
                     MyMol                            mymol);

    /*!
     * \brief Handle molecular dipole (calcDeviation)
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param mymol     the molecule
     * @param qelec     pointer to Elec properties
     * @param qcalc     pointer to calc properties
     */
    void handleMuCD(std::map<eRMS, FittingTarget>  *targets,
                    MyMol                           mymol,
                    QtypeProps                     *qelec,
                    QtypeProps                     *qcalc);

    /*!
     * \brief Handle molecular quadrupole (calcDeviation)
     * @param targets pointer to a map between the components of chi-squared and the fitting targets
     * @param qelec     pointer to Elec properties
     * @param qcalc     pointer to calc properties
     */
    void handleQuadCD(std::map<eRMS, FittingTarget>    *targets,
                      QtypeProps                       *qelec,
                      QtypeProps                       *qcalc);

    /*!
     * \brief Handle polarizability component (calcDeviation)
     * @param targets   pointer to a map between the components of chi-squared and the fitting targets
     * @param mymol     the molecule
     */
    void handlePolarCD(std::map<eRMS, FittingTarget>   *targets,
                       MyMol                            mymol);

public:
    //! Constructor
    OptACM() {}

    //! \return whether or not we use the full quadrupol
    bool fullQuadrupole() const { return bFullQuadrupole_; }

    //! \return whether or not we remove problematic compounds
    bool removeMol() const { return bRemoveMol_; }

    //! \return whether or not we are in verbose mode
    bool verbose() { return getConfigHandler()->verbose(); }

    //! \brief This function will store the current state of the force field
    void saveState();

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


#endif //ALEXANDRIA_TUNE_EEM_H
