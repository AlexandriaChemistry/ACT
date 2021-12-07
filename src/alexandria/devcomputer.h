#ifndef ALEXANDRIA_DEVCOMPUTER_H
#define ALEXANDRIA_DEVCOMPUTER_H


#include <cstdio>
#include <string>
#include <map>
#include <vector>

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/topology/topology.h"

#include "mymol.h"
#include "molgen.h"


namespace alexandria
{


/*!
 * Abstract class for chi-squared deviation computation
 * Each DevComputer handles a particular component of chi-squared
 */
class DevComputer
{

protected:

    //! Pointer to log file
    FILE *logfile_;
    //! Whether we are in verbose mode
    bool verbose_;

    /*! \brief Create a new DevComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     */
    DevComputer(      FILE *logfile,
                const bool  verbose)
    {
        logfile_ = logfile;
        verbose_ = verbose;
    }

public:

    /*! \brief Computes a component of the chi-squared deviation and place it into the appropriate FittingTarget
     * @param mymol     pointer to the molecule to compute the deviation for
     * @param targets   map between different components of chi-squared and their fitting targets
     * @param poldata   pointer to Poldata structure
     * @param param     the current force field parameter vector
     * @param commrec   pointer to communications record
     */
    virtual void calcDeviation(      MyMol                             *mymol,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     Poldata                           *poldata,
                               const std::vector<double>               &param,
                                     t_commrec                         *commrec) = 0;

};


/*!
 * DevComputer that penalizes parameters out of bounds -> eRMS::BOUNDS
 */
class BoundsDevComputer : public DevComputer
{

private:

    //! Information about each force field parameter
    std::vector<OptimizationIndex> *optIndex_;

    /*! \brief Compute penalty for variables that are out of bounds
     * \param[in] x       The actual value
     * \param[in] min     The minimum allowed value
     * \param[in] max     The maximum allowed value
     * \param[in] label   String to print if verbose
     * \return 0 when in bounds, square deviation from bounds otherwise.
     */
    double l2_regularizer(      double          x,
                                double          min,
                                double          max,
                          const std::string    &label);

public:

    /*! \brief Create a new BoundsDevComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     * @param optIndex  pointer to vector containing information about each force field parameter
     */
    BoundsDevComputer(      FILE                           *logfile,
                      const bool                            verbose,
                            std::vector<OptimizationIndex> *optIndex)
    : DevComputer(logfile, verbose)
    {
        optIndex_ = optIndex;
    }

    virtual void calcDeviation(      MyMol                             *mymol,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     Poldata                           *poldata,
                               const std::vector<double>               &param,
                                     t_commrec                         *commrec);

};

/*!
 * DevComputer that computes the deviation of the (CM5) charge
 * -> eRMS::CHARGE & eRMS::CM5
 */
class ChargeCM5DevComputer : public DevComputer
{

public:

    /*! \brief Create a new DevComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     */
    ChargeCM5DevComputer(      FILE *logfile,
                         const bool  verbose)
    : DevComputer(logfile, verbose)
    {}

    virtual void calcDeviation(      MyMol                             *mymol,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     Poldata                           *poldata,
                               const std::vector<double>               &param,
                                     t_commrec                         *commrec);

};

/*!
 * DevComputer that computes the deviation of the electrostatic potential -> eRMS::ESP
 */
class EspDevComputer : public DevComputer
{

private:

    //! Whether we fit zeta parameters
    bool fit_;

    /*! \brief Dump charges to a file
     * Debugging routine
     * \param[in] fp   The file pointer to print to
     * \param[in] mol  The molecule to read from
     * \param[in] info Additional debugging information
     */
    static void dumpQX(FILE *fp, MyMol *mol, const std::string &info)
    {
        if (false && fp)
        {
            std::string label = mol->getMolname() + "-" + info;
            fprintf(fp, "%s q:", label.c_str());
            t_mdatoms *md = mol->getMdatoms();
            auto myatoms = mol->atomsConst();
            for (int i = 0; i < myatoms.nr; i++)
            {
                fprintf(fp, " %g (%g)", myatoms.atom[i].q,
                        md->chargeA[i]);
            }
            fprintf(fp, "\n");
            fprintf(fp, "%s alpha", label.c_str());
            int ft = F_POLARIZATION;
            for(int i = 0; i < mol->ltop_->idef.il[ft].nr; i += interaction_function[ft].nratoms+1)
            {
                auto tp = mol->ltop_->idef.il[ft].iatoms[i];
                fprintf(fp, " %g", mol->ltop_->idef.iparams[tp].polarize.alpha);
            }
            fprintf(fp, "\n");
            pr_rvecs(fp, 0, label.c_str(), mol->x().rvec_array(), myatoms.nr);
        }
    }

public:

    /*! \brief Create a new DevComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     * @param fit       whether we fit zeta parameters
     */
    EspDevComputer(      FILE *logfile,
                   const bool  verbose,
                   const bool  fit)
    : DevComputer(logfile, verbose)
    {
        fit_ = fit;
    }

    virtual void calcDeviation(      MyMol                             *mymol,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     Poldata                           *poldata,
                               const std::vector<double>               &param,
                                     t_commrec                         *commrec);

};

/*!
 * DevComputer that computes the deviation of the polarizability components -> eRMS::POLAR
 */
class PolarDevComputer : public DevComputer
{

private:

    //! Whether or not to use off-diagonal elements of the quadrupole for fitting
    bool bFullQuadrupole_;

public:

    /*! \brief Create a new BoundsDevComputer
     * @param logfile           pointer to log file
     * @param verbose           whether we are in verbose mode
     * @param bFullQuadrupole   whether or not to use off-diagonal elements of the quadrupole for fitting
     */
    PolarDevComputer(    FILE  *logfile,
                     const bool verbose,
                     const bool bFullQuadrupole)
    : DevComputer(logfile, verbose)
    {
        bFullQuadrupole_ = bFullQuadrupole;
    }

    virtual void calcDeviation(      MyMol                             *mymol,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     Poldata                           *poldata,
                               const std::vector<double>               &param,
                                     t_commrec                         *commrec);

};

/*!
 * DevComputer that computes the deviation of the molecular quadrupole -> eRMS::QUAD
 */
class QuadDevComputer : public DevComputer
{

private:

    //! Whether or not to use off-diagonal elements of the quadrupole for fitting
    bool bFullQuadrupole_;

public:

    /*! \brief Create a new QuadDevComputer
     * @param logfile           pointer to log file
     * @param verbose           whether we are in verbose mode
     * @param bFullQuadrupole   whether or not to use off-diagonal elements of the quadrupole for fitting
     */
    QuadDevComputer(      FILE *logfile,
                    const bool  verbose,
                    const bool  bFullQuadrupole)
    : DevComputer(logfile, verbose)
    {
        bFullQuadrupole_ = bFullQuadrupole;
    }

    virtual void calcDeviation(      MyMol                             *mymol,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     Poldata                           *poldata,
                               const std::vector<double>               &param,
                                     t_commrec                         *commrec);

};

/*!
 * DevComputer the computes the deviation of the molecular dipole -> eRMS::MU
 */
class MuDevComputer : public DevComputer
{

private:
    //! Are we using QM only?
    bool bQM_;

public:

    /*! \brief Create a new MuComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     * @param bQM       are we using QM only?
     */
    MuDevComputer(      FILE *logfile,
                  const bool  verbose,
                  const bool  bQM)
    : DevComputer(logfile, verbose)
    {
        bQM_ = bQM;
    }

    virtual void calcDeviation(      MyMol                             *mymol,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     Poldata                           *poldata,
                               const std::vector<double>               &param,
                                     t_commrec                         *commrec);

};

/*!
 * DevComputer the computes the deviation of the molecular energy -> eRMS::EPOT
 */
class EnergyDevComputer : public DevComputer
{

public:

    /*! \brief Create a new EnergyDevComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     */
    EnergyDevComputer(      FILE *logfile,
                      const bool  verbose)
    : DevComputer(logfile, verbose)
    {
    }

    virtual void calcDeviation(      MyMol                             *mymol,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     Poldata                           *poldata,
                               const std::vector<double>               &param,
                                     t_commrec                         *commrec);

};


} // namespace alexandria


#endif //ALEXANDRIA_DEVCOMPUTER_H
