/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


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
    //! Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
    bool verbose_;

    /*! \brief Create a new DevComputer
     * @param logfile   pointer to log file
     * @param verbose   Flush output immediately rather than letting the OS buffer it. Don't use for production simulations.
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
                               const t_commrec                         *commrec) = 0;

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
                               const t_commrec                         *commrec);

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
                               const t_commrec                         *commrec);

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
     * \TODO move to cpp file
     * \param[in] fp   The file pointer to print to
     * \param[in] mol  The molecule to read from
     * \param[in] info Additional debugging information
     */
    void dumpQX(FILE *fp, MyMol *mol, const std::string &info);
   
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
                               const t_commrec                         *commrec);

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
                               const t_commrec                         *commrec);

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
                               const t_commrec                         *commrec);

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
                               const t_commrec                         *commrec);

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
                               const t_commrec                         *commrec);

};


} // namespace alexandria


#endif //ALEXANDRIA_DEVCOMPUTER_H
