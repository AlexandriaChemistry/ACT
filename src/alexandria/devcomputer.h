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

#include "act/utility/communicationrecord.h"
#include "mymol.h"
#include "molgen.h"
#include "molhandler.h"

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
    //! Whether to print stuff in the logfile
    bool verbose_;

    /*! \brief Create a new DevComputer
     * @param logfile   pointer to log file
     * @param verbose   Whether to print stuff in the logfile
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
    virtual void calcDeviation(      MyMol                         *mymol,
                                     std::map<eRMS, FittingTarget> *targets,
                                     Poldata                       *poldata,
                               const std::vector<double>           &param,
                               const CommunicationRecord           *commrec) = 0;

};


/*!
 * DevComputer that penalizes parameters out of bounds -> eRMS::BOUNDS
 */
class BoundsDevComputer : public DevComputer
{

private:

    //! Information about each force field parameter
    std::vector<OptimizationIndex> *optIndex_;

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

    virtual void calcDeviation(      MyMol                         *mymol,
                                     std::map<eRMS, FittingTarget> *targets,
                                     Poldata                       *poldata,
                               const std::vector<double>           &param,
                               const CommunicationRecord           *commrec);

};

/*!
 * DevComputer that penalizes frequencies or intensities out of bounds ->
 * eRMS::FREQUENCY or eRMS::INTENSITY
 */
class HarmonicsDevComputer : public DevComputer
{
private:
    //! The MolPropObservable, e.g. quadrupole
    MolPropObservable mpo_;
    //! Molecule handler for running frequency calculations
    MolHandler        handler_;
public:

    /*! \brief Create a new FrequencyDevComputer
     * @param logfile pointer to log file
     * @param verbose whether we are in verbose mode
     * @param mpo     indicating which of the two observables are being used
     */
    HarmonicsDevComputer(      FILE              *logfile,
                         const bool               verbose,
                               MolPropObservable  mpo);

    virtual void calcDeviation(      MyMol                         *mymol,
                                     std::map<eRMS, FittingTarget> *targets,
                                     Poldata                       *poldata,
                               const std::vector<double>           &param,
                               const CommunicationRecord           *commrec);
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

    virtual void calcDeviation(      MyMol                         *mymol,
                                     std::map<eRMS, FittingTarget> *targets,
                                     Poldata                       *poldata,
                               const std::vector<double>           &param,
                               const CommunicationRecord           *commrec);

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

    virtual void calcDeviation(      MyMol                         *mymol,
                                     std::map<eRMS, FittingTarget> *targets,
                                     Poldata                       *poldata,
                               const std::vector<double>           &param,
                               const CommunicationRecord           *commrec);

};

/*!
 * DevComputer that computes the deviation of the polarizability components -> eRMS::POLAR
 */
class PolarDevComputer : public DevComputer
{
private:
    //! Conversion from internal to external units
    //! to make the deviation more comprehensible.
    double convert_ = 1;
public:

    /*! \brief Create a new PolarDevComputer
     * @param logfile           pointer to log file
     * @param verbose           whether we are in verbose mode
     */
    PolarDevComputer(    FILE  *logfile,
                    const bool verbose);

    virtual void calcDeviation(      MyMol                         *mymol,
                                     std::map<eRMS, FittingTarget> *targets,
                                     Poldata                       *poldata,
                               const std::vector<double>           &param,
                               const CommunicationRecord           *commrec);

};

/*!
 * DevComputer that computes the deviation of the molecular multipole -> 
 * eRMS::DIPOLE, eRMS::QUAD etc.
 */
class MultiPoleDevComputer : public DevComputer
{
private:
    //! The MolPropObservable, e.g. quadrupole
    MolPropObservable mpo_;
public:

    /*! \brief Create a new MultiPoleDevComputer
     * @param logfile  Pointer to log file
     * @param verbose  Whether we are in verbose mode
     * @param mpo      The MolPropObservable, e.g. quadrupole
     */
    MultiPoleDevComputer(FILE             *logfile,
                         const bool        verbose,
                         MolPropObservable mpo)
        : DevComputer(logfile, verbose), mpo_(mpo)
    {
    }

    virtual void calcDeviation(      MyMol                         *mymol,
                                     std::map<eRMS, FittingTarget> *targets,
                                     Poldata                       *poldata,
                               const std::vector<double>           &param,
                               const CommunicationRecord                         *commrec);

};

/*!
 * DevComputer the computes the deviation of the molecular energy -> eRMS::EPOT
 * and of the forces, depending on command line options.
 */
class ForceEnergyDevComputer : public DevComputer
{

public:

    /*! \brief Create a new ForceEnergyDevComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     */
    ForceEnergyDevComputer(      FILE *logfile,
                      const bool  verbose)
    : DevComputer(logfile, verbose)
    {
    }

    virtual void calcDeviation(      MyMol                         *mymol,
                                     std::map<eRMS, FittingTarget> *targets,
                                     Poldata                       *poldata,
                               const std::vector<double>           &param,
                               const CommunicationRecord           *commrec);

};

} // namespace alexandria

#endif //ALEXANDRIA_DEVCOMPUTER_H
