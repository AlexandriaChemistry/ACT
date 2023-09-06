/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
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

#include "act/utility/communicationrecord.h"
#include "confighandler.h"
#include "actmol.h"
#include "staticindividualinfo.h"
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

    /*! \brief Computes a component of the chi-squared deviation.
     * Put the value into the appropriate FittingTarget.
     * @param forceComputer The utility to compute forces and energies
     * @param actmol     pointer to the molecule to compute the deviation for
     * @param coords    the coordinates
     * @param targets   map between different components of chi-squared and their fitting targets
     * @param forcefield   pointer to ForceField structure, that should contain the newest parameters
     */
    virtual void calcDeviation(const ForceComputer           *forceComputer,
                               ACTMol                         *actmol,
                               std::vector<gmx::RVec>        *coords,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField                 *forcefield) = 0;

};


/*!
 * DevComputer that penalizes parameters out of bounds -> eRMS::BOUNDS
 */
class BoundsDevComputer : public DevComputer
{

private:

    //! Information about each force field parameter
    std::vector<OptimizationIndex> *optIndex_;
    //! Difference between shell and core zeta to contribute to the unphysical
    double zetaDiff_ = 2;
public:

    /*! \brief Create a new BoundsDevComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     * @param optIndex  pointer to vector containing information about each force field parameter
     */
    BoundsDevComputer(      FILE                           *logfile,
                      const bool                            verbose,
                            std::vector<OptimizationIndex> *optIndex,
                            double                          zetaDiff)
    : DevComputer(logfile, verbose)
    {
        optIndex_ = optIndex;
        zetaDiff_ = zetaDiff;
    }

    virtual void calcDeviation(const ForceComputer           *forceComputer,
                               ACTMol                         *actmol,
                               std::vector<gmx::RVec>        *coords,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField                 *forcefield);

};

/*!
 * DevComputer that penalizes frequencies or intensities out of bounds ->
 * eRMS::FREQUENCY or eRMS::INTENSITY
 */
class HarmonicsDevComputer : public DevComputer
{
private:
    //! The MolPropObservable, e.g. quadrupole
    MolPropObservable       mpo_;
    //! Molecule handler for running frequency calculations
    MolHandler              handler_;
    //! Config handler for options
    SimulationConfigHandler simConfig_;
public:

    /*! \brief Create a new FrequencyDevComputer
     * @param logfile pointer to log file
     * @param verbose whether we are in verbose mode
     * @param mpo     indicating which of the two observables are being used
     */
    HarmonicsDevComputer(      FILE              *logfile,
                         const bool               verbose,
                               MolPropObservable  mpo);

    virtual void calcDeviation(const ForceComputer                 *forceComputer,
                                     ACTMol                         *actmol,
                                     std::vector<gmx::RVec>        *coords,
                                     std::map<eRMS, FittingTarget> *targets,
                               const ForceField                       *forcefield);
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

    virtual void calcDeviation(const ForceComputer                 *forceComputer,
                                     ACTMol                         *actmol,
                                     std::vector<gmx::RVec>        *coords,
                                     std::map<eRMS, FittingTarget> *targets,
                               const ForceField                       *forcefield);

};

/*!
 * DevComputer that computes the deviation of the electrostatic potential -> eRMS::ESP
 */
class EspDevComputer : public DevComputer
{

private:

    //! Whether we fit zeta parameters
    bool fitZeta_;

    /*! \brief Dump charges to a file
     * Debugging routine
     * \TODO move to cpp file
     * \param[in] fp   The file pointer to print to
     * \param[in] mol  The molecule to read from
     * \param[in] coords The atomic coordinates
     * \param[in] info Additional debugging information
     */
    void dumpQX(FILE                         *fp,
                const ACTMol                  *mol,
                const std::vector<gmx::RVec> &coords,
                const std::string            &info);
   
public:

    /*! \brief Create a new DevComputer
     * @param logfile   pointer to log file
     * @param verbose   whether we are in verbose mode
     * @param fit       whether we fit zeta parameters
     */
    EspDevComputer(      FILE *logfile,
                   const bool  verbose,
                   const bool  fitZeta)
    : DevComputer(logfile, verbose)
    {
        fitZeta_ = fitZeta;
    }

    virtual void calcDeviation(const ForceComputer                 *forceComputer,
                                     ACTMol                         *actmol,
                                     std::vector<gmx::RVec>        *coords,
                                     std::map<eRMS, FittingTarget> *targets,
                               const ForceField                       *forcefield);

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

    virtual void calcDeviation(const ForceComputer                 *forceComputer,
                                     ACTMol                         *actmol,
                                     std::vector<gmx::RVec>        *coords,
                                     std::map<eRMS, FittingTarget> *targets,
                               const ForceField                       *forcefield);

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

    virtual void calcDeviation(const ForceComputer                 *forceComputer,
                                     ACTMol                         *actmol,
                                     std::vector<gmx::RVec>        *coords,
                                     std::map<eRMS, FittingTarget> *targets,
                               const ForceField                    *forcefield);

};

/*!
 * DevComputer the computes the deviation of the molecular energy -> eRMS::EPOT
 * and of the forces, depending on command line options.
 */
class ForceEnergyDevComputer : public DevComputer
{
private:
    //! Whether or not to use Boltzmann weighting of energies and forces
    std::map<eRMS, double> boltzmannTemperature_;

    /*! \brief Compute Boltzmann weighting factor for this eRMS
     * \param[in] ermsi Which energy or force
     * \return The corresponding beta or zero if no temperature was given
     */
    double computeBeta(eRMS ermsi);
public:

    /*! \brief Create a new ForceEnergyDevComputer
     * \param[in] logfile              pointer to log file
     * \param[in] verbose              whether we are in verbose mode
     * \param[in] boltzmannTemperature Map from the different eRMS terms to a number
     *                                 indicating whether (>0) or not (0) Boltzmann weighting
     *                                 of energies and forces should be done.
     */
    ForceEnergyDevComputer(      FILE                   *logfile,
                           const bool                    verbose,
                                 std::map<eRMS, double>  boltzmannTemperature)
    : DevComputer(logfile, verbose)
    {
        boltzmannTemperature_ = boltzmannTemperature;
    }

    virtual void calcDeviation(const ForceComputer                 *forceComputer,
                                     ACTMol                         *actmol,
                                     std::vector<gmx::RVec>        *coords,
                                     std::map<eRMS, FittingTarget> *targets,
                               const ForceField                       *forcefield);

};

} // namespace alexandria

#endif //ALEXANDRIA_DEVCOMPUTER_H
