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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
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

namespace gmx
{
class TextWriter;
}

namespace alexandria
{

/*!
 * Abstract class for chi-squared deviation computation
 * Each DevComputer handles a particular component of chi-squared
 */
class DevComputer
{

protected:
    //! My name
    std::string name_;
    /*! \brief Create a new DevComputer
     * @param name My name
     */
    DevComputer(const std::string  name) : name_(name)
    {}

public:

    /*! \brief Computes a component of the chi-squared deviation.
     * Put the value into the appropriate FittingTarget.
     * @param msghandler    For status and output
     * @param forceComputer The utility to compute forces and energies
     * @param actmol        Pointer to the molecule to compute the deviation for
     * @param coords        The coordinates
     * @param targets       Map between different components of chi-squared and their fitting targets
     * @param forcefield    Pointer to ForceField structure, that should contain the newest parameters
     */
    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               ACTMol                        *actmol,
                               std::vector<gmx::RVec>        *coords,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield) = 0;
    //! \brief Return my name
    const std::string name() const { return name_; }
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
     * @param tw        Pointer to TextWriter
     * @param verbose   whether we are in verbose mode
     * @param optIndex  pointer to vector containing information about each force field parameter
     */
    BoundsDevComputer(std::vector<OptimizationIndex> *optIndex,
                      double                          zetaDiff)
        : DevComputer("Bounds")
    {
        optIndex_ = optIndex;
        zetaDiff_ = zetaDiff;
    }

    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               ACTMol                        *actmol,
                               std::vector<gmx::RVec>        *coords,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield);

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
     * @param mpo     indicating which of the two observables are being used
     */
    HarmonicsDevComputer(MolPropObservable  mpo);

    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               ACTMol                        *actmol,
                               std::vector<gmx::RVec>        *coords,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield);
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
     * \param[in] tw   The textwriter
     * \param[in] mol  The molecule to read from
     * \param[in] info Additional debugging information
     */
    void dumpQX(gmx::TextWriter              *tw,
                const ACTMol                 *mol,
                const std::string            &info);
   
public:

    /*! \brief Create a new DevComputer
     * @param fit       whether we fit zeta parameters
     */
    EspDevComputer(const bool       fitZeta)
        : DevComputer("ESP")
    {
        fitZeta_ = fitZeta;
    }

    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               ACTMol                        *actmol,
                               std::vector<gmx::RVec>        *coords,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield);

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
     */
    PolarDevComputer();

    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               ACTMol                        *actmol,
                               std::vector<gmx::RVec>        *coords,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield);

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
     * @param tw       Pointer to TextWriter
     * @param verbose  Whether we are in verbose mode
     * @param mpo      The MolPropObservable, e.g. quadrupole
     */
    MultiPoleDevComputer(MolPropObservable mpo)
        : DevComputer(mpo_name(mpo)), mpo_(mpo)
    {
    }

    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               ACTMol                        *actmol,
                               std::vector<gmx::RVec>        *coords,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield);

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

    //! Whether to split off induction correction from induction
    bool        separateInductionCorrection_ = true;

    /*! \brief Compute Boltzmann weighting factor for this eRMS
     * \param[in] ermsi Which energy or force
     * \return The corresponding beta or zero if no temperature was given
     */
    double computeBeta(eRMS ermsi);
public:

    /*! \brief Set internal variable
     * \param[in] sepIndCorr The new value
     */
    void setSeparateInductionCorrection(bool sepIndCorr) { separateInductionCorrection_ = sepIndCorr; }

    //! \return Internal variable
    bool separateInductionCorrection() const { return separateInductionCorrection_; }

    /*! \brief Create a new ForceEnergyDevComputer
     * \param[in] boltzmannTemperature Map from the different eRMS terms to a number
     *                                 indicating whether (>0) or not (0) Boltzmann weighting
     *                                 of energies and forces should be done.
     */
    ForceEnergyDevComputer(std::map<eRMS, double>  boltzmannTemperature);

    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               ACTMol                        *actmol,
                               std::vector<gmx::RVec>        *coords,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield);

};

} // namespace alexandria

#endif //ALEXANDRIA_DEVCOMPUTER_H
