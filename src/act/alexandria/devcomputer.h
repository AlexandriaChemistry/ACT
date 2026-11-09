/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#include "act/alexandria/confighandler.h"
#include "act/alexandria/loss_function.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/utility/communicationrecord.h"

namespace gmx
{
class TextWriter;
}

namespace alexandria
{

class ActMol;
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
    //! \brief Destructor
    virtual ~DevComputer() = default;

    /*! \brief Computes a component of the chi-squared deviation.
     * Put the value into the appropriate FittingTarget.
     * @param msghandler    For status and output
     * @param forceComputer The utility to compute forces and energies
     * @param actmols       The molecules to compute the deviation for
     * @param targets       Map between different components of chi-squared and their fitting targets
     * @param forcefield    Pointer to ForceField structure, that should contain the newest parameters
     * @param lossFunction  What functional form to use for the loss function
     */
    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               const std::vector<ACTMol *>   &actmols,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield,
                               LossFunction                   lossFunction) = 0;
    //! \return my name
    const std::string name() const { return name_; }
};


/*!
 * DevComputer that penalizes parameters out of bounds -> eRMS::BOUNDS
 */
class BoundsDevComputer : public DevComputer
{

private:

    //! Information about each force field parameter
    const std::vector<OptimizationIndex> *optIndex_;
public:

    /*! \brief Create a new BoundsDevComputer
     * @param optIndex  pointer to vector containing information about each force field parameter
     */
    BoundsDevComputer(const std::vector<OptimizationIndex> *optIndex)
        : DevComputer("Bounds")
    {
        optIndex_ = optIndex;
    }

    /*! \brief Add the penalty for total atomic charges (including shells and vsites)
     * outside the bounds of their ACM charge parameter to the eRMS::BOUNDS target
     * \param[in]    msghandler    Not used
     * \param[in]    forceComputer Not used
     * \param[in]    actmols       The molecules to compute the deviation for
     * \param[inout] targets       Map of fitting targets, the eRMS::BOUNDS entry is increased if present
     * \param[in]    forcefield    The force field containing the newest parameters
     * \param[in]    lossFunction  Functional form of the loss function
     */
    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               const std::vector<ACTMol *>   &actmols,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield,
                               LossFunction                   lossFunction);

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

    /*! \brief Minimize each molecule with an optimized reference structure,
     * perform a normal mode analysis and add the deviation of the
     * frequencies or intensities from the reference values to
     * the eRMS::FREQUENCY or eRMS::INTENSITY target
     * \param[in]    msghandler    For status and output
     * \param[in]    forceComputer The utility to compute forces and energies
     * \param[in]    actmols       The molecules to compute the deviation for
     * \param[inout] targets       Map of fitting targets, the frequency or intensity entry is increased
     * \param[in]    forcefield    The force field containing the newest parameters
     * \param[in]    lossFunction  Functional form of the loss function
     */
    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               const std::vector<ACTMol *>   &actmols,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield,
                               LossFunction                   lossFunction);
};

/*!
 * DevComputer that computes the deviation of the electrostatic potential -> eRMS::ESP
 */
class EspDevComputer : public DevComputer
{

private:

    /*! \brief Dump charges to a file
     * Debugging routine
     * \todo move to cpp file
     * \param[in] tw   The textwriter
     * \param[in] mol  The molecule to read from
     * \param[in] info Additional debugging information
     */
    void dumpQX(gmx::TextWriter              *tw,
                const ACTMol                 *mol,
                const std::string            &info);
   
public:

    //! Create an ESP deviation computer with current charge-distribution widths.
    EspDevComputer()
        : DevComputer("ESP")
    {
    }

    /*! \brief Compute the electrostatic potential for all molecules and add
     * the squared RMS deviation from the reference ESP to the eRMS::ESP target
     * \param[in]    msghandler    For status and debug output
     * \param[in]    forceComputer The utility to compute forces and energies, used to relax shells and vsites
     * \param[in]    actmols       The molecules to compute the deviation for
     * \param[inout] targets       Map of fitting targets, the eRMS::ESP entry is increased
     * \param[in]    forcefield    The force field containing the newest parameters
     * \param[in]    lossFunction  Not used
     */
    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               const std::vector<ACTMol *>   &actmols,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield,
                               LossFunction                   lossFunction);

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

    /*! \brief Compute the polarizability tensor for all molecules and add
     * the deviation from the reference tensor to the eRMS::Polar target
     * \param[in]    msghandler    For status and output
     * \param[in]    forceComputer The utility to compute forces and energies
     * \param[in]    actmols       The molecules to compute the deviation for
     * \param[inout] targets       Map of fitting targets, the eRMS::Polar entry is increased
     * \param[in]    forcefield    The force field containing the newest parameters
     * \param[in]    lossFunction  Functional form of the loss function
     */
    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               const std::vector<ACTMol *>   &actmols,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield,
                               LossFunction                   lossFunction);

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
     * @param mpo      The MolPropObservable, e.g. quadrupole
     */
    MultiPoleDevComputer(MolPropObservable mpo)
        : DevComputer(mpo_name(mpo)), mpo_(mpo)
    {
    }

    /*! \brief Compute the multipole moment for all molecules and add the
     * deviation from the reference to the corresponding target
     * (eRMS::MU, eRMS::QUAD, eRMS::OCT or eRMS::HEXADEC)
     * \param[in]    msghandler    For status and output
     * \param[in]    forceComputer The utility to compute forces and energies, used to relax shells and vsites
     * \param[in]    actmols       The molecules to compute the deviation for
     * \param[inout] targets       Map of fitting targets, the multipole entry is increased
     * \param[in]    forcefield    The force field containing the newest parameters
     * \param[in]    lossFunction  Functional form of the loss function
     */
    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               const std::vector<ACTMol *>   &actmols,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield,
                               LossFunction                   lossFunction);

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

    //! Whether calcDeviation computes the energies, otherwise the caller did, see energyTask
    bool        computeEnergies_ = true;

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

    //! The energy calculation that calcDeviation needs
    struct EnergyTask
    {
        //! Whether energies are needed at all
        bool compute                     = false;
        //! Whether the induction correction is separate, see ForceComputer::computeRelaxedEnergies
        bool separateInductionCorrection = true;
        //! Whether forces are needed
        bool forces                      = false;
    };

    /*! \brief The energy calculation that calcDeviation needs for these targets
     * Such that the caller can compute the energies of several sets of
     * molecules at once, see setComputeEnergies.
     * \param[in] targets The fitting targets
     * \return what to compute
     */
    EnergyTask energyTask(const std::map<eRMS, FittingTarget> &targets) const;

    /*! \brief The deviations of the energies from the energies as the GPU gives them
     * As calcDeviation after ForceComputer::computeRelaxedEnergies on the
     * GPU, without the energies of the conformation data: from
     * ForceComputer::gpuTypeEnergies, with the post-processing of the
     * energies done here. The losses are added to the targets in the same
     * order, such that the results are the same.
     * \param[in]    forceComputer The force computer, with the energies kept, see
     *                             ForceComputer::setKeepGpuEnergies
     * \param[in]    actmols       The molecules
     * \param[inout] targets       The fitting targets
     * \param[in]    lossFunction  The loss function
     * \return false if this cannot be done for these targets or molecules,
     *         then the targets may have changed
     */
    bool calcDeviationCompact(const ForceComputer           *forceComputer,
                              const std::vector<ACTMol *>   &actmols,
                              std::map<eRMS, FittingTarget> *targets,
                              LossFunction                   lossFunction);

    /*! \brief Whether calcDeviation computes the energies
     * \param[in] compute If false, the caller computed them as energyTask
     *                    says, with ForceComputer::computeRelaxedEnergies
     */
    void setComputeEnergies(bool compute) { computeEnergies_ = compute; }

    /*! \brief Create a new ForceEnergyDevComputer
     * \param[in] boltzmannTemperature Map from the different eRMS terms to a number
     *                                 indicating whether (>0) or not (0) Boltzmann weighting
     *                                 of energies and forces should be done.
     */
    ForceEnergyDevComputer(std::map<eRMS, double>  boltzmannTemperature);

    /*! \brief Compute the relaxed energies, and forces if needed, for all
     * molecules and add the deviation from the reference to the
     * eRMS::EPOT, eRMS::Force2 and interaction energy component targets
     * that have a non-zero weight
     * \param[in]    msghandler    For status and output
     * \param[in]    forceComputer The utility to compute forces and energies
     * \param[in]    actmols       The molecules to compute the deviation for
     * \param[inout] targets       Map of fitting targets, entries with non-zero weight are increased
     * \param[in]    forcefield    The force field containing the newest parameters
     * \param[in]    lossFunction  Functional form of the loss function
     */
    virtual void calcDeviation(MsgHandler                    *msghandler,
                               const ForceComputer           *forceComputer,
                               const std::vector<ACTMol *>   &actmols,
                               std::map<eRMS, FittingTarget> *targets,
                               const ForceField              *forcefield,
                               LossFunction                   lossFunction);

};

} // namespace alexandria

#endif //ALEXANDRIA_DEVCOMPUTER_H
