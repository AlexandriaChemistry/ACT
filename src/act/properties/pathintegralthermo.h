/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ACT_PATHINTEGRALTHERMO_H
#define ACT_PATHINTEGRALTHERMO_H

#include <cstddef>
#include <functional>
#include <vector>

#include "act/properties/pathintegralmd.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class ACTMol;
class ForceComputer;
class ForceField;
class JsonTree;
class MsgHandler;
class SimulationConfigHandler;

/*! \brief The points and weights of Gauss-Legendre quadrature between zero and one
 * \param[in]  n       The number of points
 * \param[out] nodes   The points, in increasing order
 * \param[out] weights Their weights, which add up to one
 * \throws gmx::InvalidInputError if n is less than one
 */
void gaussLegendre(int n, std::vector<double> *nodes, std::vector<double> *weights);

/*! \brief The energy of harmonic oscillators, as ring polymers of beads give it
 *
 * For P beads the energy of an oscillator with the frequency \f$\nu\f$ is
 * \f[ E_P = k_BT\sum_{k=0}^{P-1} r_k,\quad
 * r_k = \frac{(gh\nu)^2}{(gh\nu)^2 + (4\pi P k_BT\sin(\pi k/P))^2} \f]
 * with the constant of Planck scaled by g, half of it kinetic and half
 * potential, and its heat capacity is \f$k_B\sum_k r_k(2r_k-1)\f$. For
 * infinitely many beads these are the energy
 * \f$\frac{gh\nu}{2}\coth x\f$ and the heat capacity
 * \f$k_B(x/\sinh x)^2\f$ of the quantum oscillator, with
 * \f$x = gh\nu/2k_BT\f$. One bead gives the classical oscillator.
 * \param[in]  frequencies  The frequencies in 1/ps, those that are not positive are left out
 * \param[in]  temperature  The temperature in K
 * \param[in]  nbeads       The number of beads, zero for infinitely many
 * \param[in]  planckScale  The factor g for the constant of Planck
 * \param[out] heatCapacity The heat capacity in kJ/mol K, may be nullptr
 * \return the sum of the energies in kJ/mol
 */
double harmonicBeadEnergy(const std::vector<double> &frequencies,
                          double                     temperature,
                          int                        nbeads,
                          double                     planckScale,
                          double                    *heatCapacity = nullptr);

/*! \brief The quantum part of the free energy of harmonic oscillators
 *
 * \f[ k_BT\sum_i \ln\frac{\sinh x_i}{x_i},\quad x_i = \frac{h\nu_i}{2k_BT} \f]
 * the free energy minus the classical one.
 * \param[in] frequencies The frequencies in 1/ps
 * \param[in] temperature The temperature in K
 * \return the free energy in kJ/mol
 */
double harmonicQuantumFreeEnergy(const std::vector<double> &frequencies,
                                 double                     temperature);

/*! \brief What the normal mode analysis of the minimized compound gives
 *
 * The thermochemistry of the rigid rotor and harmonic oscillators, see
 * ThermoChemistry, which is what the simulations are compared to and
 * where the integrals for the entropy start.
 */
struct HarmonicReference
{
    //! The potential energy of the minimized structure in kJ/mol
    double              epot               = 0;
    //! The frequencies of the vibrations in 1/ps
    std::vector<double> frequencies;
    //! The number of atoms
    size_t              natom              = 0;
    //! Whether the compound is linear
    bool                linear             = false;
    //! The rotational symmetry number
    int                 symmetryNumber     = 1;
    //! The entropy of translation at 1 bar in J/mol K
    double              entropyTranslation = 0;
    //! The entropy of rotation of the classical rigid rotor in J/mol K
    double              entropyRotation    = 0;
    //! The zero point energy in kJ/mol
    double              zeroPointEnergy    = 0;
    //! The internal energy in kJ/mol relative to the minimum, with the zero point energy
    double              internalEnergy     = 0;
    //! The heat capacity at constant volume in J/mol K
    double              heatCapacity       = 0;
    //! The standard entropy in J/mol K
    double              entropy            = 0;

    //! \return the number of rotations: none for an atom, two for a linear compound
    int nrotation() const { return natom < 2 ? 0 : (linear ? 2 : 3); }

    //! \return the number of vibrations
    int nvibration() const { return 3*static_cast<int>(natom) - 3 - nrotation(); }

    /*! \brief The internal energy of the classical rigid rotor and harmonic oscillators
     * \param[in] temperature The temperature in K
     * \return the energy in kJ/mol relative to the minimum, with the translation
     */
    double classicalEnergy(double temperature) const;

    /*! \brief The entropy of the classical rigid rotor and harmonic oscillators
     * \param[in] temperature The temperature in K, the one of the entropy
     *                        of translation and rotation
     * \return the standard entropy in J/mol K
     */
    double classicalEntropy(double temperature) const;
};

/*! \brief Normal mode analysis of a minimized compound for the thermochemistry
 * \param[in] msghandler  For output and warnings
 * \param[in] pd          The force field
 * \param[in] actmol      The compound
 * \param[in] forceComp   The force computer
 * \param[in] xmin        The minimized coordinates
 * \param[in] temperature The temperature in K
 * \return the reference
 * \throws gmx::InvalidInputError if the structure is not a minimum, that
 * is if a vibration has no frequency
 */
HarmonicReference harmonicReference(MsgHandler                   *msghandler,
                                    const ForceField             *pd,
                                    const ACTMol                 *actmol,
                                    const ForceComputer          *forceComp,
                                    const std::vector<gmx::RVec> &xmin,
                                    double                        temperature);

/*! \brief The simulations of the two integrals that give the entropy
 *
 * Both integrals use the same points and weights between zero and one.
 */
struct ThermoIntegration
{
    //! The points of the quadrature
    std::vector<double>       nodes;
    //! Their weights
    std::vector<double>       weights;
    //! The potential energy in kJ/mol of classical simulations, with one
    //! bead, at the temperature times each point
    std::vector<AverageError> classicalEpot;
    //! The kinetic energy in kJ/mol of simulations with beads at the
    //! temperature, with the constant of Planck times each point
    std::vector<AverageError> quantumEkin;
    //! The number of beads of each of these simulations
    std::vector<int>          quantumBeads;
};

//! The thermochemistry of a compound from simulations with beads, see pathIntegralThermochemistry
struct PathIntegralThermo
{
    //! The temperature in K
    double            temperature     = 0;
    //! The number of beads
    int               nbeads          = 0;
    //! The harmonic values
    HarmonicReference harmonic;
    //! The internal energy in kJ/mol relative to the minimum
    AverageError      internalEnergy;
    //! The enthalpy in kJ/mol relative to the minimum, the internal energy plus k_B T
    AverageError      enthalpy;
    //! The heat capacity at constant volume in J/mol K
    AverageError      heatCapacity;
    //! Whether there is an entropy
    bool              haveEntropy     = false;
    //! The standard entropy in J/mol K
    AverageError      entropy;
    //! The Gibbs energy in kJ/mol relative to the minimum, the enthalpy minus T times the entropy
    AverageError      gibbsEnergy;
    //! The change of the free energy in kJ/mol due to the anharmonicity of the classical compound
    AverageError      anharmonic;
    //! The free energy in kJ/mol of the quantum compound minus that of the classical one
    AverageError      quantum;
    //! The same for harmonic oscillators with the frequencies of the normal modes
    double            quantumHarmonic = 0;
    //! The correction in kJ/mol of the internal energy for the finite number of beads, which is included
    double            beadEnergy       = 0;
    //! The correction in J/mol K of the heat capacity for the finite number of beads, which is included
    double            beadHeatCapacity = 0;
};

/*! \brief Thermochemistry from the averages of simulations with beads
 *
 * The internal energy \f$E\f$ and the heat capacity are those of the
 * simulation at the temperature \f$T\f$, corrected for the finite number
 * of beads: the difference between harmonic oscillators with infinitely
 * many beads and with the beads of the simulation is added, see
 * harmonicBeadEnergy. For water at room temperature with 32 beads this is
 * 1.8 kJ/mol and -11 J/mol K, what remains is the effect of the number of
 * beads on the anharmonicity. The entropy is
 * \f[ S = S_{cl}^{RRHO} + \frac{E - E_{cl}^{RRHO} - \Delta F_{anh} - \Delta F_q}{T} \f]
 * where RRHO is the classical rigid rotor with harmonic oscillators,
 * \f[ \Delta F_{anh} = -T\int_0^T \frac{E_{cl}(T') - E_{cl}^{RRHO}(T')}{T'^2}dT' \f]
 * is the change due to the anharmonicity of the classical compound, from
 * classical simulations at lower temperatures, where the compound becomes
 * harmonic, and
 * \f[ \Delta F_q = \int_0^1 \frac{2}{g}\left(\langle K\rangle_g - \frac{3N}{2}k_BT\right)dg \f]
 * is the change from the classical to the quantum compound, from
 * simulations in which the constant of Planck is \f$g\f$ times its value,
 * which is the same as atoms that are heavier by \f$1/g^2\f$ and needs
 * fewer beads by the factor \f$g\f$. Of the second integrand the value
 * for harmonic oscillators with the same beads is taken off and the
 * integral for infinitely many beads, which is known, added, such that
 * few points are needed and the number of beads is corrected for here as
 * well.
 * \param[in] harmonic    The normal mode analysis
 * \param[in] main        The simulation at the temperature with all beads
 * \param[in] integration The simulations of the integrals, nullptr for no entropy
 * \return the thermochemistry, with standard errors from those of the simulations
 * \throws gmx::InvalidInputError if the simulations of the integrals are
 * not there for all points
 */
PathIntegralThermo computeThermochemistry(const HarmonicReference   &harmonic,
                                          const PathIntegralResults &main,
                                          const ThermoIntegration   *integration);

/*! \brief Thermochemistry of a compound from simulations with beads
 *
 * Does the simulations of the two integrals of the entropy, see
 * computeThermochemistry, without files. The classical simulations at a
 * lower temperature are longer by the ratio of the temperatures, at most
 * sixteen times, since their energies are divided by that temperature
 * squared. A simulation in which the constant of Planck is g times its
 * value has g times the beads, rounded down, and is left out if that is
 * one bead, which is classical.
 * \param[in] msghandler For output and warnings
 * \param[in] pd         The force field
 * \param[in] actmol     The compound, one fragment
 * \param[in] xmin       The minimized coordinates of all particles
 * \param[in] forceComp  The force computer
 * \param[in] simConfig  The settings of the simulation with beads, with
 *                       the number of simulations for each integral
 * \param[in] harmonic   The normal mode analysis of the minimized compound
 * \param[in] main       The simulation at the temperature with all beads
 * \param[in] progress   Called after each simulation with the number done
 *                       and the total number, may be empty
 * \return the thermochemistry
 */
PathIntegralThermo pathIntegralThermochemistry(MsgHandler                         *msghandler,
                                               const ForceField                   *pd,
                                               const ACTMol                       *actmol,
                                               const std::vector<gmx::RVec>       &xmin,
                                               const ForceComputer                *forceComp,
                                               const SimulationConfigHandler      &simConfig,
                                               const HarmonicReference            &harmonic,
                                               const PathIntegralResults          &main,
                                               const std::function<void(int, int)> &progress = {});

/*! \brief Write the thermochemistry from simulations with beads
 * \param[in]    msghandler For the log file
 * \param[in]    thermo     The thermochemistry
 * \param[inout] jtree      It is added to it, may be nullptr
 */
void reportThermochemistry(MsgHandler               *msghandler,
                           const PathIntegralThermo &thermo,
                           JsonTree                 *jtree);

} // namespace alexandria

#endif
