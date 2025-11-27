/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
 * \author Julian Marrades <julian@marrad.es>
 */

#ifndef ACT_MOLHANDLER_H
#define ACT_MOLHANDLER_H

#include <vector>

#include "act/forces/forcecomputer.h"
#include "act/utility/regression.h"
#include "actmol.h"
#include "confighandler.h"

namespace alexandria
{

//! \brief enum to distinguish result of minimizer
enum class eMinimizeStatus {
    OK, TooManySteps, Solver, NoMinimum
};

//! Convert enum to string
const std::string &eMinimizeStatusToString(eMinimizeStatus e);

/*! \brief Handles molecules by performing algorithms on them
 * For example, energy minimization and hessian computation.
 * This class can access private members of ACTMol as it is a friend class
 */
class MolHandler
{
public:

    /*! \brief Compute the second derivative matrix of the potential energy
     *
     * \param[in]  pd          Pointer to force field structure
     * \param[in]  mol        Molecule to get the hessian for
     * \param[in]  forceComp  Force Computer utility
     * \param[inout] coords   Atomic coordinates to operate on
     * \param[in]  atomIndex  Vector containing the indices of the real 
     *                        atoms, not shells or vsites.
     *                        FIXME: Create another method without this argument, then
     *                        get it and call this one
     * \param[out] hessian    MatrixWrapper object that must be pre-
     *                        allocated to NxN where N = 3*atomIndex.size()
     * \param[out] forceZero  The forces on the atoms in the input structure,
     *                        that is, not on the shells or vsites. Will be cleared
     *                        and overwritten
     * \param[out] energyZero The energies corresponding to the input structure.
     * \param[out] dpdq       Derivative of dipole moment with respect to atomic coordinates.
     *                        If nullptr it will not be used. This can be used to compute
     *                        infrared intensities from a NMA. See Henschel et al.
     *                        J. Chem. Theory Comput. 16 (2020) 3307-3315.
     */
    void computeHessian(const ForceField                  *pd,
                        const ACTMol                      *mol,
                        const ForceComputer               *forceComp,
                        std::vector<gmx::RVec>            *coords,
                        const std::vector<int>            &atomIndex,
                        MatrixWrapper                     *hessian,
                        std::vector<double>               *forceZero,
                        std::map<InteractionType, double> *energyZero,
                        std::vector<gmx::RVec>            *dpdq = nullptr) const;

    /*! \brief Perform normal-mode analysis on a molecule.
     * Computes vibrational frequencies and intensities and 
     * optionally prints them.
     *
     * Also prints eigenvalues and eigenvectors of the mass-weighted 
     * hessian matrix to the debug file, if not nullptr.
     * 
     * \param[in]  pd          Pointer to force field structure
     * \param[in]  mol          The molecule to analyze
     * \param[in]  forceComp    Force Computer utility
     * \param[in]  coords       Coordinates for a minimized structure
     * \param[out] frequencies  The normal mode frequencies (in cm^-1)
     * \param[out] intensities  The normal mode intensities
     * \param[out] output       Vector of string to write information such as frequencies to, may be nullptr (default)
     * \param[in]  debugNMA     Will provide excessive printing statements
     */
    void nma(const ForceField         *pd,
             const ACTMol             *mol,
             const ForceComputer      *forceComp,
             std::vector<gmx::RVec>   *coords,
             std::vector<double>      *frequencies,
             std::vector<double>      *intensities,
             std::vector<std::string> *output = nullptr,
             bool                      debugNMA = false) const;

    /*! \brief
     * The routine will energy minimize the atomic coordinates of a molecule while
     * relaxing the shells. The minimized coordinates will be stored in the 
     * actmol object. If the msForceToler = 0, the value will be derived from the
     * tolerance inside the forceComp, accoring to:
     * toler = 10*forceComp->rmsForce()^2
     *
     * \param[in] msghandler   For output and debugging
     * \param[in] pd           Pointer to force field structure
     * \param[in] mol          The molecule object (will be modified)
     * \param[in] forceComp    Force Computer utility
     * \param[in] simConfig    Configuration options
     * \param[inout] coords    The coordinates to be minimized
     * \param[in] energies     The energies per interaction type
     * \param[in] freeze       List of atoms (not shells) that will not be
     *                         moved during the minimization.
     * \return Status flag
     */
    eMinimizeStatus minimizeCoordinates(MsgHandler                        *msghandler,
                                        const ForceField                  *pd,
                                        const ACTMol                      *mol,
                                        const ForceComputer               *forceComp,
                                        const SimulationConfigHandler     &simConfig,
                                        std::vector<gmx::RVec>            *coords,
                                        std::map<InteractionType, double> *energies,
                                        const std::vector<int>            &freeze) const;
    /*! \brief
     * The routine will perform a MD simulation of a molecule or multiple
     * molecules, while relaxing shells if present.
     *
     * \param[in] msghandler     MsgHandler
     * \param[in] pd             Pointer to force field structure
     * \param[in] mol            The molecule object (will be modified)
     * \param[in] forceComp      Force Computer utility
     * \param[in] simConfig      Simulation configuration handler
     * \param[in] trajectoryFile Filename for writing coordinates
     * \param[in] energyFile     Filename for writing energies
     * \param[in] oenv           GROMACS output environment
     */
    void simulate(MsgHandler                    *msghandler,
                  const ForceField              *pd,
                  ACTMol                        *mol,
                  const ForceComputer           *forceComp,
                  const SimulationConfigHandler &simConfig,
                  const char                    *trajectoryFile,
                  const char                    *energyFile,
                  const gmx_output_env_t        *oenv) const;

    /*! \brief
     * The routine will compute the RMSD between two sets of coordinates
     * for a given actmol object.
     *
     * \param[in]    mol  The molecule object
     * \param[in]    xref The reference coordinate set.
     * \param[inout] xfit The (minimized) coordinate set before alignment (input) 
     *                    respectively after (outptu)
     * \return Root mean square atomic deviation of atomic
     *         coordinates after superposition.
     */
    double coordinateRmsd(const ACTMol                  *mol,
                          const std::vector<gmx::RVec> &xref,
                          std::vector<gmx::RVec>       *xfit) const;

};

} // namespace alexandria


#endif // ACT_MOLHANDLER_H
