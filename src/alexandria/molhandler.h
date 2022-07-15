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
 * \author Julian Marrades <julian.marrades@hotmail.es>
 */

#ifndef ACT_MOLHANDLER_H
#define ACT_MOLHANDLER_H

#include <vector>

#include "act/forces/forcecomputer.h"
#include "act/utility/regression.h"
#include "mymol.h"
#include "confighandler.h"
//#include "gromacs/mdtypes/commrec.h"

namespace alexandria
{

/*! \brief Handles molecules by performing algorithms on them
 * For example, energy minimization and hessian computation.
 * This class can access private members of MyMol as it is a friend class
 */
class MolHandler
{

public:

    /*! \brief Compute the second derivative matrix of the potential energy
     *
     * \param[in]  mol       Molecule to get the hessian for
     * \param[in]  forceComp Force Computer utility
     * \param[in]  crtmp     Temporary communication record for one core.
     *                       FIXME: another method without this
     * \param[in]  atomIndex Vector containing the indices of the real 
     *                       atoms, not shells or vsites.
     *                       FIXME: Create another method without this argument, then
     *                       get it and call this one
     * \param[out] hessian   MatrixWrapper object that must be pre-
     *                       allocated to NxN where N = 3*atomIndex.size()
     * \param[out] forceZero The forces on the atoms in the input structure,
     *                       that is, not on the shells or vsites. Will be cleared
     *                       and overwritten
     * \return the potential energy of the input structure
     */
    double computeHessian(      MyMol               *mol,
                          const ForceComputer       *forceComp,
                          const std::vector<int>    &atomIndex,
                                MatrixWrapper       *hessian,
                                std::vector<double> *forceZero) const;

    /*! \brief Perform normal-mode analysis on a molecule.
     * Computes vibrational frequencies and intensities and 
     * optionally prints them.
     *
     * Also prints eigenvalues and eigenvectors of the mass-weighted 
     * hessian matrix to the debug file, if not nullptr.
     * 
     * \param[in]  mol          The molecule to analyze
     * \param[in]  forceComp    Force Computer utility
     * \param[out] frequencies  The normal mode frequencies (in cm^-1)
     * \param[out] intensities  The normal mode intensities
     * \param[in]  fp           File to write frequencies to, may be nullptr (default)
     */
    void nma(MyMol               *mol,
             const ForceComputer *forceComp,
             std::vector<double> *frequencies,
             std::vector<double> *intensities,
             FILE                *fp = nullptr) const;

    /*! \brief
     * The routine will energy minimize the atomic coordinates of a molecule while
     * relaxing the shells. The minimized coordinates will be stored in the 
     * mymol object.
     *
     * \param[in] mol          The molecule object (will be modified)
     * \param[in] forceComp    Force Computer utility
     * \param[in] logFile      File to write some info to, may be a nullptr
     * \param[in] maxIter      Maximum number of iterations, 0 means until convergence
     * \param[in] overRelax    Factor to use for overrelaxation of the step size
     * \param[in] msForceToler Tolerance in the mean square force for convergence
     * \return Number of iterations used
     */
    int minimizeCoordinates(MyMol               *mol,
                            const ForceComputer *forceComp,
                            FILE                *logFile,
                            int                  maxIter,
                            double               overRelax,
                            double               msForceToler) const;

    /*! \brief
     * The routine will perform a MD simulation of a molecule or multiple
     * molecules, while relaxing shells if present.
     *
     * \param[in] mol            The molecule object (will be modified)
     * \param[in] forceComp      Force Computer utility
     * \param[in] simConfig      Simulation configuration handler
     * \param[in] logFile        File for logging information
     * \param[in] trajectoryFile Filename for writing coordinates
     * \param[in] energyFile     Filename for writing energies
     * \param[in] oenv           GROMACS output environment
     */
    void simulate(MyMol                         *mol,
                  const ForceComputer           *forceComp,
                  const SimulationConfigHandler &simConfig,
                  FILE                          *logFile,
                  const char                    *trajectoryFile,
                  const char                    *energyFile,
                  const gmx_output_env_t        *oenv) const;

    /*! \brief
     * The routine will compute the RMSD between the minimized coordinates
     * and the original ones for a given mymol object.
     * This routine should be called only after minimizing the coordinates,
     * since otherwise there is no minimized structure.
     *
     * \param[in]  mol  the molecule object (will be modified)
     * \param[out] x    The two coordinate sets after alignment. Map will be cleared first.
     * \return Root mean square atomic deviation of atomic
     *         coordinates after minimization.
     */
    double coordinateRmsd(MyMol                                       *mol,
                          std::map<coordSet, std::vector<gmx::RVec> > *x) const;

};

} // namespace alexandria


#endif // ACT_MOLHANDLER_H
