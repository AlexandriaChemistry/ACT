/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2023
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
#ifndef ACT_FORCECOMPUTER_H
#define ACT_FORCECOMPUTER_H
    
#include <vector>

#include "act/alexandria/topology.h"
#include "act/poldata/poldata.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class QtypeProps;

/*! \brief Class to compute all the forces in a molecule or complex
 */
class ForceComputer
{
private:
    //! Convergence criterium for minimizing shells: mean square force
    double         msForce_;
    //! Maximum number of iterations to spend on minimizing shells
    int            maxiter_;
    //! Electric field to be optionally applied
    gmx::RVec      field_ = { 0.0, 0.0, 0.0 };
    /*! Do one actual computations.
     * Will do one force/energy computation.
     * \param[in]  pd          The force field structure
     * \param[in]  top         The molecular topology
     * \param[in]  coordinates The atomic coordinates. Coordinates of
     *                         shell particles may be changed.
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     * \param[in]  field       Optional electric field to be applied
     */
    void computeOnce(const Poldata                     *pd,
                     const Topology                    *top,
                     std::vector<gmx::RVec>            *coordinates,
                     std::vector<gmx::RVec>            *forces,
                     std::map<InteractionType, double> *energies,
                     const gmx::RVec                   &field) const;

 public:
    /*! \brief Constructor
     * \param[in] msForce The tolerance for the mean square force on shells
     * \param[in] maxiter The maximum number of iterations for shell minimization
     */
    ForceComputer(double   msForce = 1e-6,
                  int      maxiter = 25) : msForce_(msForce), maxiter_(maxiter) {}
    
    /*! Do complete energy/force computation.
     * If shells are present their positions will be minimized.
     * \param[in]  pd          Pointer to force field structure
     * \param[in]  top         The molecular topology
     * \param[in]  charges     The charges for all particles
     * \param[in]  coordinates The atomic coordinates. Coordinates of
     *                         shell particles may be changed.
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     * \param[in]  field       Optional electric field to be applied
     * \return The mean square force on the shells, or zero if not present.
     */
    double compute(const Poldata                     *pd,
                   const Topology                    *top,
                   std::vector<gmx::RVec>            *coordinates,
                   std::vector<gmx::RVec>            *forces,
                   std::map<InteractionType, double> *energies,
                   const gmx::RVec                   &field = { 0.0, 0.0, 0.0 }) const;
                 
    /*! \brief Return the gromacs type used
     * In practice this converts the InteractionType to the ftype
     * used within the force computer.
     * \param[in] pd      Pointer to force field structure
     * \param[in] itype The interaction type
     * \return the force type
     */
    int ftype(const Poldata  *pd,
              InteractionType itype) const;
    
    //! \return the force tolerance
    double forceTolerance() const { return msForce_; }

    /*! \brief Set the shell minimization RMS tolerance
     * \param[in] toler The new tolerance
     */
    void setForceTolerance(double toler) { msForce_ = toler; }
    
    /*! \brief Plot the potential functions
     * This plots the potential functions corresponding to
     * InteractionType.
     * \param[in] pd    Pointer to force field structure
     * \param[in] itype The interaction type
     */
    void plot(const Poldata   *pd,
              InteractionType  itype) const;

    /*! \brief Compute the polarizability tensor
     * \param[in]  pd          Pointer to force field structure
     * \param[in]  top         The topology
     * \param[in]  coordinates The coordinates
     * \param[out] qtp         The charge type properties
     */
    void calcPolarizability(const Poldata          *pd,
                            const Topology         *top,
                            std::vector<gmx::RVec> *coordinates,
                            QtypeProps             *qtp) const;

};

} // namespace alexandria

#endif
