/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024-2026
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
#ifndef ACT_ALEXANDRIA_B2DATA_H
#define ACT_ALEXANDRIA_B2DATA_H

#include <vector>

#include "act/basics/msg_handler.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

/*! Temporary arrays for weighted properties.
 * For the two dimensional arrays, the first is the temperature index
 */
class B2Data
{
private:
    //! Mayer-weighted energy.
    std::vector<std::vector<double>>    exp_U12_;
    //! Forces on two compounds
    std::vector<std::vector<double>>    exp_F2_[2];
    //! Torque on two compounds
    std::vector<std::vector<gmx::RVec>> exp_tau_[2];
    //! Sum of the weights of the entries per bin
    std::vector<std::vector<double>>    n_U12_;
    //! Mayer functions
    std::vector<std::vector<double>>    mayer_;
    //! Distance scale
    std::vector<double>                 dist_;
    //! Copy of temperature array
    std::vector<double>                 temperatures_;
public:
    /*! \brief Constructor
     * \param[in] nbins        Number of bins
     * \param[in] temperatures Ts at which to compute the second virial
     */
    B2Data(int                        nbins,
           const std::vector<double> &temperatures);

    /*! \brief Dump some info using a MsgHandler
     * \param[in] msghandler The message handler
     */
    void dump(MsgHandler *msghandler) const;

    /*! \brief Aggregate information from helpers to tha main node
     * \param[in] cr Communication information 
     */
    void aggregate(CommunicationRecord *cr);

    /*! \brief Add one data point.
     * Adds all values times the weight, and the weight to the internal
     * counter, such that the averages are weighted.
     * \param[in] iTemp    Temperature index
     * \param[in] index    Distance index
     * \param[in] dist     The distance between molecules
     * \param[in] exp_U12  Exp of the energy
     * \param[in] exp_F0   Exp weighted Force on molecule 0
     * \param[in] exp_F1   Exp weighted Force on molecule 1
     * \param[in] exp_tau0 Exp weighted Torque on molecule 0
     * \param[in] exp_tau1 Exp weighted Torque on molecule 1
     * \param[in] weight   Weight of the orientation, see DimerGenerator
     */
    void addData(size_t iTemp, size_t index, double dist,
                 double exp_U12, double exp_F0, double exp_F1,
                 const gmx::RVec exp_tau0, const gmx::RVec exp_tau1,
                 double weight = 1);

    /*! \brief Add the dimers of one orientation at all temperatures.
     * Gives the same sums as addData for each dimer and temperature, with
     * the Boltzmann factor g = exp(-beta U) of the energy U of each dimer,
     * g-1 as Mayer function, g times the squared forces and g times the
     * squared torque components, see Gray and Gubbins Eqns. 3.261, 3.272
     * and 3.282. Temperatures of zero are skipped.
     * \param[in] index  Distance index of each dimer
     * \param[in] dist   Distance between the molecules of each dimer
     * \param[in] energy Energy of each dimer
     * \param[in] force2 Squared force on each molecule of each dimer
     * \param[in] torque Torque on each molecule of each dimer
     * \param[in] weight Weight of each dimer, see DimerGenerator
     */
    void addOrientation(const std::vector<size_t>    &index,
                        const std::vector<double>    &dist,
                        const std::vector<double>    &energy,
                        const std::vector<double>     force2[2],
                        const std::vector<gmx::RVec>  torque[2],
                        const std::vector<double>    &weight);

    /*! \brief Do the integration of the exponent weight functions
     * Compute classical B2 and quantum correction. The bins hold equidistant
     * distances, starting at the shortest distance of the dimers. Below it
     * the Mayer function is assumed to be -1, which gives a hard core
     * contribution to the classical B2.
     * \param[in]  msghandler For debugging
     * \param[in]  iTemp      Temperature index
     * \param[in]  beta       Boltzmann factor
     * \param[in]  mass       The atomic masses
     * \param[in]  inertia    Moments of inertia of the two compounds
     * \param[out] Bclass     Classical B2
     * \param[out] BqmForce   Force contribution to B2 (first order only)
     * \param[out] BqmTorque1 Contribution due to torque on compound 1
     * \param[out] BqmTorque2 Contribution due to torque on compound 2
     */
    void integrate(MsgHandler                   *msghandler,
                   int                           iTemp,
                   double                        beta,
                   const std::vector<double>    &mass,
                   const std::vector<gmx::RVec> &inertia,
                   double                       *Bclass,
                   double                       *BqmForce,
                   double                       *BqmTorque1,
                   double                       *BqmTorque2);

    /*! \brief Generate plot with Mayer functions for all temperatures
     * \param[in] fmayer The output file name
     * \param[in] oenv   GROMACS output environment
     */
    void plotMayer(const char                *fmayer,
                   gmx_output_env_t          *oenv);
};

} // namespace

#endif
