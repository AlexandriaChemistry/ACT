/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024
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

#include <cstdio>
    
#include <vector>

#include "act/utility/communicationrecord.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

/*! \brief Compute integral over sphere section
 * \param[in] r1   The radius to start at
 * \param[in] r2   The radius to stop at
 * \param[in] val1 The function value at r1
 * \param[in] val2 The function value at r2
 * \return The integral of 4 pi r^2 f(r) from r1 to r2
 */
double sphereIntegrator(double r1, double r2, double val1, double val2);

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
    //! Number of entries per bin
    std::vector<std::vector<int>>       n_U12_;
    //! Mayer functions
    std::vector<std::vector<double>>    mayer_;
    //! Distance scale
    std::vector<double>                 dist_;
    //! Copy of temperature array
    std::vector<double>                 temperatures_;
public:
    B2Data(int                        nbins,
           double                     binwidth,
           const std::vector<double> &temperatures);

    /*! \brief Dump some info to a file
     * \param[in] fp The file pointer
     */
    void dump(FILE *fp) const;

    /*! \brief Aggregate information from helpers to tha main node
     * \param[in] cr Communication information 
     */
    void aggregate(CommunicationRecord *cr);

    /*! \brief Add one data point.
     * Adds all values and increments internal counter.
     * \param[in] iTemp    Temperature index
     * \param[in] index    Distance index
     * \param[in] exp_U12  Exp of the energy
     * \param[in] exp_F0   Exp weighted Force on molecule 0
     * \param[in] exp_F1   Exp weighted Force on molecule 1
     * \param[in] exp_tau0 Exp weighted Torque on molecule 0
     * \param[in] exp_tau1 Exp weighted Torque on molecule 1
     */
    void addData(size_t iTemp, size_t index,
                 double exp_U12, double exp_F0, double exp_F1,
                 gmx::RVec exp_tau0, gmx::RVec exp_tau1);

    /*! \brief  Fill the Mayer function with -1 values
     * \param[in] iTemp    Temperature index
     * \param[in] xmin     First non-zero energy
     * \param[in] binWidth The bin width in the arrays
     */
    void fillToXmin(int iTemp, double xmin, double binWidth);

    /*! \brief Do the integration of the exponent weight functions
     * Compute classical B2 and quantum correction.
     * \param[in]  iTemp    Temperature index
     * \param[in]  binWidth The bin width in the arrays
     * \param[in]  beta     Boltzmann factor
     * \param[in]  inertia  Moments of inertia of the two compounds
     * \param[out] Bclass     Classical B2
     * \param[out] BqmForce   Force contribution to B2 (first order only)
     * \param[out] BqmTorque1 Contribution due to torque on compound 1
     * \param[out] BqmTorque2 Contribution due to torque on compound 2
     */
    void integrate(int iTemp, double binWidth, double beta,
                   const std::vector<double>    &mass,
                   const std::vector<gmx::RVec> &inertia,
                   double *Bclass, double *BqmForce,
                   double *BqmTorque1, double *BqmTorque2);

    /*! \brief Generate plot with Mayer functions for all temperatures
     * \param[in] fmayer       The output file name
     * \param[in] temperatures The T used
     */
    void plotMayer(const char                *fmayer,
                   gmx_output_env_t          *oenv);
};

} // namespace

#endif
