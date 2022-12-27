/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

#ifndef ACT_SIMULATE_H
#define ACT_SIMULATE_H

#include <vector>

#include "act/forces/forcecomputer.h"
#include "act/poldata/poldata.h"
#include "act/utility/jsontree.h"
#include "alexandria/mymol.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/statistics/statistics.h"

namespace alexandria
{

class DimerGenerator
{
private:
    //! Number of dimers to generate
    int maxdimers_ = 1000;
    //! Number of distances to use
    int ndist_      = 20;
    //! Minimum com-com distance
    double mindist_ = 0.2;
    //! Maximum com-com distance
    double maxdist_ = 2.0;
    //! Random number seed
    int    seed_    = 1993;    
public:
    //! Constructor
    DimerGenerator() {}
    
    /*! \brief Add my options
     * \param[inout] pa The command line options
     */
    void addOptions(std::vector<t_pargs> *pa);
    
    //! Return the number of distancea
    int ndist() const { return ndist_; }
    
    /*! \brief Do the actual generation
     * \param[in] mymol The description of the two fragments
     * \param[out] mps  The generated dimers
     */
    void generate(const MyMol          *mymol,
                  std::vector<MolProp> *mps);
};


void forceFieldSummary(JsonTree      *jtree,
                       const Poldata *pd);

/*! \brief Compute integral over sphere section
 * \param[in] r1   The radius to start at
 * \param[in] r2   The radius to stop at
 * \param[in] val1 The function value at r1
 * \param[in] val2 The function value at r2
 * \return The integral of 4 pi r^2 f(r) from r1 to r2
 */
double sphereIntegrator(double r1, double r2, double val1, double val2);

void computeB2(FILE                         *logFile,
               const char                   *ehisto,
               gmx_stats                     edist,
               gmx_output_env_t             *oenv,
               const std::vector<double>    &Temperature,
               double                        mass,
               const gmx::RVec              &inertia,
               const std::vector<gmx::RVec> &force1,
               const std::vector<gmx::RVec> &torque1,
               std::vector<double>          *b2t);

void do_rerun(FILE                      *logFile,
              const Poldata             *pd,
              const MyMol               *mymol,
              ForceComputer             *forceComp,
              DimerGenerator            *gendimers,
              const char                *trajname,
              const char                *ehisto,
              const char                *b2file,
              bool                       eInter,
              double                     qtot,
              gmx_output_env_t          *oenv,
              const std::vector<double> &Temperature);
              
} // namespace alexandria

#endif // ACT_SIMULATE_H
