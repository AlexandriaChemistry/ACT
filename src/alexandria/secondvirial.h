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

/*! \brief Print a force field summary
 * \param[input] jtree JSon tree
 * \param[in]    pd    The force field
 */
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
     * \param[in]  logFile For debugging info, may be a nullptr
     * \param[in]  mymol   The description of the two fragments
     * \param[out] coords  The coordinate sets
     */
    void generate(FILE                                *logFile,
                  const MyMol                         *mymol,
                  std::vector<std::vector<gmx::RVec>> *coords);
};

class ReRunner
{
private:
    //! The force computer
    ForceComputer      *forceComp_ = nullptr;
    //! The dimer generator
    DimerGenerator     *gendimers_ = nullptr;
    //! GROMACS output stuff
    gmx_output_env_t   *oenv_      = nullptr;
               
    //! Trajectory name
    const char         *trajname_  = "";
    //! Temperature to start
    double              T1_        = 300;
    //! Final temperature
    double              T2_        = 400;
    //! Temperature step
    double              deltaT_    = 10;
    //! Whether to compute interaction energies and potentially B2
    bool                eInter_    = true;
    //! The temperature array
    std::vector<double> Temperatures_;
    //! Second virial as a function of T (see function temperatures)
    std::vector<double> b2t_;

    //! Generate temperature series if needed and return it
    const std::vector<double> &temperatures();
    /*! \brief Generate plot with Mayer functions for all temperatures
     * \param[in] ehisto The output file name
     * \param[in] mayer  The curves
     */
    void plotMayer(const char                             *ehisto,
                   const std::vector<std::vector<double>> &mayer);

    /*! \brief Generate plot with Second virial as a function of temperature
     * \param[in] b2file The B2(T) output file name
     */
    void plotB2temp(const char *b2file);


public:
    ReRunner() {}
    
    /*! \brief Make copies of utilities
     * \param[in] forceComp The force computer
     * \param[in] gendimers Dimer generator
     * \param[in] oenv      GROMACS plotting stuff
     */
    void setFunctions(ForceComputer    *forceComp,
                      DimerGenerator   *gendimers,
                      gmx_output_env_t *oenv)
    {
        forceComp_ = forceComp;
        gendimers_ = gendimers;
        oenv_      = oenv;
    }
    
    /*! \brief Add command line options
     * \param[inout] pargs  Regular flags
     * \param[inout] filenm File options
     */
    void addOptions(std::vector<t_pargs>  *pargs,
                    std::vector<t_filenm> *filenm);
    
    //! \return whether or not we will compute interaction energies                
    bool eInteraction() const { return eInter_; }
    
    //! \brief Manually set the temperatures
    void setTemperatures(double T1, double T2, double deltaT) { T1_ = T1; T2_ = T2; deltaT_ = deltaT; }
    
    //! \brief Manually set a temperature array
    void setTemperatures(const std::vector<double> &T) { Temperatures_ = T; }
    
    //! \brief Reset temperatures
    void resetTemperatures() { Temperatures_.clear(); }
    
    /*! \brief Set the interaction energy flag
     * \param[in] eInter The value
     */
    void setEInteraction(bool eInter) { eInter_ = eInter; }
    
    /*! \brief Do the rerunning with different options
     * \param[in] logFile File pointer to print info
     * \param[in] pd      Force field structure
     * \param[in] mymol   Structure with molecule info
     * \param[in] qtot    The total charge for when reading from a trajectory
     * \param[in] verbose Whether or not to print a lot
     * \param[in] fnm     The filenames
     */
    void rerun(FILE                        *logFile,
               const Poldata               *pd,
               const MyMol                 *mymol,
               double                       qtot,
               bool                         verbose,
               const std::vector<t_filenm> &fnm);

    /*! \brief Compute the second virial coefficient including QM corrections
     * \param[in] logFile   Output file for printing
     * \param[in] edist     Statistics for interaction energies
     * \param[in] inertia   The moments of inertia of the molecules
     * \param[in] force1    The interaction forces on molecule 1
     * \param[in] torqueMol The torque on both molecules
     * \param[in] fnm       The filenames
     */
    void computeB2(FILE                                      *logFile,
                   gmx_stats                                  edist,
                   double                                     mass,
                   const gmx::RVec                            inertia[2],
                   const std::vector<gmx::RVec>              &force1,
                   const std::vector<std::vector<gmx::RVec>> &torqueMol,
                   const std::vector<t_filenm>               &fnm);
                   
    //! \return the second virial as a function of T.
    const std::vector<double> &b2Temp() const { return b2t_; }

};

} // namespace alexandria

#endif // ACT_SIMULATE_H
