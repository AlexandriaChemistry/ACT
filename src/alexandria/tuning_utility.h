/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
 */

#ifndef TUNING_UTILITY_H
#define TUNING_UTILITY_H

#include <cstdio>

#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"

#include "mymol.h"
#include "molhandler.h"
#include "act/poldata/poldata.h"

/*! \brief Utility function to merge command line arguments
 * \param[inout] pargs The complete list of arguments
 * \param[in]    npa   The nmber of new elements
 * \param[in]    pa    The new elements
 */
void doAddOptions(std::vector<t_pargs> *pargs, size_t npa, t_pargs pa[]);

namespace alexandria
{

using qtStats = std::map<qType, gmx_stats>;


    /*! \brief Class to managa output from force field tuning
     */
    class TuneForceFieldPrinter
    {
    private:
        //! Can make computations for molecules, such as the Hessian and energy minimization (for now).
        MolHandler molHandler_;
        //! Tolerance (kJ/mol e) for marking ESP as an outlier in the log file
        real esp_toler_           = 30;
        //! Tolerance (Debye) for marking dipole as an outlier in the log file
        real dip_toler_           = 0.5;
        //! Tolerance (Buckingham) for marking quadrupole as an outlier in the log file
        real quad_toler_          = 5;
        //! Tolerance (Debye A^2) for marking octupole as an outlier in the log file
        real oct_toler_          = 5;
        //! Tolerance (Debye A^3) for marking hexadecapole as an outlier in the log file
        real hex_toler_          = 5;
        //! Tolerance (A^3) for marking diagonal elements of the polarizability tensor as an outlier in the log file
        real alpha_toler_         = 3;
        //! Tolerance (A^3) for marking isotropic polarizability as an outlier in the log file
        real isopol_toler_        = 2;
        //! Fit regression analysis of results to y = ax+b instead of y = ax
        bool useOffset_           = false;
        //! Perform energy minimization and compute vibrational frequencies for each molecule (after optimizing the force field if -optimize is enabled)
        bool calcFrequencies_     = true;

        //! \brief Analyse polarizability, add to statistics and print
        void analysePolarisability(FILE              *fp,
                                   alexandria::MyMol *mol,
                                   qtStats           *lsq_isoPol,
                                   qtStats           *lsq_anisoPol,
                                   qtStats           *lsq_alpha,
                                   real               efield);

        //! \brief And the atoms.
        void printAtoms(FILE              *fp,
                        alexandria::MyMol *mol);

        //! \brief do part of the printing, add to statistics
        void printEnergyForces(FILE                   *fp,
                               alexandria::MyMol      *mol,
                               const std::vector<int> &ePlot,
                               gmx_stats              *lsq_rmsf,
                               qtStats                *lsq_epot,
                               gmx_stats              *lsq_freq);
    public:
        TuneForceFieldPrinter() {}
    
        /*! \brief Add my options to the list of command line arguments
         * \param[out] pargs The vector to add to
         */
        void addOptions(std::vector<t_pargs> *pargs);
        
        /*! \brief Add my files to the list of command line arguments
         * \param[out] pargs The vector to add to
         */
        void addFileOptions(std::vector<t_filenm> *filenm);
    
        void print(FILE                           *fp,
                   std::vector<alexandria::MyMol> *mymol,
                   const Poldata                  *pd,
                   const gmx::MDLogger            &fplog,
                   const gmx_output_env_t         *oenv,
                   const CommunicationRecord      *cr,
                   real                            efield,
                   const std::vector<t_filenm>    &filenm);
    };
    
    /*! \brief Print header and command line arguments
     *
     * \param[in] fp    File pointer, if nullptr the function returns 
     *                  without doing anything
     * \param[in] pargs The command line arguments
     */
    void print_header(FILE                       *fp, 
                      const std::vector<t_pargs> &pargs);

} // namespace alexandria

#endif
