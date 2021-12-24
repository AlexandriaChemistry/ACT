/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021
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
 * In particular, this is the core of the bastat program which
 * does an analysis of bond and angle statistics.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ALLBONDEDS_H
#define ALLBONDEDS_H

#include <map>

#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/real.h"

#include "mymol.h"
#include "poldata.h"

namespace alexandria
{

    /*! \brief Base class for managing bonds, angles and dihedrals
     */
    class OneBonded
    {
    private:
        //! Atoms and bondorder
        Identifier id_;
        //! LSQ fitting structure
        gmx_stats  lsq_;
    public:
        /*! Constructor
         * \param[in] id The identifier
         */
        OneBonded(const Identifier id) : id_(id) {}
    
        //! Return my id
        const Identifier &id() const { return id_; }
    
        /*! Add data point to histogram
         * \param[in] x The point to add
         */
        void addPoint(double x);
    
        /*! Write a histogram
         * \param[in] fn
         * \param[in] xaxis
         * \param[in] oenv
         * \param[in] spacing
         */
        void writeHistogram(const char             *fn_prefix,
                            const char             *xaxis,
                            const gmx_output_env_t *oenv,
                            double                  binwidth);
        /*! \brief Extract statistics from this parameter
         * TODO Implement support for extracting the median value
         * \param[out] average The average value    
         * \param[out] sigma   The standard deviation
         * \param[out] N       The number of data points
         * \return status from statistics library
         */
        eStats getAverageSigmaN(real   *average,
                                real   *sigma,
                                size_t *N);
    };

    /*! Class to analyze structures and derive geometric statistics
     */
    class AllBondeds
    {
    private:
        //! Map from interaction type to list of bondeds
        std::map<InteractionType, std::vector<OneBonded> > bondeds_;
        //! Default dissociation energy
        real                                               Dm_        = 400;
        //! Default steepness for Morse potentials
        real                                               beta_      = 20;
        //! Default force constant for harmonic angles
        real                                               kt_        = 300;
        //! Default force constant for harmonic linear angles
        real                                               klin_      = 150000;
        //! Default force constant for dihedral angles
        real                                               kp_        = 5;
        //! Default force constant for improper dihedrals
        real                                               kimp_      = 1;
        //! Default force constant for Urey-Bradley
        real                                               kub_       = 30000;
        //! Tolerance for warning about large sigma in bond-lengths (pm)
        real                                               bond_tol_  = 5;
        //! Tolerance for warning about large sigma in angles (degrees)
        real                                               angle_tol_ = 5;
        //! Scaling factor for setting min and max for force parameters
        real                                               factor_    = 0.8;
        //! Spacing for bond histograms in pm
        real                                               bspacing_  = 1;
        //! Spacing for angle histograms in degrees
        real                                               aspacing_  = 0.5;
        //! Relative number for linear angles
        real                                               laspacing_ = 0.000001;
        //! Spacing for dihedrals in degrees
        real                                               dspacing_  = 1;
        /*! \brief Add bonds etc. for one molecule and list of atoms
         * \param[in] fplog   File to print information
         * \param[in] pd      Force field structure
         * \param[in] iType   InteractionType
         * \param[in] mmi     Molecule structure
         * \param[in] atomid  List of atoms involved in the interaction
         * \param[in] canSwap Whether or not the atom order can be swapped
         */
        void addBonded(FILE                           *fplog, 
                       const Poldata                  &pd,
                       InteractionType                 iType,
                       const MyMol                    &mmi,
                       const std::vector<int>         &atomid,
                       CanSwap                         canSwap);

    public:
        //! Constructor
        AllBondeds() {}
    
        /*! \brief Add command line options for class variables
         * \param[inout] pargs The command line parameters
         */
        void addOptions(std::vector<t_pargs> *pargs);

        /*! \brief Write histograms for bond length distributions etc.
         * \param[in] oenv GROMACS output environment structure
         */
        void writeHistogram(const gmx_output_env_t *oenv);
            
        /*! \brief Store the bond lengths etc. in the force field
         * \param[in]  fp  A file to print information to
         * \param[out] pd  The force field structure to update
         */
        void updatePoldata(FILE             *fp,
                           Poldata          *pd);
             
        /*! \brief Extract bond lengths, angles etc. from molecules
         * \param[in]  fp     File pointer for information
         * \param[in]  mp     MolProp array
         * \param[out] mymols MyMol array will be filled here
         * \param[in]  pd     Force field structure
         * \param[in]  gms    Selection of compounds
         * \param[in]  method QM method
         * \param[in]  basis  QM basis set
         * \param[in]  strict Whether to strictly adhere to the QM level of theory
         */                          
        void extractGeometries(FILE                       *fp,
                               const std::vector<MolProp> &mp,
                               std::vector<MyMol>         *mymols,
                               const Poldata              &pd,
                               const MolSelect            &gms,
                               const std::string          &method,
                               const std::string          &basis,
                               bool                        strict);

        /*! \brief Write how many bonds etc. were found
         * \param[in] fp File to write to
         */
        void writeSummary(FILE *fp);

    };
    
}

#endif
