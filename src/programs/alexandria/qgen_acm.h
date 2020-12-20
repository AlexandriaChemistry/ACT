/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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
 
 
#ifndef QGEN_ACM_H
#define QGEN_ACM_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state.h"

#include "molprop.h"
#include "poldata.h"

struct t_atoms;

enum class eQgen {
    OK, 
    NOTCONVERGED, 
    NOSUPPORT, 
    ERROR
};

namespace alexandria
{

class QgenAcm
{
    public:
    
        /*! \brief Constructor
         * 
         * \param in pd     Force field information
         * \param in atoms  Atoms data
         * \param in qtotal Total charge for the compound
         */
        QgenAcm(const Poldata           *pd,
                t_atoms                 *atoms,
                int                      qtotal);

        /*! \brief Routine that computes the charges
         * 
         * \param in  fp      File for logging information
         * \param in  molname Molecule name for output
         * \param in  pd      Force field information
         * \param out atoms   Atoms structure, charges will be updated in this
         * \param in  x       Atomic coordinates    
         * \param in  bonds   List of bonds in this compound
         */
        eQgen generateCharges(FILE              *fp,
                              const std::string &molname,
                              const Poldata     *pd,
                              t_atoms           *atoms,
                              const gmx::HostVector<gmx::RVec> x,
                              const std::vector<Bond> &bonds);
                            
        const char *message() const;
        
        int getRow(int atom);

        double getQ(int atom);

        void checkSupport(const Poldata *pd);

        double getZeta(int atom);

        void dump(FILE *fp, t_atoms *atoms);

    private:
        eQgen                            eQGEN_       = eQgen::OK;
        gmx_bool                         bWarned_     = false;
        double                           qtotal_      = 0;
        gmx_bool                         bHaveShell_  = false;
        ChargeType                       ChargeType_  = ChargeType::Point;
        int                              natom_       = 0;
        std::vector<int>                 atomnr_;
        std::vector<double>              chi0_, rhs_, jaa_;
        std::vector<gmx::RVec>           x_;
        // Store ids locally, first for chargedistristribution
        std::vector<Identifier>          qdist_id_;
        // Store ids locally, second for EEM and bond charge corrections
        std::vector<Identifier>          acm_id_;
        //! The atoms/shells to optimize charges for
        std::vector<int>                 nonFixed_;
        //! The atoms/shells not to optimize charges for
        std::vector<int>                 fixed_;
        //! Reverse mapping of the charges
        std::map<int, int>               nfToGromacs_;
        //! Mapping from nonFixed particles to shells
        std::map<int, int>               myShell_;
        std::vector<int>                 row_;       
        std::vector<double>              q_, zeta_, qsave_, zetasave_;
        std::vector<std::vector<double>> Jcc_;

        /*! \brief Re-read the EEM parameters from the FF
         *
         * Update the parameters for the Alexandria Charge model.
         * This includes, chi, JAA, zeta. In case a split charge
         * equilibration algorithm is used also the bond charge
         * correction parameters will be updated.
         * \param[in] pd  Force field database
         */
        void updateParameters(const Poldata *pd);

        double calcJ(rvec   xI, 
                     rvec   xJ,
                     double zetaI,
                     double zetaJ,
                     int    rowI,
                     int    rowJ,
                     double epsilonr);

        /*! \brief Return delta chi and bond hardness for atom pair
         *
         * \param[in]  pd  Poldata structure
         * \param[in]  ai  Atom id i
         * \param[in]  aj  Atom id j
         * \param[in]  bondorder The bond order for this bond
         * \param[out] deltachi the electronegativity correction
         * \param[out] hardness the bond hardness
         */
        void getBccParams(const Poldata *pd,
                          int            ai,
                          int            aj,
                          int            bondorder,
                          double        *deltachi,
                          double        *hardness);
    
        void copyChargesToAtoms(t_atoms *atoms);
        
        /*! \brief Compute the Jcc matrix
         *
         * Only half the matrix is computed for efficiency,
         * but the symmetric values are filled in anyway. That
         * means all the interactions are in there twice. However,
         * the off-diagonal numbers are multiplied by 0.5 such that
         * the total electrostatic potetential is correct nevertheless.
         * The diagonal is filled with the atomic hardness values.
         * \param[in] epsilonr Relative  dielectric constant
         * \param[in] bYang    Whether or not the Yang and Sharp model is used
         * \param[in] bRappe   Whether or not the Rappe and Goddard model is used
         */
        void calcJcc(double epsilonr,
                     bool   bYang,
                     bool   bRappe);
                     
        /*! \brief Compute shell potential at atom position
         * This takes into account all the shells in the molecule.
         * \param[in] top_ndx  Atom number
         * \param[in] epsilonr Relative  dielectric constant
         * \return The potential
         */
        double calcJcs(int      top_ndx,
                       double   epsilonr);

        void solveEEM(FILE *fp);
        
        /*! \brief Perform the split charge equilibration algorithm
         *
         * \param[in] fp    File for logging
         * \param[in] pd    Force field information
         * \param[in] bonds List of bonds in the compound
         */
        void solveSQE(FILE                    *fp,
                      const Poldata           *pd,
                      const std::vector<Bond> &bonds);
        
        void updatePositions(gmx::HostVector<gmx::RVec> x, t_atoms *atoms);

        double calcSij(int i, int j);

        void calcRhs(double epsilonr);
};
}
#endif
