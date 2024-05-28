/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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
 */
 
#ifndef QGEN_ACM_H
#define QGEN_ACM_H

#include <cstdio>

#include <vector>

#include "act/molprop/molprop.h"
#include "act/forcefield/forcefield.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state.h"

/*! \brief Charge generation status
 */
enum class eQgen {
    //! Charge generation status
    OK,
    //! Did not converge
    NOTCONVERGED,
    //! No support for all atom types
    NOSUPPORT,
    //! Problem solving the matrix equation
    MATRIXSOLVER,
    //! Incorrect molecule
    INCORRECTMOL,
    //! Unknown error
    ERROR
};

namespace alexandria
{

/*! \brief Class governing charge generation
 *
 * The QgenACM supports two algorithms:
 * Electronegativity equalization method (EEM)
 * Split-charge equilibration (SQE)
 * Which one of the two is used is decided by the
 * presence of interaction of 
 * InteractionType::BONDCORRECTIONS in the ForceField
 * file and structure.
 */
class QgenAcm
{
public:
    
    /*! \brief Constructor
     * 
     * \param[in] pd     Force field information
     * \param[in] atoms  Atoms data
     * \param[in] qtotal The total charge in this compound
     */
    QgenAcm(const ForceField              *pd,
            const std::vector<ActAtom> &atoms,
            int                         qtotal);
    
    /*! \brief Routine that computes the charges
     * 
     * \param[in]  fp      File for logging information
     * \param[in]  molname Molecule name for output
     * \param[in]  pd      Force field information
     * \param[out] atoms   Atoms structure, charges will be updated in this
     * \param[in]  x       Atomic coordinates    
     * \param[in]  bonds   List of bonds in this compound
     */
    eQgen generateCharges(FILE                         *fp,
                          const std::string            &molname,
                          const ForceField             *pd,
                          std::vector<ActAtom>         *atoms,
                          const std::vector<gmx::RVec> &x,
                          const std::vector<Bond>      &bonds);
                          
    /*! \brief Return a status message
     * \return A string
     */
    const char *status() const;
    /*! \brief Set the total charge
     * \param[in] qtotal The new total charge
     */
    void setQtotal(int qtotal) { qtotal_ = qtotal; }
    
    /*! \brief Return the charge corresponding to the atom
     * \param[in] atom index of the atom in the compound
     * \return the charge
     */
    double getQ(int atom);
    
    /*! \brief Return the charge distribution width of the atom
     * \param[in] atom index of the atom in the compound
     * \return the charge distribution width
     */
    double getZeta(int atom);
    
    /*! \brief Return the row corresponding to the atom
     * \param[in] atom index of the atom in the compound
     * \return the corresponding row in the periodic table
     */
    int getRow(int atom);
    
    
private:
    /*! \brief Write debug information
     * \param[in] fp The file pointer to write to
     * \param[in] atoms Information about the atoms
     */
    void dump(FILE                       *fp, 
              const std::vector<ActAtom> *atoms) const;
    
    /*! \brief Check whether the compound has support in the force field
     * \param[in] pd The force field data structure
     */
    void checkSupport(const ForceField *pd);
    
    //! Status of the algorithm
    eQgen                            eQGEN_       = eQgen::OK;
    //! Total charge
    int                              qtotal_      = 0;
    //! Whether or not shells are present
    gmx_bool                         bHaveShell_  = false;
    //! The charge type used in the force field
    ChargeType                       ChargeType_  = ChargeType::Point;
    //! Number of particles in the compound
    int                              natom_       = 0;
    //! Epsilon R
    double                           epsilonr_    = 1;
    //! Atomic number for each of the atoms
    std::vector<int>                 atomicNumber_;
    //! Parameters for EEM algorithm, chi is electronegativity
    std::vector<double>              chi0_;
    //! Parameters for EEM algorithm, eta is atomic hardness
    std::vector<double>              eta_;
    //! Right/hand side in the matrix equation
    std::vector<double>              rhs_;
    //! Atomic coordinates
    std::vector<gmx::RVec>           x_;
    // Store ids locally, first for chargedistristribution
    std::vector<Identifier>          qdist_id_;
    // Store ids locally, second for EEM and bond charge corrections
    std::vector<Identifier>          acm_id_;
    //! The atoms/shells to optimize charges for
    std::vector<int>                 nonFixed_;
    //! The atoms/shells not to optimize charges for
    std::vector<int>                 fixed_;
    //! The row number for each of the atoms
    std::vector<int>                 row_;
    //! The atomic charges
    std::vector<double>              q_;
    //! The distribution widths
    std::vector<double>              zeta_;
    //! The Coulomb matrix
    std::vector<std::vector<double>> Jcc_;
    
    /*! \brief Re-read the EEM parameters from the FF
     *
     * Update the parameters for the Alexandria Charge model.
     * This includes, chi, eta, zeta. In case a split charge
     * equilibration algorithm is used also the bond charge
     * correction parameters will be updated.
     * \param[in] pd    Force field database
     * \param[in] atoms Atoms data structure
     */
    void updateParameters(const ForceField              *pd,
                          const std::vector<ActAtom> &atoms);
    
    /*! \brief Compute Coulomb interaction
     * \param[in] xI       Coordinates for atom I
     * \param[in] xJ       Coordinates for atom J
     * \param[in] zetaI    zeta for atom I
     * \param[in] zetaJ    zeta for atom J
     * \param[in] rowI     row for atom I
     * \param[in] rowJ     row for atom J
     * \param[in] epsilonr Relative dielectric constant
     * \return the coulomb interaction
     */
    double calcJ(rvec   xI, 
                 rvec   xJ,
                 double zetaI,
                 double zetaJ,
                 int    rowI,
                 int    rowJ,
                 double epsilonr);
    
    /*! \brief Return delta_chi and delta_eta for atom pair
     *
     * \param[in]  pd  ForceField structure
     * \param[in]  ai  Atom id i
     * \param[in]  aj  Atom id j
     * \param[in]  bondorder The bond order for this bond
     * \param[out] delta_chi the electronegativity correction
     * \param[out] delta_eta the bond hardness
     */
    void getBccParams(const ForceField *pd,
                      int               ai,
                      int               aj,
                      double            bondorder,
                      double           *delta_chi,
                      double           *delta_eta);
    
    /*! \brief Store the atoms in their destination structure
     * \param[inout] atoms The array with atom properties
     */
    void copyChargesToAtoms(std::vector<ActAtom> *atoms);

    
    /*! \brief Compute the Jcc matrix
     *
     * Only half the matrix is computed for efficiency,
     * but the symmetric values are filled in anyway. That
     * means all the interactions are in there twice. However,
     * the off-diagonal numbers are multiplied by 0.5 such that
     * the total electrostatic potetential is correct nevertheless.
     * The diagonal is filled with the atomic hardness values.
     * \param[in] epsilonr Relative  dielectric constant
     */
    void calcJcc(double epsilonr);
    
    /*! \brief Compute shell potential at atom position
     * This takes into account all the shells in the molecule.
     * \param[in] top_ndx  Atom number
     * \param[in] atoms    Atom information
     * \param[in] epsilonr Relative  dielectric constant
     * \return The potential
     */
    double calcJcs(int                         top_ndx,
                   const std::vector<ActAtom> &atoms,
                   double                      epsilonr);
    
    /*! \brief Solve the matrix equation to determine charges
     * \param[in] fp File pointer for optional outptu
     * If something goes wrong, eQGEN_ will be set.
     */
    void solveEEM(FILE *fp);
    
    /*! \brief Perform the split charge equilibration algorithm
     *
     * \param[in] fp    File for logging
     * \param[in] pd    Force field information
     * \param[in] bonds List of bonds in the compound
     * If something goes wrong, eQGEN_ will be set.
     */
    void solveSQE(FILE                    *fp,
                  const ForceField        *pd,
                  const std::vector<Bond> &bonds);
    
    /*! \brief update the positions
     * \param[in] x The new coordinates
     */        
    void updatePositions(const std::vector<gmx::RVec> &x);
    
    /*! \brief Compute the right/hand side of the matrix equation
     * \param[in] atoms    Atom information
     * \param[in] epsilonr The relative dielectric constant
     */
    void calcRhs(const std::vector<ActAtom> &atoms,
                 double                      epsilonr);
};
}
#endif
