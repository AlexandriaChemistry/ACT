/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#ifndef ACT_QGEN_RESP_H
#define ACT_QGEN_RESP_H

#include <cstdio>

#include <vector>

#include "act/basics/chargemodel.h"
#include "gromacs/math/vectypes.h"

struct gmx_output_env_t;
struct t_symtab;

namespace alexandria
{

class ActAtom;
class ForceField;
class ForceFieldParameter;
class MsgHandler;

/*! \brief Class to store one grid point and it's potential.
 * The structure store both the reference and the calculated
 * electrostatic potential.
 * \todo Make clear what unit is used.
 */
class EspPoint
{
public:
    /*! \brief Constructor
     * \param[in] esp Coordinates of the point
     * \param[in] v   The potential in this point
     */
    EspPoint(const gmx::RVec &esp, double v) : esp_(esp), v_(v) {}
        
    //! Return coordinates of the grid point
    const gmx::RVec &esp() const { return esp_; }

    //! Return the reference potential at this grid point
    double v() const { return v_; }

    /*! \brief Set the reference potential at this grid point
     * \param[in] v The new potential
     */
    void setV(double v) { v_ = v; }
    
    //! Return the calculated potential at this grid point
    double vCalc() const { return vCalc_; }
    
    /*! \brief Set the calculated potential at this grid point
     * \param[in] vcalc The new potential
     */
    void setVCalc(double vcalc) { vCalc_ = vcalc; }

    //! Return the electron density in this grid point
    double rho() const { return rho_; }
    
    /*! \brief Set the electron density in this point
     * \param[in] rho The new electron density
     */
    void setRho(double rho) { rho_ = rho; }
    
private:
    //! The coordinates of a point
    gmx::RVec esp_   = { 0, 0,  0};
    //! The measured potential
    double    v_     = 0;
    //! The calculated potential
    double    vCalc_ = 0;
    //! The electron density in the point
    double    rho_   = 0;
};

/*! \brief Class to fit and evaluate charges against an electrostatic potential
 * Holds the atoms (charges, distribution widths, coordinates) and a set of
 * grid points with reference and calculated electrostatic potential. It
 * can fit charges to the potential (RESP-like, with total charge and
 * symmetry constraints), compute statistics and write cube files.
 */
class QgenResp
{
    public:
        //! Default constructor
        QgenResp() = default;

        //! \return the charge distribution type
        ChargeDistributionType chargeType() const { return ChargeDistributionType_; }

        /*! \brief Set option for ESP charge generation
         *
         * \param[in] qd Charge distribution type
         */
        void setChargeDistributionType(ChargeDistributionType qd)
        { ChargeDistributionType_ = qd; }

        //! \return the total charge of the compound (e)
        real getMolecularCharge() const { return qtot_; }

        //! \return the number of ESP grid points
        size_t nEsp() const { return ep_.size(); }

        //! \return all the ESP grid points
        const std::vector<EspPoint> &espPoints() const {return ep_; }

        /*! \brief Return one ESP grid point
         * \param[in] index The index of the point
         * \return the point
         */
        const EspPoint &espPoint(size_t index) const {return ep_[index]; }

        /*! \brief Print the number of ESP points and atoms and the symmetric atoms
         * \param[in] msg_handler For writing the information, nothing is done if nullptr
         */
        void summary(MsgHandler *msg_handler);
        
        /*! \brief Set the inforamtion about atoms
         * The size of arrays in atoms and x is checked and compared 
         * to what was there previously if anything.
         * \param[in] msg_handler For writing messages and warnings
         * \param[in] atoms  The ACT atoms structure
         * \param[in] pd     The force field
         * \param[in] qtotal Total charge of the compound, needed when
         *                   generating charges
         */
        void setAtomInfo(MsgHandler                   *msg_handler,
                         const std::vector<ActAtom>   &atoms,
                         const ForceField             *pd,
                         const int                     qtotal);

        //! \return the number of atoms
        size_t natoms() const { return nAtom_; }
        
        /*! \brief Copy input coordinates to internal ones
         * \param[in] x Input coordinates
         */
        void updateAtomCoords(const std::vector<gmx::RVec> &x);

        //! \return internal coordinates
        const std::vector<gmx::RVec> &coords() const { return x_; }

        /*! \brief Update the charges
         * \param[in] msghandler For error messages
         * \param[in] atoms Vector containing atoms with new charges
         */
        void updateAtomCharges(MsgHandler                 *msghandler,
                               const std::vector<ActAtom> &atoms);

        /*! \brief Set the charge symmetry of the atoms
         * Also counts the number of charges to fit and the number of unique ones.
         * Must be called after setAtomInfo.
         * \param[in] symmetricAtoms For each atom the index of the atom whose
         *                           charge it must equal, or empty for no symmetry
         */
        void setAtomSymmetry(const std::vector<int> &symmetricAtoms);

        /*! \brief Add a grid point with a reference potential
         * \param[in] x Coordinate X (nm)
         * \param[in] y Coordinate Y (nm)
         * \param[in] z Coordinate Z (nm)
         * \param[in] V Reference electrostatic potential (kJ/mol e)
         */
        void addEspPoint(double x,
                         double y,
                         double z,
                         double V);

        /*! \brief Make a grid for potential calculations
         *
         * Generate a grid, all numbers should be in GROMACS units,
         * that is nm.
         * \param[in] spacing Spacing between grid planes
         * \param[in] border  Distance between atoms and gird edge
         * \param[in] x       Atomic coordinates
         */
        void makeGrid(real                          spacing, 
                      real                          border,
                      const std::vector<gmx::RVec> &x);

        /*! \brief Copy the grid definition and the ESP points from another object
         * \param[in] src The object to copy from
         */
        void copyGrid(const QgenResp &src);

        /*! \brief Compute and return statistics of calculated versus reference potential
         * \param[in]  msg_handler For debug output
         * \param[out] rrms        Relative RMS deviation
         * \param[out] cosangle    Cosine of the angle between the reference and
         *                         calculated potential vectors
         * \param[out] mae         Mean absolute error (kJ/mol e)
         * \param[out] mse         Mean signed error, reference minus calculated (kJ/mol e)
         * \return the RMS deviation (Hartree/e)
         */
        double getStatistics(MsgHandler *msg_handler,
                             double     *rrms,
                             double     *cosangle,
                             double     *mae,
                             double     *mse);

        //! \brief Compute the charge density of the atoms and shells at the grid points
        void calcRho();

        /*! \brief Compute the electrostatic potential based on charges.
         * \param[in] msg_handler For warnings etc.
         * \param[in] epsilonr    Relative dielectric constant
         */
        void calcPot(MsgHandler *msg_handler,
                     double      epsilonr);

        /*! \brief Read the grid and the potential from a cube file
         * \param[in] fn       The file name, with as many atoms as this
         * \param[in] bESPonly Only read the coordinates of the atoms, not
         *                     their atomic numbers and charges
         */
        void readCube(const std::string &fn,
                      bool               bESPonly);

        /*! \brief Write a cube file
         * \param[in] fn    The file name to write to
         * \param[in] title Title to add in the file
         * \param[in] oenv  Gromacs output environment
         */
        void writeCube(const std::string      &fn,
                       const std::string      &title,
                       const gmx_output_env_t *oenv);

        /*! \brief Write a cube file containing electron density
         * \param[in] fn    The file name to write to
         * \param[in] title Title to add in the file
         * \param[in] oenv  Gromacs output environment
         */
        void writeRho(const std::string      &fn,
                      const std::string      &title,
                      const gmx_output_env_t *oenv);

        /*! \brief Write a cube file of the calculated potential minus a reference
         * Optionally also writes the deviation as a function of distance to
         * the nearest atom and, next to it, a file correlating reference and
         * calculated potential.
         * \param[in] src    Object with the reference potential, the reference is
         *                   not subtracted if it has no ESP points
         * \param[in] cubeFn The cube file name, nothing is written if empty
         * \param[in] histFn The xvg file name for the deviation as a function
         *                   of distance, nothing is written if empty
         * \param[in] title  Title to add in the files
         * \param[in] oenv   Gromacs output environment
         * \param[in] rho    If non-zero, write the charge density instead of
         *                   the potential
         */
        void writeDiffCube(const QgenResp             *src,
                           const std::string          &cubeFn,
                           const std::string          &histFn,
                           const std::string          &title,
                           const gmx_output_env_t     *oenv,
                           int                         rho);

        /*! \brief Write a histogram of the calculated potential to an xvg file
         * \param[in] fn    The file name, nothing is written if empty
         * \param[in] title Title to add in the file
         * \param[in] oenv  Gromacs output environment
         */
        void writeHisto(const std::string      &fn,
                        const std::string      &title,
                        const gmx_output_env_t *oenv);

        /*!  brief Do the ESP optimization
         *
         * Optimizes the charges using matrix inversion. No restraints are
         * taken into account, except total charge and charge symmetries.
         * \param[in] msg_handler For debugging and info
         * \param[in] epsilonr    Relative dielectric constant.
         */
        void optimizeCharges(MsgHandler *msg_handler,
                             double      epsilonr);

        /*! \brief Write reference versus calculated potential to an xvg file
         * Values are in Hartree/e, the first natoms points are a separate data set.
         * \param[in] potcomp The file name, nothing is written if nullptr
         * \param[in] oenv    Gromacs output environment
         */
        void potcomp(const char             *potcomp,
                     const gmx_output_env_t *oenv);

        /*! \brief Generate a pdb file containing atoms and grid points
         * The grid points are colored according to the ESP (b-factor field).
         * Three copies of the grid are given, the QM grid, the ACT grid
         * and the difference. The QM grid also has the atoms.
         * \param[in] atoms      The atom name information
         * \param[in] pdbdiff    The filename
         */
        void writePdbComparison(const std::vector<ActAtom>   &atoms,
                                const std::string            &pdbdiff);

        /*! \brief Update the internal copies of zeta
         * \param[in] atoms The atom info
         * \param[in] pd    The force field
         */
        void updateZeta(const std::vector<ActAtom> &atoms,
                        const ForceField              *pd);
                        
        //! Return the charge for one particle
        double getCharge(int atom) const { return q_[atom]; }


    private:
        /*! \brief Compute RMS, relative RMS, MAE, MSE and cosine of the potentials
         * \param[in] msg_handler For debug output
         */
        void calcStatistics(MsgHandler *msg_handler);

        /*! \brief Make sure the total charge is correct and that symmetry is obeyed
         * \param[in] msg_handler For debugging and info
         */
        void regularizeCharges(MsgHandler *msg_handler);

        /*! \brief Set the charge of one atom
         * \param[in] atom The atom index
         * \param[in] q    The charge (e)
         */
        void setCharge(int atom, double q) { q_[atom] = q; }

        /*! \brief Set the charge distribution width of one atom
         * \param[in] atom The atom index
         * \param[in] zeta The distribution width (1/nm)
         */
        void setZeta(int atom, double zeta) { zeta_[atom] = zeta; }

        //! The charge distribution type
        ChargeDistributionType    ChargeDistributionType_  = ChargeDistributionType::Point;
        //! Total charge of the compound (e)
        int                       qtot_        = 0;
        //! Sum of the fixed charges, i.e. non-fitted atom charges and second charges (e)
        double                    qshell_      = 0;
        //! RMS deviation of calculated from reference potential (Hartree/e)
        double                    rms_         = 0;
        //! Relative RMS deviation of the potential
        double                    rrms_        = 0;
        //! Mean signed error, reference minus calculated potential (kJ/mol e)
        double                    mse_         = 0;
        //! Mean absolute error of the potential (kJ/mol e)
        double                    mae_         = 0;
        //! Cosine of the angle between reference and calculated potential vectors
        double                    cosangle_    = 0;
        //! Origin of the grid (nm)
        dvec                      origin_      = { 0, 0, 0 };
        //! Grid spacing in each dimension (nm)
        dvec                      space_       = { 0, 0, 0 };
        //! Number of grid points in each dimension
        ivec                      nxyz_        = { 0, 0, 0 };
        //! Number of unique fitted charges, taking symmetry into account
        int                       uniqueQ_     = 0;
        //! Number of charges to fit
        size_t                    fitQ_        = 0;
        //! Number of atoms
        size_t                    nAtom_       = 0;
        //! Number of atoms with a fixed charge
        size_t                    nFixed_      = 0;

        //! The charge of each atom (e)
        std::vector<double>      q_;
        //! The force field charge parameter of each atom
        std::vector<const ForceFieldParameter *> charge_;
        //! The charge distribution width of each atom (1/nm)
        std::vector<double>      zeta_;
        //! The atomic number of each atom
        std::vector<int>         atomnumber_;
        //! The row in the periodic table of each atom, for Slater distributions
        std::vector<int>         row_;
        //! Whether the charge of each atom is fitted
        std::vector<bool>        mutable_;
        //! The second charge of each atom, at its position, see ActAtom::chargeVS
        std::vector<double>      qVS_;
        //! The width of the distribution of the second charges
        std::vector<double>      zetaVS_;
        //! The row of the second charges, for Slater distributions
        std::vector<int>         rowVS_;
        //! The atomic coordinates (nm)
        std::vector<gmx::RVec>   x_;
        //! The ESP grid points
        std::vector<EspPoint>    ep_;
        //! For each atom the index of the atom whose charge it must equal
        std::vector<size_t>      symmetricAtoms_;
};

} // namespace

#endif
