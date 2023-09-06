/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include "act/alexandria/topology.h"
#include "act/basics/chargemodel.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/statistics/statistics.h"

struct gmx_output_env_t;
struct t_symtab;

namespace alexandria
{

class ForceField;

/*! \brief Class to store one grid point and it's potential.
 * The structure store both the reference and the calculated
 * electrostatic potential.
 * TODO Make clear what unit is used.
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
     * \param[in] vCalc The new potential
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

class QgenResp
{
    public:
        QgenResp() {}

        ChargeType chargeType() const { return ChargeType_; }

        /*! \brief Set option for ESP charge generation
         *
         * \param[in] qd Charge distribution type
         */
        void setChargeType(ChargeType qd)
        { ChargeType_ = qd; }

        real getMolecularCharge() const { return qtot_; }

        size_t nEsp() const { return ep_.size(); }

        const std::vector<EspPoint> &espPoints() const {return ep_; }

        const EspPoint &espPoint(size_t index) const {return ep_[index]; }

        void summary(FILE *gp);
        
        /*! \brief Set the inforamtion about atoms
         * The size of arrays in atoms and x is checked and compared 
         * to what was there previously if anything. 
         * \param[in] atoms  The ACT atoms structure
         * \param[in] pd     The force field
         * \param[in] qtotal Total charge of the compound, needed when
         *                   generating charges
         */
        void setAtomInfo(const std::vector<ActAtom>   &atoms,
                         const ForceField             *pd,
                         const int                     qtotal);

        int natoms() const { return nAtom_; }
        
        /*! \brief Copy input coordinates to internal ones
         * \param[in] x Input coordinates
         */
        void updateAtomCoords(const std::vector<gmx::RVec> &x);

        //! \return internal coordinates
        const std::vector<gmx::RVec> &coords() const { return x_; }

        /*! \brief Update the charges
         * \param[in] q Vector containing new charges
         */
        void updateAtomCharges(const std::vector<ActAtom> &atoms);

        /*! \brief Update the charges
         * \param[in] q Vector containing new charges
         */
        void updateAtomCharges(const std::vector<double> &q);

        const std::string &getStoichiometry() const { return stoichiometry_; }

        void setAtomSymmetry(const std::vector<int> &symmetricAtoms);

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

        void copyGrid(QgenResp &src);

        void calcStatistics();

        real getStatistics(real *rrms, real *cosangle, real *mse, real *mae);

        void plotLsq(const gmx_output_env_t *oenv,
                     const char             *ESPcorr);

        void calcRho();

        /*! \brief Compute the electrostatic potential based on charges.
         *
         * \param[in] epsilonr  Relative dielectric constant
         */
        void calcPot(double epsilonr);

        void calcVShell();

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

        void writeDiffCube(QgenResp                   *src,
                           const std::string          &cubeFn,
                           const std::string          &histFn,
                           const std::string          &title,
                           const gmx_output_env_t     *oenv,
                           int                         rho);

        void writeHisto(const std::string      &fn,
                        const std::string      &title,
                        const gmx_output_env_t *oenv);

        /*!  brief Do the ESP optimization
         *
         * Optimizes the charges using matrix inversion. No restraints are
         * taken into account, except total charge and charge symmetries.
         * \param[in] epsilonr Relative dielectric constant.
         */
        void optimizeCharges(double epsilonr);

        // Make sure the total charge is correct and that symmetry is obeyed
        void regularizeCharges();

        void potcomp(const char                 *potcomp,
                     const std::vector<ActAtom> &atoms,
                     const rvec                 *x,
                     const char                 *pdbdiff,
                     const gmx_output_env_t     *oenv);

        real myWeight(int iatom) const;
    
        /*! \brief Update the internal copies of zeta
         * \param[in] atoms The atom info
         * \param[in] pd    The force field
         */
        void updateZeta(const std::vector<ActAtom> &atoms,
                        const ForceField              *pd);
                        
        //! Return the charge for one particle
        double getCharge(int atom) const { return q_[atom]; }

        double getZeta(int atom) const { return zeta_[atom]; }


    private:
        void setCharge(int atom, double q) { q_[atom] = q; }

        void setZeta(int atom, double zeta) { zeta_[atom] = zeta; }

        ChargeType                ChargeType_  = ChargeType::Point;
        int                       qtot_        = 0;
        double                    qshell_      = 0;
        double                    rms_         = 0;
        double                    rrms_        = 0;
        double                    mse_         = 0;
        double                    mae_         = 0;
        double                    cosangle_    = 0;
        dvec                      origin_      = { 0, 0, 0 };
        dvec                      space_       = { 0, 0, 0 };
        ivec                      nxyz_        = { 0, 0, 0 };
        int                       uniqueQ_     = 0;
        int                       fitQ_        = 0;
        int                       nAtom_       = 0;
        int                       nParticle_   = 0;
        int                       nFixed_      = 0;

        //! Total number of parameters
        std::vector<double>      q_;
        std::vector<double>      zeta_;
        std::vector<int>         atomnumber_;
        std::vector<int>         row_;
        std::vector<bool>        mutable_;
        std::vector<gmx::RVec>   x_;
        std::vector<std::string> dzatoms_;
        std::string              stoichiometry_;
        std::vector<EspPoint>    ep_;
        std::vector<int>         symmetricAtoms_;
};

} // namespace

#endif
