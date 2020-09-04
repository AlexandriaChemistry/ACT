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
 
 
#ifndef POLDATA_EEMPROPS_H
#define POLDATA_EEMPROPS_H

#include "gmxpre.h"

#include <string>
#include <vector>

#include "gromacs/coulombintegrals/coulombintegrals.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "communication.h"
#include "plistwrapper.h"

namespace alexandria
{
/*! \brief
 * Contains information needed for electrostatic interactions.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class RowZetaQ
{
    public:

        RowZetaQ () {}

        /*! \brief
         * RowZetaQ constructor
         *
         * \param[in] row      The row in the periodic table for each of the charge components
         * \param[in] zeta     Inverse screening length of each of the components
         * \param[in] q        Charge of each of the components
         */

        RowZetaQ(int row, double zeta, double q);

        /*! \brief
         * Return the row in the periodic table
         * for each of the charge components
         */
        int row() const { return row_; };

        /*! \brief
         * Set the row in the periodic table
         * for each of the charge components
         *
         * \param[in] row The row in the periodic table for each of the charge components
         */
        void setRow(int row) { row_ = row; }

        /*! \brief
         * Return the charge of each of the components
         */
        double q() const { return q_; }

        /*! \brief
         * Set the charge of each of the components
         *
         * \param[in] q  Charge of each of the components
         */
        void setQ(double q)
        {
            //GMX_RELEASE_ASSERT(!fixedQ_, "Trying to modify a fixed charge");
            q_ = q;
        }

        /*! \brief
         * Return the inverse screening length of each of the components
         */
        double zeta() const { return zeta_; }

        /*! \brief
         * Set the inverse screening length of each of the components
         *
         * \param[in] z  the inverse screening length of each of the components
         */
        void setZeta(double z) { zeta_ = z; }

        /*! \brief
         * Return reference (starting) value for zeta
         */
        double zetaRef() const { return zetaRef_; }

        /*! \brief
         * Set reference (starting) value for zeta
         *
         * \param[in] z  Reference value for zeta
         */
        void setZetaRef(double z) { zetaRef_ = z; }

        /*! \brief
         * Return parameter optimization index
         */
        int zIndex() const { return zindex_; }

        /*! \brief
         * Set parameter optimization index
         *
         * \param[in] zi optimization index
         */
        void setZindex(int zi) { zindex_ = zi; }

        /*! \brief
         * Return true if the charge is fixed
         */
        bool fixedQ() const { return fixedQ_; }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:
        int    row_;
        double zeta_;
        double q_;
        double zetaRef_;
        int    zindex_;
        bool   fixedQ_;
};

using RowZetaQIterator      = typename std::vector<RowZetaQ>::iterator;
using RowZetaQConstIterator = typename std::vector<RowZetaQ>::const_iterator;

class Eemprops
{
    public:

        Eemprops () {}

        Eemprops(const std::string      &name,
                 const std::string      &rowstr,
                 const std::string      &zetastr,
                 const std::string      &zeta_sigma,
                 const std::string      &qstr,
                 double                  qref,
                 double                  J0,
                 double                  J0_sigma,
                 double                  chi0,
                 double                  chi0_sigma);

        int getNzeta() const { return rzq_.size(); }

        const char *getName() const { return name_.c_str(); }

        const char *getZetastr() const { return zetastr_.c_str(); }

        const char *getZeta_sigma() const { return zeta_sigma_.c_str(); }

        const char *getQstr() const { return qstr_.c_str(); }

        const char *getRowstr() const { return rowstr_.c_str(); }
        
        double getQref() const { return qref_; }

        double getJ0() const { return J0_; }

        double getJ0_sigma() const { return J0_sigma_; }

        double getChi0() const { return chi0_; }

        double getChi0_sigma() const { return chi0_sigma_; }

        void setName(const std::string &name) { name_ = name; }

        void setRowZetaQ(const std::string &rowstr,
                         const std::string &zetastr,
                         const std::string &qstr);

        void setZetastr(const std::string &zetastr) {zetastr_ = zetastr; }

        void setZeta_sigma(const std::string &zeta_sigma) {zeta_sigma_ = zeta_sigma; }
        
        void setQref(double qref) { qref_ = qref;} 

        void setJ0(double J0) { J0_ = J0; }

        void setJ0_sigma(double J0_sigma) { J0_sigma_ = J0_sigma; }

        void setChi0(double chi0) { chi0_ = chi0; }

        void setChi0_sigma(double chi0_sigma) { chi0_sigma_ = chi0_sigma; }

        double getZeta(int index) const { return rzq_[index].zeta(); }

        double getQ(int index) const { return rzq_[index].q(); }

        int getRow(int index) const { return rzq_[index].row(); }

        void setZeta(int index, double zeta) { rzq_[index].setZeta(zeta); }

        void setQ(int index, double q) { rzq_[index].setQ(q); }

        void setRow(int index, int row) { rzq_[index].setRow(row); }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:
        std::string             name_;
        std::string             rowstr_;
        std::string             zetastr_;
        std::string             zeta_sigma_;
        std::string             qstr_;
        double                  qref_;
        double                  J0_;
        double                  J0_sigma_;
        double                  chi0_;
        double                  chi0_sigma_;
        std::vector<RowZetaQ>   rzq_;
};

using EempropsIterator      = typename std::vector<Eemprops>::iterator;
using EempropsConstIterator = typename std::vector<Eemprops>::const_iterator;

class BondCorrection
{
 public:
    /*! \brief Empty constructor 
     */
    BondCorrection() {}
    /*! \brief Constructor
     * \param name                    Bond type
     * \param hardness                Bond Hardness
     * \param hardness_sigma          Bond Hardness sigma
     * \param electronegativity       Delta Electronegativity
     * \param electronegativity_sigma Delta Electronegativity sigma
     */
    BondCorrection(const std::string &name,
                   double             hardness,
                   double             hardness_sigma,
                   double             electronegativity,
                   double             electronegativity_sigma) :
        name_(name),
        hardness_(hardness), 
        hardness_sigma_(hardness_sigma), 
        electronegativity_(electronegativity),
        electronegativity_sigma_(electronegativity_sigma)
        {}
    
    //! \brief return the name of the bondtype
    const std::string &name() const { return name_; }

    //! \brief return the hardness
    double hardness() const { return hardness_; }
    
    //! \brief return the hardness sigma
    double hardness_sigma() const { return hardness_sigma_; }
    
    //! \brief return the electronegativity
    double electronegativity() const { return electronegativity_; }
    
    //! \brief return the electronegativity sigma
    double electronegativity_sigma() const { return electronegativity_sigma_; }
    
    CommunicationStatus Send(const t_commrec *cr, int dest);

    CommunicationStatus Receive(const t_commrec *cr, int src);

 private:
    //! \brief Bond type name
    std::string name_;
    //! \brief Bond Hardness
    double      hardness_;
    //! \brief Bond Hardness sigma
    double      hardness_sigma_;
    //! \brief Delta Electronegativity
    double      electronegativity_;
    //! \brief Delta Electronegativity sigma
    double      electronegativity_sigma_;
    
};

using BondCorrectionIterator      = typename std::vector<BondCorrection>::iterator;
using BondCorrectionConstIterator = typename std::vector<BondCorrection>::const_iterator;

} // namespace aleaxndria
#endif
