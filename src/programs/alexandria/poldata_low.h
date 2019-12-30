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
 
 
#ifndef POLDATA_LOW_H
#define POLDATA_LOW_H

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

enum VsiteType
{
    evtLINEAR       = 0,
    evtPLANAR       = 1,
    evtRING_PLANAR  = 2,
    evtIN_PLANE     = 3,
    evtOUT_OF_PLANE = 4,
    evtALL          = 5,
    evtNR           = 6
};

/*! \brief
 * Convert interaction type to string
 */
const char *iType2string(InteractionType iType);

/*! \brief
 * Convert string to interaction type
 */
InteractionType string2iType(const char *string);

/*! \brief
 * Convert virtual site type to string
 */
const char *vsiteType2string(VsiteType vType);

/*! \brief
 * Convert string to virtual site type
 */
VsiteType string2vsiteType(const char *string);


/*! \brief
 * Contains all the information realted to
 * alexandria force field atom types.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Ffatype
{
    public:

        /*! \brief
         * Fftype constructor
         *
         * \param[in] desc        Description
         * \param[in] type        Atom type
         * \param[in] ptype       Polarizability type
         * \param[in] btype       Bond type
         * \param[in] ztype       Zeta type
         * \param[in] elem        Element name
         * \param[in] fixed       Fix the Van der Waals parameters?
         * \param[in] vdwparams   Van der Waals parameters
         * \param[in] refEnthalpy Reference Enthalpy of Formation
         */
        Ffatype(const std::string &desc,
                const std::string &type,
                const std::string &ptype,
                const std::string &btype,
                const std::string &ztype,
                const std::string &elem,
                bool               fixed,
                const std::string &vdwparams,
                const std::string &refEnthalpy);

        /*! \brief
         * Fftype default constructor
         */
        Ffatype () {}

        /*! \brief
         * Return the decription of atoms
         */
        const std::string &getDesc() const { return desc_; }

        /*! \brief
         * Return the type of atoms
         */
        const std::string &getType() const { return type_; }

        /*! \brief
         * Return the polarizability type of atoms
         */
        const std::string &getPtype() const { return ptype_; }

        /*! \brief
         * Return the bond type of atoms
         */
        const std::string &getBtype() const { return btype_; }

        /*! \brief
         * Return the zeta type of atoms
         */
        const std::string &getZtype() const { return ztype_; }

        /*! \brief
         * Return the element name of atoms
         */
        const std::string &getElem() const { return elem_; }

        /*! \brief
         * Return the Van der Waals parameters of atoms
         */
        const std::string &getVdwparams() const { return vdwparams_; }

        /*! \brief
         * Return whether the Van der Waals parameters are fixed
         */
        bool fixed() const { return fixed_; }

        /*! \brief
         * Determine whether the Van der Waals parameters are fixed
         */
        void setFixed(bool fixed) { fixed_ = fixed; }

        /*! \brief
         * Set the Van der Waals parameters.
         * \param[in] param The string of parameters
         */ 
        void setVdwparams(const std::string &param) 
        {
            if (!fixed_)
            {
                vdwparams_ = param;
            }
            else
            {
                GMX_THROW(gmx::InvalidInputError("It is not allowed to modify Van der Waals parameters"));
            }
        }

        /*! \brief
         * Set the modified flag
         * \param[in] modified The flag, true or false
         */
        void setModified(bool modified) { modified_ = modified; }
    
        /*! \brief
         * Returns the modified flag
         */
        bool modified() const { return modified_; }

        /*! \brief
         * Return the reference enthalpy of formation of atoms
         */
        const std::string &getRefEnthalpy() const { return refEnthalpy_; }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:
        std::string desc_;
        std::string type_;
        std::string ptype_;
        std::string btype_;
        std::string ztype_;
        std::string elem_;
        bool        fixed_;
        std::string vdwparams_;
        std::string refEnthalpy_;
        bool        modified_;
};

using FfatypeIterator      = typename std::vector<Ffatype>::iterator;
using FfatypeConstIterator = typename std::vector<Ffatype>::const_iterator;


/*! \brief
 * Contains all the information realted to
 * polarizability.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Ptype
{
    public:

        Ptype () {}

        /*! \brief
         * Ptype constructor
         *
         * \param[in] ptype           Polarizability type
         * \param[in] miller          Miller equivalent of the polarizability type
         * \param[in] bosque          Bosque equivalent of the polarizability type
         * \param[in] polarizability  Polarizability value
         * \param[in] sigPol          Uncertainty in the calulated polarizability
         */

        Ptype(const std::string &ptype,
              const std::string &miller,
              const std::string &bosque,
              double             polarizability,
              double             sigPol);

        /*! \brief
         * Return the polarizability type
         */
        const std::string &getType() const { return type_; }

        /*! \brief
         * Return Miller equivalent of the polarizability type
         */
        const std::string &getMiller() const { return miller_; }

        /*! \brief
         * Return Bosque equivalent of the polarizability type
         */
        const std::string &getBosque() const { return bosque_; }

        /*! \brief
         * Set the polarizability value
         *
         * \param[in] polarizability  Polarizability value
         */
        void setPolarizability(double polarizability) { polarizability_ = polarizability; }

        /*! \brief
         * Return the polarizability value
         *
         */
        double getPolarizability() const { return polarizability_; }

        /*! \brief
         * Set the uncertainty in the calulated polarizability
         *
         * \param[in] sigPol  uncertainty
         */
        void setSigPol(double sigPol) { sigPol_ = sigPol; }

        /*! \brief
         * Return the uncertainty in the calulated polarizability
         *
         */
        double getSigPol() const { return sigPol_; }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:
        std::string type_;
        std::string miller_;
        std::string bosque_;
        double      polarizability_;
        double      sigPol_;


};

using PtypeIterator      = typename std::vector<Ptype>::iterator;
using PtypeConstIterator = typename std::vector<Ptype>::const_iterator;

class Vsite
{
    public:

        Vsite () {};

        Vsite(const std::string &atype,
              const std::string &type,
              int                number,
              double             distance,
              double             angle,
              int                ncontrolatoms);

        /*! \brief
         * Return the atom type to which the vsites are connected.
         */
        const std::string &atype() const { return atype_; }

        /*! \brief
         * Return the type pf virtual site.
         */
        const VsiteType &type() const { return type_; }

        /*! \brief
         * Return the distance between the atom and the virtual site.
         */
        double distance() const {return distance_; }

        /*! \brief
         * Return the angle between the atom and the virtual site.
         */
        double angle() const {return angle_; }

        /*! \brief
         * Return the number of virtual sites connected to the atom.
         */
        int nvsite() const {return number_; }

        /*! \brief
         * Return the number of atoms needed to locate the vsite.
         */
        int ncontrolatoms() const {return ncontrolatoms_; }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);


    private:
        std::string atype_;
        VsiteType   type_;
        int         number_;
        double      distance_;
        double      angle_;
        int         ncontrolatoms_;
};

using VsiteIterator      = typename std::vector<Vsite>::iterator;
using VsiteConstIterator = typename std::vector<Vsite>::const_iterator;

/*! \brief
 * Contains the listed force including
 * bonds or angles, or dihedrals
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class ListedForce
{

    public:

        ListedForce () {}

        /*! \brief
         * ListedForce constructor
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] fixed       Set whether the parameters are mutable
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         */
        ListedForce(const std::vector<std::string> &atoms,
                    const std::string              &params,
                    bool                            fixed,
                    double                          refValue,
                    double                          sigma,
                    size_t                          ntrain);
                    
        /*! \brief
         * ListedForce constructor
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] fixed       Set whether the parameters are mutable
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         * \param[in] bondOrder   The order of the covalent bond
         */
        ListedForce(const std::vector<std::string> &atoms,
                    const std::string              &params,
                    bool                            fixed,
                    double                          refValue,
                    double                          sigma,
                    size_t                          ntrain,
                    size_t                          bondOrder);

        /*! \brief
         * Return the vector of atom bond types involved in the listed force
         */
        const std::vector<std::string> &atoms() const { return atoms_; }

        /*! \brief
         * Return a condensed atom string involved in the listed force
         */
        const std::string &condensed_atoms() const { return condensed_atoms_; }
        
        /*! \brief
         * Return a condensed atom string involved in the listed force
         * in reversed order.
         */
        const std::string &reverse_condensed_atoms() const { return reverse_condensed_atoms_; }

        /*! \brief
         * Return whether the bonded parameters are fixed
         */
        bool fixed() const { return fixed_; }

        /*! \brief
         * Determine whether the bonded parameters are fixed
         */
        void setFixed(bool fixed) { fixed_ = fixed; }

        /*! \brief
         * Set the listed force parameters
         */
        void setParams(const std::string &params)
        {
            if (fixed_)
            {
                GMX_THROW(gmx::InvalidInputError("It is not allowed to modify this bonded parameter"));
            }
            else
            {
                params_ = params;
            }
        }

        /*! \brief
         * Return the listed force parameters
         */
        const std::string &params() const { return params_; }

        /*! \brief
         * Set the modified flag
         * \param[in] modified The flag, true or false
         */
        void setModified(bool modified) { modified_ = modified; }
    
        /*! \brief
         * Returns the modified flag
         */
        bool modified() const { return modified_; }

        /*! \brief
         * Set reference value
         */
        void setRefValue(double value) { refValue_ = value; }

        /*! \brief
         * Return reference value
         */
        double refValue() const { return refValue_; }

        /*! \brief
         * Set uncertainty in reference value
         */
        void setSigma(double sigma) { sigma_ = sigma; }

        /*! \brief
         * Return uncertainty in reference value
         */
        double sigma() const { return sigma_; }

        /*! \brief
         * Set the number of molecules in the training set
         */
        void setNtrain(int ntrain) { ntrain_ = ntrain; }

        /*! \brief
         * Return the number of molecules in the training set
         */
        int ntrain() const { return ntrain_; }
        
        /*! \brief
         * Set the bond order for a specific covalent bond
         */
        void setBondOrder(int bondOrder) {bondOrder_ = bondOrder;}
        
        /*! \brief
         * Return the bond order for a specific covalent bond
         */
        size_t bondOrder() const {return bondOrder_;}

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:

        const std::vector<std::string> atoms_;
        std::string                    condensed_atoms_;
        std::string                    reverse_condensed_atoms_;
        std::string                    params_;
        bool                           fixed_;
        double                         refValue_;
        double                         sigma_;
        size_t                         ntrain_;
        size_t                         bondOrder_;
        bool                           modified_;
        // Generate condensed atom strings.
        void MakeCondensed();
};

using ListedForceIterator      = typename std::vector<ListedForce>::iterator;
using ListedForceConstIterator = typename std::vector<ListedForce>::const_iterator;


/*! \brief
 * Contains the listed forces including
 * bonds and angles, and dihedrals
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class ListedForces
{

    public:

        ListedForces () {}

        /*! \brief
         * ListedForces constructor
         *
         * \param[in] iType       Interaction type
         * \param[in] function    The function used to calculate the force corresponding to the iType
         * \param[in] Unit        The unit of the reference value
         */
        ListedForces(const std::string      iType,
                     const std::string     &function,
                     const std::string     &unit);

        /* \brief
         * Return the unit of the reference value
         */
        const std::string &unit() const {return unit_; }

        /* \brief
         * Return the function name
         */
        const std::string &function() const {return function_; }

        /* \brief
         * Return the function type
         */
        const unsigned int &fType() const {return fType_; }

        /* \brief
         * Return the interaction type
         */
        const InteractionType &iType() const {return iType_; }

        /* \brief
         * Find a listed force in the force_ vector
         *
         * \param[in] atoms   A vector of atom bond types invloved in the listed force
         */
        ListedForceIterator findForce(const std::vector<std::string> &atoms);

        /* \brief
         * Find a listed force in the force_ vector, ignore order of atoms
         *
         * \param[in] atoms   A vector of atom bond types invloved in the listed force
         */
        ListedForceConstIterator findForceRandomOrder(const std::vector<std::string> &atoms) const;

        /* \brief
         * Find a listed force in the force_ vector
         *
         * \param[in] atoms   A vector of atom bond types invloved in the listed force
         */
        ListedForceConstIterator findForce(const std::vector<std::string> &atoms) const;
                                           
        /* \brief
         * Find a listed force in the force_ vector
         *
         * \param[in] atoms     A vector of atom bond types invloved in the listed force
         * \param[in] bondOrder bond order of the covalent bond
         */
        ListedForceIterator findForce(const std::vector<std::string> &atoms, 
                                      size_t                          bondOrder);

        /* \brief
         * Find a listed force in the force_ vector
         *
         * \param[in] atoms   A vector of atom bond types invloved in the listed force
         * \param[in] bondOrder bond order of the covalent bond
         */
        ListedForceConstIterator findForce(const std::vector<std::string> &atoms,
                                           size_t                          bondOrder) const;

        /* \brief
         * Set the listed force paramters
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         */
        bool setForceParams(const std::vector<std::string> &atoms,
                            const std::string              &params,
                            double                          refValue,
                            double                          sigma,
                            size_t                          ntrain);
                            
        /* \brief
         * Set the listed force paramters
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         * \param[in] bondOrder bond order of the covalent bond
         */
        bool setForceParams(const std::vector<std::string> &atoms,
                            const std::string              &params,
                            double                          refValue,
                            double                          sigma,
                            size_t                          ntrain,
                            size_t                          bondOrder);

        /* \brief
         * Generate a ListedForce object and
         * add it to the force_ vector, if it does not exist already
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] fixed       Determines whether this parameter is mutable
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         */
        void addForce(const std::vector<std::string> &atoms,
                      const std::string              &params,
                      bool                            fixed,
                      double                          refValue,
                      double                          sigma,
                      size_t                          ntrain);
                      
        /* \brief
         * Generate a ListedForce object and
         * add it to the force_ vector, if it does not exist already
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] fixed       Determines whether this parameter is mutable
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         * \param[in] bondOrder bond order of the covalent bond
         */
        void addForce(const std::vector<std::string> &atoms,
                      const std::string              &params,
                      bool                            fixed,
                      double                          refValue,
                      double                          sigma,
                      size_t                          ntrain,
                      size_t                          bondOrder);

        /* \brief
         * Search a listed force in the force_ vector
         *
         * \param[in]  atoms       A vector of atom bond types invloved in the listed force
         * \param[in]  params      Listed force paramters
         * \param[out] refvalue    The reference value such as the reference bond length
         * \param[out] sigma       Uncertainty in the calculated reference value
         * \param[out] ntrain      Number of molecules in the training set
         * \param[in]  randomAtomOrder If true a random order of atoms is accepted
         */
        bool searchForce(std::vector<std::string> &atoms,
                         std::string              &params,
                         double                   *refValue,
                         double                   *sigma,
                         size_t                   *ntrain,
                         bool                      randomAtomOrder) const;
                         
        /* \brief
         * Search a listed force in the force_ vector
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         * \param[in] bondOrder bond order of the covalent bond
         */
        bool searchForceBondOrder(std::vector<std::string> &atoms,
                                  std::string              &params,
                                  double                   *refValue,
                                  double                   *sigma,
                                  size_t                   *ntrain,
                                  size_t                    bondOrder) const;


        /*! \brief
         * Erase the force_ vector
         */
        void eraseListedForce()
        {
            force_.clear();
        }

        ListedForceIterator forceBegin() {return force_.begin(); }

        ListedForceConstIterator forceBegin() const {return force_.begin(); }

        ListedForceIterator forceEnd() {return force_.end(); }

        ListedForceConstIterator forceEnd() const {return force_.end(); }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:

        InteractionType          iType_;
        const std::string        function_;
        const std::string        unit_;
        std::vector<ListedForce> force_;
        unsigned int             fType_;
};

using ListedForcesIterator      = typename std::vector<ListedForces>::iterator;
using ListedForcesConstIterator = typename std::vector<ListedForces>::const_iterator;


/*! \brief
 * Contains Bosque polarizability.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Bosque
{
    public:

        Bosque () {}

        /*! \brief
         * Bosque constructor
         *
         * \param[in] bosque          Bosque atom type name
         * \param[in] polarizability  Polarizability value
         */
        Bosque(const std::string &bosque, double polarizability);

        /*! \brief
         * Return Bosque equivalent of the polarizability type
         */
        const std::string &getBosque() const { return bosque_; }

        /*! \brief
         * Return polarizability value
         */
        double getPolarizability() const { return polarizability_; }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:
        std::string bosque_;
        double      polarizability_;

};

using BosqueIterator      = typename std::vector<Bosque>::iterator;
using BosqueConstIterator = typename std::vector<Bosque>::const_iterator;


/*! \brief
 * Contains Miller polarizability.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Miller
{
    public:

        Miller () {}

        /*! \brief
         * Miller constructor
         *
         * \param[in] miller            Miller atom type name
         * \param[in] atomnumber        Atomic number
         * \param[in] tauAhc            Polarizability description tau
         * \param[in] alphaAhp          Polarizability description alpha
         * \param[in] alexandria_equiv  Alexandria type
         */
        Miller(const std::string &miller,
               int                atomnumber,
               double             tauAhc,
               double             alphaAhp,
               const std::string &alexandria_equiv);

        /*! \brief
         * Return Miller atom type name
         */
        const std::string &getMiller() const { return miller_; }

        /*! \brief
         * Return atomic number
         */
        int getAtomnumber() const { return atomnumber_; }

        /*! \brief
         * Return polarizability description tau
         */
        double getTauAhc() const { return tauAhc_; }

        /*! \brief
         * Return polarizability description alpha
         */
        double getAlphaAhp() const { return alphaAhp_; }

        /*! \brief
         * Return Alexandria type
         */
        const std::string &getAlexandriaEquiv() const { return alexandria_equiv_; }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:
        std::string miller_;
        int         atomnumber_;
        double      tauAhc_;
        double      alphaAhp_;
        std::string alexandria_equiv_;
};

using MillerIterator      = typename std::vector<Miller>::iterator;
using MillerConstIterator = typename std::vector<Miller>::const_iterator;

class Symcharges
{
    public:

        Symcharges () {}

        Symcharges(const std::string &central,
                   const std::string &attached,
                   int                numattach);

        const std::string &getCentral() const { return central_; }

        const std::string &getAttached() const { return attached_; }

        int getNumattach() const { return numattach_; }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:
        const std::string central_;
        const std::string attached_;
        int               numattach_;
};

using SymchargesIterator      = typename std::vector<Symcharges>::iterator;
using SymchargesConstIterator = typename std::vector<Symcharges>::const_iterator;

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
        double                  J0_;
        double                  J0_sigma_;
        double                  chi0_;
        double                  chi0_sigma_;
        std::vector<RowZetaQ>   rzq_;
};

using EempropsIterator      = typename std::vector<Eemprops>::iterator;
using EempropsConstIterator = typename std::vector<Eemprops>::const_iterator;

} // namespace aleaxndria
#endif
