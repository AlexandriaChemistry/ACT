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

#include <map>
#include <string>
#include <vector>

#include "gromacs/coulombintegrals/coulombintegrals.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "communication.h"
#include "identifier.h"
#include "mutability.h"
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
         * \param[in] acmtype     Alexandria Charge Model type
         * \param[in] elem        Element name
         * \param[in] mass        Mass of the particle
         * \param[in] atomnumber  Atomic number
         * \param[in] charge      The charge on this particle
         * \param[in] row         The row in the periodic table
         * \param[in] mutability  Whether or not this charge can be changed
         * \param[in] refEnthalpy Reference Enthalpy of Formation
         */
        Ffatype(const std::string &desc,
                const std::string &type,
                const std::string &ptype,
                const std::string &btype,
                const std::string &ztype,
                const std::string &acmtype,
                const std::string &elem,
                double             mass,
                int                atomnumber,
                double             charge,
                int                row,
                Mutability         mutability,
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
         * Return an identifier for an interaction type
         * \throws if iType is not present
         */
        Identifier id(InteractionType iType) const;

        /*! \brief Return mass
         * \return mass of the particle
         */
        double mass() const { return mass_; }
        
        /*! \brief Return atomnumber
         * \return atomnumber of the particle
         */
        int atomnumber() const { return atomnumber_; }
    
        /*! \brief Return charge
         * \return charge of the particle
         */
        double charge() const { return charge_; }
        
        /*! \brief Return row in periodic table
         * \return row
         */
        int row() const { return row_; }
    
        /*! \brief
         * Return whether an identifier for an interaction type is present
         */
        bool hasId(InteractionType iType) const
        { 
            return subType_.find(iType) != subType_.end();
        }

        /*! \brief
         * Return whether or not this charge can be changed
         */
        Mutability mutability() const { return mutability_; }
        /*! \brief
         * Return the element name of atoms
         */
        const std::string &getElem() const { return elem_; }

        /*! \brief
         * Return the reference enthalpy of formation of atoms
         */
        const std::string &getRefEnthalpy() const { return refEnthalpy_; }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:
        std::string desc_;
        std::string type_;
        std::string elem_;
        std::string refEnthalpy_;
        double      mass_;
        int         atomnumber_;
        double      charge_;
        int         row_;
        Mutability  mutability_;
        std::map<InteractionType, const std::string> subType_;
};

using FfatypeIterator      = typename std::vector<Ffatype>::iterator;
using FfatypeConstIterator = typename std::vector<Ffatype>::const_iterator;

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

} // namespace aleaxndria
#endif
