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
#ifndef POLDATA_LOW_H
#define POLDATA_LOW_H

#include <map>
#include <string>
#include <vector>

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "communication.h"
#include "identifier.h"
#include "mutability.h"
#include "plistwrapper.h"

namespace alexandria
{

//! Class to distguish virtual site types
enum class VsiteType
{
    //! Linear virtual site (gv_linear) dependent on two atoms
    LINEAR,
    //! Planar virtual site (gv_planar) dependent on three atoms
    PLANAR,
    //! Planar virtual site (gv_ringplanar) dependent on three atoms in a ring
    RING_PLANAR,
    //! Planar virtual site (gv_inplane) dependent on three atoms
    IN_PLANE,
    //! Out-of-plane virtual site (gv_outplane) dependent on three atoms
    OUT_OF_PLANE,
    //! All possible virtual sites
    ALL         
};

/*! \brief
 * Convert virtual site type to string
 * \param[in] vType The virtual site type
 * \return VsiteType name
 */
const char *vsiteType2string(VsiteType vType);

/*! \brief
 * Convert string to virtual site type
 * \param[in] string The string to convert
 * \returns The vsite type
 * \throws if the string is not recognized
 */
VsiteType string2vsiteType(const char *string);

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
