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

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "act/utility/communicationrecord.h"

struct t_excls;
struct t_symtab;

namespace alexandria
{


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

        CommunicationStatus Send(const CommunicationRecord *cr, int dest);

        CommunicationStatus Receive(const CommunicationRecord *cr, int src);

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

        CommunicationStatus Send(const CommunicationRecord *cr, int dest);

        CommunicationStatus Receive(const CommunicationRecord *cr, int src);

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

        CommunicationStatus Send(const CommunicationRecord *cr, int dest);

        CommunicationStatus Receive(const CommunicationRecord *cr, int src);

    private:
        std::string central_;
        std::string attached_;
        int         numattach_;
};

using SymchargesIterator      = typename std::vector<Symcharges>::iterator;
using SymchargesConstIterator = typename std::vector<Symcharges>::const_iterator;

} // namespace aleaxndria
#endif
