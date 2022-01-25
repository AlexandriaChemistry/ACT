/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#ifndef COMPOSITION_H
#define COMPOSITION_H

#include <string.h>

#include <map>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

#include "communication.h"
#include "experiment.h"

namespace alexandria
{

/*! \brief
 * Specifies the name of an atom type and the number in a molecular composition
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class AtomNum
{
    private:
        /*! \brief
         * Atom name
         */
        std::string catom_;
        /*! \brief
         * Atom number
         */
        int         cnumber_;
    public:
        //! Default constructor
        AtomNum() {};

        /*! \brief
         * Creates a new AtomNum object.
         *
         * \param[in] catom   Atom name
         * \param[in] cnumber Number of copies of this atom
         */
        AtomNum(const char *catom, int cnumber) { SetAtom(catom); SetNumber(cnumber); }

        /*! \brief
         * Creates a new AtomNum object.
         *
         * \param[in] catom   Atom name
         * \param[in] cnumber Number of copies of this atom
         */
        AtomNum(const std::string &catom, int cnumber)
        {
            SetAtom(catom);
            SetNumber(cnumber);
        }

        /*! \brief
         * Return the name of the atom for this AtomNum
         */
        const std::string &getAtom() const { return catom_; }

        /*! \brief
         * Set the name of the atom for this AtomNum
         */
        void SetAtom(const std::string &catom) { catom_ = catom; }

        /*! \brief
         * Set the name of the atom for this AtomNum
         */
        void SetAtom(const char *catom) { catom_.assign(catom); }

        /*! \brief
         * Return the number of atoms for this AtomNum
         */
        int getNumber() const { return cnumber_; }

        /*! \brief
         * Set the number of atoms for this AtomNum
         */
        void SetNumber(int cnumber) { cnumber_ = cnumber; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] cr   GROMACS data structure for MPI communication
         * \param[in] dest Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(const CommunicationRecord *cr,
                                 int                        dest) const;

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] cr  GROMACS data structure for MPI communication
         * \param[in] src Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(const CommunicationRecord *cr,
                                    int                        src);
};

//! Iterator over a vector of AtomNum
using  AtomNumIterator      = typename std::vector<AtomNum>::iterator;
//! Const iterator over a vector of AtomNum
using  AtomNumConstIterator = typename std::vector<AtomNum>::const_iterator;

/*! \brief
 * Contains the molecular composition in terms of atoms
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularComposition
{
    private:
        /*! \brief
         * Composition name
         */
        std::string          compname_;
        /*! \brief
         * A vector of AtomNum object
         */
        std::vector<AtomNum> atomnum_;
    public:
        //! Defult constructor
        MolecularComposition() {}

        /*! \brief
         * Creates a new MolecularComposition object.
         *
         * \param[in] compname  Name of the composition type
         */
        MolecularComposition(const char *compname) { compname_.assign(compname); }

        /*! \brief
         * Creates a new MolecularComposition object.
         *
         * \param[in] compname  Name of the composition type
         */
        MolecularComposition(const std::string &compname)
        {
            compname_ = compname;
        }

        /*! \brief
         * Return the composition name
         */
        const std::string &getCompName() const { return compname_; }

        /*! \brief
         * Set the composition name
         */
        void SetCompName(const std::string &compname)
        {
            compname_ = compname;
        }

        /*! \brief
         * Set the composition name
         */
        void SetCompName(char *compname) { compname_.assign(compname); }

        /*! \brief
         * Add an AtomNum object to the composition
         *
         * \param[in] an  Atom number
         */
        void AddAtom(AtomNum an);

        /*! \brief
         * Remove the atom with name catom from the composition
         *
         * \param[in] catom Atom name
         */
        void DeleteAtom(const std::string &catom);

        /*! \brief
         * Replace the oldatom by newatom
         *
         * \param[in] oldatom   Name of the old atom
         * \param[in] newatom   Name of the new atom
         */
        void ReplaceAtom(const std::string &oldatom,
                         const std::string &newatom);

        //! Return const vector of AtomNum
        const std::vector<AtomNum> atomNumConst() const { return atomnum_; }
    
        /*! \brief
         * Return const iterator pointing to a specific atom or EndAtomNum if not found
         *
         * \param[in] an Atom number
         */
        AtomNumConstIterator searchAtomConst(const std::string &an) const;

        /*! \brief
         * Return iterator pointing to a specific atom or EndAtomNum if not found
         *
         * \param[in] an Atom number
         */
        AtomNumIterator searchAtom(const std::string &an);

        /*! \brief
         * Return the number of atoms of a certain type
         *
         **\param[in] atom Atom name
         */
        int CountAtoms(const std::string &atom) const;

        /*! \brief
         * Return the total number of atoms
         */
        int CountAtoms() const;

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] cr   GROMACS data structure for MPI communication
         * \param[in] dest Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(const CommunicationRecord *cr,
                                 int                        dest) const;

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] cr  GROMACS data structure for MPI communication
         * \param[in] src Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(const CommunicationRecord *cr,
                                    int                        src);
};
//! Iterator over MolecularComposition items
using MolecularCompositionIterator      = typename std::vector<MolecularComposition>::iterator;
//! Const iterator over MolecularComposition items
using MolecularCompositionConstIterator = typename std::vector<MolecularComposition>::const_iterator;

} // namespace alexandria

#endif
