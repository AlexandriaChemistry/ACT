/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019 
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
 
 
#ifndef GMX_RESP_ATOM_H
#define GMX_RESP_ATOM_H

#include <vector>

#include "gromacs/math/vectypes.h"

#include "poldata.h"

namespace alexandria
{
class RespAtomType
{
    public:
        RespAtomType(int                             atype,
                     int                             particleType,
                     bool                            hasShell,
                     const char                     *atomtype,
                     const Poldata                  *pd,
                     const std::vector<std::string> &dzatoms,
                     double                          zeta,
                     double                          q);

        size_t getNZeta() const { return rz_.size(); }

        const std::string &getAtomtype() const { return atomtype_; }

        void setAtomtype(const std::string &type) { atomtype_ = type; }

        bool getBRestrained() const { return bRestrained_; }

        bool hasShell() const { return bHasShell_; }
        
        int getAtype() const { return atype_; }
        
        int ptype() const { return ptype_; }

        void setAtype(int i) { atype_ = i; }

        RowZetaQIterator beginRZ() { return rz_.begin(); }

        RowZetaQIterator endRZ() { return rz_.end(); }

        RowZetaQConstIterator beginRZ() const { return rz_.begin(); }

        RowZetaQConstIterator endRZ() const { return rz_.end(); }

    private:
        //! Atom type index
        int                   atype_;
        //! particleType
        int                   ptype_;
        //! Signify whether this charge should be restrained during fitting
        bool                  bRestrained_;
        //! Signify whether this particle has a shell
        bool                  bHasShell_;
        //! String corresponding to atom type
        std::string           atomtype_;
        //! Arrays of charge components
        std::vector<RowZetaQ> rz_;
};

using RespAtomTypeIterator      = typename std::vector<RespAtomType>::iterator;
using RespAtomTypeConstIterator = typename std::vector<RespAtomType>::const_iterator;

class RespAtom
{
    public:
 RespAtom(int atomnumber, int atype, int ptype, double q, double qref, gmx::RVec x)
         : atomnumber_(atomnumber), atype_(atype), pType_(ptype), q_(q), qref_(qref), x_(x)
        {
            qindex_ = -1;
            fixedQ_ = (pType_ != eptAtom);
        }

        //! Return the atom type
        int atype() const { return atype_; }

        //! Set the atom type
        void setAtype(int i) { atype_ = i; }

        //! Return the atomic number
        int atomnumber() const { return atomnumber_; }

        //! Set the atomic number
        void setAtomnumber(int i) { atomnumber_ = i; }

        //! Return whether this charge is fixed
        bool fixedQ() const { return fixedQ_; }

        //! Return the charge
        double q() const { return q_; }

        //! Set the charge
        void setQ(double value)
        {
            GMX_RELEASE_ASSERT(!fixedQ_, "Trying to modify a fixed charge");
            q_ = value;
        }
        
        //! Return the reference charge (the non-flexible part)
        double qRef() const { return qref_; }

        //! Return the coordinates
        const gmx::RVec &x() const { return x_; }

        //! Set the coordinates
        void setX(const gmx::RVec &x) { x_ = x; }

        //! Return the charge index in parameterization array
        int qIndex() const { return qindex_; }

        //! Set the charge index
        void setQindex(int qi) { qindex_ = qi; }
    private:
        //! Atomic number
        int       atomnumber_;
        //! Atom type
        int       atype_;
        //! Particle type
        int       pType_;
        //! Total charge of the atom (which is optimized in resp)
        double    q_;
        //! Reference charge
        double    qref_;
        //! Coordinates for this atom
        gmx::RVec x_;
        //! Index in parameterization array
        int       qindex_;
        //! Tells us whether this charge is fixed
        bool      fixedQ_;
};

using RespAtomIterator      = typename std::vector<RespAtom>::iterator;
using RespAtomConstIterator = typename std::vector<RespAtom>::const_iterator;

} // namespace

#endif
