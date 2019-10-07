/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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
 
 
#ifndef PLISTWRAPPER_H
#define PLISTWRAPPER_H

//#include <algorithm>
#include <vector>

#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/hackblock.h"

namespace alexandria
{
//! Interaction type
enum InteractionType
{
    eitBONDS              = 0,
    eitANGLES             = 1,
    eitLINEAR_ANGLES      = 2,
    eitPROPER_DIHEDRALS   = 3,
    eitIMPROPER_DIHEDRALS = 4,
    eitVDW                = 5,
    eitLJ14               = 6,
    eitPOLARIZATION       = 7,
    eitCONSTR             = 8,
    eitVSITE2             = 9,
    eitVSITE3FAD          = 10,
    eitVSITE3OUT          = 11,
    eitNR                 = 12
};

using ParamIterator      = typename std::vector<t_param>::iterator;
using ConstParamIterator = typename std::vector<t_param>::const_iterator;

using BondOrderIterator      = typename std::vector<size_t>::iterator;
using ConstBondOrderIterator = typename std::vector<size_t>::const_iterator;

//! Cleaner version of plist array
class PlistWrapper
{
    public:
        //! Constructor
        PlistWrapper(InteractionType itype,
                     int             ftype) 
                     : 
                         ftype_(ftype), 
                         itype_(itype)     
                     {}

        //! Add one parameter
        void addParam(t_param p) { p_.push_back(p); }

        //! Return the function type
        int getFtype() const { return ftype_; }

        //! Return the interaction type
        InteractionType getItype() const { return itype_; }

        //! Update the function type
        void setFtype(int ftype) { ftype_ = ftype; }

        //! Loop over parameters
        ConstParamIterator beginParam() const { return p_.begin(); }

        //! Loop over parameters
        ConstParamIterator endParam() const { return p_.end(); }

        //! Loop over parameters
        ParamIterator beginParam() { return p_.begin(); }

        //! Loop over parameters
        ParamIterator endParam() { return p_.end(); }

        //! Remove one parameter from the array and return array for next
        ParamIterator eraseParam(ParamIterator p) { return p_.erase(p); }

        //! Remove all parameters
        void eraseParams() { p_.clear(); }

        //! Return number of parameters
        unsigned int nParam() const { return p_.size(); }
        
        //! Set bond order
        void addBondOrder (size_t bondOrder) {bondOrder_.push_back(bondOrder);}
        
        //! Return bond order vector
        std::vector<size_t> bondOrder () const {return bondOrder_;}
        
        //! Return bond order for bond j
        size_t bondOrder (int j) const {return bondOrder_[j];}
        
        //! Loop over parameters
        ConstBondOrderIterator beginBondOrder() const { return bondOrder_.begin(); }

        //! Loop over parameters
        ConstBondOrderIterator endBondOrder() const { return bondOrder_.end(); }

        //! Loop over parameters
        BondOrderIterator beginBondOrder() { return bondOrder_.begin(); }

        //! Loop over parameters
        BondOrderIterator endBondOrder() { return bondOrder_.end(); }
    
    private:
        //! Function type
        int                  ftype_;
        //! Interaction type
        InteractionType      itype_;
        //! Array of parameters
        std::vector<t_param> p_;
        //! Bond order 
        std::vector<size_t>  bondOrder_;
};

//! Another utility typedef for a looper
using  PlistWrapperIterator      = typename std::vector<PlistWrapper>::iterator;
//! Another utility typedef for a looper
using  ConstPlistWrapperIterator = typename std::vector<PlistWrapper>::const_iterator;

ConstPlistWrapperIterator SearchPlist(const std::vector<PlistWrapper> &plist, int ftype);

 ConstPlistWrapperIterator SearchPlist(const std::vector<PlistWrapper> &plist, InteractionType itype);

 PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, int ftype);

 PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, InteractionType itype);

 unsigned int CountPlist(const std::vector<PlistWrapper> &plist, int ftype);

void delete_params(std::vector<PlistWrapper> &plist_,
                   const int                  ftype,
                   const int                  alist[]);
                   

void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        int                        ftype,
                        InteractionType            itype,
                        const t_param             &p);
                        
void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        int                        ftype,
                        InteractionType            itype,
                        const t_param             &p,
                        size_t                     bondOrder);
}

#endif
