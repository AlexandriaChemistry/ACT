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
 
#ifndef PLISTWRAPPER_H
#define PLISTWRAPPER_H

#include <vector>

#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/hackblock.h"

#include "interactiontype.h"
#include "mymol_low.h"
#include "poldata.h"

namespace alexandria
{
//! Shortcut for vector iterator
using ParamIterator      = typename std::vector<t_param>::iterator;

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
        void addParam(t_param p,
                      const std::vector<double> bondOrder);

        //! Return the function type
        int getFtype() const { return ftype_; }

        //! Return the interaction type
        InteractionType getItype() const { return itype_; }

        //! Update the function type
        void setFtype(int ftype) { ftype_ = ftype; }

        //! Return  the parameter array for editing
        std::vector<t_param> *params() { return &p_; }
        
        //! Return  the parameter array for reading only
        const std::vector<t_param> &paramsConst() const { return p_; }
        
        //! Remove one parameter from the array and return array for next
        ParamIterator eraseParam(ParamIterator p) { return p_.erase(p); }

        //! Remove all parameters
        void eraseParams() { p_.clear(); }

        //! Return number of parameters
        unsigned int nParam() const { return p_.size(); }
        
        //! Return bond order vector
        const std::vector<double> &bondOrder () const {return bondOrder_;}
        
        //! Return bond order for bond j
        double bondOrder (int j) const {return bondOrder_[j];}
        
    private:
        //! Function type
        int                  ftype_;
        //! Interaction type
        InteractionType      itype_;
        //! Array of parameters
        std::vector<t_param> p_;
        //! Bond order 
        std::vector<double>  bondOrder_;
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

void delete_params(std::vector<PlistWrapper> *plist_,
                   const int                  ftype,
                   const int                  alist[]);
                   
void add_param_to_plist(std::vector<PlistWrapper> *plist,
                        int                        ftype,
                        InteractionType            itype,
                        const t_param             &p,
                        const std::vector<double> &bondOrder);

void cp_plist(t_params                   plist[],
              int                        ftype,
              InteractionType            itype,
              std::vector<PlistWrapper> &plist_);

/*! \brief Update GROMACS force field parameters
 *
 * The force field parameters in gromacs idef-style structures are updated from the
 * Poldata structure. This can be a heavy calculation for a large system, but it will 
 * only handle the  bonds etc. in the system (a molecule typically).
 * \param[in]  pd      Poldata structure
 * \param[out] plist   The parameter vector
 * \param[in]  atoms   GROMACS atom structure
 * \param[in]  missing Determining whether parameters can be missing
 * \param[in]  molname Molecule name
 * \param[out] errors  Error messages are appended to this vector.
 * \return Warning status
 */ 
immStatus updatePlist(const Poldata             *pd,
                      std::vector<PlistWrapper> *plist,
                      const t_atoms             *atoms,
                      missingParameters          missing,
                      const std::string         &molname,
                      std::vector<std::string>  *errors);



void print_top_header(FILE                    *fp,
                      const Poldata           *pd,
                      bool                     bPol,
                      std::vector<std::string> commercials,
                      bool                     bItp);

}
#endif
