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
#ifndef MOLSELECT_H
#define MOLSELECT_H

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "ga/Dataset.h"

//#include "act/molprop/molprop.h"
//#include "act/poldata/poldata.h"

namespace alexandria
{

class IMolSelect 
{
    private:
        std::string iupac_;
        iMolSelect  status_;
        int         index_;
        
    public:
        IMolSelect(const std::string &iupac, iMolSelect status, int index)
            :
                iupac_(iupac), 
                status_(status), 
                index_(index) 
            {}

        const std::string &iupac() const { return iupac_; }
        
        iMolSelect status() const { return status_; }
        
        int index() const { return index_; }
};

class MolSelect
{
    private:
        //! Vector of iMolSelect entries
        std::vector<IMolSelect> ims_;

    public:
    
        MolSelect() {};

        void read(const char *filename);
        
        size_t nMol() const { return ims_.size(); }
    
        /*! \brief Get data set for iupac
         * \param[in]  iupac The molecule name
         * \param[out] ims   The data set
         * \return true if found, false otherwise
         */    
        bool status(const std::string &iupac, iMolSelect *ims) const;
        
        /*! \brief Get index for iupac
         * \param[in]  iupac The molecule name
         * \param[out] index The index
         * \return true if found, false otherwise
         */    
        bool index(const std::string &iupac, int *index) const;
        
        int count(iMolSelect ims) const
        {
            return std::count_if(ims_.begin(), ims_.end(),
                                 [ims](IMolSelect const &i)
                                 { return i.status() == ims; });
        }
};

} // namespace

#endif
