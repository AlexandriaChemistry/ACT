/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024,2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
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

#include "act/basics/dataset.h"
#include "act/utility/communicationrecord.h"

namespace alexandria
{

/*! \brief One entry in a molecule selection: a compound name,
 * the data set it belongs to and its index in the selection
 */
class IMolSelect 
{
    private:
        //! The molecule name
        std::string iupac_;
        //! The data set this molecule belongs to
        iMolSelect  status_;
        //! Sequence number of the molecule in the selection
        int         index_;
        
    public:
        //! Constructor
        IMolSelect(const std::string &iupac, iMolSelect status, int index)
            :
                iupac_(iupac), 
                status_(status), 
                index_(index) 
            {}

        //! \return the molecule name
        const std::string &iupac() const { return iupac_; }
        
        //! \return the data set this molecule belongs to
        iMolSelect status() const { return status_; }
        
        //! \return the sequence number of the molecule in the selection
        int index() const { return index_; }

        //! \return true if this is in fact a dimer
        bool isDimer() const { return iupac_.find("#") != std::string::npos; };
};

/*! \brief Selection of molecules, each assigned to a data set
 * (e.g. training or test), typically read from a selection file
 */
class MolSelect
{
private:
    //! Vector of iMolSelect entries
    std::vector<IMolSelect> ims_;
    //! Look up iupac
    const IMolSelect *findIupac(const std::string &iupac) const;
public:
    
    //! Constructor
    MolSelect() = default;
    /*! \brief Add one entry
     * \param[in] iupac The molecule name
     * \param[in] index A number
     * \param[in] ims   The selection group
     */
    void addOne(const std::string &iupac,
                int                index,
                iMolSelect         ims);
    /*! \brief Read selection from a file
     * \param[in] filename
     */
    void read(const char *filename);
    
    //! \return the number of molecules in the selection
    size_t nMol() const { return ims_.size(); }

    /*! \brief Get data set for iupac
     * \param[in]  iupac The molecule name
     * \param[out] ims   The data set
     * \return true if found, false otherwise
     */    
    bool status(const std::string &iupac, iMolSelect *ims) const;

    //! \return the vector of selected compounds
    const std::vector<IMolSelect> &imolSelect() const { return ims_; }
    /*! \brief Get index for iupac
     * \param[in]  iupac The molecule name
     * \param[out] index The index
     * \return true if found, false otherwise
     */    
    bool index(const std::string &iupac, int *index) const;
    
    /*! \brief Count the molecules in a data set
     * \param[in] ims The data set
     * \return the number of molecules in that data set
     */
    int count(iMolSelect ims) const
    {
        return std::count_if(ims_.begin(), ims_.end(),
                             [ims](IMolSelect const &i)
                                 { return i.status() == ims; });
    }

    //! \return number of dimers in the data selection
    size_t countDimers() const
    {
        return std::count_if(ims_.begin(), ims_.end(),
                             [](IMolSelect const &i)
                                 { return i.isDimer(); });
    }

    //! Broadcast my data
    void bcast(const CommunicationRecord *cr);
};

} // namespace

#endif
