/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GA_DATASET_H
#define GA_DATASET_H

#include <map>
#include <string>

//! Class to distinguish data sets
enum class iMolSelect {
    //! Training data set
    Train,
    //! Testing data set
    Test,
    //! Data set to ignore
    Ignore
};

/*! \brief Return string corresponding to data set
 * \param[in] ims The data set
 * \return a string
 */
const char *iMolSelectName(iMolSelect ims);

//! \return map of all names of data sets
const std::map<iMolSelect, const char *> &iMolSelectNames();

/*! \brief Look up data set name
 * \param[in]  name The data set name
 * \param[out] ims  The type of data set
 * \return true if found, false otherwise
 */
bool name2molselect(const std::string &name, iMolSelect *ims);

#endif
