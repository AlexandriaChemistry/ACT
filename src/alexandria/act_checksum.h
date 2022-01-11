/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ACT_CHECKSUM_H
#define ACT_CHECKSUM_H

#include <string>

namespace alexandria
{

    class Poldata;
    
    /*! \brief Generate an MD5 checksum of a file
     * \param[in] filename The name of the file
     */
    std::string computeCheckSum(const std::string &filename);

    /*! \brief Will generate an MD5 checksum of the Poldata structure.
     * The Poldata structure will be written to a temporary file with
     * empty internal checksum and timestamp fields. Then the checksum
     * will be computed using the above routine, and the original 
     * checksum and timestamp will be restored. 
     * The integrity of a force field file can then
     * be checked by comparing pd->version() to poldataCheckSum(pd).
     * \param[in] pd Poldata structure. Even though this is mutable it
     *                should be unchanged after the call to this function.
     */
    std::string poldataCheckSum(Poldata *pd);


} // namespace alexandria

#endif
