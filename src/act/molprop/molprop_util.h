/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#ifndef MOLPROP_UTIL_H
#define MOLPROP_UTIL_H

#include <string>
#include <vector>

#include "act/alexandria/molselect.h"
#include "act/molprop/molprop.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/arrayref.h"

/*! \brief
 * Enumerated type for MolPropSort function
 *
 * \ingroup group_molprop
 */
enum MolPropSortAlgorithm {
    MPSA_MOLNAME     = 0,
    MPSA_FORMULA     = 1,
    MPSA_COMPOSITION = 2,
    MPSA_SELECTION   = 3,
    MPSA_NR          = 4
};

namespace alexandria
{
/*! \brief Reset the molprop index number
 * \param[in] mp The molprops
 */
void generate_index(std::vector<MolProp> *mp);

/*! \brief Will merge multiple molprop files.
 * \param[in] msg_handler For info and debugging
 * \param[in] infiles     MolProp filenames
 * \param[in] mp          The molprop vector
 */
void merge_xml(MsgHandler                             *msg_handler,
               const gmx::ArrayRef<const std::string> &infiles,
               std::vector<alexandria::MolProp>       *mp);

/*! \brief
 * Sorts a vector of molprops
 *
 * Function that uses the std::sort routine and can apply different sorting
 * keys.
 *
 * \param[in]     msg_handler For info and debugging
 * \param[inout]  mp          The vector of MolProp
 * \param[in]     mpsa        The algorithm used for sorting
 * \param[in]     apt         Database of atom properties
 * \param[in]     gms         Optional structure containing selection criteria
 * \ingroup group_molprop
 */
void MolPropSort(MsgHandler           *msg_handler,
                 std::vector<MolProp> *mp,
                 MolPropSortAlgorithm  mpsa,
                 gmx_atomprop_t        apt,
                 const MolSelect      &gms);

} // namespace alexandria

/*! \brief Utility to split a user-provided lot
 *
 * \param[in]  lot    Level of theory
 * \param[out] method QM method
 * \param[out] basis  QM basis set
 */
void splitLot(const char  *lot,
              std::string *method,
              std::string *basis);

#endif
