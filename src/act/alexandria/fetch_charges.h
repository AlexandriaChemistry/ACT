/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023
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
#ifndef FETCH_CHARGES_H
#define FETCH_CHARGES_H

#include <map>
#include <string>
#include <vector>

#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/molprop/molprop.h"
#include "act/utility/communicationrecord.h"

namespace alexandria
{

typedef std::map<std::string, std::vector<std::pair<Identifier, double> > > chargeMap;

/*! \brief Generate charges for all compounds in a molprop file.
 *
 * Will create a map containing the charges on all atoms in "molecules" containing a single fragment.
 * For polarizable models the charge of the shell are added explicitly in the list, and the same
 * goes for virtual sites. The fragment id (InChi) is used as the string in the returned map.
 *
 * \param[in] pd        The force field structure
 * \param[in] forceComp A force computer
 * \param[in] charge_fn The name of a molprop file
 * \param[in] qt        Charge type, by default ACM charges will be generated.
 * \return the map.
 */
chargeMap fetchChargeMap(ForceField    *pd,
                         ForceComputer *forceComp,
                         const char    *charge_fn,
                         qType          qt = qType::ACM);

/*! \brief Generate charges for all compounds in a molprop file.
 * \param[in] pd        The force field structure
 * \param[in] forceComp A force computer
 * \param[in] mps       Vector of molprops
 * \param[in] qt        Charge type, by default ACM charges will be generated.
 * \return the map, see above.
 */
chargeMap fetchChargeMap(ForceField                 *pd,
                         ForceComputer              *forceComp,
                         const std::vector<MolProp> &mps,
                         qType                       qt = qType::ACM);

/*! \brief Broadcast a charge map to processors
 * \param[in]    cr   The communication data structure
 * \param[inout] qmap The charge map
 */
void broadcastChargeMap(const CommunicationRecord *cr,
                        chargeMap                 *qmap);

} // namespace alexandria

#endif
