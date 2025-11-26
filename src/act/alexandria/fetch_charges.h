/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2025
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
#ifndef FETCH_CHARGES_H
#define FETCH_CHARGES_H

#include <map>
#include <set>
#include <string>
#include <vector>

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/molprop/molprop.h"
#include "act/utility/communicationrecord.h"

namespace alexandria
{

class ACTMol;

/*! \brief chargeMap definition.
 * First string is the molecule identifier, which should be the InChi.
 * Vector contains pairs of identifier and charge.
 */
typedef std::map<std::string, std::vector<std::pair<Identifier, double> > > ChargeMap;

/*! \brief Generate charges for all compounds in a molprop file.
 *
 * Will create a map containing the charges on all atoms in "molecules" containing a single fragment.
 * For polarizable models the charge of the shell are added explicitly in the list, and the same
 * goes for virtual sites. The fragment id (InChi) is used as the string in the returned map.
 *
 * \param[in]  msghandler The message handler
 * \param[in]  pd         The force field structure
 * \param[in]  forceComp  A force computer
 * \param[in]  charge_fn  The name of a molprop file
 * \param[out] mols       The molecules
 * \param[in]  lookup     Set of compounds to look up. If empty charges for all compounds will be determined.
 * \param[in]  algorithm  Method to generate charges
 * \param[in]  qread      Charge type if read
 * \return the map.
 */
ChargeMap fetchChargeMap(MsgHandler                  *msghandler,
                         ForceField                  *pd,
                         const ForceComputer         *forceComp,
                         const char                  *charge_fn,
                         std::vector<ACTMol>         *mols,
                         const std::set<std::string> &lookup,
                         ChargeGenerationAlgorithm    algorithm,
                         const char                  *qread);

/*! \brief Broadcast a charge map to processors
 * \param[in]    cr   The communication data structure
 * \param[inout] qmap The charge map
 */
void broadcastChargeMap(const CommunicationRecord *cr,
                        ChargeMap                 *qmap);

} // namespace alexandria

#endif
