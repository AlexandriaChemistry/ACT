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
    
#include "fetch_charges.h"    

#include "act/alexandria/actmol.h"
#include "act/basics/allmols.h"
#include "act/basics/msg_handler.h"
#include "act/molprop/molprop_xml.h"

namespace alexandria
{

ChargeMap fetchChargeMap(MsgHandler                  *msghandler,
                         ForceField                  *pd,
                         const ForceComputer         *forceComp,
                         const std::vector<MolProp>  &mps,
                         const std::set<std::string> &lookup,
                         ChargeGenerationAlgorithm    algorithm,
                         const char                  *qread)
{
    AlexandriaMols amols;
    ChargeMap      qmap;
    for(auto mp = mps.begin(); mp < mps.end(); ++mp)
    {
        auto &frags = mp->fragments();
        // Only do monomers
        if (frags.size() != 1)
        {
            continue;
        }
        bool addThisMol = false;
        // Check lookup table of compounds
        if (!lookup.empty())
        {
            addThisMol = lookup.find(mp->getIupac()) != lookup.end();
            if (!addThisMol)
            {
                // Look for synonyms
                auto amol = amols.findMol(mp->getMolname());
                if (amol)
                {
                    for(const auto &syn : amol->synonyms)
                    {
                        addThisMol = addThisMol || (lookup.find(syn) != lookup.end());
                    }
                }
            }
        }
        if (!addThisMol)
        {
            continue;
        }
        alexandria::ACTMol actmol;
        actmol.Merge(&(*mp));
        actmol.GenerateTopology(msghandler, pd, missingParameters::Error);
        // Add charges according to the algorithm selected
        std::vector<std::pair<Identifier, double>> newq;
        if (!msghandler->ok())
        {
            msghandler->resetStatus();
            continue;
        }
        if (algorithm == ChargeGenerationAlgorithm::Read)
        {
            actmol.setCharges(msghandler, pd, qread);
        }
        else
        {
            if (algorithm == ChargeGenerationAlgorithm::ESP)
            {
                std::map<MolPropObservable, iqmType> iqm = {
                    { MolPropObservable::POTENTIAL, iqmType::QM }
                };
                actmol.getExpProps(msghandler, pd, iqm, 0.0, 100);
            }
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            std::vector<gmx::RVec> coords = actmol.xOriginal();
            actmol.generateCharges(msghandler, pd, forceComp, algorithm,
                                   &coords, &forces);
        }
        if (msghandler->ok())
        {
            for(auto &atom: actmol.atomsConst())
            {
                newq.push_back({atom.id(), atom.charge()});
            }
            qmap.insert( { frags[0].inchi(), newq } );
        }
        else
        {
            msghandler->resetStatus();
        }
    }
    return qmap;
}

ChargeMap fetchChargeMap(MsgHandler                  *msghandler,
                         ForceField                  *pd,
                         const ForceComputer         *forceComp,
                         const char                  *charge_fn,
                         const std::set<std::string> &lookup,
                         ChargeGenerationAlgorithm    algorithm,
                         const char                  *qread)
{
    std::vector<MolProp> mps;
    MolPropRead(msghandler, charge_fn, &mps);
    return fetchChargeMap(msghandler, pd, forceComp, mps, lookup, algorithm, qread);
}

void broadcastChargeMap(const CommunicationRecord *cr,
                        ChargeMap                 *qmap)
{
    auto mycomm = MPI_COMM_WORLD;
    int nqmap   = qmap->size();
    cr->bcast(&nqmap, mycomm);
    if (cr->isMaster())
    {
        for(auto &q : *qmap)
        {
            std::string cmp = q.first;
            cr->bcast(&cmp, mycomm);
            int natom = q.second.size();
            cr->bcast(&natom, mycomm);
            for(size_t i = 0; i < q.second.size(); i++)
            {
                auto &mypair = q.second[i];
                mypair.first.BroadCast(cr, 0, mycomm);
                cr->bcast(&mypair.second, mycomm);
            }
        }
    }
    else
    {
        for(int i = 0; i < nqmap; i++)
        {
            std::string cmp;
            cr->bcast(&cmp, mycomm);
            int natom;
            cr->bcast(&natom, mycomm);
            std::vector<std::pair<Identifier, double>> newpairs;
            for(int i = 0; i < natom; i++)
            {
                std::pair<Identifier, double> mypair;
                mypair.first.BroadCast(cr, 0, mycomm);
                cr->bcast(&mypair.second, mycomm);
                newpairs.push_back(mypair);
            }
            qmap->insert({cmp, newpairs});
        }
    }
}

} // namespace alexandria
