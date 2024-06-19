/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023,2024
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
    
#include "fetch_charges.h"    

#include "act/alexandria/actmol.h"
#include "act/molprop/molprop_xml.h"

namespace alexandria
{

chargeMap fetchChargeMap(ForceField                 *pd,
                         ForceComputer              *forceComp,
                         const std::vector<MolProp> &mps,
                         qType                       qt)
{
    chargeMap qmap;
    for(auto mp = mps.begin(); mp < mps.end(); mp++)
    {
        alexandria::ACTMol actmol;
        actmol.Merge(&(*mp));
        auto imm = actmol.GenerateTopology(nullptr, pd, missingParameters::Error);
        if (immStatus::OK != imm)
        {
            continue;
        }
        std::vector<gmx::RVec> coords = actmol.xOriginal();
        std::map<MolPropObservable, iqmType> iqm = {
            { MolPropObservable::CHARGE, iqmType::QM }
        };
        actmol.getExpProps(pd, iqm, 0.0, 100);
        auto fhandler = actmol.fragmentHandler();
        if (fhandler->topologies().size() == 1)
        {
            std::vector<double> dummy;
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            auto alg = pd->chargeGenerationAlgorithm();
            if (qType::ACM != qt)
            {
                alg = ChargeGenerationAlgorithm::Read;
            }
            imm = actmol.GenerateCharges(pd, forceComp, alg,
                                         qt, dummy, &coords, &forces);
            if (immStatus::OK == imm)
            {
                // Add ACM charges
                std::vector<std::pair<Identifier, double>> newq;
                for(auto atom: actmol.atomsConst())
                {
                    newq.push_back({atom.id(), atom.charge()});
                }
                qmap.insert({fhandler->ids()[0], newq});
            }
        }
    }
    return qmap;
}

chargeMap fetchChargeMap(ForceField    *pd,
                         ForceComputer *forceComp,
                         const char    *charge_fn,
                         qType          qt)
{
    std::vector<MolProp> mps;
    MolPropRead(charge_fn, &mps);
    return fetchChargeMap(pd, forceComp, mps, qt);
}

void broadcastChargeMap(const CommunicationRecord *cr,
                        chargeMap                 *qmap)
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
