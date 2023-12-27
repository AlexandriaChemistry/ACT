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
    
#include "fetch_charges.h"    

#include "act/alexandria/actmol.h"
#include "act/molprop/molprop_xml.h"

namespace alexandria
{

std::map<std::string, std::vector<double> > fetchChargeMap(const ForceField           *pd,
                                                           ForceComputer              *forceComp,
                                                           const std::vector<MolProp> &mps)
{
    std::map<std::string, std::vector<double> > qmap;
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
        actmol.getExpProps(pd, iqm, 0.0, 0.0, 100);
        auto fhandler = actmol.fragmentHandler();
        if (fhandler->topologies().size() == 1)
        {
            std::vector<double> dummy;
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            imm = actmol.GenerateCharges(pd, forceComp, pd->chargeGenerationAlgorithm(),
                                         qType::ACM, dummy, &coords, &forces);
            if (immStatus::OK == imm)
            {
                // Add ACM charges
                std::vector<double> newq;
                for(auto atom: actmol.atomsConst())
                {
                    newq.push_back(atom.charge());
                }
                qmap.insert({fhandler->ids()[0], newq});
            }
        }
    }
    return qmap;
}

std::map<std::string, std::vector<double> > fetchChargeMap(const ForceField *pd,
                                                           ForceComputer    *forceComp,
                                                           const char       *charge_fn)
{
    std::vector<MolProp> mps;
    MolPropRead(charge_fn, &mps);
    return fetchChargeMap(pd, forceComp, mps);
}

void broadcastChargeMap(const CommunicationRecord                   *cr,
                        std::map<std::string, std::vector<double> > *qmap)
{
    auto mycomm = MPI_COMM_WORLD;
    int nqmap   = qmap->size();
    cr->bcast(&nqmap, mycomm);
    if (cr->isMaster())
    {
        for(auto &q : *qmap)
        {
            std::string tmp = q.first;
            cr->bcast(&tmp, mycomm);
            cr->bcast(&q.second, mycomm);
        }
    }
    else
    {
        for(int i = 0; i < nqmap; i++)
        {
            std::string id;
            cr->bcast(&id, mycomm);
            std::vector<double> q;
            cr->bcast(&q, mycomm);
            qmap->insert({id, q});
        }
    }
}

} // namespace alexandria
