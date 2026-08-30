/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2026
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
/*! \internal \file
 * \brief
 * Implements test of bonded force routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_act_tests
 */
#include "actpre.h"

#include "actmol_util.h"
#include "act/import/import.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/topology/actmol.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

void initACTMol(MsgHandler          *msghandler,
                const std::string   &dataName,
                ForceField          *pd,
                ForceComputer       *fcomp,
                std::vector<ACTMol> *mps)
{
    if (!msghandler)
    {
        fprintf(stderr, "Cannot initACTMol without MsgHandler\n");
        return;
    }
    const char   *conf     = (char *)"minimum";
    //    std::string   dataName = gmx::test::TestFileManager::getInputFilePath(molname);
    std::vector<alexandria::MolProp> molprops;
    double        qtot     = 0;
    // Charge gen params
    auto alg = pd->chargeGenerationAlgorithm();
    std::vector<double> qcustom;
    bool userqtot = false;
    importFile(msghandler, pd, dataName.c_str(), &molprops,
               conf, JobType::OPT, userqtot,
               &qtot, true);
    if (!msghandler->ok())
    {
        msghandler->fatal("Could not import file");
    }

    for(auto &molprop : molprops)
    {
        ACTMol mm;
        mm.Merge(&molprop);
        mm.GenerateTopology(msghandler, pd,
                            missingParameters::Ignore);
        if (!msghandler->ok())
        {
            msghandler->fatal("Could not generate topology");
        }
        std::map<MolPropObservable, iqmType> iqmMap = 
            {
                { MolPropObservable::DELTAE0,           iqmType::QM },
                { MolPropObservable::POTENTIAL,         iqmType::QM },
                { MolPropObservable::INTERACTIONENERGY, iqmType::QM },
                { MolPropObservable::DIPOLE,            iqmType::QM },
                { MolPropObservable::QUADRUPOLE,        iqmType::QM },
                { MolPropObservable::OCTUPOLE,          iqmType::QM },
                { MolPropObservable::HEXADECAPOLE,      iqmType::QM },
                { MolPropObservable::POLARIZABILITY,    iqmType::QM },
                { MolPropObservable::CHARGE,            iqmType::QM }
            };
        mm.getExpProps(msghandler, pd, iqmMap, 0);
        if (msghandler->ok())
        {
            std::vector<gmx::RVec> forces(mm.atomsConst().size());
            std::vector<gmx::RVec> coords = mm.xOriginal();
            mm.generateCharges(msghandler, pd, fcomp, alg, &coords, &forces, true);
        }
        if (msghandler->ok())
        {
            mps->push_back(std::move(mm));
        }
    }
}

}
