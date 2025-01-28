/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2025
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
 * \ingroup module_listed-forces
 */
#include "actpre.h"
#include "actmol_util.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/atype_mapping.h"
#include "act/alexandria/actmol.h"
#include "act/basics/msg_handler.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

void initACTMol(const char          *molname, 
                ForceField          *pd,
                ForceComputer       *fcomp,
                std::vector<ACTMol> *mps)
    {
        int           maxpot   = 100;
        int           nsymm    = 0;
        const char   *conf     = (char *)"minimum";
        std::string   method, basis;
        const char   *jobtype  = (char *)"Opt";
        
        std::string   dataName = gmx::test::TestFileManager::getInputFilePath(molname);
        std::vector<alexandria::MolProp> molprops;
        double        qtot     = 0;
        // Charge gen params
        auto alg = ChargeGenerationAlgorithm::NONE;
        std::vector<double> qcustom;
        bool userqtot = false;
        matrix box;
        bool readOK = readBabel(pd, dataName.c_str(), &molprops, molname, molname,
                                conf, &method, &basis,
                                maxpot, nsymm, jobtype, userqtot,
                                &qtot, false, box, true);
        EXPECT_TRUE(readOK);
        if (readOK)
        {
            for(auto &molprop : molprops)
            {
                ACTMol mm;
                mm.Merge(&molprop);
                auto imm = mm.GenerateTopology(stdout, pd,
                                               missingParameters::Ignore);
                EXPECT_TRUE(ACTMessage::OK ==imm);
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
                imm = mm.getExpProps(pd, iqmMap, 0);
                std::vector<gmx::RVec> forces(mm.atomsConst().size());
                std::vector<gmx::RVec> coords = mm.xOriginal();
                mm.GenerateCharges(pd, fcomp, alg, qType::Calc, qcustom, &coords, &forces, true);
                mps->push_back(mm);
            }
        }
    }    
}
