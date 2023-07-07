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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cctype>
#include <cstdlib>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/actmol.h"
#include "act/alexandria/tuning_utility.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

int min_complex(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "alexandria min_complex read a molprop file, and for each complex",
        "(more than 1 fragment) in the file it will take the lowest energy",
        "structure, minimize the coordinates and store them."
    };
    std::vector<t_filenm>     fnm = {
        { efXML, "-ff", "aff",      ffREAD  },
        { efXML, "-mp", "molprop",  ffREAD  },
        { efLOG, "-g",  "gen_scan", ffWRITE }
    };
    gmx_output_env_t *oenv;
    double shellToler = 1e-6;
    bool   verbose    = false;
    double minfrac    = 0.9;
    double maxfrac    = 1.5;
    int    nfrac      = 5;
    std::vector<t_pargs>      pa = {
        { "-minfrac",   FALSE, etREAL, {&minfrac},
          "Fraction of center of mass distance to start the scan on" },
        { "-maxfrac",   FALSE, etREAL, {&maxfrac},
          "Fraction of center of mass distance to end the scan on" },
        { "-nfrac", FALSE, etINT, {&nfrac},
          "Number of distances within the range from minfrac to maxfrac" },
        { "-v", FALSE, etBOOL, {&verbose},
          "Print more information to the log file." },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)" }
    };
    SimulationConfigHandler  sch;
    sch.add_pargs(&pa, true);

    int status = 0;
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        status = 1;
        return status;
    }
    sch.check_pargs();

    ForceField        pd;
    try
    {
        readForceField(opt2fn("-ff", fnm.size(), fnm.data()), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    (void) pd.verifyCheckSum(stderr);
    const char *logFileName = opt2fn("-g", fnm.size(),fnm.data());
    FILE *logFile   = gmx_ffopen(logFileName, "w");
    if (shellToler >= sch.forceTolerance())
    {
        shellToler = sch.forceTolerance()/10;
        fprintf(logFile, "Shell tolerance larger than atom tolerance, changing it to %g\n", shellToler);
    }
    auto  forceComp = new ForceComputer(shellToler, 100);
    print_header(logFile, pa, fnm);
    
    std::vector<MolProp> mps;
    MolPropRead(opt2fn("-mp", fnm.size(), fnm.data()), &mps);
    for(auto mp = mps.begin(); mp < mps.end(); ++mp)
    {
        if (mp->fragments().size() != 2)
        {
            fprintf(logFile, "Ignoring '%s' with %zu fragments\n",
                    mp->getMolname().c_str(), mp->fragments().size());
            continue;
        }
        ACTMol actmol;
        matrix box;
        clear_mat(box);
        actmol.Merge(&(*mp));
    
        immStatus imm = immStatus::OK;
        imm = actmol.GenerateTopology(logFile, &pd, missingParameters::Error);
        std::vector<gmx::RVec> coords = actmol.xOriginal();
        if (immStatus::OK == imm)
        {
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            
            std::vector<double> myq;
            auto alg   = pd.chargeGenerationAlgorithm();
            auto qtype = qType::Calc;
            imm        = actmol.GenerateCharges(&pd, forceComp, alg, qtype, myq, &coords, &forces);
        }
        if (immStatus::OK == imm)
        {
            auto eMin = eMinimizeStatus::OK;
            /* Generate output file for debugging if requested */
            MolHandler molhandler;
            std::vector<gmx::RVec> xmin   = coords;
            std::map<InteractionType, double> energies;
            eMin = molhandler.minimizeCoordinates(&pd, &actmol, forceComp, sch,
                                                  &xmin, &energies, 
                                                  verbose ? logFile : nullptr, {});
    
            if (eMinimizeStatus::OK == eMin)
            {
                std::vector<gmx::RVec> interactionForces;
                std::map<InteractionType, double> einter;
                actmol.calculateInteractionEnergy(&pd, forceComp, &einter, &interactionForces, &xmin);
                auto rmsd = molhandler.coordinateRmsd(&actmol, coords, &xmin);
                fprintf(logFile, "%s final energy: %g. Interaction energy: %g. RMSD wrt original structure %g nm.\n",
                        actmol.getMolname().c_str(),
                        energies[InteractionType::EPOT], 
                        einter[InteractionType::EPOT], rmsd);
                std::string confout = gmx::formatString("%s.pdb",
                                                        actmol.getMolname().c_str());
                actmol.PrintConformation(confout.c_str(), xmin, false, box);
            }
            else
            {
                fprintf(stderr, "Minimization failed for %s: %s, check log file %s\n",
                        actmol.getMolname().c_str(),
                        eMinimizeStatusToString(eMin).c_str(),
                        logFileName);
            }
        }
    }
    gmx_ffclose(logFile);
    return status;
}

} // namespace alexandria
