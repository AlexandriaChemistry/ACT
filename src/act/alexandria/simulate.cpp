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

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cctype>
#include <cstdlib>

#include "act/alexandria/actmol.h"
#include "act/alexandria/alex_modules.h"
#include "act/alexandria/compound_reader.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/fragmenthandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/secondvirial.h"
#include "act/alexandria/train_utility.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/utility/jsontree.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

int simulate(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "alexandria simulate performs a proof-of-principle MD simulation, using",
        "a force field derived by the Alexandria Chemistry Toolkit.", 
        "The program can perform energy minimization and/or do a",
        "constant energy molecular dynamics simulation in vacuum.",
        "In addition, a series of conformations (trajectory) may be",
        "submitted after which the energy per conformation is printed.[PAR]",
        "The input is given by a coordinate file, a force field file and",
        "command line options. During the simulation an energy file,",
        "a trajectory file and a log file are generated.[PAR]",
        "During minimization a user select group of atoms can be frozen.",
        "To do so, supply an index file with atom numbers (first atom is 1)",
        "and numbering should disregard shells and virtual sites if present",
        "in the model."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff",      "aff",        ffREAD  },
        { efPDB, "-o",       "trajectory", ffWRITE },
        { efSTO, "-c",       "confout",    ffOPTWR },
        { efXVG, "-e",       "energy",     ffWRITE },
        { efNDX, "-freeze",  "freeze",     ffOPTRD }
    };
    gmx_output_env_t         *oenv;
    double                    shellToler = 1e-6;
    bool                      json       = false;
    std::vector<t_pargs>      pa = {
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)" },
        { "-json", FALSE, etBOOL, {&json},
          "Print part of the output in json format" }
    };
    SimulationConfigHandler  sch;
    sch.add_options(&pa, &fnm);
    sch.add_MD_options(&pa);
    DimerGenerator           gendimers;

    gendimers.addOptions(&pa, &fnm, &desc);
    ReRunner                 rerun(false);
    rerun.addOptions(&pa, &fnm);
    CompoundReader compR;
    compR.addOptions(&pa, &fnm, &desc);
    MsgHandler msghandler;
    msghandler.addOptions(&pa, &fnm, "simulate");
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 1;
    }
    CommunicationRecord cr;
    cr.init(cr.size());
    msghandler.optionsFinished(fnm, &cr);

    sch.check_pargs(&msghandler);
    compR.optionsOK(&msghandler, fnm);
    if (!msghandler.ok())
    {
        return 1;
    }
    gendimers.finishOptions(fnm);

    print_header(msghandler.tw(), pa, fnm);

    if (shellToler >= sch.forceTolerance())
    {
        shellToler = sch.forceTolerance()/10;
        msghandler.msg(ACTStatus::Warning,
                       gmx::formatString("\nShell tolerance larger than atom tolerance, changing it to %g\n", shellToler));
    }
    ForceField        pd;
    try
    {
        readForceField(opt2fn("-ff", fnm.size(), fnm.data()), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    (void) pd.verifyCheckSum(stderr);

    auto forceComp = new ForceComputer(shellToler, sch.maxIter());
    std::vector<ACTMol> actmols = compR.read(&msghandler, pd, forceComp);
    if (actmols.empty())
    {
        fprintf(stderr, "Could not read or process molecules. Please check the log file for more information.\n");
        return 1;
    }
    auto &actmol = actmols[0];
    JsonTree jtree("simulate");
    if (msghandler.info())
    {
        forceFieldSummary(&jtree, &pd);
    }

    if (pd.polarizable())
    {
        // Make a copy since it maybe changed
        auto xx     = actmol.xOriginal();
        auto qprops = actmol.qProps();
        for(auto qp = qprops->begin(); qp < qprops->end(); ++qp)
        {
            auto qCalc = qp->qPact();
            qCalc->initializeMoments();
            qCalc->calcPolarizability(&pd, actmol.topology(), forceComp);
            auto alpha = qCalc->polarizabilityTensor();
            std::string unit("A^3");
            double fac = convertFromGromacs(1, unit);
            JsonTree poltree("Polarizability");

            poltree.addValueUnit("XX", gmx_ftoa(fac*alpha[XX][XX]), unit);
            poltree.addValueUnit("YY", gmx_ftoa(fac*alpha[YY][YY]), unit);
            poltree.addValueUnit("ZZ", gmx_ftoa(fac*alpha[ZZ][ZZ]), unit);
            poltree.addValueUnit("Average", gmx_ftoa(fac*qCalc->isotropicPolarizability()), unit);
            jtree.addObject(poltree);
        }
    }

    if (debug)
    {
        actmol.topology()->dump(debug);
    }
    auto eMin = eMinimizeStatus::OK;
    /* Generate output file for debugging if requested */
    if (gendimers.hasTrajectory())
    {
        rerun.setFunctions(forceComp, &gendimers, oenv);
        rerun.setEInteraction(actmol.fragmentHandler()->topologies().size() > 1);
        rerun.rerun(&msghandler, &pd, &actmol, msghandler.verbose());
    }
    else
    {
        MolHandler molhandler;
        std::vector<gmx::RVec> coords = actmol.xOriginal();
        std::vector<gmx::RVec> xmin   = coords;
        auto nfrag = actmol.fragmentHandler()->topologies().size();
        msghandler.msg(ACTStatus::Info, gmx::formatString("There are %lu fragments.", nfrag));
        if (sch.singlePoint())
        {
            std::map<InteractionType, double> energies;
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            (void) forceComp->compute(&pd, actmol.topology(), &xmin, &forces, &energies);
            JsonTree jtener("Energies before");
            std::string unit("kJ/mol");
            for (const auto &ener : energies)
            {
                auto val = gmx::formatString("%.4f", ener.second);
                jtener.addValueUnit(interactionTypeToString(ener.first),
                                    val.c_str(), unit);
            }
            jtree.addObject(jtener);
            if (nfrag == 2)
            {
                // Make a copy to keep the correct shell coordinates when saving the structure.
                auto newxmin = xmin;
                std::map<InteractionType, double> einter;
                std::vector<gmx::RVec>            interactionForces;
                actmol.calculateInteractionEnergy(&msghandler, &pd, forceComp, &einter,
                                                  &interactionForces, &newxmin, true);
                for(const auto &ei : einter)
                {
                    msghandler.msg(ACTStatus::Info,
                                   gmx::formatString("Interaction energy %s: %g",
                                                     interactionTypeToString(ei.first).c_str(), ei.second));
                }
            }
        }
        else if (sch.minimize())
        {
            std::vector<int> freeze;
            auto freezeName = opt2fn_null("-freeze", fnm.size(),fnm.data());
            if (nullptr != freezeName)
            {
                int isize;
                int *myindex;
                char *grpnames;
                rd_index(freezeName, 1, &isize, &myindex, &grpnames);
                if (isize > 0)
                {
                    printf("Will freeze %d atoms %s\n", isize, grpnames);
                    for(int ii = 0; ii < isize; ii++)
                    {
                        freeze.push_back(myindex[ii]);
                    }
                }
            }

            std::map<InteractionType, double> energies;
            {
                std::vector<gmx::RVec> forces(actmol.atomsConst().size());
                (void) forceComp->compute(&pd, actmol.topology(), &xmin, &forces, &energies);
                JsonTree jtener("Energies before");
                std::string unit("kJ/mol");
                for (const auto &ener : energies)
                {
                    auto val = gmx::formatString("%.4f", ener.second);
                    jtener.addValueUnit(interactionTypeToString(ener.first),
                                        val.c_str(), unit);
                }
                jtree.addObject(jtener);
            }
            eMin = molhandler.minimizeCoordinates(&msghandler, &pd, &actmol, forceComp, sch,
                                                  &xmin, &energies, freeze);
            if (eMinimizeStatus::OK == eMin)
            {
                auto rmsd = molhandler.coordinateRmsd(&actmol, coords, &xmin);
                msghandler.msg(ACTStatus::Info,
                               gmx::formatString("Final energy: %g RMSD wrt original structure %g nm.",
                                                 energies[InteractionType::EPOT], rmsd));
                if (nfrag == 2)
                {
                    // Make a copy to keep the correct shell coordinates when saving the structure.
                    auto newxmin = xmin;
                    std::map<InteractionType, double> einter;
                    std::vector<gmx::RVec>            interactionForces;
                    actmol.calculateInteractionEnergy(&msghandler, &pd, forceComp, &einter,
                                                      &interactionForces, &newxmin, true);
                    for(const auto &ei : einter)
                    {
                        msghandler.msg(ACTStatus::Info,
                                       gmx::formatString("Interaction energy %s: %g",
                                                         interactionTypeToString(ei.first).c_str(), ei.second));
                    }
                }
                JsonTree jtener("Energies after");
                std::string unit("kJ/mol");
                for (const auto &ener : energies)
                {
                    auto val = gmx::formatString("%.4f", ener.second);
                    jtener.addValueUnit(interactionTypeToString(ener.first),
                                        val.c_str(), unit);
                }
                jtree.addObject(jtener);
            }
        }
        if (eMinimizeStatus::OK == eMin && sch.nsteps() > 0)
        {
            molhandler.simulate(&msghandler, &pd, &actmol, forceComp, sch,
                                opt2fn("-o", fnm.size(),fnm.data()),
                                opt2fn("-e", fnm.size(),fnm.data()),
                                oenv);
        }
        auto confout = opt2fn_null("-c", fnm.size(),fnm.data());
        if (confout && eMinimizeStatus::OK == eMin)
        {
            matrix box = { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } };
            actmol.PrintConformation(confout, xmin, sch.writeShells(), box);
        }
    }

    auto imm = ACTMessage::OK;
    if (eMinimizeStatus::OK != eMin)
    {
        fprintf(stderr, "Minimization failed: %s, check log file %s\n",
                eMinimizeStatusToString(eMin).c_str(),
                msghandler.filename().c_str());
        status = 1;
    }
    else if (ACTMessage::OK != imm)
    {
        msghandler.msg(ACTStatus::Error, ACTMessage::MinimizationFailed,
                       gmx::formatString("Please check the log file %s for error messages.\n",
                                         msghandler.filename().c_str()));
        status = 1;
    }
    if (json)
    {
        jtree.write("simulate.json", json);
    }
    else
    {
        int indent = 0;
        msghandler.write(jtree.writeString(json, &indent));
    }
    if (!msghandler.ok())
    {
        printf("Simulate failed with an error:\n");
    }
    return status;
}

} // namespace alexandria
