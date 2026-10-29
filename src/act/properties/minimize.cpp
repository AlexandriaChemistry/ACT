/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
#include "minimize.h"

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/train_utility.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/import/compound_reader.h"
#include "act/molprop/molprop_xml.h"
#include "act/properties/secondvirial.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_low.h"
#include "act/topology/fragmenthandler.h"
#include "act/utility/jsontree.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

//! \brief Add energies to a json tree
static void addEnergies(JsonTree                                *jtree,
                        const char                              *name,
                        const std::map<InteractionType, double> &energies)
{
    JsonTree    jtener(name);
    std::string unit("kJ/mol");
    for (const auto &ener : energies)
    {
        auto val = gmx::formatString("%.4f", ener.second);
        jtener.addValueUnit(interactionTypeToString(ener.first),
                            val.c_str(), unit);
    }
    jtree->addObject(jtener);
}

eMinimizeStatus minimizeMolecule(MsgHandler                    *msghandler,
                                 const ForceField              *pd,
                                 const ForceComputer           *forceComp,
                                 const SimulationConfigHandler &sch,
                                 const char                    *freezeFile,
                                 ACTMol                        *actmol,
                                 std::vector<gmx::RVec>        *xmin,
                                 JsonTree                      *jtree)
{
    std::vector<int> freeze;
    if (nullptr != freezeFile)
    {
        int isize;
        int *myindex;
        char *grpnames;
        rd_index(freezeFile, 1, &isize, &myindex, &grpnames);
        if (isize > 0)
        {
            printf("Will freeze %d atoms %s\n", isize, grpnames);
            for(int ii = 0; ii < isize; ii++)
            {
                freeze.push_back(myindex[ii]);
            }
        }
    }

    auto coords = *xmin;
    std::map<InteractionType, double> energies;
    {
        std::vector<gmx::RVec> forces(actmol->atomsConst().size());
        forceComp->compute(msghandler, pd, actmol->topology(), xmin, &forces, &energies);
        addEnergies(jtree, "Energies before", energies);
        // Multipoles
        msghandler->write("Multipole before minimization");
        analyse_multipoles(msghandler, actmol, {}, pd, forceComp, nullptr);
    }
    // Now do the minimization.
    msghandler->write(gmx::formatString("Starting minimization of '%s' using %s algorithm.\n",
                                        actmol->getMolname().c_str(),
                                        eMinimizeAlgorithmToString(sch.minAlg()).c_str()));

    MolHandler molhandler;
    auto eMin = molhandler.minimizeCoordinates(msghandler, pd, actmol->topology(), forceComp, sch,
                                               xmin, &energies, freeze);
    if (eMinimizeStatus::OK == eMin)
    {
        auto rmsd = molhandler.coordinateRmsd(actmol->atomsConst(), coords, xmin);
        msghandler->msg(ACTStatus::Info,
                        gmx::formatString("Final energy: %g RMSD wrt original structure %g nm.",
                                          energies[InteractionType::EPOT], rmsd));
        if (actmol->fragmentHandler()->topologies().size() == 2)
        {
            // Make a copy to keep the correct shell coordinates when saving the structure.
            auto newxmin = *xmin;
            std::map<InteractionType, double> einter;
            std::vector<gmx::RVec>            interactionForces;
            actmol->calculateInteractionEnergy(msghandler, pd, forceComp, &einter,
                                               &interactionForces, &newxmin, true);
            for(const auto &ei : einter)
            {
                msghandler->msg(ACTStatus::Info,
                                gmx::formatString("Interaction energy %s: %g",
                                                  interactionTypeToString(ei.first).c_str(), ei.second));
            }
        }
        // Multipoles
        msghandler->write("Multipole after minimization");
        // This is a hack for printing.
        for(auto &qp : (*actmol->qProps()))
        {
            qp.qPact()->setX(*xmin);
        }
        analyse_multipoles(msghandler, actmol, {}, pd, forceComp, nullptr);
        addEnergies(jtree, "Energies after", energies);
    }
    return eMin;
}

/*! \brief Tool to minimize the energy of a molecule
 * \param[in] argc Number of arguments on the cmd line
 * \param[in] argv The actual arguments
 * \return 0 if ok, 1 otherwise
 */
int minimize(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "alexandria minimize minimizes the energy of a molecule or a complex",
        "with respect to the coordinates, using a force field derived by the",
        "Alexandria Chemistry Toolkit. The input is given by a coordinate",
        "file, a force field file and command line options. The minimized",
        "structure and a log file are written.[PAR]",
        "A user selected group of atoms can be frozen.",
        "To do so, supply an index file with atom numbers (first atom is 1)",
        "and numbering should disregard shells and virtual sites if present",
        "in the model. Bond lengths, or bond lengths and angles, can be",
        "constrained with the -constrain option.[PAR]",
        "To obtain information on the used model, the final structure can be",
        "stored as a molprop file. The same minimization can be done with",
        "alexandria simulate -minimize, which can also run a simulation",
        "starting from the minimized structure."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff",      "aff",        ffREAD  },
        { efSTO, "-c",       "confout",    ffWRITE },
        { efNDX, "-freeze",  "freeze",     ffOPTRD },
        { efXML, "-omp",     "molprop",    ffOPTWR }
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
    sch.add_minimize_options(&pa, true);
    CompoundReader compR;
    compR.addOptions(&pa, &fnm, &desc);
    MsgHandler msghandler;
    msghandler.addOptions(&pa, &fnm, "minimize");
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
    sch.setMinimize(true);
    compR.optionsFinished(&msghandler, fnm);
    if (!msghandler.ok())
    {
        return 1;
    }
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
        readForceField(opt2fn("-ff", fnm.size(), fnm.data()), &pd, &msghandler);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    (void) pd.verifyCheckSum(&msghandler);
    for(const auto &ss : pd.info())
    {
        msghandler.write(ss);
    }

    ForceComputer forceComp;
    double epsilonr;
    if (!ffOption(pd, InteractionType::ELECTROSTATICS, "epsilonr", &epsilonr))
    {
        epsilonr = 1;
    }
    forceComp.init(shellToler, sch.maxIter(), forceComp.maxShellDistance(), epsilonr);

    std::vector<ACTMol> actmols;
    try
    {
        compR.read(&msghandler, pd, &forceComp, &actmols);
    }
    catch (const std::exception &ex)
    {
        // E.g. a file that does not exist
        msghandler.msg(ACTStatus::Error,
                       gmx::formatString("Could not read molecules: %s", ex.what()));
        actmols.clear();
    }
    if (actmols.empty())
    {
        fprintf(stderr, "Could not read or process molecules. Please check the log file for more information.\n");
        return 1;
    }
    auto &actmol = actmols[0];
    // Check whether we need to add position restraints
    if (sch.kFBPR() > 0)
    {
        actmol.topologyPtr()->addFBPR(&msghandler, sch.kFBPR(), sch.r0FBPR());
    }
    JsonTree jtree("minimize");
    if (msghandler.info())
    {
        forceFieldSummary(&jtree, &pd);
    }
    if (msghandler.debug())
    {
        actmol.topology()->dump(msghandler.twDebug());
    }
    msghandler.msg(ACTStatus::Info,
                   gmx::formatString("There are %lu fragments.",
                                     actmol.fragmentHandler()->topologies().size()));
    auto xmin = actmol.xOriginal();
    auto eMin = minimizeMolecule(&msghandler, &pd, &forceComp, sch,
                                 opt2fn_null("-freeze", fnm.size(), fnm.data()),
                                 &actmol, &xmin, &jtree);
    int status = 0;
    if (eMinimizeStatus::OK == eMin)
    {
        matrix box = { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } };
        actmol.PrintConformation(opt2fn("-c", fnm.size(), fnm.data()), xmin,
                                 sch.writeShells(), box);
        auto omp = opt2fn_null("-omp", fnm.size(),fnm.data());
        if (omp)
        {
            std::vector<MolProp> mpt;
            // We cannot use the actmol after this move command.
            mpt.push_back(std::move(*actmol.molProp()));
            MolPropWrite(omp, mpt, false);
        }
    }
    else
    {
        fprintf(stderr, "Minimization failed: %s, check log file %s\n",
                eMinimizeStatusToString(eMin).c_str(),
                msghandler.filename().c_str());
        status = 1;
    }
    if (json)
    {
        jtree.write("minimize.json", json);
    }
    else
    {
        int indent = 0;
        msghandler.write(jtree.writeString(json, &indent));
    }
    if (!msghandler.ok())
    {
        printf("Minimize failed with an error:\n");
        status = 1;
    }
    return status;
}

} // namespace alexandria
