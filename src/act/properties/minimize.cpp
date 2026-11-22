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
#include <cmath>
#include <cstdio>
#include <exception>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/train_utility.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/import/compound_reader.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "act/properties/secondvirial.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_json.h"
#include "act/topology/actmol_low.h"
#include "act/topology/fragmenthandler.h"
#include "act/utility/jsontree.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "minimize.h"

namespace alexandria
{

//! \brief Add energies to a json tree
static void addEnergies(JsonTree                                *jtree,
                        const char                              *name,
                        const std::map<InteractionType, double> &energies)
{
    JsonTree    jtener(name);
    const std::string unit("kJ/mol");
    for (const auto &ener : energies)
    {
        const auto val = gmx::formatString("%.4f", ener.second);
        jtener.addValueUnit(interactionTypeToString(ener.first),
                            val, unit);
    }
    jtree->addObject(jtener);
}

/*! \brief Move coordinates such that the center of mass is in the origin
 * \param[in] atoms  The particles with their masses
 * \param[in] coords The coordinates
 * \return the moved coordinates
 */
static std::vector<gmx::RVec> centeredCoordinates(const std::vector<ActAtom>   &atoms,
                                                  const std::vector<gmx::RVec> &coords)
{
    double    tmass = 0;
    gmx::RVec com   = { 0, 0, 0 };
    for (size_t i = 0; i < coords.size(); i++)
    {
        const double mass = atoms[i].mass();
        tmass += mass;
        for (int m = 0; m < DIM; m++)
        {
            com[m] += mass*coords[i][m];
        }
    }
    auto centered = coords;
    if (tmass > 0)
    {
        for (auto &x : centered)
        {
            for (int m = 0; m < DIM; m++)
            {
                x[m] -= com[m]/tmass;
            }
        }
    }
    return centered;
}

MinimizeResult minimizeStructure(MsgHandler                    *msghandler,
                                 const ForceField              *pd,
                                 const ForceComputer           *forceComp,
                                 const SimulationConfigHandler &sch,
                                 const std::vector<int>        &freeze,
                                 ACTMol                        *actmol,
                                 std::vector<gmx::RVec>        *xmin)
{
    // The minimizer counts the atoms only
    int natom = 0;
    for (const auto &atom : actmol->atomsConst())
    {
        natom += (atom.pType() == ActParticle::Atom) ? 1 : 0;
    }
    for (int index : freeze)
    {
        if (index < 0 || index >= natom)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("There is no atom %d to freeze, the molecule has %d atoms, numbered from zero",
                                                               index, natom)));
        }
    }
    MinimizeResult result;
    const auto coords = *xmin;
    std::map<InteractionType, double> energies;
    {
        std::vector<gmx::RVec> forces(actmol->atomsConst().size());
        forceComp->compute(msghandler, pd, actmol->topology(), xmin, &forces, &energies);
        result.energiesBefore = energies;
        // Multipoles
        msghandler->write("Multipole before minimization");
        analyse_multipoles(msghandler, actmol, {}, pd, forceComp, nullptr);
    }
    // Now do the minimization.
    msghandler->write(gmx::formatString("Starting minimization of '%s' using %s algorithm.\n",
                                        actmol->getMolname().c_str(),
                                        eMinimizeAlgorithmToString(sch.minAlg()).c_str()));

    const MolHandler molhandler;
    auto eMin = molhandler.minimizeCoordinates(msghandler, pd, actmol->topology(), forceComp, sch,
                                               xmin, &energies, freeze);
    result.status        = eMin;
    result.energiesAfter = energies;
    if (minimizationDidNotConverge(eMin))
    {
        msghandler->msg(ACTStatus::Warning,
                        gmx::formatString("Minimization did not converge: %s. The results are those of the structure with the lowest energy.",
                                          minimizeStatusMessage(eMin, energies[InteractionType::EPOT]).c_str()));
    }
    if (eMinimizeStatus::OK == eMin || minimizationDidNotConverge(eMin))
    {
        result.usable = true;
        auto rmsd = molhandler.coordinateRmsd(actmol->atomsConst(), coords, xmin);
        result.rmsd = rmsd;
        msghandler->msg(ACTStatus::Info,
                        gmx::formatString("Final energy: %g RMSD wrt original structure %g nm.",
                                          energies[InteractionType::EPOT], rmsd));
        if (actmol->fragmentHandler()->topologies().size() == 2)
        {
            // Make a copy to keep the correct shell coordinates when saving the structure.
            std::vector<std::vector<gmx::RVec>>            newxmin = { *xmin };
            std::vector<std::map<InteractionType, double>> einter;
            std::vector<double>                            epot;
            std::vector<std::vector<gmx::RVec>>            interactionForces;
            forceComp->computeInteraction(msghandler, pd, actmol, &newxmin, true,
                                          &einter, &epot, &interactionForces);
            for(const auto &ec : einter)
            {
                for(const auto &ei : ec)
                {
                    msghandler->msg(ACTStatus::Info,
                                    gmx::formatString("Interaction energy %s: %g",
                                                      interactionTypeToString(ei.first).c_str(), ei.second));
                }
            }
            if (!einter.empty())
            {
                result.interactionEnergies = einter[0];
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
    }
    return result;
}

eMinimizeStatus minimizeMolecule(MsgHandler                    *msghandler,
                                 const ForceField              *pd,
                                 const ForceComputer           *forceComp,
                                 const SimulationConfigHandler &sch,
                                 const char                    *freezeFile,
                                 ACTMol                        *actmol,
                                 std::vector<gmx::RVec>        *xmin,
                                 JsonTree                      *jtree,
                                 double                        *epot)
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
    const auto result = minimizeStructure(msghandler, pd, forceComp, sch, freeze, actmol, xmin);
    if (epot)
    {
        *epot = result.energiesAfter.count(InteractionType::EPOT) > 0 ?
            result.energiesAfter.at(InteractionType::EPOT) : 0.0;
    }
    addEnergies(jtree, "Energies before", result.energiesBefore);
    // The status for programs that read the output
    jtree->addObject("status", eMinimizeStatusToString(result.status));
    jtree->addBool("converged", eMinimizeStatus::OK == result.status);
    if (result.usable)
    {
        jtree->addValueUnit("rmsd", result.rmsd, "nm");
        // For programs that read the output, such as a scan of the distance
        if (!result.interactionEnergies.empty())
        {
            addEnergies(jtree, "Interaction energies", result.interactionEnergies);
        }
        addEnergies(jtree, "Energies after", result.energiesAfter);
    }
    return result.status;
}

void MinimizeOptions::check() const
{
    if (forceFieldFile.empty())
    {
        GMX_THROW(gmx::InvalidInputError("There is no force field file (-ff)"));
    }
    if (std::string(compounds.fileName()).empty() && std::string(compounds.smiles()).empty() &&
        std::string(compounds.dbName()).empty())
    {
        GMX_THROW(gmx::InvalidInputError("There is no molecule: give a structure file (-f), a SMILES string (-smiles) or a compound of the charge file (-db)"));
    }
    checkSettings();
}

void MinimizeOptions::checkSettings() const
{
    if (!std::isfinite(shellTolerance) || shellTolerance <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -shelltoler should be more than zero",
                                                           shellTolerance)));
    }
    simulation.validate();
}

int runMinimize(MsgHandler *msghandler, const MinimizeOptions &options)
{
    options.check();
    // The minimization is what this does, whatever the settings say
    auto sch = options.simulation;
    sch.setMinimize(true);
    const double shellToler = sch.shellTolerance(options.shellTolerance, msghandler);
    ForceField   pd;
    readForceField(options.forceFieldFile, &pd, msghandler);
    (void) pd.verifyCheckSum(msghandler);
    for(const auto &ss : pd.info())
    {
        msghandler->write(ss);
    }

    ForceComputer forceComp;
    forceComp.init(shellToler, sch.shellMaxIter(), pd.epsilonr());

    std::vector<ACTMol> actmols;
    // Reading changes the reader
    auto                compR = options.compounds;
    try
    {
        compR.read(msghandler, pd, &forceComp, &actmols);
    }
    catch (const std::exception &ex)
    {
        // E.g. a file that does not exist
        msghandler->msg(ACTStatus::Error,
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
        actmol.topologyPtr()->addFBPR(msghandler, sch.kFBPR(), sch.r0FBPR(), actmol.xOriginal());
    }
    JsonTree jtree("minimize");
    if (msghandler->info())
    {
        forceFieldSummary(&jtree, &pd);
    }
    if (msghandler->debug())
    {
        actmol.topology()->dump(msghandler->twDebug());
    }
    msghandler->msg(ACTStatus::Info,
                    gmx::formatString("There are %lu fragments.",
                                      actmol.fragmentHandler()->topologies().size()));
    auto xmin = actmol.xOriginal();
    // The description of the molecule has to be made before the molprop
    // is moved out of it below.
    const JsonTree jmolecule = moleculeToJson(&actmol);
    // Both structures are given with the center of mass in the origin, and
    // the final structure is superimposed on the initial one when the root
    // mean square deviation is computed, such that they can be shown together.
    JsonTree       jcoords("coordinates");
    jcoords.addObject("unit", "nm");
    jcoords.addObject("frame", "center of mass in the origin, after is superimposed on before");
    jcoords.addObject(coordinatesToJson("before", centeredCoordinates(actmol.atomsConst(), xmin)));
    double epot = 0;
    auto   eMin = minimizeMolecule(msghandler, &pd, &forceComp, sch,
                                   options.freezeFile.empty() ? nullptr : options.freezeFile.c_str(),
                                   &actmol, &xmin, &jtree, &epot);
    int status = 0;
    // Also without convergence the structure with the lowest energy is written
    if (eMinimizeStatus::OK == eMin || minimizationDidNotConverge(eMin))
    {
        if (!options.confout.empty())
        {
            matrix box = { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } };
            actmol.PrintConformation(options.confout.c_str(), xmin, sch.writeShells(), box);
        }
        if (!options.molpropFile.empty())
        {
            std::vector<MolProp> mpt;
            // We cannot use the actmol after this move command.
            mpt.push_back(std::move(*actmol.molProp()));
            MolPropWrite(options.molpropFile.c_str(), mpt, false);
        }
    }
    if (minimizationDidNotConverge(eMin))
    {
        // Without a lower minimum the structure is usable, see minimumUsable
        fprintf(stderr, "Minimization did not converge: %s, the structure with the lowest energy is written, check log file %s\n",
                minimizeStatusMessage(eMin, epot).c_str(),
                msghandler->filename().c_str());
        status = minimumUsable(eMin) ? 0 : 1;
    }
    else if (eMinimizeStatus::OK != eMin)
    {
        fprintf(stderr, "Minimization failed: %s, check log file %s\n",
                eMinimizeStatusToString(eMin).c_str(),
                msghandler->filename().c_str());
        status = 1;
    }
    const bool json = !options.jsonFile.empty();
    if (json)
    {
        jtree.addObject(jmolecule);
        if (eMinimizeStatus::OK == eMin || minimizationDidNotConverge(eMin))
        {
            jcoords.addObject(coordinatesToJson("after", centeredCoordinates(actmol.atomsConst(), xmin)));
        }
        jtree.addObject(jcoords);
        jtree.write(options.jsonFile, json);
    }
    else
    {
        int indent = 0;
        msghandler->write(jtree.writeString(json, &indent));
    }
    if (!msghandler->ok())
    {
        printf("Minimize failed with an error:\n");
        status = 1;
    }
    return status;
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
        "starting from the minimized structure.[PAR]",
        "With [TT]-json[tt], the results are written in json format to the",
        "given file, minimize.json by default: the status of the",
        "minimization, the energies before and after, the root mean square",
        "deviation between the two structures, the particles and bonds of",
        "the molecule and the coordinates of both structures."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff",      "aff",        ffREAD  },
        { efSTO, "-c",       "confout",    ffWRITE },
        { efNDX, "-freeze",  "freeze",     ffOPTRD },
        { efXML, "-omp",     "molprop",    ffOPTWR },
        { efJSON, "-json",   "minimize",   ffOPTWR }
    };
    gmx_output_env_t         *oenv;
    MinimizeOptions           options;
    real                      shellToler = options.shellTolerance;
    std::vector<t_pargs>      pa = {
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)" }
    };
    options.simulation.add_minimize_options(&pa, true);
    options.compounds.addOptions(&pa, &fnm, &desc);
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
    try
    {
        // What the command line asks for
        options.simulation.check_pargs(&msghandler);
        options.compounds.optionsFinished(&msghandler, fnm);
        if (!msghandler.ok())
        {
            return 1;
        }
        print_header(msghandler.tw(), pa, fnm);
        options.shellTolerance = shellToler;
        options.forceFieldFile = opt2fn("-ff", fnm.size(), fnm.data());
        options.confout        = opt2fn("-c", fnm.size(), fnm.data());
        const auto optional    = [&fnm](const char *option)
        {
            const char *name = opt2fn_null(option, fnm.size(), fnm.data());
            return std::string(name ? name : "");
        };
        options.freezeFile  = optional("-freeze");
        options.molpropFile = optional("-omp");
        options.jsonFile    = optional("-json");
        return runMinimize(&msghandler, options);
    }
    catch (const gmx::InvalidInputError &ex)
    {
        fprintf(stderr, "%s\n", ex.what());
        return 1;
    }
}

} // namespace alexandria
