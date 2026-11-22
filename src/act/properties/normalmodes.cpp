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

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <map>
#include <string>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/train_utility.h"
#include "act/basics/atomization_energy.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forces/forcecomputer.h"
#include "act/import/compound_reader.h"
#include "act/import/import.h"
#include "act/molprop/molprop_xml.h"
#include "act/properties/secondvirial.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_json.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/jsontree.h"
#include "act/utility/output_environment.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "normalmodes.h"

namespace alexandria
{

void NmaOptions::check() const
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

void NmaOptions::checkSettings() const
{
    if (!std::isfinite(lineWidth) || lineWidth <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -linewidth should be more than zero",
                                                           lineWidth)));
    }
    if (!std::isfinite(shellTolerance) || shellTolerance <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -shelltoler should be more than zero",
                                                           shellTolerance)));
    }
    simulation.validate();
}

int runNma(MsgHandler *msghandler, const NmaOptions &options, gmx_output_env_t *oenv)
{
    options.check();
    const auto &sch  = options.simulation;
    const bool  json = !options.jsonFile.empty();
    ForceField  pd;
    readForceField(options.forceFieldFile, &pd, msghandler);
    (void) pd.verifyCheckSum(msghandler);
    const double  shellToler = sch.shellTolerance(options.shellTolerance, msghandler);
    ForceComputer forceComp;
    forceComp.init(shellToler, sch.shellMaxIter(), pd.epsilonr());

    JsonTree jtree("nma");
    if (msghandler->verbose())
    {
        forceFieldSummary(&jtree, &pd);
    }

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
        return 1;
    }
    auto &actmol = actmols[0];

    if (msghandler->debug())
    {
        actmol.topology()->dump(msghandler->twDebug());
    }
    auto eMin = eMinimizeStatus::OK;
    const MolHandler molhandler;
    const std::vector<gmx::RVec> coords = actmol.xOriginal();
    std::vector<gmx::RVec> xmin   = coords;
    if (sch.minimize())
    {
        std::map<InteractionType, double> energies;
        eMin = molhandler.minimizeCoordinates(msghandler, &pd, actmol.topology(), &forceComp, sch,
                                              &xmin, &energies, {});
        // The status of the minimization for programs that read the output
        jtree.addObject("status", eMinimizeStatusToString(eMin));
        jtree.addBool("converged", eMinimizeStatus::OK == eMin);
        if (eMinimizeStatus::OK == eMin)
        {
            auto rmsd = molhandler.coordinateRmsd(actmol.atomsConst(), coords, &xmin);
            msghandler->msg(ACTStatus::Info,
                            gmx::formatString("Final energy: %g. RMSD wrt original structure %g nm.",
                                              energies[InteractionType::EPOT], rmsd));
            jtree.addValueUnit("rmsd", rmsd, "nm");
            JsonTree jtener("Energies");
            const std::string unit("kJ/mol");
            for (const auto &ener : energies)
            {
                jtener.addValueUnit(interactionTypeToString(ener.first),
                                    ener.second, unit);
            }
            jtree.addObject(jtener);

            if (!options.confout.empty())
            {
                matrix box = { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } };
                actmol.PrintConformation(options.confout.c_str(), xmin, false, box);
            }
        }
    }
    else
    {
        fprintf(stderr, "Running NMA without prior minimization, check your output!\n");
    }
    int status = 0;
    if (eMinimizeStatus::OK == eMin)
    {
        AtomizationEnergy atomenergy;
        atomenergy.read(msghandler);
        // A program without a command line has no output environment for
        // the spectrum
        const OutputEnvironment output(oenv);
        doFrequencyAnalysis(msghandler, &pd, &actmol, molhandler, &forceComp, &xmin,
                            atomenergy, nullptr, &jtree,
                            options.spectrumFile.empty() ? nullptr : options.spectrumFile.c_str(),
                            options.lineWidth, output.get(), msghandler->verbose(), json);
        if (json)
        {
            // The molecule and the structure the normal modes belong to
            jtree.addObject(moleculeToJson(&actmol));
            JsonTree jcoords("coordinates");
            jcoords.addObject("unit", "nm");
            jcoords.addObject(coordinatesToJson("minimized", xmin));
            jtree.addObject(jcoords);
        }
    }

    if (eMinimizeStatus::OK != eMin)
    {
        msghandler->msg(ACTStatus::Error, ACTMessage::MinimizationFailed,
                        gmx::formatString("check log file %s", msghandler->filename().c_str()));
        status = 1;
    }

    if (json)
    {
        jtree.write(options.jsonFile, json);
    }
    else
    {
        int indent = 0;
        msghandler->write(jtree.writeString(false, &indent));
    }
    return status;
}

int nma(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "alexandria nma performs a normal mode analysis, typically using",
        "a force field derived by the Alexandria Chemistry Toolkit.", 
        "The program will perform an energy minimization and a normal mode",
        "analysis including thermochemistry calculations.[PAR]",
        "The input is given by a coordinate file, a force field file,",
        "a molprop file containing charge information and",
        "command line options.[PAR]",
        "With [TT]-json[tt], the results are written in json format to the",
        "given file, nma.json by default: the energies, the frequencies and",
        "intensities, the thermochemistry, the particles and bonds of the",
        "molecule, the minimized coordinates and for each normal mode the",
        "reduced mass and the displacement of each particle, with which",
        "the vibrations can be shown."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff",      "aff",        ffREAD  },
        { efSTO, "-c",       "confout",    ffOPTWR },
        { efXVG, "-ir",      "IRspectrum", ffOPTWR },
        { efJSON, "-json",   "nma",        ffOPTWR }
    };
    gmx_output_env_t         *oenv;
    NmaOptions                options;
    real                      shellToler = options.shellTolerance;
    //! Line width (cm^-1) for a Lorentzian when computing infrared intensities and plotting an IR spectrum
    real                      linewidth  = options.lineWidth;
    std::vector<t_pargs>      pa = {
        { "-linewidth", FALSE, etREAL, {&linewidth},
          "Line width (cm^-1) for a Lorentzian when computing infrared intensities and plotting an IR spectrum" },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)" }
    };
    options.simulation.add_options(&pa, &fnm);
    MsgHandler msghandler;
    options.compounds.addOptions(&pa, &fnm, &desc);
    msghandler.addOptions(&pa, &fnm, "nma");
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
        options.lineWidth      = linewidth;
        options.forceFieldFile = opt2fn("-ff", fnm.size(), fnm.data());
        const auto optional    = [&fnm](const char *option)
        {
            const char *name = opt2fn_null(option, fnm.size(), fnm.data());
            return std::string(name ? name : "");
        };
        options.confout      = optional("-c");
        options.spectrumFile = optional("-ir");
        options.jsonFile     = optional("-json");
        return runNma(&msghandler, options, oenv);
    }
    catch (const gmx::InvalidInputError &ex)
    {
        fprintf(stderr, "%s\n", ex.what());
        return 1;
    }
}

} // namespace alexandria
