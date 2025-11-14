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

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/compound_reader.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/actmol.h"
#include "act/alexandria/secondvirial.h"
#include "act/alexandria/train_utility.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/jsontree.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

int nma(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "alexandria nma performs a normal mode analysis, typically using",
        "a force field derived by the Alexandria Chemistry Toolkit.", 
        "The program will perform an energy minimization and a normal mode",
        "analysis including thermochemistry calculations.[PAR]",
        "The input is given by a coordinate file, a force field file,",
        "a molprop file containing charge information and",
        "command line options."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff",      "aff",        ffREAD  },
        { efSTO, "-c",       "confout",    ffOPTWR },
        { efXVG, "-ir",      "IRspectrum", ffOPTWR }
    };
    gmx_output_env_t         *oenv;
    double                    shellToler = 1e-6;
    bool                      json       = false;
    //! Line width (cm^-1) for a Lorentzian when computing infrared intensities and plotting an IR spectrum
    double                    linewidth  = 24;
    std::vector<t_pargs>      pa = {
        { "-linewidth", FALSE, etREAL, {&linewidth},
          "Line width (cm^-1) for a Lorentzian when computing infrared intensities and plotting an IR spectrum" },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)" },
        { "-json", FALSE, etBOOL, {&json},
          "Print part of the output in json format" }
    };
    SimulationConfigHandler  sch;
    sch.setMinimize(true);
    sch.add_options(&pa, &fnm);
    MsgHandler msghandler;
    CompoundReader compR;
    compR.addOptions(&pa, &fnm, &desc);
    msghandler.addOptions(&pa, &fnm, "nma");
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
    compR.optionsFinished(&msghandler, fnm);
    if (!msghandler.ok())
    {
        return 1;
    }

    ForceField        pd;
    try
    {
        readForceField(opt2fn("-ff", fnm.size(), fnm.data()), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    (void) pd.verifyCheckSum(stderr);
    if (shellToler >= sch.forceTolerance())
    {
        shellToler = sch.forceTolerance()/10;
        msghandler.msg(ACTStatus::Warning,
                       gmx::formatString("Shell tolerance larger than atom tolerance, changing it to %g",
                                         shellToler));
    }
    ForceComputer forceComp;
    forceComp.init(shellToler, 100);
    print_header(msghandler.tw(), pa, fnm);
    
    JsonTree jtree("simulate");
    if (msghandler.verbose())
    {
        forceFieldSummary(&jtree, &pd);
    }

    std::vector<ACTMol> actmols = compR.read(&msghandler, pd, &forceComp);
    if (actmols.empty())
    {
        return 1;
    }
    auto &actmol = actmols[0];

    if (debug)
    {
        actmol.topology()->dump(debug);
    }
    auto eMin = eMinimizeStatus::OK;
    /* Generate output file for debugging if requested */
    MolHandler molhandler;
    std::vector<gmx::RVec> coords = actmol.xOriginal();
    std::vector<gmx::RVec> xmin   = coords;
    if (sch.minimize())
    {
        std::map<InteractionType, double> energies;
        eMin = molhandler.minimizeCoordinates(&msghandler, &pd, &actmol, &forceComp, sch,
                                              &xmin, &energies, {});
        if (eMinimizeStatus::OK == eMin)
        {
            auto rmsd = molhandler.coordinateRmsd(&actmol, coords, &xmin);
            msghandler.msg(ACTStatus::Info,
                           gmx::formatString("Final energy: %g. RMSD wrt original structure %g nm.",
                                             energies[InteractionType::EPOT], rmsd));
            JsonTree jtener("Energies");
            std::string unit("kJ/mol");
            for (const auto &ener : energies)
            {
                jtener.addValueUnit(interactionTypeToString(ener.first),
                                    gmx_ftoa(ener.second), unit);
            }
            jtree.addObject(jtener);

            auto confout = opt2fn_null("-c", fnm.size(),fnm.data());
            if (confout)
            {
                matrix box;
                clear_mat(box);
                // TODO This will crash
                write_sto_conf(confout, actmol.getMolname().c_str(),
                               nullptr,
                               as_rvec_array(xmin.data()), nullptr,
                               epbcNONE, box);
            }
        }
    }
    else
    {
        fprintf(stderr, "Running NMA with prior minimization, check your output!\n");
    }
    if (eMinimizeStatus::OK == eMin)
    {
        AtomizationEnergy atomenergy;
        atomenergy.read(&msghandler);
        doFrequencyAnalysis(&pd, &actmol, molhandler, &forceComp, &xmin,
                            atomenergy, nullptr, &jtree,
                            opt2fn_null("-ir", fnm.size(), fnm.data()),
                            linewidth, oenv, msghandler.verbose());
    }

    if (eMinimizeStatus::OK != eMin)
    {
        msghandler.msg(ACTStatus::Error, ACTMessage::MinimizationFailed,
                       gmx::formatString("check log file %s", msghandler.filename().c_str()));
        status = 1;
    }

    if (json)
    {
        jtree.write("nma.json", json);
    }
    else
    {
        int indent = 0;
        msghandler.write(jtree.writeString(false, &indent));
    }
    return status;
}

} // namespace alexandria
