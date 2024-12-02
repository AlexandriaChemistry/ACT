/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022,2023,2024
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
        { efXML, "-charges", "charges",    ffOPTRD },
        { efSTO, "-c",       "confout",    ffOPTWR },
        { efLOG, "-g",       "nma",        ffWRITE },
        { efXVG, "-ir",      "IRspectrum", ffOPTWR }
    };
    gmx_output_env_t         *oenv;
    static char              *filename   = (char *)"";
    static char              *molnm      = (char *)"";
    static char              *qqm        = (char *)"";
    double                    qtot       = 0;
    double                    shellToler = 1e-6;
    bool                      verbose    = false;
    bool                      json       = false;
    //! Line width (cm^-1) for a Lorentzian when computing infrared intensities and plotting an IR spectrum
    double                    linewidth  = 24;
    std::vector<t_pargs>      pa = {
        { "-f",      FALSE, etSTR,  {&filename},
          "Molecular structure file in e.g. pdb format" },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Combined charge of the molecule(s). This will be taken from the input file by default, but that is not always reliable." },
        { "-qqm",    FALSE, etSTR,  {&qqm},
          "Use a method from quantum mechanics that needs to be present in the input file. Either ESP, Hirshfeld, CM5 or Mulliken may be available." },
        { "-linewidth", FALSE, etREAL, {&linewidth},
          "Line width (cm^-1) for a Lorentzian when computing infrared intensities and plotting an IR spectrum" },
        { "-v", FALSE, etBOOL, {&verbose},
          "Print more information to the log file." },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)" },
        { "-json", FALSE, etBOOL, {&json},
          "Print part of the output in json format" }
    };
    SimulationConfigHandler  sch;
    sch.setMinimize(true);
    sch.add_options(&pa, &fnm);
    CompoundReader compR;
    compR.addOptions(&pa, &fnm, &desc);
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 1;
    }
    sch.check_pargs();
    if (!compR.optionsOK(fnm))
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
    const char *logFileName = opt2fn("-g", fnm.size(),fnm.data());
    FILE *logFile   = gmx_ffopen(logFileName, "w");
    if (shellToler >= sch.forceTolerance())
    {
        shellToler = sch.forceTolerance()/10;
        printf("Shell tolerance larger than atom tolerance, changing it to %g\n", shellToler);
    }
    auto  forceComp = new ForceComputer(shellToler, 100);
    print_header(logFile, pa, fnm);
    compR.setLogfile(logFile);
    
    JsonTree jtree("simulate");
    if (verbose)
    {
        forceFieldSummary(&jtree, &pd);
    }

    std::vector<ACTMol> actmols = compR.read(pd, forceComp);
    if (actmols.empty())
    {
        return 1;
    }
    auto &actmol = actmols[0];
    std::vector<gmx::RVec> coords = actmol.xOriginal();

    if (debug)
    {
        actmol.topology()->dump(debug);
    }
    auto eMin = eMinimizeStatus::OK;
    /* Generate output file for debugging if requested */
    if (actmol.errors().empty())
    {
        MolHandler molhandler;
        std::vector<gmx::RVec> coords = actmol.xOriginal();
        std::vector<gmx::RVec> xmin   = coords;
        if (sch.minimize())
        {
            std::map<InteractionType, double> energies;
            eMin = molhandler.minimizeCoordinates(&pd, &actmol, forceComp, sch,
                                                  &xmin, &energies, logFile, {});
            if (eMinimizeStatus::OK == eMin)
            {
                auto rmsd = molhandler.coordinateRmsd(&actmol, coords, &xmin);
                fprintf(logFile, "Final energy: %g. RMSD wrt original structure %g nm.\n",
                        energies[InteractionType::EPOT], rmsd);
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
            AtomizationEnergy        atomenergy;
            doFrequencyAnalysis(&pd, &actmol, molhandler, forceComp, &xmin,
                                atomenergy, nullptr, &jtree,
                                opt2fn_null("-ir", fnm.size(), fnm.data()),
                                    linewidth, oenv, verbose);
        }
    }

    if (eMinimizeStatus::OK != eMin)
    {
        fprintf(stderr, "Minimization failed: %s, check log file %s\n",
                eMinimizeStatusToString(eMin).c_str(),
                logFileName);
        status = 1;
    }

    if (json)
    {
        jtree.write("nma.json", json);
    }
    else
    {
        jtree.fwrite(logFile, json);
    }
    gmx_ffclose(logFile);
    return status;
}

} // namespace alexandria
