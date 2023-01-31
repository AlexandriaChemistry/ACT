/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022,2023
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
#include "act/alexandria/atype_mapping.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/mymol.h"
#include "act/alexandria/secondvirial.h"
#include "act/alexandria/tuning_utility.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/poldata/poldata_xml.h"
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
        "The input is given by a coordinate file, a force field file and",
        "command line options."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff", "aff",        ffREAD  },
        { efSTO, "-c",  "confout",    ffOPTWR },
        { efLOG, "-g",  "simulation", ffWRITE },
        { efXVG, "-ir", "IRspectrum", ffOPTWR }
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
    sch.add_pargs(&pa);
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        status = 1;
        return status;
    }
    sch.check_pargs();

    Poldata        pd;
    try
    {
        readPoldata(opt2fn("-ff", fnm.size(), fnm.data()), &pd);
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
    print_header(logFile, pa);
    
    JsonTree jtree("simulate");
    if (verbose)
    {
        forceFieldSummary(&jtree, &pd);
    }

    MyMol                mymol;
    {
        std::vector<MolProp> mps;
        double               qtot_babel = qtot;
        std::string          method, basis;
        int                  maxpot = 100;
        int                  nsymm  = 1;
        if (!readBabel(filename, &mps, molnm, molnm, "", &method,
                       &basis, maxpot, nsymm, "Opt", &qtot_babel,
                       false))
        {
            fprintf(logFile, "Reading %s failed.\n", filename);
            status = 1;
        }
        else
        {
            std::map<std::string, std::string> g2a;
            gaffToAlexandria("", &g2a);
            if (!g2a.empty())
            {
                if (!renameAtomTypes(&mps[0], g2a))
                {
                    status = 1;
                }
            }
        }
            
        if (status == 0)
        {
            if (mps.size() > 1)
            {
                fprintf(stderr, "Warning: will only use the first compound (out of %zu) in %s\n", mps.size(), filename);
            }
            mymol.Merge(&mps[0]);
        }
    }
    immStatus imm = immStatus::OK;
    if (status == 0)
    {
        imm = mymol.GenerateTopology(logFile, &pd, missingParameters::Error, false);
    }
    std::vector<gmx::RVec> coords = mymol.xOriginal();
    if (immStatus::OK == imm && status == 0)
    {
        CommunicationRecord cr;
        gmx::MDLogger  mdlog {};
        std::vector<gmx::RVec> forces(mymol.atomsConst().size());

        std::vector<double> myq;
        auto alg   = pd.chargeGenerationAlgorithm();
        auto qtype = qType::Calc;
        if (strlen(qqm) > 0)
        {
            alg   = ChargeGenerationAlgorithm::Read;
            qtype = stringToQtype(qqm);
        }
        imm    = mymol.GenerateCharges(&pd, forceComp, mdlog, &cr, alg, qtype, myq, &coords, &forces);
    }
    if (immStatus::OK == imm && status == 0)
    {
        if (debug)
        {
            mymol.topology()->dump(debug);
        }
        auto eMin = eMinimizeStatus::OK;
        /* Generate output file for debugging if requested */
        if (mymol.errors().empty())
        {
            MolHandler molhandler;
            std::vector<gmx::RVec> coords = mymol.xOriginal();
            std::vector<gmx::RVec> xmin   = coords;
            if (sch.minimize())
            {
                std::map<InteractionType, double> energies;
                eMin = molhandler.minimizeCoordinates(&pd, &mymol, forceComp, sch,
                                                      &xmin, &energies, logFile);
                if (eMinimizeStatus::OK == eMin)
                {
                    auto rmsd = molhandler.coordinateRmsd(&mymol, coords, &xmin);
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
                        write_sto_conf(confout, mymol.getMolname().c_str(),
                                       mymol.gmxAtomsConst(),
                                       as_rvec_array(xmin.data()), nullptr,
                                       epbcNONE, box);
                    }
                    
                    AtomizationEnergy        atomenergy;
                    doFrequencyAnalysis(&pd, &mymol, molhandler, forceComp, &coords,
                                        atomenergy, nullptr, &jtree,
                                        opt2fn_null("-ir", fnm.size(), fnm.data()),
                                        linewidth, oenv, sch.lapack(), verbose);
                }
            }
        }
        
        if (eMinimizeStatus::OK != eMin)
        {
            fprintf(stderr, "Minimization failed: %s, check log file %s\n",
                    eMinimizeStatusToString(eMin).c_str(),
                    logFileName);
            status = 1;
        }
        else if (immStatus::OK != imm)
        {
            fprintf(stderr, "\nFatal Error. Please check the log file %s for error messages.\n", logFileName);
            fprintf(logFile, "%s\n", immsg(imm));
            for(const auto &err: mymol.errors())
            {
                fprintf(logFile, "%s\n", err.c_str());
            }
            status = 1;
        }
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
