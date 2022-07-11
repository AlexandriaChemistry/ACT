/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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

#include <ctype.h>
#include <stdlib.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/utility/futil.h"

#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/poldata/poldata_xml.h"
#include "alexandria/alex_modules.h"
#include "alexandria/atype_mapping.h"
#include "alexandria/babel_io.h"
#include "alexandria/confighandler.h"
#include "alexandria/molhandler.h"
#include "alexandria/mymol.h"
#include "alexandria/tuning_utility.h"

namespace alexandria
{

int simulate(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "alexandria simulate performs a proof-of-principle MD simulation, typically using",
        "a force field derived by the Alexandria Chemistry Toolkit.", 
        "The features are limited to",
        "constant energy simulations in vacuum on one core.",
        "The input is given by a coordinate file, a force field file and",
        "command line options. During the simulation an energy file,",
        "a trajectory file and a log file are generated."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff", "gentop",     ffREAD  },
        { efPDB, "-o",  "trajectory", ffWRITE },
        { efSTO, "-c",  "confout",    ffOPTWR },
        { efXVG, "-e",  "energy",     ffWRITE },
        { efLOG, "-g",  "simulation", ffWRITE }
    };
    gmx_output_env_t         *oenv;
    static char              *filename = (char *)"";
    static char              *molnm    = (char *)"";
    double                    qtot     = 0;
    std::vector<t_pargs>      pa = {
        { "-f",      FALSE, etSTR,  {&filename},
          "Molecular structure file in e.g. pdb format" },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Combined charge of the molecule(s). This will be taken from the input file by default, but that is not always reliable." },
    };
    SimulationConfigHandler  sch;
    sch.add_pargs(&pa);
    
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 0;
    }

    Poldata        pd;
    try
    {
        readPoldata(opt2fn("-ff", fnm.size(), fnm.data()), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    (void) pd.verifyCheckSum(stderr);
    auto  forceComp = new ForceComputer(&pd);
    FILE *logFile   = gmx_ffopen(opt2fn("-g", fnm.size(),fnm.data()), "w");
    print_header(logFile, pa);
    
    MyMol                mymol;
    {
        MolProp     mp;
        double      qtot_babel = qtot;
        std::string method, basis;
        int         maxpot = 100;
        int         nsymm  = 1;
        if (readBabel(filename, &mp, molnm, molnm, "", &method,
                      &basis, maxpot, nsymm, "Opt", &qtot_babel,
                      false))
        {
            std::map<std::string, std::string> g2a;
            gaffToAlexandria("", &g2a);
            bool mappingOK = true;
            if (!g2a.empty())
            {
                mappingOK = renameAtomTypes(&mp, g2a);
            }
            if (mappingOK)
            {
                mymol.Merge(&mp);
            }
        }
        else
        {
            fprintf(logFile, "Reading %s failed.\n", filename);
        }
    }
    auto imm = mymol.GenerateTopology(logFile, &pd, missingParameters::Error,
                                      false);

    CommunicationRecord cr;
    gmx::MDLogger  mdlog {};
    if (immStatus::OK == imm)
    {
        std::vector<double> myq;
        auto alg = pd.chargeGenerationAlgorithm();
        imm    = mymol.GenerateCharges(&pd, forceComp, mdlog, &cr, alg, myq);
    }
    /* Generate output file for debugging if requested */
    if (immStatus::OK == imm && mymol.errors().size() == 0)
    {
        MolHandler molhandler;
        if (sch.nma() || sch.minimize())
        {
            int myIter = molhandler.minimizeCoordinates(&mymol, forceComp, logFile, 10000);
            fprintf(logFile, "Number of iterations %d, final energy %g\n",
                    myIter, mymol.potentialEnergy());
            matrix box;
            clear_mat(box);
            std::vector<gmx::RVec> xx;
            for(const auto &x1 : mymol.coordinateSet(coordSet::Minimized))
            {
                xx.push_back(x1);
            }
            write_sto_conf(opt2fn("-c", fnm.size(),fnm.data()), 
                           mymol.getMolname().c_str(),
                           mymol.gmxAtomsConst(),
                           as_rvec_array(xx.data()), nullptr,
                           epbcNONE, box);
        }
        if (immStatus::OK == imm)
        {
            if (sch.nma())
            {
                std::vector<std::string> output;
                doFrequencyAnalysis(&mymol, molhandler, forceComp, nullptr, &output);
                for(const auto &op : output)
                {
                    fprintf(logFile, "%s\n", op.c_str());
                }
            }
            else
            {
                molhandler.simulate(&mymol, forceComp, sch, logFile,
                                    opt2fn("-o", fnm.size(),fnm.data()),
                                    opt2fn("-e", fnm.size(),fnm.data()),
                                    oenv);
            }
        }
    }
    else
    {
        fprintf(stderr, "\nFatal Error. Please check the log file for error messages.\n");
        fprintf(logFile, "%s\n", immsg(imm));
        for(const auto &err: mymol.errors())
        {
            fprintf(logFile, "%s\n", err.c_str());
        }
    }
    gmx_ffclose(logFile);
    return 0;
}

} // namespace alexandria
