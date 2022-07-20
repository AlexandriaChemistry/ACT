/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

static void forceFieldSummary(FILE          *fp,
                              const Poldata *pd)
{
    fprintf(fp, "\nForce field summary.\n--------------------\n");
    fprintf(fp, "Input file:        %s\n", pd->filename().c_str());
    fprintf(fp, "Created:           %s\n", pd->timeStamp().c_str());
    fprintf(fp, "Checksum:          %s\n", pd->checkSum().c_str());
    fprintf(fp, "Polarizable:       %s\n", yesno_names[pd->polarizable()]);
    fprintf(fp, "Charge generation: %s\n", 
            chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str());
    fprintf(fp, "Using %d exclusions\n", pd->getNexcl());
    fprintf(fp, "Relative dielectric constant epsilon_r %g\n", pd->getEpsilonR());
    fprintf(fp, "There are %zu particle types\n", pd->getNatypes());
    for(const auto &fs : pd->forcesConst())
    {
        auto itype = fs.first;
        auto &ffpl = fs.second;
        if (!ffpl.parametersConst().empty())
        {
            
            if (!ffpl.function().empty())
            {
                fprintf(fp, "InteractionType %s force function %s has %zu entries.\n",
                        interactionTypeToString(itype).c_str(),
                        ffpl.function().c_str(),
                        ffpl.parametersConst().size());
            }
            else
            {
                fprintf(fp, "There are %5zu entries for %s.\n",
                        ffpl.parametersConst().size(),
                        interactionTypeToString(itype).c_str());
            }
        }
    }
    fprintf(fp, "\n");
}

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
    static char              *filename   = (char *)"";
    static char              *molnm      = (char *)"";
    double                    qtot       = 0;
    double                    forceToler = 0;
    double                    overRelax  = 1.0;
    int                       maxIter    = 0;
    bool                      verbose    = false;
    std::vector<t_pargs>      pa = {
        { "-f",      FALSE, etSTR,  {&filename},
          "Molecular structure file in e.g. pdb format" },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Combined charge of the molecule(s). This will be taken from the input file by default, but that is not always reliable." },
        { "-toler",  FALSE, etREAL, {&forceToler},
          "Convergence tolerance on the mean square force for the energy minimizer. If zero it will be determined automatically." },
        { "-maxiter",FALSE, etINT,  {&maxIter},
          "Maximum number of iterations for the energy minimizer, 0 is until convergence." },
        { "-overrelax", FALSE, etREAL, {&overRelax},
          "Apply overrelaxation (if > 1) to speed up minimization. Can be dangerous for poor energy functions." },
        { "-v", FALSE, etBOOL, {&verbose},
          "Print more information to the log file." },
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
    if (verbose)
    {
        forceFieldSummary(logFile, &pd);
    }

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
    std::vector<gmx::RVec> forces(mymol.atomsConst().size());
    if (immStatus::OK == imm)
    {
        std::vector<double> myq;
        auto alg = pd.chargeGenerationAlgorithm();
        imm    = mymol.GenerateCharges(&pd, forceComp, mdlog, &cr, alg, myq, &forces);
    }
    if (verbose && pd.polarizable())
    {
        std::vector<gmx::RVec> xx(mymol.atomsConst().size());
        for(size_t i = 0; i < xx.size(); i++)
        {
            copy_rvec(mymol.x()[i], xx[i]);
        }
        auto qCalc = mymol.qTypeProps(qType::Calc);
        forceComp->calcPolarizability(mymol.topology(), &xx, qCalc);
        auto alpha = qCalc->polarizabilityTensor();
        double fac = convertFromGromacs(1, "A^3");
        fprintf(logFile, "Alpha trace: %10g %10g %10g. Isotropic: %10g\n",
                fac*alpha[XX][XX], fac*alpha[YY][YY], fac*alpha[ZZ][ZZ], 
                fac*qCalc->isotropicPolarizability());
    }

    if (debug)
    {
        mymol.topology()->dump(debug);
    }
    /* Generate output file for debugging if requested */
    if (immStatus::OK == imm && mymol.errors().size() == 0)
    {
        MolHandler molhandler;
        if (verbose && pd.polarizable())
        {
            
        }
        if (sch.nma() || sch.minimize())
        {
            int myIter = molhandler.minimizeCoordinates(&mymol, forceComp, logFile, maxIter, overRelax, forceToler);
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
                AtomizationEnergy        atomenergy;
                std::vector<std::string> output;
                doFrequencyAnalysis(&mymol, molhandler, forceComp, atomenergy, nullptr, &output);
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
