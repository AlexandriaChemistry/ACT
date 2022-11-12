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

static void do_rerun(FILE          *logFile,
                     const Poldata *pd,
                     const MyMol   *mymol,
                     ForceComputer *forceComp,
                     const char    *trajname,
                     double         qtot)
{
    std::vector<MolProp> mps;
    std::string          method, basis;
    int                  maxpot = 100;
    int                  nsymm  = 1;
    const char          *molnm  = "";
    if (readBabel(trajname, &mps, molnm, molnm, "", &method,
                  &basis, maxpot, nsymm, "Opt", &qtot, false))
    {
        fprintf(logFile, "Doing energy calculation for %zu structures from %s\n",
                mps.size(), trajname);       
        std::map<InteractionType, double> energies;
        int mp_index = 0;
        for (auto mp : mps)
        {
            auto exper = mp.experimentConst();
            if (exper.size() == 1)
            {
                auto expx = exper[0].getCoordinates();
                std::vector<gmx::RVec> coords;
                const auto &atoms = mymol->atomsConst();
                if (expx.size() == atoms.size())
                {
                    // Assume there are shells in the input
                    coords = expx;
                }
                else
                {
                    size_t index = 0;
                    for(size_t i = 0; i < atoms.size(); i++)
                    {
                        if (index <= expx.size())
                        {
                            gmx::RVec xnm;
                            for(int m = 0; m < DIM; m++)
                            {
                                xnm[m] = expx[index][m];
                            }
                            coords.push_back(xnm);
                        }
                        else
                        {
                            GMX_THROW(gmx::InvalidInputError("Number of coordinates in trajectory does not match input file"));
                        }
                        if (atoms[i].pType() == eptAtom)
                        {
                            index++;
                        }
                    }
                }
                std::vector<gmx::RVec> forces(coords.size());
                forceComp->compute(pd, mymol->topology(),
                                   &coords, &forces, &energies);
                fprintf(logFile, "%5d", mp_index);
                for(const auto &ee : energies)
                {
                    fprintf(logFile, "  %s %8g", interactionTypeToString(ee.first).c_str(), ee.second);
                }
                fprintf(logFile, "\n");
            }
        }
    }
    else
    {
        fprintf(stderr, "Could not read compounds from %s\n", trajname);
    }
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
        { efLOG, "-g",  "simulation", ffWRITE },
        { efXVG, "-ir", "IRspectrum", ffOPTWR }
    };
    gmx_output_env_t         *oenv;
    static char              *filename   = (char *)"";
    static char              *trajname   = (char *)"";
    static char              *molnm      = (char *)"";
    double                    qtot       = 0;
    double                    shellToler = 1e-6;
    bool                      verbose    = false;
    std::vector<t_pargs>      pa = {
        { "-f",      FALSE, etSTR,  {&filename},
          "Molecular structure file in e.g. pdb format" },
        { "-traj",   FALSE, etSTR,  {&trajname},
          "Trajectory or series of structures of the same compound for which the energies will be computed. If this option is present, no simulation will be performed." },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Combined charge of the molecule(s). This will be taken from the input file by default, but that is not always reliable." },
        { "-v", FALSE, etBOOL, {&verbose},
          "Print more information to the log file." },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)" }
    };
    SimulationConfigHandler  sch;
    sch.add_pargs(&pa);
    
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 0;
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
    if (verbose)
    {
        forceFieldSummary(logFile, &pd);
    }

    MyMol                mymol;
    {
        std::vector<MolProp> mps;
        double               qtot_babel = qtot;
        std::string          method, basis;
        int                  maxpot = 100;
        int                  nsymm  = 1;
        if (readBabel(filename, &mps, molnm, molnm, "", &method,
                      &basis, maxpot, nsymm, "Opt", &qtot_babel,
                      false))
        {
            if (mps.size() > 1)
            {
                fprintf(stderr, "Warning: will only use the first compound (out of %zu) in %s\n", mps.size(), filename);
            }
            std::map<std::string, std::string> g2a;
            gaffToAlexandria("", &g2a);
            bool mappingOK = true;
            if (!g2a.empty())
            {
                mappingOK = renameAtomTypes(&mps[0], g2a);
            }
            if (mappingOK)
            {
                mymol.Merge(&mps[0]);
            }
        }
        else
        {
            fprintf(logFile, "Reading %s failed.\n", filename);
        }
    }
    auto imm = mymol.GenerateTopology(logFile, &pd, missingParameters::Error,
                                      false);
    if (immStatus::OK != imm)
    {
        fprintf(stderr, "Problem generating topology: %s\n", immsg(imm));
        return 0;
    }
    CommunicationRecord cr;
    gmx::MDLogger  mdlog {};
    std::vector<gmx::RVec> forces(mymol.atomsConst().size());
    std::vector<gmx::RVec> coords = mymol.xOriginal();
    if (immStatus::OK == imm)
    {
        std::vector<double> myq;
        auto alg = pd.chargeGenerationAlgorithm();
        imm    = mymol.GenerateCharges(&pd, forceComp, mdlog, &cr, alg, myq, &coords, &forces);
    }
    if (verbose && pd.polarizable())
    {
        // Make a copy since it maybe changed
        auto xx    = coords;
        auto qCalc = mymol.qTypeProps(qType::Calc);
        forceComp->calcPolarizability(&pd, mymol.topology(), &xx, qCalc);
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
    if (strlen(trajname) > 0)
    {
        do_rerun(logFile, &pd, &mymol, forceComp, trajname, qtot);
    }
    else if (immStatus::OK == imm && mymol.errors().empty())
    {
        MolHandler molhandler;
        std::vector<gmx::RVec> coords = mymol.xOriginal();
        std::vector<gmx::RVec> xmin   = coords;
        auto eMin = eMinimizeStatus::OK;
        if (sch.nma() || sch.minimize())
        {
            std::map<InteractionType, double> energies;
            eMin = molhandler.minimizeCoordinates(&pd, &mymol, forceComp, sch,
                                                  &xmin, &energies, logFile);
            if (eMinimizeStatus::OK == eMin)
            {
                auto rmsd = molhandler.coordinateRmsd(&mymol, coords, &xmin);
                fprintf(logFile, "Final energy: %g. RMSD wrt original structure %g nm.\n",
                        energies[InteractionType::EPOT], rmsd);
                matrix box;
                clear_mat(box);
                write_sto_conf(opt2fn("-c", fnm.size(),fnm.data()), 
                               mymol.getMolname().c_str(),
                               mymol.gmxAtomsConst(),
                               as_rvec_array(xmin.data()), nullptr,
                               epbcNONE, box);
            }
            else
            {
                fprintf(stderr, "Minimization failed: %s, check log file %s\n",
                        eMinimizeStatusToString(eMin).c_str(),
                        logFileName);
            }
        }
        if (immStatus::OK == imm && eMinimizeStatus::OK == eMin)
        {
            if (sch.nma())
            {
                if (eMinimizeStatus::OK == eMin)
                {
                    AtomizationEnergy        atomenergy;
                    std::vector<std::string> output;
                    doFrequencyAnalysis(&pd, &mymol, molhandler, forceComp, &coords,
                                        atomenergy, nullptr, &output,
                                        opt2fn_null("-ir", fnm.size(), fnm.data()),
                                        sch.lineWidth(), oenv,
                                        sch.lapack(), verbose);
                    for(const auto &op : output)
                    {
                        fprintf(logFile, "%s\n", op.c_str());
                    }
                }
                else
                {
                    fprintf(logFile, "Cannot do NMA because energy minimization failed to converge.\n");
                }
            }
            else
            {
                molhandler.simulate(&pd, &mymol, forceComp, sch, logFile,
                                    opt2fn("-o", fnm.size(),fnm.data()),
                                    opt2fn("-e", fnm.size(),fnm.data()),
                                    oenv);
            }
        }
    }
    else
    {
        fprintf(stderr, "\nFatal Error. Please check the log file %s for error messages.\n", logFileName);
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
