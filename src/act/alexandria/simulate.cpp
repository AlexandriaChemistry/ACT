/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024
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
#include "act/alexandria/babel_io.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/fetch_charges.h"
#include "act/alexandria/fragmenthandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/secondvirial.h"
#include "act/alexandria/train_utility.h"
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
        "alexandria simulate performs a proof-of-principle MD simulation, typically using",
        "a force field derived by the Alexandria Chemistry Toolkit.", 
        "The program can perform energy minimization and/or do a",
        "constant energy molecular dynamics simulation in vacuum.",
        "In addition, a series of conformations (trajectory) may be",
        "submitted after which the energy per conformation is printed.[PAR]",
        "The input is given by a coordinate file, a force field file and",
        "command line options. During the simulation an energy file,",
        "a trajectory file and a log file are generated.[PAR]",
        "It is highly recommended to provide the [TT]-charges[tt] option",
        "with a molprop file that will be used to generate charges for the",
        "file specified with the [TT]-f[tt] option.",
        "During minimization a user select group of atoms can be frozen.",
        "To do so, supply an index file with atom numbers (first atom is 1)",
        "and numbering should disregard shells if present."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff",      "aff",        ffREAD  },
        { efXML, "-charges", "molprop",    ffOPTRD },
        { efPDB, "-o",       "trajectory", ffWRITE },
        { efSTO, "-c",       "confout",    ffOPTWR },
        { efXVG, "-e",       "energy",     ffWRITE },
        { efLOG, "-g",       "simulation", ffWRITE },
        { efNDX, "-freeze",  "freeze",     ffOPTRD }
    };
    gmx_output_env_t         *oenv;
    static char              *filename   = (char *)"";
    static char              *molnm      = (char *)"";
    static char              *qqm        = (char *)"";
    double                    qtot       = 0;
    double                    shellToler = 1e-6;
    bool                      verbose    = false;
    bool                      json       = false;
    static char              *qcustom    = (char *)"";
    std::vector<t_pargs>      pa = {
        { "-f",      FALSE, etSTR,  {&filename},
          "Molecular structure file in e.g. pdb format" },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Combined charge of the molecule(s). This will be taken from the input file by default, but that is not always reliable. If the -qcustom flag is supplied, that will be used instead." },
        { "-qqm",    FALSE, etSTR,  {&qqm},
          "Use a method from quantum mechanics that needs to be present in the input file. Either ESP, Hirshfeld, CM5 or Mulliken may be available." },
        { "-qcustom", FALSE, etSTR, {&qcustom}, 
          "Here a quoted string of custom charges can be provided such that a third party source can be used. It is then possible to generate multipoles and compare the ESP to a quantum chemistry result. The number of charges provided must match the number of particles (including shells if present in the force field used)." },
        { "-v", FALSE, etBOOL, {&verbose},
          "Print more information to the log file." },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)" },
        { "-json", FALSE, etBOOL, {&json},
          "Print part of the output in json format" }
    };
    SimulationConfigHandler  sch;
    sch.add_options(&pa, &fnm);
    sch.add_MD_options(&pa);
    DimerGenerator           gendimers;
    // We do not want to see those options in simulate, just in b2.
    // gendimers.addOptions(&pa);
    ReRunner                 rerun(false);
    rerun.addOptions(&pa, &fnm);
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        status = 1;
        return status;
    }
    sch.check_pargs();

    if (strlen(filename) == 0)
    {
        fprintf(stderr, "Please provide a filename with a structure to optimize.\n");
        return 0;
    }
    const char *logFileName = opt2fn("-g", fnm.size(),fnm.data());
    FILE *logFile   = gmx_ffopen(logFileName, "w");
    print_header(logFile, pa, fnm);
    if (shellToler >= sch.forceTolerance())
    {
        shellToler = sch.forceTolerance()/10;
        fprintf(logFile, "\nShell tolerance larger than atom tolerance, changing it to %g\n", shellToler);
    }
    ForceField        pd;
    try
    {
        readForceField(opt2fn("-ff", fnm.size(), fnm.data()), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    (void) pd.verifyCheckSum(stderr);

    chargeMap qmap;
    auto forceComp = new ForceComputer(shellToler, sch.maxIter());
    bool userqtot  = false;
    auto qfn       = opt2fn_null("-charges", fnm.size(), fnm.data());
    if (qfn)
    {
        qmap = fetchChargeMap(&pd, forceComp, qfn);
        fprintf(logFile, "\nRead %lu entries into charge map from %s\n", qmap.size(), qfn);
        if (strlen(qqm) > 0)
        {
            fprintf(stderr, "WARNING: Since you provide a charge map, I will ignore the -qqm flag.\n");
        }
        if (strlen(qcustom) > 0)
        {
            fprintf(stderr, "WARNING: Since you provide a charge map, I will ignore the -qcustom flag.\n");
        }
    }
    else if (strlen(qcustom) > 0)
    {
        std::string descr("qcustom");
        double myq = 0;
        for(const auto &word : split(qcustom, ' '))
        {
            myq += my_atof(word, descr);
        }
        qtot     = static_cast<int>(std::round(myq));
        userqtot = true;
    }

    JsonTree jtree("simulate");
    if (verbose)
    {
        forceFieldSummary(&jtree, &pd);
    }

    ACTMol actmol;
    matrix box;
    clear_mat(box);
    {
        std::vector<MolProp> mps;
        double               qtot_babel = qtot;
        std::string          method, basis;
        int                  maxpot = 100;
        int                  nsymm  = 1;
        if (!readBabel(&pd, filename, &mps, molnm, molnm, "", &method,
                       &basis, maxpot, nsymm, "Opt", userqtot, &qtot_babel,
                       false, box, sch.oneH()))
        {
            fprintf(logFile, "Reading %s failed.\n", filename);
            status = 1;
        }
        if (status == 0)
        {
            if (mps.size() > 1)
            {
                fprintf(stderr, "Warning: will only use the first compound (out of %zu) in %s\n", mps.size(), filename);
            }
            if (mps.size() == 0)
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Failed to import coordinate file %s using OpenBabel", filename).c_str()));
            }
            actmol.Merge(&mps[0]);
        }
    }
    if (actmol.totalCharge() != qtot)
    {
        fprintf(logFile, "WARNING: detected total charge %d, command line says %g.\n",
                actmol.totalCharge(), qtot);
    }
    
    immStatus imm = immStatus::OK;
    if (status == 0)
    {
        imm = actmol.GenerateTopology(logFile, &pd, missingParameters::Error);
    }
    std::vector<gmx::RVec> coords = actmol.xOriginal();
    if (immStatus::OK == imm && status == 0)
    {
        auto fragments  = actmol.fragmentHandler();
        if (!qmap.empty())
        {
            if (fragments->setCharges(qmap))
            {
                // Copy charges to the high-level topology as well
                fragments->fetchCharges(actmol.atoms());
            }
            else
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Not all compounds present in the charge map %s", qfn).c_str()));
            }
        }
        else
        {
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());

            std::vector<double> myq;
            auto alg   = pd.chargeGenerationAlgorithm();
            auto qtype = qType::Calc;
            if (strlen(qcustom) > 0)
            {
                auto mycharges = gmx::splitString(qcustom);
                for(auto &q : mycharges)
                {
                    myq.push_back(my_atof(q.c_str(), "custom q"));
                }
                alg = ChargeGenerationAlgorithm::Custom;
            }
            else if (strlen(qqm) > 0)
            {
                alg   = ChargeGenerationAlgorithm::Read;
                qtype = stringToQtype(qqm);
            }
            imm    = actmol.GenerateCharges(&pd, forceComp, alg, qtype, myq, &coords, &forces);
        }
    }
    if (immStatus::OK == imm && status == 0)
    {
        if (pd.polarizable())
        {
            // Make a copy since it maybe changed
            auto xx    = coords;
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
        if (strlen(rerun.trajectoryFileName()) > 0)
        {
            rerun.setFunctions(forceComp, &gendimers, oenv);
            rerun.setEInteraction(actmol.fragmentHandler()->topologies().size() > 1);
            rerun.rerun(logFile, &pd, &actmol, userqtot, qtot, verbose, sch.oneH());
        }
        else if (actmol.errors().empty())
        {
            MolHandler molhandler;
            std::vector<gmx::RVec> coords = actmol.xOriginal();
            std::vector<gmx::RVec> xmin   = coords;
            if (sch.minimize())
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
                eMin = molhandler.minimizeCoordinates(&pd, &actmol, forceComp, sch,
                                                      &xmin, &energies, logFile, freeze);
                if (eMinimizeStatus::OK == eMin)
                {
                    auto rmsd = molhandler.coordinateRmsd(&actmol, coords, &xmin);
                    fprintf(logFile, "Final energy: %g RMSD wrt original structure %g nm.\n",
                            energies[InteractionType::EPOT], rmsd);
                    auto nfrag = actmol.fragmentHandler()->topologies().size();
                    printf("There are %lu fragments\n", nfrag);
                    if (nfrag == 2)
                    {
                        std::map<InteractionType, double> einter;
                        std::vector<gmx::RVec>            interactionForces;
                        actmol.calculateInteractionEnergy(&pd, forceComp, &einter,
                                                          &interactionForces, &xmin);
                        for(const auto &ei : einter)
                        {
                            fprintf(logFile, "Interaction energy %s: %g\n",
                                    interactionTypeToString(ei.first).c_str(), ei.second);
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
                molhandler.simulate(&pd, &actmol, forceComp, sch, logFile,
                                    opt2fn("-o", fnm.size(),fnm.data()),
                                    opt2fn("-e", fnm.size(),fnm.data()),
                                    oenv);
            }
            auto confout = opt2fn_null("-c", fnm.size(),fnm.data());
            if (confout)
            {
                actmol.PrintConformation(confout, xmin, sch.writeShells(), box);
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
            for(const auto &err: actmol.errors())
            {
                fprintf(logFile, "%s\n", err.c_str());
            }
            status = 1;
        }
    }
    if (json)
    {
        jtree.write("simulate.json", json);
    }
    else
    {
        jtree.fwrite(logFile, json);
    }
    gmx_ffclose(logFile);
    if (immStatus::OK != imm)
    {
        printf("Simulate failed with an error:\n%s\n", immsg(imm));
    }
    return status;
}

} // namespace alexandria
