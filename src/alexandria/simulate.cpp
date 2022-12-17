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
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/futil.h"

#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/poldata/poldata_xml.h"
#include "act/utility/jsontree.h"
#include "act/utility/stringutil.h"
#include "alexandria/alex_modules.h"
#include "alexandria/atype_mapping.h"
#include "alexandria/babel_io.h"
#include "alexandria/confighandler.h"
#include "alexandria/molhandler.h"
#include "alexandria/mymol.h"
#include "alexandria/tuning_utility.h"

namespace alexandria
{

static void forceFieldSummary(JsonTree      *jtree,
                              const Poldata *pd)
{
    jtree->addObject(JsonTree("Force field file", pd->filename()));
    jtree->addObject(JsonTree("Created", pd->timeStamp()));
    jtree->addObject(JsonTree("Checksum", pd->checkSum()));
    jtree->addObject(JsonTree("Polarizable", yesno_names[pd->polarizable()]));
    jtree->addObject(JsonTree("Charge generation", 
                              chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str()));
    jtree->addObject(JsonTree("# exclusions", gmx_itoa(pd->getNexcl())));
    jtree->addObject(JsonTree("Relative dielectric constant epsilon_r",
                              gmx_ftoa(pd->getEpsilonR())));
    jtree->addObject(JsonTree("# particle types", gmx_itoa(pd->getNatypes())));
    
    JsonTree ftree("InteractionTypes");
    for(const auto &fs : pd->forcesConst())
    {
        auto itype = fs.first;
        auto &ffpl = fs.second;
        if (!ffpl.parametersConst().empty())
        {
            JsonTree fftree(interactionTypeToString(itype));
            if (!ffpl.function().empty())
            {
                fftree.addObject(JsonTree("Force function", ffpl.function()));
            }
            fftree.addObject(JsonTree("# entries",
                                      gmx_itoa(ffpl.parametersConst().size())));
            ftree.addObject(fftree);
        }
    }
    if (!ftree.empty())
    {
        jtree->addObject(ftree);
    }
}

static double sphere_int(double r1, double r2, double val1, double val2)
{
    // Approximate trapezium by y = ax + b (a == slope)
    // then integrate that multiplied by x^2 to get
    // a/4 x^4 + b/3 x^3
    // insert old and new point
    double a        = (val2-val1)/(r2-r1);
    double b        = val1 - a*r1;
    double integral = ((a/4)*(r2*r2*r2*r2 - r1*r1*r1*r1) + 
                       (b/3)*(r2*r2*r2 - r1*r1*r1));
    return 4*M_PI*integral;
}

static void computeB2(FILE                         *logFile,
                      const char                   *ehisto,
                      const char                   *b2file,
                      gmx_stats                     edist,
                      gmx_output_env_t             *oenv,
                      const std::vector<double>    &Temperature,
                      double                        mass,
                      const gmx::RVec              &inertia,
                      const std::vector<gmx::RVec> &force1,
                      const std::vector<gmx::RVec> &torque1)
{
    auto                      N = edist.get_npoints();
    const std::vector<double> x = edist.getX();
    const std::vector<double> y = edist.getY();
    double xmin                 = *std::min_element(x.begin(), x.end());
    double xmax                 = *std::max_element(x.begin(), x.end());
    if (N > 2 && xmax > xmin)
    {
        FILE *fp = nullptr;
        if (ehisto)
        {
            fp = xvgropen(ehisto, "Mayer function", "r (nm)",
                          "< exp[-U12/kBT]-1 >", oenv);
            std::vector<std::string> label;
            for(auto T : Temperature)
            {
                if (T != 0)
                {
                    label.push_back(gmx::formatString("T = %g K", T));
                }
            }
            xvgrLegend(fp, label, oenv);
        }
        FILE *b2p = nullptr;
        if (b2file)
        {
            b2p = xvgropen(b2file, "Second virial coefficient",
                           "Temperature (K)", "B2(T) cm^3/mol", oenv);
        }
        real   binwidth = 0.01; // nm
        size_t nbins    = 1+std::round((xmax-xmin)/binwidth);
        binwidth        = (xmax-xmin)/(nbins-1);
        for(auto T : Temperature)
        {
            if (T == 0)
            {
                fprintf(stderr, "Please provide a finite temperature to compute second virial.\n");
                continue;
            }
            std::vector<double>    exp_U12(nbins, 0.0);
            std::vector<double>    exp_F2(nbins, 0.0);
            std::vector<gmx::RVec> exp_tau2(nbins, { 0.0, 0.0, 0.0 });
            std::vector<int>       n_U12(nbins, 0);
            double beta = 1.0/(BOLTZ*T);
            for(size_t ii = 0; ii < x.size(); ii++)
            {
                double rindex = (x[ii]-xmin)/binwidth;
                size_t index  = rindex;
                // Gray and Gubbins Eqn. 3.261
                double g0_12 = std::exp(-y[ii]*beta);
                // Gray and Gubbins Eqn. 3.272
                exp_U12[index] += g0_12-1;
                // Gray and Gubbins Eqn. 3.281
                exp_F2[index]  += g0_12*iprod(force1[ii], force1[ii]);
                for(int m = 0; m < DIM; m++)
                {
                    // Gray and Gubbins Eqn. 3.282
                    exp_tau2[index][m] += g0_12*torque1[ii][m]*torque1[ii][m];
                }
                n_U12[index]   += 1;
            }
            double    Bclass    =  0;
            double    BqmForce  =  0;
            double    BqmTorque =  0;
            // We start in the origin even if there is no data.
            double    r1        =  0;
            double    Uprev     = -1;
            double    Fprev     =  0;
            gmx::RVec Tprev     = { 0, 0, 0 };
            if (fp)
            {
                fprintf(fp, "@ type xy\n");
                fprintf(fp, "%10g  %10g\n", r1, Uprev);
            }
            double hbarfac  = gmx::square(PLANCK/(M_PI*beta))/24;
            for(size_t ii = 0; ii < nbins; ii++)
            {
                double r2 = xmin+ii*binwidth;
                if (n_U12[ii] > 0)
                {
                    double Unew = exp_U12[ii]/n_U12[ii];
                    if (fp)
                    {
                        fprintf(fp, "%10g  %10g\n", r2, Unew);
                    }
                    auto dB       = sphere_int(r1, r2, Uprev, Unew);
                    // TODO: Should really be factor 0.5 here!
                    Bclass       -= dB;
                    Uprev         = Unew;
                    double Fnew   = exp_F2[ii]/(mass*n_U12[ii]);
                    BqmForce     += hbarfac*sphere_int(r1, r2, Fprev, Fnew);
                    Fprev         = Fnew;
                    for(int m = 0; m < DIM; m++)
                    {
                        if (inertia[m] > 0)
                        {
                            double Tnew  = exp_tau2[ii][m]/inertia[m];
                            BqmTorque   += hbarfac*sphere_int(r1, r2, Tprev[m], Tnew);
                            Tprev[m]     = Tnew;
                        }
                    }
                    r1            = r2;
                }
            }
            if (fp)
            {
                fprintf(fp, "&\n");
            }
            double Btot = (Bclass + BqmForce + BqmTorque)*AVOGADRO*1e-21;
            fprintf(logFile, "T = %g K. Classical second virial coefficient B2cl %g BqmForce %g BqmTorque %g nm^3 Total %g cm^3/mol\n", T,
                    Bclass, BqmForce, BqmTorque, Btot);
            if (b2p)
            {
                fprintf(b2p, "%10g  %10g\n", T, Btot);
            }
        }
        if (fp)
        {
            xvgrclose(fp);
        }
        if (b2p)
        {
            xvgrclose(b2p);
        }
    }
}

static void do_rerun(FILE                      *logFile,
                     const Poldata             *pd,
                     const MyMol               *mymol,
                     ForceComputer             *forceComp,
                     const char                *trajname,
                     const char                *ehisto,
                     const char                *b2file,
                     bool                       eInter,
                     double                     qtot,
                     gmx_output_env_t          *oenv,
                     const std::vector<double> &Temperature)
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
        gmx_stats edist;
        std::vector<gmx::RVec> force1;
        std::vector<gmx::RVec> torque1;
        gmx::RVec              inertia  = { 0, 0, 0 };
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
                fprintf(logFile, "%5d", mp_index);
                if (eInter)
                {
                    auto EE         = mymol->calculateInteractionEnergy(pd, forceComp, &forces, &coords);
                    auto atomStart  = mymol->fragmentHandler()->atomStart();
                    if (atomStart.size() != 3)
                    {
                        GMX_THROW(gmx::InvalidInputError(gmx::formatString("This is not a dimer, there are %zu fragments instead of 2", atomStart.size()-1).c_str()));
                    }
                    gmx::RVec f[2]    = { { 0, 0, 0 }, { 0, 0, 0 } };
                    gmx::RVec com[2]  = { { 0, 0, 0 }, { 0, 0, 0 } };
                    double    mtot[2] = { 0, 0 };
                    auto      tops  = mymol->fragmentHandler()->topologies();
                    for(int kk = 0; kk < 2; kk++)
                    {
                        auto atoms = tops[kk].atoms();
                        for(size_t i = atomStart[kk]; i < atomStart[kk+1]; i++)
                        {
                            gmx::RVec mr1;
                            auto      mi = atoms[i-atomStart[kk]].mass();
                            svmul(mi, coords[i], mr1);
                            mtot[kk] += mi;
                            rvec_inc(com[kk], mr1);
                            rvec_inc(f[kk], forces[i]);
                        }
                        GMX_RELEASE_ASSERT(mtot[kk] > 0, "Zero mass");
                        for(size_t m = 0; m < DIM; m++)
                        {
                            com[kk][m] /= mtot[kk];
                        }
                    }
                    force1.push_back(f[0]);
                    gmx::RVec torque = { 0, 0, 0 };
                    for(size_t i = atomStart[0]; i < atomStart[1]; i++)
                    {
                        gmx::RVec ri;
                        rvec_sub(coords[i], com[0], ri);
                        gmx::RVec ti;
                        cprod(ri, forces[i], ti);
                        rvec_inc(torque, ti);
                        for(int m = 0; m < DIM; m++)
                        {
                            inertia[m] += ri[m]*ri[m]*atoms[i].mass();
                        }
                    }
                    torque1.push_back(torque);
                    gmx::RVec dcom;
                    rvec_sub(com[0], com[1], dcom);
                    double rcom = norm(dcom);
                    fprintf(logFile, " r %g Einter %g Force %g %g %g Torque %g %g %g",
                            rcom, EE, f[0][XX], f[0][YY], f[0][ZZ], torque[XX], torque[YY], torque[ZZ]);
                    edist.add_point(rcom, EE, 0, 0);
                }
                else
                {
                    forceComp->compute(pd, mymol->topology(),
                                       &coords, &forces, &energies);
                    for(const auto &ee : energies)
                    {
                        fprintf(logFile, "  %s %8g", interactionTypeToString(ee.first).c_str(), ee.second);
                    }
                }
                fprintf(logFile, "\n");
            }
            mp_index += 1;
        }
        if (eInter)
        {
            for(int m = 0; m < DIM; m++)
            {
                inertia[m] /= edist.get_npoints();
            }
            computeB2(logFile, ehisto, b2file, edist, oenv, Temperature,
                      mymol->fragmentHandler()->topologies()[0].mass(),
                      inertia, force1, torque1);
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
        "The program can perform energy minimization and normal mode",
        "analysis including thermochemistry calculations, or do a",
        "constant energy molecular dynamics simulation in vacuum.",
        "In addition, a series of conformations (trajectory) may be",
        "submitted after which the energy per conformation is printed.",
        "If the input trajectory consists of dimers, the second virial",
        "coefficient can be estimated as well.[PAR]",
        "The input is given by a coordinate file, a force field file and",
        "command line options. During the simulation an energy file,",
        "a trajectory file and a log file are generated. If a trajectory",
        "of dimers is presented as input for energy calculations, the",
        "corresponding molecule file, used for generating the topology",
        "needs to be a molprop (xml) file and contain information about",
        "the compounds in the dimer. Based on dimer energies the second",
        "virial coefficient will be estimated. Mayer curves can optionally",
        "be plotted and the second viriail can be plotted as a function",
        "of temperature."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff", "gentop",     ffREAD  },
        { efXML, "-mp", "molprop",    ffOPTRD },
        { efPDB, "-o",  "trajectory", ffWRITE },
        { efSTO, "-c",  "confout",    ffOPTWR },
        { efXVG, "-e",  "energy",     ffWRITE },
        { efXVG, "-eh", "mayer",      ffOPTWR },
        { efXVG, "-b2", "B2T",        ffOPTWR },
        { efLOG, "-g",  "simulation", ffWRITE },
        { efXVG, "-ir", "IRspectrum", ffOPTWR }
    };
    gmx_output_env_t         *oenv;
    static char              *filename   = (char *)"";
    static char              *trajname   = (char *)"";
    static char              *molnm      = (char *)"";
    static char              *qqm        = (char *)"";
    double                    qtot       = 0;
    double                    shellToler = 1e-6;
    double                    T1         = 0;
    double                    T2         = 0;
    double                    deltaT     = 10;
    bool                      verbose    = false;
    bool                      eInter     = false;
    bool                      json       = false;
    std::vector<t_pargs>      pa = {
        { "-f",      FALSE, etSTR,  {&filename},
          "Molecular structure file in e.g. pdb format" },
        { "-traj",   FALSE, etSTR,  {&trajname},
          "Trajectory or series of structures of the same compound for which the energies will be computed. If this option is present, no simulation will be performed." },
        { "-T1",     FALSE, etREAL, {&T1},
          "Starting temperature for second virial calculations" },
        { "-T2",     FALSE, etREAL, {&T2},
          "Starting temperature for second virial calculations" },
        { "-dT",     FALSE, etREAL, {&deltaT},
          "Temperature increment for calculation of second virial" },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Combined charge of the molecule(s). This will be taken from the input file by default, but that is not always reliable." },
        { "-qqm",    FALSE, etSTR,  {&qqm},
          "Use a method from quantum mechanics that needs to be present in the input file. Either ESP, Hirshfeld, CM5 or Mulliken may be available." },
        { "-einter", FALSE, etBOOL, {&eInter},
          "Compute dimer interaction energies when doing a rerun" },
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

    if (opt2bSet("-mp", fnm.size(), fnm.data()) && strlen(filename) > 0)
    {
        fprintf(stderr, "Please supply either a molprop file (-mp option) or an input filename (-f option), but not both.\n");
        status = 1;
        return status;
    }
    else if (!opt2bSet("-mp", fnm.size(), fnm.data()) && strlen(filename) == 0)
    {
        fprintf(stderr, "Please supply either a molprop file (-mp option) or an input filename (-f option)\n");
        status = 1;
        return status;
    }
    
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
        if (opt2bSet("-mp", fnm.size(), fnm.data()))
        {
            MolPropRead(opt2fn("-mp", fnm.size(), fnm.data()), &mps);
        }
        else
        {
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
        if (pd.polarizable())
        {
            // Make a copy since it maybe changed
            auto xx    = coords;
            auto qCalc = mymol.qTypeProps(qType::Calc);
            qCalc->initializeMoments();
            forceComp->calcPolarizability(&pd, mymol.topology(), &xx, qCalc);
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

        if (debug)
        {
            mymol.topology()->dump(debug);
        }
        auto eMin = eMinimizeStatus::OK;
        /* Generate output file for debugging if requested */
        if (strlen(trajname) > 0)
        {
            std::vector<double> Temperature;
            Temperature.push_back(T1);
            if (deltaT > 0)
            {
                double T = T1+deltaT;
                while (T <= T2)
                {
                    Temperature.push_back(T);
                    T += deltaT;
                }
            }
            else
            {
                Temperature.push_back(T2);
            }
            do_rerun(logFile, &pd, &mymol, forceComp, trajname,
                     opt2fn_null("-eh", fnm.size(), fnm.data()),
                     opt2fn_null("-b2", fnm.size(), fnm.data()),
                     eInter, qtot, oenv, Temperature);
        }
        else if (mymol.errors().empty())
        {
            MolHandler molhandler;
            std::vector<gmx::RVec> coords = mymol.xOriginal();
            std::vector<gmx::RVec> xmin   = coords;
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
                    JsonTree jtener("Energies");
                    std::string unit("kJ/mol");
                    for (const auto &ener : energies)
                    {
                        jtener.addValueUnit(interactionTypeToString(ener.first),
                                            gmx_ftoa(ener.second), unit);
                    }
                    jtree.addObject(jtener);
                    matrix box;
                    clear_mat(box);
                    write_sto_conf(opt2fn("-c", fnm.size(),fnm.data()), 
                                   mymol.getMolname().c_str(),
                                   mymol.gmxAtomsConst(),
                                   as_rvec_array(xmin.data()), nullptr,
                                   epbcNONE, box);
                    
                    if (sch.nma())
                    {
                        AtomizationEnergy        atomenergy;
                        doFrequencyAnalysis(&pd, &mymol, molhandler, forceComp, &coords,
                                            atomenergy, nullptr, &jtree,
                                            opt2fn_null("-ir", fnm.size(), fnm.data()),
                                            sch.lineWidth(), oenv,
                                            sch.lapack(), verbose);
                    }
                }
            }
            if (!sch.nma() && eMinimizeStatus::OK == eMin)
            {
                molhandler.simulate(&pd, &mymol, forceComp, sch, logFile,
                                    opt2fn("-o", fnm.size(),fnm.data()),
                                    opt2fn("-e", fnm.size(),fnm.data()),
                                    oenv);
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
        jtree.write("simulate.json", json);
    }
    else
    {
        jtree.fwrite(logFile, json);
    }
    gmx_ffclose(logFile);
    return status;
}

} // namespace alexandria
