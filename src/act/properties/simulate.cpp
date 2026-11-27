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
#include <memory>
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
#include "act/forces/forcecomputerstatistics.h"
#include "act/forces/simd_information.h"
#include "act/import/compound_reader.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "act/properties/minimize.h"
#include "act/properties/pathintegralmd.h"
#include "act/properties/pathintegralthermo.h"
#include "act/properties/secondvirial.h"
#include "act/properties/velocityhandler.h"
#include "act/topology/actmol.h"
#include "act/topology/fragmenthandler.h"
#include "act/utility/jsontree.h"
#include "act/utility/output_environment.h"
#include "act/utility/stringutil.h"
#include "act/utility/timings.h"
#include "act/utility/units.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "simulate.h"

namespace alexandria
{

/*! \brief Tool to perform simple simulations and minimations
 * \param[in] argc Number of arguments on the cmd line
 * \param[in] argv The actual arguments
 * \return 0 if ok, 1 otherwise
 */
namespace
{

//! \return the name of a file for the functions that take nullptr for none
const char *fileOrNull(const std::string &fileName)
{
    return fileName.empty() ? nullptr : fileName.c_str();
}

/*! \brief Simulate at constant temperature with beads, see pathIntegralMD
 *
 * Simulates the compounds together and, on request, each compound alone,
 * and writes the averages and the binding energy.
 * \param[in]    msghandler For output and warnings
 * \param[in]    pd         The force field
 * \param[in]    actmol     The compounds
 * \param[in]    x          The starting coordinates
 * \param[in]    forceComp  The force computer
 * \param[in]    sch        The settings
 * \param[in]    options    The files
 * \param[in]    oenv       The output environment
 * \param[inout] jtree      The averages are added
 */
void simulateWithBeads(MsgHandler                    *msghandler,
                       const ForceField              *pd,
                       const ACTMol                  *actmol,
                       const std::vector<gmx::RVec>  &x,
                       const ForceComputer           *forceComp,
                       const SimulationConfigHandler &sch,
                       const SimulateOptions         &options,
                       const gmx_output_env_t        *oenv,
                       JsonTree                      *jtree)
{
    const auto        fragments = actmol->fragmentHandler();
    const auto       &tops      = fragments->topologies();
    // The first particle of each compound, and the number of particles
    std::vector<size_t> atomStart = { 0 };
    for(const auto &top : tops)
    {
        atomStart.push_back(atomStart.back() + top.atoms().size());
    }
    if (atomStart.back() != actmol->atomsConst().size() || atomStart.back() != x.size())
    {
        GMX_THROW(gmx::InternalError("The compounds do not have the particles of the molecule"));
    }
    PathIntegralFiles files;
    files.trajectory = fileOrNull(options.trajectoryFile);
    files.energy     = fileOrNull(options.energyFile);
    files.progress   = fileOrNull(options.progressFile);
    files.histogram  = fileOrNull(options.histogramFile);
    files.oenv       = oenv;
    // The thermochemistry is that of one compound, and starts from its
    // normal modes: a structure that is no minimum is refused before the
    // simulations
    const bool        thermo = sch.thermochemistry() && tops.size() == 1;
    HarmonicReference harmonic;
    if (sch.thermochemistry() && !thermo)
    {
        msghandler->msg(ACTStatus::Warning,
                        gmx::formatString("There are %zu compounds, the thermochemistry (-thermo) is that of one compound and is not computed.",
                                          tops.size()));
    }
    if (thermo)
    {
        harmonic = harmonicReference(msghandler, pd, actmol, forceComp, x, sch.temperature());
    }
    const auto all   = pathIntegralMD(msghandler, pd, actmol->topology(), atomStart, x,
                                      forceComp, sch, sch.seed(), files);
    JsonTree   tree("Simulation with beads");
    reportPathIntegral(msghandler, tops.size() > 1 ? "all compounds" : "the compound", all, &tree);
    if (thermo && all.etot.n > 0)
    {
        // The simulations for the entropy write no files, a program that
        // follows the progress gets a line after each
        const auto progress = [&files](int done, int total)
        {
            if (files.progress)
            {
                FILE *fp = gmx_ffopen(files.progress, "a");
                fprintf(fp, "{\"entropy\": {\"done\": %d, \"total\": %d}}\n", done, total);
                gmx_ffclose(fp);
            }
        };
        const auto result = pathIntegralThermochemistry(msghandler, pd, actmol, x, forceComp, sch,
                                                        harmonic, all, progress);
        reportThermochemistry(msghandler, result, &tree);
    }
    if (sch.monomers() && tops.size() < 2)
    {
        msghandler->msg(ACTStatus::Warning, "There is one compound, hence no binding energy (-monomers).");
    }
    else if (sch.monomers())
    {
        // Each compound alone, without files or wall, and with other random numbers
        std::vector<PathIntegralResults> alone;
        auto                             schAlone = sch;
        schAlone.setWall(0, sch.kwall());
        for(size_t f = 0; f < tops.size(); f++)
        {
            const std::vector<gmx::RVec> xf(x.begin() + atomStart[f], x.begin() + atomStart[f + 1]);
            const std::vector<size_t>    one  = { 0, xf.size() };
            const int                    seed = sch.seed() == 0 ? 0 : sch.seed() + 1 + static_cast<int>(f);
            msghandler->write(gmx::formatString("\nCompound %zu alone.", f + 1));
            alone.push_back(pathIntegralMD(msghandler, pd, &tops[f], one, xf, forceComp, schAlone, seed,
                                           PathIntegralFiles()));
            reportPathIntegral(msghandler, gmx::formatString("compound %zu", f + 1), alone.back(), &tree);
        }
        AverageError deltaH;
        const auto   deltaU = bindingEnergy(all, alone, &deltaH);
        msghandler->write(gmx::formatString("\nBinding at %g K with %d bead%s: the total energy of all compounds minus those of the compounds alone.",
                                            all.temperature, all.nbeads, all.nbeads == 1 ? "" : "s"));
        msghandler->write(gmx::formatString("  %-22s %14.4f +/- %.4f kJ/mol", "Energy", deltaU.average, deltaU.error));
        msghandler->write(gmx::formatString("  %-22s %14.4f +/- %.4f kJ/mol", "Enthalpy", deltaH.average, deltaH.error));
        JsonTree btree("Binding");
        btree.addObject("Unit", std::string("kJ/mol"));
        btree.addObject("Energy", deltaU.average);
        btree.addObject("Enthalpy", deltaH.average);
        btree.addObject("Standard error", deltaU.error);
        tree.addObject(btree);
        if (all.movedApart || all.wallFraction > 0)
        {
            msghandler->msg(ACTStatus::Warning,
                            gmx::formatString("The compounds %s, they are not bound all the time at %g K. The binding energy then depends on the volume that the wall leaves them, and is not that of a bound cluster.",
                                              all.movedApart ? "moved apart" : "reached the wall", all.temperature));
        }
    }
    jtree->addObject(tree);
}

} // namespace

void SimulateOptions::check() const
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

void SimulateOptions::checkSettings() const
{
    if (!std::isfinite(shellTolerance) || shellTolerance <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -shelltoler should be more than zero",
                                                           shellTolerance)));
    }
    simulation.validate();
}

int runSimulate(MsgHandler                *msghandler,
                const SimulateOptions     &options,
                gmx_output_env_t          *oenv,
                const CommunicationRecord *cr)
{
    options.check();
    // Time the phases, see printTimings
    PhaseTimer phase("Startup");
    // The settings may change below
    auto       sch    = options.simulation;
    const bool json   = !options.jsonFile.empty();
    int        status = 0;
    // A program without a command line has no output environment for the
    // energies
    const OutputEnvironment output(oenv);
    const double shellToler = sch.shellTolerance(options.shellTolerance, msghandler);
    ForceField        pd;
    readForceField(options.forceFieldFile, &pd, msghandler);
    (void) pd.verifyCheckSum(msghandler);
    for(const auto &ss : pd.info())
    {
        msghandler->write(ss);
    }

    ForceComputer forceComp;
    double epsilonr;
    epsilonr = pd.epsilonr();
    forceComp.init(shellToler, sch.shellMaxIter(), epsilonr);

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
        // The restraints move the center of mass and rotate around it,
        // removing that motion would counteract them
        if (COMMRemoval::None != sch.commRemoval())
        {
            msghandler->msg(ACTStatus::Warning,
                           "The motion of the center of mass is not removed, since there are position restraints.");
            sch.setCommRemoval(COMMRemoval::None);
        }
    }
    JsonTree jtree("simulate");
    if (msghandler->info())
    {
        forceFieldSummary(&jtree, &pd);
    }

    if (pd.polarizable())
    {
        phase.next("Polarizability");
        // Make a copy since it maybe changed
        const auto xx     = actmol.xOriginal();
        auto qprops = actmol.qProps();
        for(auto qp = qprops->begin(); qp < qprops->end(); ++qp)
        {
            auto qCalc = qp->qPact();
            qCalc->initializeMoments();
            qCalc->calcPolarizability(msghandler, &pd, actmol.topology(), &forceComp);
            auto alpha = qCalc->polarizabilityTensor();
            const std::string unit("A^3");
            double fac = convertFromGromacs(1, unit);
            JsonTree poltree("Polarizability");

            poltree.addValueUnit("XX", fac*alpha[XX][XX], unit);
            poltree.addValueUnit("YY", fac*alpha[YY][YY], unit);
            poltree.addValueUnit("ZZ", fac*alpha[ZZ][ZZ], unit);
            poltree.addValueUnit("Average", fac*qCalc->isotropicPolarizability(), unit);
            jtree.addObject(poltree);
        }
    }

    if (msghandler->debug())
    {
        actmol.topology()->dump(msghandler->twDebug());
    }
    auto   eMin = eMinimizeStatus::OK;
    // The potential energy after minimizing
    double epot = 0;
    /* Generate output file for debugging if requested */
    if (options.rerun && options.rerun->gendimersConst().hasTrajectory())
    {
        // The rerun keeps the force computer and the output environment
        auto rerun = *options.rerun;
        phase.next("Rerun");
        if (sch.nbeads() > 0)
        {
            msghandler->msg(ACTStatus::Warning,
                            "The energies of the structures of a trajectory are computed, there is no simulation with beads (-nbeads).");
        }
        rerun.setFunctions(&forceComp, output.get());
        rerun.setEInteraction(actmol.fragmentHandler()->topologies().size() > 1);
        rerun.rerun(msghandler, &pd, &actmol, msghandler->verbose());
    }
    else
    {
        const MolHandler molhandler;
        const std::vector<gmx::RVec> coords = actmol.xOriginal();
        std::vector<gmx::RVec> xmin   = coords;
        auto nfrag = actmol.fragmentHandler()->topologies().size();
        msghandler->msg(ACTStatus::Info, gmx::formatString("There are %lu fragments.", nfrag));
        if (sch.singlePoint())
        {
            phase.next("Single point");
            std::map<InteractionType, double> energies;
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            forceComp.compute(msghandler, &pd, actmol.topology(), &xmin, &forces, &energies);
            JsonTree jtener("Energies before");
            const std::string unit("kJ/mol");
            for (const auto &ener : energies)
            {
                const auto val = gmx::formatString("%.4f", ener.second);
                jtener.addValueUnit(interactionTypeToString(ener.first),
                                    val, unit);
            }
            jtree.addObject(jtener);
            if (nfrag == 2)
            {
                // Make a copy to keep the correct shell coordinates when saving the structure.
                std::vector<std::vector<gmx::RVec>>            newxmin = { xmin };
                std::vector<std::map<InteractionType, double>> einter;
                std::vector<double>                            epot;
                std::vector<std::vector<gmx::RVec>>            interactionForces;
                bool                                           separateInductionCorrection = false;
                forceComp.computeInteraction(msghandler, &pd, &actmol, &newxmin,
                                             separateInductionCorrection,
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
            }
        }
        else if (sch.minimize() || (sch.thermochemistry() && sch.nbeads() > 0 && sch.nsteps() > 0))
        {
            // The thermochemistry is relative to the minimum
            phase.next("Energy minimization");
            eMin = minimizeMolecule(msghandler, &pd, &forceComp, sch,
                                    fileOrNull(options.freezeFile),
                                    &actmol, &xmin, &jtree, &epot);
        }
        // Without a lower minimum the structure is usable, see minimumUsable
        if (minimumUsable(eMin) && sch.nsteps() > 0)
        {
            phase.next("Molecular dynamics");
            // Start from the minimized coordinates, if we minimized
            if (sch.nbeads() > 0)
            {
                simulateWithBeads(msghandler, &pd, &actmol, xmin, &forceComp, sch, options,
                                  output.get(), &jtree);
            }
            else
            {
                molhandler.simulate(msghandler, &pd, &actmol, xmin, &forceComp, sch,
                                    fileOrNull(options.trajectoryFile),
                                    fileOrNull(options.energyFile),
                                    output.get(),
                                    fileOrNull(options.progressFile));
            }
        }
        phase.next("Output");
        auto confout = fileOrNull(options.confout);
        if (confout && minimumUsable(eMin))
        {
            matrix box = { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } };
            actmol.PrintConformation(confout, xmin, sch.writeShells(), box);
        }
    }
    phase.next("Output");
    auto omp = fileOrNull(options.molpropFile);
    if (omp)
    {
        std::vector<MolProp> mpt;
        // We cannot use the actmol after this move command.
        mpt.push_back(std::move(*actmol.molProp()));
        MolPropWrite(omp, mpt, false);
    }

    auto imm = ACTMessage::OK;
    if (eMinimizeStatus::NoMinimum == eMin)
    {
        fprintf(stderr, "Minimization did not converge: %s, check log file %s\n",
                minimizeStatusMessage(eMin, epot).c_str(),
                msghandler->filename().c_str());
    }
    else if (eMinimizeStatus::OK != eMin)
    {
        fprintf(stderr, "Minimization failed: %s, check log file %s\n",
                minimizeStatusMessage(eMin, epot).c_str(),
                msghandler->filename().c_str());
        status = 1;
    }
    else if (ACTMessage::OK != imm)
    {
        msghandler->msg(ACTStatus::Error, ACTMessage::MinimizationFailed,
                       gmx::formatString("Please check the log file %s for error messages.\n",
                                         msghandler->filename().c_str()));
        status = 1;
    }
    if (pd.polarizable() && cr)
    {
        // The number of shell iterations per force calculation
        ForceComputerStatistics fcStats;
        msghandler->write(fcStats.statistics(cr, &forceComp, 0, true));
    }
    if (!json)
    {
        int indent = 0;
        msghandler->write(jtree.writeString(json, &indent));
    }
    phase.stop();
    // The timings go into the json file, or else only into the log file
    if (cr)
    {
        const auto peak = theoreticalPeakPerCore();
        printTimings(cr, msghandler->tw(), phase.wallTime(), json ? &jtree : nullptr, &peak);
    }
    if (json)
    {
        jtree.write(options.jsonFile, json);
    }
    if (!msghandler->ok())
    {
        printf("Simulate failed with an error:\n");
    }
    return status;
}


int simulate(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "alexandria simulate performs a proof-of-principle MD simulation, using",
        "a force field derived by the Alexandria Chemistry Toolkit.", 
        "The program can perform energy minimization and/or do a",
        "constant energy molecular dynamics simulation in vacuum.",
        "With [TT]-nbeads[tt] the simulation is at the constant temperature",
        "[TT]-temp[tt] instead, see below.",
        "In addition, a series of conformations (trajectory) may be",
        "submitted after which the energy per conformation is printed.[PAR]",
        "The input is given by a coordinate file, a force field file and",
        "command line options. During the simulation an energy file,",
        "a trajectory file and a log file are generated.[PAR]",
        "The minimization is the same as in alexandria minimize.",
        "During minimization a user select group of atoms can be frozen.",
        "To do so, supply an index file with atom numbers (first atom is 1)",
        "and numbering should disregard shells and virtual sites if present",
        "in the model.[PAR]",
        "To obtain information on the used model, the final structure can be",
        "store as a molprop file.[PAR]",
        "With [TT]-json[tt], part of the output, such as the timings, is",
        "written in json format to the given file, simulate.json by default.[PAR]",
        "With [TT]-progress[tt], a line with a json object is added to the",
        "given file each time the energies are written, with the step, the",
        "time, the temperature and the energies. The energy, trajectory and",
        "progress files are flushed a few times per second, such that",
        "another program can follow the simulation while it runs.[PAR]",
        "With [TT]-nbeads[tt] 1 the simulation is classical dynamics with a",
        "Langevin thermostat. With more beads it is path integral molecular",
        "dynamics: each atom is a ring polymer of beads, which gives the",
        "nuclear quantum effects on the energies and the structure, such as",
        "the zero point energy. Shells are relaxed and virtual sites are",
        "constructed for each bead. The charges are those of the input",
        "structure. The center of mass does not move, and [TT]-comm[tt] is not",
        "used. Constraints cannot be used. The averages over the steps after",
        "[TT]-equil[tt] are written to the log file, with standard errors from",
        "eight blocks. Since the thermostat changes the dynamics, only the",
        "averages have a meaning, not the time dependence.[PAR]",
        "The energy file then has the energy terms averaged over the beads,",
        "the kinetic energy from the centroid virial estimator (EKIN), which",
        "includes the translation of the center of mass, their sum (ETOT), the",
        "temperature of the momenta of the beads, the spread of the beads",
        "around their centroids (RGYR, nm) and, for more than one compound, the",
        "distance between the centers of mass of the compounds (DISTANCE, nm).",
        "The trajectory has the centroids of the beads.[PAR]",
        "In vacuum, compounds that move apart do not come back. A wall",
        "([TT]-rwall[tt], [TT]-kwall[tt]) keeps them together. It acts on the",
        "centers of mass and does not change the structure or orientation of",
        "the compounds. For two compounds [TT]-hist[tt] gives the distribution",
        "of the distance and the potential of mean force along it.[PAR]",
        "With [TT]-monomers[tt] each compound is simulated alone as well, and",
        "the binding energy is printed: the total energy of all compounds",
        "together minus those of the compounds alone, and the enthalpy, which",
        "has -kT in addition for each compound that binds. At low temperature",
        "it is the binding energy with zero point motion, -D0. This has a",
        "meaning only as long as the compounds stay together without the",
        "wall: the program warns if they reach it. Please check the convergence",
        "with the number of beads, which has to increase when the temperature",
        "decreases.[PAR]",
        "With [TT]-thermo[tt] the thermochemistry of one compound is computed",
        "from the simulations: the internal energy, the enthalpy, the heat",
        "capacity, the standard entropy and the Gibbs energy at 1 bar, printed",
        "next to those of the rigid rotor with harmonic oscillators of the",
        "normal mode analysis ([TT]alexandria nma[tt]). The compound is",
        "minimized first, and the energies are relative to the minimum. The",
        "values include a correction for the finite number of beads. The",
        "entropy is computed by thermodynamic integration, which takes",
        "2 x [TT]-nquad[tt] more simulations: classical ones at lower",
        "temperatures for the anharmonicity, and ones with beads in which",
        "the constant of Planck is scaled down for the quantum effects."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff",      "aff",        ffREAD  },
        { efPDB, "-o",       "trajectory", ffWRITE },
        { efSTO, "-c",       "confout",    ffOPTWR },
        { efXVG, "-e",       "energy",     ffWRITE },
        { efNDX, "-freeze",  "freeze",     ffOPTRD },
        { efXML, "-omp",     "molprop",    ffOPTWR },
        { efJSON, "-json",   "simulate",   ffOPTWR },
        { efJSONL, "-progress", "progress", ffOPTWR },
        { efXVG, "-hist",    "distance",   ffOPTWR }
    };
    gmx_output_env_t         *oenv;
    SimulateOptions           options;
    real                      shellToler = options.shellTolerance;
    std::vector<t_pargs>      pa = {
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)" }
    };
    options.simulation.add_options(&pa, &fnm);
    options.simulation.add_MD_options(&pa);
    options.simulation.add_bead_options(&pa);
    options.rerun.emplace(false);
    options.rerun->addOptions(&pa, &fnm, &desc);
    options.compounds.addOptions(&pa, &fnm, &desc);
    MsgHandler msghandler;
    msghandler.addOptions(&pa, &fnm, "simulate");
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
        options.rerun->finishOptions(fnm);
        print_header(msghandler.tw(), pa, fnm);
        options.shellTolerance = shellToler;
        options.forceFieldFile = opt2fn("-ff", fnm.size(), fnm.data());
        options.trajectoryFile = opt2fn("-o", fnm.size(), fnm.data());
        options.energyFile     = opt2fn("-e", fnm.size(), fnm.data());
        const auto optional    = [&fnm](const char *option)
        {
            const char *name = opt2fn_null(option, fnm.size(), fnm.data());
            return std::string(name ? name : "");
        };
        options.confout      = optional("-c");
        options.freezeFile   = optional("-freeze");
        options.molpropFile  = optional("-omp");
        options.jsonFile     = optional("-json");
        options.progressFile  = optional("-progress");
        options.histogramFile = optional("-hist");
        return runSimulate(&msghandler, options, oenv, &cr);
    }
    catch (const gmx::InvalidInputError &ex)
    {
        fprintf(stderr, "%s\n", ex.what());
        return 1;
    }
}

} // namespace alexandria
