/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <map>
#include <nlohmann/json.hpp>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/pdbwriter.h"
#include "act/basics/act_particle.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/properties/fragmentwall.h"
#include "act/properties/ringpolymer.h"
#include "act/topology/actatom.h"
#include "act/topology/topology.h"
#include "act/utility/jsontree.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"
#include "pathintegralmd.h"

namespace alexandria
{

namespace
{

//! The fragments moved apart if their distance grew by this much (nm), without a wall
const double apartDistance = 1.0;

//! The fragments are not bound if they are at the wall in more than this fraction of the steps
const double wallFractionLimit = 0.01;

//! The largest number of bins of the histogram of the distance
const size_t maxBins = 50000;

/*! \brief The names of the columns of the energy file
 * \param[in] terms    The energy terms
 * \param[in] wall     Whether there is a wall
 * \param[in] distance Whether there is a distance between fragments
 * \return the names
 */
std::vector<std::string> energyNames(const std::map<InteractionType, double> &terms,
                                     bool                                     wall,
                                     bool                                     distance)
{
    std::vector<std::string> names;
    for(const auto &t : terms)
    {
        names.push_back(interactionTypeToString(t.first));
    }
    if (wall)
    {
        names.emplace_back("WALL");
    }
    for(const char *name : { "EKIN", "ETOT", "TEMPERATURE", "RGYR" })
    {
        names.emplace_back(name);
    }
    if (distance)
    {
        names.emplace_back("DISTANCE");
    }
    return names;
}

//! What is known of a configuration of the ring polymers
struct Sample
{
    //! The step
    int    step        = 0;
    //! The potential energy of the force field
    double epot        = 0;
    //! The energy of the wall
    double ewall       = 0;
    //! The kinetic energy
    double ekin        = 0;
    //! The temperature of the momenta
    double temperature = 0;
    //! The spread of the beads
    double gyration    = 0;
    //! The distance between the fragments
    double distance    = 0;
};

/*! \brief Write the progress of a simulation as one line with a JSON object
 * As MolHandler::simulate does, with the number of beads, the spread of
 * the beads and the distance between the fragments in addition.
 * \param[in] fp       The file
 * \param[in] sample   The configuration
 * \param[in] nsteps   The number of steps of the simulation
 * \param[in] deltat   The time step (ps)
 * \param[in] nbeads   The number of beads
 * \param[in] terms    The energy terms
 * \param[in] wall     Whether there is a wall
 * \param[in] distance Whether there is a distance between fragments
 */
void printProgress(FILE                                    *fp,
                   const Sample                            &sample,
                   int                                      nsteps,
                   double                                   deltat,
                   int                                      nbeads,
                   const std::map<InteractionType, double> &terms,
                   bool                                     wall,
                   bool                                     distance)
{
    nlohmann::ordered_json js;
    js["step"]        = sample.step;
    js["nsteps"]      = nsteps;
    js["time"]        = { { "value", sample.step*deltat }, { "unit", "ps" } };
    js["temperature"] = { { "value", sample.temperature }, { "unit", "K" } };
    nlohmann::ordered_json jener;
    jener["unit"] = "kJ/mol";
    for (const auto &e : terms)
    {
        jener[interactionTypeToString(e.first)] = e.second;
    }
    if (wall)
    {
        jener["WALL"] = sample.ewall;
    }
    jener["EKIN"]  = sample.ekin;
    jener["ETOT"]  = sample.ekin + sample.epot;
    js["energies"] = jener;
    js["beads"]    = nbeads;
    js["gyration"] = { { "value", sample.gyration }, { "unit", "nm" } };
    if (distance)
    {
        js["distance"] = { { "value", sample.distance }, { "unit", "nm" } };
    }
    // On one line, replacing what is not a valid character or number
    fprintf(fp, "%s\n", js.dump(-1, ' ', false, nlohmann::json::error_handler_t::replace).c_str());
}

/*! \brief Write the histogram of the distance between two fragments
 * With the potential of mean force along the distance, from the
 * probability density divided by the squared distance, which is zero at
 * its minimum.
 * \param[in] file    The file
 * \param[in] oenv    The output environment
 * \param[in] results The results with the histogram
 */
void writeHistogram(const char                *file,
                    const gmx_output_env_t    *oenv,
                    const PathIntegralResults &results)
{
    const auto &histo = results.histogram;
    double      total = 0;
    for(auto h : histo)
    {
        total += h;
    }
    if (total <= 0)
    {
        return;
    }
    double              kT = BOLTZ*results.temperature;
    double              dr = results.binWidth;
    std::vector<double> pmf(histo.size(), 0.0);
    double              pmfMin = 0;
    bool                first  = true;
    for(size_t i = 0; i < histo.size(); i++)
    {
        if (histo[i] > 0)
        {
            double r = (i + 0.5)*dr;
            pmf[i]   = -kT*std::log(histo[i]/(total*dr*r*r));
            if (first || pmf[i] < pmfMin)
            {
                pmfMin = pmf[i];
                first  = false;
            }
        }
    }
    FILE *fp = xvgropen(file, "Distance between the centers of mass", "R (nm)", "", oenv);
    std::vector<const char *> legend = { "Probability density (1/nm)", "Potential of mean force (kJ/mol)" };
    xvgr_legend(fp, legend.size(), legend.data(), oenv);
    for(size_t i = 0; i < histo.size(); i++)
    {
        // The potential of mean force is only known where the fragments were
        if (histo[i] > 0)
        {
            fprintf(fp, "%10g  %12g  %12g\n", (i + 0.5)*dr, histo[i]/(total*dr), pmf[i] - pmfMin);
        }
    }
    xvgrclose(fp);
}

} // namespace

BlockAverager::BlockAverager(size_t nexpected, size_t nblocks) :
    nexpected_(nexpected), sum_(std::max<size_t>(nblocks, 1), 0.0), count_(std::max<size_t>(nblocks, 1), 0)
{
}

void BlockAverager::add(double value)
{
    size_t nblocks = sum_.size();
    size_t block   = nblocks - 1;
    if (n_ < nexpected_)
    {
        block = (n_*nblocks)/nexpected_;
    }
    sum_[block]   += value;
    count_[block] += 1;
    total_        += value;
    n_            += 1;
}

AverageError BlockAverager::result() const
{
    AverageError ae;
    ae.n = n_;
    if (0 == n_)
    {
        return ae;
    }
    ae.average     = total_/n_;
    // The blocks that got all their values
    size_t nblocks = sum_.size();
    size_t ndone   = nblocks;
    if (n_ < nexpected_)
    {
        ndone = (n_*nblocks)/nexpected_;
    }
    std::vector<double> means;
    for(size_t b = 0; b < ndone; b++)
    {
        if (count_[b] > 0)
        {
            means.push_back(sum_[b]/count_[b]);
        }
    }
    if (means.size() >= 2)
    {
        double mean = 0, var = 0;
        for(auto m : means)
        {
            mean += m/means.size();
        }
        for(auto m : means)
        {
            var += (m - mean)*(m - mean);
        }
        ae.error = std::sqrt(var/(means.size()*(means.size() - 1)));
    }
    return ae;
}

CovarianceAverager::CovarianceAverager(size_t nexpected, size_t nblocks) :
    nexpected_(nexpected), block_(std::max<size_t>(nblocks, 1))
{
}

void CovarianceAverager::add(double a, double b)
{
    size_t nblocks = block_.size();
    size_t block   = nblocks - 1;
    if (total_.n < nexpected_)
    {
        block = (total_.n*nblocks)/nexpected_;
    }
    for(Sums *sums : { &total_, &block_[block] })
    {
        sums->n  += 1;
        sums->a  += a;
        sums->b  += b;
        sums->ab += a*b;
    }
}

AverageError CovarianceAverager::result() const
{
    AverageError ae;
    ae.n = total_.n;
    if (0 == total_.n)
    {
        return ae;
    }
    ae.average     = total_.covariance();
    // The blocks that got all their values
    size_t nblocks = block_.size();
    size_t ndone   = nblocks;
    if (total_.n < nexpected_)
    {
        ndone = (total_.n*nblocks)/nexpected_;
    }
    std::vector<double> cov;
    for(size_t b = 0; b < ndone; b++)
    {
        if (block_[b].n > 1)
        {
            cov.push_back(block_[b].covariance());
        }
    }
    if (cov.size() >= 2)
    {
        double mean = 0, var = 0;
        for(auto c : cov)
        {
            mean += c/cov.size();
        }
        for(auto c : cov)
        {
            var += (c - mean)*(c - mean);
        }
        ae.error = std::sqrt(var/(cov.size()*(cov.size() - 1)));
    }
    return ae;
}

PathIntegralResults pathIntegralMD(MsgHandler                    *msghandler,
                                   const ForceField              *pd,
                                   const Topology                *top,
                                   const std::vector<size_t>     &atomStart,
                                   const std::vector<gmx::RVec>  &x,
                                   const ForceComputer           *forceComp,
                                   const SimulationConfigHandler &simConfig,
                                   int                            seed,
                                   const PathIntegralFiles       &files,
                                   const SimulationObserver      *observer,
                                   double                         planckScale)
{
    PathIntegralResults results;
    const int nsteps  = simConfig.nsteps();
    const int nbeads  = simConfig.nbeads();
    if (nsteps <= 0 || nbeads <= 0)
    {
        return results;
    }
    const int    nequil  = simConfig.nequil();
    const int    nstener = simConfig.nstener();
    const int    nstxout = simConfig.nstxout();
    const double deltat  = simConfig.deltat();
    RingPolymerSettings settings;
    settings.nbeads      = nbeads;
    settings.temperature = simConfig.temperature();
    settings.timeStep    = deltat;
    settings.centroidTau = simConfig.tau();
    settings.seed        = seed;
    settings.planckScale = planckScale;
    RingPolymerSampler sampler(pd, top, forceComp, settings);
    // The fragments, and the wall that keeps them together
    const auto         &atoms = top->atoms();
    std::vector<double> mass;
    for(const auto &atom : atoms)
    {
        mass.push_back(atom.mass());
    }
    const FragmentWall wall(atomStart, mass, simConfig.rwall(), simConfig.kwall());
    const bool         haveDistance = wall.nfragments() > 1;
    if (wall.active())
    {
        sampler.setCentroidPotential([&wall](const std::vector<gmx::DVec> &xc, std::vector<gmx::DVec> *fc)
                                     { return wall.compute(xc, fc); });
    }
    results.nbeads      = nbeads;
    results.temperature = settings.temperature;
    results.nequil      = nequil;
    results.nfragments  = wall.nfragments();
    results.wall        = wall.active();

    // A program that observes the simulation may do without files
    FILE *traj     = files.trajectory ? gmx_ffopen(files.trajectory, "w") : nullptr;
    FILE *exvg     = files.energy ? xvgropen(files.energy, "ACT Energies", "Time (ps)",
                                             "Energy (kJ/mol)", files.oenv) : nullptr;
    FILE *progress = files.progress ? gmx_ffopen(files.progress, "w") : nullptr;
    // The files are flushed a few times per second, such that another
    // program can follow the simulation
    auto       lastFlush = std::chrono::steady_clock::now();
    const auto flush     = [&](bool now)
    {
        const auto time = std::chrono::steady_clock::now();
        if (now || time - lastFlush >= std::chrono::milliseconds(200))
        {
            for(FILE *fp : { traj, exvg, progress })
            {
                if (fp)
                {
                    fflush(fp);
                }
            }
            lastFlush = time;
        }
    };
    // For trajectory writing
    matrix           box   = { { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 } };
    const char      *title = "ACT trajectory of centroids";
    std::vector<int> trajIndex;
    for(size_t i = 0; i < atoms.size(); i++)
    {
        if (simConfig.writeShells() || ActParticle::Atom == atoms[i].pType())
        {
            trajIndex.push_back(i);
        }
    }
    // The averages over the steps after the equilibration
    const size_t  nsample = nsteps - nequil;
    size_t        nterm   = 0;
    for(int step = nequil; nstener > 0 && step < nsteps; step++)
    {
        nterm += (step % nstener == 0) ? 1 : 0;
    }
    BlockAverager epot(nsample), ekin(nsample), etot(nsample), gyration(nsample), distance(nsample);
    BlockAverager ekinPrimitive(nsample);
    // The fluctuations of the energy, relative to the first value since
    // the differences are small compared to the energy itself
    CovarianceAverager fluctuation(nsample);
    double             energy0     = 0;
    bool               haveEnergy0 = false;
    std::map<InteractionType, BlockAverager> terms;
    size_t                 nwall = 0;
    std::vector<gmx::DVec> xc;
    std::vector<gmx::RVec> xr(atoms.size());
    double                 distance0  = 0;
    bool                   legendDone = false;
    bool                   keepGoing  = true;

    msghandler->write(gmx::formatString("\nWill start a simulation of %d steps at %g K with %d bead%s: %s.", nsteps,
                                        settings.temperature, nbeads, nbeads == 1 ? "" : "s",
                                        nbeads == 1 ? "classical dynamics with a thermostat" : "path integral molecular dynamics"));
    if (wall.active())
    {
        msghandler->write(gmx::formatString("A wall with the radius %g nm and the force constant %g kJ/mol nm^2 keeps the %zu compounds together.",
                                            simConfig.rwall(), simConfig.kwall(), wall.nfragments()));
    }
    else if (simConfig.rwall() > 0)
    {
        msghandler->msg(ACTStatus::Warning, "There is one compound, the wall (-rwall) is not used.");
    }
    sampler.start(msghandler, x, nstener > 0);
    int step = 0;
    for (; keepGoing && step < nsteps; step++)
    {
        // The configuration at this time, before integrating
        Sample sample;
        sample.step        = step;
        sample.epot        = sampler.potentialEnergy();
        sample.ewall       = sampler.centroidEnergy();
        sample.ekin        = sampler.kineticEnergy();
        sample.temperature = sampler.beadTemperature()/nbeads;
        sample.gyration    = std::sqrt(sampler.gyration2());
        const bool writeEnergies = nstener > 0 && step % nstener == 0;
        const bool writeX        = nstxout > 0 && step % nstxout == 0;
        if (haveDistance || writeX)
        {
            sampler.centroids(&xc);
        }
        if (haveDistance)
        {
            sample.distance = wall.distance(xc);
            if (0 == step)
            {
                distance0 = sample.distance;
            }
            if (!wall.active() && !results.movedApart && sample.distance > distance0 + apartDistance)
            {
                results.movedApart = true;
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("The compounds moved apart: at step %d their distance is %.3f nm, it was %.3f nm at the start. In vacuum they will not come back. A wall (-rwall) keeps them together.",
                                                  step, sample.distance, distance0));
            }
        }
        if (step >= nequil)
        {
            epot.add(sample.epot);
            ekin.add(sample.ekin);
            etot.add(sample.epot + sample.ekin);
            const double kprim = sampler.primitiveKineticEnergy();
            ekinPrimitive.add(kprim);
            if (!haveEnergy0)
            {
                energy0     = sample.epot + sample.ekin;
                haveEnergy0 = true;
            }
            fluctuation.add(sample.epot + sample.ekin - energy0, sample.epot + kprim - energy0);
            gyration.add(sample.gyration);
            if (haveDistance)
            {
                distance.add(sample.distance);
                nwall += (sample.ewall > 0) ? 1 : 0;
            }
            if (2 == wall.nfragments())
            {
                size_t bin = std::min(static_cast<size_t>(sample.distance/results.binWidth), maxBins - 1);
                if (bin >= results.histogram.size())
                {
                    results.histogram.resize(bin + 1, 0.0);
                }
                results.histogram[bin] += 1;
            }
        }
        if (writeEnergies)
        {
            const auto &eterms = sampler.energyTerms();
            if (step >= nequil)
            {
                for(const auto &t : eterms)
                {
                    auto tt = terms.find(t.first);
                    if (terms.end() == tt)
                    {
                        tt = terms.insert({ t.first, BlockAverager(nterm) }).first;
                    }
                    tt->second.add(t.second);
                }
            }
            if (exvg)
            {
                if (!legendDone)
                {
                    const auto                names = energyNames(eterms, wall.active(), haveDistance);
                    std::vector<const char *> cnames;
                    for(const auto &name : names)
                    {
                        cnames.push_back(name.c_str());
                    }
                    xvgr_legend(exvg, cnames.size(), cnames.data(), files.oenv);
                    legendDone = true;
                }
                fprintf(exvg, "%10g", step*deltat);
                for(const auto &t : eterms)
                {
                    fprintf(exvg, "  %.5f", t.second);
                }
                if (wall.active())
                {
                    fprintf(exvg, "  %.5f", sample.ewall);
                }
                fprintf(exvg, " %.5f %.5f %.5f %.6f", sample.ekin, sample.epot + sample.ekin,
                        sample.temperature, sample.gyration);
                if (haveDistance)
                {
                    fprintf(exvg, " %.6f", sample.distance);
                }
                fprintf(exvg, "\n");
            }
            if (progress)
            {
                printProgress(progress, sample, nsteps, deltat, nbeads, eterms, wall.active(), haveDistance);
            }
            if (observer && observer->energies)
            {
                keepGoing = observer->energies(step, step*deltat, sample.temperature, sample.ekin, eterms) && keepGoing;
            }
            flush(false);
        }
        // The centroids of the same time, numbering the frames from 1
        if (writeX)
        {
            for(size_t i = 0; i < atoms.size(); i++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    xr[i][m] = xc[i][m];
                }
            }
            if (traj)
            {
                pdbWriter(traj, title, atoms, xr, top->residueNames(), epbcNONE,
                          box, ' ', step/nstxout + 1, trajIndex, {}, true);
            }
            if (observer && observer->coordinates)
            {
                keepGoing = observer->coordinates(step, step*deltat, xr) && keepGoing;
            }
            flush(false);
        }
        if (!keepGoing)
        {
            // What was due for this step is out, and no further step is taken
            break;
        }
        // The energy terms are computed for the steps in which they are written
        sampler.run(msghandler, 1, nstener > 0 && (step + 1) % nstener == 0);
    }
    flush(true);
    if (traj)
    {
        gmx_ffclose(traj);
    }
    if (exvg)
    {
        xvgrclose(exvg);
    }
    if (progress)
    {
        gmx_ffclose(progress);
    }
    results.nsteps   = step;
    results.finished = keepGoing;
    results.epot     = epot.result();
    results.ekin     = ekin.result();
    results.etot     = etot.result();
    results.ekinPrimitive = ekinPrimitive.result();
    {
        // The heat capacity from the fluctuations of the energy
        const double kT   = BOLTZ*settings.temperature;
        const auto   cov  = fluctuation.result();
        results.heatCapacity.n       = cov.n;
        results.heatCapacity.average = BOLTZ*(cov.average/(kT*kT) + 0.5*sampler.ncoordinates());
        results.heatCapacity.error   = BOLTZ*cov.error/(kT*kT);
    }
    results.gyration = gyration.result();
    results.distance = distance.result();
    for(const auto &t : terms)
    {
        results.terms.insert({ t.first, t.second.result() });
    }
    if (results.wall && results.distance.n > 0)
    {
        results.wallFraction = static_cast<double>(nwall)/results.distance.n;
        if (results.wallFraction > wallFractionLimit)
        {
            msghandler->msg(ACTStatus::Warning,
                            gmx::formatString("The compounds were at the wall in %.1f %% of the steps. They are not bound at %g K, and the energies depend on the radius of the wall.",
                                              100*results.wallFraction, results.temperature));
        }
    }
    if (files.histogram && !results.histogram.empty())
    {
        writeHistogram(files.histogram, files.oenv, results);
    }
    msghandler->write(keepGoing ? "Simulation finished correctly." : "Simulation ended on request.");
    return results;
}

void reportPathIntegral(MsgHandler                *msghandler,
                        const std::string         &title,
                        const PathIntegralResults &results,
                        JsonTree                  *jtree)
{
    if (0 == results.epot.n)
    {
        msghandler->write(gmx::formatString("\nNo averages for %s: the simulation ended within the %d steps of equilibration.",
                                            title.c_str(), results.nequil));
        return;
    }
    msghandler->write(gmx::formatString("\nAverages for %s at %g K with %d bead%s, over %zu steps after %d steps of equilibration.",
                                        title.c_str(), results.temperature, results.nbeads,
                                        results.nbeads == 1 ? "" : "s", results.epot.n, results.nequil));
    msghandler->write("The standard errors are from eight blocks, they are too small if the simulation is short.");
    JsonTree tree(title);
    tree.addObject("Beads", results.nbeads);
    tree.addValueUnit("Temperature", results.temperature, "K");
    tree.addObject("Steps", static_cast<int>(results.epot.n));
    const auto add = [&](const std::string &name, const AverageError &ae, const char *unit, int decimals)
    {
        msghandler->write(gmx::formatString("  %-22s %14.*f +/- %.*f %s", name.c_str(), decimals, ae.average,
                                            decimals, ae.error, unit));
        JsonTree sub(name);
        sub.addObject("Average", ae.average);
        sub.addObject("Standard error", ae.error);
        sub.addObject("Unit", std::string(unit));
        tree.addObject(sub);
    };
    add("Potential energy", results.epot, "kJ/mol", 4);
    add("Kinetic energy", results.ekin, "kJ/mol", 4);
    add("Total energy", results.etot, "kJ/mol", 4);
    add("Spread of the beads", results.gyration, "nm", 5);
    if (results.nfragments > 1)
    {
        add(results.nfragments == 2 ? "Distance" : "Largest distance", results.distance, "nm", 4);
    }
    if (results.wall)
    {
        msghandler->write(gmx::formatString("  The compounds were at the wall in %.2f %% of the steps.", 100*results.wallFraction));
        tree.addObject("Fraction at the wall", results.wallFraction);
    }
    if (results.movedApart)
    {
        tree.addBool("Moved apart", true);
    }
    if (!results.terms.empty())
    {
        msghandler->write(gmx::formatString("Energy terms, averages over the beads and over the %zu steps of which the energies were written:",
                                            results.terms.begin()->second.n));
        JsonTree termTree("Energy terms");
        for(const auto &t : results.terms)
        {
            const auto name = interactionTypeToString(t.first);
            msghandler->write(gmx::formatString("  %-22s %14.4f +/- %.4f kJ/mol", name.c_str(),
                                                t.second.average, t.second.error));
            JsonTree sub(name);
            sub.addObject("Average", t.second.average);
            sub.addObject("Standard error", t.second.error);
            termTree.addObject(sub);
        }
        tree.addObject(termTree);
    }
    if (jtree)
    {
        jtree->addObject(tree);
    }
}

AverageError bindingEnergy(const PathIntegralResults              &cluster,
                           const std::vector<PathIntegralResults> &fragments,
                           AverageError                           *deltaH)
{
    AverageError de;
    de.average  = cluster.etot.average;
    de.n        = cluster.etot.n;
    double var  = cluster.etot.error*cluster.etot.error;
    for(const auto &f : fragments)
    {
        de.average -= f.etot.average;
        var        += f.etot.error*f.etot.error;
    }
    de.error = std::sqrt(var);
    if (deltaH)
    {
        // Each fragment that binds takes k_B T of pV of an ideal gas away
        *deltaH = de;
        if (!fragments.empty())
        {
            deltaH->average -= (static_cast<double>(fragments.size()) - 1)*BOLTZ*cluster.temperature;
        }
    }
    return de;
}

} // namespace alexandria
