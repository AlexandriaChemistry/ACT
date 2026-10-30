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
#include "secondvirial.h"

#include <cctype>
#include <cmath>
#include <cstdlib>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/train_utility.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/import/compound_reader.h"
#include "act/import/import.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/properties/b2data.h"
#include "act/properties/princ.h"
#include "act/topology/actmol.h"
#include "act/utility/jsontree.h"
#include "act/utility/timings.h"
#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

//! Map b2Type to string
const std::map<b2Type, std::string> b2Type2str = {
    { b2Type::Classical, "B2(Classical)" },
    { b2Type::Force, "B2(Force)" },
    { b2Type::Torque1, "B2(Torque1)" },
    { b2Type::Torque2, "B2(Torque2)" },
    { b2Type::Total, "B2(Total)" }
};                    

const std::string &b2TypeToString(b2Type b2t)
{
    return b2Type2str.find(b2t)->second;
}

void forceFieldSummary(JsonTree      *jtree,
                       const ForceField *pd)
{
    jtree->addObject(JsonTree("Force field file", pd->filename()));
    jtree->addObject(JsonTree("Created", pd->timeStamp()));
    jtree->addObject(JsonTree("Checksum", pd->checkSum()));
    jtree->addObject(JsonTree("Polarizable", yesno_names[pd->polarizable()]));
    jtree->addObject(JsonTree("Charge generation", 
                              chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str()));
    int nexclvdw;
    if (!ffOption(*pd, InteractionType::VDW, 
                  "nexcl", &nexclvdw))
    {
        nexclvdw = 0;
    }
    jtree->addObject(JsonTree("# vanderwaals exclusions", gmx_itoa(nexclvdw)));
    int nexclqq;
    if (!ffOption(*pd, InteractionType::ELECTROSTATICS, 
                  "nexcl", &nexclqq))
    {
        nexclqq = 0;
    }
    jtree->addObject(JsonTree("# coulomb exclusions", gmx_itoa(nexclqq)));
    double epsilonr;
    if (!ffOption(*pd, InteractionType::ELECTROSTATICS, 
                  "epsilonr", &epsilonr))
    {
        epsilonr = 1;
    }
    jtree->addObject(JsonTree("Relative dielectric constant epsilon_r",
                              gmx_ftoa(epsilonr)));
    jtree->addObject(JsonTree("# particle types", gmx_itoa(pd->getNatypes())));
    
    JsonTree ftree("InteractionTypes");
    for(const auto &fs : pd->forcesConst())
    {
        auto itype = fs.first;
        auto &ffpl = fs.second;
        if (!ffpl.parametersConst().empty())
        {
            JsonTree fftree(interactionTypeToString(itype));
            auto function = potentialToString(ffpl.potential());
            if (!function.empty())
            {
                fftree.addObject(JsonTree("Force function", function));
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

void ReRunner::addOptions(std::vector<t_pargs>      *pargs,
                          std::vector<t_filenm>     *filenm,
                          std::vector<const char *> *desc)
{
    std::vector<t_pargs> pa = {
        { "-T1",     FALSE, etREAL, {&T1_},
          "Starting temperature for second virial calculations." },
        { "-T2",     FALSE, etREAL, {&T2_},
          "Starting temperature for second virial calculations." },
        { "-dT",     FALSE, etREAL, {&deltaT_},
          "Temperature increment for calculation of second virial." },
        { "-totalOnly", FALSE, etBOOL, {&totalOnly_},
          "Plot the total B2(T) only and not the components" }
    };
    pargs->push_back(pa[0]);
    if (computeB2_)
    {
        for(size_t i = 1; i < pa.size(); i++)
        {
            pargs->push_back(pa[i]);
        }
    
        std::vector<t_filenm> fnm = {
            { efXVG, "-eh", "mayer",      ffOPTWR },
            { efXVG, "-b2", "B2T",        ffOPTWR }
        };
        for(const auto &fn : fnm)
        {
            filenm->push_back(fn);
        }
    }
    gendimers_.addOptions(pargs, filenm, desc);
}

const std::vector<double> &ReRunner::temperatures()
{
    if (Temperatures_.empty())
    {
        Temperatures_.push_back(T1_);
        if (deltaT_ > 0)
        {
            double T = T1_+deltaT_;
            while (T <= T2_)
            {
                Temperatures_.push_back(T);
                T += deltaT_;
            }
        }
        else
        {
            Temperatures_.push_back(T2_);
        }
    }
    return Temperatures_;
}


void ReRunner::plotB2temp(const char *b2file)
{
    if (!b2file || strlen(b2file) == 0)
    {
        return;
    }
    auto T = temperatures();
    if (b2t_.empty() || b2t_[b2Type::Total].size() != T.size())
    {
        fprintf(stderr, "Internal inconsistency. There are %zu temperatures and %zu B2 values.\n",
                T.size(), b2t_.size());
        return;
    }
    FILE *b2p = xvgropen(b2file, "Second virial coefficient",
                         "Temperature (K)", "B2(T) cm^3/mol", oenv_);
    if (totalOnly_)
    {
        for(size_t ii = 0; ii < T.size(); ii++)
        {
            fprintf(b2p, "%10g  %10g\n", T[ii], b2t_[b2Type::Total][ii]);
        }
    }
    else
    {
        std::vector<std::string> legend = {
            "Total", "Classical", "Force", "Torque1", "Torque2"
        };
        xvgrLegend(b2p, legend, oenv_);
        for(size_t ii = 0; ii < T.size(); ii++)
        {
            fprintf(b2p, "%10g  %10g  %10g  %10g  %10g  %10g\n", T[ii], 
                    b2t_[b2Type::Total][ii], b2t_[b2Type::Classical][ii],
                    b2t_[b2Type::Force][ii], b2t_[b2Type::Torque1][ii],
                    b2t_[b2Type::Torque2][ii]);
        }
    }
    xvgrclose(b2p);
}

/*! \brief Make the coordinates of all particles from those of the atoms
 *
 * If there are as many input coordinates as particles, they are used as
 * they are. Otherwise the input has the atoms only, and the other particles
 * get the coordinates of their parent atom.
 * \param[in] atoms The particles
 * \param[in] x     The input coordinates
 * \return the coordinates of all particles
 */
static std::vector<gmx::RVec> allParticleCoords(const std::vector<ActAtom>   &atoms,
                                                const std::vector<gmx::RVec> &x)
{
    if (x.size() == atoms.size())
    {
        return x;
    }
    std::vector<gmx::RVec> coords;
    size_t index = 0;
    for(size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].pType() == ActParticle::Atom)
        {
            if (index >= x.size())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Number of (generated) coordinates in trajectory (%zu) does not match molecule file (%zu)", x.size(), atoms.size()).c_str()));
            }
            coords.push_back(x[index++]);
        }
        else if (!atoms[i].parents().empty() && atoms[i].parents()[0] < static_cast<int>(i))
        {
            coords.push_back(coords[atoms[i].parents()[0]]);
        }
        else
        {
            // Will be constructed, e.g. a virtual site
            coords.push_back({ 0, 0, 0 });
        }
    }
    if (index != x.size())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Number of (generated) coordinates in trajectory (%zu) does not match molecule file (%zu)", x.size(), atoms.size()).c_str()));
    }
    return coords;
}

void ReRunner::rerun(MsgHandler       *msghandler,
                     const ForceField *pd,
                     ACTMol           *actmol,
                     bool              verbose)
{
    std::vector<std::vector<gmx::RVec>> dimers;
    gendimers_.read(&dimers);
    msghandler->write(gmx::formatString("Doing energy calculation for %zu structures from %s\n",
                                        dimers.size(), gendimers_.trajname()));

    if (verbose && msghandler)
    {
        msghandler->writeDebug(memory_usage());
    }
    gmx_stats edist;
    // Coordinates of all particles of each structure
    const auto &atoms = actmol->atomsConst();
    std::vector<std::vector<gmx::RVec>> allCoords;
    for (const auto &x : dimers)
    {
        allCoords.push_back(allParticleCoords(atoms, x));
    }
    // Compute the energies of all structures at once
    std::vector<std::map<InteractionType, double>> allEnergies;
    {
        ConformationSet set;
        set.top = actmol->topology();
        set.setCoordinates(allCoords);
        forceComp_->computeTotal(msghandler, pd, &set, ShellOptions(), &allEnergies);
        for(size_t c = 0; c < allCoords.size(); c++)
        {
            allCoords[c] = set.confs[c].x;
        }
    }
    for (size_t idim = 0; idim < allCoords.size(); idim++)
    {
        const auto &coords   = allCoords[idim];
        auto       &energies = allEnergies[idim];
        std::string out = gmx::formatString("%5zu", idim);
        if (eInter_)
        {
            auto atomStart  = actmol->fragmentHandler()->atomStart();
            std::vector<gmx::RVec> com        = { { 0, 0, 0 }, { 0, 0, 0 } };
            std::vector<double>    mtot       = { 0, 0 };
            const auto &tops  = actmol->fragmentHandler()->topologies();
            for(int kk = 0; kk < 2; kk++)
            {
                auto natom = tops[kk].atoms().size();
                for(size_t i = atomStart[kk]; i < atomStart[kk]+natom; i++)
                {
                    gmx::RVec mr1;
                    auto      mi = atoms[i].mass();
                    svmul(mi, coords[i], mr1);
                    mtot[kk] += mi;
                    // Compute center of mass of compound kk
                    rvec_inc(com[kk], mr1);
                }
                GMX_RELEASE_ASSERT(mtot[kk] > 0, "Zero mass");
                for(size_t m = 0; m < DIM; m++)
                {
                    // Normalize
                    com[kk][m] /= mtot[kk];
                }
            }
            gmx::RVec dcom;
            rvec_sub(com[0], com[1], dcom);
            double rcom = norm(dcom);
            out += gmx::formatString(" r %g", rcom);
            for (auto &EE: energies)
            {
                out += gmx::formatString(" %s %g", interactionTypeToString(EE.first).c_str(), EE.second);
            }
            edist.add_point(rcom, energies[InteractionType::EPOT], 0, 0);
        }
        else
        {
            for(const auto &ee : energies)
            {
                out += gmx::formatString("  %s %8g", interactionTypeToString(ee.first).c_str(), ee.second);
            }
        }
        if (msghandler->info())
        {
            msghandler->write(out);
        }
    }
    if (verbose && msghandler)
    {
        msghandler->writeDebug(memory_usage());
    }
}

void ReRunner::runB2(CommunicationRecord         *cr,
                     MsgHandler                  *msghandler,
                     const ForceField            *pd,
                     ACTMol                      *actmol,
                     int                          maxdimer,
                     const std::vector<t_filenm> &fnm)
{
    // Compute the relative masses
    std::vector<double> masses = {
        actmol->fragmentHandler()->topologies()[0].mass(),
        actmol->fragmentHandler()->topologies()[1].mass()
    };
    std::vector<dimer_coords> dimers;
    // The orientations are computed in rounds of 1024, 2048, 4096, ...
    // until maxdimer, each round adding to the data of the previous ones,
    // such that the convergence with the number of orientations can be
    // followed. All orientations together are the same as in one round.
    // A grid of orientations is only complete in one round.
    std::vector<int> rounds;
    if (RandAlg::Gauss == gendimers_.randAlg())
    {
        maxdimer = gendimers_.numberOfOrientations(maxdimer);
    }
    else
    {
        for(int n = 1024; n < maxdimer; n *= 2)
        {
            rounds.push_back(n);
        }
    }
    rounds.push_back(maxdimer);
    // The master generates the random numbers of all orientations and hands
    // them out in each round. Different nodes get different random numbers.
    // For each dimer, 2 x 3 = 6 random numbers are needed.
    std::vector<std::vector<double>> allRand;
    if (cr->isMaster())
    {
        static const size_t timer = Timings::instance().index("Random numbers");
        ScopedTimer         scoped(timer);
        gendimers_.generateRandomNumbers(maxdimer);
        allRand = gendimers_.allRandom();
        gendimers_.clearAllRandom();
    }
    /* Distribute the orientations from start to end over the nodes, the
     * master takes the first ones, including the rest of the division.
     * The random numbers are added to those of the previous rounds.
     * Returns the number of orientations of this node.
     */
    auto distribute = [&](int start, int end)
    {
        static const size_t timer = Timings::instance().index("Distributing orientations");
        ScopedTimer         scoped(timer);
        int nper  = (end - start) / cr->size();
        int nrest = (end - start) % cr->size();
        int nmine = 0;
        if (cr->isMaster())
        {
            for(int dst = 1; dst < cr->size(); dst++)
            {
                cr->send(dst, nper);
                int nn = start + nrest + dst*nper;
                for(int j = 0; j < nper; j++)
                {
                    cr->send(dst, allRand[nn + j]);
                }
            }
            nmine = nper + nrest;
            for(int j = 0; j < nmine; j++)
            {
                gendimers_.addRandomNumbers(allRand[start + j]);
            }
            msghandler->write(gmx::formatString("Will generate %d dimers on helpers and %d on master.", nper, nmine));
        }
        else
        {
            int src = cr->superior();
            cr->recv(src, &nmine);
            for(int i = 0; i < nmine; i++)
            {
                std::vector<double> q;
                cr->recv(src, &q);
                gendimers_.addRandomNumbers(q);
            }
        }
        return nmine;
    };
    // Will be used to obtain a seed for the random number engine for bootstrapping
    std::random_device                 bsRand;
    // Standard mersenne_twister_engine seeded with rd()
    std::mt19937                       bsGen(bsRand());
    std::uniform_int_distribution<int> bsDistr(0, gendimers_.ndist());

    // Temperature array.
    auto   Temperature = temperatures();
    for(const auto &b2b : b2Type2str)
    {
        b2t_.insert({b2b.first, {}});
        b2t_[b2b.first].resize(Temperature.size());
    }
    // Temporary arrays for weighted properties.
    B2Data b2data(gendimers_.ndist(), Temperature);
    
    gmx::RVec rvnul = { 0, 0, 0 };
    // Moments of inertia, averaged over all dimers.
    std::vector<gmx::RVec> inertia    = { rvnul, rvnul };
    std::vector<gmx::RVec> inertiaSum = { rvnul, rvnul };
    size_t                 nInertia   = 0;
    // The other parts are timed where they are done
    auto               &timings      = Timings::instance();
    static const size_t prepareTimer = timings.index("Dimer generation");
    static const size_t energyTimer  = timings.index("Interaction energies");
    static const size_t torqueTimer  = timings.index("Forces and torques");
    static const size_t addTimer     = timings.index("B2 accumulation");
    static const size_t sumTimer     = timings.index("B2 aggregation");
    static const size_t b2Timer      = timings.index("B2 integration");
    static const size_t outputTimer  = timings.index("B2 output");
    // First time around before generating
    timings.start(prepareTimer);
    gendimers_.prepare(msghandler, pd, actmol);
    timings.stop(prepareTimer);
    // Collect all generated dimer coordinates for -ox output (only when requested)
    std::vector<std::vector<gmx::RVec>> allCoords;
    // The data of all nodes of all rounds so far, on the master
    B2Data b2total(gendimers_.ndist(), Temperature);
    // The total B2 of each round, for the convergence table
    std::vector<std::vector<double>> b2rounds;
    // The molecules of the dimers, the same in all dimers
    const auto &atoms     = actmol->atomsConst();
    const auto &atomStart = actmol->fragmentHandler()->atomStart();
    const auto &tops      = actmol->fragmentHandler()->topologies();
    // Whether the monomers are flexible, otherwise they are rigid
    bool        flexible  = gendimers_.flexible();
    size_t                 natom[2];
    std::vector<int>       molIndex[2];
    std::vector<real>      molMass[2];
    // Coordinates relative to the center of mass
    std::vector<gmx::RVec> xCom[2];
    // Principal axes and moments of inertia of the current dimer
    matrix                 trans[2];
    gmx::RVec              inertia1[2];
    for(int kk = 0; kk < 2; kk++)
    {
        natom[kk] = tops[kk].atoms().size();
        for(size_t i = atomStart[kk]; i < atomStart[kk]+natom[kk]; i++)
        {
            // Store atom index relative to molecule start
            molIndex[kk].push_back(i-atomStart[kk]);
            molMass[kk].push_back(atoms[i].mass());
        }
        xCom[kk].resize(natom[kk]);
    }
    // Squared forces and torques on the molecules, energies, distances,
    // distance indices and weights of the dimers of an orientation
    std::vector<double>    force2Mol[2];
    std::vector<gmx::RVec> torqueMol[2];
    std::vector<double>    dimEnergy, dimDist, dimWeight;
    std::vector<size_t>    dimIndex;
    int ndone = 0;
    for(auto nround : rounds)
    {
        int ndimer = distribute(ndone, nround);
        //! Loop over my dimers, completely independently
        for(int iRotation = 0; iRotation < ndimer; iRotation++)
        {
            // Generate a new set of dimers for all distances
            if (!gendimers_.hasTrajectory())
            {
                ScopedTimer scoped(prepareTimer);
                dimers = gendimers_.generateDimers(msghandler);
                if (gendimers_.outCoords())
                {
                    for (auto &d : dimers)
                    {
                        allCoords.push_back(d.coords);
                    }
                }
            }
            // Coordinates of all particles of the dimers
            timings.start(energyTimer);
            std::vector<std::vector<gmx::RVec>> allCoordsShells;
            for (const auto &d : dimers)
            {
                allCoordsShells.push_back(allParticleCoords(atoms, d.coords));
            }
            // Compute the interaction energies and forces of all dimers at once
            std::vector<std::map<InteractionType, double>> allEinter;
            std::vector<double>                            allEpot;
            std::vector<std::vector<gmx::RVec>>            allForces;
            forceComp_->computeInteraction(msghandler, pd, actmol, &allCoordsShells, true,
                                           &allEinter, &allEpot, &allForces);
            timings.stop(energyTimer);
            // Loop over dimers
            timings.start(torqueTimer);
            size_t ndim = dimers.size();
            for(int kk = 0; kk < 2; kk++)
            {
                force2Mol[kk].resize(ndim);
                torqueMol[kk].resize(ndim);
            }
            dimIndex.resize(ndim);
            dimDist.resize(ndim);
            dimEnergy.resize(ndim);
            dimWeight.resize(ndim);
            for (size_t idim = 0; idim < ndim; idim++)
            {
                const auto &coords = allCoordsShells[idim];
                const auto &forces = allForces[idim];
                auto       &einter = allEinter[idim];
                auto        epot   = allEpot[idim];

                gmx::RVec f[2]         = { rvnul, rvnul };
                gmx::RVec com[2]       = { rvnul, rvnul };
                double    mtot[2]      = { 0, 0 };
                gmx::RVec torque[2]    = { rvnul, rvnul };
                gmx::RVec torqueRot[2] = { rvnul, rvnul };
                for(int kk = 0; kk < 2; kk++)
                {
                    for(size_t i = atomStart[kk]; i < atomStart[kk]+natom[kk]; i++)
                    {
                        gmx::RVec mr1;
                        auto      mi = atoms[i].mass();
                        svmul(mi, coords[i], mr1);
                        mtot[kk] += mi;
                        // Compute center of mass of compound kk
                        rvec_inc(com[kk], mr1);
                        // Compute total force on compound kk
                        rvec_inc(f[kk], forces[i]);
                    }
                    GMX_RELEASE_ASSERT(mtot[kk] > 0, "Zero mass");
                    for(size_t m = 0; m < DIM; m++)
                    {
                        // Normalize
                        com[kk][m] /= mtot[kk];
                    }
                    // Compute the coordinates relative to the center of mass
                    auto &xc = xCom[kk];
                    for(size_t i = atomStart[kk]; i < atomStart[kk]+natom[kk]; i++)
                    {
                        rvec_sub(coords[i], com[kk], xc[i-atomStart[kk]]);
                    }
                    // Compute moments of inertia and transformation matrix.
                    // Rigid molecules are only translated between the dimers
                    // of an orientation, so their principal axes are the same.
                    if (flexible || 0 == idim)
                    {
                        clear_rvec(inertia1[kk]);
                        clear_mat(trans[kk]);
                        principal_comp(molIndex[kk], molMass[kk], xc, &trans[kk], &inertia1[kk]);
                    }

                    // Move to inertial frame (only well-defined for
                    // rigid molecules).
                    // The trans matrix should convert that coordinate to the inertial frame,
                    // but what about the force on the atoms? It likely has to be rotated in the
                    // same manner. After that, the torque can be computed.
                    // Before doing the rotations, the sum of the torque vectors is zero.
                    // Since the rotations are different for both molecules, this does not
                    // hold after the rotations are done.
                    //! \todo: write out the math.
                    for(size_t i = atomStart[kk]; i < atomStart[kk]+natom[kk]; i++)
                    {
                        gmx::RVec ri, fi, ti;
                        cprod(xc[i-atomStart[kk]], forces[i], ti);
                        rvec_inc(torque[kk], ti);
                        // Rotate coordinates
                        mvmul(trans[kk], xc[i-atomStart[kk]], ri);
                        // Rotate force vector
                        mvmul(trans[kk], forces[i], fi);
                        // Compute torque on this atom
                        cprod(ri, fi, ti);
                        // Update total torque
                        rvec_inc(torqueRot[kk], ti);
                    }
                    rvec_inc(inertiaSum[kk], inertia1[kk]);
                    force2Mol[kk][idim] = iprod(f[kk], f[kk]);
                    torqueMol[kk][idim] = torqueRot[kk];
                }
                nInertia += 1;
                // Check consistency from the rotation algorithm.
                rvec dx;
                rvec_sub(com[0], com[1], dx);
                double rcom = norm(dx);
                if (std::abs(rcom - dimers[idim].dist) > 1e-4)
                {
                    msghandler->fatal(gmx::formatString(" rcom calculated %g from input %g",
                                                        rcom, dimers[idim].dist));
                }
                rcom = dimers[idim].dist;
                if (msghandler->verbose())
                {
                    std::string out = gmx::formatString(" r %g epot %g", rcom, epot);
                    for (auto &EE: einter)
                    {
                        out += gmx::formatString(" %s %g", interactionTypeToString(EE.first).c_str(), EE.second);
                    }
                    out += gmx::formatString(" Force %g %g %g Torque[0] %g %g %g Torque[1] %g %g %g Rotated Torque[0] %g %g %g Rotated Torque[1] %g %g %g",
                                             f[0][XX], f[0][YY], f[0][ZZ],
                                             torque[0][XX], torque[0][YY], torque[0][ZZ],
                                             torque[1][XX], torque[1][YY], torque[1][ZZ],
                                             torqueRot[0][XX], torqueRot[0][YY], torqueRot[0][ZZ],
                                             torqueRot[1][XX], torqueRot[1][YY], torqueRot[1][ZZ]);
                    msghandler->write(out);
                }
                dimIndex[idim]  = dimers[idim].index;
                dimDist[idim]   = rcom;
                dimEnergy[idim] = flexible ? epot : einter[InteractionType::EPOT];
                dimWeight[idim] = dimers[idim].weight;
                if (msghandler->verbose())
                {
                    msghandler->write(gmx::formatString("epot %g einter[InteractionType::EPOT] %g diff %g", epot, einter[InteractionType::EPOT], epot-einter[InteractionType::EPOT]));
                }
            }
            if (msghandler->verbose())
            {
                msghandler->write("");
            }
            timings.stop(torqueTimer);
            // Now we have computed energies, forces and torques for all
            // dimers of this orientation. Time to update the intermediate structure
            ScopedTimer scoped(addTimer);
            b2data.addOrientation(dimIndex, dimDist, dimEnergy, force2Mol, torqueMol, dimWeight);
        }
        ndone = nround;
        for(int kk = 0; kk < 2 && nInertia > 0; kk++)
        {
            for(int m = 0; m < DIM; m++)
            {
                inertia[kk][m] = inertiaSum[kk][m]/nInertia;
            }
        }
        // Sum b2data over processors, keeping the sums of this node
        // for the next rounds
        timings.start(sumTimer);
        b2total = b2data;
        b2total.aggregate(cr);
        timings.stop(sumTimer);

        if (cr->isMaster())
        {
            ScopedTimer scoped(b2Timer);
            msghandler->write(gmx::formatString("\nSecond virial coefficient (cm^3/mol) from %d dimer orientations", nround));
            auto header = gmx::formatString("Temperature");
            for(const auto &b2b : b2Type2str)
            {
                header += gmx::formatString("  %16s", b2b.second.c_str());
            }
            msghandler->write(header);
            std::vector<double> btotal;
            for(size_t iTemp = 0; iTemp < Temperature.size(); iTemp++)
            {
                auto T      = Temperature[iTemp];
                double beta = 1.0/(BOLTZ*T);
                // Now compute the components
                double Bclass, BqmForce, BqmTorque1, BqmTorque2;
                b2total.integrate(msghandler, iTemp, beta, masses, inertia,
                                  &Bclass, &BqmForce, &BqmTorque1, &BqmTorque2);
                if (msghandler->ok())
                {
                    // Conversion to regular units cm^3/mol.
                    double fac  = AVOGADRO*1e-21;
                    //! \todo: Fix the torque contribution
                    double bqt  = (BqmTorque1+BqmTorque2)*0.5;
                    double Btot = (Bclass + BqmForce + bqt)*fac;
                    // Add to bootstrapping statistics
                    b2t_[b2Type::Classical][iTemp] = Bclass*fac;
                    b2t_[b2Type::Force][iTemp]     = BqmForce*fac;
                    b2t_[b2Type::Torque1][iTemp]   = BqmTorque1*fac;
                    b2t_[b2Type::Torque2][iTemp]   = BqmTorque2*fac;
                    b2t_[b2Type::Total][iTemp]     = Btot;

                    std::string out = gmx::formatString("%11.2f", T);
                    for(const auto &b2b : b2Type2str)
                    {
                        out += gmx::formatString("  %16.2f", b2t_[b2b.first][iTemp]);
                    }
                    msghandler->write(out);
                }
                btotal.push_back(b2t_[b2Type::Total][iTemp]);
            }
            b2rounds.push_back(btotal);
        }
    }

    ScopedTimer scoped(outputTimer);
    // Gather generated dimer coordinates from all nodes and write to -ox file
    if (gendimers_.outCoords() && !gendimers_.hasTrajectory())
    {
        if (!cr->isMaster())
        {
            int dest = cr->superior();
            int n    = static_cast<int>(allCoords.size());
            cr->send(dest, n);
            for (auto &d : allCoords)
            {
                cr->send(dest, d);
            }
        }
        else
        {
            for (int src = 1; src < cr->size(); src++)
            {
                int n = 0;
                cr->recv(src, &n);
                for (int i = 0; i < n; i++)
                {
                    std::vector<gmx::RVec> d;
                    cr->recv(src, &d);
                    allCoords.push_back(std::move(d));
                }
            }
            gendimers_.writeCoords(actmol, allCoords);
        }
    }

    if (cr->isMaster())
    {
        // Convergence of the total B2 with the number of orientations
        msghandler->write("\nConvergence of the total second virial coefficient (cm^3/mol) with the number of dimer orientations");
        auto header = gmx::formatString("%10s", "Dimers");
        for(auto T : Temperature)
        {
            header += gmx::formatString("  %10.2f", T);
        }
        msghandler->write(header);
        for(size_t r = 0; r < b2rounds.size(); r++)
        {
            auto out = gmx::formatString("%10d", rounds[r]);
            for(auto b : b2rounds[r])
            {
                out += gmx::formatString("  %10.2f", b);
            }
            msghandler->write(out);
        }
        if (!fnm.empty())
        {
            // Store Mayer functions if requested.
            b2total.plotMayer(opt2fn_null("-eh", fnm.size(), fnm.data()), oenv_);
            // Print B2(T) if requested
            plotB2temp(opt2fn("-b2", fnm.size(), fnm.data()));
        }
    }
}

/*! \brief Tool to compute second virial
 * \param[in] argc Number of arguments on the cmd line
 * \param[in] argv The actual arguments
 * \return 0 if ok, 1 otherwise
 */
int b2(int argc, char *argv[])
{
    // Time the phases of the program, see printTimings.
    PhaseTimer phase("Startup");
    std::vector<const char *> desc = {
        "alexandria b2 will generate conformations of dimers. The",
        "corresponding charges (molprop) file, used for generating the topology",
        "needs to contain information about",
        "the compounds in the dimer. Based on dimer energies the second",
        "virial coefficient will be estimated. Mayer curves can optionally",
        "be plotted and the second virial can be plotted as a function",
        "of temperature.[PAR]",
        "To run the program you need to provide a pdb file with two monomeric",
        "compounds that have support in the force field and that are present",
        "in the charges file."
    };

    std::vector<t_filenm>     fnm = {
        { efXML, "-ff",      "aff",     ffREAD  }
    };
    gmx_output_env_t         *oenv;
    static char              *molnm            = (char *)"";
    double                    shellToler       = 1e-6;
    int                       maxdimers        = 1024;
    int                       shellMaxIter     = 25;
    bool                      json             = false;
    std::vector<t_pargs>      pa = {
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule." },
        { "-maxdimer", FALSE, etINT, {&maxdimers},
          "Number of dimer orientations to generate if you do not provide a trajectory. For each of these a distance scan will be performed." },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)." },
        { "-shellMaxIter", FALSE, etINT, {&shellMaxIter},
          "Max number of iterations for minimizing shell positions" },
        { "-json", FALSE, etBOOL, {&json},
          "Print part of the output in json format to b2.json." }
    };
    CommunicationRecord cr;
    cr.init(cr.size());
    ReRunner            rerun(true);
    rerun.addOptions(&pa, &fnm, &desc);
    CompoundReader      compR;
    compR.addOptions(&pa, &fnm, &desc);
    MsgHandler          msghandler;
    msghandler.addOptions(&pa, &fnm, "b2");
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        status = 1;
        return status;
    }
    msghandler.optionsFinished(fnm, &cr);
    if (cr.isMaster())
    {
        print_header(msghandler.tw(), pa, fnm);
    }
    rerun.gendimers()->finishOptions(fnm);
    compR.optionsFinished(&msghandler, fnm);
    if (!msghandler.ok())
    {
        return 1;
    }
    ForceField        pd;
    try
    {
        readForceField(opt2fn("-ff", fnm.size(), fnm.data()), &pd, &msghandler);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    (void) pd.verifyCheckSum(&msghandler);
    double epsilonr = 1;
    if (!ffOption(pd, InteractionType::ELECTROSTATICS, 
                  "epsilonr", &epsilonr))
    {
        epsilonr = 1;
    }

    ForceComputer forceComp(shellToler, shellMaxIter, epsilonr);
    
    JsonTree jtree("SecondVirialCoefficient");
    if (msghandler.info())
    {
        forceFieldSummary(&jtree, &pd);
    }

    std::vector<ACTMol> actmols;
    compR.read(&msghandler, pd, &forceComp, &actmols);
    auto &actmol = actmols[0];

    if (msghandler.ok() && status == 0)
    {
        if (msghandler.debug())
        {
            actmol.topology()->dump(msghandler.twDebug());
        }
        phase.next("Second virial");
        rerun.setFunctions(&forceComp, oenv);
        rerun.runB2(&cr, &msghandler, &pd, &actmol, maxdimers, fnm);
    }
    phase.next("Output");
    bool writeJson = json && cr.isMaster();
    if (!writeJson)
    {
        int indent = 0;
        msghandler.write(jtree.writeString(json, &indent));
    }
    phase.stop();
    // The timings go into the json file, or else only into the log file
    printTimings(&cr, cr.isMaster() ? msghandler.tw() : nullptr, phase.wallTime(),
                 writeJson ? &jtree : nullptr);
    if (writeJson)
    {
        jtree.write("b2.json", json);
    }
    return status;
}

} // namespace alexandria
