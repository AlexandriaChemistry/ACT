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
#include "act/forces/forcecomputerstatistics.h"
#include "act/forces/simd_information.h"
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
                              chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm())));
    // The same for van der Waals and Coulomb
    int nexcl = std::max(0, pd->nexcl());
    jtree->addObject(JsonTree("# vanderwaals exclusions", gmx_itoa(nexcl)));
    jtree->addObject(JsonTree("# coulomb exclusions", gmx_itoa(nexcl)));
    double epsilonr;
    epsilonr = pd->epsilonr();
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
            const auto function = potentialToString(ffpl.potential());
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
          "Plot the total B2(T) only and not the components" },
        { "-block", FALSE, etINT, {&blockSize_},
          "Largest number of distances of an orientation that are computed together, the distances are divided into blocks of equal size. 0 means all of them. Blocks of a few hundred distances fit in the cache of a core, which is faster than computing all distances at once, in particular when many processes share a processor." },
        { "-emin", FALSE, etREAL, {&eMin_},
          "Dimers with a lower interaction energy (kJ/mol) are treated as overlapping atoms, with a Mayer function of -1. In models without repulsion on all charged sites, such as TIP5P, a charge can come arbitrarily close to an opposite charge of the other compound, with an energy that goes to minus infinity, and the second virial would diverge. The number of these dimers is reported." }
    };
    pargs->push_back(pa[0]);
    if (computeB2_)
    {
        for(size_t i = 1; i < pa.size(); i++)
        {
            pargs->push_back(pa[i]);
        }
    
        const std::vector<t_filenm> fnm = {
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


void ReRunner::addB2(JsonTree *jtree)
{
    auto &T = temperatures();
    if (b2t_.empty() || b2t_[b2Type::Total].size() != T.size())
    {
        return;
    }
    JsonTree b2tree("Second virial coefficients");
    b2tree.addObject(JsonTree("Unit", "cm^3/mol"));
    for(size_t ii = 0; ii < T.size(); ii++)
    {
        JsonTree ttree(gmx::formatString("T = %g K", T[ii]));
        for(const auto &b2b : b2Type2str)
        {
            ttree.addObject(JsonTree(b2b.second, gmx::formatString("%.17g", b2t_[b2b.first][ii])));
        }
        b2tree.addObject(ttree);
    }
    jtree->addObject(b2tree);
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
        const std::vector<std::string> legend = {
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
 * \param[in]  atoms  The particles
 * \param[in]  x      The input coordinates
 * \param[out] coords The coordinates of all particles, its memory is used again
 */
static void allParticleCoords(const std::vector<ActAtom>   &atoms,
                              const std::vector<gmx::RVec> &x,
                              std::vector<gmx::RVec>       *coords)
{
    if (x.size() == atoms.size())
    {
        coords->assign(x.begin(), x.end());
        return;
    }
    coords->clear();
    size_t index = 0;
    for(size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].pType() == ActParticle::Atom)
        {
            if (index >= x.size())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Number of (generated) coordinates in trajectory (%zu) does not match molecule file (%zu)", x.size(), atoms.size()).c_str()));
            }
            coords->push_back(x[index++]);
        }
        else if (!atoms[i].parents().empty() && atoms[i].parents()[0] < static_cast<int>(i))
        {
            coords->push_back((*coords)[atoms[i].parents()[0]]);
        }
        else
        {
            // Will be constructed, e.g. a virtual site
            coords->emplace_back(0, 0, 0);
        }
    }
    if (index != x.size())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Number of (generated) coordinates in trajectory (%zu) does not match molecule file (%zu)", x.size(), atoms.size()).c_str()));
    }
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
    std::vector<std::vector<gmx::RVec>> allCoords(dimers.size());
    for (size_t i = 0; i < dimers.size(); i++)
    {
        allParticleCoords(atoms, dimers[i], &allCoords[i]);
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
            allCoords[c] = set.data.x(c).toRVec();
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
    const std::vector<double> masses = {
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
        const ScopedTimer   scoped(timer);
        gendimers_.generateRandomNumbers(maxdimer);
        allRand = gendimers_.allRandom();
        gendimers_.clearAllRandom();
    }
    /* Distribute the orientations from start to end over the nodes, the
     * master takes the first ones, including the rest of the division.
     * The random numbers are added to those of the previous rounds.
     * Returns the number of orientations of this node.
     */
    const auto distribute = [&](int start, int end)
    {
        static const size_t timer = Timings::instance().index("Distributing orientations");
        const ScopedTimer   scoped(timer);
        int nper  = (end - start) / cr->size();
        int nrest = (end - start) % cr->size();
        int nmine = 0;
        if (cr->isMaster())
        {
            // The random numbers of each node in one message, with
            // the number of them for each orientation.
            std::vector<int>    counts;
            std::vector<double> numbers;
            for(int dst = 1; dst < cr->size(); dst++)
            {
                counts.clear();
                numbers.clear();
                int nn = start + nrest + dst*nper;
                for(int j = 0; j < nper; j++)
                {
                    counts.push_back(allRand[nn + j].size());
                    numbers.insert(numbers.end(), allRand[nn + j].begin(), allRand[nn + j].end());
                }
                cr->send(dst, counts);
                cr->send(dst, numbers);
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
            int                 src = cr->superior();
            std::vector<int>    counts;
            std::vector<double> numbers;
            cr->recv(src, &counts);
            cr->recv(src, &numbers);
            nmine = counts.size();
            auto q = numbers.begin();
            for(int count : counts)
            {
                gendimers_.addRandomNumbers(std::vector<double>(q, q + count));
                q += count;
            }
        }
        return nmine;
    };
    // Will be used to obtain a seed for the random number engine for bootstrapping
    std::random_device                 bsRand;
    // Standard mersenne_twister_engine seeded with rd()
    const std::mt19937                 bsGen(bsRand());
    const std::uniform_int_distribution<int> bsDistr(0, gendimers_.ndist());

    // Temperature array.
    auto   Temperature = temperatures();
    for(const auto &b2b : b2Type2str)
    {
        b2t_.insert({b2b.first, {}});
        b2t_[b2b.first].resize(Temperature.size());
    }
    // Temporary arrays for weighted properties.
    B2Data b2data(gendimers_.ndist(), Temperature, eMin_);
    
    const gmx::RVec rvnul = { 0, 0, 0 };
    // Moments of inertia, averaged over all dimers. The sum is over many
    // dimers, in double precision.
    std::vector<gmx::RVec> inertia    = { rvnul, rvnul };
    std::vector<gmx::DVec> inertiaSum = { { 0, 0, 0 }, { 0, 0, 0 } };
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
    B2Data b2total(gendimers_.ndist(), Temperature, eMin_);
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
    // Coordinates of all particles, interaction energies and forces of the
    // dimers of an orientation, kept such that their memory is used again
    std::vector<std::vector<gmx::RVec>> allCoordsShells;
    DimerInteractions                   inter;
    // Per molecule the center of mass, the total force, the torques and the
    // principal axes of the dimers of an orientation, with the dimers in the
    // inner loop like in inter, and the coordinates relative to the center
    // of mass of one molecule
    std::vector<double> comDim[2], forceDim[2], torqueDim[2], torqueRotDim[2], transDim[2];
    std::vector<double> xcDim;
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
                const ScopedTimer scoped(prepareTimer);
                gendimers_.generateDimers(msghandler, &dimers);
                if (gendimers_.outCoords())
                {
                    for (auto &d : dimers)
                    {
                        allCoords.push_back(d.coords);
                    }
                }
            }
            size_t ndimAll = dimers.size();
            for(int kk = 0; kk < 2; kk++)
            {
                force2Mol[kk].resize(ndimAll);
                torqueMol[kk].resize(ndimAll);
            }
            dimIndex.resize(ndimAll);
            dimDist.resize(ndimAll);
            dimEnergy.resize(ndimAll);
            dimWeight.resize(ndimAll);
            // The dimers are computed in blocks of distances of equal size, see -block
            size_t block = std::max<size_t>(ndimAll, 1);
            if (blockSize_ > 0)
            {
                size_t bmax   = static_cast<size_t>(blockSize_);
                size_t nblock = std::max<size_t>((ndimAll + bmax - 1)/bmax, 1);
                block         = std::max<size_t>((ndimAll + nblock - 1)/nblock, 1);
            }
            for(size_t b0 = 0; b0 < ndimAll; b0 += block)
            {
                size_t ndim = std::min(block, ndimAll - b0);
                // Coordinates of all particles of the dimers
                timings.start(energyTimer);
                allCoordsShells.resize(ndim);
                for (size_t i = 0; i < ndim; i++)
                {
                    allParticleCoords(atoms, dimers[b0+i].coords, &allCoordsShells[i]);
                }
                // Compute the interaction energies and forces of the dimers at once
                forceComp_->computeInteraction(msghandler, pd, actmol, allCoordsShells, true, &inter);
                timings.stop(energyTimer);
                // Loop over dimers
                timings.start(torqueTimer);
                // The dimers are in the inner loop of the arrays, such that the
                // loops over them can be vectorized
                size_t        st      = inter.stride;
                const double *xs      = inter.x.data();
                const double *fs      = inter.f.data();
                bool          verbose = msghandler->verbose();
                for(int kk = 0; kk < 2; kk++)
                {
                    double mtot = 0;
                    for(size_t i = atomStart[kk]; i < atomStart[kk]+natom[kk]; i++)
                    {
                        mtot += atoms[i].mass();
                    }
                    GMX_RELEASE_ASSERT(mtot > 0, "Zero mass");
                    // Compute center of mass and total force of compound kk
                    auto &com  = comDim[kk];
                    auto &ftot = forceDim[kk];
                    com.assign(DIM*st, 0.0);
                    ftot.assign(DIM*st, 0.0);
                    for(size_t i = atomStart[kk]; i < atomStart[kk]+natom[kk]; i++)
                    {
                        double mi = atoms[i].mass();
                        for(int m = 0; m < DIM; m++)
                        {
                            const double *x  = xs + (DIM*i+m)*st;
                            const double *f  = fs + (DIM*i+m)*st;
                            double       *cm = com.data() + m*st;
                            double       *fm = ftot.data() + m*st;
                            for(size_t c = 0; c < ndim; c++)
                            {
                                double mr = mi*x[c];
                                cm[c] += mr;
                                fm[c] += f[c];
                            }
                        }
                    }
                    for(size_t k = 0; k < DIM*ndim; k++)
                    {
                        // Normalize, the rows are ndim apart here
                        com[(k/ndim)*st + k%ndim] /= mtot;
                    }
                    // Compute the coordinates relative to the center of mass
                    xcDim.resize(DIM*natom[kk]*st);
                    for(size_t j = 0; j < natom[kk]; j++)
                    {
                        for(int m = 0; m < DIM; m++)
                        {
                            const double *x  = xs + (DIM*(atomStart[kk]+j)+m)*st;
                            const double *cm = com.data() + m*st;
                            double       *xc = xcDim.data() + (DIM*j+m)*st;
                            for(size_t c = 0; c < ndim; c++)
                            {
                                xc[c] = x[c] - cm[c];
                            }
                        }
                    }
                    // Compute moments of inertia and transformation matrix.
                    // Rigid molecules are only translated between the dimers
                    // of an orientation, so their principal axes are the same.
                    auto &transk = transDim[kk];
                    transk.resize(DIM*DIM*st);
                    for(size_t c = 0; c < ndim; c++)
                    {
                        if (flexible || 0 == b0 + c)
                        {
                            auto &xc = xCom[kk];
                            for(size_t j = 0; j < natom[kk]; j++)
                            {
                                for(int m = 0; m < DIM; m++)
                                {
                                    xc[j][m] = xcDim[(DIM*j+m)*st + c];
                                }
                            }
                            clear_rvec(inertia1[kk]);
                            clear_mat(trans[kk]);
                            principal_comp(molIndex[kk], molMass[kk], xc, &trans[kk], &inertia1[kk]);
                        }
                        for(int m = 0; m < DIM; m++)
                        {
                            inertiaSum[kk][m] += inertia1[kk][m];
                        }
                        for(int a = 0; a < DIM; a++)
                        {
                            for(int b = 0; b < DIM; b++)
                            {
                                transk[(DIM*a+b)*st + c] = trans[kk][a][b];
                            }
                        }
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
                    auto &trot = torqueRotDim[kk];
                    trot.assign(DIM*st, 0.0);
                    for(size_t j = 0; j < natom[kk]; j++)
                    {
                        const double *xc = xcDim.data() + DIM*j*st;
                        const double *f  = fs + DIM*(atomStart[kk]+j)*st;
                        for(size_t c = 0; c < ndim; c++)
                        {
                            const dvec xi = { xc[c], xc[st+c], xc[2*st+c] };
                            const dvec fc = { f[c], f[st+c], f[2*st+c] };
                            // Rotate coordinates and force vector, in
                            // double precision as the input
                            dvec ri = { 0, 0, 0 }, fi = { 0, 0, 0 }, ti;
                            for(int a = 0; a < DIM; a++)
                            {
                                for(int b = 0; b < DIM; b++)
                                {
                                    double tab = transk[(DIM*a+b)*st + c];
                                    ri[a] += tab*xi[b];
                                    fi[a] += tab*fc[b];
                                }
                            }
                            // Compute torque on this atom
                            dcprod(ri, fi, ti);
                            // Update total torque
                            for(int m = 0; m < DIM; m++)
                            {
                                trot[m*st + c] += ti[m];
                            }
                        }
                    }
                    // The torque in the original frame is only printed
                    if (verbose)
                    {
                        auto &tq = torqueDim[kk];
                        tq.assign(DIM*st, 0.0);
                        for(size_t j = 0; j < natom[kk]; j++)
                        {
                            const double *xc = xcDim.data() + DIM*j*st;
                            const double *f  = fs + DIM*(atomStart[kk]+j)*st;
                            for(size_t c = 0; c < ndim; c++)
                            {
                                dvec xi = { xc[c], xc[st+c], xc[2*st+c] };
                                dvec fc = { f[c], f[st+c], f[2*st+c] };
                                dvec ti;
                                dcprod(xi, fc, ti);
                                for(int m = 0; m < DIM; m++)
                                {
                                    tq[m*st + c] += ti[m];
                                }
                            }
                        }
                    }
                }
                nInertia += ndim;
                for (size_t idim = 0; idim < ndim; idim++)
                {
                    auto      einter       = inter.einter[idim];
                    auto      epot         = inter.epot[idim];
                    gmx::RVec f[2], com[2], torque[2], torqueRot[2];
                    for(int kk = 0; kk < 2; kk++)
                    {
                        for(int m = 0; m < DIM; m++)
                        {
                            f[kk][m]         = forceDim[kk][m*st + idim];
                            com[kk][m]       = comDim[kk][m*st + idim];
                            torque[kk][m]    = verbose ? torqueDim[kk][m*st + idim] : 0;
                            torqueRot[kk][m] = torqueRotDim[kk][m*st + idim];
                        }
                        force2Mol[kk][b0+idim] = iprod(f[kk], f[kk]);
                        torqueMol[kk][b0+idim] = torqueRot[kk];
                    }
                    // Check consistency from the rotation algorithm.
                    rvec dx;
                    rvec_sub(com[0], com[1], dx);
                    double rcom = norm(dx);
                    if (std::abs(rcom - dimers[b0+idim].dist) > 1e-4)
                    {
                        msghandler->fatal(gmx::formatString(" rcom calculated %g from input %g",
                                                            rcom, dimers[b0+idim].dist));
                    }
                    rcom = dimers[b0+idim].dist;
                    if (msghandler->verbose())
                    {
                        std::string out = gmx::formatString(" r %g epot %g einter %g", rcom, epot, einter);
                        out += gmx::formatString(" Force %g %g %g Torque[0] %g %g %g Torque[1] %g %g %g Rotated Torque[0] %g %g %g Rotated Torque[1] %g %g %g",
                                                 f[0][XX], f[0][YY], f[0][ZZ],
                                                 torque[0][XX], torque[0][YY], torque[0][ZZ],
                                                 torque[1][XX], torque[1][YY], torque[1][ZZ],
                                                 torqueRot[0][XX], torqueRot[0][YY], torqueRot[0][ZZ],
                                                 torqueRot[1][XX], torqueRot[1][YY], torqueRot[1][ZZ]);
                        msghandler->write(out);
                    }
                    dimIndex[b0+idim]  = dimers[b0+idim].index;
                    dimDist[b0+idim]   = rcom;
                    dimEnergy[b0+idim] = flexible ? epot : einter;
                    dimWeight[b0+idim] = dimers[b0+idim].weight;
                    if (msghandler->verbose())
                    {
                        msghandler->write(gmx::formatString("epot %g einter %g diff %g", epot, einter, epot-einter));
                    }
                }
                // Per particle the center of mass and the force 10, the coordinates
                // relative to it 3 and the torques 54, per molecule the normalization
                // with three divisions, the moments of inertia and the squared
                // force 23, and the distance between the centers of mass with its
                // square root 13. The principal axes are not counted.
                timings.addFlops(ndim*(67.0*(natom[0] + natom[1]) + 2*23 + 13));
                timings.stop(torqueTimer);
            }
            if (msghandler->verbose())
            {
                msghandler->write("");
            }
            // Now we have computed energies, forces and torques for all
            // dimers of this orientation. Time to update the intermediate structure
            const ScopedTimer scoped(addTimer);
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
            const ScopedTimer scoped(b2Timer);
            msghandler->write(gmx::formatString("\nSecond virial coefficient (cm^3/mol) from %d dimer orientations", nround));
            if (b2total.nOverlap() > 0)
            {
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("%.0f of %.0f dimers have an energy below %g kJ/mol and are treated as overlapping atoms, see -emin.",
                                                  b2total.nOverlap(), b2total.nDimer(), eMin_));
            }
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

    const ScopedTimer scoped(outputTimer);
    // Gather generated dimer coordinates from all nodes and write to -ox file
    if (gendimers_.outCoords() && !gendimers_.hasTrajectory())
    {
        // The coordinates of each node in one message, with the number
        // of them for each dimer
        std::vector<int>       counts;
        std::vector<gmx::RVec> coords;
        if (!cr->isMaster())
        {
            int dest = cr->superior();
            for (auto &d : allCoords)
            {
                counts.push_back(d.size());
                coords.insert(coords.end(), d.begin(), d.end());
            }
            cr->send(dest, counts);
            cr->send(dest, coords);
        }
        else
        {
            for (int src = 1; src < cr->size(); src++)
            {
                cr->recv(src, &counts);
                cr->recv(src, &coords);
                auto x = coords.begin();
                for (int count : counts)
                {
                    allCoords.emplace_back(x, x + count);
                    x += count;
                }
            }
            gendimers_.writeCoords(actmol, allCoords);
        }
    }

    if (cr->isMaster())
    {
        // Below the shortest distance the Mayer function is taken as -1,
        // which is only correct if it is -1 at that distance
        size_t nsoft = 0;
        double fmax  = -1, Tmax = 0;
        for(size_t iTemp = 0; iTemp < Temperature.size(); iTemp++)
        {
            double f = b2total.mayerAtShortestDistance(iTemp);
            if (f > -1 + 1e-6)
            {
                nsoft += 1;
                if (f > fmax)
                {
                    fmax = f;
                    Tmax = Temperature[iTemp];
                }
            }
        }
        if (nsoft > 0)
        {
            msghandler->msg(ACTStatus::Warning,
                            gmx::formatString("The Mayer function at the shortest distance %g nm is larger than -1 at %zu temperature(s), up to %g at %g K. Since it is taken as -1 below that distance, please decrease -mindist.",
                                              gendimers_.mindist(), nsoft, fmax, Tmax));
        }
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
        "in the charges file.[PAR]",
        "With [TT]-json[tt], part of the output, including the second virial",
        "coefficients with all digits, is written in json format to the given",
        "file, b2.json by default."
    };

    std::vector<t_filenm>     fnm = {
        { efXML,  "-ff",     "aff",     ffREAD  },
        { efJSON, "-json",   "b2",      ffOPTWR }
    };
    gmx_output_env_t         *oenv;
    char                     *molnm            = (char *)"";
    real                      shellToler       = 1e-6;
    int                       maxdimers        = 1024;
    int                       shellMaxIter     = 25;
    std::vector<t_pargs>      pa = {
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule." },
        { "-maxdimer", FALSE, etINT, {&maxdimers},
          "Number of dimer orientations to generate if you do not provide a trajectory. For each of these a distance scan will be performed." },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)." },
        { "-shellMaxIter", FALSE, etINT, {&shellMaxIter},
          "Max number of iterations for minimizing shell positions" }
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
    // Part of the output in json format, including the second virial
    // coefficients with all digits, in the file of -json if it is given
    const char *jsonFile = opt2fn_null("-json", fnm.size(), fnm.data());
    bool        json     = jsonFile != nullptr;
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
    epsilonr = pd.epsilonr();

    ForceComputer forceComp(shellToler, shellMaxIter, epsilonr);
    
    JsonTree jtree("SecondVirialCoefficient");
    if (msghandler.info())
    {
        forceFieldSummary(&jtree, &pd);
    }

    std::vector<ACTMol> actmols;
    compR.read(&msghandler, pd, &forceComp, &actmols);
    if (actmols.empty())
    {
        fprintf(stderr, "Could not read or process molecules. Please check the log file for more information.\n");
        return 1;
    }
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
        if (json && cr.isMaster())
        {
            rerun.addB2(&jtree);
        }
    }
    if (pd.polarizable())
    {
        // The number of shell iterations per force calculation, summed
        // over the processes
        ForceComputerStatistics fcStats;
        const auto stats = fcStats.statistics(&cr, &forceComp, 0, true);
        if (cr.isMaster())
        {
            msghandler.write(stats);
        }
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
    const auto peak = theoreticalPeakPerCore();
    printTimings(&cr, cr.isMaster() ? msghandler.tw() : nullptr, phase.wallTime(),
                 writeJson ? &jtree : nullptr, &peak);
    if (writeJson)
    {
        jtree.write(jsonFile, json);
    }
    return status;
}

} // namespace alexandria
