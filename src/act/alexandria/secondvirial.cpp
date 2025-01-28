/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2025
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
#include <cstdlib>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/atype_mapping.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/b2data.h"
#include "act/alexandria/compound_reader.h"
#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/actmol.h"
#include "act/alexandria/princ.h"
#include "act/alexandria/train_utility.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/jsontree.h"
#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

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

void ReRunner::addOptions(std::vector<t_pargs>  *pargs,
                          std::vector<t_filenm> *filenm)
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

void ReRunner::rerun(FILE             *logFile,
                     const ForceField *pd,
                     const ACTMol     *actmol,
                     bool              verbose)
{
    std::vector<std::vector<gmx::RVec>> dimers;
    gendimers_->read(&dimers);
    if (logFile)
    {
        fprintf(logFile, "Doing energy calculation for %zu structures from %s\n",
                dimers.size(), gendimers_->trajname());
        fflush(logFile);
    }
    if (verbose && debug)
    {
        print_memory_usage(debug);
    }
    std::map<InteractionType, double> energies;
    int mp_index = 0;
    gmx_stats edist;
    // Loop over molecules
    const auto &atoms = actmol->atomsConst();
    for (size_t idim = 0; idim < dimers.size(); idim++)
    {
        std::vector<gmx::RVec> coords;
        if (dimers[idim].size() == atoms.size())
        {
            // Assume there are shells in the input
            coords = dimers[idim];
        }
        else
        {
            size_t index = 0;
            for(size_t i = 0; i < atoms.size(); i++)
            {
                if (index < dimers[idim].size())
                {
                    coords.push_back(dimers[idim][index]);
                }
                else
                {
                    GMX_THROW(gmx::InvalidInputError(gmx::formatString("Number of (generated) coordinates in trajectory (%zu) does not match molecule file (%zu)", dimers[idim].size(), atoms.size()).c_str()));
                }
                if (atoms[i].pType() == ActParticle::Atom)
                {
                    index++;
                }
            }
        }
        std::vector<gmx::RVec> forces(coords.size());
        if (verbose)
        {
            fprintf(logFile, "%5d", mp_index);
        }
        if (eInter_)
        {
            std::map<InteractionType, double> einter;
            actmol->calculateInteractionEnergy(pd, forceComp_, &einter, &forces, &coords, true);

            auto atomStart  = actmol->fragmentHandler()->atomStart();
            std::vector<gmx::RVec> com        = { { 0, 0, 0 }, { 0, 0, 0 } };
            std::vector<double>    mtot       = { 0, 0 };
            auto      tops  = actmol->fragmentHandler()->topologies();
            for(int kk = 0; kk < 2; kk++)
            {
                auto natom = tops[kk]->atoms().size();
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
            if (logFile && verbose)
            {
                fprintf(logFile, " r %g", rcom);
                for (auto &EE: einter)
                {
                    fprintf(logFile, " %s %g", interactionTypeToString(EE.first).c_str(), EE.second);
                }
            }
            edist.add_point(rcom, einter[InteractionType::EPOT], 0, 0);
        }
        else
        {
            forceComp_->compute(pd, actmol->topology(),
                                &coords, &forces, &energies);
            for(const auto &ee : energies)
            {
                fprintf(logFile, "  %s %8g", interactionTypeToString(ee.first).c_str(), ee.second);
            }
        }
        if (verbose)
        {
            fprintf(logFile, "\n");
        }
        mp_index++;
    }
    if (verbose && debug)
    {
        print_memory_usage(debug);
    }
}

void ReRunner::runB2(CommunicationRecord         *cr,
                     FILE                        *logFile,
                     const ForceField            *pd,
                     const ACTMol                *actmol,
                     int                          maxdimer,
                     bool                         verbose,
                     const std::vector<t_filenm> &fnm)
{
    // Compute the relative masses
    std::vector<double> masses = {
        actmol->fragmentHandler()->topologies()[0]->mass(),
        actmol->fragmentHandler()->topologies()[1]->mass()
    };
    int ndimer = 0;
    int nrest  = 0;
    std::vector<std::vector<gmx::RVec>> dimers;
    if (gendimers_->hasTrajectory())
    {
        if (cr->isMaster())
        {
            gendimers_->read(&dimers);
        }
        maxdimer = dimers.size();
        ndimer   = 1;
    }
    else
    {
        // Do this in parallel and with little memory
        ndimer = maxdimer / cr->size();
        nrest  = maxdimer % cr->size();

        // Make sure that different nodes have different random number generator seed.
        // For each dimer, 2 x 3 = 6 random numbers are needed. That means, each node
        // uses ndimer x 6. However (see below), if nrest != 0, the master will get that
        // number. 
        if (cr->isMaster())
        {
            gendimers_->generateRandomNumbers(maxdimer);
            auto allRand = gendimers_->allRandom();
            int ioffset = 0;
            if (nrest > 0)
            {
                ioffset = 1;
            }
            for(int dst = 1; dst < cr->size(); dst++)
            {
                cr->send(dst, ndimer);
                int nn = nrest + (dst-ioffset)*ndimer;
                for(int j = 0; j < ndimer; j++)
                {
                    cr->send(dst, allRand[nn + j]);
                }
            }
            if (nrest > 0)
            {
                ndimer = nrest;
            }
            gendimers_->resetAllRandom(ndimer);
            if (logFile)
            {
                fprintf(logFile, "Will generate or read %d dimers on helpers and %d on master.\n", ndimer, nrest);
            }
        }
        else
        {
            int src = cr->superior();
            cr->recv(src, &ndimer);
            for(int i = 0; i < ndimer; i++)
            {
                std::vector<double> q;
                cr->recv(src, &q);
                gendimers_->addRandomNumbers(q);
            }
        }
    }
    if (logFile)
    {
        fflush(logFile);
    }
    // Will be used to obtain a seed for the random number engine for bootstrapping
    std::random_device                 bsRand;
    // Standard mersenne_twister_engine seeded with rd()
    std::mt19937                       bsGen(bsRand());
    std::uniform_int_distribution<int> bsDistr(0, gendimers_->ndist());

    // Bins start from zero for proper integration
    size_t nbins    = 1+std::round(gendimers_->maxdist()/gendimers_->binwidth());

    // Temperature array.
    auto   Temperature = temperatures();
    for(const auto &b2b : b2Type2str)
    {
        b2t_.insert({b2b.first, {}});
        b2t_[b2b.first].resize(Temperature.size());
    }
    // Temporary arrays for weighted properties.
    B2Data b2data(nbins, gendimers_->binwidth(), Temperature);
    
    double xmin = 0;
    std::vector<gmx::RVec> inertia = { { 0, 0, 0 }, { 0, 0, 0 } };
    //! Loop over my dimers, completely independently
    for(int idimer = 0; idimer < ndimer; idimer++)
    {
        // Generate a new set of dimers for all distances
        if (!gendimers_->hasTrajectory())
        {
            dimers = gendimers_->generateDimers(debug, actmol);
        }
        // Structures to store energies, forces and torques
        gmx_stats                           edist;
        std::map<InteractionType, double>   energies;
        std::vector<std::vector<gmx::RVec>> forceMol;
        std::vector<std::vector<gmx::RVec>> torqueMol;
        forceMol.resize(2);
        torqueMol.resize(2);
        for(int kk = 0; kk < 2; kk++)
        {
            torqueMol[kk].resize(dimers.size());
            forceMol[kk].resize(dimers.size());
        }

        // Loop over molecules
        const auto &atoms = actmol->atomsConst();
        for (size_t idim = 0; idim < dimers.size(); idim++)
        {
            std::vector<gmx::RVec> coords;
            if (dimers[idim].size() == atoms.size())
            {
                // Assume there are shells in the input
                coords = dimers[idim];
            }
            else
            {
                size_t index = 0;
                for(size_t i = 0; i < atoms.size(); i++)
                {
                    if (index < dimers[idim].size())
                    {
                        coords.push_back(dimers[idim][index]);
                    }
                    else
                    {
                        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Number of (generated) coordinates in trajectory (%zu) does not match molecule file (%zu)", dimers[idim].size(), atoms.size()).c_str()));
                    }
                    if (atoms[i].pType() == ActParticle::Atom)
                    {
                        index++;
                    }
                }
            }
            std::vector<gmx::RVec> forces(coords.size());
            std::map<InteractionType, double> einter;
            actmol->calculateInteractionEnergy(pd, forceComp_, &einter, &forces, &coords, true);

            auto atomStart  = actmol->fragmentHandler()->atomStart();
            std::vector<gmx::RVec> f          = { { 0, 0, 0 }, { 0, 0, 0 } };
            std::vector<gmx::RVec> com        = { { 0, 0, 0 }, { 0, 0, 0 } };
            std::vector<double>    mtot       = { 0, 0 };
            std::vector<gmx::RVec> torque     = { { 0, 0, 0 }, { 0, 0, 0 } };
            std::vector<gmx::RVec> torqueRot  = { { 0, 0, 0 }, { 0, 0, 0 } };
            auto      tops  = actmol->fragmentHandler()->topologies();
            for(int kk = 0; kk < 2; kk++)
            {
                auto natom = tops[kk]->atoms().size();
                for(size_t i = atomStart[kk]; i < atomStart[kk]+natom; i++)
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
                copy_rvec(f[kk], forceMol[kk][idim]);

                // Compute the coordinates relative to the center of mass
                std::vector<real>      mass;
                std::vector<int>       index;
                std::vector<gmx::RVec> x_com;
                for(size_t i = atomStart[kk]; i < atomStart[kk]+natom; i++)
                {
                    // Store atom index relative to molecule start
                    index.push_back(i-atomStart[kk]);
                    // Store mass
                    mass.push_back(atoms[i].mass());
                    // Subtract COM and store
                    gmx::RVec ri;
                    rvec_sub(coords[i], com[kk], ri);
                    x_com.push_back(ri);
                }
                // Compute moments of inertia and transformation matrix
                gmx::RVec inertia1;
                clear_rvec(inertia1);
                matrix trans;
                principal_comp(index, mass, x_com, &trans, &inertia1);

                // Move to inertial frame (only well-defined for
                // rigid molecules).
                // The trans matrix should convert that coordinate to the inertial frame,
                // but what about the force on the atoms? It likely has to be rotated in the
                // same manner. After that, the torque can be computed.
                // Before doing the rotations, the sum of the torque vectors is zero.
                // Since the rotations are different for both molecules, this does not
                // hold after the rotations are done.
                // TODO: write out the math.
                for(size_t i = atomStart[kk]; i < atomStart[kk]+natom; i++)
                {
                    gmx::RVec ri, fi, ti;
                    cprod(x_com[i-atomStart[kk]], forces[i], ti);
                    rvec_inc(torque[kk], ti);
                    // Rotate coordinates
                    mvmul(trans, x_com[i-atomStart[kk]], ri);
                    // Rotate force vector
                    mvmul(trans, forces[i], fi);
                    // Compute torque on this atom
                    cprod(ri, fi, ti);
                    // Update total torque
                    rvec_inc(torqueRot[kk], ti);
                }
                rvec_inc(inertia[kk], inertia1);
                torqueMol[kk][idim] = torqueRot[kk];
            }
            gmx::RVec dcom;
            rvec_sub(com[0], com[1], dcom);
            double rcom = norm(dcom);
            if (verbose && logFile)
            {
                fprintf(logFile, " r %g", rcom);
                for (auto &EE: einter)
                {
                    fprintf(logFile, " %s %g", interactionTypeToString(EE.first).c_str(), EE.second);
                }
                fprintf(logFile, " Force %g %g %g Torque[0] %g %g %g Torque[1] %g %g %g Rotated Torque[0] %g %g %g Rotated Torque[1] %g %g %g",
                        f[0][XX], f[0][YY], f[0][ZZ],
                        torque[0][XX], torque[0][YY], torque[0][ZZ],
                        torque[1][XX], torque[1][YY], torque[1][ZZ],
                        torqueRot[0][XX], torqueRot[0][YY], torqueRot[0][ZZ],
                        torqueRot[1][XX], torqueRot[1][YY], torqueRot[1][ZZ]);
                fprintf(logFile, "\n");
            }
            edist.add_point(rcom, einter[InteractionType::EPOT], 0, 0);
        }
        if (verbose && logFile)
        {
            fprintf(logFile, "\n");
        }
        for(int kk = 0; kk < 2; kk++)
        {
            for(int m = 0; m < DIM; m++)
            {
                inertia[kk][m] /= edist.get_npoints();
            }
        }
        // Now we have computed energies, forces and torques for one dimer.
        // Time to update the intermediate structure
        auto &x = edist.getX();
        auto &y = edist.getY();
        xmin = *std::min_element(x.begin(), x.end());
        for(size_t iTemp = 0; iTemp < Temperature.size(); iTemp++)
        {
            auto T      = Temperature[iTemp];
            if (T == 0)
            {
                fprintf(stderr, "Please provide a finite temperature to compute second virial.\n");
                continue;
            }
            double beta = 1.0/(BOLTZ*T);
            for(size_t jj = 0; jj < x.size(); jj++)
            {
                size_t ii = jj;
                size_t index = x[ii]/gendimers_->binwidth();
                // Gray and Gubbins Eqn. 3.261
                double g0_12 = std::exp(-y[ii]*beta);
                gmx::RVec tau[2];
                for(int kk = 0; kk < 2; kk++)
                {
                    clear_rvec(tau[kk]);
                    for(int m = 0; m < DIM; m++)
                    {
                        // Gray and Gubbins Eqn. 3.282
                        tau[kk][m] = g0_12*torqueMol[kk][ii][m]*torqueMol[kk][ii][m];
                    }
                }
                // Gray and Gubbins Eqn. 3.272
                b2data.addData(iTemp, index, g0_12-1,
                               // Gray and Gubbins Eqn. 3.281
                               g0_12*iprod(forceMol[0][ii], forceMol[0][ii]),
                               g0_12*iprod(forceMol[1][ii], forceMol[1][ii]),
                               tau[0], tau[1]);
            }
        }
    }
    // Sum b2data over processors
    b2data.aggregate(cr);

    if (cr->isMaster())
    {
        // Starting energy, all values until first data entry
        
        for(size_t iTemp = 0; iTemp < Temperature.size(); iTemp++)
        {
            auto T      = Temperature[iTemp];
            double beta = 1.0/(BOLTZ*T);
            b2data.fillToXmin(iTemp, xmin, gendimers_->binwidth());
            // Now compute the components
            double Bclass, BqmForce, BqmTorque1, BqmTorque2; 
            b2data.integrate(iTemp, gendimers_->binwidth(), beta, masses, inertia,
                             &Bclass, &BqmForce, &BqmTorque1, &BqmTorque2);

            // Conversion to regular units cm^3/mol.
            double fac  = AVOGADRO*1e-21;
            // TODO: Fix the torque contribution
            double bqt  = (BqmTorque1+BqmTorque2)*0.5;
            double Btot = (Bclass + BqmForce + bqt)*fac;
            // Add to bootstrapping statistics
            b2t_[b2Type::Classical][iTemp] = Bclass*fac;
            b2t_[b2Type::Force][iTemp]     = BqmForce*fac;
            b2t_[b2Type::Torque1][iTemp]   = BqmTorque1*fac;
            b2t_[b2Type::Torque2][iTemp]   = BqmTorque2*fac;
            b2t_[b2Type::Total][iTemp]     = Btot;
            
            if (logFile)
            {
                fprintf(logFile, "T = %g K. ", T);
                for(const auto &b2b : b2Type2str)
                {
                    fprintf(logFile, " %s %8.1f", b2b.second.c_str(),
                            b2t_[b2b.first][iTemp]);
                }
                fprintf(logFile, "\n");
            }
        }
        if (!fnm.empty())
        {
            // Store Mayer functions if requested.
            b2data.plotMayer(opt2fn_null("-eh", fnm.size(), fnm.data()), oenv_);
            // Print B2(T) if requested
            plotB2temp(opt2fn("-b2", fnm.size(), fnm.data()));
        }
    }
}

int b2(int argc, char *argv[])
{
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
        { efXML, "-ff",      "aff",     ffREAD  },
        { efLOG, "-g",       "b2",      ffWRITE } 
    };
    gmx_output_env_t         *oenv;
    static char              *molnm      = (char *)"";
    double                    shellToler = 1e-6;
    int                       maxdimers  = 1024;
    bool                      verbose    = false;
    bool                      oneH       = false;
    bool                      json       = false;
    std::vector<t_pargs>      pa = {
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule." },
        { "-v", FALSE, etBOOL, {&verbose},
          "Print more information to the log file." },
        { "-oneH", FALSE, etBOOL, {&oneH},
          "Map all different hydrogen atom types back to H, mainly for debugging." },
        { "-maxdimer", FALSE, etINT, {&maxdimers},
          "Number of dimer orientations to generate if you do not provide a trajectory. For each of these a distance scan will be performed." },
        { "-shelltoler", FALSE, etREAL, {&shellToler},
          "Tolerance for shell force optimization (mean square force)." },
        { "-json", FALSE, etBOOL, {&json},
          "Print part of the output in json format." }
    };
    CommunicationRecord cr;
    cr.init(cr.size());
    DimerGenerator      gendimers;
    gendimers.addOptions(&pa, &fnm, &desc);
    ReRunner            rerun(true);
    rerun.addOptions(&pa, &fnm);
    CompoundReader      compR;
    compR.addOptions(&pa, &fnm, &desc);
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, 
                           fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        status = 1;
        return status;
    }
    gendimers.finishOptions(fnm);
    if (!compR.optionsOK(fnm))
    {
        return 1;
    }
    ForceField        pd;
    try
    {
        readForceField(opt2fn("-ff", fnm.size(), fnm.data()), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    (void) pd.verifyCheckSum(stderr);
    FILE *logFile = nullptr;
    if (cr.isMaster())
    {
        const char *logFileName = opt2fn("-g", fnm.size(),fnm.data());
        logFile = gmx_ffopen(logFileName, "w");
        print_header(logFile, pa, fnm);
    }
    compR.setLogfile(logFile);
    auto  forceComp = new ForceComputer(shellToler, 100);
    
    JsonTree jtree("SecondVirialCoefficient");
    if (verbose)
    {
        forceFieldSummary(&jtree, &pd);
    }

    std::vector<ACTMol> actmols = compR.read(pd, forceComp);
    auto &actmol = actmols[0];
    ACTMessage imm = ACTMessage::OK;
    std::vector<gmx::RVec> coords = actmol.xOriginal();

    if (ACTMessage::OK == imm && status == 0)
    {
        if (debug)
        {
            actmol.topology()->dump(debug);
        }
        rerun.setFunctions(forceComp, &gendimers, oenv);
        rerun.runB2(&cr, logFile, &pd, &actmol, maxdimers, verbose, fnm);
    }
    if (json && cr.isMaster())
    {
        jtree.write("simulate.json", json);
    }
    else if (logFile)
    {
        jtree.fwrite(logFile, json);
    }
    if (logFile)
    {
        gmx_ffclose(logFile);
    }
    return status;
}

} // namespace alexandria
