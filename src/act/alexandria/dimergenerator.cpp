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
#include "dimergenerator.h"

#include <cctype>
#include <cmath>
#include <cstdlib>

#include "act/alexandria/rotator.h"
#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "external/quasirandom_sequences/sobol.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

void DimerGenerator::addOptions(std::vector<t_pargs>  *pa,
                                std::vector<t_filenm> *fnm)
{
    std::vector<t_pargs> mypa = {
        { "-ndist", FALSE, etINT, {&ndist_},
          "Number of equidistant distances to use for computing interaction energies and forces. Total number of dimers is the product of maxdimer and ndist. If left to zero, the first 0.5 nm will use 2.5 pm, and beyond that 7.5 pm." },
        { "-mindist", FALSE, etREAL, {&mindist_},
          "Minimum com-com distance to generate dimers for." },
        { "-maxdist", FALSE, etREAL, {&maxdist_},
          "Maximum com-com distance to generate dimers for." },
        { "-seed", FALSE, etINT, {&seed_},
          "Random number seed to generate monomer orientations, applied if seed is larger than 0. If not, the built-in default will be used. For the Sobol quasi-random sequence it is advised to leave the seed at 0." },
        { "-rotalg", FALSE, etSTR, {&rotalg_},
          "Rotation algorithm should be either Cartesian, Polar or Sobol. Default is Cartesian and the other two algorithms are experimental. Please verify your output when using those." },
        { "-dbgGD", FALSE, etBOOL, {&debugGD_},
          "Low-level debugging of routines. Gives complete information only when run on a single processor." }
    };
    for(auto &pp : mypa)
    {
        pa->push_back(pp);
    }
    std::vector<t_filenm>  myfnm = {
        { efSTX, "-ox", "dimers",     ffOPTWR  }
    };
    for(auto &ff : myfnm)
    {
        fnm->push_back(ff);
    }
}

DimerGenerator::~DimerGenerator()
{
    if (rot_)
    {
        rot_->printAngleHisto();
        if (debugGD_)
        {
            rot_->printAverageMatrix(stdout);
        }
        delete rot_;
    }
}

void DimerGenerator::setSeed(int seed)
{
    seed_      = seed;
    sobolSeed_ = seed_;
}

void DimerGenerator::finishOptions()
{
    if (seed_ > 0)
    {
        gen_.seed(seed_);
    }
    sobolSeed_ = seed_;
    rot_ = new Rotator(rotalg_, debugGD_);
    if (0 == ndist_)
    {
        ndist_ = 1+std::round((maxdist_-mindist_)/binWidth_);
    }
    else
    {
        binWidth_ = (maxdist_-mindist_)/(ndist_);
    }
}

static void dump_coords(const char                                *outcoords,
                        const ACTMol                              *actmol,
                        const std::vector<std::vector<gmx::RVec>> &coords)
{
    std::string outxyz(outcoords);
    auto pos = outxyz.rfind(".");
    outxyz = outxyz.substr(0, pos) + ".xyz";
    FILE *fp = gmx_ffopen(outxyz.c_str(), "w");
    auto atoms = actmol->atomsConst();
    for(size_t ix = 0; ix < coords.size(); ix++)
    {
        fprintf(fp, "%5lu\n", atoms.size());
        fprintf(fp, "Conformation %10zu\n", ix);
        for(size_t iy = 0; iy < atoms.size(); iy++)
        {
            fprintf(fp, "%5s  %10g  %10g  %10g\n", atoms[iy].name().c_str(),
                    10*coords[ix][iy][XX], 10*coords[ix][iy][YY], 10*coords[ix][iy][ZZ]);
        }
    }
    gmx_ffclose(fp);
}

void DimerGenerator::generate(FILE                                *logFile,
                              const ACTMol                        *actmol,
                              int                                  maxdimer,
                              std::vector<std::vector<gmx::RVec>> *coords,
                              const char                          *outcoords)
{
    size_t nmp  = ndist_*maxdimer;
    size_t mem  = (nmp*actmol->xOriginal().size()*sizeof(double)*DIM)/(1024*1024);
    auto   info = gmx::formatString("Will generate %d dimer configurations at %d distances using %s algorithm. Memory usage: %zu Mb",
                                    maxdimer, ndist_, rotalgToString(rot_->rotalg()).c_str(), mem);
    if (logFile)
    {
        fprintf(logFile, "%s\n", info.c_str());
    }
    printf("%s\n", info.c_str());
    if (!rot_)
    {
        GMX_THROW(gmx::InternalError("Forgot to call DimerGenerator::finishOptions"));
    }
    // Loop over the dimers
    for(int ndim = 0; ndim < maxdimer; ndim++)
    {
        for(auto &newx : generateDimers(actmol))
        {
            coords->push_back(newx);
        }
    }
    if (nullptr != outcoords && strlen(outcoords) > 0)
    {
        dump_coords(outcoords, actmol, *coords);
    }

}

std::vector<std::vector<gmx::RVec>> DimerGenerator::generateDimers(const ACTMol *actmol)
{
    auto fraghandler = actmol->fragmentHandlerConst();
    if (fraghandler.topologies().size() != 2)
    {
        gmx_fatal(FARGS, "Cannot generate dimers, numer of compounds is %lu",
                  fraghandler.topologies().size());
    }

    // Copy original coordinates
    auto xorig     = actmol->xOriginal();
    // Split the coordinates into two fragments
    auto atomStart = fraghandler.atomStart();
    // Topologies
    auto tops = fraghandler.topologies();
    std::vector<gmx::RVec> xmOrig[2];
    for(int m = 0; m < 2; m++)
    {
        auto   atoms   = tops[m].atoms();
        for(size_t j = atomStart[m]; j < atomStart[m]+atoms.size(); j++)
        {
            xmOrig[m].push_back(xorig[j]);
        } 
    }
    // Move molecules to their respective COM
    gmx::RVec com[2];
    for(int m = 0; m < 2; m++)
    {
        // Compute center of mass
        clear_rvec(com[m]);
        auto   atoms   = tops[m].atoms();
        double totmass = 0;
        for(size_t j = 0; j < atoms.size(); j++)
        {
            gmx::RVec mx;
            svmul(atoms[j].mass(), xmOrig[m][j], mx);
            rvec_inc(com[m], mx);
            totmass += atoms[j].mass();
        }
        for(int n = 0; n < DIM; n++)
        {
            com[m][n] /= totmass;
        }
        // Subtract center of mass
        for(size_t j = 0; j < atoms.size(); j++)
        {
            rvec_sub(xmOrig[m][j], com[m], xmOrig[m][j]);
        }
    }
    // Loop over orientations
    // Then initiate the big array
    std::vector<std::vector<gmx::RVec>> coords(ndist_);
    if (debugGD_ && debug)
    {
        print_memory_usage(debug);
    }
    // Copy the original coordinates
    std::vector<gmx::RVec> xrand[2];
    for(int m = 0; m < 2; m++)
    {
        xrand[m] = xmOrig[m];
    }
    std::vector<double> q(2*DIM, 0.0);
    if (RotationAlgorithm::Sobol == rot_->rotalg())
    {
        // Quasi random numbers
        i8_sobol(2*DIM, &sobolSeed_, q.data());
    }
    else
    {
        for(size_t j = 0; j < 2*DIM; j++)
        {
            // "True" random numbers
            q[j] = dis_(gen_);
        }
    }
    // Rotate the coordinates
    for(int m = 0; m < 2; m++)
    {
        // Random rotation, using the pre-calculated random numbers
        size_t j0 = m*DIM;
        xrand[m] = rot_->random(q[j0], q[j0+1], q[j0+2], xrand[m]);
    }
    // Loop over distances from mindist to maxdist
    for(int idist = 0; idist < ndist_; idist++)
    {
        double    dist  = mindist_;
        if (ndist_ > 1)
        {
            dist += idist*binWidth_;
        }
        gmx::RVec trans = { 0, 0, dist };
        auto      atoms = tops[1].atoms();
        for(size_t j = 0; j < atoms.size(); j++)
        {
            rvec_inc(xrand[1][j], trans);
        }
        coords[idist].resize(xorig.size());
        for(int m = 0; m < 2; m++)
        {
            auto   atoms   = tops[m].atoms();
            for(size_t j = atomStart[m]; j < atomStart[m]+atoms.size(); j++)
            {
                auto jindex = j - atomStart[m];
                copy_rvec(xrand[m][jindex], coords[idist][j]);
            }
        }
            // Put the coordinates back!
        for(size_t j = 0; j < atoms.size(); j++)
        {
            rvec_dec(xrand[1][j], trans);
        }
    }
    if (debugGD_ && debug)
    {
        print_memory_usage(debug);
    }
    return coords;
}

} // namespace alexandria
