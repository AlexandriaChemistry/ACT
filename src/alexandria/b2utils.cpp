/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022,2023
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

#include <ctype.h>
#include <stdlib.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"

#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "alexandria/princ.h"

namespace alexandria
{

class Rotator
{
public:
    Rotator() {}
    
    void random(double                  ralpha,
                double                  rbeta,
                double                  rgamma,
                std::vector<gmx::RVec> *coords)
    {
        // Distribution is 0-M_PI, multiply by two to get to 2*M_PI
        double alpha = ralpha * 2 * M_PI;
        double beta  = rbeta  * 2 * M_PI; 
        double gamma = std::acos(2*rgamma-1);
        double cosa  = std::cos(alpha);
        double sina  = std::sin(alpha);
        double cosb  = std::cos(beta);
        double sinb  = std::sin(beta);
        double cosg  = std::cos(gamma);
        double sing  = std::sin(gamma);
        
        matrix A;
        A[0][0] = cosb * cosg;
        A[0][1] =-cosb * sing;
        A[0][2] = sinb;
        
        A[1][0] = sina * sinb * cosg + cosa * sing;
        A[1][1] =-sina * sinb * sing + cosa * cosg;
        A[1][2] =-sina * cosb;
        
        A[2][0] =-cosa * sinb * cosg + sina * sing;
        A[2][1] = cosa * sinb * sing + sina * cosg;
        A[2][2] = cosa * cosb;
        
        auto oldX = *coords;
        for(size_t i = 0; i < oldX.size(); i++)
        {
            gmx::RVec newx;
            mvmul(A, oldX[i], newx);
            copy_rvec(newx, (*coords)[i]);
        }
    }
};

void DimerGenerator::addOptions(std::vector<t_pargs>  *pa,
                                std::vector<t_filenm> *fnm)
{
    std::vector<t_pargs> mypa = {
        { "-maxdimer", FALSE, etINT, {&maxdimers_},
          "Number of dimer orientations to generate if you do not provide a trajectory. For each of these a distance scan will be performed." },
        { "-ndist", FALSE, etINT, {&ndist_},
          "Number of equidistant distances to use for computing interaction energies and forces. Total number of dimers is the product of maxdimer and ndist. If left to zero, the first 0.5 nm will use 2.5 pm, and beyond that 7.5 pm." },
        { "-mindist", FALSE, etREAL, {&mindist_},
          "Minimum com-com distance to generate dimers for." },
        { "-maxdist", FALSE, etREAL, {&maxdist_},
          "Maximum com-com distance to generate dimers for." },
        { "-seed", FALSE, etINT, {&seed_},
          "Random number seed to generate monomer orientations, applied if seed is larger than 0. If not, the built-in default will be used." }
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

void DimerGenerator::finishOptions()
{
    if (seed_ > 0)
    {
        gen_.seed(seed_);
    }
}
    
void DimerGenerator::generate(FILE                                *logFile,
                              const MyMol                         *mymol,
                              std::vector<std::vector<gmx::RVec>> *coords,
                              gmx_unused const char               *outcoords)
{
    auto fragptr = mymol->fragmentHandler();
    if (fragptr->topologies().size() == 2)
    {
        // Random number generation
        Rotator rot;
        
        // Copy original coordinates
        auto xorig     = mymol->xOriginal();
        // Split the coordinates into two fragments
        auto atomStart = fragptr->atomStart();
        std::vector<gmx::RVec> xmOrig[2];
        for(int m = 0; m < 2; m++)
        {
            for(size_t j = atomStart[m]; j < atomStart[m+1]; j++)
            {
                xmOrig[m].push_back(xorig[j]);
            } 
        }
        // Topologies
        auto tops = fragptr->topologies();
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
        size_t nmp = maxdimers_*ndist_;
        // Initiate all the MolProps with a copy of the input molecule
        MolProp tmp = *mymol;
        auto exper  = tmp.experiment();
        // Do some cleaning
        exper->clear();
        tmp.clearCategory();
        // Then initiate the big array
        coords->resize(nmp);
        if (logFile)
        {
            print_memory_usage(logFile);
        }
        // Shortcut to the current molecules
        // MolProp *mp = &((*mps)[mp_index++]);
        // exper = mp->experiment();
        for(int ndim = 0; ndim < maxdimers_; ndim++)
        {
            // Copy the coordinates and rotate them
            std::vector<gmx::RVec> xrand[2];
            for(int m = 0; m < 2; m++)
            {
                xrand[m] = xmOrig[m];
                // Random rotation
                rot.random(dis_(gen_), dis_(gen_), dis_(gen_), &xrand[m]);
            }
            // Loop over distances from mindist to maxdist
            double range = maxdist_ - mindist_;
            for(int idist = 0; idist < ndist_; idist++)
            {
                double    dist  = mindist_ + range*((1.0*idist)/(ndist_-1));
                gmx::RVec trans = { 0, 0, dist };
                auto      atoms = tops[1].atoms();
                for(size_t j = 0; j < atoms.size(); j++)
                {
                    rvec_inc(xrand[1][j], trans);
                }
                size_t idim = ndim*ndist_+idist;
                (*coords)[idim].resize(atomStart[2]-atomStart[0]);
                for(int m = 0; m < 2; m++)
                {
                    for(size_t j = atomStart[m]; j < atomStart[m+1]; j++)
                    {
                        auto jindex = j - atomStart[m];
                        copy_rvec(xrand[m][jindex], (*coords)[idim][j]);
                    }
                }
                // Put the coordinates back!
                for(size_t j = 0; j < atoms.size(); j++)
                {
                    rvec_dec(xrand[1][j], trans);
                }
            }
            if (logFile)
            {
                print_memory_usage(logFile);
            }
        }
    }
}

} // namespace alexandria
