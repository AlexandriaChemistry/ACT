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
#include "b2utils.h"

#include <cctype>
#include <cmath>
#include <cstdlib>

#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "alexandria/princ.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/vec.h"

namespace alexandria
{

enum class RotationAlgorithm 
    { 
        Cartesian, Polar
    };

std::map<std::string, RotationAlgorithm> stringToRotationAlgorithm = {
    { "Cartesian", RotationAlgorithm::Cartesian },
    { "Polar", RotationAlgorithm::Polar },
    { "cartesian", RotationAlgorithm::Cartesian },
    { "polar", RotationAlgorithm::Polar }
};

class Rotator
{
private:
    //! The rotation matrix
    matrix            A_;
    //! The average matrix after many calls to rotate
    matrix            Average_;
    //! The number of matrices added
    size_t            naver_  = 0;
    //! Rotation algorithm to use
    RotationAlgorithm rotalg_ = RotationAlgorithm::Cartesian; 
    //! \brief Reset the matrix to a unity matrix
    void resetMatrix()
    {
        clear_mat(A_);
        A_[XX][XX] = A_[YY][YY] = A_[ZZ][ZZ] = 1;
        clear_mat(Average_);
    }
    
    /*! \brief Do the actual rotation of input coordinates
     * \param[in] coords Input coordinates
     * \return the rotated coordinates
     */
    std::vector<gmx::RVec> rotate(const std::vector<gmx::RVec> &coords)
    {
        std::vector<gmx::RVec> newcoords;
        for(size_t i = 0; i < coords.size(); i++)
        {
            gmx::RVec newx;
            mvmul(A_, coords[i], newx);
            newcoords.push_back(newx);
        }
        m_add(A_, Average_, Average_);
        naver_ += 1;
        return newcoords;
    }
    
    std::vector<gmx::RVec> cartesian(double                        ralpha,
                                     double                        rbeta,
                                     double                        rgamma,
                                     const std::vector<gmx::RVec> &coords)
    {
        // Distribution is 0-1, multiply by two to get to 2*M_PI
        double alpha = ralpha * 2 * M_PI;
        double beta  = rbeta  * 2 * M_PI;
        double gamma = rgamma * 2 * M_PI;
        double cosa  = std::cos(alpha);
        double sina  = std::sin(alpha);
        double cosb  = std::cos(beta);
        double sinb  = std::sin(beta);
        double cosg  = std::cos(gamma);
        double sing  = std::sin(gamma);
        
        A_[0][0] = cosb * cosg;
        A_[0][1] =-cosb * sing;
        A_[0][2] = sinb;
        
        A_[1][0] = sina * sinb * cosg + cosa * sing;
        A_[1][1] =-sina * sinb * sing + cosa * cosg;
        A_[1][2] =-sina * cosb;
        
        A_[2][0] =-cosa * sinb * cosg + sina * sing;
        A_[2][1] = cosa * sinb * sing + sina * cosg;
        A_[2][2] = cosa * cosb;
        
        return rotate(coords);
    }
    
    std::vector<gmx::RVec>  polar(double                        rtheta,
                                  double                        rphi,
                                  double                        rgamma,
                                  const std::vector<gmx::RVec> &coords)
    {
        // Distribution is 0-1, multiply by two to get to 2*M_PI
        double    theta = std::acos(2*rtheta-1);
        double    phi   = rphi  * 2 * M_PI;
        // Create random vector
        // https://stackoverflow.com/questions/20769011/converting-3d-polar-coordinates-to-cartesian-coordinates
        gmx::RVec u     = { 
            std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi),
            std::cos(theta)
        };
        // Now create rotation matrix corresponding to rotation about this vector
        // https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
        // Confusing notation with two theta angles, Wikipedia is strange
        double gamma = rgamma * 2 * M_PI;
        double costh  = std::cos(gamma);
        double sinth  = std::sin(gamma);
        
        matrix B = {
            { costh + u[XX]*u[XX]*(1-costh),
              u[XX]*u[YY]*(1-costh) - u[ZZ]*sinth,
              u[XX]*u[ZZ]*(1-costh) + u[YY]*sinth },
            { u[YY]*u[XX]*(1-costh) + u[ZZ]*sinth,
              costh + u[YY]*u[YY]*(1-costh),
              u[YY]*u[ZZ]*(1-costh) - u[XX]*sinth },
            { u[ZZ]*u[XX]*(1-costh) - u[YY]*sinth,
              u[ZZ]*u[YY]*(1-costh) + u[XX]*sinth,
              costh + u[ZZ]*u[ZZ]*(1-costh) }
        };
        copy_mat(B, A_);
        return rotate(coords);
    }
    
public:
    Rotator(const std::string &rotalg)
    {
        resetMatrix();
        if (stringToRotationAlgorithm.end() != stringToRotationAlgorithm.find(rotalg))
        {
            rotalg_ = stringToRotationAlgorithm[rotalg];
        }
    }
    
    std::vector<gmx::RVec> random(double                        rtheta,
                                  double                        rphi,
                                  double                        rgamma,
                                  const std::vector<gmx::RVec> &coords)
    {
        std::vector<gmx::RVec> rx;
        switch(rotalg_)
        {
        case RotationAlgorithm::Cartesian:
            rx = cartesian(rtheta, rphi, rgamma, coords);
            break;
        case RotationAlgorithm::Polar:
            rx = polar(rtheta, rphi, rgamma, coords);
            break;
        }
        return rx;
    }
    void checkMatrix(FILE *fp)
    {
        fprintf(fp, "Norms of rows: %g %g %g\n",
                norm(A_[XX]), norm(A_[YY]), norm(A_[ZZ]));
        matrix B;
        transpose(A_, B);
        fprintf(fp, "Norms of columns: %g %g %g\n",
                norm(B[XX]), norm(B[YY]), norm(B[ZZ]));
    }
        
    void printAverageMatrix(FILE *fp)
    {
        if (fp && naver_ > 0)
        {
            fprintf(fp, "Average Matrix (n=%zu)\n", naver_);
            for(int m = 0; m < DIM; m++)
            {
                for(int n = 0; n < DIM; n++)
                {
                    fprintf(fp, "  %10g", Average_[m][n]/naver_);
                }
                fprintf(fp, "\n");
            }
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
          "Random number seed to generate monomer orientations, applied if seed is larger than 0. If not, the built-in default will be used." },
        { "-rotalg", FALSE, etSTR, {&rotalg_},
          "Rotation algorithm should be either Cartesian or Polar" }
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
    if (fragptr->topologies().size() != 2)
    {
        gmx_fatal(FARGS, "Cannot generate dimers, numer of compounds is %lu",
                  fragptr->topologies().size());
    }
    // Random number generation
    Rotator rot(rotalg_);
    
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
    // Then initiate the big array
    coords->resize(nmp);
    if (debug)
    {
        print_memory_usage(debug);
    }
    // Loop over the dimers
    for(int ndim = 0; ndim < maxdimers_; ndim++)
    {
        // Copy the coordinates and rotate them
        std::vector<gmx::RVec> xrand[2];
        for(int m = 0; m < 2; m++)
        {
            // Random rotation
            xrand[m] = rot.random(dis_(gen_), dis_(gen_), dis_(gen_), xmOrig[m]);
            if (ndim == 0)
            {
                rot.checkMatrix(logFile);
            }
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
        if (debug)
        {
            print_memory_usage(debug);
        }
    }
    rot.printAverageMatrix(logFile);
}

} // namespace alexandria
