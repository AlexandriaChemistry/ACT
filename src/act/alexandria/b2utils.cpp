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
#include "external/quasirandom_sequences/sobol.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

enum class RotationAlgorithm 
    { 
        Cartesian, Polar, Sobol
    };

std::map<std::string, RotationAlgorithm> stringToRotationAlgorithm = {
    { "Cartesian", RotationAlgorithm::Cartesian },
    { "Polar", RotationAlgorithm::Polar },
    { "Sobol", RotationAlgorithm::Sobol },
    { "cartesian", RotationAlgorithm::Cartesian },
    { "polar", RotationAlgorithm::Polar },
    { "sobol", RotationAlgorithm::Sobol },
};

static const std::string &rotalgToString(RotationAlgorithm rotalg)
{
    for(auto &stra : stringToRotationAlgorithm)
    {
        if (stra.second == rotalg)
        {
            return stra.first;
        }
    }
    return stringToRotationAlgorithm.begin()->first;
}

class Rotator
{
private:
    //! The rotation matrix
    matrix            A_;
    //! The average matrix after many calls to rotate
    matrix            Average_;
    //! The number of matrices added
    size_t            naver_  = 0;
    //! Sobol seed, only use this for debugging
    long long int     sobolSeed_ = 0;
    //! Rotation algorithm to use
    RotationAlgorithm rotalg_ = RotationAlgorithm::Cartesian;
    //! Debug angles?
    bool              debugAngles_ = false;
    //! Statistics of angles used
    gmx_stats         alpha_, beta_, gamma_;
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
    
    /*! \brief Store the angles generated if requested
     * \param[in] alpha First angle, unit radians
     * \param[in] beta  Second angle
     * \param[in] gamma Third angle
     */
    void storeAngles(double alpha, double beta, double gamma)
    {
        if (debugAngles_)
        {
            alpha_.add_point(RAD2DEG*alpha);
            beta_.add_point(RAD2DEG*beta);
            gamma_.add_point(RAD2DEG*gamma);
        }
    }
    
    std::vector<gmx::RVec> cartesian(double                        alpha,
                                     double                        beta,
                                     double                        gamma,
                                     const std::vector<gmx::RVec> &coords)
    {
        storeAngles(alpha, beta, gamma);
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
    
    std::vector<gmx::RVec> polar(double                        phi,
                                 double                        theta,
                                 double                        gamma,
                                 const std::vector<gmx::RVec> &coords)
    {
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
        double costh  = std::cos(gamma);
        double sinth  = std::sin(gamma);
        storeAngles(phi, theta, gamma);
        
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

    /*! \brief Rotate using angles generated by a Sobol sequence
     * \param[in] coords Input coordinates 
     * \returns the rotated coordinates 
     */  
    std::vector<gmx::RVec> sobol(double                        alpha,
                                 double                        beta,
                                 double                        gamma,
                                 const std::vector<gmx::RVec> &coords)
    {
        // Orientation described by Euler angles
        storeAngles(alpha, beta, gamma);
        double cosa = std::cos(alpha);
        double sina = std::sin(alpha);
        double cosb = std::cos(beta);
        double sinb = std::sin(beta);
        double cosc = std::cos(gamma);
        double sinc = std::sin(gamma);
        A_[XX][XX] =  cosa*cosb*cosc-sina*sinc;
        A_[YY][XX] =  sina*cosb*cosc+cosa*sinc;
        A_[ZZ][XX] = -sinb*cosc;
        A_[XX][YY] = -cosa*cosb*sinc-sina*cosc;
        A_[YY][YY] = -sina*cosb*sinc+cosa*cosc;
        A_[ZZ][YY] =  sinb*sinc;
        A_[XX][ZZ] =  cosa*sinb;
        A_[YY][ZZ] =  sina*sinb;
        A_[ZZ][ZZ] =  cosb;

        return rotate(coords);
    }
    /*! \brief Print a histogram of an angle
     * \param[in] angle The statistics container
     * \param[in] file  The filename to print to 
     */
    void printOneAngleHisto(gmx_stats angle, const char *file)
    {
        if (angle.get_npoints() == 0)
        {
            return;
        }
        real binwidth   = 2;
        int  nbins      = 0;
        bool normalized = true;
        std::vector<double> xx, yy;
        if (eStats::OK == angle.make_histogram(binwidth, &nbins, eHisto::Y,
                                               normalized, &xx, &yy))
        {
            FILE *fp = gmx_ffopen(file, "w");
            for(size_t i = 0; i < yy.size(); i++)
            {
                fprintf(fp, "%10g  %10g\n", xx[i], yy[i]);
            }
            gmx_ffclose(fp);
        }
    }
    
public:
    /*! \brief Constructor setting up algorithm
     * \param[in] rotalg      The rotation algorithm string
     * \param[in] debugAngles Whether or not to print histograms of angles
     * \param[in] sobolSeed   Starting index in the Sobol sequence. Keep at 
     *                        zero except for debugging.
     */
    Rotator(const std::string &rotalg, bool debugAngles, long long int sobolSeed = 0)
    {
        resetMatrix();
        if (stringToRotationAlgorithm.end() != stringToRotationAlgorithm.find(rotalg))
        {
            rotalg_ = stringToRotationAlgorithm[rotalg];
        }
        debugAngles_ = debugAngles;
        sobolSeed_   = sobolSeed;
    }
    
    //! \return the rotation algorithm selected
    RotationAlgorithm rotalg() const { return rotalg_; }
    
    /*! \brief Do a (quasi) random rotation
     */
    std::vector<gmx::RVec> random(double                        r1,
                                  double                        r2,
                                  double                        r3,
                                  const std::vector<gmx::RVec> &coords)
    {
        // Distribution is 0-1, multiply by two to get to 2*M_PI
        double alpha = r1 * 2 * M_PI;
        double gamma = r2 * 2 * M_PI;
        std::vector<gmx::RVec> rx;
        switch(rotalg_)
        {
        case RotationAlgorithm::Cartesian:
            {
                double beta  = r3 * 2 * M_PI;
                rx = cartesian(alpha, beta, gamma, coords);
            }
            break;
        case RotationAlgorithm::Polar:
            {
                double beta  = std::acos(2*r3-1);
                rx = sobol(alpha, beta, gamma, coords);
            }
            break;
        case RotationAlgorithm::Sobol:
            {
                double beta  = std::acos(2*r3-1);
                rx = polar(alpha, beta, gamma, coords);
            }
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
    void printAngleHisto()
    {
        printOneAngleHisto(alpha_, "alpha.xvg");
        printOneAngleHisto(beta_, "beta.xvg");
        printOneAngleHisto(gamma_, "gamma.xvg");
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
          "Random number seed to generate monomer orientations, applied if seed is larger than 0. If not, the built-in default will be used. For the Sobol quasi-random sequence it is advised to leave the seed at 0." },
        { "-rotalg", FALSE, etSTR, {&rotalg_},
          "Rotation algorithm should be either Cartesian, Polar or Sobol" },
        { "-dbgGD", FALSE, etBOOL, {&debugGD_},
          "Low-level debugging of routines" }
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
                              const ACTMol                         *actmol,
                              std::vector<std::vector<gmx::RVec>> *coords,
                              gmx_unused const char               *outcoords)
{
    auto fragptr = actmol->fragmentHandler();
    if (fragptr->topologies().size() != 2)
    {
        gmx_fatal(FARGS, "Cannot generate dimers, numer of compounds is %lu",
                  fragptr->topologies().size());
    }
    // Random number generation
    long long int sobolSeed = seed_;
    Rotator rot(rotalg_, debugGD_, sobolSeed);
    size_t nmp = maxdimers_*ndist_;
    
    auto info = gmx::formatString("Will generate %zu dimer configurations at %d distances using %s algorithm.",
                                  nmp, ndist_, rotalgToString(rot.rotalg()).c_str());
    if (logFile)
    {
        fprintf(logFile, "%s\n", info.c_str());
    }
    printf("%s\n", info.c_str());
    
    // Copy original coordinates
    auto xorig     = actmol->xOriginal();
    // Split the coordinates into two fragments
    auto atomStart = fragptr->atomStart();
    // Topologies
    auto tops = fragptr->topologies();
    std::vector<gmx::RVec> xmOrig[2];
    for(int m = 0; m < 2; m++)
    {
        auto   atoms   = tops[m]->atoms();
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
        auto   atoms   = tops[m]->atoms();
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
    coords->resize(nmp);
    if (debugGD_ && debug)
    {
        print_memory_usage(debug);
    }
    // Loop over the dimers
    for(int ndim = 0; ndim < maxdimers_; ndim++)
    {
        // Copy the original coordinates
        std::vector<gmx::RVec> xrand[2];
        for(int m = 0; m < 2; m++)
        {
            xrand[m] = xmOrig[m];
        }
        std::vector<double> q(2*DIM, 0.0);
        if (RotationAlgorithm::Sobol == rot.rotalg())
        {
            // Quasi random numbers
            i8_sobol(2*DIM, &sobolSeed, q.data());
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
            xrand[m] = rot.random(q[j0], q[j0+1], q[j0+2], xrand[m]);
            if (debugGD_ && ndim == 0)
            {
                rot.checkMatrix(logFile);
            }
        }
        // Loop over distances from mindist to maxdist
        double range = maxdist_ - mindist_;
        for(int idist = 0; idist < ndist_; idist++)
        {
            double    dist  = mindist_;
            if (ndist_ > 1)
            {
                dist += range*((1.0*idist)/(ndist_-1));
            }
            gmx::RVec trans = { 0, 0, dist };
            auto      atoms = tops[1]->atoms();
            for(size_t j = 0; j < atoms.size(); j++)
            {
                rvec_inc(xrand[1][j], trans);
            }
            size_t idim = ndim*ndist_+idist;
            (*coords)[idim].resize(xorig.size());
            for(int m = 0; m < 2; m++)
            {
                auto   atoms   = tops[m]->atoms();
                for(size_t j = atomStart[m]; j < atomStart[m]+atoms.size(); j++)
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
        if (debugGD_ && debug)
        {
            print_memory_usage(debug);
        }
    }
    rot.printAngleHisto();
    if (debugGD_)
    {
        rot.printAverageMatrix(logFile);
    }
    if (nullptr != outcoords && strlen(outcoords) > 0)
    {
        dump_coords(outcoords, actmol, *coords);
    }
}

} // namespace alexandria
