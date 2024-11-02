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

std::map<std::string, RotationAlgorithm> stringToRotationAlgorithm = {
    { "Cartesian", RotationAlgorithm::Cartesian },
    { "Polar", RotationAlgorithm::Polar },
    { "Sobol", RotationAlgorithm::Sobol },
    { "cartesian", RotationAlgorithm::Cartesian },
    { "polar", RotationAlgorithm::Polar },
    { "sobol", RotationAlgorithm::Sobol },
};

const std::string &rotalgToString(RotationAlgorithm rotalg)
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

void Rotator::resetMatrix()
{
    clear_mat(A_);
    A_[XX][XX] = A_[YY][YY] = A_[ZZ][ZZ] = 1;
    clear_mat(Average_);
}
    
std::vector<gmx::RVec> Rotator::rotate(const std::vector<gmx::RVec> &coords)
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
    
void Rotator::storeAngles(double alpha, double beta, double gamma)
{
    if (debugAngles_)
    {
        alpha_.add_point(RAD2DEG*alpha);
        beta_.add_point(RAD2DEG*beta);
        gamma_.add_point(RAD2DEG*gamma);
    }
}
    
std::vector<gmx::RVec> Rotator::cartesian(double                        alpha,
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
    
std::vector<gmx::RVec> Rotator::polar(double                        phi,
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

std::vector<gmx::RVec> Rotator::sobol(double                        alpha,
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

void Rotator::printOneAngleHisto(gmx_stats angle, const char *file)
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
    
Rotator::Rotator(const std::string &rotalg, bool debugAngles)
{
    resetMatrix();
    if (stringToRotationAlgorithm.end() != stringToRotationAlgorithm.find(rotalg))
    {
        rotalg_ = stringToRotationAlgorithm[rotalg];
    }
    debugAngles_ = debugAngles;
}
    
std::vector<gmx::RVec> Rotator::random(double                        r1,
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
            rx = polar(alpha, beta, gamma, coords);
        }
        break;
    case RotationAlgorithm::Sobol:
        {
            double beta  = std::acos(2*r3-1);
            rx = sobol(alpha, beta, gamma, coords);
        }
        break;
    }
    return rx;
}
    
void Rotator::checkMatrix(FILE *fp)
{
    fprintf(fp, "Norms of rows: %g %g %g\n",
            norm(A_[XX]), norm(A_[YY]), norm(A_[ZZ]));
    matrix B;
    transpose(A_, B);
    fprintf(fp, "Norms of columns: %g %g %g\n",
            norm(B[XX]), norm(B[YY]), norm(B[ZZ]));
}
        
void Rotator::printAverageMatrix(FILE *fp)
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

void Rotator::printAngleHisto()
{
    printOneAngleHisto(alpha_, "alpha.xvg");
    printOneAngleHisto(beta_, "beta.xvg");
    printOneAngleHisto(gamma_, "gamma.xvg");
}

} // namespace alexandria

