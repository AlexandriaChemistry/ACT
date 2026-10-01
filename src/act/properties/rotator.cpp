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
#include "rotator.h"

#include <cctype>
#include <cmath>
#include <cstdlib>

#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "external/quasirandom_sequences/sobol.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

void Rotator::resetMatrix()
{
    clear_mat(A_);
    A_[XX][XX] = A_[YY][YY] = A_[ZZ][ZZ] = 1;
    clear_mat(Average_);
}
    
std::vector<gmx::RVec> Rotator::doRotate(const std::vector<gmx::RVec> &coords)
{
    std::vector<gmx::RVec> newcoords(coords.size());
    for(size_t i = 0; i < coords.size(); i++)
    {
        mvmul(A_, coords[i], newcoords[i]);
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
    
Rotator::Rotator(bool debugAngles)
{
    resetMatrix();
    debugAngles_ = debugAngles;
}
    
std::vector<gmx::RVec> Rotator::randomRotate(double                        r1,
                                             double                        r2,
                                             double                        r3,
                                             const std::vector<gmx::RVec> &coords)
{
    // Distribution is 0-1, multiply by two to get to 2*M_PI
    double alpha = r1 * 2 * M_PI;
    double gamma = r3 * 2 * M_PI;
    // Azimuthal angle to generate even sampling on a sphere
    double beta  = std::acos(2*r2-1);

    // Orientation is described by Euler angles
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
    
    return doRotate(coords);
}
    
void Rotator::checkMatrix(MsgHandler *msghandler)
{
    if (msghandler)
    {
        msghandler->writeDebug(gmx::formatString("Norms of rows: %g %g %g",
                                                  norm(A_[XX]), norm(A_[YY]), norm(A_[ZZ])));
        matrix B;
        transpose(A_, B);
        msghandler->writeDebug(gmx::formatString("Norms of columns: %g %g %g",
                                                  norm(B[XX]), norm(B[YY]), norm(B[ZZ])));
    }
}
        
void Rotator::printAverageMatrix(MsgHandler *msghandler)
{
    if (msghandler && naver_ > 0)
    {
        msghandler->writeDebug(gmx::formatString("Average Matrix (n=%zu)", naver_));
        double inv = 1.0 / naver_;
        for(int m = 0; m < DIM; m++)
        {
            msghandler->writeDebug(gmx::formatString("  %10g  %10g  %10g",
                                                      Average_[m][0]*inv,
                                                      Average_[m][1]*inv,
                                                      Average_[m][2]*inv));
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

