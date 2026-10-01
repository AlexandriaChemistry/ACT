/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024,2026
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

#ifndef ACT_ROTATOR_H
#define ACT_ROTATOR_H

#include <random>
#include <vector>

#include "act/basics/msg_handler.h"
#include "act/statistics/statistics.h"
#include "gromacs/math/vec.h"

namespace alexandria
{

class Rotator
{
private:
    //! The rotation matrix
    matrix            A_;
    //! The average matrix after many calls to rotate
    matrix            Average_;
    //! The number of matrices added
    size_t            naver_  = 0;
    //! Debug angles?
    bool              debugAngles_ = false;
    //! Statistics of angles used
    gmx_stats         alpha_, beta_, gamma_;
    //! \brief Reset the matrix to a unity matrix
    void resetMatrix();
    
    /*! \brief Do the actual rotation of input coordinates
     * \param[in] coords Input coordinates
     * \return the rotated coordinates
     */
    std::vector<gmx::RVec> doRotate(const std::vector<gmx::RVec> &coords);
    
    /*! \brief Store the angles generated if requested
     * \param[in] alpha First angle, unit radians
     * \param[in] beta  Second angle
     * \param[in] gamma Third angle
     */
    void storeAngles(double alpha, double beta, double gamma);
    
    /*! \brief Print a histogram of an angle
     * \param[in] angle The statistics container
     * \param[in] file  The filename to print to 
     */
    void printOneAngleHisto(gmx_stats angle, const char *file);

public:
    /*! \brief Constructor setting up algorithm
     * \param[in] debugAngles Whether or not to print histograms of angles
     */
    Rotator(bool debugAngles);
    
    /*! \brief Do a (quasi) random rotation
     * All random numbers should be between 0 and 1.
     * Rotation is about the origin, so if molecules are not centered
     * in their center of mass (c.o.m.) the c.o.m. will move as well.
     * \param[in] r1     Random number corresponding to first angle
     * \param[in] r2     Second
     * \param[in] r3     Third
     * \param[in] coords Input coordinates 
     * \returns the rotated coordinates 
     */
    std::vector<gmx::RVec> randomRotate(double                        r1,
                                        double                        r2,
                                        double                        r3,
                                        const std::vector<gmx::RVec> &coords);

    void checkMatrix(MsgHandler *msghandler);

    void printAverageMatrix(MsgHandler *msghandler);

    void printAngleHisto();
};

}

#endif // ACT_ROTATOR_H
