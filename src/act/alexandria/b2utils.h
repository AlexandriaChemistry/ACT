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

#ifndef ACT_B2UTILS_H
#define ACT_B2UTILS_H

#include <random>
#include <vector>

#include "act/alexandria/mymol.h"
#include "act/forces/forcecomputer.h"
#include "act/poldata/poldata.h"
#include "act/utility/jsontree.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/statistics/statistics.h"

namespace alexandria
{

class DimerGenerator
{
private:
    //! Number of dimers to generate
    int maxdimers_      = 1000;
    //! Number of distances to use. See alexandria b2 -h for an explanation
    int ndist_          =    0;
    //! Minimum com-com distance (nm)
    double mindist_     =  0.1;
    //! Maximum com-com distance (nm)
    double maxdist_     =  4.0;
    //! (Quasi-) random number seed
    int    seed_        =    0;
    //! Low-level debugging
    bool   debugGD_     = false;
    //! Rotation algorithm
    const char *rotalg_ = "";
    //! Random device
    std::random_device                     rd_;
    //! Generator
    std::mt19937                           gen_;
    //! Distribution
    std::uniform_real_distribution<double> dis_;
public:
    //! Constructor
    DimerGenerator() : gen_(rd_()), dis_(std::uniform_real_distribution<double>(0.0, 1.0)) {}
    
    /*! \brief Add my options
     * \param[inout] pa  The command line options
     * \param[inout] fnm File names for the command line
     */
    void addOptions(std::vector<t_pargs>  *pa,
                    std::vector<t_filenm> *fnm);

    //! \brief Process the options
    void finishOptions();
    
    //! Return the number of distancea
    int ndist() const { return ndist_; }
    
    /*! \brief Do the actual generation
     * \param[in]  logFile   For debugging info, may be a nullptr
     * \param[in]  mymol     The description of the two fragments
     * \param[out] coords    The coordinate sets
     * \param[in]  outcoords Output file name for the generated coordinates.
     *                       If nullptr, they will not be written.
     */
    void generate(FILE                                *logFile,
                  const MyMol                         *mymol,
                  std::vector<std::vector<gmx::RVec>> *coords,
                  const char                          *outcoords);
};

} // namespace alexandria

#endif // ACT_B2UTILS_H
