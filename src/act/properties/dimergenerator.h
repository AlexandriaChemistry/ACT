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

#ifndef ACT_DIMERGENERATOR_H
#define ACT_DIMERGENERATOR_H

#include <random>
#include <vector>

#include "act/alexandria/actmol.h"
#include "act/forces/forcecomputer.h"
#include "act/forcefield/forcefield.h"
#include "act/properties/rotator.h"
#include "act/utility/jsontree.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class DimerGenerator
{
private:
    //! Number of distances to use. See alexandria b2 -h for an explanation
    int ndist_          =  0;
    //! The bin width (nm)
    double binWidth_    =  0.0025;
    //! Minimum com-com distance (nm)
    double mindist_     =  0.1;
    //! Maximum com-com distance (nm)
    double maxdist_     =  4.0;
    //! Random number seed, not applied to the Sobol algorithm
    int    dimerseed_   =  1993;
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
    //! Rotator
    Rotator                               *rot_ = nullptr;
    //! Trajectory name
    const char                            *trajname_ = nullptr;
    //! The random numbers
    std::vector<std::vector<double>>       allRandom_;
    //! Random number index
    size_t                                 randIndex_ = 0;
public:
    //! Constructor
    DimerGenerator() : gen_(rd_()), dis_(std::uniform_real_distribution<double>(0.0, 1.0)) {}

    //! Destructor
    ~DimerGenerator();
  
    /*! \brief Add my options
     * \param[inout] pa   The command line options
     * \param[inout] fnm  File names for the command line
     * \param[inout] desc Descriptive text
     */
    void addOptions(std::vector<t_pargs>      *pa,
                    std::vector<t_filenm>     *fnm,
                    std::vector<const char *> *desc);

    /*! \brief Process the options
     * \param[in] fnm The file names after processing
     */
    void finishOptions(const std::vector<t_filenm> &fnm);

    //! Get the original seed
    int seed() const { return dimerseed_; }

    //! Set a new seed
    void setSeed(int seed);

    /*! Generate the random numbers for all rotations
     * \param[in] ndimers The number of dimers
     */
    void generateRandomNumbers(int ndimers);

    //! \return The list of random numbers
    const std::vector<std::vector<double>> &allRandom() const { return allRandom_; }

    //! Empty random number array
    void clearAllRandom() { allRandom_.clear(); }

    //! Empty random number array
    void resetAllRandom(size_t ndimer) { allRandom_.resize(ndimer); }

    /*! \brief Set the random numbers
     * \param[in] newRand The new random numbers
     */
    void addRandomNumbers(const std::vector<double> &newRand)
    {
        allRandom_.push_back(newRand);
    }

    //! Return the number of distancea
    int ndist() const { return ndist_; }

    //! Return trajectory name
    const char *trajname() const { return trajname_; }

    //! Return max distance
    double maxdist() const { return maxdist_; }

    //! Return bin width
    double binwidth() const { return binWidth_; }

    //! Return whether or not there is a trajectory to read
    bool hasTrajectory() const { return trajname_ && strlen(trajname_) > 0; }

    /*! \brief Do the actual generation for one pair. Two molecules
     * will be oriented and a distance series will be generated.
     * \param[in]  logFile   For debugging info, may be a nullptr
     * \param[in]  actmol  The description of the two fragments
     * \return The coordinate sets
     */
    std::vector<std::vector<gmx::RVec>> generateDimers(FILE         *logFile,
                                                       const ACTMol *actmol);
    /*! \brief Do the actual generation for many dimers.
     * Note that the memory usage can be significant (many Gb).
     * 
     * \param[in]  logFile   For debugging info, may be a nullptr
     * \param[in]  actmol    The description of the two fragments
     * \param[in]  maxdimer  The number of different dimers to generate
     * \param[out] coords    The coordinate sets
     * \param[in]  outcoords Output file name for the generated coordinates.
     *                       If nullptr, they will not be written.
     */
    void generate(FILE                                *logFile,
                  const ACTMol                        *actmol,
                  int                                  maxdimer,
                  std::vector<std::vector<gmx::RVec>> *coords,
                  const char                          *outcoords);

    /*! \brief Read all the dimers at once from a file. 
     * \param[out]   coords   The coordinates. If empty reading failed or the variable was empty
     */
    void read(std::vector<std::vector<gmx::RVec>> *coords);
};

} // namespace alexandria

#endif // ACT_DIMERGENERATOR_H
