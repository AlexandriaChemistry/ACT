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

#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/properties/rotator.h"
#include "act/topology/actatom.h"
#include "act/topology/actmol.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

namespace alexandria
{

/*! \brief To distinguish algorithms for spherical sampling
 * Sobol and SobolOwen use quasi-random numbers, Pseudo pseudo-random
 * numbers. Gauss uses a grid of Euler angles, as in the program Orient,
 * see DimerGenerator::generateRandomNumbers.
 */
enum class RandAlg { Sobol, SobolOwen, Pseudo, Gauss };

/*! Convert string to corresponding algorithm
 * \param[in] str The string from the commandline
 * \return a RandAlg, or RandAlg::Sobol in case of incorrect string. Will print a message in that case.
 */
RandAlg stringToRandAlg(const char *str);

//! The coordinates of a dimer at one distance
struct dimer_coords
{
    //! Dstance between center of mass of compounds
    double                 dist;
    //! Index in order
    size_t                 index;
    //! Number of particles in compounds
    size_t                 natom[2];
    //! The actual coordinates as a single array for the whole dimer
    std::vector<gmx::RVec> coords;
    //! The weight of the orientation in averages over orientations
    double                 weight = 1;
};

/*! \brief Generate dimer coordinates from two monomers by sampling
 * relative orientations (see RandAlg) and series of center of mass
 * distances, or read them from a trajectory
 */
class DimerGenerator
{
private:
    //! Number of distances to use. See alexandria b2 -h for an explanation
    int ndist_          =  0;
    //! The bin width (nm)
    double binWidth_    =  0.0025;
    //! Minimum com-com distance (nm)
    real   mindist_     =  0.1;
    //! Maximum com-com distance (nm)
    real   maxdist_     =  4.0;
    //! Random number seed, not applied to the Sobol algorithm
    int    dimerseed_   =  1993;
    //! Whether to use flexible monomers
    bool   flexible_    = false;
    //! Whether to use minimize the monomers before generating dimers
    bool   minimize_    = true;
    //! Low-level debugging
    bool   debugGD_     = false;
    //! The random algorithm
    RandAlg                                randAlg_ = RandAlg::Sobol;
    //! Number of points for alpha in the Gauss algorithm
    int                                    nalpha_  = 7;
    //! Number of points for beta in the Gauss algorithm
    int                                    nbeta_   = 7;
    //! Number of points for gamma in the Gauss algorithm
    int                                    ngamma_  = 7;
    //! Generator, seeded from a random device unless a seed is given
    std::mt19937                           gen_;
    //! Distribution
    std::uniform_real_distribution<double> dis_;
    //! Rotator
    Rotator                                rot_ { false };
    //! Trajectory to read the dimers from, none if empty (-traj)
    std::string                            trajname_;
    //! The random numbers
    std::vector<std::vector<double>>       allRandom_;
    //! Random number index
    size_t                                 randIndex_ = 0;
    //! Coordinates for the two monomers
    std::vector<gmx::RVec>                 xmOrig_[2];
    //! Atom info for both monomers
    std::vector<ActAtom>                   atoms_[2];
    //! Output coordinate file name, none if empty (-ox)
    std::string                            outcoords_;
    /*! \brief The choice of -randalg as the command line parser keeps it
     * The first element points to the choice, nullptr before a command
     * line was parsed. Each generator has its own.
     */
    std::array<const char *, 6>            randAlgChoice_ = { nullptr, "Sobol", "SobolOwen", "Pseudo", "Gauss", nullptr };
public:
    //! Constructor
    DimerGenerator() : gen_(std::random_device{}()), dis_(std::uniform_real_distribution<double>(0.0, 1.0)) {}
  
    /*! \brief Add my options
     * \param[inout] pa   The command line options
     * \param[inout] fnm  File names for the command line
     * \param[inout] desc Descriptive text
     */
    void addOptions(std::vector<t_pargs>      *pa,
                    std::vector<t_filenm>     *fnm,
                    std::vector<const char *> *desc);

    /*! \brief Process the options
     * Takes the file names and the choice of the algorithm from the
     * command line, then see finishSettings.
     * \param[in] fnm The file names after processing
     * \throws gmx::InvalidInputError for settings that make no sense
     */
    void finishOptions(const std::vector<t_filenm> &fnm);

    /*! \brief Check the settings and prepare for generating dimers
     * Seeds the generator, and divides the range of distances in the number
     * of distances or in bins of the bin width. finishOptions calls it; a
     * program without a command line calls it after the setters.
     * \throws gmx::InvalidInputError that names the option of a setting
     * that is out of range
     */
    void finishSettings();

    /*! \brief The distances of the centers of mass, as -mindist, -maxdist and -ndist
     * \param[in] mindist The shortest distance in nm
     * \param[in] maxdist The longest distance in nm
     * \param[in] ndist   The number of distances, less than two for bins of 2.5 pm
     */
    void setDistances(double mindist, double maxdist, int ndist)
    {
        mindist_ = mindist;
        maxdist_ = maxdist;
        ndist_   = ndist;
    }

    /*! \brief Set the algorithm for the orientations, as -randalg
     * \param[in] randAlg The algorithm
     */
    void setRandAlg(RandAlg randAlg) { randAlg_ = randAlg; }

    /*! \brief Read the dimers from a trajectory instead of generating them, as -traj
     * \param[in] trajname The file, none if empty
     */
    void setTrajectory(const std::string &trajname) { trajname_ = trajname; }

    /*! \brief Write the generated dimers to a file, as -ox
     * \param[in] outcoords The file, none if empty
     */
    void setOutputCoordinates(const std::string &outcoords) { outcoords_ = outcoords; }

    /*! \brief Print the statistics of the orientations, with -dbgGD only
     * \param[in] msghandler For the debug output, may be nullptr
     */
    void printStatistics(MsgHandler *msghandler);

    //! Get the original seed
    int seed() const { return dimerseed_; }

    //! \return whether flexibility is used
    bool flexible() const { return flexible_; }

    //! Set a new seed
    void setSeed(int seed);

    /*! Generate the random numbers for all rotations
     * For each orientation there are 6 numbers between 0 and 1, three for
     * the Euler angles of each molecule, see Rotator::randomRotate. With
     * the Gauss algorithm, the orientations are the points of the grid
     * and there is a seventh number, the weight of the orientation.
     * \param[in] ndimers The number of dimers, not used for the Gauss algorithm
     */
    void generateRandomNumbers(int ndimers);

    //! \return The algorithm for sampling orientations
    RandAlg randAlg() const { return randAlg_; }

    /*! \brief Use the Gauss algorithm with a grid of orientations
     * \param[in] nalpha Number of points for alpha of the second molecule
     * \param[in] nbeta  Number of points for beta of each molecule
     * \param[in] ngamma Number of points for gamma of each molecule
     */
    void setGaussGrid(int nalpha, int nbeta, int ngamma)
    {
        randAlg_ = RandAlg::Gauss;
        nalpha_  = nalpha;
        nbeta_   = nbeta;
        ngamma_  = ngamma;
    }

    /*! \brief The number of orientations
     * \param[in] ndimers The number of orientations asked for
     * \return The number of points of the grid for the Gauss algorithm,
     *         ndimers otherwise
     */
    int numberOfOrientations(int ndimers) const;

    //! \return The list of random numbers
    const std::vector<std::vector<double>> &allRandom() const { return allRandom_; }

    //! Empty random number array
    void clearAllRandom() { allRandom_.clear(); }

    /*! \brief Set the random numbers
     * \param[in] newRand The new random numbers
     */
    void addRandomNumbers(const std::vector<double> &newRand)
    {
        allRandom_.push_back(newRand);
    }

    //! Return the number of distancea
    int ndist() const { return ndist_; }

    //! \return the shortest distance between the centers of mass
    double mindist() const { return mindist_; }

    //! Return trajectory name
    const char *trajname() const { return trajname_.empty() ? nullptr : trajname_.c_str(); }

    //! Return max distance
    double maxdist() const { return maxdist_; }

    //! Return bin width
    double binwidth() const { return binWidth_; }

    //! Return output coordinates filename (nullptr if not set)
    const char *outCoords() const { return outcoords_.empty() ? nullptr : outcoords_.c_str(); }

    /*! \brief Write all generated coordinates to file
     * \param[in] actmol  The molecule description (for atom names)
     * \param[in] coords  All the generated coordinate sets to write
     */
    void writeCoords(const ACTMol                              *actmol,
                     const std::vector<std::vector<gmx::RVec>> &coords) const;

    //! Return whether or not there is a trajectory to read
    bool hasTrajectory() const { return !trajname_.empty(); }

    /*! \brief Prepare for generating, this must be called once
     * This will perform monomer minimization if requested.
     * \param[in] msghandler MsgHandler for output and debugging, may be nullptr
     * \param[in] pd         Force field for minimizing input coordinates
     * \param[in] actmol Information about the dimer
     */
    void prepare(MsgHandler       *msghandler,
                 const ForceField *pd,
                 const ACTMol     *actmol);

    /*! \brief Do the actual generation for one pair. Two molecules
     * will be oriented and a distance series will be generated.
     * The memory of the coordinate sets of the previous call is used
     * again, such that it does not have to be allocated for every orientation.
     * \param[in]  msghandler MsgHandler for output and debugging, may be nullptr
     * \param[out] coords     The coordinate sets
     */
    void generateDimers(MsgHandler                *msghandler,
                        std::vector<dimer_coords> *coords);

    /*! \brief Read all the dimers at once from a file. 
     * \param[out]   coords   The coordinates. If empty reading failed or the variable was empty
     */
    void read(std::vector<std::vector<gmx::RVec>> *coords);
};

} // namespace alexandria

#endif // ACT_DIMERGENERATOR_H
