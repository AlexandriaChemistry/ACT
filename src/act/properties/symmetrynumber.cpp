/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
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
 * Implements the computation of the rotational symmetry number.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
// This include has to come first, to prevent the GROMACS definition of
// "real" from messing up the library.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-dtor"
#pragma GCC diagnostic ignored "-Wdeprecated-enum-enum-conversion"
#include <Eigen/Dense>
#include <Eigen/SVD>
#pragma GCC diagnostic pop

#include "symmetrynumber.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

namespace
{

//! Search for distance preserving permutations
class SymmetrySearch
{
private:
    //! Number of atoms
    size_t                            natom_;
    //! Atomic numbers
    const std::vector<int>           &atomicNumbers_;
    //! Coordinates with respect to the center of mass
    std::vector<Eigen::Vector3d>      x_;
    //! Interatomic distances
    std::vector<std::vector<double>>  dist_;
    //! Tolerance
    double                            tolerance_;
    //! Maximum number of steps
    size_t                            maxNodes_;
    //! Number of steps taken
    size_t                            nodes_ = 0;
    //! Order in which atoms are assigned
    std::vector<size_t>               order_;
    //! Candidate images for each atom
    std::vector<std::vector<size_t>>  candidates_;
    //! Current permutation
    std::vector<size_t>               perm_;
    //! Which atoms are already used as image
    std::vector<bool>                 used_;
    //! Number of proper rotations found
    int                               nRotations_ = 0;

    //! \return whether the complete permutation is a proper rotation
    bool isProperRotation() const
    {
        // Kabsch algorithm, restricted to rotations with determinant +1
        Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
        for(size_t i = 0; i < natom_; i++)
        {
            H += x_[i]*x_[perm_[i]].transpose();
        }
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d U = svd.matrixU();
        Eigen::Matrix3d V = svd.matrixV();
        Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
        D(2, 2)           = (V*U.transpose()).determinant() < 0 ? -1 : 1;
        Eigen::Matrix3d R = V*D*U.transpose();
        double msd = 0;
        for(size_t i = 0; i < natom_; i++)
        {
            msd += (R*x_[i] - x_[perm_[i]]).squaredNorm();
        }
        return std::sqrt(msd/natom_) < tolerance_;
    }

    //! Assign an image to the atom at position k in the order
    void search(size_t k)
    {
        if (nodes_ > maxNodes_)
        {
            return;
        }
        if (k == natom_)
        {
            if (isProperRotation())
            {
                nRotations_ += 1;
            }
            return;
        }
        size_t i = order_[k];
        for(auto j : candidates_[i])
        {
            if (used_[j])
            {
                continue;
            }
            nodes_ += 1;
            // All distances to atoms that were assigned before should be preserved
            bool ok = true;
            for(size_t l = 0; l < k && ok; l++)
            {
                size_t m = order_[l];
                ok = std::abs(dist_[i][m] - dist_[j][perm_[m]]) < tolerance_;
            }
            if (ok)
            {
                perm_[i] = j;
                used_[j] = true;
                search(k+1);
                used_[j] = false;
            }
        }
    }

public:
    SymmetrySearch(const std::vector<int>       &atomicNumbers,
                   const std::vector<double>    &masses,
                   const std::vector<gmx::RVec> &coords,
                   double                        tolerance,
                   size_t                        maxNodes)
        : natom_(coords.size()), atomicNumbers_(atomicNumbers),
          tolerance_(tolerance), maxNodes_(maxNodes)
    {
        // Center of mass
        Eigen::Vector3d xcm  = Eigen::Vector3d::Zero();
        double          mtot = 0;
        for(size_t i = 0; i < natom_; i++)
        {
            xcm  += masses[i]*Eigen::Vector3d(coords[i][XX], coords[i][YY], coords[i][ZZ]);
            mtot += masses[i];
        }
        xcm /= mtot;
        for(size_t i = 0; i < natom_; i++)
        {
            x_.push_back(Eigen::Vector3d(coords[i][XX], coords[i][YY], coords[i][ZZ]) - xcm);
        }
        dist_.resize(natom_, std::vector<double>(natom_, 0.0));
        for(size_t i = 0; i < natom_; i++)
        {
            for(size_t j = 0; j < natom_; j++)
            {
                dist_[i][j] = (x_[i] - x_[j]).norm();
            }
        }
        // A rotation preserves the element and the distance to the center of mass
        candidates_.resize(natom_);
        for(size_t i = 0; i < natom_; i++)
        {
            for(size_t j = 0; j < natom_; j++)
            {
                if (atomicNumbers_[i] == atomicNumbers_[j] &&
                    std::abs(x_[i].norm() - x_[j].norm()) < tolerance_)
                {
                    candidates_[i].push_back(j);
                }
            }
        }
        // Assign the atoms with the fewest candidates first, then those
        // furthest from the center since they constrain the rotation most.
        order_.resize(natom_);
        std::iota(order_.begin(), order_.end(), 0);
        std::stable_sort(order_.begin(), order_.end(), [this](size_t a, size_t b)
                         {
                             if (candidates_[a].size() != candidates_[b].size())
                             {
                                 return candidates_[a].size() < candidates_[b].size();
                             }
                             return x_[a].norm() > x_[b].norm();
                         });
        perm_.resize(natom_, 0);
        used_.resize(natom_, false);
    }

    //! \return the symmetry number, or 0 if the search did not finish
    int run()
    {
        search(0);
        return nodes_ > maxNodes_ ? 0 : nRotations_;
    }
};

} // namespace

int rotationalSymmetryNumber(const std::vector<int>       &atomicNumbers,
                             const std::vector<double>    &masses,
                             const std::vector<gmx::RVec> &coords,
                             double                        tolerance,
                             size_t                        maxNodes)
{
    GMX_RELEASE_ASSERT(atomicNumbers.size() == coords.size() &&
                       masses.size() == coords.size(),
                       "Inconsistent input for symmetry number");
    if (coords.size() <= 1)
    {
        return 1;
    }
    SymmetrySearch search(atomicNumbers, masses, coords, tolerance, maxNodes);
    return search.run();
}

} // namespace alexandria
