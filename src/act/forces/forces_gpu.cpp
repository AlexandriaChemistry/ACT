/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
 * Implements the potentials on a GPU of forces_gpu.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "forces_gpu.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <type_traits>

#include "act/utility/timings.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

static_assert(std::is_same_v<Coordinate, int64_t>, "The GPU code expects 64-bit integer coordinates");

//! Default of gpuMinimumWork
static const size_t defaultGpuMinimumWork = 100000;

/*! \brief Width of the combined Gaussian, as in Coulomb_GG, zero for point charges
 * \param[in] zi Width of the first charge, zero for a point charge
 * \param[in] zj Width of the second charge, zero for a point charge
 * \return The combined width
 */
static double combinedWidth(double zi, double zj)
{
    if (zi == 0 || zj == 0)
    {
        return zi + zj;
    }
    return zi*zj/std::sqrt(zi*zi + zj*zj);
}

/*! \brief The Coulomb interactions for the GPU, as Coulomb::prepare in forces_explicit_simd.cpp
 * \param[in]  entries  The interactions
 * \param[in]  atoms    The atoms
 * \param[in]  epsilonr Relative dielectric constant
 * \param[out] pairs    The interactions for the GPU
 */
static void prepareCoulomb(TopologyEntryRange              entries,
                           const std::vector<ActAtom>     &atoms,
                           double                          epsilonr,
                           std::vector<gpu::CoulombPair>  *pairs)
{
    pairs->clear();
    pairs->reserve(entries.size());
    for (const auto &b : entries)
    {
        const auto        &ind   = b->atomIndices();
        const auto        &atomi = atoms[ind[0]];
        const auto        &atomj = atoms[ind[1]];
        const double       qi[2] = { atomi.charge(), atomi.chargeVS() };
        const double       qj[2] = { atomj.charge(), atomj.chargeVS() };
        const double       zi[2] = { atomi.zeta(), atomi.zetaVS() };
        const double       zj[2] = { atomj.zeta(), atomj.zetaVS() };
        gpu::CoulombPair   pair;
        pair.i = static_cast<uint32_t>(ind[0]);
        pair.j = static_cast<uint32_t>(ind[1]);
        for (int a = 0; a < 2; a++)
        {
            for (int c = 0; c < 2; c++)
            {
                int k        = 2*a + c;
                pair.qq[k]   = static_cast<float>(ONE_4PI_EPS0*qi[a]*qj[c]/epsilonr);
                pair.zeta[k] = static_cast<float>(combinedWidth(zi[a], zj[c]));
            }
        }
        pairs->push_back(pair);
    }
}

//! The least work for the GPU, see gpuMinimumWork
static size_t &minimumWork()
{
    static size_t work = []()
    {
        const char *env = std::getenv("ACT_GPU_MIN_WORK");
        return env ? static_cast<size_t>(std::strtoull(env, nullptr, 10)) : defaultGpuMinimumWork;
    }();
    return work;
}

size_t gpuMinimumWork()
{
    return minimumWork();
}

size_t setGpuMinimumWork(size_t work)
{
    size_t before = minimumWork();
    minimumWork() = work;
    return before;
}

bool GpuBatch::canCompute(Potential mypot)
{
    return (gpu::useGpu() &&
            (Potential::COULOMB_POINT == mypot || Potential::COULOMB_GAUSSIAN == mypot));
}

void GpuBatch::clear()
{
    ntask_    = 0;
    computed_ = false;
}

void GpuBatch::add(Potential                   mypot,
                   TopologyEntryRange          entries,
                   const std::vector<ActAtom> &atoms,
                   double                      epsilonr,
                   size_t                      stride,
                   size_t                      nslot,
                   const Coordinate           *x)
{
    GMX_RELEASE_ASSERT(canCompute(mypot), "The GPU cannot compute this potential");
    if (ntask_ == tasks_.size())
    {
        tasks_.emplace_back();
    }
    auto &t  = tasks_[ntask_++];
    t.natom  = atoms.size();
    t.stride = stride;
    t.nslot  = nslot;
    t.x      = x;
    prepareCoulomb(entries, atoms, epsilonr, &t.pairs);
    computed_ = false;
}

bool GpuBatch::compute(bool forces)
{
    computed_ = false;
    forces_   = forces;
    size_t work = 0;
    for (size_t k = 0; k < ntask_; k++)
    {
        work += tasks_[k].pairs.size()*tasks_[k].nslot;
    }
    if (ntask_ == 0 || work < gpuMinimumWork())
    {
        return false;
    }
    static const size_t timer = Timings::instance().index("GPU pair potentials");
    ScopedTimer         scoped(timer);
    std::vector<gpu::CoulombBlock> blocks(ntask_);
    for (size_t k = 0; k < ntask_; k++)
    {
        auto &t = tasks_[k];
        t.energies.assign(t.nslot, 0.0);
        if (forces)
        {
            t.f.assign(3*t.natom*t.stride, 0.0f);
        }
        blocks[k] = { t.natom, t.stride, t.nslot, t.x, t.pairs.data(), t.pairs.size(),
                      forces ? t.f.data() : nullptr, t.energies.data() };
    }
    computed_ = gpu::computeCoulombBlocks(blocks.data(), blocks.size(), forces);
    return computed_;
}

void GpuBatch::addResults(size_t k, float *f, double *energies) const
{
    GMX_RELEASE_ASSERT(computed_ && k < ntask_, "No results for this task");
    const auto &t = tasks_[k];
    for (size_t c = 0; c < t.nslot; c++)
    {
        energies[c] += t.energies[c];
    }
    if (forces_)
    {
        for (size_t i = 0; i < 3*t.natom; i++)
        {
            for (size_t c = 0; c < t.nslot; c++)
            {
                f[i*t.stride + c] += t.f[i*t.stride + c];
            }
        }
    }
}

bool computeGpu(Potential                   mypot,
                TopologyEntryRange          entries,
                const std::vector<ActAtom> &atoms,
                double                      epsilonr,
                size_t                      stride,
                size_t                      nslot,
                const Coordinate           *x,
                float                      *f,
                double                     *energies,
                bool                        forces)
{
    if (!gpu::useGpu() || entries.size()*nslot < gpuMinimumWork() ||
        !(Potential::COULOMB_POINT == mypot || Potential::COULOMB_GAUSSIAN == mypot))
    {
        return false;
    }
    // Kept between calls, such that the memory is used again
    static thread_local std::vector<gpu::CoulombPair> pairs;
    prepareCoulomb(entries, atoms, epsilonr, &pairs);
    // The Coulomb energy is the first energy term
    return gpu::computeCoulomb(atoms.size(), stride, nslot, x, pairs.data(), pairs.size(),
                               f, energies, forces);
}

} // namespace alexandria
