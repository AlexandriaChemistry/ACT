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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <type_traits>

#include "act/forcefield/forcefield_parametername.h"
#include "act/gpu/pairrows.h"
#include "act/utility/timings.h"
#include "gromacs/math/functions.h"
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

bool gpuPairKernel(Potential mypot, gpu::PairKernel *kernel)
{
    switch (mypot)
    {
    case Potential::COULOMB_POINT:
    case Potential::COULOMB_GAUSSIAN: *kernel = gpu::PairKernel::Coulomb; return true;
    case Potential::POLARIZATION:     *kernel = gpu::PairKernel::Polarization; return true;
    case Potential::LJ12_6:           *kernel = gpu::PairKernel::LJ12_6; return true;
    case Potential::LJ14_7:           *kernel = gpu::PairKernel::LJ14_7; return true;
    case Potential::WANG_BUCKINGHAM:  *kernel = gpu::PairKernel::WangBuckingham; return true;
    case Potential::MACDANIEL_SCHMIDT: *kernel = gpu::PairKernel::DoubleExponential; return true;
    case Potential::HARMONIC_BONDS:   *kernel = gpu::PairKernel::HarmonicBond; return true;
    case Potential::CUBIC_BONDS:      *kernel = gpu::PairKernel::CubicBond; return true;
    case Potential::MORSE_BONDS:      *kernel = gpu::PairKernel::MorseBond; return true;
    case Potential::HUA_BONDS:        *kernel = gpu::PairKernel::HuaBond; return true;
    default:                          return false;
    }
}

/*! \brief The interaction of a pair for the GPU, see gpuPairEntry
 * \tparam E TopologyEntry or PairView
 */
template <typename E>
static bool gpuPairEntryT(gpu::PairKernel             kernel,
                          const E                    &b,
                          const std::vector<ActAtom> &atoms,
                          double                      epsilonr,
                          gpu::PairEntry             *pe)
{
    const auto &ind    = b.atomIndices();
    pe->i = static_cast<uint32_t>(ind[0]);
    pe->j = static_cast<uint32_t>(ind[1]);
    std::fill(pe->p, pe->p + 8, 0.0f);
    switch (kernel)
    {
    case gpu::PairKernel::Coulomb:
    {
        const auto  &atomi = atoms[ind[0]];
        const auto  &atomj = atoms[ind[1]];
        const double qi[2] = { atomi.charge(), atomi.chargeVS() };
        const double qj[2] = { atomj.charge(), atomj.chargeVS() };
        const double zi[2] = { atomi.zeta(), atomi.zetaVS() };
        const double zj[2] = { atomj.zeta(), atomj.zetaVS() };
        for (int a = 0; a < 2; a++)
        {
            for (int c = 0; c < 2; c++)
            {
                int k       = 2*a + c;
                pe->p[k]    = static_cast<float>(ONE_4PI_EPS0*qi[a]*qj[c]/epsilonr);
                pe->p[4+k]  = static_cast<float>(combinedWidth(zi[a], zj[c]));
            }
        }
        return true;
    }
    default:
    {
        // The other kernels only need the parameters, see pairrows.h
        double qshell = gpu::PairKernel::Polarization == kernel ? atoms[ind[1]].charge() : 0.0;
        return gpu::rows::pairRow(kernel, b.params().data(), qshell, ONE_4PI_EPS0, pe->p);
    }
    }
}

bool gpuPairEntry(gpu::PairKernel             kernel,
                  const TopologyEntry        &b,
                  const std::vector<ActAtom> &atoms,
                  double                      epsilonr,
                  gpu::PairEntry             *pe)
{
    return gpuPairEntryT(kernel, b, atoms, epsilonr, pe);
}

bool gpuPairEntry(gpu::PairKernel             kernel,
                  const PairView             &b,
                  const std::vector<ActAtom> &atoms,
                  double                      epsilonr,
                  gpu::PairEntry             *pe)
{
    return gpuPairEntryT(kernel, b, atoms, epsilonr, pe);
}

bool gpuBondedKernel(Potential mypot, gpu::BondedKernel *kernel)
{
    switch (mypot)
    {
    case Potential::HARMONIC_ANGLES:    *kernel = gpu::BondedKernel::HarmonicAngle; return true;
    case Potential::LINEAR_ANGLES:      *kernel = gpu::BondedKernel::LinearAngle; return true;
    case Potential::FOURIER_DIHEDRALS:  *kernel = gpu::BondedKernel::FourierDihedral; return true;
    case Potential::HARMONIC_DIHEDRALS: *kernel = gpu::BondedKernel::HarmonicDihedral; return true;
    default:                            return false;
    }
}

void gpuBondedEntry(gpu::BondedKernel kernel, const TopologyEntry &b, gpu::BondedEntry *be)
{
    const auto &ind    = b.atomIndices();
    const auto &params = b.params();
    for (int k = 0; k < 4; k++)
    {
        be->a[k] = static_cast<uint32_t>(k < static_cast<int>(ind.size()) ? ind[k] : ind[0]);
    }
    gpu::rows::bondedRow(kernel, params.data(), params.size(), be->p);
}

bool gpuVsiteParameters(InteractionType itype, const TopologyEntry &b, int index[3], double factor[3])
{
    for (int k = 0; k < 3; k++)
    {
        index[k]  = -1;
        factor[k] = 1;
    }
    const auto use = [index, factor](int k, int i, double f = 1)
    {
        index[k]  = i;
        factor[k] = f;
    };
    // The sign of the out-of-plane distance of the site, see Vsite3OUT
    const auto sign = [&b]() { return static_cast<const Vsite3OUT *>(b.self())->sign(); };
    switch (itype)
    {
    case InteractionType::VSITE1:
        return true;
    case InteractionType::VSITE2:
        use(0, vsite2A);
        return true;
    case InteractionType::VSITE2FD:
        use(0, vsite2fdA);
        return true;
    case InteractionType::VSITE3:
        use(0, vsite3A);
        use(1, vsite3B);
        return true;
    case InteractionType::VSITE3S:
        use(0, vsite3sA);
        use(1, vsite3sA);
        return true;
    case InteractionType::VSITE3FD:
        use(0, vsite3fdA);
        use(1, vsite3fdB);
        return true;
    case InteractionType::VSITE3OUT:
        use(0, vsite3outA);
        use(1, vsite3outB);
        use(2, vsite3outC, sign());
        return true;
    case InteractionType::VSITE3OUTS:
        use(0, vsite3outsA);
        use(1, vsite3outsA);
        use(2, vsite3outsC, sign());
        return true;
    case InteractionType::VSITE4:
        use(0, vsite4A);
        use(1, vsite4B);
        use(2, vsite4C);
        return true;
    case InteractionType::VSITE4S:
        use(0, vsite4sA);
        use(1, vsite4sB);
        use(2, vsite4sB);
        return true;
    case InteractionType::VSITE4S3:
        use(0, vsite4s3A);
        use(1, vsite4s3A);
        use(2, vsite4s3A);
        return true;
    default:
        return false;
    }
}

bool gpuVsiteEntry(InteractionType itype, const TopologyEntry &b, gpu::VsiteEntry *ve)
{
    const auto &ai     = b.atomIndices();
    const auto &params = b.params();
    if (ai.empty() || params.empty())
    {
        return false;
    }
    std::fill(ve->p, ve->p + 3, 0.0f);
    std::fill(ve->parent, ve->parent + 4, 0u);
    ve->site    = static_cast<uint32_t>(ai.back());
    ve->nparent = static_cast<uint32_t>(ai.size() - 1);
    for (size_t k = 0; k + 1 < ai.size(); k++)
    {
        ve->parent[k] = static_cast<uint32_t>(ai[k]);
    }
    // The second parent first, for the routines that get it first
    const auto secondFirst = [ve]() { std::swap(ve->parent[0], ve->parent[1]); };
    switch (itype)
    {
    case InteractionType::VSITE1:     ve->kind = gpu::VsiteKind::Copy; break;
    case InteractionType::VSITE2:     ve->kind = gpu::VsiteKind::Linear2; break;
    case InteractionType::VSITE2FD:   ve->kind = gpu::VsiteKind::Fd2; break;
    case InteractionType::VSITE3:
    case InteractionType::VSITE3S:    ve->kind = gpu::VsiteKind::Linear3; break;
    case InteractionType::VSITE3FD:   ve->kind = gpu::VsiteKind::Fd3; break;
    case InteractionType::VSITE3OUT:
    case InteractionType::VSITE3OUTS: ve->kind = gpu::VsiteKind::Out3; secondFirst(); break;
    case InteractionType::VSITE4:
    case InteractionType::VSITE4S:
    case InteractionType::VSITE4S3:   ve->kind = gpu::VsiteKind::Linear4; break;
    default:                          return false;
    }
    if (InteractionType::VSITE4S3 == itype)
    {
        secondFirst();
    }
    // The parameters, as the GPU makes them from the recipes, see gpuVsiteParameters
    int    index[3];
    double factor[3];
    (void) gpuVsiteParameters(itype, b, index, factor);
    for (int k = 0; k < 3; k++)
    {
        if (index[k] >= 0)
        {
            ve->p[k] = static_cast<float>(factor[k]*params[index[k]]);
        }
    }
    return true;
}

/*! \brief The interactions with prepared parameters for the GPU, see gpuPairEntry
 * \param[in]  kernel   The kernel
 * \param[in]  entries  The interactions
 * \param[in]  atoms    The atoms
 * \param[in]  epsilonr Relative dielectric constant
 * \param[out] pairs    The interactions for the GPU, without those that have none
 */
template <typename Range>
static void preparePairs(gpu::PairKernel              kernel,
                         const Range                 &entries,
                         const std::vector<ActAtom>  &atoms,
                         double                       epsilonr,
                         std::vector<gpu::PairEntry> *pairs)
{
    pairs->clear();
    pairs->reserve(entries.size());
    for (const auto &b : entries)
    {
        gpu::PairEntry e;
        if (gpuPairEntry(kernel, *b.operator->(), atoms, epsilonr, &e))
        {
            pairs->push_back(e);
        }
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
    // The GPU has kernels for more potentials, see gpuPairKernel, but each
    // potential of each molecule returns its own forces, which are copied
    // and added on the host. Except for Coulomb, which has the most
    // interactions, that takes longer than the processor needs to compute
    // the potentials, as measured for training with 94 compounds.
    return (gpu::useGpu() &&
            (Potential::COULOMB_POINT == mypot || Potential::COULOMB_GAUSSIAN == mypot));
}

void GpuBatch::clear()
{
    ntask_    = 0;
    computed_ = false;
}

template <typename Range>
void GpuBatch::addT(Potential                   mypot,
                    const Range                &entries,
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
    (void) gpuPairKernel(mypot, &t.kernel);
    preparePairs(t.kernel, entries, atoms, epsilonr, &t.pairs);
    computed_ = false;
}

void GpuBatch::add(Potential mypot, TopologyEntryRange entries, const std::vector<ActAtom> &atoms,
                   double epsilonr, size_t stride, size_t nslot, const Coordinate *x)
{
    addT(mypot, entries, atoms, epsilonr, stride, nslot, x);
}

void GpuBatch::add(Potential mypot, const PairRange &entries, const std::vector<ActAtom> &atoms,
                   double epsilonr, size_t stride, size_t nslot, const Coordinate *x)
{
    addT(mypot, entries, atoms, epsilonr, stride, nslot, x);
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
    const ScopedTimer   scoped(timer);
    std::vector<gpu::PairBatchBlock> blocks(ntask_);
    for (size_t k = 0; k < ntask_; k++)
    {
        auto &t = tasks_[k];
        t.energies.assign(gpu::pairEnergyTerms*t.stride, 0.0);
        if (forces)
        {
            t.f.assign(3*t.natom*t.stride, 0.0f);
        }
        blocks[k] = { t.natom, t.stride, t.nslot, t.x, t.kernel, t.pairs.data(), t.pairs.size(),
                      forces ? t.f.data() : nullptr, t.energies.data() };
    }
    computed_ = gpu::computePairBlocks(blocks.data(), blocks.size(), forces);
    return computed_;
}

void GpuBatch::addResults(size_t k, float *f, double *energies) const
{
    GMX_RELEASE_ASSERT(computed_ && k < ntask_, "No results for this task");
    const auto &t = tasks_[k];
    for (int term = 0; term < gpu::pairEnergyTerms; term++)
    {
        for (size_t c = 0; c < t.nslot; c++)
        {
            energies[term*t.stride + c] += t.energies[term*t.stride + c];
        }
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

/*! \brief Compute a potential for all conformations on the GPU, see computeGpu
 * \tparam Range TopologyEntryRange or PairRange
 */
template <typename Range>
static bool computeGpuT(Potential                   mypot,
                const Range                &entries,
                const std::vector<ActAtom> &atoms,
                double                      epsilonr,
                size_t                      stride,
                size_t                      nslot,
                const Coordinate           *x,
                float                      *f,
                double                     *energies,
                bool                        forces)
{
    gpu::PairKernel kernel;
    if (!gpu::useGpu() || entries.size()*nslot < gpuMinimumWork() || !gpuPairKernel(mypot, &kernel))
    {
        return false;
    }
    // Kept between calls, such that the memory is used again
    static thread_local std::vector<gpu::PairEntry> pairs;
    preparePairs(kernel, entries, atoms, epsilonr, &pairs);
    const gpu::PairBatchBlock block = { atoms.size(), stride, nslot, x, kernel, pairs.data(), pairs.size(),
                                  forces ? f : nullptr, energies };
    return gpu::computePairBlocks(&block, 1, forces);
}

bool computeGpu(Potential mypot, TopologyEntryRange entries, const std::vector<ActAtom> &atoms,
                double epsilonr, size_t stride, size_t nslot, const Coordinate *x, float *f,
                double *energies, bool forces)
{
    return computeGpuT(mypot, entries, atoms, epsilonr, stride, nslot, x, f, energies, forces);
}

bool computeGpu(Potential mypot, const PairRange &entries, const std::vector<ActAtom> &atoms,
                double epsilonr, size_t stride, size_t nslot, const Coordinate *x, float *f,
                double *energies, bool forces)
{
    return computeGpuT(mypot, entries, atoms, epsilonr, stride, nslot, x, f, energies, forces);
}

} // namespace alexandria
