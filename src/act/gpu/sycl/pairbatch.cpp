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
 * Pair interactions of many conformations on the GPU with SYCL, see
 * computePairBlocks in gpu.h.
 * The results are the same in every run and do not depend on which other
 * conformations are computed, as on the processor, such that for instance
 * relaxing shells gives the same structures. Therefore no atomic
 * operations are used, whose order is unknown, but three kernels:
 * the first computes the force and energy of each interaction in each
 * conformation, the second sums those of each atom in the order of the
 * interactions, the third the energies of the atoms of each conformation.
 * The conformation is the fastest index, such that neighboring work-items
 * read neighboring memory. The molecules of all blocks are computed by
 * the same kernels, each work-item finds its block by bisection.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gpu.h"

#include <algorithm>
#include <utility>
#include <vector>

#include <sycl/sycl.hpp>

#include "sycl/context.h"
#include "sycl/pairkernels.h"

namespace alexandria
{

namespace gpu
{

namespace
{

//! Where the data of a block are in the device arrays, see PairBatchBlock
struct DeviceBlock
{
    //! Start of the coordinates
    uint64_t xOffset;
    //! Start of the interactions
    uint64_t pairOffset;
    //! Start of the interactions times slots
    uint64_t workOffset;
    //! Start of the atoms
    uint64_t atomOffset;
    //! Start of the atoms times slots
    uint64_t atomWorkOffset;
    //! Start of the slots
    uint64_t slotOffset;
    //! Number of atoms
    uint32_t natom;
    //! Distance between the rows of the coordinates
    uint32_t stride;
    //! Number of slots
    uint32_t nslot;
    //! The potential
    PairKernel kernel;
};

//! The device memory for computePairBlocks and the data on the host to copy
struct PairBatchMemory
{
    //! The input of the kernels, see computePairBlocks
    DeviceArray<char>   in;
    //! The output of the kernels, see computePairBlocks
    DeviceArray<char>   out;
    //! The input on the host
    HostArray<char>     inHost;
    //! The output on the host
    HostArray<char>     outHost;
    //! The force of each interaction in each slot, 3 rows of all of these
    DeviceArray<float>  fpair;
    //! The energy of each interaction in each slot
    DeviceArray<float>  epair;
    //! The energy of the interactions of each atom in each slot, as the first atom
    DeviceArray<double> eatom;
};

//! The device memory, never deleted, see gpuQueue
PairBatchMemory &pairBatchMemory()
{
    static PairBatchMemory *memory = new PairBatchMemory;
    return *memory;
}

/*! \brief The block of a work-item, by bisection of the starts of the blocks
 * \param[in] blocks The blocks
 * \param[in] nblock Number of blocks
 * \param[in] index  The index of the work-item
 * \param[in] start  Member of DeviceBlock with the start of the work of each block
 * \return the block whose work holds index
 */
inline uint32_t findBlock(const DeviceBlock *blocks, uint32_t nblock, uint64_t index,
                          uint64_t DeviceBlock::*start)
{
    uint32_t lo = 0, hi = nblock;
    while (hi - lo > 1)
    {
        uint32_t mid = (lo + hi)/2;
        if (blocks[mid].*start <= index)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }
    return lo;
}

/*! \brief The kernels for computePairBlocks, see the description of this file
 * \tparam Forces Whether to compute forces
 * \param[in]  queue     The queue
 * \param[in]  blocks    The blocks
 * \param[in]  nblock    Number of blocks
 * \param[in]  nwork     Interactions times slots of all blocks
 * \param[in]  natomWork Atoms times slots of all blocks
 * \param[in]  nslotAll  Slots of all blocks
 * \param[in]  x         The coordinates
 * \param[in]  pairs     The interactions
 * \param[in]  atomStart Start of the interactions of each atom in atomPairs
 * \param[in]  atomPairs Twice the interactions of each atom, plus one if it is the second
 * \param[out] fpair     The force of each interaction, 3*nwork
 * \param[out] epair     The energy terms of each interaction, pairEnergyTerms*nwork
 * \param[out] eatom     The energy terms of each atom, pairEnergyTerms*natomWork
 * \param[out] f         The forces, 3*natomWork
 * \param[out] e         The energy terms, pairEnergyTerms*nslotAll
 */
template <bool Forces>
void pairKernels(sycl::queue       &queue,
                    const DeviceBlock *blocks,
                    uint32_t           nblock,
                    size_t             nwork,
                    size_t             natomWork,
                    size_t             nslotAll,
                    const int64_t     *x,
                    const PairEntry *pairs,
                    const uint64_t    *atomStart,
                    const uint32_t    *atomPairs,
                    float             *fpair,
                    float             *epair,
                    double            *eatom,
                    float             *f,
                    double            *e)
{
    // Each interaction in each slot
    queue.parallel_for(sycl::range<1>(nwork), [=](sycl::id<1> id)
    {
        uint64_t           w     = id[0];
        const DeviceBlock &b     = blocks[findBlock(blocks, nblock, w, &DeviceBlock::workOffset)];
        uint64_t           local = w - b.workOffset;
        uint64_t           p     = local / b.nslot;
        uint64_t           c     = local % b.nslot;
        const PairEntry &pair  = pairs[b.pairOffset + p];
        uint64_t           xi    = b.xOffset + 3*uint64_t(pair.i)*b.stride + c;
        uint64_t           xj    = b.xOffset + 3*uint64_t(pair.j)*b.stride + c;
        float              dx[3];
        float              rsq   = 0;
        for (int m = 0; m < 3; m++)
        {
            // The exact difference of the integers, rounded once
            dx[m] = static_cast<float>(x[xi + m*b.stride] - x[xj + m*b.stride])*
                static_cast<float>(coordinateUnitNm);
            rsq  += dx[m]*dx[m];
        }
        float e[pairEnergyTerms];
        float fscal = pairInteraction(b.kernel, pair.p, rsq, e);
        for (int k = 0; k < pairEnergyTerms; k++)
        {
            epair[k*nwork + w] = e[k];
        }
        if constexpr (Forces)
        {
            for (int m = 0; m < 3; m++)
            {
                fpair[m*nwork + w] = rsq > 0 ? fscal*dx[m] : 0.0f;
            }
        }
    });
    // The sums for each atom in the order of the interactions, the
    // energies of the interactions in which it is the first atom
    queue.parallel_for(sycl::range<1>(natomWork), [=](sycl::id<1> id)
    {
        uint64_t           w     = id[0];
        const DeviceBlock &b     = blocks[findBlock(blocks, nblock, w, &DeviceBlock::atomWorkOffset)];
        uint64_t           local = w - b.atomWorkOffset;
        uint64_t           a     = local / b.nslot;
        uint64_t           c     = local % b.nslot;
        double             ea[pairEnergyTerms] = { 0, 0, 0 };
        float              fa[3] = { 0, 0, 0 };
        for (uint64_t k = atomStart[b.atomOffset + a]; k < atomStart[b.atomOffset + a + 1]; k++)
        {
            uint64_t p      = atomPairs[k] / 2;
            bool     second = (atomPairs[k] % 2) == 1;
            uint64_t wp     = b.workOffset + p*b.nslot + c;
            if (!second)
            {
                for (int t = 0; t < pairEnergyTerms; t++)
                {
                    ea[t] += epair[t*nwork + wp];
                }
            }
            if constexpr (Forces)
            {
                for (int m = 0; m < 3; m++)
                {
                    float fp = fpair[m*nwork + wp];
                    fa[m]   += second ? -fp : fp;
                }
            }
        }
        for (int t = 0; t < pairEnergyTerms; t++)
        {
            eatom[t*natomWork + w] = ea[t];
        }
        if constexpr (Forces)
        {
            for (int m = 0; m < 3; m++)
            {
                f[3*b.atomWorkOffset + (3*a + m)*b.nslot + c] = fa[m];
            }
        }
    });
    // The energy of each slot in the order of the atoms
    queue.parallel_for(sycl::range<1>(nslotAll), [=](sycl::id<1> id)
    {
        uint64_t           s  = id[0];
        const DeviceBlock &b  = blocks[findBlock(blocks, nblock, s, &DeviceBlock::slotOffset)];
        uint64_t           c  = s - b.slotOffset;
        for (int t = 0; t < pairEnergyTerms; t++)
        {
            double ec = 0;
            for (uint64_t a = 0; a < b.natom; a++)
            {
                ec += eatom[t*natomWork + b.atomWorkOffset + a*b.nslot + c];
            }
            e[t*nslotAll + s] = ec;
        }
    });
}

} // namespace

bool computePairBlocks(const PairBatchBlock *blocks,
                       size_t                nblock,
                       bool                  forces)
{
    sycl::queue *queue = gpuQueue();
    if (!queue)
    {
        return false;
    }
    auto &mem = pairBatchMemory();
    // The blocks with work, the coordinates of each array once
    std::vector<size_t>      used;
    std::vector<DeviceBlock> dblocks;
    std::vector<char>        firstX;
    std::vector<std::pair<const int64_t *, uint64_t>> xKnown;
    uint64_t nx = 0, npair = 0, natom = 0, nwork = 0, natomWork = 0, nslotAll = 0;
    for (size_t k = 0; k < nblock; k++)
    {
        const auto &cb = blocks[k];
        if (cb.npair == 0 || cb.nslot == 0)
        {
            continue;
        }
        used.push_back(k);
        DeviceBlock b;
        b.natom  = static_cast<uint32_t>(cb.natom);
        b.stride = static_cast<uint32_t>(cb.stride);
        b.nslot  = static_cast<uint32_t>(cb.nslot);
        b.kernel = cb.kernel;
        auto xk  = std::find_if(xKnown.begin(), xKnown.end(),
                                [&cb](const auto &xo) { return xo.first == cb.x; });
        firstX.push_back(xk == xKnown.end());
        if (xk != xKnown.end())
        {
            b.xOffset = xk->second;
        }
        else
        {
            b.xOffset = nx;
            xKnown.push_back({ cb.x, nx });
            nx += 3*cb.natom*cb.stride;
        }
        b.pairOffset     = npair;
        b.workOffset     = nwork;
        b.atomOffset     = natom;
        b.atomWorkOffset = natomWork;
        b.slotOffset     = nslotAll;
        npair           += cb.npair;
        natom           += cb.natom;
        nwork           += cb.npair*cb.nslot;
        natomWork       += cb.natom*cb.nslot;
        nslotAll        += cb.nslot;
        dblocks.push_back(b);
    }
    if (used.empty())
    {
        return true;
    }
    // All input in one buffer and all output in another, such that each
    // is copied at once
    PackedArrays in;
    size_t       oBlocks    = in.add<DeviceBlock>(dblocks.size());
    size_t       oX         = in.add<int64_t>(nx);
    size_t       oPairs     = in.add<PairEntry>(npair);
    size_t       oAtomStart = in.add<uint64_t>(natom + 1);
    size_t       oAtomPairs = in.add<uint32_t>(2*npair);
    PackedArrays out;
    size_t       oF         = out.add<float>(forces ? 3*natomWork : 0);
    size_t       oE         = out.add<double>(pairEnergyTerms*nslotAll);
    char        *hin        = mem.inHost.reserve(*queue, in.size());
    char        *hout       = mem.outHost.reserve(*queue, out.size());
    auto        *hBlocks    = reinterpret_cast<DeviceBlock *>(hin + oBlocks);
    auto        *hX         = reinterpret_cast<int64_t *>(hin + oX);
    auto        *hPairs     = reinterpret_cast<PairEntry *>(hin + oPairs);
    auto        *hAtomStart = reinterpret_cast<uint64_t *>(hin + oAtomStart);
    auto        *hAtomPairs = reinterpret_cast<uint32_t *>(hin + oAtomPairs);
    std::copy(dblocks.begin(), dblocks.end(), hBlocks);
    // The interactions of each atom in the order of the interactions
    hAtomStart[0] = 0;
    for (size_t u = 0; u < used.size(); u++)
    {
        const auto &cb = blocks[used[u]];
        const auto &b  = dblocks[u];
        if (firstX[u])
        {
            std::copy(cb.x, cb.x + 3*cb.natom*cb.stride, hX + b.xOffset);
        }
        std::copy(cb.pairs, cb.pairs + cb.npair, hPairs + b.pairOffset);
        uint64_t *start = hAtomStart + b.atomOffset;
        for (size_t a = 1; a <= cb.natom; a++)
        {
            start[a] = 0;
        }
        for (size_t p = 0; p < cb.npair; p++)
        {
            start[cb.pairs[p].i + 1] += 1;
            start[cb.pairs[p].j + 1] += 1;
        }
        // start[0] is the end of the previous block
        for (size_t a = 1; a <= cb.natom; a++)
        {
            start[a] += start[a-1];
        }
        std::vector<uint64_t> next(start, start + cb.natom);
        for (size_t p = 0; p < cb.npair; p++)
        {
            hAtomPairs[next[cb.pairs[p].i]++] = 2*p;
            hAtomPairs[next[cb.pairs[p].j]++] = 2*p + 1;
        }
    }
    char    *din  = mem.in.reserve(*queue, in.size());
    char    *dout = mem.out.reserve(*queue, out.size());
    auto     dep  = mem.epair.reserve(*queue, pairEnergyTerms*nwork);
    auto     dea  = mem.eatom.reserve(*queue, pairEnergyTerms*natomWork);
    auto     dfp  = forces ? mem.fpair.reserve(*queue, 3*nwork) : nullptr;
    uint32_t nb   = static_cast<uint32_t>(dblocks.size());
    auto     db   = reinterpret_cast<const DeviceBlock *>(din + oBlocks);
    auto     dx   = reinterpret_cast<const int64_t *>(din + oX);
    auto     dp   = reinterpret_cast<const PairEntry *>(din + oPairs);
    auto     das  = reinterpret_cast<const uint64_t *>(din + oAtomStart);
    auto     dap  = reinterpret_cast<const uint32_t *>(din + oAtomPairs);
    auto     df   = forces ? reinterpret_cast<float *>(dout + oF) : nullptr;
    auto     de   = reinterpret_cast<double *>(dout + oE);
    queue->memcpy(din, hin, in.size());
    if (forces)
    {
        pairKernels<true>(*queue, db, nb, nwork, natomWork, nslotAll, dx, dp, das, dap, dfp, dep, dea, df, de);
    }
    else
    {
        pairKernels<false>(*queue, db, nb, nwork, natomWork, nslotAll, dx, dp, das, dap, dfp, dep, dea, df, de);
    }
    queue->memcpy(hout, dout, out.size());
    queue->wait_and_throw();
    const auto *hE = reinterpret_cast<const double *>(hout + oE);
    const auto *hF = reinterpret_cast<const float *>(hout + oF);
    for (size_t u = 0; u < used.size(); u++)
    {
        const auto &cb = blocks[used[u]];
        const auto &b  = dblocks[u];
        for (int t = 0; t < pairEnergyTerms; t++)
        {
            for (size_t c = 0; c < cb.nslot; c++)
            {
                cb.energies[t*cb.stride + c] += hE[t*nslotAll + b.slotOffset + c];
            }
        }
        if (forces)
        {
            const float *fb = hF + 3*b.atomWorkOffset;
            for (size_t i = 0; i < 3*cb.natom; i++)
            {
                for (size_t c = 0; c < cb.nslot; c++)
                {
                    cb.f[i*cb.stride + c] += fb[i*cb.nslot + c];
                }
            }
        }
    }
    return true;
}

} // namespace gpu

} // namespace alexandria
