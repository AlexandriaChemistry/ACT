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
 * Coulomb interactions of many conformations on the GPU with SYCL.
 * The results are the same in every run and do not depend on which other
 * conformations are computed, as on the processor, such that for instance
 * relaxing shells gives the same structures. Therefore no atomic
 * operations are used, whose order is unknown, but three kernels:
 * the first computes the force and energy of each interaction in each
 * conformation, the second sums those of each atom in the order of the
 * interactions, the third the energies of the atoms of each conformation.
 * The conformation is the fastest index, such that neighboring work-items
 * read neighboring memory.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gpu.h"

#include <algorithm>
#include <vector>

#include <sycl/sycl.hpp>

#include "sycl/context.h"

namespace alexandria
{

namespace gpu
{

namespace
{

//! The unit of the coordinates in nm, 2^-56, see coordinate.h
constexpr double coordinateUnitNm = 1.0/72057594037927936.0;

//! 1/sqrt(pi)
constexpr float invSqrtPi = 0.564189583547756287f;

/*! \brief Scalar force divided by r of two Gaussian charges, without charges
 * As gaussianForceOverR in forces_simd.cpp in single precision: below
 * x^2 = 0.16 the series is more accurate than the difference of the terms.
 * \param[in] zeta Width of the combined distribution
 * \param[in] rinv Inverse distance
 * \param[in] x    zeta times the distance
 * \param[in] erfx The error function of x
 * \return the scalar force divided by r
 */
inline float gaussianForceOverR(float zeta, float rinv, float x, float erfx)
{
    float x2 = x*x;
    if (x2 < 0.16f)
    {
        return zeta*zeta*zeta*2*invSqrtPi*
            (2.0f/3.0f + x2*(-2.0f/5.0f + x2*(1.0f/7.0f + x2*(-1.0f/27.0f + x2*(1.0f/132.0f)))));
    }
    return (erfx*rinv - 2*invSqrtPi*zeta*sycl::exp(-x2))*rinv*rinv;
}

/*! \brief Energy and scalar force divided by r of one term of a Coulomb interaction
 * As Coulomb::compute in forces_explicit_simd.cpp: particles on top of each
 * other only have an energy with Gaussian charges, and no force.
 * \param[in]  qq    The charge product
 * \param[in]  zeta  The combined width, zero for point charges
 * \param[in]  rsq   The distance squared
 * \param[in]  rinv  The inverse distance, not used if rsq is zero
 * \param[out] fscal The scalar force divided by r, added to
 * \return the energy
 */
inline float coulombTerm(float qq, float zeta, float rsq, float rinv, float *fscal)
{
    if (zeta == 0)
    {
        if (rsq > 0)
        {
            *fscal += qq*rinv*rinv*rinv;
            return qq*rinv;
        }
        return 0;
    }
    if (rsq > 0)
    {
        float rz   = zeta*rsq*rinv;
        float erfz = sycl::erf(rz);
        *fscal    += qq*gaussianForceOverR(zeta, rinv, rz, erfz);
        return qq*erfz*rinv;
    }
    return qq*2*zeta*invSqrtPi;
}

//! The device memory for computeCoulomb
struct CoulombMemory
{
    //! The coordinates
    DeviceArray<int64_t>     x;
    //! The forces of the atoms
    DeviceArray<float>       f;
    //! The energy of each conformation
    DeviceArray<double>      e;
    //! The interactions
    DeviceArray<CoulombPair> pairs;
    //! The force of each interaction in each conformation
    DeviceArray<float>       fpair;
    //! The energy of each interaction in each conformation
    DeviceArray<float>       epair;
    //! The energy of the interactions of each atom as the first atom
    DeviceArray<double>      eatom;
    //! Start of the interactions of each atom in atomPairs
    DeviceArray<uint32_t>    atomStart;
    //! Twice the interactions of each atom, plus one if it is the second atom
    DeviceArray<uint32_t>    atomPairs;
    //! On the host, as above
    std::vector<uint32_t>    atomStartHost;
    //! On the host, as above
    std::vector<uint32_t>    atomPairsHost;
    //! The forces on the host
    std::vector<float>       fHost;
    //! The energies on the host
    std::vector<double>      eHost;
};

//! The device memory, never deleted, see gpuQueue
CoulombMemory &coulombMemory()
{
    static CoulombMemory *memory = new CoulombMemory;
    return *memory;
}

/*! \brief The kernels for computeCoulomb, see the description of this file
 * \tparam Forces Whether to compute forces
 * \param[in]  queue     The queue
 * \param[in]  natom     Number of atoms
 * \param[in]  stride    Distance between the rows of x and f
 * \param[in]  nslot     Number of conformations
 * \param[in]  x         The coordinates
 * \param[in]  pairs     The interactions
 * \param[in]  npair     Number of interactions
 * \param[in]  atomStart Start of the interactions of each atom in atomPairs, natom + 1
 * \param[in]  atomPairs Twice the interactions of each atom, plus one if it is the second
 * \param[out] fpair     The force of each interaction, 3*npair*nslot
 * \param[out] epair     The energy of each interaction, npair*nslot
 * \param[out] eatom     The energy of each atom, natom*nslot
 * \param[out] f         The forces, the rows of x
 * \param[out] e         The energies, nslot
 */
template <bool Forces>
void coulombKernels(sycl::queue       &queue,
                    size_t             natom,
                    size_t             stride,
                    size_t             nslot,
                    const int64_t     *x,
                    const CoulombPair *pairs,
                    size_t             npair,
                    const uint32_t    *atomStart,
                    const uint32_t    *atomPairs,
                    float             *fpair,
                    float             *epair,
                    double            *eatom,
                    float             *f,
                    double            *e)
{
    // Each interaction in each conformation
    queue.parallel_for(sycl::range<1>(npair*nslot), [=](sycl::id<1> id)
    {
        size_t             p    = id[0] / nslot;
        size_t             c    = id[0] % nslot;
        const CoulombPair &pair = pairs[p];
        size_t             xi   = 3*pair.i*stride + c;
        size_t             xj   = 3*pair.j*stride + c;
        float              dx[3];
        float              rsq  = 0;
        for (int m = 0; m < 3; m++)
        {
            // The exact difference of the integers, rounded once
            dx[m] = static_cast<float>(x[xi + m*stride] - x[xj + m*stride])*
                static_cast<float>(coordinateUnitNm);
            rsq  += dx[m]*dx[m];
        }
        float rinv  = rsq > 0 ? sycl::rsqrt(rsq) : 0.0f;
        float fscal = 0;
        float ep    = coulombTerm(pair.qq[0], pair.zeta[0], rsq, rinv, &fscal);
        for (int k = 1; k < 4; k++)
        {
            if (pair.qq[k] != 0)
            {
                ep += coulombTerm(pair.qq[k], pair.zeta[k], rsq, rinv, &fscal);
            }
        }
        epair[p*nslot + c] = ep;
        if constexpr (Forces)
        {
            for (int m = 0; m < 3; m++)
            {
                fpair[(3*p + m)*nslot + c] = rsq > 0 ? fscal*dx[m] : 0.0f;
            }
        }
    });
    // The sums for each atom in the order of the interactions, the
    // energies of the interactions in which it is the first atom
    queue.parallel_for(sycl::range<1>(natom*nslot), [=](sycl::id<1> id)
    {
        size_t a     = id[0] / nslot;
        size_t c     = id[0] % nslot;
        double ea    = 0;
        float  fa[3] = { 0, 0, 0 };
        for (uint32_t k = atomStart[a]; k < atomStart[a+1]; k++)
        {
            uint32_t p      = atomPairs[k] / 2;
            bool     second = (atomPairs[k] % 2) == 1;
            if (!second)
            {
                ea += epair[p*nslot + c];
            }
            if constexpr (Forces)
            {
                for (int m = 0; m < 3; m++)
                {
                    float fp = fpair[(3*p + m)*nslot + c];
                    fa[m]   += second ? -fp : fp;
                }
            }
        }
        eatom[a*nslot + c] = ea;
        if constexpr (Forces)
        {
            for (int m = 0; m < 3; m++)
            {
                f[(3*a + m)*stride + c] = fa[m];
            }
        }
    });
    // The energy of each conformation in the order of the atoms
    queue.parallel_for(sycl::range<1>(nslot), [=](sycl::id<1> id)
    {
        size_t c  = id[0];
        double ec = 0;
        for (size_t a = 0; a < natom; a++)
        {
            ec += eatom[a*nslot + c];
        }
        e[c] = ec;
    });
}

/*! \brief The interactions of each atom, in their order
 * \param[in]  natom     Number of atoms
 * \param[in]  pairs     The interactions
 * \param[in]  npair     Number of interactions
 * \param[out] atomStart Start of the interactions of each atom, natom + 1
 * \param[out] atomPairs Twice the interaction, plus one for the second atom
 */
void atomInteractions(size_t                 natom,
                      const CoulombPair     *pairs,
                      size_t                 npair,
                      std::vector<uint32_t> *atomStart,
                      std::vector<uint32_t> *atomPairs)
{
    atomStart->assign(natom + 1, 0);
    for (size_t p = 0; p < npair; p++)
    {
        (*atomStart)[pairs[p].i + 1] += 1;
        (*atomStart)[pairs[p].j + 1] += 1;
    }
    for (size_t a = 0; a < natom; a++)
    {
        (*atomStart)[a+1] += (*atomStart)[a];
    }
    atomPairs->resize(2*npair);
    std::vector<uint32_t> next(atomStart->begin(), atomStart->end() - 1);
    for (size_t p = 0; p < npair; p++)
    {
        (*atomPairs)[next[pairs[p].i]++] = 2*p;
        (*atomPairs)[next[pairs[p].j]++] = 2*p + 1;
    }
}

} // namespace

bool computeCoulomb(size_t             natom,
                    size_t             stride,
                    size_t             nslot,
                    const int64_t     *x,
                    const CoulombPair *pairs,
                    size_t             npair,
                    float             *f,
                    double            *energies,
                    bool               forces)
{
    sycl::queue *queue = gpuQueue();
    if (!queue)
    {
        return false;
    }
    if (npair == 0 || nslot == 0)
    {
        return true;
    }
    auto  &mem = coulombMemory();
    atomInteractions(natom, pairs, npair, &mem.atomStartHost, &mem.atomPairsHost);
    size_t nx  = 3*natom*stride;
    auto   dx  = mem.x.reserve(*queue, nx);
    auto   dp  = mem.pairs.reserve(*queue, npair);
    auto   de  = mem.e.reserve(*queue, nslot);
    auto   das = mem.atomStart.reserve(*queue, natom + 1);
    auto   dap = mem.atomPairs.reserve(*queue, 2*npair);
    auto   dep = mem.epair.reserve(*queue, npair*nslot);
    auto   dea = mem.eatom.reserve(*queue, natom*nslot);
    auto   dfp = forces ? mem.fpair.reserve(*queue, 3*npair*nslot) : nullptr;
    auto   df  = forces ? mem.f.reserve(*queue, nx) : nullptr;
    queue->memcpy(dx, x, nx*sizeof(int64_t));
    queue->memcpy(dp, pairs, npair*sizeof(CoulombPair));
    queue->memcpy(das, mem.atomStartHost.data(), (natom + 1)*sizeof(uint32_t));
    queue->memcpy(dap, mem.atomPairsHost.data(), 2*npair*sizeof(uint32_t));
    if (forces)
    {
        coulombKernels<true>(*queue, natom, stride, nslot, dx, dp, npair, das, dap, dfp, dep, dea, df, de);
    }
    else
    {
        coulombKernels<false>(*queue, natom, stride, nslot, dx, dp, npair, das, dap, dfp, dep, dea, df, de);
    }
    mem.eHost.resize(nslot);
    queue->memcpy(mem.eHost.data(), de, nslot*sizeof(double));
    if (forces)
    {
        mem.fHost.resize(nx);
        queue->memcpy(mem.fHost.data(), df, nx*sizeof(float));
    }
    queue->wait_and_throw();
    for (size_t c = 0; c < nslot; c++)
    {
        energies[c] += mem.eHost[c];
    }
    if (forces)
    {
        for (size_t i = 0; i < 3*natom; i++)
        {
            for (size_t c = 0; c < nslot; c++)
            {
                f[i*stride + c] += mem.fHost[i*stride + c];
            }
        }
    }
    return true;
}

} // namespace gpu

} // namespace alexandria
