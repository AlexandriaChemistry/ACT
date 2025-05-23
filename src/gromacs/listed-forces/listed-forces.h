/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \defgroup module_listed-forces Interactions between lists of particles
 * \ingroup group_mdrun
 *
 * \brief Handles computing energies and forces for listed
 * interactions.
 *
 * Located here is the the code for
 * - computing energies and forces for interactions between a small
     number of particles, e.g bonds, position restraints and listed
     non-bonded interactions (e.g. 1-4).
 * - high-level functions used by mdrun for computing a set of such
     quantities
 * - managing thread-wise decomposition, thread-local buffer output,
     and reduction of output data across threads.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 */
/*! \libinternal \file
 *
 * \brief This file contains declarations of high-level functions used
 * by mdrun to compute energies and forces for listed interactions.
 *
 * Clients of libactgromacs that want to evaluate listed interactions
 * should call functions declared here.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_listed-forces
 */
#ifndef GMX_LISTED_FORCES_LISTED_FORCES_H
#define GMX_LISTED_FORCES_LISTED_FORCES_H

#include "gromacs/gpu_utils/gpu_macros.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_enerdata_t;
struct gmx_grppairener_t;
struct gmx_multisim_t;
struct gmx_ffparams_t;
struct GpuBondedLists;
class history_t;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_idef;
struct t_graph;
struct t_lambda;
struct t_mdatoms;
struct t_nrnb;
class t_state;

namespace gmx
{
class ForceWithVirial;
}

//! Type of CPU function to compute a bonded interaction.
using BondedFunction = real(*)(int nbonds, const t_iatom iatoms[],
                               const t_iparams iparams[],
                               const rvec x[], rvec4 f[], rvec fshift[],
                               const t_pbc *pbc, const t_graph *g,
                               real lambda, real *dvdlambda,
                               const t_mdatoms *md, t_fcdata *fcd,
                               int *ddgatindex);

//! Getter for finding a callable CPU function to compute an \c ftype interaction.
BondedFunction bondedFunction(int ftype);

/*! \brief Calculates all listed force interactions.
 *
 * Note that pbc_full is used only for position restraints, and is
 * not initialized if there are none. */
void calc_listed(const t_commrec *cr,
                 const gmx_multisim_t *ms,
                 struct gmx_wallcycle *wcycle,
                 const t_idef *idef,
                 const rvec x[], history_t *hist,
                 rvec f[],
                 gmx::ForceWithVirial *forceWithVirial,
                 const t_forcerec *fr,
                 const struct t_pbc *pbc, const struct t_pbc *pbc_full,
                 const struct t_graph *g,
                 gmx_enerdata_t *enerd, t_nrnb *nrnb, const real *lambda,
                 const t_mdatoms *md,
                 struct t_fcdata *fcd, int *ddgatindex,
                 int force_flags);

/*! \brief As calc_listed(), but only determines the potential energy
 * for the perturbed interactions.
 *
 * The shift forces in fr are not affected. */
void calc_listed_lambda(const t_idef *idef,
                        const rvec x[],
                        const t_forcerec *fr,
                        const struct t_pbc *pbc, const struct t_graph *g,
                        gmx_grppairener_t *grpp, real *epot, t_nrnb *nrnb,
                        const real *lambda,
                        const t_mdatoms *md,
                        struct t_fcdata *fcd, int *global_atom_index);

/*! \brief Do all aspects of energy and force calculations for mdrun
 * on the set of listed interactions */
void
do_force_listed(struct gmx_wallcycle           *wcycle,
                matrix                          box,
                const t_lambda                 *fepvals,
                const t_commrec                *cr,
                const gmx_multisim_t           *ms,
                const t_idef                   *idef,
                const rvec                      x[],
                history_t                      *hist,
                rvec                           *forceForUseWithShiftForces,
                gmx::ForceWithVirial           *forceWithVirial,
                const t_forcerec               *fr,
                const struct t_pbc             *pbc,
                const struct t_graph           *graph,
                gmx_enerdata_t                 *enerd,
                t_nrnb                         *nrnb,
                const real                     *lambda,
                const t_mdatoms                *md,
                struct t_fcdata                *fcd,
                int                            *global_atom_index,
                int                             flags);

/*! \brief Initializes the GPU bonded setup */
CUDA_FUNC_QUALIFIER
void
init_gpu_bonded(GpuBondedLists gmx_unused       *gpuBondedLists,
                const gmx_ffparams_t gmx_unused &ffparams,
                void gmx_unused                 *streamPtr) CUDA_FUNC_TERM

/*! \brief Updates the bonded work to run on a GPU
 *
 * Intended to be called after each domain decomposition stage. */
CUDA_FUNC_QUALIFIER
void update_gpu_bonded(GpuBondedLists gmx_unused *gpuBondedLists) CUDA_FUNC_TERM

/*! \brief Launches bonded kernels on a GPU */
CUDA_FUNC_QUALIFIER
void do_bonded_gpu(t_forcerec gmx_unused   *fr,
                   int gmx_unused           forceFlags,
                   void gmx_unused         *xqDevicePtr,
                   const matrix gmx_unused  box,
                   void gmx_unused         *forceDevicePtr,
                   rvec gmx_unused         *fshiftDevicePtr) CUDA_FUNC_TERM

/*! \brief Copies back the bonded energies */
CUDA_FUNC_QUALIFIER
void bonded_gpu_get_energies(t_forcerec gmx_unused     *fr,
                             gmx_enerdata_t gmx_unused *enerd) CUDA_FUNC_TERM

/*! \brief Clears the device side energy buffer */
CUDA_FUNC_QUALIFIER
void bonded_gpu_clear_energies(GpuBondedLists gmx_unused *gpuBondedLists) CUDA_FUNC_TERM

#endif
