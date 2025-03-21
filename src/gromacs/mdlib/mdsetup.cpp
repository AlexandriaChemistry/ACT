/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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
#include "actpre.h"

#include "mdsetup.h"

#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

/* TODO: Add a routine that collects the initial setup of the algorithms.
 *
 * The final solution should be an MD algorithm base class with methods
 * for initialization and atom-data setup.
 */

void mdAlgorithmsSetupAtomData(const t_inputrec  *ir,
                               const gmx_mtop_t  *top_global,
                               gmx_localtop_t    *top,
                               t_forcerec        *fr,
                               t_graph          **graph,
                               gmx::MDAtoms      *mdAtoms,
                               gmx_vsite_t       *vsite,
                               gmx_shellfc_t     *shellfc)
{
    int   numAtomIndex, numHomeAtoms;
    int  *atomIndex;

    {
        numAtomIndex = -1;
        atomIndex    = nullptr;
        numHomeAtoms = top_global->natoms;
    }
    atoms2md(top_global, ir, numAtomIndex, atomIndex, numHomeAtoms, mdAtoms);

    auto mdatoms = mdAtoms->mdatoms();
    {
        /* Currently gmx_generate_local_top allocates and returns a pointer.
         * We should implement a more elegant solution.
         */
        gmx_localtop_t *tmpTop;

        tmpTop = gmx_mtop_generate_local_top(top_global, ir->efep != efepNO);
        *top   = *tmpTop;
        sfree(tmpTop);
    }

    if (vsite)
    {
        set_vsite_top(vsite, top, mdatoms);
    }

    if (ir->ePBC != epbcNONE && !fr->bMolPBC)
    {
        GMX_ASSERT(graph != nullptr, "We use a graph with PBC (no periodic mols) and without DD");

        *graph = mk_graph(nullptr, &(top->idef), 0, top_global->natoms, FALSE, FALSE);
    }

    /* Note that with DD only flexible constraints, not shells, are supported
     * and these don't require setup in make_local_shells().
     */
    if (shellfc)
    {
        make_local_shells(mdatoms, shellfc);
    }

    setup_bonded_threading(fr->bondedThreading,
                           fr->natoms_force,
                           fr->gpuBondedLists != nullptr,
                           top->idef);

}
