/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "basenetwork.h"

#include "config.h"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"

bool gmx_mpi_initialized()
{
#if !GMX_MPI
    return 0;
#else
    int n;
    MPI_Initialized(&n);

    return n != 0;
#endif
}

int gmx_node_num()
{
#if !GMX_MPI
    return 1;
#else
#if GMX_THREAD_MPI
    if (!gmx_mpi_initialized())
    {
        return 1;
    }
#endif
    int i;
    (void) MPI_Comm_size(MPI_COMM_WORLD, &i);
    return i;
#endif
}

int gmx_node_rank()
{
#if !GMX_MPI
    return 0;
#else
#if GMX_THREAD_MPI
    if (!gmx_mpi_initialized())
    {
        return 0;
    }
#endif
    int i;
    (void) MPI_Comm_rank(MPI_COMM_WORLD, &i);
    return i;
#endif
}

static int mpi_hostname_hash()
{
    int hash_int;

#if GMX_LIB_MPI
    int  resultlen;
    char mpi_hostname[MPI_MAX_PROCESSOR_NAME];

    /* This procedure can only differentiate nodes with different names.
     * Architectures where different physical nodes have identical names,
     * such as IBM Blue Gene, should use an architecture specific solution.
     */
    MPI_Get_processor_name(mpi_hostname, &resultlen);

    /* The string hash function returns an unsigned int. We cast to an int.
     * Negative numbers are converted to positive by setting the sign bit to 0.
     * This makes the hash one bit smaller.
     * A 63-bit hash (with 64-bit int) should be enough for unique node hashes,
     * even on a million node machine. 31 bits might not be enough though!
     */
    hash_int = static_cast<int>(gmx_string_fullhash_func(mpi_hostname, gmx_string_hash_init));
    if (hash_int < 0)
    {
        hash_int -= INT_MIN;
    }
#else

    /* thread-MPI currently puts the thread number in the process name,
     * we might want to change this, as this is inconsistent with what
     * most MPI implementations would do when running on a single node.
     */
    hash_int = 0;
#endif

    return hash_int;
}

int gmx_physicalnode_id_hash()
{
    int hash = 0;

    if (GMX_MPI)
    {
        hash = mpi_hostname_hash();
    }

    if (debug)
    {
        fprintf(debug, "In gmx_physicalnode_id_hash: hash %d\n", hash);
    }

    return hash;
}

void gmx_broadcast_world(int size, void *buffer)
{
#if GMX_MPI
    MPI_Bcast(buffer, size, MPI_BYTE, 0, MPI_COMM_WORLD);
#else
    GMX_UNUSED_VALUE(size);
    GMX_UNUSED_VALUE(buffer);
#endif
}

#if GMX_LIB_MPI
void gmx_abort(int errorno)
{
    MPI_Abort(MPI_COMM_WORLD, errorno);
    std::abort();
}
#endif
