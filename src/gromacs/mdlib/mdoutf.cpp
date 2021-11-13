/*
 * This file is part of the GROMACS molecular simulation package.
 *
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

#include "mdoutf.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/trajectory_writing.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

struct gmx_mdoutf {
    t_fileio               *fp_trn;
    ener_file_t             fp_ene;
    int                     eIntegrator;
    gmx_bool                bExpanded;
    int                     elamstats;
    int                     simulation_part;
    FILE                   *fp_dhdl;
    int                     natoms_global;
    int                     natoms_x_compressed;
    gmx_groups_t           *groups; /* for compressed position writing */
    gmx_wallcycle_t         wcycle;
    rvec                   *f_global;
};


gmx_mdoutf_t init_mdoutf(int nfile, const t_filenm fnm[],
                         const MdrunOptions &mdrunOptions,
                         const t_commrec *cr,
                         const t_inputrec *ir, gmx_mtop_t *top_global,
                         const gmx_output_env_t *oenv, gmx_wallcycle_t wcycle)
{
    gmx_mdoutf_t   of;
    const char    *appendMode = "a+", *writeMode = "w+", *filemode;
    gmx_bool       bAppendFiles;
    int            i;

    snew(of, 1);

    of->fp_trn       = nullptr;
    of->fp_ene       = nullptr;
    of->fp_dhdl      = nullptr;

    of->eIntegrator             = ir->eI;
    of->bExpanded               = ir->bExpanded;
    of->elamstats               = ir->expandedvals->elamstats;
    of->simulation_part         = ir->simulation_part;
    of->wcycle                  = wcycle;
    of->f_global                = nullptr;

    if (MASTER(cr))
    {
        bAppendFiles = mdrunOptions.continuationOptions.appendFiles;

        filemode = bAppendFiles ? appendMode : writeMode;

        if ((EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI)) &&
            (!GMX_FAHCORE &&
             !(EI_DYNAMICS(ir->eI) &&
               ir->nstxout == 0 &&
               ir->nstvout == 0 &&
               ir->nstfout == 0)
            )
            )
        {
            const char *filename;
            filename = ftp2fn(efTRN, nfile, fnm);
            switch (fn2ftp(filename))
            {
                case efTRR:
                case efTRN:
                    /* If there is no uncompressed coordinate output and
                       there is compressed TNG output write forces
                       and/or velocities to the TNG file instead. */
                    if (ir->nstxout != 0)
                    {
                        of->fp_trn = gmx_trr_open(filename, filemode);
                    }
                    break;
                default:
                    gmx_incons("Invalid full precision file format");
            }
        }
        if (EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI))
        {
            of->fp_ene = open_enx(ftp2fn(efEDR, nfile, fnm), filemode);
        }

        if ((ir->efep != efepNO || ir->bSimTemp) && ir->fepvals->nstdhdl > 0 &&
            (ir->fepvals->separate_dhdl_file == esepdhdlfileYES ) &&
            EI_DYNAMICS(ir->eI))
        {
            if (bAppendFiles)
            {
                of->fp_dhdl = gmx_fio_fopen(opt2fn("-dhdl", nfile, fnm), filemode);
            }
            else
            {
                of->fp_dhdl = open_dhdl(opt2fn("-dhdl", nfile, fnm), ir, oenv);
            }
        }

        /* Set up atom counts so they can be passed to actual
           trajectory-writing routines later. Also, XTC writing needs
           to know what (and how many) atoms might be in the XTC
           groups, and how to look up later which ones they are. */
        of->natoms_global       = top_global->natoms;
        of->groups              = &top_global->groups;
        of->natoms_x_compressed = 0;
        for (i = 0; (i < top_global->natoms); i++)
        {
            if (getGroupType(*of->groups, egcCompressedX, i) == 0)
            {
                of->natoms_x_compressed++;
            }
        }

        if (ir->nstfout && DOMAINDECOMP(cr))
        {
            snew(of->f_global, top_global->natoms);
        }
    }

    return of;
}

ener_file_t mdoutf_get_fp_ene(gmx_mdoutf_t of)
{
    return of->fp_ene;
}

FILE *mdoutf_get_fp_dhdl(gmx_mdoutf_t of)
{
    return of->fp_dhdl;
}

gmx_wallcycle_t mdoutf_get_wcycle(gmx_mdoutf_t of)
{
    return of->wcycle;
}

void mdoutf_write_to_trajectory_files(FILE *fplog, const t_commrec *cr,
                                      gmx_mdoutf_t of,
                                      int mdof_flags,
                                      gmx_mtop_t *top_global,
                                      int64_t step, double t,
                                      t_state *state_local, t_state *state_global,
                                      ObservablesHistory *observablesHistory,
                                      gmx::ArrayRef<gmx::RVec> f_local)
{
    rvec *f_global;

    if (DOMAINDECOMP(cr))
    {
        {
            if (mdof_flags & (MDOF_X | MDOF_X_COMPRESSED))
            {
                gmx::ArrayRef<gmx::RVec> globalXRef = MASTER(cr) ? makeArrayRef(state_global->x) : gmx::EmptyArrayRef();
                dd_collect_vec(cr->dd, state_local, makeArrayRef(state_local->x), globalXRef);
            }
            if (mdof_flags & MDOF_V)
            {
                gmx::ArrayRef<gmx::RVec> globalVRef = MASTER(cr) ? makeArrayRef(state_global->v) : gmx::EmptyArrayRef();
                dd_collect_vec(cr->dd, state_local, makeArrayRef(state_local->v), globalVRef);
            }
        }
        f_global = of->f_global;
        if (mdof_flags & MDOF_F)
        {
            dd_collect_vec(cr->dd, state_local, f_local, gmx::arrayRefFromArray(reinterpret_cast<gmx::RVec *>(f_global), f_local.size()));
        }
    }
    else
    {
        /* We have the whole state locally: copy the local state pointer */
        state_global = state_local;

        f_global     = as_rvec_array(f_local.data());
    }

    if (MASTER(cr))
    {
        if (mdof_flags & (MDOF_X | MDOF_V | MDOF_F))
        {
            const rvec *x = (mdof_flags & MDOF_X) ? state_global->x.rvec_array() : nullptr;
            const rvec *v = (mdof_flags & MDOF_V) ? state_global->v.rvec_array() : nullptr;
            const rvec *f = (mdof_flags & MDOF_F) ? f_global : nullptr;

            if (of->fp_trn)
            {
                gmx_trr_write_frame(of->fp_trn, step, t, state_local->lambda[efptFEP],
                                    state_local->box, top_global->natoms,
                                    x, v, f);
                if (gmx_fio_flush(of->fp_trn) != 0)
                {
                    gmx_file("Cannot write trajectory; maybe you are out of disk space?");
                }
            }
        }
    }
}

void done_mdoutf(gmx_mdoutf_t of)
{
    if (of->fp_ene != nullptr)
    {
        close_enx(of->fp_ene);
    }
    if (of->fp_trn)
    {
        gmx_trr_close(of->fp_trn);
    }
    if (of->fp_dhdl != nullptr)
    {
        gmx_fio_fclose(of->fp_dhdl);
    }
    if (of->f_global != nullptr)
    {
        sfree(of->f_global);
    }

    sfree(of);
}

