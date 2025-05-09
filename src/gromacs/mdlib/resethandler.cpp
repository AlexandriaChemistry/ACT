/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Defines the reset handler class.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */

#include "actpre.h"

#include "resethandler.h"

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{
ResetHandler::ResetHandler(
        compat::not_null<SimulationSignal*> signal,
        bool                                simulationsShareState,
        int64_t                             nsteps,
        bool                                isMaster,
        bool                                resetHalfway,
        real                                maximumHoursToRun,
        const MDLogger                     &mdlog,
        gmx_wallcycle_t                     wcycle,
        gmx_walltime_accounting_t           walltime_accounting) :
    signal_(*signal),
    rankCanSetSignal_(false),
    simulationNeedsReset_(false),
    maximumHoursToRun_(maximumHoursToRun)
{
    if (simulationsShareState)
    {
        signal_.isLocal = false;
    }
    if (resetHalfway)
    {
        GMX_LOG(mdlog.info).asParagraph().
            appendText(
                "The -resethway functionality is deprecated, and may be removed in a future version.");
        if (nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle, nsteps / 2);
        }
        simulationNeedsReset_ = true;

        if (isMaster && (maximumHoursToRun > 0))
        {
            rankCanSetSignal_ = true;
        }
    }
    else if (wcycle_get_reset_counters(wcycle) > 0)
    {
        simulationNeedsReset_ = true;
    }
    else
    {
        // if no reset is happening, this will always be valid
        walltime_accounting_set_valid_finish(walltime_accounting);
    }
}

bool ResetHandler::setSignalImpl(gmx_walltime_accounting_t walltime_accounting)
{
    const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);
    if (secondsSinceStart > maximumHoursToRun_ * 60.0 * 60.0 * 0.495)
    {
        /* Set flag that will communicate the signal to all ranks in the simulation */
        signal_.sig = static_cast<signed char>(ResetSignal::doResetCounters);
        /* Let handler know that we did signal a reset */
        return true;
    }
    /* Let handler know that we did not signal a reset */
    return false;
}

bool ResetHandler::resetCountersImpl(
        int64_t                     step,
        int64_t                     step_rel,
        const MDLogger             &mdlog,
        FILE                       *fplog,
        const t_commrec            *cr,
        t_nrnb                     *nrnb,
        gmx_wallcycle_t             wcycle,
        gmx_walltime_accounting_t   walltime_accounting)
{
    /* Reset either if signal has been passed,  */
    if (static_cast<ResetSignal>(signal_.set) == ResetSignal::doResetCounters ||
        step_rel == wcycle_get_reset_counters(wcycle))
    {
        char sbuf[STEPSTRSIZE];

        /* Reset all the counters related to performance over the run */
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "step %s: resetting all time and cycle counters",
                gmx_step_str(step, sbuf));

        wallcycle_stop(wcycle, ewcRUN);
        wallcycle_reset_all(wcycle);
        init_nrnb(nrnb);
        wallcycle_start(wcycle, ewcRUN);
        walltime_accounting_reset_time(walltime_accounting, step);
        print_date_and_time(fplog, cr->nodeid, "Restarted time", gmx_gettime());

        wcycle_set_reset_counters(wcycle, -1);
        /* Reset can only happen once, so clear the triggering flag. */
        signal_.set = static_cast<signed char>(ResetSignal::noSignal);
        /* We have done a reset, so the finish will be valid. */
        walltime_accounting_set_valid_finish(walltime_accounting);
        /* Let handler know that we handled a reset */
        return true;
    }

    /* Let handler know that we did not handle a reset */
    return false;
}

}  // namespace gmx
