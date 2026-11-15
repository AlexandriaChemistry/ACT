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
 * How evenly the middlemen of train_ff share the individuals of a
 * generation when the master hands them out from a queue, see
 * HybridGAMC and -lockstep.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_ALEXANDRIA_QUEUEBALANCE_H
#define ACT_ALEXANDRIA_QUEUEBALANCE_H

#include <algorithm>
#include <cstddef>

namespace alexandria
{

//! How the groups of individuals of a generation are spread over the middlemen
struct QueueBalance
{
    //! The number of groups of at most -lockstep individuals
    size_t groups = 0;
    //! The number of groups that the busiest middleman computes in a generation
    size_t rounds = 0;
    //! The fraction of the time of a generation that the middlemen wait, on average
    double idle = 0;
    //! The largest number of middlemen below the given one with no waiting, 0 if none
    size_t fewer = 0;
    //! The smallest number of middlemen above the given one with no waiting
    size_t more = 0;
};

/*! \brief How evenly the middlemen share the individuals of a generation
 * The master hands out groups of lockstep individuals, a next group to the
 * first middleman that is done. If the groups take equally long, a
 * generation takes as many group times as the busiest middleman computes
 * groups, while the others wait for it.
 * \param[in] population The number of individuals
 * \param[in] lockstep   The number of individuals in a group
 * \param[in] nmiddlemen The number of middlemen
 * \return the balance
 */
inline QueueBalance queueBalance(size_t population, size_t lockstep, size_t nmiddlemen)
{
    QueueBalance qb;
    if (population == 0 || lockstep == 0 || nmiddlemen == 0)
    {
        return qb;
    }
    qb.groups = (population + lockstep - 1)/lockstep;
    qb.rounds = (qb.groups + nmiddlemen - 1)/nmiddlemen;
    qb.idle   = 1.0 - static_cast<double>(qb.groups)/(qb.rounds*nmiddlemen);
    // With m middlemen nobody waits if m divides the number of groups, or
    // if there are at least as many middlemen as groups
    for (size_t m = nmiddlemen - 1; m > 0 && qb.fewer == 0; m--)
    {
        if (qb.groups % m == 0)
        {
            qb.fewer = m;
        }
    }
    for (size_t m = nmiddlemen + 1; qb.more == 0; m++)
    {
        if (qb.groups % m == 0 || m >= qb.groups)
        {
            qb.more = m;
        }
    }
    return qb;
}

//! The settings that -lockstep 0 chooses, see planLockstep
struct LockstepPlan
{
    //! The number of individuals that a middleman computes at once
    size_t lockstep   = 1;
    //! The number of middlemen
    size_t middlemen  = 0;
    //! The number of processes that work, the middlemen and the master
    size_t processes  = 0;
    //! The number of processes that are left idle
    size_t idle       = 0;
    //! The number of groups that each middleman computes in a generation
    size_t rounds     = 0;
    //! The number of processes with which each middleman would compute one group, 0 if that is the plan
    size_t better     = 0;
};

//! The largest number of individuals that a middleman computes at once with -lockstep 0
constexpr size_t maxAutoLockstep = 16;

//! The number of middlemen per GPU that -lockstep 0 aims for at least
constexpr size_t minMiddlemenPerGpu = 2;

//! The number of middlemen per GPU that -lockstep 0 never exceeds, below the 48 processes that the NVIDIA multi-process service takes
constexpr size_t maxMiddlemenPerGpu = 40;

/*! \brief Choose -lockstep and the number of middlemen for training on GPUs
 * The middlemen aim for at least minMiddlemenPerGpu per GPU, and as many as
 * needed for groups of at most maxAutoLockstep individuals, but at most
 * maxMiddlemenPerGpu per GPU and at most the processes there are minus the
 * master. The individuals are divided into groups of equal size, and the
 * middlemen are the fewest that compute them in the same number of rounds;
 * the processes left over are idle. On the computers of the benchmarks in
 * the manual this gives the settings that were the fastest there, for
 * instance for a population of 128: groups of 16 on 8 middlemen with 4
 * GPUs, and of 8 on 16 middlemen with 8 GPUs.
 * \param[in] population The number of individuals
 * \param[in] processes  The number of processes
 * \param[in] gpus       The number of GPUs of all nodes
 * \return the plan, with lockstep 1 and all processes working if there is
 *         no GPU or no more than one process
 */
inline LockstepPlan planLockstep(size_t population, size_t processes, size_t gpus)
{
    LockstepPlan plan;
    plan.processes = processes;
    plan.middlemen = processes > 1 ? processes - 1 : 0;
    if (population == 0 || processes < 2 || gpus == 0)
    {
        return plan;
    }
    // The middlemen of a plan with at most m middlemen
    const auto make = [population](size_t m)
    {
        LockstepPlan p;
        p.lockstep  = std::min(maxAutoLockstep, (population + m - 1)/m);
        size_t groups = (population + p.lockstep - 1)/p.lockstep;
        p.rounds    = (groups + m - 1)/m;
        p.middlemen = (groups + p.rounds - 1)/p.rounds;
        p.processes = p.middlemen + 1;
        return p;
    };
    size_t wanted = std::max((population + maxAutoLockstep - 1)/maxAutoLockstep,
                             minMiddlemenPerGpu*gpus);
    wanted        = std::min(wanted, maxMiddlemenPerGpu*gpus);
    plan          = make(std::min(wanted, processes - 1));
    plan.idle     = processes - plan.processes;
    if (wanted > processes - 1)
    {
        const auto ideal = make(wanted);
        if (ideal.rounds < plan.rounds)
        {
            plan.better = ideal.processes;
        }
    }
    return plan;
}

} // namespace alexandria

#endif
