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

} // namespace alexandria

#endif
