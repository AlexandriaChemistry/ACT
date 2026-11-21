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
 * A limit on the wall time of a program.
 *
 * A calculation that may take longer than the time it is given, for
 * instance in a queue, asks now and then whether its time is up and
 * then ends in an orderly way, with its results written. The time counts
 * from the start of the program, such that every process of a parallel
 * run has the same limit without talking to the others.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_UTILITY_TIME_LIMIT_H
#define ACT_UTILITY_TIME_LIMIT_H

#include <chrono>

namespace alexandria
{

/*! \brief A limit on the wall time since the program started
 * Without a limit expired() never holds, and costs nothing but a comparison.
 */
class TimeLimit
{
public:
    //! The clock, which does not jump when the time of the machine is set
    using Clock = std::chrono::steady_clock;

    /*! \brief A limit counted from the start of the program
     * \param[in] hours The time in hours after which expired() holds, no limit if zero or less
     */
    explicit TimeLimit(double hours = 0);

    /*! \brief A limit counted from a given moment
     * \param[in] hours The time in hours after which expired() holds, no limit if zero or less
     * \param[in] start The moment the time counts from
     */
    TimeLimit(double hours, Clock::time_point start);

    //! \return Whether there is a limit
    bool limited() const { return seconds_ > 0; }

    //! \return The limit in hours, zero if there is none
    double hours() const { return limited() ? seconds_/3600 : 0; }

    //! \return Whether there is a limit and the time is up
    bool expired() const;

    //! \return The time in hours since the moment the limit counts from
    double elapsedHours() const;

    /*! \brief When the program started
     * \return The moment this library was loaded, before main
     */
    static Clock::time_point programStart();

private:
    //! The limit in seconds, none if zero or less
    double            seconds_ = 0;
    //! The moment the time counts from
    Clock::time_point start_;
};

} // namespace alexandria

#endif
