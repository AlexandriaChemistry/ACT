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
 * Implements a limit on the wall time of a program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "time_limit.h"

namespace alexandria
{

namespace
{

//! The moment this file was initialized, which is before main
const TimeLimit::Clock::time_point startOfProgram = TimeLimit::Clock::now();

} // namespace

TimeLimit::TimeLimit(double hours) : TimeLimit(hours, startOfProgram)
{
}

TimeLimit::TimeLimit(double hours, Clock::time_point start)
    : seconds_(hours > 0 ? 3600*hours : 0), start_(start)
{
}

bool TimeLimit::expired() const
{
    if (!limited())
    {
        return false;
    }
    return std::chrono::duration<double>(Clock::now() - start_).count() >= seconds_;
}

double TimeLimit::elapsedHours() const
{
    return std::chrono::duration<double>(Clock::now() - start_).count()/3600;
}

TimeLimit::Clock::time_point TimeLimit::programStart()
{
    return startOfProgram;
}

} // namespace alexandria
