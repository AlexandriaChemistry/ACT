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
 * Integer coordinates.
 *
 * The coordinates of the conformations are stored as integers in units of
 * coordinateUnit, 2^-56 nm or about 1.4e-17 nm, such that differences of
 * coordinates are exact, independent of where the particles are, and
 * converted to floating point only once, with one rounding. Since the
 * unit is a power of two, scaling a converted difference is exact. The
 * resolution is that of double precision coordinates of about 1 nm, such
 * that also very stiff shells can be relaxed as far as in double precision.
 * A 64-bit integer covers coordinates up to 128 nm from the origin.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ACT_BASICS_COORDINATE_H
#define ACT_BASICS_COORDINATE_H

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

//! A coordinate component in units of coordinateUnit
using Coordinate = std::int64_t;

//! The number of units per nm, a power of two
constexpr double coordinateScale = 72057594037927936.0; // 2^56

//! The unit of the integer coordinates in nm
constexpr double coordinateUnit  = 1/coordinateScale;

//! The largest magnitude of a coordinate in nm, 2^7
constexpr double coordinateMax   = 128.0;

/*! \brief Report a coordinate beyond the range, kept out of line
 * \param[in] x The coordinate in nm
 * \throws gmx::RangeError always
 */
[[noreturn, gnu::cold, gnu::noinline]] inline void coordinateRangeError(double x)
{
    GMX_THROW(gmx::RangeError(gmx::formatString("Coordinate %g nm is not within %g nm from the origin",
                                                x, coordinateMax)));
}

/*! \brief Convert a coordinate in nm to integer units
 * The coordinate is rounded to the nearest unit, halfway cases to even,
 * with the default rounding mode, which std::rint uses without the cost
 * of std::nearbyint.
 * \param[in] x The coordinate in nm
 * \return the coordinate in units of coordinateUnit
 * \throws gmx::RangeError if the coordinate is not finite or too far from the origin
 */
inline Coordinate toCoordinate(double x)
{
    double units = std::rint(x*coordinateScale);
    // Below 2^63, the largest magnitude of the integers
    if (!(std::abs(units) < coordinateMax*coordinateScale)) [[unlikely]]
    {
        coordinateRangeError(x);
    }
    return static_cast<Coordinate>(units);
}

/*! \brief Convert a coordinate in nm to integer units, limited to the range
 * As toCoordinate, but a coordinate beyond the range becomes the largest
 * one, for moves of particles that may run away, such as shells.
 * \param[in] x The coordinate in nm, not NaN
 * \return the coordinate in units of coordinateUnit
 */
inline Coordinate toCoordinateClamped(double x)
{
    const double limit = 0.999*coordinateMax;
    return toCoordinate(std::clamp(x, -limit, limit));
}

/*! \brief Convert an integer coordinate to nm
 * \param[in] c The coordinate in units of coordinateUnit
 * \return the coordinate in nm
 */
inline double fromCoordinate(Coordinate c)
{
    return c*coordinateUnit;
}

/*! \brief Convert an integer coordinate to nm in the precision of T
 * The integer is converted to T once, the unit is a power of two, such
 * that the result is that of fromCoordinate rounded to T. In single
 * precision this saves converting to double first.
 * \param[in] c The coordinate in units of coordinateUnit
 * \return the coordinate in nm
 */
template <typename T>
inline T fromCoordinateAs(Coordinate c)
{
    return static_cast<T>(c)*static_cast<T>(coordinateUnit);
}

/*! \brief The difference of two coordinates
 * Computed with wrap-around, such that it is correct if the coordinates
 * are less than 2^63 units apart, wherever they are.
 * \param[in] a The first coordinate
 * \param[in] b The second coordinate
 * \return a - b
 */
inline Coordinate coordinateDifference(Coordinate a, Coordinate b)
{
    return static_cast<Coordinate>(static_cast<std::uint64_t>(a) - static_cast<std::uint64_t>(b));
}

/*! \brief The sum of two coordinates, or of a coordinate and a difference
 * Computed with wrap-around, such that it is the inverse of
 * coordinateDifference.
 * \param[in] a The first coordinate
 * \param[in] b The second coordinate
 * \return a + b
 */
inline Coordinate coordinateSum(Coordinate a, Coordinate b)
{
    return static_cast<Coordinate>(static_cast<std::uint64_t>(a) + static_cast<std::uint64_t>(b));
}

/*! \brief Move a coordinate, limited to the range
 * The displacement is rounded and added to the integer, which gives the
 * same as rounding the moved coordinate, without converting the
 * coordinate. A move beyond the range ends at the range, for particles
 * that may run away, such as shells.
 * \param[in] c The coordinate
 * \param[in] d The displacement in nm
 * \return the moved coordinate
 * \throws gmx::RangeError if the displacement is not a number
 */
inline Coordinate moveCoordinate(Coordinate c, double d)
{
    // Displacements of less than a nm, as when relaxing shells, fit in the
    // integers without checking the range
    Coordinate step = (std::abs(d) < 1) ? static_cast<Coordinate>(std::rint(d*coordinateScale)) : toCoordinateClamped(d);
    Coordinate moved;
    if (__builtin_add_overflow(c, step, &moved))
    {
        moved = (d > 0) ? toCoordinate(0.999*coordinateMax) : toCoordinate(-0.999*coordinateMax);
    }
    return moved;
}

} // namespace alexandria

#endif
