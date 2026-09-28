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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#ifndef ALEXANDRIA_SIMPSON_H
#define ALEXANDRIA_SIMPSON_H

#include <vector>

#include "act/basics/msg_handler.h"

/*! \brief Use Simpson's rule to integrate a function
 * https://en.wikipedia.org/wiki/Simpson%27s_rule
 * There are two prerequisites for this to work:
 * 1) the spacing in the x-coordinates must be uniform
 * 2) the size of the x and y array should be odd, or zeroPadding should be true
 * \param[in] msghandler  For error handling
 * \param[in] x           X coordinate
 * \param[in] y           Y coordinate
 * \param[in] zeroPadding allow to add a zero value at the end of the vector if the number of point is even
 * \return the integral
 */
double simpsonIntegrate(alexandria::MsgHandler    *msghandler,
                        const std::vector<double> &x,
                        const std::vector<double> &y,
                        bool                       zeroPadding = true);
                        
#endif
