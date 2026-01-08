/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2025,2026
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
#include "actpre.h" 

#include "version.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/stringutil.h"

std::string act_welcome()
{
    return gmx::formatString("\n     Welcome to the Alexandria Chemistry Toolkit\n\n"
                             "       Version %s\n\n"
                             "              Copyright (c) 2014-2026\n\n"
                             "Mohammad M. Ghahremanpour, Paul J. van Maaren and David van der Spoel\n\n"
                             "See https://github.com/AlexandriaChemistry/ACT for details.\n\n"
                             "The Alexandria Chemistry Toolkit (ACT) is free software under the GNU Public License v 2.\n"
                             "Read more at https://www.gnu.org/licenses/gpl-2.0.html.\n\n"
                             "Parts of ACT are under the GNU Lesser Public License 2.1 as indicated\n"
                             "in the individual files.\n\n", act_version());
}

std::string act_goodbye()
{
    return gmx::formatString("Thanks for using the Alexandria Chemistry Toolkit.\n"
                             "    https://github.com/AlexandriaChemistry/ACT\n");
}
