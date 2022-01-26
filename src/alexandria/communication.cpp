/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */


#include "communication.h"

#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

const char *cs_name(CommunicationStatus cs)
{
    switch (cs)
    {
        case CS_OK:
            static const char *ok = "Communication OK";
            return ok;
        case CS_ERROR:
            static const char *err = "Communication Error";
            return err;
        case CS_SEND_DATA:
            static const char *sd = "Communication sent data";
            return sd;
        case CS_RECV_DATA:
            static const char *rd = "Communication OK";
            return rd;
        default:
            gmx_fatal(FARGS, "Unknown communication status %d", (int) cs);
    }
    return nullptr;
};


}
