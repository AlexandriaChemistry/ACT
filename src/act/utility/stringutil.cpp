/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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
 
#include "stringutil.h"

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <map>
#include <sstream>
#include <string>

#include "gromacs/utility/fatalerror.h"

int get_option(const char **opts)
{
    int val = 0;

    if (!opts)
    {
        return -1;
    }
    if (opts[val] != nullptr)
    {
        for (val = 1; (opts[val] != nullptr); val++)
        {
            if (strcasecmp(opts[0], opts[val]) == 0)
            {
                break;
            }
        }
    }
    if (opts[val] == nullptr)
    {
        val = 0;
    }
    else
    {
        val--;
    }

    return val;
}

std::vector<std::string> split(const std::string        &s,
                               char                      delim,
                               std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string       item;
    while (std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s,
                               char               delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::string gmx_ftoa(double f)
{
    char buf[32];

    if (fabs(f) < 100)
    {
        sprintf(buf, "%.3f", f);
    }
    else
    {
        sprintf(buf, "%g", f);
    }
    return std::string(buf);
}

std::string gmx_itoa(int f)
{
    char a[32];

    sprintf(a, "%d", f);

    return std::string(a);
}

double my_atof(const char *str, const char *description)
{
    char   *ptr = nullptr;
    double  d   = std::strtod(str, &ptr);
    if (ptr == nullptr)
    {
        fprintf(stderr, "Could not read double precision number %s from '%s' found %f\n",
                description ? description : "", str, d);
        d = -1;
    }
    return d;
}

int my_atoi(const char *str, const char *description)
{
    char *ptr = nullptr;
    int   d   = std::strtol(str, &ptr, 10);
    if (ptr == nullptr)
    {
        fprintf(stderr, "Could not read double precision number %s from '%s' found %d\n",
                description ? description : "", str, d);
        d = -1;
    }
    return d;
}

bool stringEqual(const std::string &a, const std::string &b)
{
    size_t sz = a.size();
    if (b.size() != sz)
    {
        return false;
    }
    for (size_t i = 0; i < sz; ++i)
    {
        if (tolower(a[i]) != tolower(b[i]))
        {
            return false;
        }
    }
    return true;
}

