/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "xml_util.h"

#include <string.h>

#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

void add_xml_int(xmlNodePtr ptr, const std::string &name, int val)
{
    auto buf = gmx::formatString("%d", val);
    if (xmlSetProp(ptr, (xmlChar *)name.c_str(), (xmlChar *)buf.c_str()) == 0)
    {
        gmx_fatal(FARGS, "XML problem setting %s to %d", name.c_str(), val);
    }
}

void add_xml_double(xmlNodePtr ptr, const std::string &name, double val)
{
    std::string buf = gmx::formatString("%.13g", val);
    if (xmlSetProp(ptr, (xmlChar *)name.c_str(), (xmlChar *)buf.c_str()) == 0)
    {
        gmx_fatal(FARGS, "XML problem setting %s to %g", name.c_str(), val);
    }
}

void add_xml_char(xmlNodePtr ptr, const std::string &name, const char *val)
{
    if (xmlSetProp(ptr, (xmlChar *)name.c_str(), (xmlChar *)val) == 0)
    {
        gmx_fatal(FARGS, "XML problem setting %s to %s", name.c_str(), val);
    }
}

void add_xml_string(xmlNodePtr ptr, const std::string &name, const std::string &val)
{
    if (xmlSetProp(ptr, xmlCharStrdup(name.c_str()), xmlCharStrdup(val.c_str())) == 0)
    {
        gmx_fatal(FARGS, "XML problem setting %s to %s", name.c_str(), val.c_str());
    }
}

xmlNodePtr add_xml_child(xmlNodePtr parent, const std::string &type)
{
    xmlNodePtr child;

    if ((child = xmlNewChild(parent, NULL, (xmlChar *)type.c_str(), NULL)) == NULL)
    {
        gmx_fatal(FARGS, "Creating element %s", type.c_str());
    }

    return child;
}

xmlNodePtr add_xml_child_val(xmlNodePtr parent, const std::string &type, const char *value)
{
    xmlNodePtr child;

    if ((child = xmlNewChild(parent, NULL, (xmlChar *)type.c_str(), (xmlChar *)value)) == NULL)
    {
        gmx_fatal(FARGS, "Creating element %s with value %s", (char *)type.c_str(), value);
    }

    return child;
}

xmlNodePtr add_xml_comment(xmlDocPtr   doc,
                           xmlNodePtr  prev,
                           const char *comment)
{
    xmlNodePtr comm, ptr;

    if ((comm = xmlNewComment((const xmlChar *)comment)) == NULL)
    {
        gmx_fatal(FARGS, "Creating doc comment element %s", "");
    }
    ptr = prev;
    while (ptr->next != NULL)
    {
        ptr = ptr->next;
    }
    ptr->next    = comm;
    comm->prev   = ptr;
    comm->doc    = doc;

    return comm;
}
