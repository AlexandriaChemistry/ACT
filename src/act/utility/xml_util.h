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
 
#ifndef XML_UTIL_H
#define XML_UTIL_H

#include <stdio.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

/*! Add integer to xmlNode
 * \param[in] ptr  The XML structure
 * \paran[in] name The name of the variable
 * \param[in] val  The value of the variable
 * Will crash with fatal error if operation fails
 */
void add_xml_int(xmlNodePtr ptr, const std::string &name, int val);

/*! Add double to xmlNode
 * \param[in] ptr  The XML structure
 * \paran[in] name The name of the variable
 * \param[in] val  The value of the variable
 * Will crash with fatal error if operation fails
 */
void add_xml_double(xmlNodePtr ptr, const std::string &name, double val);

/*! Add integer to xmlNode
 * \param[in] ptr  The XML structure
 * \paran[in] name The name of the variable
 * \param[in] val  The value of the variable
 * Will crash with fatal error if operation fails
 */
void add_xml_char(xmlNodePtr ptr, const std::string &name, const char *val);

/*! Add std::string to xmlNode
 * \param[in] ptr  The XML structure
 * \paran[in] name The name of the variable
 * \param[in] val  The value of the variable
 * Will crash with fatal error if operation fails
 */
void add_xml_string(xmlNodePtr ptr, const std::string &name, const std::string &val);

/*! Add a child to the xml tree
 * \param[in] parent The parent in the tree
 * \param[in] type   The name of the child
 * \returns the child
 * Will crash with fatal error if it was not possible to add a child
 */
xmlNodePtr add_xml_child(xmlNodePtr parent, const std::string &type);

/*! Add a child with a value to the xml tree
 * \param[in] parent The parent in the tree
 * \param[in] type   The name of the child
 * \param[in] value  The string value
 * \returns the child
 * Will crash with fatal error if it was not possible to add a child
 */
xmlNodePtr add_xml_child_val(xmlNodePtr parent, const std::string &type, const char *value);

/*! Add a comment to an XML file
 * \param[in] doc     The document level pointer
 * \param[in] prev    The node where to add the comment to, at the bottom of the tree
 * \param[in] comment The text of the comment
 * \returns the comment pointer
 * Will crash with fatal error if it was not possible to add a child
 */
xmlNodePtr add_xml_comment(xmlDocPtr   doc,
                           xmlNodePtr  prev,
                           const char *comment);

#endif
