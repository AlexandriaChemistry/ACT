/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

#include "jsontree.h"

#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"
    
namespace alexandria
{

static void add_space(std::string *str, int n)
{
    for(int i = 0; i < n; i++)
    {
        *str += " ";
    }
}

const std::string JsonTree::writeString(bool json, int *indent) const
{
    std::string str;
    add_space(&str, *indent);
    if (json)
    {
        str += gmx::formatString("{\"%s\":", key_.c_str());
        *indent += 2;
        if (!value_.empty())
        {
            str += gmx::formatString("\"%s\"", value_.c_str());
        }
        else
        {
            str += "[\n";
            size_t count = 0;
            for (const auto &obj : objects_)
            {
                str += obj.writeString(json, indent);
                if (count < objects_.size()-1)
                {
                    str += ",";
                }
                str += "\n";
                count++;
            }
            add_space(&str, *indent);
            str += "]";
        }
        *indent -= 2;
        str += "}";
    }
    else
    {
        str += gmx::formatString("%s:", key_.c_str());
        *indent += 2;
        if (!value_.empty())
        {
            str += gmx::formatString(" %s", value_.c_str());
        }
        else
        {
            str += "\n";
            size_t count = 0;
            for (const auto &obj : objects_)
            {
                str += obj.writeString(json, indent);
                if (count < objects_.size()-1)
                {
                    str += ",";
                }
                str += "\n";
                count++;
            }
            add_space(&str, *indent);
        }
        *indent -= 2;
    }
    return str;
}
        
void JsonTree::fwrite(FILE *fp,
                      bool  json)
{
    int indent = 0;
    fprintf(fp, "%s\n", writeString(json, &indent).c_str());
}
        
void JsonTree::write(const std::string &fileName,
                     bool               json)
{
    FILE *fp = gmx_ffopen(fileName.c_str(), "w");
    fwrite(fp, json);
    gmx_ffclose(fp);
}
        
} // namespace alexandria
