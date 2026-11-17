/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#include <algorithm>
#include <cstdlib>
#include <string>

#include "categories.h"



namespace alexandria
{

/*! \brief Compare two CategoryListElements for sorting
 * \param[in] ca First
 * \param[in] cb Second
 * \return whether first is less than second
 */
static bool CompareCategoryListElements(CategoryListElement ca,
                                        CategoryListElement cb)
{
    return (ca.getName().compare(cb.getName()) < 0);
}

/*! \brief Compare two strings for sorting
 * \param[in] ca First
 * \param[in] cb Second
 * \return whether first is less than second
 */
static bool CompareStrings(std::string ca, std::string cb)
{
    return (ca.compare(cb) < 0);
}

void CategoryListElement::addMolecule(const std::string &molecule)
{
    for (const auto &i : molecule_)
    {
        if (i == molecule)
        {
            return;
        }
    }
    molecule_.push_back(molecule);
}

void CategoryListElement::sortMolecules()
{
    std::sort(molecule_.begin(), molecule_.end(), CompareStrings);
}

bool CategoryListElement::hasMolecule(const std::string &molecule) const
{
    return std::find(molecule_.begin(), molecule_.end(), molecule) != molecule_.end();
}

void CategoryList::addCategory(const std::string &catname,
                               const std::string &molecule)
{
    bool foundCategory = false;
    for (auto &i : catli_)
    {
        if (i.getName().compare(catname) == 0)
        {
            i.addMolecule(molecule);
            foundCategory = true;
            break;
        }
    }
    if (!foundCategory)
    {
        catli_.emplace_back(catname, molecule);
    }
}

void CategoryList::sortCategories()
{
    std::sort(catli_.begin(), catli_.end(), CompareCategoryListElements);
    for (auto &i : catli_)
    {
        i.sortMolecules();
    }
}

}
