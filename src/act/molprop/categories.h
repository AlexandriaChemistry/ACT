/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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


#ifndef CATEGORIES_H
#define CATEGORIES_H

#include <string.h>

#include <string>
#include <vector>

#include "gromacs/utility/real.h"

#include "alexandria/molselect.h"

/*! \brief
 * Contains classes related to alexandria force field tools
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
namespace alexandria
{
class MolProp;

class CategoryListElement
{
    private:
        //! The name of the category
        std::string              cat_;
        //! Molecules in this category
        std::vector<std::string> molecule_;
    public:
        /*! \brief Constructor
         * \param[in] cat      The category
         * \param[in] molecule The molecule
         */
        CategoryListElement(const std::string &cat, const std::string &molecule)
            : cat_(cat)
        {
            addMolecule(molecule);
        }
        /*! \brief Add a molecule to the list
         * \param[in] molecule The molecule to add
         */
        void addMolecule(const std::string &molecule);

        //! \return the number of molecules
        int nMolecule() const { return molecule_.size(); }

        /*! \brief Check whether molecule is present
         * \param[in] molecule The molecule to look for
         * \return true if found, false otherwise
         */
        bool hasMolecule(const std::string &molecule) const;

        //! \return the name of the category
        const std::string &getName() const { return cat_; }

        //! \brief Sort the molecules in this category
        void sortMolecules();

        //! \return The molecule names
        const std::vector<std::string> &molecules() const { return molecule_; }
};

/*! \brief List of compound categories and the compounds associated with that
 */
class CategoryList
{
    private:
        //! Vector of CategoryListElement items
        std::vector<CategoryListElement> catli_;
    public:
        //! Constructor
        CategoryList() {};

        /*! \brief Add one category and associated molecule
         * \param[in] catname The category name
         * \param[in] molecule The molecule name
         */
        void addCategory(const std::string &catname,
                         const std::string &molecule);

        //! Sort the categories alphabetically
        void sortCategories();

        //! The number of categories
        int nCategories() { return catli_.size(); }

        //! Non-mutable array of items
        const std::vector<CategoryListElement> &elementsConst() const { return catli_; }

        //! Mutable list of items
         std::vector<CategoryListElement> *elements() { return &catli_; }
};

/*! \brief Make the entire category list
 * \param[out] cList The final list
 * \param[in]  mp    The molecule properties
 * \param[in]  gms   The selection structure
 * \param[in]  ims   The data set selected
 */
void makeCategoryList(CategoryList               *cList,
                      const std::vector<MolProp> &mp,
                      const MolSelect            &gms,
                      iMolSelect                  ims);

}

#endif
