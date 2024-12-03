/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023
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

#ifndef ACT_ALLMOLS_H
#define ACT_ALLMOLS_H

#include <map>
#include <set>
#include <string>
#include <vector>

namespace alexandria
{
    class AlexandriaMol
    {
    public:
        //! IUPAC name
        std::string iupac;
        //! Formula
        std::string formula;
        //! Charge
        int         charge;
        //! Multiplicity
        int         mult;
        //! Mass
        double      mass;
        //! CAS number
        std::string cas;
        //! ChemSpider ID
        int         csid;
        //! PubChem ID
        int         pubid;
        //! InChi
        std::string inchi;
        //! InChiKey
        std::string inchikey;
        //! Class identifiers
        std::set<std::string> classid;
        //! Synonyms
        std::set<std::string> synonyms;

        /*! \brief Constructor
         * \param[in] line Read from csv file and interpreted.
         */
        AlexandriaMol(const std::vector<std::string> &line);
    };
    
    class AlexandriaMols
    {
    private:
        //! The actual data
        std::map<std::string, AlexandriaMol> mols_;
        //! Translation map from synonym to InChi
        std::map<std::string, std::string>   nameToInChi_;
    public:
        //! Constructor, reads input file and stores it
        AlexandriaMols();
      
        /*! \brief Look up an AlexandriaMol
         * \param[in] inchi The Standard InChi of the compound
         * \return pointer to AlexandriaMol or nullptr if not found
         */
        const AlexandriaMol *findInChi(const std::string &inchi) const;

        /*! \brief Look up an AlexandriaMol
         * \param[in] name A molecule name
         * \return pointer to AlexandriaMol or nullptr if not found
         */
        const AlexandriaMol *findMol(const std::string &name) const;
    };

} // namespace

#endif
