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

#ifndef ACT_ATOMIZATION_ENERGY_H
#define ACT_ATOMIZATION_ENERGY_H

#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "msg_handler.h"

namespace alexandria
{

    class AtomizationEnergyTerm;
    
    class AtomizationEnergy
    {
    public:
        //! Constructor
        AtomizationEnergy() {}

        /*! Reads input file and stores it
         * \param[in] msg_handler Will be used for information if not nullptr
         */
        void read(MsgHandler *msg_handler = nullptr);

        //! Destructor to clean up memory
        ~AtomizationEnergy();
        
        /*! \brief Return a term according to specifications
         * \param[in] elem   The element or atom
         * \param[in] charge The charge of the element/atom
         * \param[in] source The source of the data
         * \param[in] prop   The property to extract
         * \param[in] T      The temperature (multiples of 100K or 298.15K only)
         * \param[out] unit  The unit of the value (may be nullptr)
         * \param[out] ref   The reference (may be nullptr)
         */
        double term(const std::string &elem,
                    int                charge,
                    const std::string &source,
                    const std::string &prop,
                    double             T,
                    std::string       *unit,
                    std::string       *ref) const;

        /*! \brief Dump the data to a new file
         * \param[in] filenm The name of the file to dump to
         */
        void dump(const std::string &filenm);

    private:
        //! All the terms in one array
        std::vector<AtomizationEnergyTerm *> terms_;
};

}

#endif
