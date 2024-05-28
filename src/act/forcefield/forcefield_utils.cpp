/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018
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
 * Implements test utilities for the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "forcefield_utils.h"

#include <map>

#include "gromacs/utility/exceptions.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/symcharges.h"
#include "act/forcefield/forcefield_xml.h"

#include "testutils/testfilemanager.h"

static std::map<std::string, alexandria::ForceField> pdTest;

alexandria::ForceField *getForceField(const std::string &qdist)
{
    if (pdTest.count(qdist) == 0)
    {
        std::string baseName = gmx::formatString("%s.xml", qdist.c_str());
        std::string dataName = gmx::test::TestFileManager::getInputFilePath(baseName);
        try
        {
            alexandria::ForceField pd;
            alexandria::readForceField(dataName, &pd);
            pdTest.insert(std::pair<std::string, alexandria::ForceField>(qdist, std::move(pd)));
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    return &(pdTest[qdist]);
}

