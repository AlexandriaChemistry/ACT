/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2020
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include "alexandria/identifier.h"

#include <cstdio>

#include <map>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

typedef std::map<const Identifier, int> idmap;

class IdentifierMapTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;

        IdentifierMapTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        static void TearDownTestCase()
        {
        }

        void runTest(const idmap &idmap)
        {
            for (const auto &id : idmap)
            {
                auto mapfind = idmap.find(id.first);
                if (mapfind == idmap.end())
                {
                    fprintf(stderr, "Could not find %s in map. Shame on you.\n", id.first.id().c_str());
                }
                else
                {
                    checker_.checkInteger(mapfind->second, id.first.id().c_str());
                }
            }
        }

};

TEST_F(IdentifierMapTest, FindInMap) {
    Identifier PH2({"P", "H"}, 2, CanSwap::Yes);
    Identifier BP2({"B", "P"}, 2, CanSwap::Yes);
    Identifier HC1({"H", "C"}, 1, CanSwap::Yes);
    Identifier NC1({"N", "C"}, 1, CanSwap::Yes);
    Identifier NC2({"N", "C"}, 2, CanSwap::Yes);
    Identifier NO2({"N", "O"}, 2, CanSwap::No);
    Identifier OC2({"O", "C"}, 2, CanSwap::Yes);
    Identifier CO1({"O", "C"}, 1, CanSwap::Yes);
    Identifier HO1({"H", "O"}, 1, CanSwap::Yes);
    std::map<const Identifier, int> idmap = {
        { PH2, 1 },
        { BP2, 2 },
        { HC1, 3 },
        { NC2, 4 },
        { CO1, 5 },
        { HO1, 6 },
        { OC2, 7 },
        { NC1, 8 },
        { NO2, 9 }
    };
    runTest(idmap);
}


}

}
