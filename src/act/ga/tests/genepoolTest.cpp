/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2024
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
#include <cmath>

#include <map>

#include <gtest/gtest.h>

#include "act/ga/gene_pool.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"


namespace ga
{

namespace
{

class GenePoolTest : public gmx::test::CommandLineTestBase
{
protected:
 
    GenePoolTest ()
    {
    }

    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }

    static void TearDownTestCase()
    {
    }

    void compare(const GenePool &gp)
    {
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-8);
        checker.setDefaultTolerance(tolerance);
        int gene = 0;
        for(const auto &g : gp.genePool())
        {
            gene++;
            std::string buf = "Gene " + std::to_string(gene);
            checker.checkSequence(g.bases().begin(), g.bases().end(),
                                  buf.c_str());
        }
    }

    void write()
    {
        GenePool gp;
        std::vector<std::vector<double> > bases = {
            { 3.4, 5.6, 0.1, 2.3 },
            { 6.7, 7.8, 0.9, 2.4 },
            { 1.2, 3.4, 5.4,-2.1 }
        };
        FitnessMap fm = { { iMolSelect::Train, 0.1 } };
        for(const auto &b : bases)
        {
            Genome g(b, fm);
            gp.addGenome(g);
        }
        std::string fileName("gptest.csv");
        gp.write(fileName);
        GenePool gp2;
        gp2.read(fileName);
        compare(gp2);
    }

    void read()
    {
        GenePool gp;
        gp.read(gmx::test::TestFileManager::getInputFilePath("gp.csv"));
        compare(gp);
    }
};

TEST_F (GenePoolTest, Write)
{
    write();
}

TEST_F (GenePoolTest, Read)
{
    read();
}

}

} 
