/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2026
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
 * Implements test of the penalizer algorithgm
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cmath>

#include <map>

#include <gtest/gtest.h>

#include "act/alexandria/acminitializer.h"
#include "act/ga/gene_pool.h"
#include "act/ga/penalizer.h"
#include "ga_test_helper.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"


namespace ga
{

namespace
{

class PenalizerTest : public gmx::test::CommandLineTestBase
{
protected:
    
 
    PenalizerTest ()
    {
    }

    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }

    static void TearDownTestCase()
    {
    }

    void checkGenePool(gmx::test::TestReferenceChecker *checker_,
                       const GenePool                  *pool,
                       const std::string                label,
                       double                           volume)
    {
        auto checker = checker_->checkCompound("GenePool", label);
        int i = 0;
        for(const auto &gene : pool->genePool())
        {
            auto gcheck = checker.checkCompound("Gene", std::to_string(i));
            gcheck.checkSequence(gene.bases().begin(),
                                 gene.bases().end(),
                                 label.c_str());
            std::string flabel("Fitness_");
            flabel += std::to_string(i);
            gcheck.checkFloat(gene.fitness(iMolSelect::Train), flabel.c_str());
            i += 1;
        }
        checker.checkFloat(volume, "Volume");
    }

    void vfpTest(bool logVolume, double volFracLimit, double popFrac)
    {
        // Checker
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-4);
        checker.setDefaultTolerance(tolerance);
        // Helper structure
        alexandria::GaTestHelper gth(1, {"sigma"}, {alexandria::eRMS::Interaction});
        // Make some genes. The values have to be between 0.1 and 0.5 to match the FF.
        std::vector<std::vector<double>> genes(4);
        double totalVolume = 1;
        for(size_t i = 0; i < gth.sii->nParam(); i++)
        {
            double lb = gth.sii->lowerBoundAtIndex(i);
            double ub = gth.sii->upperBoundAtIndex(i);
            totalVolume *= (ub-lb);
            for(size_t j = 0; j < genes.size(); j++)
            {
                genes[j].push_back(lb + 0.1*(i+j+1)*(ub-lb)/gth.sii->nParam());
            }
        };
        checker.checkFloat(totalVolume, "totalVolume");
        GenePool gp(genes.size());
        for(size_t i = 0; i < genes.size(); i++)
        {
            gp.addGenome(genes[i], i);
        }
        VolumeFractionPenalizer vfp(nullptr, logVolume, totalVolume,
                                    volFracLimit, popFrac,
                                    static_cast<Initializer *>(gth.initializer));
        double oldVolume = vfp.getPoolVolume(gp);
        if (logVolume)
        {
            oldVolume = std::exp(oldVolume);
        }
        checkGenePool(&checker, &gp, "Before", oldVolume);
        bool penalized = vfp.penalize(nullptr, &gp, 1);
        EXPECT_TRUE(penalized);
        double newVolume =  vfp.getPoolVolume(gp);
        if (logVolume)
        {
            newVolume = std::exp(newVolume);
        }
        checkGenePool(&checker, &gp, "After", newVolume);
        EXPECT_TRUE(newVolume > oldVolume);
    }
};

TEST_F (PenalizerTest, Volume1)
{
    bool   logVolume    = false;
    double volFracLimit = 0.2;
    double popFrac      = 0.25;
    vfpTest(logVolume, volFracLimit, popFrac);
}

TEST_F (PenalizerTest, LogVolume1)
{
    bool   logVolume    = true;
    double volFracLimit = 0.2;
    double popFrac      = 0.25;
    vfpTest(logVolume, volFracLimit, popFrac);
}

TEST_F (PenalizerTest, Volume2)
{
    bool   logVolume    = false;
    double volFracLimit = 0.5;
    double popFrac      = 0.4;
    vfpTest(logVolume, volFracLimit, popFrac);
}

TEST_F (PenalizerTest, LogVolume2)
{
    bool   logVolume    = true;
    double volFracLimit = 0.5;
    double popFrac      = 0.4;
    vfpTest(logVolume, volFracLimit, popFrac);
}

}

} 
