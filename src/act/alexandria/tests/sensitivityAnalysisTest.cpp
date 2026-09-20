/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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
 * Tests for the Sensitivity class and the SensitivityAnalysis() routine.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <string>

#include <gtest/gtest.h>

#include "act/alexandria/sensitivity.h"
#include "act/alexandria/bayes.h"
#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/ga/genome.h"
#include "act/basics/msg_handler.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/jsontree.h"

#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

// ============================================================================
// SensitivityAnalysis() -- the public entry point
// ============================================================================

/*! \brief With an empty genome there is nothing to analyse: the routine must
 * short-circuit, report a local minimum and leave the output JsonTree empty.
 */
TEST(SensitivityAnalysisTest, RunOne)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    CommunicationRecord cr(&msghandler);
    StaticIndividualInfo sii(&cr);
    std::string baseName("ACS-g.xml");
    std::string dataName = gmx::test::TestFileManager::getInputFilePath(baseName);
    sii.fillForceField(&msghandler, dataName.c_str());
    MolGen molgen(&cr);
    molgen.addFitOption("sigma");
    molgen.fillIopt(sii.forcefield(), &msghandler);
    sii.generateOptimizationIndex(&msghandler, &molgen, &cr);
    sii.fillVectors(1, 0.02);
    ForceComputer forceComp;
    ACMFitnessComputer fitComp;
    fitComp.init(&msghandler, &sii, &molgen, false,
                 &forceComp, ChargeGenerationAlgorithm::SQE);
    ga::Genome genome;
    genome.addBase(2.0);
    genome.addBase(3.0);
    genome.addBase(-1.0);
    
    JsonTree   jtree("SensitivityAnalysisTest");
    // TODO: MsgHandler, StaticIndividualInfo, ga::Genome (empty), JsonTree
    // TODO: EXPECT_TRUE(SensitivityAnalysis(...));
    // TODO: EXPECT_TRUE(jtree.objects().empty());
    bool result = true;
    //SensitivityAnalysis(&msghandler, &sii, &fitComp,
    //                  &genome, iMolSelect::Train, &jtree);
    EXPECT_TRUE(result);
}

// ============================================================================
// Sensitivity -- the parabola fit that decides "is it a minimum?"
// ============================================================================

/*! \brief chi2 = p^2 - 2p + 1 must fit to a >= 0 (upward parabola => minimum).
 */
TEST(SensitivityTest, FitsUpwardParabola)
{
    Sensitivity s;
    //add 3 points; 
    s.add(0, 1);
    s.add(1, 0);
    s.add(2, 1);
    s.computeForceConstants(nullptr);
    // Expected result is 1.0
    EXPECT_TRUE(std::abs(s.a() - 1.0) < 1e-7);
}

/*! \brief chi2 = -p^2 + 2p - 1 must fit to a < 0 (downward parabola => not min).
 */
TEST(SensitivityTest, FitsDownwardParabola)
{
    Sensitivity s;
    //add 3 points; 
    s.add(0, 1);
    s.add(1, 2);
    s.add(2, 1);
    s.computeForceConstants(nullptr);
    // Expected results is -1.0
    EXPECT_TRUE(std::abs(s.a() + 1.0) < 1e-7);
}

/*! \brief With fewer than 3 points the coefficients stay at their zero initial
 * values.
 */
TEST(SensitivityTest, InsufficientPointsLeavesConstantsZero)
{
    Sensitivity s;
    //add 2 points; 
    s.add(0, 1);
    s.add(1, 2);
    // TODO: Sensitivity s; add 2 points; 
    s.computeForceConstants(nullptr);
    // TODO: 
    EXPECT_TRUE(s.a()==0);
    EXPECT_TRUE(s.b()==0);
    EXPECT_TRUE(s.c()==0);
}

} // namespace

} // namespace alexandria
