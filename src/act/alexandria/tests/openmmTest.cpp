/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2025
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

#include "actmol_util.h"

#include "act/alexandria/openmm_xml.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "gromacs/mdtypes/inputrec.h"

#include "testutils/cmdlinetest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace alexandria
{

namespace
{

class OpenMMXmlTest : public gmx::test::CommandLineTestBase
{
    protected:
    void testOpenMM(const std::string &forceField,
                    const std::string &fileName,
                    bool               numberAtypes,
                    double             mDrude = 0.1)
    {
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        auto ifm = gmx::test::TextFileMatch(gmx::test::ExactTextMatch()).createFileMatcher();
        gmx::test::TestFileManager tfm;
        MsgHandler                 msghandler; 
        msghandler.setPrintLevel(ACTStatus::Warning);
        auto ff         = getForceField(forceField);
        double rmsToler = 0.0000001;
        ForceComputer fcomp(rmsToler, 25);
        std::vector<ACTMol> mps;
        std::string molfn("mols/");
        molfn += fileName;
        initACTMol(molfn.c_str(), ff, &fcomp, &mps);

        auto tmpFile  = tfm.getTemporaryFilePath("xml");
        writeOpenMM(&msghandler, tmpFile, ff, mps, mDrude, numberAtypes);
        ifm->checkFile(tmpFile, &checker);
    }
};

TEST_F(OpenMMXmlTest, Sulfate)
{
    testOpenMM("ACS-g", "sulfate.sdf", true);
}

TEST_F(OpenMMXmlTest, Furan)
{
    testOpenMM("ACS-g", "furan.sdf", true);
}

TEST_F(OpenMMXmlTest, Acetone)
{
    testOpenMM("ACS-g", "acetone.sdf", true);
}

TEST_F(OpenMMXmlTest, AcetonePolMDrude)
{
    testOpenMM("ACS-pg", "acetone.sdf", true, 0.2);
}

TEST_F(OpenMMXmlTest, SulfatePol)
{
    testOpenMM("ACS-pg", "sulfate.sdf", true);
}

TEST_F(OpenMMXmlTest, FuranPol)
{
    testOpenMM("ACS-pg", "furan.sdf", true);
}

TEST_F(OpenMMXmlTest, FuranPolNoNum)
{
    testOpenMM("ACS-pg", "furan.sdf", false);
}

TEST_F(OpenMMXmlTest, AcetonePol)
{
    testOpenMM("ACS-pg", "acetone.sdf", true);
}

TEST_F(OpenMMXmlTest, Merged)
{
    testOpenMM("ACS-g", "merged.pdb", true);
}



}

}
