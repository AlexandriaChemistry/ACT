/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2025
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

#include "actpre.h"

#include <cmath>

#include <cstdlib>

#include <gtest/gtest.h>

#include "act/alexandria/compound_reader.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_utils.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

class CompoundReaderTest : public gmx::test::CommandLineTestBase
{
protected:
    MsgHandler             msghandler_;
    std::vector<Bond>      bonds_;
    std::vector<gmx::RVec> x_, y_;
    gmx::test::TestReferenceChecker checker_;

    //! Constructor that does initialization
    CompoundReaderTest() : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
    }

    //! Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    //! Do the actual testing
    void testCompoundReader ()
    {
        CompoundReader            compR;
        ForceComputer             forceComp;
        std::vector<ACTMol>       mols;
        std::vector<t_filenm>     filenm;
        std::vector<t_pargs>      pargs;
        std::vector<const char *> desc;

        auto pd      = getForceField("ff-hh");
        auto charges = gmx::test::TestFileManager::getInputFilePath("mp2-hh.xml");
        compR.addOptions(&pargs, &filenm, &desc);
        std::string qstr("-charges");
        for(auto &fn : filenm)
        {
            if (qstr.compare(fn.opt) == 0)
            {
                fn.fn    = charges.c_str();
                fn.filenames.push_back(charges);
                fn.flag &= ffSET;
            }
        }
        compR.optionsFinished(&msghandler_, filenm);
        checker_.checkString(compR.qread(), "Charge type to read");
        checker_.checkString(chargeGenerationAlgorithmName(compR.algorithm()).c_str(), "Algorithm");
        compR.read(&msghandler_, *pd, &forceComp, &mols);
        checker_.checkInt64(mols.size(), "Number of molecules");
        for(const auto &mol : mols)
        {
            const auto &atoms = mol.atomsConst();
            int i = 1;
            for(const auto &atom : atoms)
            {
                std::string aname = gmx::formatString("%s-%s-%d", mol.formula().c_str(), atom.name().c_str(), i);
                checker_.checkDouble(atom.charge(), aname.c_str());
            }
        }
    }

    //! Clean the test data.
    static void TearDownTestCase()
    {
    }
    
};

TEST_F (CompoundReaderTest, Count){
    testCompoundReader();
}

} // namespace alexandria
