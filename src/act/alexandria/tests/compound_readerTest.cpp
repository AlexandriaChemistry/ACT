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
    MsgHandler                       msghandler_;
    CompoundReader                   compR_;
    ForceComputer                    forceComp_;
    ForceField                      *pd;
    std::string                      qopt_;
    std::string                      fopt_;
    std::string                      charges_;
    gmx::test::TestReferenceChecker  checker_;

    //! Constructor that does initialization
    CompoundReaderTest() : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        pd = getForceField("ff-hh");
        qopt_.assign("-charges");
        fopt_.assign("-f");
        charges_ = gmx::test::TestFileManager::getInputFilePath("mp2-hh.xml");
    }

    //! Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    void hack_fnm(std::vector<t_filenm> *filenm)
    {
        for(auto &fn : *filenm)
        {
            if (qopt_.compare(fn.opt) == 0)
            {
                fn.fn    = charges_.c_str();
                fn.filenames.push_back(charges_);
                fn.flag &= ffSET;
            }
        }
    }
    void hack_pargs(std::vector<t_pargs> *pargs, std::string &infile)
    {
        for(auto &pa : *pargs)
        {
            if (fopt_.compare(pa.option) == 0 && infile.size() > 0)
            {
                *pa.u.c = infile.c_str();
                pa.bSet = true;
            }
        }
    }
    //! Do the actual testing
    void testCompoundReader ()
    {
        std::vector<ACTMol>       mols;
        std::vector<t_filenm>     filenm;
        std::vector<t_pargs>      pargs;
        std::vector<const char *> desc;
        CompoundReader            compR_;

        compR_.addOptions(&pargs, &filenm, &desc);
        hack_fnm(&filenm);
        compR_.optionsFinished(&msghandler_, filenm);
        checker_.checkString(compR_.qread(), "Charge type to read");
        checker_.checkString(chargeGenerationAlgorithmName(compR_.algorithm()).c_str(), "Algorithm");
        compR_.read(&msghandler_, *pd, &forceComp_, &mols);
        auto &qMap = compR_.chargeMapConst();
        checker_.checkInt64(qMap.size(), "Number of molecules");
        int i = 1;
        for(const auto &mol : qMap)
        {
            for(const auto &atom : mol.second)
            {
                std::string aname = gmx::formatString("%s-%d",
                                                      atom.first.id().c_str(), i);
                checker_.checkDouble(atom.second, aname.c_str());
                i += 1;
            }
        }
    }

    //! Do the actual testing
    void testCompoundReaderFile (const std::string &input, bool hackFNM)
    {
        std::vector<ACTMol>       mols;
        std::vector<t_filenm>     filenm;
        std::vector<t_pargs>      pargs;
        std::vector<const char *> desc;
        CompoundReader            compR_;

        std::string infile = gmx::test::TestFileManager::getInputFilePath(input);
        compR_.addOptions(&pargs, &filenm, &desc);
        hack_pargs(&pargs, infile);
        if (hackFNM)
        {
            hack_fnm(&filenm);
        }
        compR_.optionsFinished(&msghandler_, filenm);
        checker_.checkString(compR_.qread(), "Charge type to read");
        checker_.checkString(chargeGenerationAlgorithmName(compR_.algorithm()).c_str(), "Algorithm");
        compR_.read(&msghandler_, *pd, &forceComp_, &mols);
        checker_.checkInt64(mols.size(), "Number of molecules");
        for(const auto &mol : mols)
        {
            const auto &atoms = mol.atomsConst();
            int i = 1;
            for(const auto &atom : atoms)
            {
                std::string aname = gmx::formatString("%s-%s-%d", mol.formula().c_str(), atom.name().c_str(), i);
                checker_.checkDouble(atom.charge(), aname.c_str());
                i += 1;
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

TEST_F (CompoundReaderTest, HCl){
    testCompoundReaderFile("hydrogen-chloride.sdf", false);
}

TEST_F (CompoundReaderTest, Dimer){
    testCompoundReaderFile("hydrogen-bromide#hydrogen-fluoride.pdb", false);
}

TEST_F (CompoundReaderTest, Fluorane){
    testCompoundReaderFile("fluorane.sdf", false);
}

TEST_F (CompoundReaderTest, HCl_Charges){
    testCompoundReaderFile("hydrogen-chloride.sdf", true);
}

TEST_F (CompoundReaderTest, Dimer_Charges){
    testCompoundReaderFile("hydrogen-bromide#hydrogen-fluoride.pdb", true);
}

TEST_F (CompoundReaderTest, Fluorane_Charges){
    testCompoundReaderFile("fluorane.sdf", true);
}

} // namespace alexandria
