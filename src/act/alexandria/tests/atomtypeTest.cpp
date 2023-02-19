/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2019-2023
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cmath>

#include <map>

#include <gtest/gtest.h>

#include "act/alexandria/babel_io.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/actmol.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/qgen/qgen_acm.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class AtomtypeTest : public gmx::test::CommandLineTestBase
{
private:
    ForceField *pd;

public:   
    AtomtypeTest()
    {
        pd = getForceField("ACS-g");
    }
    
protected:
    void testAtype(const char *molname, bool bondTest=false)
        {
            gmx::test::TestReferenceChecker checker_(this->rootChecker());
            int                             maxpot   = 100;
            int                             nsymm    = 0;
            const char                     *conf     = (char *)"minimum";
            std::string                     basis, method;
            const char                     *jobtype  = (char *)"Opt";

            std::string                     dataName;
            std::vector<alexandria::MolProp> molprops;

            dataName    = gmx::test::TestFileManager::getInputFilePath(molname);
            double qtot = 0.0;
            bool readOK = readBabel(pd, dataName.c_str(), &molprops, molname, molname,
                                    conf, &method, &basis, maxpot, nsymm, jobtype, &qtot, false, true);
            EXPECT_TRUE(readOK);
            for(auto &molprop: molprops)
            {
                std::vector<std::string> atypes;
                auto                     exper = molprop.experimentConst().begin();
                for (auto &ca : exper->calcAtomConst())
                {
                    atypes.push_back(ca.getObtype());
                }
                if (!bondTest)
                {
                    checker_.checkInteger(static_cast<int>(atypes.size()), molname);
                    checker_.checkSequence(atypes.begin(), atypes.end(), "atomtypes");
                }
                else
                {
                    std::vector<std::string> bondorder;
                    for (auto &bond : molprop.bondsConst())
                    {
                        auto ai = bond.aI();
                        auto aj = bond.aJ();
                        bondorder.push_back(gmx::formatString("%s-%d %s-%d: %g",
                                                              atypes[ai].c_str(), ai,
                                                              atypes[aj].c_str(), aj,
                                                              bond.bondOrder()));
                    }
                    checker_.checkInteger(static_cast<int>(bondorder.size()), molname);
                    checker_.checkSequence(bondorder.begin(), bondorder.end(), "bondorder");
                }
            }
        }
};

class BondtypeTest : public AtomtypeTest
{
    protected:
    void testAtype(const char *molname)
    {
        AtomtypeTest::testAtype(molname, true);
    }
};

#include "atom_bond_test_include.cpp"

}

}
