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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cmath>

#include <map>

#include <gtest/gtest.h>

#include "act/constraints/shake.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forces/forcecomputer.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_util.h"
#include "act/topology/topology.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{


namespace
{

class ConstraintsTest : public gmx::test::CommandLineTestBase
{
protected:
    gmx::test::TestReferenceChecker checker_;
    
    ConstraintsTest () : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-6);
        checker_.setDefaultTolerance(tolerance);
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    void testConstrainer(const std::string &molname,
                         int                maxiter,
                         double             toler)
    {
        auto               *pd = getForceField("PC+GV-elec");
        Constrainer         ccc(maxiter, toler);
        ForceComputer       forceComp;
        std::vector<ACTMol> mps;
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        std::string dirmolname("../../alexandria/tests/mols/");
        dirmolname     += molname;
        initACTMol(dirmolname.c_str(), pd, &forceComp, &mps);
        gmx::RVec           fZero = { 0, 0, 0 };
        for(auto &mol: mps)
        {
            checker_.checkString(mol.getMolname(), "molname");
            std::vector<gmx::RVec> coords = mol.xOriginal();
            std::vector<gmx::RVec> forces(coords.size(), fZero);
            int error = ccc.shake(&msghandler,
                                  mol.topology()->entry(InteractionType::BONDS),
                                  mol.topology()->atoms(),
                                  &coords,
                                  &forces);
            EXPECT_TRUE(error == 0);
            if (error == 0)
            {
                checker_.checkSequence(coords.begin(), coords.end(),
                                       "coordinates");
                checker_.checkSequence(forces.begin(), forces.end(),
                                       "forces");
            }
        }
    }
    
    static void TearDownTestCase()
    {
    }
};

TEST_F (ConstraintsTest, Water) 
{
    testConstrainer("water-3-oep.log.pdb",100, 1e-4);
}

TEST_F (ConstraintsTest, Acetate) 
{
    testConstrainer("acetate.sdf", 100, 1e-4);
}

}

} // namespace
