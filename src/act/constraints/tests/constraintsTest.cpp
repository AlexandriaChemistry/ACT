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
                         double             toler,
                         double             invdt)
    {
        auto               *pd = getForceField("PC+GV-elec");
        Constrainer         ccc(maxiter, toler, invdt);
        ForceComputer       forceComp;
        std::vector<ACTMol> mps;
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        std::string dirmolname("../../alexandria/tests/mols/");
        dirmolname += molname;
        auto dataName = gmx::test::TestFileManager::getInputFilePath(dirmolname);
        initACTMol(&msghandler, dataName, pd, &forceComp, &mps);
        gmx::RVec           vZero = { 0, 0, 0 };
        int isign                 = 0;
        std::vector<double> signs = { 0.3, 0.7, -0.6, -0.9, 0.5 };
        for(auto &mol: mps)
        {
            checker_.checkString(mol.getMolname(), "molname");
            std::vector<gmx::RVec> coords = mol.xOriginal();
            std::vector<gmx::RVec> velocities(coords.size(), vZero);
            auto atoms = mol.topology()->atoms();
            for(size_t i = 0; i < atoms.size(); i++)
            {
                if (atoms[i].pType() == ActParticle::Atom)
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        coords[i][m]    += 0.002*signs[isign];
                        velocities[i][m] = 0.1*signs[isign];
                        isign            = (isign+1) % signs.size();
                    }
                }
            }
            checker_.checkSequence(coords.begin(), coords.end(), "coordinates before");
            auto pot = pd->findForcesConst(InteractionType::BONDS).potential();
            int error = ccc.shake(&msghandler,
                                  pot,
                                  mol.topology()->entry(InteractionType::BONDS),
                                  atoms,
                                  &coords);
            EXPECT_TRUE(error == 0);
            if (error == 0)
            {
                checker_.checkSequence(coords.begin(), coords.end(), "coordinates after");
            }
            checker_.checkSequence(velocities.begin(), velocities.end(), "velocities before");
            error = ccc.rattle(&msghandler,
                               pot,
                               mol.topology()->entry(InteractionType::BONDS),
                               mol.topology()->atoms(),
                               coords, &velocities);
            EXPECT_TRUE(error == 0);
            if (error == 0)
            {
                checker_.checkSequence(velocities.begin(), velocities.end(), "velocities after");
            }
        }
    }
    
    static void TearDownTestCase()
    {
    }
};

TEST_F (ConstraintsTest, Water) 
{
    testConstrainer("water-3-oep.log.pdb",100, 1e-4, 1000);
}

TEST_F (ConstraintsTest, Acetate) 
{
    testConstrainer("acetate.sdf", 1000, 1e-4, 1000);
}

}

} // namespace
