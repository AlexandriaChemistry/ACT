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
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forces/forcecomputer.h"
#include "act/forces/forcecomputerutils.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_util.h"
#include "act/topology/topology.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/exceptions.h"
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
                         ConstraintType     ctype,
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
            ccc.makeConstraints(&msghandler, ctype, pd, mol.topology(), coords);
            if (ConstraintType::Bonds != ctype)
            {
                checker_.checkInt64(ccc.numBondConstraints(), "bond constraints");
                checker_.checkInt64(ccc.numAngleConstraints(), "angle constraints");
                checker_.checkInt64(ccc.numRedundantAngles(), "redundant angles");
            }
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
            int error = ccc.shake(&msghandler, &coords);
            EXPECT_TRUE(error == 0);
            if (error == 0)
            {
                checker_.checkSequence(coords.begin(), coords.end(), "coordinates after");
            }
            // All distances should be at their constrained values
            auto &pairs = ccc.pairs();
            for(size_t ll = 0; ll < pairs.size(); ll++)
            {
                rvec dx;
                rvec_sub(coords[pairs[ll][0]], coords[pairs[ll][1]], dx);
                EXPECT_NEAR(ccc.length(ll), norm(dx), 2*toler*ccc.length(ll)) << "constraint " << ll;
            }
            if (ConstraintType::Angles == ctype && ccc.numRedundantAngles() == 0)
            {
                // Without redundant constraints all angles are at their reference values
                for(const auto &a : mol.topology()->entry(InteractionType::ANGLES))
                {
                    auto &ind = a->atomIndices();
                    rvec  r_ij, r_kj;
                    real  costh;
                    auto  theta = bond_angle(coords[ind[0]], coords[ind[1]], coords[ind[2]],
                                             r_ij, r_kj, &costh);
                    EXPECT_NEAR(a->params()[angleANGLE], theta, 1e-3);
                }
            }
            checker_.checkSequence(velocities.begin(), velocities.end(), "velocities before");
            error = ccc.rattle(&msghandler, coords, &velocities);
            EXPECT_TRUE(error == 0);
            if (error == 0)
            {
                checker_.checkSequence(velocities.begin(), velocities.end(), "velocities after");
            }
            // No relative velocity along the constraints, within the
            // tolerance of rattle, which is relative to the change of
            // the distance in one time step.
            for(size_t ll = 0; ll < pairs.size(); ll++)
            {
                rvec dx, dv;
                rvec_sub(coords[pairs[ll][0]], coords[pairs[ll][1]], dx);
                rvec_sub(velocities[pairs[ll][0]], velocities[pairs[ll][1]], dv);
                EXPECT_NEAR(0, iprod(dx, dv)/norm(dx), 2*toler*ccc.length(ll)*invdt) << "constraint " << ll;
            }
        }
    }
    
    static void TearDownTestCase()
    {
    }
};

TEST_F (ConstraintsTest, Water) 
{
    testConstrainer("water-3-oep.log.pdb", ConstraintType::Bonds, 100, 1e-4, 1000);
}

TEST_F (ConstraintsTest, Acetate) 
{
    testConstrainer("acetate.sdf", ConstraintType::Bonds, 1000, 1e-4, 1000);
}

TEST_F (ConstraintsTest, WaterAngles) 
{
    testConstrainer("water-3-oep.log.pdb", ConstraintType::Angles, 100, 1e-4, 1000);
}

TEST_F (ConstraintsTest, AcetateAngles) 
{
    testConstrainer("acetate.sdf", ConstraintType::Angles, 1000, 1e-4, 1000);
}

TEST_F (ConstraintsTest, MethanethiolAngles) 
{
    testConstrainer("methanethiol.sdf", ConstraintType::Angles, 1000, 1e-4, 1000);
}

TEST_F (ConstraintsTest, NoConstraints) 
{
    testConstrainer("acetate.sdf", ConstraintType::No, 1000, 1e-4, 1000);
}

//! \brief Load a molecule with a force field for the constraint tests
static void loadMolecule(ForceField          *pd,
                         const std::string   &molname,
                         std::vector<ACTMol> *mps)
{
    ForceComputer forceComp;
    MsgHandler    msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    auto dataName = gmx::test::TestFileManager::getInputFilePath("../../alexandria/tests/mols/" + molname);
    initACTMol(&msghandler, dataName, pd, &forceComp, mps);
}

TEST(ConstrainerTest, ProjectForces)
{
    auto               *pd  = getForceField("PC+GV-elec");
    std::vector<ACTMol> mps;
    loadMolecule(pd, "acetate.sdf", &mps);
    ASSERT_EQ(1u, mps.size());
    auto               &mol = mps[0];
    Constrainer         ccc(100, 1e-8, 1);
    auto        coords = mol.xOriginal();
    ccc.makeConstraints(nullptr, ConstraintType::Angles, pd, mol.topology(), coords);
    ASSERT_EQ(0, ccc.shake(nullptr, &coords));
    std::vector<gmx::RVec> forces(coords.size());
    for(size_t i = 0; i < forces.size(); i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            forces[i][m] = std::sin(1.0 + 3*i + m);
        }
    }
    ASSERT_EQ(0, ccc.projectForces(nullptr, coords, &forces));
    // The accelerations do not change the constrained distances
    auto &atoms = mol.topology()->atoms();
    for(size_t c = 0; c < ccc.numConstraints(); c++)
    {
        auto  &p = ccc.pairs()[c];
        rvec   dx, ai, aj, da;
        rvec_sub(coords[p[0]], coords[p[1]], dx);
        svmul(atoms[p[0]].invmass(), forces[p[0]], ai);
        svmul(atoms[p[1]].invmass(), forces[p[1]], aj);
        rvec_sub(ai, aj, da);
        EXPECT_NEAR(0, iprod(dx, da), 1e-10) << "constraint " << c;
    }
    // Projecting again does not change anything
    auto forces2 = forces;
    ASSERT_EQ(0, ccc.projectForces(nullptr, coords, &forces2));
    for(size_t i = 0; i < forces.size(); i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            EXPECT_NEAR(forces[i][m], forces2[i][m], 1e-10);
        }
    }
}

TEST(ConstrainerTest, Frozen)
{
    auto               *pd  = getForceField("PC+GV-elec");
    std::vector<ACTMol> mps;
    loadMolecule(pd, "acetate.sdf", &mps);
    ASSERT_EQ(1u, mps.size());
    auto               &mol = mps[0];
    Constrainer         ccc(100, 1e-8, 1);
    auto        coords = mol.xOriginal();
    ccc.makeConstraints(nullptr, ConstraintType::Bonds, pd, mol.topology(), coords);
    // Freeze the two carbons, the bond between them is not constrained anymore
    ccc.setFrozen({ 0, 1 });
    EXPECT_EQ(5u, ccc.numConstraints());
    auto xfrozen = coords;
    ASSERT_EQ(0, ccc.shake(nullptr, &coords));
    for(int i : { 0, 1 })
    {
        for(int m = 0; m < DIM; m++)
        {
            EXPECT_EQ(xfrozen[i][m], coords[i][m]);
        }
    }
    for(size_t c = 0; c < ccc.numConstraints(); c++)
    {
        auto &p = ccc.pairs()[c];
        rvec  dx;
        rvec_sub(coords[p[0]], coords[p[1]], dx);
        EXPECT_NEAR(ccc.length(c), norm(dx), 1e-8);
    }
}

TEST(ConstraintTypeTest, Strings)
{
    for(auto ct : { ConstraintType::No, ConstraintType::Bonds, ConstraintType::Angles })
    {
        EXPECT_EQ(ct, stringToConstraintType(constraintTypeToString(ct)));
    }
    EXPECT_EQ("angles", constraintTypeToString(ConstraintType::Angles));
    EXPECT_THROW(stringToConstraintType("dihedrals"), gmx::InvalidInputError);
}

}

} // namespace
