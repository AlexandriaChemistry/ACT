/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018
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
#include <math.h>

#include <map>

#include <gtest/gtest.h>

#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_low.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "gromacs/topology/atoms.h"

namespace alexandria
{

namespace
{

class ForceFieldTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
        static std::vector<std::string> atomNames;
        static std::string              atomName;

        ForceFieldTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
            ForceField *mypd     = getForceField("ACM-g");
            for (const auto &iter : mypd->particleTypesConst())
            {
                atomNames.push_back(iter.second.id().id());
            }
        }

        static void TearDownTestCase()
        {
        }

        void addAtype(alexandria::ForceField *pd)
        {
            std::string          atype("U");
            std::string          elem("U");
            std::string          desc("temporary test atom");
            Identifier   atpId("U");
            ParticleType atp(atpId, desc, eptAtom);
            atp.setOption("poltype", "p_U");
            atp.setOption("zetatype", "z_U");
            atp.setOption("bondtype", "b_U");
            atp.setOption("acmtype", "z_U");
            atp.setOption("element", elem);
            atp.setOption("atomnumber", "92");
            ForceFieldParameter mm("Da", 238.29, 0.0,  1, 230, 240, Mutability::Free, true, false);
            atp.addForceFieldParameter("mass", mm);
            ForceFieldParameter rr("kJ/mol", 1000.0, 0.0,  1, 990, 1010, Mutability::Free, true, true);
            atp.addForceFieldParameter("ref_enthalpy", rr);
            pd->addParticleType(atp);
        }
};

std::vector<std::string> ForceFieldTest::atomNames;
std::string              ForceFieldTest::atomName;

TEST_F (ForceFieldTest, getAtype){
    auto aType = getForceField("ACM-g")->findParticleType("h");

    checker_.checkString(aType->element(), "element");
    checker_.checkString(aType->description(), "description");
    checker_.checkString(aType->id().id(), "type");
    auto itp = InteractionType::POLARIZATION;
    if (aType->hasInteractionType(itp))
    {
        checker_.checkString(aType->interactionTypeToIdentifier(itp).id().c_str(), "poltype");
    }
    checker_.checkString(aType->interactionTypeToIdentifier(InteractionType::BONDS).id().c_str(), "bondtype");
    checker_.checkString(aType->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION).id().c_str(), "zetatype");
}

TEST_F(ForceFieldTest, addingParticleTypeTwiceThrows)
{
    alexandria::ForceField *pd = getForceField("ACM-g");
    addAtype(pd);
    EXPECT_THROW(addAtype(pd), gmx::InvalidInputError);
}

TEST_F(ForceFieldTest, addAtype){
    alexandria::ForceField *pd = getForceField("ACM-g");
    std::string          atype("U");
    // For some reason this test "inherits" the same ForceField as the previous test
    //addAtype(pd);

    if (pd->hasParticleType(atype))
    {
        auto fa = pd->findParticleType(atype);
        /* Test if the extractions where correct */
        checker_.checkString(fa->element().c_str(), "element");
        checker_.checkString(fa->description().c_str(), "description");
        checker_.checkString(fa->interactionTypeToIdentifier(InteractionType::POLARIZATION).id(), "poltype");
        checker_.checkString(fa->interactionTypeToIdentifier(InteractionType::BONDS).id(), "bondtype");
        checker_.checkString(fa->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION).id(), "zetatype");
    }
}

TEST_F (ForceFieldTest, chi)
{
    std::vector<double>      chi0s;
    std::string              atoms[]  = { "h", "f", "p2", "br" };
    std::vector<std::string> eqd = { "ACM-g", "ACM-pg" };

    for (auto &atom : atoms)
    {
        for (auto model : eqd)
        {
            auto pd       = getForceField(model);
            auto fa       = pd->findParticleType(atom)->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            auto eep      = pd->findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            auto p        = eep.findParameterTypeConst(fa, "chi");
            chi0s.push_back(p.value());
        }
    }
    checker_.checkSequence(chi0s.begin(), chi0s.end(), "chi");
}

TEST_F (ForceFieldTest, alpha)
{
    std::vector<double>      alphas;
    std::string              atoms[]  = { "h", "c2", "f", "p2", "br" };
    std::vector<std::string> eqd = { "ACM-pg" };

    for (auto &atom : atoms)
    {
        for (auto model : eqd)
        {
            auto pd       = getForceField(model);
            auto fa       = pd->findParticleType(atom)->interactionTypeToIdentifier(InteractionType::POLARIZATION);
            auto eep      = pd->findForcesConst(InteractionType::POLARIZATION);
            auto p        = eep.findParameterTypeConst(fa, "alpha");
            alphas.push_back(p.value());
        }
    }
    checker_.checkSequence(alphas.begin(), alphas.end(), "alpha");
}

TEST_F (ForceFieldTest, row)
{
    std::vector<int> rows;
    std::string      atoms[]  = { "h", "n2", "s2", "br" };
    std::string      models[] = { "ACM-g", "ACM-pg" };

    for (auto &atom : atoms)
    {
        for (auto &model : models)
        {
            auto pd  = getForceField(model);
            auto fa  = pd->findParticleType(atom);
            auto p   = fa->row();
            rows.push_back(p);
        }
    }
    checker_.checkSequence(rows.begin(), rows.end(), "row");
}

TEST_F (ForceFieldTest, zeta)
{
    std::vector<double>      zetas;
    std::string              atoms[]  = { "h", "c3", "s2", "br" };
    std::vector<std::string> eqd = { "ACM-pg" };

    for (auto &atom : atoms)
    {
        for (auto &model : eqd)
        {
            auto pd  = getForceField(model);
            auto fa  = pd->findParticleType(atom);
            auto ztp = fa->interactionTypeToIdentifier(InteractionType::COULOMB);
            auto eep = pd->findForcesConst(InteractionType::COULOMB);
            auto p   = eep.findParameterTypeConst(ztp, "zeta");
            zetas.push_back(p.value());
        }
    }
    checker_.checkSequence(zetas.begin(), zetas.end(), "zeta");
}

}

} // namespace
