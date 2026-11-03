/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2026
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

#include <fstream>
#include <functional>
#include <map>
#include <sstream>

#include <gtest/gtest.h>

#include "act/forcefield/forcefield.h"
#include "act/forcefield/symcharges.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

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
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
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
            ParticleType atp(atpId, desc, ActParticle::Atom);
            atp.setOption("poltype", "p_U");
            atp.setOption("zetatype", "z_U");
            atp.setOption("bondtype", "b_U");
            atp.setOption("acmtype", "z_U");
            atp.setOption("element", elem);
            atp.setOption("atomnumber", "92");
            ForceFieldParameter mm("Da", 238.29, 0.0,  1, 230, 240, Mutability::Bounded, true, false);
            atp.addForceFieldParameter("mass", mm);
            ForceFieldParameter rr("kJ/mol", 1000.0, 0.0,  1, 990, 1010, Mutability::Bounded, true, true);
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
            // The width of the charge distribution of the particle type
            auto fa  = pd->findParticleType(atom);
            zetas.push_back(fa->parameterConst(zetaName).value());
        }
    }
    checker_.checkSequence(zetas.begin(), zetas.end(), "zeta");
}

// Test: A force field with a mutability that is not supported, such as
// the former Free, cannot be read
TEST(ForceFieldSimpleTest, unknownMutabilityThrows)
{
    std::string dataName = gmx::test::TestFileManager::getInputFilePath("ACM-g.xml");
    std::ifstream in(dataName);
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string xml = buffer.str();
    auto pos = xml.find("mutability=\"Bounded\"");
    ASSERT_NE(std::string::npos, pos);
    xml.replace(pos, std::string("mutability=\"Bounded\"").size(), "mutability=\"Free\"");
    gmx::test::TestFileManager fileManager;
    std::string freeName = fileManager.getTemporaryFilePath("ACM-g-free.xml");
    {
        std::ofstream out(freeName);
        out << xml;
    }
    ForceField pd;
    EXPECT_THROW(readForceField(freeName, &pd), gmx::InvalidInputError);
}

//! \return the file name of ACM-g.xml with a change made by edit
static std::string modifiedForceField(gmx::test::TestFileManager       *fileManager,
                                      const std::string                &name,
                                      std::function<void(std::string *)> edit,
                                      const std::string                &base = "ACM-g.xml")
{
    std::ifstream     in(gmx::test::TestFileManager::getInputFilePath(base));
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string xml = buffer.str();
    edit(&xml);
    auto fileName = fileManager->getTemporaryFilePath(name);
    std::ofstream out(fileName);
    out << xml;
    return fileName;
}

// Test: A particle type of which the ACM type has no parameters cannot be
// read, unless without the check of consistency
TEST(ForceFieldSimpleTest, missingAcmParametersThrows)
{
    gmx::test::TestFileManager fileManager;
    auto fileName = modifiedForceField(&fileManager, "no-iodine.xml", [](std::string *xml)
    {
        auto eem   = xml->find("interaction type=\"ELECTRONEGATIVITYEQUALIZATION\"");
        auto begin = xml->find("    <parameterlist identifier=\"i_z\">", eem);
        auto end   = xml->find("</parameterlist>", begin);
        ASSERT_NE(std::string::npos, begin);
        xml->erase(begin, end + std::string("</parameterlist>\n").size() - begin);
    });
    ForceField pd;
    EXPECT_THROW(readForceField(fileName, &pd), gmx::InvalidInputError);
    ForceField incomplete;
    readForceField(fileName, &incomplete, nullptr, true, false);
    EXPECT_FALSE(incomplete.findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION).parameterExists(Identifier("i_z")));
}

// Test: A particle type with a charge from ACM but without ACM type cannot be read
TEST(ForceFieldSimpleTest, missingAcmTypeThrows)
{
    gmx::test::TestFileManager fileManager;
    auto fileName = modifiedForceField(&fileManager, "no-acmtype.xml", [](std::string *xml)
    {
        auto        ptype  = xml->find("<particletype identifier=\"i\"");
        std::string option("      <option key=\"acmtype\" value=\"i_z\"/>\n");
        auto        pos    = xml->find(option, ptype);
        ASSERT_NE(std::string::npos, ptype);
        ASSERT_NE(std::string::npos, pos);
        xml->erase(pos, option.size());
    });
    ForceField pd;
    EXPECT_THROW(readForceField(fileName, &pd), gmx::InvalidInputError);
    ForceField incomplete;
    readForceField(fileName, &incomplete, nullptr, true, false);
    EXPECT_FALSE(incomplete.findParticleType("i")->hasInteractionType(InteractionType::ELECTRONEGATIVITYEQUALIZATION));
}

// Test: Point charges with a width cannot be read
TEST(ForceFieldSimpleTest, pointChargesWithZetaThrow)
{
    gmx::test::TestFileManager fileManager;
    auto fileName = modifiedForceField(&fileManager, "point.xml", [](std::string *xml)
    {
        std::string gaussian("<option key=\"chargedistribution\" value=\"Gaussian\"/>");
        auto        pos = xml->find(gaussian);
        ASSERT_NE(std::string::npos, pos);
        xml->replace(pos, gaussian.size(), "<option key=\"chargedistribution\" value=\"Point\"/>");
    });
    ForceField pd;
    EXPECT_THROW(readForceField(fileName, &pd), gmx::InvalidInputError);
}

TEST(ForceFieldSimpleTest, changeChargeMinMaxForFixedParticle)
{
    // Load a fresh copy to avoid mutating the shared singleton
    std::string dataName = gmx::test::TestFileManager::getInputFilePath("ACM-g.xml");
    ForceField pd;
    readForceField(dataName, &pd);

    auto particle = pd.findParticleType("Be2+");
    ASSERT_NE(nullptr, particle);

    auto charge = particle->parameter("charge");
    ASSERT_NE(nullptr, charge);

    // The charge parameter in particletypes is Fixed
    EXPECT_FALSE(charge->isMutable());

    double origValue = charge->value();

    // Changing min/max must work regardless of mutability (Fixed),
    // since setMinimum/setMaximum are independent of mutability.
    EXPECT_TRUE(charge->setMinimum(2.0));
    charge->setMaximum(4.0);

    EXPECT_DOUBLE_EQ(2.0, charge->minimum());
    EXPECT_DOUBLE_EQ(4.0, charge->maximum());
    // Value must remain unchanged (Fixed)
    EXPECT_DOUBLE_EQ(origValue, charge->value());
}

}

} // namespace
