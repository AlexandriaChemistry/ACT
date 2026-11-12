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
            const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
            ForceField *mypd     = getForceField("ACM-g1");
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
            const std::string    atype("U");
            const std::string    elem("U");
            const std::string    desc("temporary test atom");
            const Identifier   atpId("U");
            ParticleType atp(atpId, desc, ActParticle::Atom);
            atp.setOption("poltype", "p_U");
            atp.setOption("zetatype", "z_U");
            atp.setOption("bondtype", "b_U");
            atp.setOption("acmtype", "z_U");
            atp.setOption("element", elem);
            atp.setOption("atomnumber", "92");
            const ForceFieldParameter mm("Da", 238.29, 0.0,  1, 230, 240, Mutability::Bounded, true, false);
            atp.addForceFieldParameter("mass", mm);
            const ForceFieldParameter rr("kJ/mol", 1000.0, 0.0,  1, 990, 1010, Mutability::Bounded, true, true);
            atp.addForceFieldParameter("ref_enthalpy", rr);
            pd->addParticleType(atp);
        }
};

std::vector<std::string> ForceFieldTest::atomNames;
std::string              ForceFieldTest::atomName;

TEST_F (ForceFieldTest, getAtype){
    auto aType = getForceField("ACM-g1")->findParticleType("h");

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
    alexandria::ForceField *pd = getForceField("ACM-g1");
    addAtype(pd);
    EXPECT_THROW(addAtype(pd), gmx::InvalidInputError);
}

TEST_F(ForceFieldTest, addAtype){
    alexandria::ForceField *pd = getForceField("ACM-g1");
    const std::string    atype("U");
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
    const std::vector<std::string> eqd = { "ACM-g1", "ACM-pg1" };

    for (auto &atom : atoms)
    {
        for (const auto& model : eqd)
        {
            auto pd       = getForceField(model);
            const auto fa       = pd->findParticleType(atom)->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            const auto eep      = pd->findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            const auto p        = eep.findParameterTypeConst(fa, "chi");
            chi0s.push_back(p.value());
        }
    }
    checker_.checkSequence(chi0s.begin(), chi0s.end(), "chi");
}

TEST_F (ForceFieldTest, alpha)
{
    std::vector<double>      alphas;
    std::string              atoms[]  = { "h", "c2", "f", "p2", "br" };
    const std::vector<std::string> eqd = { "ACM-pg1" };

    for (auto &atom : atoms)
    {
        for (const auto& model : eqd)
        {
            auto pd       = getForceField(model);
            const auto fa       = pd->findParticleType(atom)->interactionTypeToIdentifier(InteractionType::POLARIZATION);
            const auto eep      = pd->findForcesConst(InteractionType::POLARIZATION);
            const auto p        = eep.findParameterTypeConst(fa, "alpha");
            alphas.push_back(p.value());
        }
    }
    checker_.checkSequence(alphas.begin(), alphas.end(), "alpha");
}

TEST_F (ForceFieldTest, row)
{
    std::vector<int> rows;
    std::string      atoms[]  = { "h", "n2", "s2", "br" };
    std::string      models[] = { "ACM-g1", "ACM-pg1" };

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
    const std::vector<std::string> eqd = { "ACM-pg1" };

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
    const std::string dataName = gmx::test::TestFileManager::getInputFilePath("forcefields/ACM-g1.xml");
    const std::ifstream in(dataName);
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string xml = buffer.str();
    auto pos = xml.find("mutability=\"Bounded\"");
    ASSERT_NE(std::string::npos, pos);
    xml.replace(pos, std::string("mutability=\"Bounded\"").size(), "mutability=\"Free\"");
    gmx::test::TestFileManager fileManager;
    const std::string freeName = fileManager.getTemporaryFilePath("ACM-g-free.xml");
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
                                      const std::string                &base = "forcefields/ACM-g1.xml")
{
    const std::ifstream     in(gmx::test::TestFileManager::getInputFilePath(base));
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string xml = buffer.str();
    edit(&xml);
    auto fileName = fileManager->getTemporaryFilePath(name);
    std::ofstream out(fileName);
    out << xml;
    return fileName;
}

// Test: A force field with eta, chi, delta_eta or delta_chi in another unit
// than eV, or zeta in another unit than 1/nm, cannot be read, since the
// charge generation uses these values in those units as they are. Other
// parameters are converted with ForceFieldParameter::internalValue.
TEST(ForceFieldSimpleTest, unitsOfChargeGenerationChecked)
{
    gmx::test::TestFileManager fileManager;
    // The start of the parameter up to its unit, and another unit
    for (const auto &change : std::vector<std::pair<std::string, std::string>>{
             { "<parameter type=\"chi\" unit=\"", "kJ/mol e" },
             { "<parameter type=\"zeta\" unit=\"", "" } })
    {
        SCOPED_TRACE(change.first + change.second);
        const auto fileName = modifiedForceField(&fileManager, "units.xml", [&change](std::string *xml)
              {
                  auto pos = xml->find(change.first);
                  ASSERT_NE(std::string::npos, pos);
                  auto start = pos + change.first.size();
                  auto end   = xml->find("\"", start);
                  xml->replace(start, end - start, change.second);
              });
        ForceField pd;
        try
        {
            readForceField(fileName, &pd);
            ADD_FAILURE() << "The force field was read";
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find("cannot be used")) << ex.what();
        }
    }
    // Unchanged it can be read
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACM-g1.xml"), &pd);
}

// Test: A force field in which a type of virtual site has the function of
// another type cannot be read, since the virtual sites are constructed by
// their type with the parameters of their function.
TEST(ForceFieldSimpleTest, vsiteFunctionOfOtherTypeThrows)
{
    gmx::test::TestFileManager fileManager;
    const std::string          base = "forcefields/ACS-pg-vs3.xml";
    const auto fileName = modifiedForceField(&fileManager, "vsites.xml", [](std::string *xml)
          {
              const std::string from = "<interaction type=\"VSITE4S\" function=\"VSITE4S\"";
              auto              pos  = xml->find(from);
              ASSERT_NE(std::string::npos, pos);
              xml->replace(pos, from.size(), "<interaction type=\"VSITE4S\" function=\"VSITE4S3\"");
          }, base);
    ForceField pd;
    try
    {
        readForceField(fileName, &pd);
        ADD_FAILURE() << "The force field was read";
    }
    catch (const gmx::InvalidInputError &ex)
    {
        EXPECT_NE(std::string::npos, std::string(ex.what()).find("VSITE4S have function VSITE4S3")) << ex.what();
    }
    // Unchanged it can be read
    ForceField unchanged;
    readForceField(gmx::test::TestFileManager::getInputFilePath(base), &unchanged);
}

// Test: A particle type of which the ACM type has no parameters cannot be
// read, unless without the check of consistency
TEST(ForceFieldSimpleTest, missingAcmParametersThrows)
{
    gmx::test::TestFileManager fileManager;
    const auto fileName = modifiedForceField(&fileManager, "no-iodine.xml", [](std::string *xml)
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
    const auto fileName = modifiedForceField(&fileManager, "no-acmtype.xml", [](std::string *xml)
          {
              auto        ptype  = xml->find("<particletype identifier=\"i\"");
              const std::string option("      <option key=\"acmtype\" value=\"i_z\"/>\n");
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
    const auto fileName = modifiedForceField(&fileManager, "point.xml", [](std::string *xml)
          {
              const std::string gaussian("<option key=\"chargedistribution\" value=\"Gaussian\"/>");
              auto        pos = xml->find(gaussian);
              ASSERT_NE(std::string::npos, pos);
              xml->replace(pos, gaussian.size(), "<option key=\"chargedistribution\" value=\"Point\"/>");
          });
    ForceField pd;
    EXPECT_THROW(readForceField(fileName, &pd), gmx::InvalidInputError);
}

// Test: The same for a force field file in which the charge distribution
// is the function of the interaction COULOMB
TEST(ForceFieldSimpleTest, pointChargesWithZetaThrowOldStyle)
{
    gmx::test::TestFileManager fileManager;
    const auto fileName = modifiedForceField(&fileManager, "point-old.xml", [](std::string *xml)
          {
              const std::string gaussian("<interaction type=\"COULOMB\" function=\"COULOMB_GAUSSIAN\"");
              auto        pos = xml->find(gaussian);
              ASSERT_NE(std::string::npos, pos);
              xml->replace(pos, gaussian.size(), "<interaction type=\"COULOMB\" function=\"COULOMB_POINT\"");
          }, "forcefields/ACM-g2.xml");
    ForceField pd;
    EXPECT_THROW(readForceField(fileName, &pd), gmx::InvalidInputError);
}

//! A force field file in which the charge distributions are those of COULOMB
static const char oldStyle[] = "forcefields/ACM-g2.xml";

//! Replace the first occurrence of a text, which has to be present
static void replaceFirst(std::string *xml, const std::string &from, const std::string &to)
{
    auto pos = xml->find(from);
    ASSERT_NE(std::string::npos, pos) << from;
    xml->replace(pos, from.size(), to);
}

// Test: In an older force field file the width of the charge distribution is
// a parameter of COULOMB to which the zetatype of a particle type refers, and
// the number of exclusions and the charge distribution are options of COULOMB.
// They become a parameter of each particle type and options of the particle
// types, and COULOMB and the zetatypes are removed.
TEST(ForceFieldSimpleTest, OldStyleCoulomb)
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath(oldStyle), &pd);
    EXPECT_FALSE(pd.interactionPresent(InteractionType::ELECTROSTATICS));
    EXPECT_EQ(ChargeDistributionType::Gaussian, pd.chargeDistribution());
    EXPECT_EQ(Potential::COULOMB_GAUSSIAN, pd.potential(InteractionType::ELECTROSTATICS));
    EXPECT_EQ(0, pd.nexcl());
    EXPECT_DOUBLE_EQ(1, pd.epsilonr());
    EXPECT_FALSE(pd.findForcesConst(InteractionType::VDW).optionExists("nexcl"));
    for (const auto &pt : pd.particleTypesConst())
    {
        EXPECT_FALSE(pt.second.hasOption("zetatype")) << pt.first.id();
    }
    auto be = pd.findParticleType("Be2+");
    ASSERT_TRUE(be->hasParameter(zetaName));
    const auto &zeta = be->parameterConst(zetaName);
    EXPECT_DOUBLE_EQ(13, zeta.value());
    EXPECT_DOUBLE_EQ(6, zeta.minimum());
    EXPECT_DOUBLE_EQ(20, zeta.maximum());
    EXPECT_EQ(Mutability::Bounded, zeta.mutability());
    EXPECT_EQ(1u, zeta.ntrain());
    EXPECT_DOUBLE_EQ(12, pd.findParticleType("c3")->zeta());
}

// Test: Particle types that share a zetatype get a copy of its zeta each
TEST(ForceFieldSimpleTest, OldStyleSharedZetatype)
{
    gmx::test::TestFileManager fileManager;
    const auto fileName = modifiedForceField(&fileManager, "shared.xml", [](std::string *xml)
          {
              replaceFirst(xml, "<option key=\"zetatype\" value=\"Be2+_z\"/>", "<option key=\"zetatype\" value=\"c3_z\"/>");
          }, oldStyle);
    ForceField pd;
    readForceField(fileName, &pd);
    auto be = pd.findParticleType("Be2+");
    auto c3 = pd.findParticleType("c3");
    EXPECT_DOUBLE_EQ(12, be->zeta());
    EXPECT_DOUBLE_EQ(12, c3->zeta());
    be->parameters()->find(zetaName)->second.setValue(15);
    EXPECT_DOUBLE_EQ(15, be->zeta());
    EXPECT_DOUBLE_EQ(12, c3->zeta());
}

// Test: An older force field file with a zeta both in COULOMB and in the
// particle type, with different values, cannot be read
TEST(ForceFieldSimpleTest, OldStyleConflictingZetaThrows)
{
    gmx::test::TestFileManager fileManager;
    const auto fileName = modifiedForceField(&fileManager, "conflict.xml", [](std::string *xml)
          {
              const std::string zetatype("<option key=\"zetatype\" value=\"Be2+_z\"/>");
              replaceFirst(xml, zetatype, zetatype + "\n      <parameter type=\"zeta\" unit=\"1/nm\" value=\"14\" uncertainty=\"0\" minimum=\"6\" maximum=\"20\" ntrain=\"1\" mutability=\"Bounded\" nonnegative=\"yes\"/>");
          }, oldStyle);
    ForceField pd;
    EXPECT_THROW(readForceField(fileName, &pd), gmx::InvalidInputError);
}

// Test: An older force field file with different numbers of exclusions for
// COULOMB and VANDERWAALS cannot be read
TEST(ForceFieldSimpleTest, OldStyleDifferentExclusionsThrow)
{
    gmx::test::TestFileManager fileManager;
    const auto fileName = modifiedForceField(&fileManager, "nexcl.xml", [](std::string *xml)
          {
              // The first is that of VANDERWAALS, which precedes COULOMB
              replaceFirst(xml, "<option key=\"nexcl\" value=\"0\"/>", "<option key=\"nexcl\" value=\"1\"/>");
          }, oldStyle);
    ForceField pd;
    EXPECT_THROW(readForceField(fileName, &pd), gmx::InvalidInputError);
}

// Test: A force field written after reading an older file has the new form,
// and reading it gives the same parameters and options
TEST(ForceFieldSimpleTest, OldStyleRoundTrip)
{
    gmx::test::TestFileManager fileManager;
    ForceField old, roundTrip;
    readForceField(gmx::test::TestFileManager::getInputFilePath(oldStyle), &old);
    const auto fileName = fileManager.getTemporaryFilePath("roundtrip.xml");
    writeForceField(fileName, &old, false);
    const std::ifstream     in(fileName);
    std::stringstream buffer;
    buffer << in.rdbuf();
    EXPECT_EQ(std::string::npos, buffer.str().find("zetatype"));
    EXPECT_EQ(std::string::npos, buffer.str().find("type=\"COULOMB\""));
    readForceField(fileName, &roundTrip);
    EXPECT_EQ(old.particleTypesOptions(), roundTrip.particleTypesOptions());
    ASSERT_EQ(old.particleTypesConst().size(), roundTrip.particleTypesConst().size());
    for (const auto &pt : old.particleTypesConst())
    {
        auto rt = roundTrip.findParticleType(pt.first.id());
        EXPECT_EQ(pt.second.optionsConst(), rt->optionsConst()) << pt.first.id();
        for (const auto &p : pt.second.parametersConst())
        {
            ASSERT_TRUE(rt->hasParameter(p.first)) << pt.first.id() << " " << p.first;
            const auto &q = rt->parameterConst(p.first);
            EXPECT_EQ(p.second.value(), q.value()) << pt.first.id() << " " << p.first;
            EXPECT_EQ(p.second.minimum(), q.minimum()) << pt.first.id() << " " << p.first;
            EXPECT_EQ(p.second.maximum(), q.maximum()) << pt.first.id() << " " << p.first;
            EXPECT_EQ(p.second.mutability(), q.mutability()) << pt.first.id() << " " << p.first;
            EXPECT_EQ(p.second.ntrain(), q.ntrain()) << pt.first.id() << " " << p.first;
        }
    }
}

// Test: In a force field file of the new form a particle type without zeta
// has point charges, and without the options of the particle types there are
// point charges, epsilon_r is one and the number of exclusions is unknown
TEST(ForceFieldSimpleTest, NewStyleDefaults)
{
    gmx::test::TestFileManager fileManager;
    const auto noZeta = modifiedForceField(&fileManager, "nozeta.xml", [](std::string *xml)
          {
              auto be  = xml->find("identifier=\"Be2+\"");
              ASSERT_NE(std::string::npos, be);
              auto pos = xml->find("<parameter type=\"zeta\"", be);
              ASSERT_LT(pos, xml->find("</particletype>", be));
              xml->erase(pos, xml->find("/>", pos) + 2 - pos);
          });
    ForceField pd;
    readForceField(noZeta, &pd);
    EXPECT_FALSE(pd.findParticleType("Be2+")->hasParameter(zetaName));
    EXPECT_EQ(0, pd.findParticleType("Be2+")->zeta());
    EXPECT_EQ(ChargeDistributionType::Gaussian, pd.chargeDistribution());
    ForceField none;
    none.setFilename("none.xml");
    EXPECT_EQ(ChargeDistributionType::Point, none.chargeDistribution());
    EXPECT_DOUBLE_EQ(1, none.epsilonr());
    EXPECT_EQ(-1, none.nexcl());
}

// Test: An unknown option of the particle types or an invalid value cannot be read
TEST(ForceFieldSimpleTest, NewStyleInvalidOptionThrows)
{
    gmx::test::TestFileManager fileManager;
    for (const auto &opt : { std::pair<std::string, std::string>("unknown", "1"),
                             { "nexcl", "two" }, { "chargedistribution", "Cubic" } })
    {
        const auto fileName = modifiedForceField(&fileManager, "option.xml", [&opt](std::string *xml)
              {
                  const std::string nexcl("<option key=\"nexcl\" value=\"0\"/>");
                  replaceFirst(xml, nexcl, nexcl + "\n    <option key=\"" + opt.first + "\" value=\"" + opt.second + "\"/>");
              });
        ForceField pd;
        EXPECT_THROW(readForceField(fileName, &pd), gmx::InvalidInputError) << opt.first;
    }
}

TEST(ForceFieldSimpleTest, changeChargeMinMaxForFixedParticle)
{
    // Load a fresh copy to avoid mutating the shared singleton
    const std::string dataName = gmx::test::TestFileManager::getInputFilePath("forcefields/ACM-g1.xml");
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
