/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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

#include "programs/alexandria/plistwrapper.h"
#include "programs/alexandria/poldata.h"
#include "programs/alexandria/poldata_low.h"
#include "programs/alexandria/poldata_xml.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "poldata_utils.h"

namespace alexandria
{

namespace
{

class PoldataTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
        static std::vector<std::string> atomNames;
        static std::string              atomName;

        PoldataTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
            Poldata *mypd     = getPoldata("ACM-g");
            for (const auto &iter : mypd->particleTypesConst())
            {
                atomNames.push_back(iter.id().id());
            }
        }

        static void TearDownTestCase()
        {
        }


};

std::vector<std::string> PoldataTest::atomNames;
std::string              PoldataTest::atomName;

TEST_F (PoldataTest, getAtype){
    auto aType = getPoldata("ACM-g")->findParticleType("h1");

    checker_.checkString(aType->element(), "element");
    checker_.checkString(aType->description(), "description");
    checker_.checkString(aType->id().id(), "type");
    checker_.checkString(aType->interactionTypeToIdentifier(InteractionType::POLARIZATION).id().c_str(), "poltype");
    checker_.checkString(aType->interactionTypeToIdentifier(InteractionType::BONDS).id().c_str(), "bondtype");
    checker_.checkString(aType->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION).id().c_str(), "zetatype");
    checker_.checkDouble(aType->refEnthalpy(), "refEnthalpy");
}

TEST_F(PoldataTest, addAtype){
    std::string          atype("U");
    std::string          elem("U");
    std::string          desc("temporary test atom");
    alexandria::Poldata *pd = getPoldata("ACM-g");
    Identifier   atpId({"U"}, CanSwap::No);
    ParticleType atp(atpId, desc, eptAtom);
    atp.setOption("poltype", "p_U");
    atp.setOption("zetatype", "z_U");
    atp.setOption("bondtype", "b_U");
    atp.setOption("acmtype", "z_U");
    atp.setOption("element", elem);
    atp.setOption("atomnumber", "92");
    ForceFieldParameter mm("Da", 238.29, 0.0,  1, 230, 240, Mutability::Free, true);
    atp.addForceFieldParameter("mass", mm);
    ForceFieldParameter rr("kJ/mol", 1000.0, 0.0,  1, 990, 1010, Mutability::Free, true);
    atp.addForceFieldParameter("ref_enthalpy", rr);
    pd->addParticleType(atp);

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

TEST_F (PoldataTest, Verstraelen)
{
    auto pd  = getPoldata("Verstraelen");
    auto bcs = pd->findForcesConst(InteractionType::BONDCORRECTIONS);
    std::vector<std::string> name;
    std::vector<double>      hardness;
    std::vector<double>      electronegativity;
    for ( auto bc : bcs.parametersConst() )
    {
        name.push_back(bc.first.id());
        hardness.push_back(bcs.findParameterTypeConst(bc.first, "hardness").value());
        electronegativity.push_back(bcs.findParameterTypeConst(bc.first, "electronegativity").value());
    }
    checker_.checkSequence(name.begin(), name.end(), "name");
    checker_.checkSequence(hardness.begin(), hardness.end(), "hardness");
    checker_.checkSequence(electronegativity.begin(), electronegativity.end(),
                           "electronegativity");
}

TEST_F (PoldataTest, chi)
{
    std::vector<double>      chi0s;
    std::string              atoms[]  = { "hp", "f", "p2", "br" };
    std::vector<std::string> eqd = { "ACM-g", "ACM-pg" };

    for (auto &atom : atoms)
    {
        for (auto model : eqd)
        {
            auto pd       = getPoldata(model);
            auto fa       = pd->findParticleType(atom)->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            auto eep      = pd->findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            auto p        = eep.findParameterTypeConst(fa, "chi");
            chi0s.push_back(p.value());
        }
    }
    checker_.checkSequence(chi0s.begin(), chi0s.end(), "chi");
}

TEST_F (PoldataTest, alpha)
{
    std::vector<double>      alphas;
    std::string              atoms[]  = { "hp", "c2", "f", "p2", "br" };
    std::vector<std::string> eqd = { "ACM-pg" };

    for (auto &atom : atoms)
    {
        for (auto model : eqd)
        {
            auto pd       = getPoldata(model);
            auto fa       = pd->findParticleType(atom)->interactionTypeToIdentifier(InteractionType::POLARIZATION);
            auto eep      = pd->findForcesConst(InteractionType::POLARIZATION);
            auto p        = eep.findParameterTypeConst(fa, "alpha");
            alphas.push_back(p.value());
        }
    }
    checker_.checkSequence(alphas.begin(), alphas.end(), "alpha");
}

TEST_F (PoldataTest, row)
{
    std::vector<int> rows;
    std::string      atoms[]  = { "ha", "n", "s", "br" };
    std::string      models[] = { "ESP-ps", "ACM-g", "Rappe" };

    for (auto &atom : atoms)
    {
        for (auto &model : models)
        {
            auto pd  = getPoldata(model);
            auto fa  = pd->findParticleType(atom);
            auto p   = fa->row();
            rows.push_back(p);
        }
    }
    checker_.checkSequence(rows.begin(), rows.end(), "row");
}

TEST_F (PoldataTest, zeta)
{
    std::vector<double>      zetas;
    std::string              atoms[]  = { "h1", "c3", "s4", "br" };
    std::vector<std::string> eqd = { "ACM-pg", "ESP-ps" };

    for (auto &atom : atoms)
    {
        for (auto &model : eqd)
        {
            auto pd  = getPoldata(model);
            auto fa  = pd->findParticleType(atom);
            auto ztp = fa->interactionTypeToIdentifier(InteractionType::CHARGEDISTRIBUTION);
            auto eep = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
            auto p   = eep.findParameterTypeConst(ztp, "zeta");
            zetas.push_back(p.value());
        }
    }
    checker_.checkSequence(zetas.begin(), zetas.end(), "zeta");
}

TEST_F (PoldataTest, chargeType)
{
    std::vector<std::string> eqd = { "ESP-pp", "ESP-pg", "ESP-ps" };

    std::vector<std::string> forces;
    for (auto &model : eqd)
    {
        auto pd = getPoldata(model);
        auto qt = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
        auto ct = name2ChargeType(qt.optionValue("chargetype"));
        forces.push_back(chargeTypeName(ct));
    }
    checker_.checkSequence(forces.begin(), forces.end(), "chargeType");
}

}

} // namespace
