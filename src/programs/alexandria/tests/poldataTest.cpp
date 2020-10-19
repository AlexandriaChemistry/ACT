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
            auto     atomName = mypd->findAtype("ha")->getType();
            for (auto iter = mypd->getAtypeBegin(); iter != mypd->getAtypeEnd(); iter++)
            {
                atomNames.push_back(iter->getType());
            }
        }

        static void TearDownTestCase()
        {
        }


};

std::vector<std::string> PoldataTest::atomNames;
std::string              PoldataTest::atomName;

TEST_F (PoldataTest, getAtype){
    alexandria::FfatypeIterator aType =  getPoldata("ACM-g")->findAtype("h1");

    checker_.checkString(aType->getElem(), "elem");
    checker_.checkString(aType->getDesc(), "desc");
    checker_.checkString(aType->getType(), "type");
    checker_.checkString(aType->id(eitPOLARIZATION).id(), "ptype");
    checker_.checkString(aType->id(eitBONDS).id(), "btype");
    checker_.checkString(aType->id(eitELECTRONEGATIVITYEQUALIZATION).id(), "ztype");
    checker_.checkString(aType->getRefEnthalpy(), "refEnthalpy");
}

TEST_F(PoldataTest, addAtype){
    const std::string        elem         = "U";
    const std::string        desc         = "temporary test atom";
    const std::string        atype        = "U";
    const std::string        ptype        = "p_U";
    const std::string        ztype        = "z_U";
    const std::string        btype        = "b_U";
    const std::string        ref_enthalpy = "1000";
    alexandria::Poldata     *pd           = getPoldata("ACM-g");
    pd->addAtype(elem,
                 desc,
                 atype,
                 ptype,
                 btype,
                 ztype,
                 ref_enthalpy);

    auto fa = pd->findAtype(atype);
    if (fa != pd->getAtypeEnd())
    {
        /* Test if the extractions where correct */
        checker_.checkString(fa->getElem(), elem.c_str());
        checker_.checkString(fa->getDesc(), desc.c_str());
        checker_.checkString(fa->getType(), atype.c_str());
        checker_.checkString(fa->id(eitPOLARIZATION).id(), "ptype");
        checker_.checkString(fa->id(eitBONDS).id(), "btype");
        checker_.checkString(fa->id(eitELECTRONEGATIVITYEQUALIZATION).id(), "ztype");
    }
}

TEST_F (PoldataTest, Verstraelen)
{
    auto pd  = getPoldata("Verstraelen");
    auto bcs = pd->findForcesConst(eitBONDCORRECTIONS);
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
            auto fa       = pd->findAtype(atom)->id(eitELECTRONEGATIVITYEQUALIZATION);
            auto eep      = pd->findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
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
            auto fa       = pd->findAtype(atom)->id(eitPOLARIZATION);
            auto eep      = pd->findForcesConst(eitPOLARIZATION);
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
            auto fa  = pd->findAtype(atom)->id(eitELECTRONEGATIVITYEQUALIZATION);
            auto eep = pd->findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
            auto p   = eep.findParameterTypeConst(fa, "row");
            rows.push_back(p.value());
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
            auto fa  = pd->findAtype(atom);
            auto ztp = fa->id(eitELECTRONEGATIVITYEQUALIZATION);
            auto eep = pd->findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
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
        auto ct = getPoldata(model)->chargeType();
        forces.push_back(chargeTypeName(ct));
    }
    checker_.checkSequence(forces.begin(), forces.end(), "chargeType");
}

}

} // namespace
