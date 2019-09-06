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
#include <map>
 
#include <math.h>

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
            Poldata *mypd = getPoldata("ACM-g");
            auto atomName = mypd->findAtype("ha")->getType();
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
    checker_.checkString(aType->getPtype(), "ptype");
    checker_.checkString(aType->getBtype(), "btype");
    checker_.checkString(aType->getZtype(), "ztype");
    checker_.checkString(aType->getVdwparams(), "vdwparams");
    checker_.checkString(aType->getRefEnthalpy(), "refEnthalpy");
}

TEST_F(PoldataTest, addAtype){
    const std::string        elem         = "U";
    const std::string        desc         = "temporary test atom";
    const std::string        atype        = "U";
    const std::string        ptype        = "p_U";
    const std::string        ztype        = "z_U";
    const std::string        btype        = "b_U";
          std::string        vdwparams    = "10.0 11.1 12.2";
    const std::string        ref_enthalpy = "1000";
    bool                     fixVdw       = true;
    alexandria::Poldata     *pd           = getPoldata("ACM-g");
    pd->addAtype(elem,
                 desc,
                 atype,
                 ptype,
                 btype,
                 ztype,
                 fixVdw,
                 vdwparams,
                 ref_enthalpy);

    auto fa = pd->findAtype(atype);
    if (fa != pd->getAtypeEnd())
    {
        /* Test if the extractions where correct */
        checker_.checkString(fa->getElem(), elem.c_str());
        checker_.checkString(fa->getDesc(), desc.c_str());
        checker_.checkString(fa->getType(), atype.c_str());
        checker_.checkString(fa->getPtype(), ptype.c_str());
        checker_.checkString(fa->getBtype(), btype.c_str());
        checker_.checkString(fa->getZtype(), ztype.c_str());
        checker_.checkString(fa->getVdwparams(), vdwparams.c_str());
    }
    EXPECT_THROW_GMX(fa->setVdwparams(vdwparams), gmx::InvalidInputError);
    fa->setFixed(false);
    fa->setVdwparams(vdwparams);
}

TEST_F (PoldataTest, Ptype)
{
    auto pd    = getPoldata("ACM-pg");
    auto ptype = pd->findPtype("p_ha");
    if (ptype != pd->getPtypeEnd())
    {
        checker_.checkString(ptype->getType(), "type");
        checker_.checkString(ptype->getMiller(), "miller");
        checker_.checkString(ptype->getBosque(), "bosque");
        checker_.checkDouble(ptype->getPolarizability(), "polarizability");
        checker_.checkDouble(ptype->getSigPol(), "sigPol");
    }
}

TEST_F (PoldataTest, Miller)
{
    auto pd    = getPoldata("ACM-g");
    alexandria::MillerIterator miller = pd->getMillerBegin();
    checker_.checkInteger(miller->getAtomnumber(), "atomnumber");
    checker_.checkDouble(miller->getTauAhc(), "tauAhc");
    checker_.checkDouble(miller->getAlphaAhp(), "alphaAhp");
}


TEST_F (PoldataTest, Bosque)
{
    auto pd    = getPoldata("ACM-g");
    alexandria::BosqueIterator bosque = pd->getBosqueBegin();
    checker_.checkString(bosque->getBosque(), "bosque");
    checker_.checkDouble(bosque->getPolarizability(), "polarizability");
}

TEST_F (PoldataTest, chi)
{
    std::vector<double>                  chi0s;
    std::vector<ChargeModel> eqd = { eqdACM_g, eqdACM_pg };

    for (auto model : eqd)
    {
        auto pd       = getPoldata(model);
        auto atomName = pd->findAtype("ha")->getType();
        chi0s.push_back(pd->getChi0(atomName));
    }
    checker_.checkSequence(chi0s.begin(), chi0s.end(), "chi");
}

TEST_F (PoldataTest, row){
    std::vector<int> rows;
    std::string atoms[] = { "ha", "s", "br" };
    std::string models[] = { "ESP-ps", "ACM-g", "Rappe" };

    for (auto &atom : atoms)
    {
        for (auto &model : models)
        {
            auto pd       = getPoldata(model);
            auto atomName = pd->findAtype(atom)->getType();
            rows.push_back(pd->getRow(atomName, 0));
        }
    }
    checker_.checkSequence(rows.begin(), rows.end(), "row");
}

TEST_F (PoldataTest, zeta)
{
    std::vector<double>                  zetas;
    std::vector<ChargeModel> eqd = { eqdACM_pg, eqdESP_ps };

    for (auto model : eqd)
    {        
        auto pd       = getPoldata(model);
        auto atomName = pd->findAtype("ha")->getType();
        for(int i=0; i < pd->getNzeta(atomName); i++)
        {
            zetas.push_back(pd->getZeta(atomName, i));
        }
    }
    checker_.checkSequence(zetas.begin(), zetas.end(), "zeta");
}

TEST_F (PoldataTest, chargeModel)
{
    std::vector<ChargeModel> eqd = { eqdESP_pp, eqdESP_pg, eqdESP_ps };

    std::vector<std::string> forces;
    for (auto model : eqd)
    {        
        auto mypd = getPoldata(model);
        forces.push_back(getEemtypeName(mypd->getChargeModel()));
    }
    checker_.checkSequence(forces.begin(), forces.end(), "chargeModel");
}


TEST_F (PoldataTest, lenghtUnit)
{
    auto fs = getPoldata("ACM-g")->findForces(alexandria::eitBONDS);
    std::string length =  fs->unit();
    checker_.checkString(length, "lenghtUnit");
}

TEST_F (PoldataTest, polarUnit)
{
    std::string polarUnit = getPoldata("ACM-g")->getPolarUnit( );
    checker_.checkString(polarUnit, "polarUnit");
}


TEST_F (PoldataTest, polarRef)
{
    std::string polarRef =  getPoldata("ACM-g")->getPolarRef( );
    checker_.checkString(polarRef, "polarRef");
}

}

}
