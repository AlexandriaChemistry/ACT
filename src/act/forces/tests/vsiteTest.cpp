/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2025
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
/*! \internal \file
 * \brief
 * Implements test of bonded force routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author A. Najla Hosseini <atiyeh.hosseini@icm.uu.se>
 * \ingroup group_forces_tests
 */
#include "actpre.h"

#include "../vsitehandler.h"

#include <cmath>
#include <memory>

#include <gtest/gtest.h>

#include "act/alexandria/actmol.h"
#include "act/alexandria/babel_io.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forces/forcecomputer.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class VsiteHandlerTest : public gmx::test::CommandLineTestBase
{
private:
    matrix box = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    real   dt  = 0.001;
protected:
    bool setup(const std::string &ff,
               const std::string &molName,
               ACTMol            *mol)
    {
        // Get forcefield
        auto pd  = getForceField(ff);
        
        std::vector<alexandria::MolProp> molprops;
        bool   userqtot   = false;
        double qtot_babel = 0;
        const char           *conf      = (char *)"minimum";
        const char           *jobtype   = (char *)"Opt";
        int                   maxpot    = 100;
        int                   nsymm     = 0;
        std::string           method;
        std::string           basis;
        std::string fileName = gmx::formatString("%s.sdf", molName.c_str());
        std::string dataName = gmx::test::TestFileManager::getInputFilePath(fileName);
        EXPECT_TRUE(readBabel(nullptr, pd, dataName.c_str(), &molprops,
                              molName.c_str(), molName.c_str(),
                              conf, &method, &basis, maxpot,
                              nsymm, jobtype, userqtot ,&qtot_babel,
                              false, box, true));
        EXPECT_TRUE(molprops.size() == 1);
        auto &molprop = molprops[0];
        alexandria::ACTMol               mp_;

        mol->Merge(&molprop);
        MsgHandler msghandler;
        // Uncomment in case of issues
        msghandler.setPrintLevel(ACTStatus::Debug);
        // Generate charges and topology
        mol->GenerateTopology(&msghandler, pd, missingParameters::Ignore);

        return msghandler.ok();
    }
    void testX(const std::string &ff,
               const std::string &molName)
    {
        ACTMol mp_;

        bool ok = setup(ff, molName, &mp_);
        EXPECT_TRUE(ok);
        if (ok)
        {
            std::vector<gmx::RVec> coords = mp_.xOriginal();

            VsiteHandler vh;
            vh.init(box, dt);
            vh.constructPositions(mp_.topology(), &coords, box);

            gmx::test::TestReferenceChecker checker_(this->rootChecker());
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
            auto  myatoms = mp_.topology()->atoms();
            const char *xyz[DIM] = { "X", "Y", "Z" };
            for(size_t i = 0; i < coords.size(); i++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    auto label = gmx::formatString("coords[%s-%zu][%s]", myatoms[i].name().c_str(), i, xyz[m]);
                    checker_.checkReal(coords[i][m], label.c_str());
                }
            }
        }
    }
    void testF(const std::string &ff,
               const std::string &molName)
    {
        ACTMol mp_;
        bool ok = setup(ff, molName, &mp_);
        EXPECT_TRUE(ok);
        if (ok)
        {
            std::vector<gmx::RVec> coords = mp_.xOriginal();

            VsiteHandler vh;
            vh.init(box, dt);
            vh.constructPositions(mp_.topology(), &coords, box);
            gmx::RVec fzero = { 0, 0, 0 };
            std::vector<gmx::RVec> forces(coords.size(), fzero);
            auto  myatoms = mp_.topology()->atoms();
            size_t i = 0;
            int    m = 0;
            for (const auto &a : myatoms)
            {
                if (a.pType() == ActParticle::Vsite)
                {
                    forces[i][m] = 1;
                    m = (m+1) % DIM;
                }
                i += 1;
            }
            vh.distributeForces(mp_.topology(), coords, &forces, box);
            gmx::test::TestReferenceChecker checker_(this->rootChecker());
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
            const char *xyz[DIM] = { "X", "Y", "Z" };

            for(size_t i = 0; i < forces.size(); i++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    auto label = gmx::formatString("forces[%s-%zu][%s]", myatoms[i].name().c_str(), i, xyz[m]);
                    checker_.checkReal(forces[i][m], label.c_str());
                }
            }
        }
    }
};

TEST_F (VsiteHandlerTest, VSite2HF)
{
    testX("ACS-pg-vs2", "hydrogen-fluoride");
}

TEST_F (VsiteHandlerTest, VSite2HBr)
{
    testX("ACS-pg-vs2", "hydrogen-bromide");
}

TEST_F (VsiteHandlerTest, VSite3H20)
{
    testX("ACS-pg-vs3", "water");
}

TEST_F (VsiteHandlerTest, VSite2HFForce)
{
    testF("ACS-pg-vs2", "hydrogen-fluoride");
}

TEST_F (VsiteHandlerTest, VSite2HBrForce)
{
    testF("ACS-pg-vs2", "hydrogen-bromide");
}

TEST_F (VsiteHandlerTest, VSite3H20Force)
{
    testF("ACS-pg-vs3", "water");
}

TEST(Vsite1, HF)
{
    std::string forcefield("ACS-pg-vs2");
    // Get forcefield
    auto pd  = getForceField(forcefield);

    auto fh = Identifier({ "f", "v1f1" }, { 9 }, CanSwap::Yes);

    auto itype = InteractionType::VSITE1;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(fh));
    auto vsite1_name = potentialToParameterName(Potential::VSITE1);
    auto param_fh = fs.findParameterType(fh, vsite1_name[vsite1A]);
    EXPECT_TRUE(param_fh->internalValue() == 1);
    auto fh2 = Identifier({ "g", "v2f1" }, { 9 }, CanSwap::Yes);
    EXPECT_FALSE(fs.parameterExists(fh2));
    
}

TEST(Vsite2, HF)
{
    std::string forcefield("ACS-pg-vs2");
    // Get forcefield
    auto pd  = getForceField(forcefield);

    auto fh = Identifier({ "f_b", "h_b", "v2" }, { 1, 9 }, CanSwap::Vsite2);
    auto hf = Identifier({ "h_b", "f_b", "v2" }, { 1, 9 }, CanSwap::Vsite2);
    // Compare identifiers
    EXPECT_TRUE(fh == hf);
    auto itype = InteractionType::VSITE2;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(fh));
    auto vsite2_name = potentialToParameterName(Potential::VSITE2);
    auto param_fh = fs.findParameterType(fh, vsite2_name[vsite2A]);
    EXPECT_TRUE(param_fh->internalValue() == -1.05);
    if (fs.parameterExists(hf))
    {
        auto param_hf = fs.findParameterType(hf, vsite2_name[vsite2A]);
        EXPECT_TRUE(param_fh->internalValue() == param_hf->internalValue());
    }
}

TEST(Vsite2FD, HF)
{
    std::string forcefield("ACS-pg-vs2");
    // Get forcefield
    auto pd  = getForceField(forcefield);

    auto fh = Identifier({ "f_b", "h_b", "v2" }, { 1, 9 }, CanSwap::Vsite2);
    auto hf = Identifier({ "h_b", "f_b", "v2" }, { 1, 9 }, CanSwap::Vsite2);
    // Compare identifiers
    EXPECT_TRUE(fh == hf);
    auto itype = InteractionType::VSITE2;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(fh));
    auto vsite2_name = potentialToParameterName(Potential::VSITE2);
    auto param_fh = fs.findParameterType(fh, vsite2_name[vsite2A]);
    EXPECT_TRUE(param_fh->internalValue() == -1.05);
    if (fs.parameterExists(hf))
    {
        auto param_hf = fs.findParameterType(hf, vsite2_name[vsite2A]);
        EXPECT_TRUE(param_fh->internalValue() == param_hf->internalValue());
    }
}

TEST(Vsite2, HFCanSwapVsite2)
{
    std::string forcefield("ACS-pg-vs2");
    // Get forcefield
    auto pd  = getForceField(forcefield);

    auto fh = Identifier({ "f_b", "h_b", "v2" }, { 1, 9 }, CanSwap::Vsite2);
    auto hf = Identifier({ "h_b", "f_b", "v2" }, { 1, 9 }, CanSwap::Vsite2);
    // Compare identifiers
    EXPECT_TRUE(fh == hf);
    auto itype = InteractionType::VSITE2;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(fh));
    auto vsite2_name = potentialToParameterName(Potential::VSITE2);
    auto param_fh = fs.findParameterType(fh, vsite2_name[vsite2A]);
    EXPECT_TRUE(param_fh->internalValue() == -1.05);
    EXPECT_TRUE(fs.parameterExists(hf));
}

TEST(Vsite3, hoh)
{
    std::string forcefield("ACS-pg-vs3");
    // Get forcefield
    auto pd  = getForceField(forcefield);

    auto hoh = Identifier({ "h_b", "o3_b", "h_b", "v3" }, { 1, 1, 9 }, CanSwap::Vsite3);
    // Compare identifiers
    EXPECT_TRUE(hoh == hoh);
    auto itype = InteractionType::VSITE3;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(hoh));
    auto vsite3_name = potentialToParameterName(Potential::VSITE3);
    auto param_hoh = fs.findParameterType(hoh, vsite3_name[vsite3A]);
    EXPECT_TRUE(param_hoh->internalValue() == -0.2);
    if (fs.parameterExists(hoh))
    {
        auto param_hoh = fs.findParameterType(hoh, vsite3_name[vsite3A]);
        EXPECT_TRUE(param_hoh->internalValue() == param_hoh->internalValue());
    }
}

TEST(VSite3, hohCanSwapVsite3)
{
    std::string forcefield("ACS-pg-vs3");
    // Get the forcefield
    auto pd = getForceField(forcefield);

    auto hoh = Identifier({ "h_b", "o3_b", "h_b", "v3" }, { 1, 1, 9 }, CanSwap::Vsite3);
    auto coh = Identifier({ "c3_b", "o3_b", "h_b", "v3" }, { 1, 1, 9 }, CanSwap::Vsite3);
    // Compare identifiers
    //9 means v-site
    EXPECT_FALSE(hoh == coh);
    auto itype = InteractionType::VSITE3;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(hoh));
    auto vsite3_name = potentialToParameterName(Potential::VSITE3);
    auto param_hoh = fs.findParameterType(hoh, vsite3_name[vsite3A]);
    EXPECT_TRUE(param_hoh->internalValue() == -0.2);
    EXPECT_FALSE(fs.parameterExists(coh));
}

TEST(VSite3fd, hohCanSwapNo)
{
    std::string forcefield("ACS-v3fds");
    // Get the forcefield
    auto pd = getForceField(forcefield);

    auto hoh = Identifier({ "h_b", "o3_b", "h_b", "v3Ow" }, { 1, 1, 9 }, CanSwap::No);
    auto coh = Identifier({ "c3_b", "o3_b", "h_b", "v3Ow" }, { 1, 1, 9 }, CanSwap::No);
    // Compare identifiers
    //9 means v-site
    EXPECT_FALSE(hoh == coh);
    auto itype = InteractionType::VSITE3FD;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(hoh));
    auto vsite3fd_name = potentialToParameterName(Potential::VSITE3FD);
    auto param_hoh = fs.findParameterType(hoh, vsite3fd_name[vsite3fdB]);
    EXPECT_TRUE(param_hoh->internalValue() == 0.15);
    EXPECT_FALSE(fs.parameterExists(coh));
}



TEST(VSite3OUT, hohCanSwapNo)
{
    std::string forcefield("ACS-pg-vs3");
    // Get the forcefield
    auto pd = getForceField(forcefield);

    auto hoh = Identifier({ "h_b", "o3_b", "h_b", "v3out" }, { 1, 1, 9 }, CanSwap::No);
    auto coh = Identifier({ "c3_b", "o3_b", "h_b", "v3out" }, { 1, 1, 9 }, CanSwap::No);
    // Compare identifiers
    EXPECT_FALSE(hoh == coh);
    auto itype = InteractionType::VSITE3OUT;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(hoh));
    auto vsite3out_name = potentialToParameterName(Potential::VSITE3OUT);
    auto param_hoh = fs.findParameterType(hoh, vsite3out_name[vsite3outA]);
    EXPECT_TRUE(param_hoh->internalValue() == -0.35);
    EXPECT_FALSE(fs.parameterExists(coh));
}

TEST(VSite4, nh3CanSwapNo)
{
    std::string forcefield("ACS-pg-vs3");
    // Get the forcefield
    auto pd = getForceField(forcefield);

    // 9 means v-site
    auto nh3 = Identifier({ "n3_b", "h_b", "h_b", "c3_b", "v4" }, { 1, 1, 1, 9 }, CanSwap::No);
    auto cnh = Identifier({ "n3_b", "h_b", "h_b", "h_b", "v4" }, { 1, 1, 1, 9 }, CanSwap::No);
    // Compare identifiers
    EXPECT_FALSE(nh3 == cnh);
    auto itype = InteractionType::VSITE4;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(nh3));
    auto vsite4_name = potentialToParameterName(Potential::VSITE4);
    EXPECT_TRUE(fs.findParameterType(nh3, vsite4_name[vsite4A])->internalValue() == -0.2);
    EXPECT_TRUE(fs.findParameterType(nh3, vsite4_name[vsite4B])->internalValue() == -0.3);
    EXPECT_TRUE(fs.findParameterType(nh3, vsite4_name[vsite4C])->internalValue() == -0.4);
    EXPECT_FALSE(fs.parameterExists(cnh));
}

TEST(VSite4S, nh3CanSwapNo)
{
    std::string forcefield("ACS-pg-vs3");
    // Get the forcefield
    auto pd = getForceField(forcefield);

    auto nh3 = Identifier({ "n3_b", "c3_b", "h_b", "h_b", "v4" }, { 1, 1, 1, 9 }, CanSwap::No);
    // Compare identifiers
    auto itype = InteractionType::VSITE4S;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(nh3));
    auto vsite4_name = potentialToParameterName(Potential::VSITE4S);
    EXPECT_TRUE(fs.findParameterType(nh3, vsite4_name[vsite4sA])->internalValue() == -0.3);
    EXPECT_TRUE(fs.findParameterType(nh3, vsite4_name[vsite4sB])->internalValue() == -0.5);
}

TEST(VSite4S3, nh3CanSwapNo)
{
    std::string forcefield("ACS-pg-vs3");
    // Get the forcefield
    auto pd = getForceField(forcefield);

    auto nh3 = Identifier({ "n3_b", "h_b", "h_b", "h_b", "v4" }, { 1, 1, 1, 9 }, CanSwap::No);
    // Compare identifiers
    auto itype = InteractionType::VSITE4S3;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(nh3));
    auto vsite4_name = potentialToParameterName(Potential::VSITE4S3);
    EXPECT_TRUE(fs.findParameterType(nh3, vsite4_name[vsite4s3A])->internalValue() == -0.2);
}

}  // namespace

}  // namespace alexandria
