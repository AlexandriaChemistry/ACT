/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2026
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

#include <set>

#include <gtest/gtest.h>

#include "act/import/import.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forces/forcecomputer.h"
#include "act/topology/actmol.h"
#include "act/topology/topology.h"

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
protected:
    bool setup(const std::string &ff,
               const std::string &molName,
               ACTMol            *mol)
    {
        std::string fileName = gmx::formatString("%s.sdf", molName.c_str());
        return setupFile(ff, gmx::test::TestFileManager::getInputFilePath(fileName), mol);
    }

    //! Read a molecule from a file and generate the topology
    bool setupFile(const std::string &ff,
                   const std::string &dataName,
                   ACTMol            *mol)
    {
        // Get forcefield
        auto pd  = getForceField(ff);
        
        std::vector<alexandria::MolProp> molprops;
        bool   userqtot   = false;
        double qtot_babel = 0;
        const char           *conf      = (char *)"minimum";
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        importFile(&msghandler, pd, dataName.c_str(), &molprops,
                   conf, JobType::OPT, userqtot ,&qtot_babel,
                   true);
        EXPECT_TRUE(msghandler.ok());
        EXPECT_TRUE(molprops.size() == 1);
        if (molprops.size() != 1)
        {
            return false;
        }
        auto &molprop = molprops[0];
        alexandria::ACTMol               mp_;

        mol->Merge(&molprop);
        // Uncomment in case of issues
        // msghandler.setPrintLevel(ACTStatus::Debug);
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
            vh.constructPositions(mp_.topology(), &coords);

            gmx::test::TestReferenceChecker checker_(this->rootChecker());
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
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
            vh.constructPositions(mp_.topology(), &coords);
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
            vh.distributeForces(mp_.topology(), coords, &forces);
            gmx::test::TestReferenceChecker checker_(this->rootChecker());
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
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

    /*! \brief Check the distribution of forces on virtual sites by finite differences
     *
     * A constant force F_v on each virtual site corresponds to the energy
     * E = -sum_v F_v . x_v(x), where the virtual site positions x_v depend on
     * the positions x of the atoms. The force on atom a is then
     * -dE/dx_a = sum_v F_v . dx_v/dx_a, which is computed numerically.
     * \param[in] ff       The force field
     * \param[in] dataName The molecule file
     * \param[in] vsites   The virtual site types that should be present
     */
    void testSpread(const std::string                 &ff,
                    const std::string                 &dataName,
                    const std::set<InteractionType>   &vsites)
    {
        ACTMol mol;
        ASSERT_TRUE(setupFile(ff, dataName, &mol));
        auto top = mol.topology();
        // The virtual site types in the topology
        std::set<InteractionType> found;
        for (const auto &entry : top->entries())
        {
            if (isVsite(entry.first) && !entry.second.empty())
            {
                found.insert(entry.first);
            }
        }
        auto names = [](const std::set<InteractionType> &its)
        {
            std::string str;
            for (auto it : its)
            {
                str += interactionTypeToString(it) + " ";
            }
            return str;
        };
        EXPECT_EQ(names(vsites), names(found));
        VsiteHandler vh;
        std::vector<gmx::RVec> coords = mol.xOriginal();
        vh.constructPositions(top, &coords);
        // A different force on each virtual site
        const auto &atoms = top->atoms();
        std::vector<gmx::RVec> fvsite(coords.size(), { 0, 0, 0 });
        int nvsite = 0;
        for (size_t i = 0; i < atoms.size(); i++)
        {
            if (atoms[i].pType() == ActParticle::Vsite)
            {
                nvsite += 1;
                fvsite[i] = { 1.0*nvsite, -0.5, 0.25*nvsite };
            }
        }
        ASSERT_GT(nvsite, 0);
        std::vector<gmx::RVec> forces = fvsite;
        vh.distributeForces(top, coords, &forces);
        // Energy as a function of the coordinates
        auto energy = [&](std::vector<gmx::RVec> x)
        {
            vh.constructPositions(top, &x);
            double e = 0;
            for (size_t i = 0; i < x.size(); i++)
            {
                e -= iprod(fvsite[i], x[i]);
            }
            return e;
        };
        double h = 1e-6;
        for (size_t i = 0; i < atoms.size(); i++)
        {
            for (int m = 0; m < DIM; m++)
            {
                if (atoms[i].pType() == ActParticle::Vsite)
                {
                    // The force has been moved to the atoms
                    EXPECT_EQ(0, forces[i][m]) << "vsite " << i;
                    continue;
                }
                auto xp = coords;
                auto xm = coords;
                xp[i][m] += h;
                xm[i][m] -= h;
                double fnum = -(energy(xp) - energy(xm))/(2*h);
                EXPECT_NEAR(fnum, forces[i][m], 1e-5) << "atom " << i << " " << atoms[i].name() << " dim " << m;
            }
        }
    }
};

TEST_F (VsiteHandlerTest, SpreadVsite2HF)
{
    testSpread("ACS-pg-vs2", gmx::test::TestFileManager::getInputFilePath("hydrogen-fluoride.sdf"),
               { InteractionType::VSITE1, InteractionType::VSITE2 });
}

TEST_F (VsiteHandlerTest, SpreadVsite3Water)
{
    testSpread("ACS-pg-vs3", gmx::test::TestFileManager::getInputFilePath("water.sdf"),
               { InteractionType::VSITE3, InteractionType::VSITE3OUT });
}

TEST_F (VsiteHandlerTest, SpreadVsite3SWater)
{
    testSpread("water_hyperpol", gmx::test::TestFileManager::getInputFilePath("water.sdf"),
               { InteractionType::VSITE3S });
}

TEST_F (VsiteHandlerTest, SpreadVsite3FDWater)
{
    testSpread("ACS-v3fds", gmx::test::TestFileManager::getInputFilePath("water.sdf"),
               { InteractionType::VSITE3FD });
}

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

/*! \brief A topology with one virtual site of each type that can be constructed
 * Five atoms and the virtual sites, each constructed from some of the
 * atoms. VSITE3FAD is left out, since it cannot be constructed yet.
 * \return the topology
 */
static Topology allVsiteTypes()
{
    Topology top;
    const int natom = 5;
    for(int i = 0; i < natom; i++)
    {
        top.addAtom(ActAtom("A", "C", "c", ActParticle::Atom, 6, 12.0, 0.0, 1));
    }
    struct Type
    {
        InteractionType     itype;
        int                 nparent;
        std::vector<double> params;
        int                 sign;
    };
    std::vector<Type> types = {
        { InteractionType::VSITE1,     1, { 1.0 },             0 },
        { InteractionType::VSITE2,     2, { 0.3 },             0 },
        { InteractionType::VSITE2FD,   2, { 0.05 },            0 },
        { InteractionType::VSITE3,     3, { 0.2, 0.3 },        0 },
        { InteractionType::VSITE3S,    3, { 0.25 },            0 },
        { InteractionType::VSITE3FD,   3, { 0.4, 0.07 },       0 },
        { InteractionType::VSITE3OUT,  3, { 0.1, 0.2, 2.0 },   1 },
        { InteractionType::VSITE3OUTS, 3, { 0.15, -1.5 },     -1 },
        { InteractionType::VSITE4,     4, { 0.1, 0.2, 0.3 },   0 },
        { InteractionType::VSITE4S,    4, { 0.15, 0.25 },      0 },
        { InteractionType::VSITE4S3,   4, { 0.2 },             0 }
    };
    int vs = natom;
    for(const auto &t : types)
    {
        top.addAtom(ActAtom("V", "X", "v", ActParticle::Vsite, 0, 0.0, 0.0, 1));
        // The parents in an order that differs per type
        std::vector<int> parents;
        for(int k = 0; k < t.nparent; k++)
        {
            parents.push_back((vs + 2*k) % natom);
        }
        TopologyEntryVector entries;
        if (t.sign != 0)
        {
            Vsite3OUT v3(parents[0], parents[1], parents[2], vs, t.sign);
            v3.setParams(t.params);
            entries.push_back(BaseContainer<TopologyEntry>(std::move(v3)));
        }
        else
        {
            TopologyEntry te;
            for(auto p : parents)
            {
                te.addAtom(p);
            }
            te.addAtom(vs);
            te.setParams(t.params);
            entries.push_back(BaseContainer<TopologyEntry>(std::move(te)));
        }
        top.addEntry(t.itype, entries);
        vs += 1;
    }
    return top;
}

// Test: The virtual sites of many conformations at once, in the layout of
// ConformationArrays, are the same as one conformation at a time, for all
// types, such that optimizing one of the routines does not change results
TEST(VsiteHandlerArrays, AllTypes)
{
    auto   top    = allVsiteTypes();
    size_t natom  = top.atoms().size();
    // More than four times the widest SIMD, and padding slots
    size_t nconf  = 37;
    size_t stride = 40;
    std::vector<std::vector<gmx::RVec>> x(nconf), f(nconf);
    unsigned int seed = 1993;
    auto random = [&seed]()
    {
        seed = 1103515245*seed + 12345;
        return ((seed >> 8) % 2001)/1000.0 - 1.0;
    };
    for(size_t c = 0; c < nconf; c++)
    {
        x[c].resize(natom);
        f[c].resize(natom);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                // Atoms about 0.1 nm apart, the virtual sites are constructed
                x[c][a][m] = 0.1*random() + (m == 0 ? 0.1*a : 0.0);
                f[c][a][m] = 100*random();
            }
        }
    }
    std::vector<double> xa(DIM*natom*stride), fa(DIM*natom*stride);
    for(size_t c = 0; c < stride; c++)
    {
        size_t cc = std::min(c, nconf-1);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                xa[(DIM*a+m)*stride + c] = x[cc][a][m];
                fa[(DIM*a+m)*stride + c] = f[cc][a][m];
            }
        }
    }
    // All types are there
    size_t ntype = 0;
    for(const auto &entry : top.entries())
    {
        ntype += isVsite(entry.first) ? 1 : 0;
    }
    EXPECT_EQ(11u, ntype);
    auto         xstart = x;
    VsiteHandler vsiteHandler;
    vsiteHandler.constructPositions(&top, xa.data(), stride, stride, nconf);
    vsiteHandler.distributeForces(&top, xa.data(), fa.data(), stride, stride, nconf);
    for(size_t c = 0; c < nconf; c++)
    {
        vsiteHandler.constructPositions(&top, &x[c]);
        vsiteHandler.distributeForces(&top, x[c], &f[c]);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                EXPECT_EQ(x[c][a][m], xa[(DIM*a+m)*stride + c]) << "conformation " << c << " particle " << a;
                EXPECT_EQ(f[c][a][m], fa[(DIM*a+m)*stride + c]) << "conformation " << c << " particle " << a;
                if (top.atoms()[a].pType() == ActParticle::Vsite)
                {
                    EXPECT_EQ(0, f[c][a][m]) << "the force on virtual site " << a << " is spread";
                    EXPECT_NE(xstart[c][a][m], x[c][a][m]) << "virtual site " << a << " is constructed";
                }
            }
        }
    }
}

}  // namespace

}  // namespace alexandria
