/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements test of bonded force routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_listed-forces
 */
#include "actpre.h"

#include "act/forces/forcecomputer.h"

#include <cmath>

#include <memory>

#include <gtest/gtest.h>

#include "act/poldata/poldata_utils.h"
#include "alexandria/atype_mapping.h"
#include "alexandria/babel_io.h"
#include "alexandria/fill_inputrec.h"
#include "alexandria/mymol.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/utility/fatalerror.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class ForceComputerTest : public ::testing::Test
{
protected:
    gmx::test::TestReferenceData    refData_;
    gmx::test::TestReferenceChecker checker_;
    
    ForceComputerTest( ) :
        checker_(refData_.rootChecker())
    {
        gmx::test::FloatingPointTolerance tolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-6));
        checker_.setDefaultTolerance(tolerance);
    }

    void test(const char *molname, const char *forcefield)
    {
        int                             maxpot   = 100;
        int                             nsymm    = 0;
        const char                     *conf     = (char *)"minimum";
        std::string                     method, basis;
        const char                     *jobtype  = (char *)"Opt";
        
        std::string                     dataName;
        auto molprop = new alexandria::MolProp;
        
        dataName = gmx::test::TestFileManager::getInputFilePath(molname);
        double qtot = 0;
        bool readOK = readBabel(dataName.c_str(), molprop, molname, molname,
                                conf, &method, &basis,
                                maxpot, nsymm, jobtype, &qtot, false);
        EXPECT_TRUE(readOK);
        if (readOK)
        {
            std::map<std::string, std::string> g2a;
            gaffToAlexandria("", &g2a);
            if (!g2a.empty())
            {
                EXPECT_TRUE(renameAtomTypes(molprop, g2a));
            }
        }
        MyMol mp_;
        mp_.Merge(molprop);
        // Generate charges and topology
        t_inputrec      inputrecInstance;
        t_inputrec     *inputrec   = &inputrecInstance;
        fill_inputrec(inputrec);
        mp_.setInputrec(inputrec);
        
        // Get poldata
        auto pd  = getPoldata(forcefield);
        auto imm = mp_.GenerateTopology(stdout, pd,
                                        missingParameters::Error);
        EXPECT_TRUE(immStatus::OK == imm);
        if (immStatus::OK != imm)
        {
            fprintf(stderr, "Could not generate topology because '%s'. Used basis %s and method %s.\n",
                    immsg(imm), basis.c_str(), method.c_str());
            return;
        }
        // Needed for GenerateCharges
        CommunicationRecord cr;
        gmx::MDLogger  mdlog {};
        auto fcomp = new ForceComputer(pd);
        auto alg = ChargeGenerationAlgorithm::NONE;
        std::vector<double> qcustom;
        bool qSymm = false;
        mp_.symmetrizeCharges(pd, qSymm, nullptr);
        mp_.GenerateCharges(pd, fcomp, mdlog, &cr, alg, qcustom);
        // For debugging: print all the interactions in the topology
        // gmx_init_debug(1, gmx::formatString("%s-%s.debug", molname, forcefield).c_str());
        // mp_.topology()->dump(stdout);
        
        auto atoms = mp_.atoms();
        std::vector<gmx::RVec>            forces, coordinates;
        std::map<InteractionType, double> energies;
        for(int i = 0; i < atoms->nr; i++)
        {
            coordinates.push_back(mp_.x()[i]);
            forces.push_back({ 0, 0, 0 });
        }
        double shellRmsf;
        t_commrec *crtmp = init_commrec();
        crtmp->nnodes = 1;
        mp_.calculateEnergyOld(crtmp, &shellRmsf);
        auto ed = mp_.enerdata();
        for(int i = 0; i < F_NRE; i++)
        {
            if (ed->term[i] != 0)
            {
                std::string label = gmx::formatString("%s_gmx", interaction_function[i].name);
                checker_.checkReal(ed->term[i], label.c_str());
            }
        }
        fcomp->compute(mp_.topology(), &coordinates, &forces, &energies);
        
        for(auto &ener: energies)
        {
            if (ener.second != 0)
            {
                // TODO remove this hack and make a real interactiontypetoftype.
                int ftype = F_EPOT;
                if (pd->interactionPresent(ener.first))
                {
                    ftype = pd->findForcesConst(ener.first).fType();
                }
                checker_.checkReal(ener.second, interaction_function[ftype].name);
                EXPECT_TRUE(std::abs(ener.second-ed->term[ftype]) < 1e-3);
            }
        }
        const char *xyz[DIM] = { "X", "Y", "Z" };
        auto mpf = mp_.f();
        for(size_t i = 0; i < forces.size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                auto label = gmx::formatString("%s-%zu_gmx f%s", 
                                               *atoms->atomtype[i],
                                               i+1, xyz[m]);
                checker_.checkReal(mpf[i][m], label.c_str());
            }
        }
        for(size_t i = 0; i < forces.size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                auto label = gmx::formatString("%s-%zu f%s", 
                                               *atoms->atomtype[i],
                                               i+1, xyz[m]);
                checker_.checkReal(forces[i][m], label.c_str());
                EXPECT_TRUE(std::abs(forces[i][m]-mpf[i][m]) < 1e-3);
            }
        }
        // gmx_stop_debug();
    }
};

TEST_F (ForceComputerTest, CarbonDioxide)
{
    test("carbon-dioxide.sdf", "ACS-g");
}

TEST_F (ForceComputerTest, HydrogenChloride)
{

    test("hydrogen-chloride.sdf", "ACS-g");
}

TEST_F (ForceComputerTest, Water)
{

    test("water-3-oep.log.pdb", "ACS-g");
}

TEST_F (ForceComputerTest, Acetone)
{
    test("acetone-3-oep.log.pdb", "ACS-g");
}

TEST_F (ForceComputerTest, Uracil)
{
    test("uracil.sdf", "ACS-g");
}

TEST_F (ForceComputerTest, CarbonDioxidePol)
{
  test("carbon-dioxide.sdf", "ACS-pg");
}

TEST_F (ForceComputerTest, HydrogenChloridePol)
{

    test("hydrogen-chloride.sdf", "ACS-pg");
}

TEST_F (ForceComputerTest, WaterPol)
{

    test("water-3-oep.log.pdb", "ACS-pg");
}

TEST_F (ForceComputerTest, AcetonePol)
{
    test("acetone-3-oep.log.pdb", "ACS-pg");
}

TEST_F (ForceComputerTest, UracilPol)
{
    test("uracil.sdf", "ACS-pg");
}

}  // namespace

}  // namespace gmx
