/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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
 * \ingroup module_listed-forces
 */
#include "actpre.h"

#include "act/forces/forcecomputer.h"

#include <cmath>

#include <memory>

#include <gtest/gtest.h>

#include "actmol_util.h"

#include "act/alexandria/atype_mapping.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/actmol.h"
#include "act/forcefield/forcefield_utils.h"
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
    
    void test(const char *molname, const char *forcefield, 
              bool testPolarizability, double stretch = 1)
    {
        // Get forcefield
        auto pd  = getForceField(forcefield);
        
        double rmsToler = 0.0000001;
        auto fcomp = new ForceComputer(rmsToler, 25);
        
        t_inputrec      inputrecInstance;
        
        // The molecule
        std::vector<ACTMol> mps;
        initACTMol(molname, pd, fcomp, &inputrecInstance, &mps);
   
        for(auto &mp : mps)
        {
            std::vector<gmx::RVec> coordinates, forces;
            auto xo = mp.xOriginal();
            for(size_t i = 0; i < xo.size(); i++)
            {
                rvec xxx;
                for(int m = 0; m < DIM; m++)
                {
                    xxx[m] = xo[i][m]*stretch;
                }
                coordinates.push_back(xxx);
                forces.push_back({ 0, 0, 0 });
            }
            if (testPolarizability)
            {
                auto qCalc = mp.qTypeProps(qType::Calc);
                qCalc->initializeMoments();
                fcomp->calcPolarizability(pd, mp.topology(), &coordinates, qCalc);
                auto alpha = qCalc->polarizabilityTensor();
                const char *xyz[DIM] = { "X", "Y", "Z" };
                
                for(int m = 0; m < DIM; m++)
                {
                    for(int n = 0; n < DIM; n++)
                    {
                        std::string label = gmx::formatString("alpha[%s][%s]", xyz[m], xyz[n]);
                        checker_.checkReal(convertFromGromacs(alpha[m][n], "A^3"),
                                           label.c_str());
                    }
                }
            }
            else
            {
                std::map<InteractionType, double> gmxEnergies, actEnergies;
                // This turns on comparison of gromacs and ACT. For debugging
                // you may want to set this to false.
                bool       strict    = true;
                double     shellRmsf;
                t_commrec *crtmp     = init_commrec();
                crtmp->nnodes = 1;
                PaddedVector<gmx::RVec> gmxforces;
                gmxforces.resizeWithPadding(mp.atomsConst().size());
                
                auto fsc = pd->forcesConst();
                mp.calculateEnergyOld(crtmp, &coordinates, &gmxforces, &gmxEnergies, &shellRmsf);
                fcomp->compute(pd, mp.topology(), &coordinates, &forces, &actEnergies);
                for(auto &ifm : gmxEnergies)
                {
                    int ftype = F_EPOT;
                    if (ifm.first != InteractionType::EPOT)
                    {
                        auto fs = fsc.find(ifm.first);
                        EXPECT_TRUE(fsc.end() != fs);
                        ftype = fs->second.gromacsType();
                    }
                    std::string label = gmx::formatString("%s", interaction_function[ftype].name);
                    if (stretch != 1)
                    {
                        label += gmx::formatString("%g", stretch);
                    }
                    auto gmxlabel = label+"_gmx";
                    checker_.checkReal(ifm.second, gmxlabel.c_str());
                    auto actlabel = label+"_act";
                    // Hack to compare GROMACS Buckingham or LJ to ACT
                    double actEner;
                    auto irep  = InteractionType::REPULSION;
                    auto idisp = InteractionType::DISPERSION;
                    if (ifm.first == InteractionType::VDW &&
                        actEnergies.find(irep) != actEnergies.end() &&
                        actEnergies.find(idisp) != actEnergies.end())
                    {
                        actEner  = (actEnergies[irep]+actEnergies[idisp]);
                    }
                    else
                    {
                        actEner  = actEnergies[ifm.first];
                    }
                    checker_.checkReal(actEner, actlabel.c_str());
                    if (strict)
                    {
                        EXPECT_TRUE(std::abs(ifm.second-actEner) < 1e-3);
                    }
                }
                auto atoms = mp.atomsConst();
                const char *xyz[DIM] = { "X", "Y", "Z" };
                for(size_t i = 0; i < forces.size(); i++)
                {
                    bool shell = atoms[i].pType() == eptShell;
                    for(int m = 0; m < DIM; m++)
                    {
                        if (!shell)
                        {
                            std::string stretchName;
                            if (stretch != 1)
                            {
                                stretchName += gmx::formatString("%g", stretch);
                            }
                            checker_.checkReal(gmxforces[i][m], gmx::formatString("%s-%zu_gmx%s f%s", 
                                                                                  atoms[i].ffType().c_str(),
                                                                                  i+1, stretchName.c_str(),
                                                                                  xyz[m]).c_str());
                            checker_.checkReal(forces[i][m], gmx::formatString("%s-%zu_act%s f%s", 
                                                                               atoms[i].ffType().c_str(),
                                                                               i+1, stretchName.c_str(),
                                                                               xyz[m]).c_str());
                        }
                        if (strict)
                        {
                            EXPECT_TRUE(std::abs(forces[i][m]-gmxforces[i][m]) < 1e-3);
                        }
                    }
                }
            }
        }
    }
    
    void testSimple(const char *molname, const char *forcefield)
    {
        // Get forcefield
        auto pd  = getForceField(forcefield);
        
        double rmsToler = 0.0000001;
        auto fcomp = new ForceComputer(rmsToler, 25);
        
        t_inputrec      inputrecInstance;
        
        // The molecule
        std::vector<ACTMol> mps;
        initACTMol(molname, pd, fcomp, &inputrecInstance, &mps);
   
        for(auto &mp : mps)
        {
            std::vector<gmx::RVec> coordinates, forces;
            auto xo = mp.xOriginal();
            for(size_t i = 0; i < xo.size(); i++)
            {
                coordinates.push_back(xo[i]);
                forces.push_back({ 0, 0, 0 });
            }
            std::map<InteractionType, double> actEnergies;
            auto fsc = pd->forcesConst();
            fcomp->compute(pd, mp.topology(), &coordinates, &forces, &actEnergies);
            for(auto &ifm : actEnergies)
            {
                int ftype = F_EPOT;
                switch (ifm.first)
                {
                case InteractionType::EPOT:
                    break;
                case InteractionType::DISPERSION:
                    ftype = F_DISPERSION;
                    break;
                case InteractionType::REPULSION:
                    ftype = F_REPULSION;
                    break;
                default:
                    {
                        auto fs = fsc.find(ifm.first);
                        EXPECT_TRUE(fsc.end() != fs);
                        ftype = fs->second.gromacsType();
                    }
                }
                std::string label = gmx::formatString("%s", interaction_function[ftype].name);
                auto actEner  = actEnergies[ifm.first];
                checker_.checkReal(actEner, label.c_str());
                auto atoms = mp.atomsConst();
                const char *xyz[DIM] = { "X", "Y", "Z" };
                for(size_t i = 0; i < forces.size(); i++)
                {
                    bool shell = atoms[i].pType() == eptShell;
                    for(int m = 0; m < DIM; m++)
                    {
                        if (!shell)
                        {
                            checker_.checkReal(forces[i][m], gmx::formatString("%s-%zu f%s", 
                                                                               atoms[i].ffType().c_str(),
                                                                               i+1, xyz[m]).c_str());
                        }
                    }
                }
            }
        }
    }
};

TEST_F (ForceComputerTest, MethaneThiol)
{
    test("methanethiol.sdf", "ACS-g", false);
}

TEST_F (ForceComputerTest, CarbonDioxide)
{
    test("carbon-dioxide.sdf", "ACS-g", false);
}

TEST_F (ForceComputerTest, HydrogenChloride)
{

    test("hydrogen-chloride.sdf", "ACS-g", false);
}

TEST_F (ForceComputerTest, HydrogenChlorideStretch)
{
    test("hydrogen-chloride.sdf", "ACS-g", false, 0.98);
    test("hydrogen-chloride.sdf", "ACS-g", false, 1);
    test("hydrogen-chloride.sdf", "ACS-g", false, 1.02);
}

TEST_F (ForceComputerTest, Water)
{
    test("water-3-oep.log.pdb", "ACS-g", false);
}

TEST_F (ForceComputerTest, WaterStretch)
{
    test("water-3-oep.log.pdb", "ACS-g", false, 0.98);
    test("water-3-oep.log.pdb", "ACS-g", false, 1.02);
}

TEST_F (ForceComputerTest, Acetone)
{
    test("acetone-3-oep.log.pdb", "ACS-g", false);
}

TEST_F (ForceComputerTest, AcetoneSdf)
{
    test("acetone.sdf", "ACS-g", false);
}

TEST_F (ForceComputerTest, AcetoneStretch)
{
    test("acetone-3-oep.log.pdb", "ACS-g", false, 0.98);
    test("acetone-3-oep.log.pdb", "ACS-g", false, 1.07);
}

TEST_F (ForceComputerTest, AcetoneNonPlanar)
{
    test("acetone-nonplanar.pdb", "ACS-g", false);
}

TEST_F (ForceComputerTest, Uracil)
{
    test("uracil.sdf", "ACS-g", false);
}

TEST_F (ForceComputerTest, UracilStretch)
{
    test("uracil.sdf", "ACS-g", false, 0.9);
    test("uracil.sdf", "ACS-g", false, 1.1);
}

TEST_F (ForceComputerTest, AcetoneDih)
{
    test("acetone-3-oep.log.pdb", "ACS-g-dih", false);
}

TEST_F (ForceComputerTest, UracilDih)
{
    test("uracil.sdf", "ACS-g-dih", false);
}

TEST_F (ForceComputerTest, WaterUB)
{

    test("water-3-oep.log.pdb", "ACS-g-dih", false);
}

TEST_F (ForceComputerTest, HydrogenChloridePolStretch)
{
    test("hydrogen-chloride.sdf", "ACS-pg", false, 0.98);
    test("hydrogen-chloride.sdf", "ACS-pg", false, 1.05);
}

TEST_F (ForceComputerTest, WaterPol)
{
    test("water-3-oep.log.pdb", "ACS-pg", false);
}

TEST_F (ForceComputerTest, WaterPolStretch)
{
    test("water-3-oep.log.pdb", "ACS-pg", false, 0.98);
    test("water-3-oep.log.pdb", "ACS-pg", false, 1.04);
}

TEST_F (ForceComputerTest, AcetonePol)
{
    test("acetone-3-oep.log.pdb", "ACS-pg", false);
}

TEST_F (ForceComputerTest, AcetonePolStretch)
{
    test("acetone-3-oep.log.pdb", "ACS-pg", false, 0.97);
    test("acetone-3-oep.log.pdb", "ACS-pg", false, 1.01);
}

TEST_F (ForceComputerTest, AcetoneNonPlanarPol)
{
    test("acetone-nonplanar.pdb", "ACS-pg", false);
}

TEST_F (ForceComputerTest, UracilPol)
{
    test("uracil.sdf", "ACS-pg", false);
}

TEST_F (ForceComputerTest, UracilPolStretch)
{
    test("uracil.sdf", "ACS-pg", false, 0.99);
    test("uracil.sdf", "ACS-pg", false, 1.06);
}

TEST_F (ForceComputerTest, CarbonDioxidePolarizability)
{
  test("carbon-dioxide.sdf", "ACS-pg", true);
}

TEST_F (ForceComputerTest, HydrogenChloridePolarizability)
{
    test("hydrogen-chloride.sdf", "ACS-pg", true);
}

TEST_F (ForceComputerTest, WaterPolarizability)
{
    test("water-3-oep.log.pdb", "ACS-pg", true);
}

TEST_F (ForceComputerTest, AcetonePolarizability)
{
    test("acetone-3-oep.log.pdb", "ACS-pg", true);
}

TEST_F (ForceComputerTest, AcetonePolGbham)
{
    testSimple("acetone-3-oep.log.pdb", "ACS-pg-gbham");
}

TEST_F (ForceComputerTest, AcetonePolarizabilityGbham)
{
    test("acetone-3-oep.log.pdb", "ACS-pg-gbham", true);
}

TEST_F (ForceComputerTest, UracilPolarizability)
{
    test("uracil.sdf", "ACS-pg", true);
}

}  // namespace

}  // namespace gmx
