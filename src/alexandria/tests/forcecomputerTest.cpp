/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
    gmx::MDLogger                   mdlog_ {};
 
    ForceComputerTest( ) :
        checker_(refData_.rootChecker())
    {
        gmx::test::FloatingPointTolerance tolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-6));
        checker_.setDefaultTolerance(tolerance);
    }
    
    void initMyMol(const char    *molname, 
                   const Poldata *pd,
                   ForceComputer *fcomp,
                   t_inputrec    *inputrec,
                   MyMol         *mp)
    {
        int           maxpot   = 100;
        int           nsymm    = 0;
        const char   *conf     = (char *)"minimum";
        std::string   method, basis;
        const char   *jobtype  = (char *)"Opt";
        
        std::string   dataName = gmx::test::TestFileManager::getInputFilePath(molname);
        auto          molprop  = new alexandria::MolProp;
        double        qtot     = 0;
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
        mp->Merge(molprop);
        // Generate charges and topology
        fill_inputrec(inputrec);
        mp->setInputrec(inputrec);
        
        auto imm = mp->GenerateTopology(stdout, pd,
                                        missingParameters::Error, true);
        EXPECT_TRUE(immStatus::OK == imm);
        
        // Needed for GenerateCharges
        CommunicationRecord cr;
        auto alg = ChargeGenerationAlgorithm::NONE;
        std::vector<double> qcustom;
        bool qSymm = false;
        mp->symmetrizeCharges(pd, qSymm, nullptr);
        std::vector<gmx::RVec> forces(mp->atomsConst().size());
        mp->GenerateCharges(pd, fcomp, mdlog_, &cr, alg, qcustom, &forces);
    }
    
    void test(const char *molname, const char *forcefield, 
              bool testPolarizability, double stretch = 1)
    {
        // Get poldata
        auto pd  = getPoldata(forcefield);
        
        double rmsToler = 0.00001;
        auto fcomp = new ForceComputer(pd, rmsToler, 25);
        
        t_inputrec      inputrecInstance;
        
        // The molecule
        MyMol mp_;
        initMyMol(molname, pd, fcomp, &inputrecInstance, &mp_);
   
        
        auto atoms = mp_.atomsConst();
        std::vector<gmx::RVec>            forces, coordinates;
        std::map<InteractionType, double> energies;
        for(size_t i = 0; i < atoms.size(); i++)
        {
            rvec xxx;
            for(int m = 0; m < DIM; m++)
            {
                xxx[m] = mp_.x()[i][m]*stretch;
            }
            coordinates.push_back(xxx);
            forces.push_back({ 0, 0, 0 });
        }
        if (testPolarizability)
        {
            auto qCalc = mp_.qTypeProps(qType::Calc);
            fcomp->calcPolarizability(mp_.topology(), &coordinates, qCalc);
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
            if (stretch != 1)
            {
                mp_.setX(coordinates);
            }
            // This turn comparison of gromacs and ACT on. For debugging
            // you may want to set this to false.
            bool       strict    = true;
            double     shellRmsf;
            t_commrec *crtmp     = init_commrec();
            crtmp->nnodes = 1;
            PaddedVector<gmx::RVec> gmxforces;
            gmxforces.resizeWithPadding(mp_.atomsConst().size());
            mp_.calculateEnergyOld(crtmp, &gmxforces, &shellRmsf);
            auto ed = mp_.enerdata();
            for(int i = 0; i < F_NRE; i++)
            {
                if (ed->term[i] != 0)
                {
                    std::string label = gmx::formatString("%s", interaction_function[i].name);
                    if (stretch != 1)
                    {
                        label += gmx::formatString("%g", stretch);
                    }
                    label += "_gmx";
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
                    std::string label(interaction_function[ftype].name);
                    if (stretch != 1)
                    {
                        label += gmx::formatString("%g", stretch);
                    }
                    label += "_act";
                    checker_.checkReal(ener.second, label.c_str());
                    if (strict)
                    {
                        EXPECT_TRUE(std::abs(ener.second-ed->term[ftype]) < 1e-3);
                    }
                }
            }
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
};

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

TEST_F (ForceComputerTest, UracilPolarizability)
{
    test("uracil.sdf", "ACS-pg", true);
}

}  // namespace

}  // namespace gmx
