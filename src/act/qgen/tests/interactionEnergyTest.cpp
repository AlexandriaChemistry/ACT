/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2025
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
#include <cmath>

#include <map>
#include <numeric>

#include <gtest/gtest.h>

#include "act/alexandria/atype_mapping.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/actmol.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/qgen/qgen_acm.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! Simple enum to distinguish file formats
enum class inputFormat {
    LOG,
    PDB,
    SDF,
    ZIP
};

//! Class to test the Alexandria Charge Model
class InteractionEnergyTest : public gmx::test::CommandLineTestBase
{
protected:
    gmx::test::TestReferenceChecker checker_;
    alexandria::ACTMol               mp_;

    //init set tolerance
    InteractionEnergyTest () : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-5);
        checker_.setDefaultTolerance(tolerance);
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    void testAcm(const std::string               &model, 
                 inputFormat                      inputformat, 
                 const std::string               &molname, 
                 std::vector<double>              qtotal,
                 const std::vector<double>       &qcustom,
                 bool                             useHF,
                 bool                             oneH = true)
    {
        int                   maxpot    = 100;
        int                   nsymm     = 0;
        const char           *conf      = (char *)"minimum";
        const char           *jobtype   = (char *)"Opt";
        std::string           method;
        std::string           basis;
        std::string           fileName(molname);
        std::vector<alexandria::MolProp> molprops;
        bool                  trustObCharge = false;
        
        if (inputformat == inputFormat::LOG)
        {
            if (useHF)
            {
                fileName.append(".log");
                method.assign("HF");
                basis.assign("3-21G");
            }
            else
            {
                fileName.append("-3-oep.log");
                method.assign("B3LYP");
                basis.assign("GEN");
            }
            trustObCharge = true;
        }
        else if (inputformat == inputFormat::ZIP)
        {
            method.assign("B3LYP");
            basis.assign("GEN");
            trustObCharge = true;
        }
        else if (inputformat == inputFormat::PDB)
        {
            fileName.append(".pdb");
        }
        else if (inputformat == inputFormat::SDF)
        {
            fileName.append(".sdf");
        }
        std::string dataName = gmx::test::TestFileManager::getInputFilePath(fileName);
        
        // Compute total charge from input
        int myqtot = 0;
        for(size_t i = 0; i < qtotal.size(); i++)
        {
            myqtot += qtotal[i];
        }
        // Get forcefield
        auto pd = getForceField(model);
        bool   userqtot   = !qcustom.empty();
        double qtot_babel = myqtot;
        matrix box;
        EXPECT_TRUE(readBabel(pd, dataName.c_str(), &molprops,
                              molname.c_str(), molname.c_str(),
                              conf, &method, &basis, maxpot, nsymm,
                              jobtype, userqtot, &qtot_babel, false, box, oneH));

        if (trustObCharge)
        {
            EXPECT_TRUE(myqtot == qtot_babel);
        }
        for(auto &molprop: molprops)
        {
            double qtot_sum = 0;
            auto fptr = molprop.fragmentPtr();
            EXPECT_TRUE(fptr->size() == qtotal.size());
            for(size_t i = 0; i < qtotal.size(); i++)
            {
                qtot_sum += qtotal[i];
                (*fptr)[i].setCharge(qtotal[i]);
            }
            mp_.Merge(&molprop);
            MsgHandler msghandler;
            msghandler.setPrintLevel(ACTStatus::Warning);
            // Generate charges and topology
            mp_.GenerateTopology(&msghandler, pd, missingParameters::Ignore);
            if (!msghandler.ok())
            {
                return;
            }
            
            // Needed for GenerateCharges
            auto forceComp = new ForceComputer();
            std::vector<gmx::RVec> forces(mp_.atomsConst().size());
            std::vector<gmx::RVec> coords = mp_.xOriginal();
            auto alg = ChargeGenerationAlgorithm::NONE;
            if (!qcustom.empty())
            {
                alg = ChargeGenerationAlgorithm::Custom;
            }
            mp_.GenerateCharges(&msghandler, pd, forceComp, alg, qType::Calc, qcustom, &coords, &forces);
            
            std::vector<double> qtotValues;
            auto myatoms = mp_.atomsConst();
            for (size_t atom = 0; atom < myatoms.size(); atom++)
            {
                qtotValues.push_back(myatoms[atom].charge());
            }
            double qtot_all = std::accumulate(qtotValues.begin(), qtotValues.end(), 0.0);
            EXPECT_TRUE(std::fabs(qtot_all - qtot_sum) < 1e-4);
            char   buf[256];
            snprintf(buf, sizeof(buf), "qtotValuesEqdAlgorithm_%s", 
                     chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str());
            checker_.checkInteger(static_cast<int>(qtotValues.size()), "qtotSize");
            checker_.checkSequence(qtotValues.begin(), qtotValues.end(), buf);
            // This vector has n+1 entries for n fragments
            auto fh        = mp_.fragmentHandler();
            if (fh)
            {
                auto atomStart = fh->atomStart();
                if (atomStart.size() >= 2)
                {
                    for(size_t f = 0; f < atomStart.size(); f++)
                    {
                        double qt    = 0;
                        auto   natom = fh->topologies()[f]->atoms().size();
                        for(size_t atom = atomStart[f]; atom < atomStart[f]+natom; atom++)
                        {
                            qt += myatoms[atom].charge();
                        }
                        auto label = gmx::formatString("Molecule %zu charge", f+1);
                        checker_.checkReal(qt, label.c_str());
                    }
                }
            }
            // Now the energies
            double rmsToler = 0.00001;
            auto fcomp = new ForceComputer(rmsToler, 25);
            if (mp_.fragmentHandler()->topologies().size() > 1)
            {
                std::vector<gmx::RVec> forces;
                std::vector<gmx::RVec> coords = mp_.xOriginal();
                std::map<InteractionType, double> einter;
                mp_.calculateInteractionEnergy(&msghandler, pd, fcomp, &einter, &forces, &coords, true);
                for(const auto &e : einter)
                {
                    checker_.checkReal(e.second, interactionTypeToString(e.first).c_str());
                }
                double esum = (einter[InteractionType::ELECTROSTATICS] + einter[InteractionType::INDUCTION] +
                               einter[InteractionType::INDUCTIONCORRECTION] +
                               einter[InteractionType::DISPERSION] + einter[InteractionType::EXCHANGE]);
                double toler = 1e-3;
                auto iepot = InteractionType::EPOT;
                EXPECT_TRUE(std::abs(esum - einter[iepot]) <= toler);
                // Check what happens when we sum induction correction into induction
                std::map<InteractionType, double> einter2;
                mp_.calculateInteractionEnergy(&msghandler, pd, fcomp, &einter2, &forces, &coords, true);
                EXPECT_TRUE(std::abs(einter2[iepot]-einter[iepot]) <= toler);
            }
        }
    }
            
    static void TearDownTestCase()
    {
    }
    
};

TEST_F (InteractionEnergyTest, WaterDimerACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::PDB, "water_dimer", {0,0}, qcustom, false);
}

TEST_F (InteractionEnergyTest, WaterIodideACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::LOG, "water_I", {0,-1}, qcustom, true);
}

TEST_F (InteractionEnergyTest, WaterDimerACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::PDB, "water_dimer", {0,0}, qcustom, false);
}

TEST_F (InteractionEnergyTest, WaterIodideACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::LOG, "water_I", {0, -1}, qcustom, true);
}

TEST_F (InteractionEnergyTest, MethanolWaterACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::SDF, "methanol-water", {0,0}, qcustom, true);
}

TEST_F (InteractionEnergyTest, MethanolWaterACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::SDF, "methanol-water", {0,0}, qcustom, true);
}

TEST_F (InteractionEnergyTest, AcetateWaterACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::SDF, "acetate-water", {-1,0}, qcustom, true);
}

TEST_F (InteractionEnergyTest, AcetateWaterACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::SDF, "acetate-water", {-1,0}, qcustom, true);
}

TEST_F (InteractionEnergyTest, HydrogenFluorideDimerACSps)
{
    std::vector<double> qcustom;
    testAcm("ACS-ps", inputFormat::SDF, "hfdimer", {0,0}, qcustom, true, false);
}

}

}
