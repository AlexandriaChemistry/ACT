/*
 * This source file is part of the Alexandria program.
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
class AcmTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
        alexandria::ACTMol               mp_;

        //init set tolerance
        AcmTest () : checker_(this->rootChecker())
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
                     bool                             qSymm,
                     std::vector<double>              qtotal,
                     const std::vector<double>       &qcustom,
                     bool                             useHF)
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
            auto pd  = getForceField(model);

            double qtot_babel = myqtot;
            matrix box;
            EXPECT_TRUE(readBabel(pd, dataName.c_str(), &molprops,
                                  molname.c_str(), molname.c_str(),
                                  conf, &method, &basis, maxpot,
                                  nsymm, jobtype, &qtot_babel, false, box));

            if (trustObCharge)
            {
                EXPECT_TRUE(myqtot == qtot_babel);
            }
            auto &molprop = molprops[0];
            
            double qtot_sum = 0;
            auto fptr = molprop.fragmentPtr();
            EXPECT_TRUE(fptr->size() == qtotal.size());
            for(size_t i = 0; i < qtotal.size(); i++)
            {
                qtot_sum += qtotal[i];
                (*fptr)[i].setCharge(qtotal[i]);
            }
            mp_.Merge(&molprop);
            // Generate charges and topology
            auto imm = mp_.GenerateTopology(stdout, pd,
                                            missingParameters::Error);
            if (immStatus::OK != imm)
            {
                fprintf(stderr, "Error generating topology: %s\n", immsg(imm));
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
            mp_.symmetrizeCharges(pd, qSymm, nullptr);
            mp_.GenerateCharges(pd, forceComp, alg, qType::Calc, qcustom, &coords, &forces);
            
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
            auto fh        = mp_.fragmentHandler();
            if (fh)
            {
                auto atomStart = fh->atomStart();
                if (atomStart.size() >= 2)
                {
                    for(size_t f = 0; f < atomStart.size(); f++)
                    {
                        auto   natom = fh->topologies()[f]->atoms().size();
                        double qt    = 0;
                        for(size_t atom = atomStart[f]; atom < atomStart[f]+natom; atom++)
                        {
                            qt += myatoms[atom].charge();
                        }
                        auto label = gmx::formatString("Molecule %zu charge", f+1);
                        checker_.checkReal(qt, label.c_str());
                    }
                }
            }
        }

        static void TearDownTestCase()
        {
        }

};

TEST_F (AcmTest, WaterDimerACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::PDB, "water_dimer", true, {0,0}, qcustom, false);
}

TEST_F (AcmTest, WaterIodideACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::LOG, "water_I", true, {0,-1}, qcustom, true);
}

TEST_F (AcmTest, AXpgLOG)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::LOG, "1-butanol", true, {0}, qcustom, false);
}

TEST_F (AcmTest, AXpgZIP)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::ZIP, "trimethylphosphate3-esp.log.gz", true, {0}, qcustom, false);
}

TEST_F (AcmTest, AXpgPDB)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::PDB, "1-butanol", true, {0}, qcustom, false);
}

TEST_F (AcmTest, AXpgNoSymmLOG)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::LOG, "1-butanol", false, {0}, qcustom, false);
}

TEST_F (AcmTest, AXpgNoSymmPDB)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::PDB, "1-butanol", false, {0}, qcustom, false);
}

TEST_F (AcmTest, AXgLOG)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::LOG, "1-butanol", true, {0}, qcustom, false);
}

TEST_F (AcmTest, AXgPDB)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::PDB, "1-butanol", true, {0}, qcustom, false);
}

TEST_F (AcmTest, AXgNoSymmLOG)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::LOG, "1-butanol", false, {0}, qcustom, false);
}

TEST_F (AcmTest, AXgNoSymmPDB)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::PDB, "1-butanol", false, {0}, qcustom, false);
}

TEST_F (AcmTest, AXgNegative)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::LOG, "acetate", false, {-1}, qcustom, false);
}

TEST_F (AcmTest, AXpgNegative)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::LOG, "acetate", true, {-1}, qcustom, false);
}

TEST_F (AcmTest, AXgPositive)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::SDF, "guanidinium", false, {1}, qcustom, false);
}

TEST_F (AcmTest, AXpgPositive)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::SDF, "guanidinium", true, {1}, qcustom, false);
}

TEST_F (AcmTest, SQEgNeutral)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::PDB, "1-butanol", false, {0}, qcustom, false);
}

TEST_F (AcmTest, SQEgNegative)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::LOG, "acetate", false, {-1}, qcustom, false);
}

TEST_F (AcmTest, SQEgPositive)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::SDF, "guanidinium", false, {1}, qcustom, false);
}

TEST_F (AcmTest, SQEpgNeutral)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::PDB, "1-butanol", false, {0}, qcustom, false);
}

TEST_F (AcmTest, SQEpgNegative)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::LOG, "acetate", true, {-1}, qcustom, false);
}

TEST_F (AcmTest, SQEpgPositive)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::SDF, "guanidinium", true, {1}, qcustom, false);
}

TEST_F (AcmTest, WaterDimerACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::PDB, "water_dimer", true, {0,0}, qcustom, false);
}

TEST_F (AcmTest, WaterIodideACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::LOG, "water_I", true, {0, -1}, qcustom, true);
}

TEST_F (AcmTest, CustomButanolPDB)
{
    std::vector<double> qcustom = { -0.18, 0.06, 0.06, 0.06, -0.12, 0.06, 0.06,  -0.12, 0.06, 0.06, -0.12, 0.06, 0.06, -0.4, 0.4};
    testAcm("ACM-g", inputFormat::PDB, "1-butanol", false, {0}, qcustom, false);
}

TEST_F (AcmTest, CustomAcetatePDB)
{
    std::vector<double> qcustom = { -0.18, 0.06, 0.06, 0.06, 0.1, -0.55, -0.55 };
    testAcm("ACM-g", inputFormat::LOG, "acetate", false, {-1}, qcustom, false);
}

TEST_F (AcmTest, CustomGuanidiniumSDF)
{
    std::vector<double> qcustom = { -0.3, -0.3, -0.3, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3 };
    testAcm("ACM-g", inputFormat::SDF, "guanidinium", false, {1}, qcustom, false);
}

TEST_F (AcmTest, MethanolWaterACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::SDF, "methanol-water", true, {0,0}, qcustom, true);
}

TEST_F (AcmTest, MethanolWaterACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::SDF, "methanol-water", true, {0,0}, qcustom, true);
}

TEST_F (AcmTest, AcetateWaterACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::SDF, "acetate-water", true, {-1,0}, qcustom, true);
}

TEST_F (AcmTest, AcetateWaterACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::SDF, "acetate-water", true, {-1,0}, qcustom, true);
}

}

}
