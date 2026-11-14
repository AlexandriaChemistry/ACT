/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2026
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

#include "act/import/import.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/qgen/qgen_acm.h"
#include "act/topology/actmol.h"

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
    PDB,
    SDF,
    ZIP
};

//! Class to test the Alexandria Charge Model
class AcmTest : public gmx::test::CommandLineTestBase
{
    protected:
        //! For checking values computed in these tests
        gmx::test::TestReferenceChecker checker_;
        //! Our testing molecule
        alexandria::ACTMol               mp_;

        //! Init set tolerance
        AcmTest () : checker_(this->rootChecker())
        {
            const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 1e-5);
            checker_.setDefaultTolerance(tolerance);
        }

        //! Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        //! The actual testing routine
        void testAcm(const std::string               &model, 
                     inputFormat                      inputformat, 
                     const std::string               &molname, 
                     std::vector<double>              qtotal,
                     const std::vector<double>       &qcustom)
        {
            const char           *conf      = (char *)"minimum";
            std::string           fileName("molecules/" + molname);
            std::vector<alexandria::MolProp> molprops;
            
            if (inputformat == inputFormat::PDB)
            {
                fileName.append(".pdb");
            }
            else if (inputformat == inputFormat::SDF)
            {
                fileName.append(".sdf");
            }
            const std::string dataName = gmx::test::TestFileManager::getInputFilePath(fileName);
            
            // Compute total charge from input
            int myqtot = 0;
            for(size_t i = 0; i < qtotal.size(); i++)
            {
                myqtot += qtotal[i];
            }
            // Get forcefield
            auto pd           = getForceField(model);
            bool   userqtot   = !qcustom.empty();
            double qtot_babel = myqtot;
            MsgHandler msghandler;
            msghandler.setPrintLevel(ACTStatus::Warning);
            importFile(&msghandler, pd, dataName, &molprops,
                       conf, JobType::OPT, userqtot ,&qtot_babel,
                       true);
            EXPECT_TRUE(msghandler.ok());
            if (!msghandler.ok())
            {
                return;
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
            // Uncomment in case of issues
            // msghandler.setACTStatus(ACTStatus::Debug);
            // Generate charges and topology
            mp_.GenerateTopology(&msghandler, pd,
                                 missingParameters::Ignore);
            EXPECT_TRUE(msghandler.ok());
            if (!msghandler.ok())
            {
                return;
            }
            
            // Needed for GenerateCharges
            auto forceComp = new ForceComputer();
            std::vector<gmx::RVec> forces(mp_.atomsConst().size());
            std::vector<gmx::RVec> coords = mp_.xOriginal();
            auto alg = pd->chargeGenerationAlgorithm();
            if (!qcustom.empty())
            {
                alg = ChargeGenerationAlgorithm::Custom;
                mp_.setCharges(qcustom);
            }
            else
            {
                mp_.generateCharges(&msghandler, pd, forceComp, alg, &coords, &forces);
            }
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
                        auto   natom = fh->topologies()[f].atoms().size();
                        double qt    = 0;
                        for(size_t atom = atomStart[f]; atom < atomStart[f]+natom; atom++)
                        {
                            qt += myatoms[atom].charge();
                        }
                        const auto label = gmx::formatString("Molecule %zu charge", f+1);
                        checker_.checkReal(qt, label.c_str());
                    }
                }
            }
        }

        //! For cleaning
        static void TearDownTestCase()
        {
        }

};

TEST_F (AcmTest, WaterDimerACSg)
{
    const std::vector<double> qcustom;
    testAcm("ACS-g4", inputFormat::PDB, "water_dimer", {0,0}, qcustom);
}

TEST_F (AcmTest, AXpgPDB)
{
    const std::vector<double> qcustom;
    testAcm("ACM-pg2", inputFormat::PDB, "1-butanol", {0}, qcustom);
}

TEST_F (AcmTest, AXpgNoSymmPDB)
{
    const std::vector<double> qcustom;
    testAcm("ACM-pg2", inputFormat::PDB, "1-butanol", {0}, qcustom);
}

TEST_F (AcmTest, AXgPDB)
{
    const std::vector<double> qcustom;
    testAcm("ACM-g2", inputFormat::PDB, "1-butanol", {0}, qcustom);
}

TEST_F (AcmTest, AXgNoSymmPDB)
{
    const std::vector<double> qcustom;
    testAcm("ACM-g2", inputFormat::PDB, "1-butanol", {0}, qcustom);
}

TEST_F (AcmTest, AXgPositive)
{
    const std::vector<double> qcustom;
    testAcm("ACM-g2", inputFormat::SDF, "guanidinium", {1}, qcustom);
}

TEST_F (AcmTest, AXgPositive_FixedQ)
{
    const std::vector<double> qcustom;
    // This call is expected to throw gmx::InvalidInputError due to invalid input.
    EXPECT_THROW(testAcm("ACM-g-fixed", inputFormat::SDF, "guanidinium", {1}, qcustom), gmx::InvalidInputError);
}

TEST_F (AcmTest, AXpgPositive)
{
    const std::vector<double> qcustom;
    testAcm("ACM-pg2", inputFormat::SDF, "guanidinium", {1}, qcustom);
}

TEST_F (AcmTest, SQEgNeutral)
{
    const std::vector<double> qcustom;
    testAcm("ACS-g4", inputFormat::PDB, "1-butanol", {0}, qcustom);
}

TEST_F (AcmTest, SQEgPositive)
{
    const std::vector<double> qcustom;
    testAcm("ACS-g4", inputFormat::SDF, "guanidinium", {1}, qcustom);
}

TEST_F (AcmTest, SQEpgNeutral)
{
    const std::vector<double> qcustom;
    testAcm("ACS-pg2", inputFormat::PDB, "1-butanol", {0}, qcustom);
}

TEST_F (AcmTest, SQEpgPositive)
{
    const std::vector<double> qcustom;
    testAcm("ACS-pg2", inputFormat::SDF, "guanidinium", {1}, qcustom);
}

TEST_F (AcmTest, WaterDimerACSpg)
{
    const std::vector<double> qcustom;
    testAcm("ACS-pg2", inputFormat::PDB, "water_dimer", {0,0}, qcustom);
}

TEST_F (AcmTest, CustomButanolPDB)
{
    const std::vector<double> qcustom = { -0.18, 0.06, 0.06, 0.06, -0.12, 0.06, 0.06,  -0.12, 0.06, 0.06, -0.12, 0.06, 0.06, -0.4, 0.4};
    testAcm("ACM-g2", inputFormat::PDB, "1-butanol", {0}, qcustom);
}

TEST_F (AcmTest, CustomGuanidiniumSDF)
{
    const std::vector<double> qcustom = { -0.3, -0.3, -0.3, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3 };
    testAcm("ACM-g2", inputFormat::SDF, "guanidinium", {1}, qcustom);
}

TEST_F (AcmTest, MethanolWaterACSg)
{
    const std::vector<double> qcustom;
    testAcm("ACS-g4", inputFormat::SDF, "methanol-water", {0,0}, qcustom);
}

TEST_F (AcmTest, MethanolWaterACSpg)
{
    const std::vector<double> qcustom;
    testAcm("ACS-pg2", inputFormat::SDF, "methanol-water", {0,0}, qcustom);
}

TEST_F (AcmTest, AcetateWaterACSg)
{
    const std::vector<double> qcustom;
    testAcm("ACS-g4", inputFormat::SDF, "acetate-water", {-1,0}, qcustom);
}

TEST_F (AcmTest, AcetateWaterACSpg)
{
    const std::vector<double> qcustom;
    testAcm("ACS-pg2", inputFormat::SDF, "acetate-water", {-1,0}, qcustom);
}


//! \return The charges of SQE for butanol, which must be supported if ok is true
static std::vector<double> butanolSQE(ForceField *pd, bool ok)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    std::vector<MolProp> molprops;
    double qtot = 0;
    importFile(&msghandler, pd, gmx::test::TestFileManager::getInputFilePath("molecules/1-butanol.pdb"),
               &molprops, "minimum", JobType::OPT, false, &qtot, true);
    EXPECT_TRUE(msghandler.ok());
    ACTMol mol;
    mol.Merge(&molprops[0]);
    mol.GenerateTopology(&msghandler, pd, missingParameters::Ignore);
    EXPECT_TRUE(msghandler.ok());
    auto    atoms = mol.atomsConst();
    const auto x     = mol.xOriginal();
    QgenAcm qgen(&msghandler, pd, atoms, mol.bondsConst(), 0);
    auto    eqgen = qgen.generateCharges(&msghandler, "1-butanol", pd, &atoms, x, mol.bondsConst());
    EXPECT_EQ(ok ? eQgen::OK : eQgen::NOSUPPORT, eqgen) << qgen.status();
    std::vector<double> q;
    for (const auto &a : atoms)
    {
        q.push_back(a.charge());
    }
    return q;
}

/* Test: SQE for a compound with a bond that has no bond charge correction,
 * which cannot be swapped, is not supported, rather than done with zero
 * corrections for all bonds.
 */
TEST(AcmBondCorrectionTest, MissingBondCorrection)
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACS-g4.xml"), &pd);
    ASSERT_EQ(ChargeGenerationAlgorithm::SQE, pd.chargeGenerationAlgorithm());
    auto &bcc = pd.findForcesConst(InteractionType::BONDCORRECTIONS);
    ASSERT_EQ(CanSwap::No, bcc.canSwap());
    butanolSQE(&pd, true);
    // Remove the corrections of the O-H bond in either order
    auto   *params  = pd.findForces(InteractionType::BONDCORRECTIONS)->parameters();
    size_t  nerased = 0;
    for (auto it = params->begin(); it != params->end(); )
    {
        auto atoms = it->first.atoms();
        if ((atoms[0] == "h_z" && atoms[1] == "o3_z") || (atoms[0] == "o3_z" && atoms[1] == "h_z"))
        {
            it = params->erase(it);
            nerased++;
        }
        else
        {
            ++it;
        }
    }
    ASSERT_GT(nerased, 0u);
    butanolSQE(&pd, false);
}
}

}
