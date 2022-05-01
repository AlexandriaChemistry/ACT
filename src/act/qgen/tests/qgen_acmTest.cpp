/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2022
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
#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "alexandria/atype_mapping.h"
#include "alexandria/babel_io.h"
#include "alexandria/fill_inputrec.h"
#include "alexandria/mymol.h"
#include "act/poldata/poldata.h"
#include "act/poldata/poldata_utils.h"
#include "act/poldata/poldata_xml.h"
#include "../qgen_acm.h"

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
        alexandria::MyMol               mp_;

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
                     std::vector<int>                 moleculeStart,
                     const std::vector<const char *>  formula,
                     bool                             useHF)
        {
            int                   maxpot    = 100;
            int                   nsymm     = 0;
            const char           *conf      = (char *)"minimum";
            const char           *jobtype   = (char *)"Opt";
            std::string           method;
            std::string           basis;
            std::string           fileName(molname);
            alexandria::MolProp   molprop;
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

            double qtot_babel = myqtot;
            if (readBabel(dataName.c_str(),
                          &molprop,
                          molname.c_str(),
                          molname.c_str(),
                          conf,
                          &method,
                          &basis,
                          maxpot,
                          nsymm,
                          jobtype,
                          &qtot_babel,
                          false))
            {
                std::map<std::string, std::string> g2a;
                gaffToAlexandria("", &g2a);
                if (!g2a.empty())
                {
                    EXPECT_TRUE(renameAtomTypes(&molprop, g2a));
                }
            }
            else
            {
                fprintf(stderr, "Error reading file %s using OpenBabel.\n",
                        dataName.c_str());
                return;
            }
            if (trustObCharge)
            {
                EXPECT_TRUE(myqtot == qtot_babel);
            }
            if (qtotal.size() != moleculeStart.size())
            {
                GMX_THROW(gmx::InternalError("Different numbers of qtotal and moleculeStart"));
            }
            molprop.clearFragments();
            molprop.generateComposition();
            double qtot_sum = 0;
            for(size_t i = 0; i < qtotal.size(); i++)
            {
                std::vector<int> atomIndices;
                size_t           moleculeEnd = molprop.NAtom();
                if (i < qtotal.size()-1)
                {
                    moleculeEnd = moleculeStart[i+1];
                }
                for(size_t k = moleculeStart[i]; k < moleculeEnd; k++)
                {
                    atomIndices.push_back(k);
                }
                molprop.addFragment(Fragment(std::to_string(i), 0, qtotal[i], 1, formula[i], atomIndices));
                qtot_sum += qtotal[i];
            }

            mp_.Merge(&molprop);
            // Generate charges and topology
            t_inputrec      inputrecInstance;
            t_inputrec     *inputrec   = &inputrecInstance;
            fill_inputrec(inputrec);
            mp_.setInputrec(inputrec);

            // Get poldata
            auto pd  = getPoldata(model);
            auto imm = mp_.GenerateTopology(stdout, pd,
                                            missingParameters::Error);
            if (immStatus::OK != imm)
            {
                fprintf(stderr, "Error generating topology: %s\n", immsg(imm));
                return;
            }

            // Needed for GenerateCharges
            CommunicationRecord cr;
            auto           pnc      = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
            gmx::MDLogger  mdlog {};
            auto alg = ChargeGenerationAlgorithm::NONE;
            if (!qcustom.empty())
            {
                alg = ChargeGenerationAlgorithm::Custom;
            }
            mp_.symmetrizeCharges(pd, qSymm, nullptr);
            mp_.GenerateCharges(pd, mdlog, &cr, alg, qcustom);
                                
            std::vector<double> qtotValues;
            auto myatoms = mp_.atomsConst();
            for (int atom = 0; atom < myatoms.nr; atom++)
            {
                qtotValues.push_back(myatoms.atom[atom].q);
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
                if (atomStart.size() > 2)
                {
                    for(size_t f = 0; f < atomStart.size()-1; f++)
                    {
                        double qt = 0;
                        for(size_t atom = atomStart[f]; atom < atomStart[f+1]; atom++)
                        {
                            qt += myatoms.atom[atom].q;
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
    testAcm("ACS-g", inputFormat::PDB, "water_dimer", true, {0,0}, qcustom, {0,3}, { "H2O", "H2O" },  false);
}

TEST_F (AcmTest, WaterDimerACSgAsMonomer)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::PDB, "water_dimer", true, {0}, qcustom, {0}, { "H4O2"},  false);
}

TEST_F (AcmTest, WaterIodideACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::LOG, "water_I", true, {0,-1}, qcustom, {0,3}, { "H2O", "I-" }, true);
}

TEST_F (AcmTest, WaterIodideACSgAsMonomer)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::LOG, "water_I", true, {-1}, qcustom, {0}, { "H2OI-" }, true);
}

TEST_F (AcmTest, AXpgLOG)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::LOG, "1-butanol", true, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, AXpgZIP)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::ZIP, "trimethylphosphate3-esp.log.gz", true, {0}, qcustom, {0}, { "C3H9PO4" }, false);
}

TEST_F (AcmTest, AXpgPDB)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::PDB, "1-butanol", true, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, AXpgNoSymmLOG)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::LOG, "1-butanol", false, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, AXpgNoSymmPDB)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::PDB, "1-butanol", false, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, AXgLOG)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::LOG, "1-butanol", true, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, AXgPDB)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::PDB, "1-butanol", true, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, AXgNoSymmLOG)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::LOG, "1-butanol", false, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, AXgNoSymmPDB)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::PDB, "1-butanol", false, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, AXgNegative)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::LOG, "acetate", false, {-1}, qcustom, {0}, { "C2H3O2" }, false);
}

TEST_F (AcmTest, AXpgNegative)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::LOG, "acetate", true, {-1}, qcustom, {0}, { "C2H3O2" }, false);
}

TEST_F (AcmTest, AXgPositive)
{
    std::vector<double> qcustom;
    testAcm("ACM-g", inputFormat::SDF, "guanidinium", false, {1}, qcustom, {0}, { "CH6N3" }, false);
}

TEST_F (AcmTest, AXpgPositive)
{
    std::vector<double> qcustom;
    testAcm("ACM-pg", inputFormat::SDF, "guanidinium", true, {1}, qcustom, {0}, { "CH6N3" }, false);
}

TEST_F (AcmTest, SQEgNeutral)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::PDB, "1-butanol", false, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, SQEgNegative)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::LOG, "acetate", false, {-1}, qcustom, {0}, { "C2H3O2" }, false);
}

TEST_F (AcmTest, SQEgPositive)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::SDF, "guanidinium", false, {1}, qcustom, {0}, { "CH6N3" }, false);
}
#define LATER
#ifdef LATER
TEST_F (AcmTest, SQEpgNeutral)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::PDB, "1-butanol", false, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, SQEpgNegative)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::LOG, "acetate", true, {-1}, qcustom, {0}, { "C2H3O2" }, false);
}

TEST_F (AcmTest, SQEpgPositive)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::SDF, "guanidinium", true, {1}, qcustom, {0}, { "CH6N3" }, false);
}

TEST_F (AcmTest, WaterDimerACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::PDB, "water_dimer", true, {0,0}, qcustom, {0,3}, { "H2O", "H2O" }, false);
}

TEST_F (AcmTest, WaterIodideACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::LOG, "water_I", true, {0, -1}, qcustom, {0,3}, { "H2O", "I-" }, true);
}

#endif

TEST_F (AcmTest, CustomButanolPDB)
{
    std::vector<double> qcustom = { -0.18, 0.06, 0.06, 0.06, -0.12, 0.06, 0.06,  -0.12, 0.06, 0.06, -0.12, 0.06, 0.06, -0.4, 0.4};
    testAcm("ACM-g", inputFormat::PDB, "1-butanol", false, {0}, qcustom, {0}, { "C4H10O" }, false);
}

TEST_F (AcmTest, CustomAcetatePDB)
{
    std::vector<double> qcustom = { -0.18, 0.06, 0.06, 0.06, 0.1, -0.55, -0.55 };
    testAcm("ACM-g", inputFormat::LOG, "acetate", false, {-1}, qcustom, {0}, { "C2H3O2" }, false);
}

TEST_F (AcmTest, CustomGuanidiniumSDF)
{
    std::vector<double> qcustom = { -0.3, -0.3, -0.3, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3 };
    testAcm("ACM-g", inputFormat::SDF, "guanidinium", false, {1}, qcustom, {0}, { "CH6N3" }, false);
}

TEST_F (AcmTest, MethanolWaterACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::SDF, "methanol-water", true, {0,0}, qcustom, {0,6}, { "CH4O", "H2O" }, true);
}

TEST_F (AcmTest, MethanolWaterACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::SDF, "methanol-water", true, {0,0}, qcustom, {0,6}, { "CH4O", "H2O" }, true);
}

TEST_F (AcmTest, AcetateWaterACSg)
{
    std::vector<double> qcustom;
    testAcm("ACS-g", inputFormat::SDF, "acetate-water", true, {-1,0}, qcustom, {0,7}, { "C2H3O2", "H2O" }, true);
}

TEST_F (AcmTest, AcetateWaterACSpg)
{
    std::vector<double> qcustom;
    testAcm("ACS-pg", inputFormat::SDF, "acetate-water", true, {-1,0}, qcustom, {0,7}, { "C2H3O2", "H2O" }, true);
}

}

}
