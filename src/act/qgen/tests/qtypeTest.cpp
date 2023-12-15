/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022,2023
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

#include <gtest/gtest.h>

#include "act/alexandria/atype_mapping.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/actmol.h"
#include "act/molprop/multipole_names.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/qgen/qtype.h"
#include "act/utility/units.h"

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
    XYZ,
    ZIP
};

//! Class to test the Alexandria Charge Model
class QtypeTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;

        //init set tolerance
        QtypeTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-5);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        void testQtype(const std::string &model, inputFormat inputformat, 
                       const std::string &molname,
                       double qtotal,
                       std::vector<double> qcustom)
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
            
            switch (inputformat)
            {
            case inputFormat::LOG:
                {
                    fileName.append("-3-oep.log");
                    method.assign("B3LYP");
                    basis.assign("GEN");
                    trustObCharge = true;
                }
                break;
            case inputFormat::ZIP:
                {
                    method.assign("B3LYP");
                    basis.assign("GEN");
                    trustObCharge = true;
                }
                break;
            case inputFormat::PDB:
                {
                    fileName.append(".pdb");
                }
                break;
            case inputFormat::SDF:
                {
                    fileName.append(".sdf");
                }
                break;
            case inputFormat::XYZ:
                {
                    fileName.append(".xyz");
                }
                break;
            }
            // Get forcefield
            auto pd  = getForceField(model);

            auto dataName = gmx::test::TestFileManager::getInputFilePath(fileName);
            double qtot_babel = qtotal;
            matrix box;
            EXPECT_TRUE(readBabel(pd, dataName.c_str(), &molprops,
                                  molname.c_str(), molname.c_str(),
                                  conf, &method, &basis, maxpot,
                                  nsymm, jobtype, &qtot_babel, false, box));

            if (trustObCharge)
            {
                EXPECT_TRUE(qtotal == qtot_babel);
            }
            auto forceComp = new ForceComputer();

            // Now loop over molprops, there may be more than one
            int mp_index = 1;
            for(const auto &molprop : molprops)
            {
                alexandria::ACTMol actmol;

                actmol.Merge(&molprop);
                // Generate charges and topology
                auto imm = actmol.GenerateTopology(stdout, pd,
                                                   missingParameters::Error);
                if (immStatus::OK != imm)
                {
                    fprintf(stderr, "Error generating topology: %s\n", immsg(imm));
                    break;
                }

                // Needed for GenerateCharges
                std::vector<gmx::RVec> forces(actmol.atomsConst().size());
                std::vector<gmx::RVec> coords = actmol.xOriginal();
                auto alg = ChargeGenerationAlgorithm::NONE;
                if (!qcustom.empty())
                {
                    alg = ChargeGenerationAlgorithm::Custom;
                }
                actmol.GenerateCharges(pd, forceComp, alg, qType::Calc, qcustom, &coords, &forces);
                
                std::vector<double> q;
                auto myatoms = actmol.atomsConst();
                for (size_t atom = 0; atom < myatoms.size(); atom++)
                {
                    q.push_back(myatoms[atom].charge());
                }
                std::string qlabel("Charge");
                if (molprops.size() > 1)
                {
                    qlabel += gmx::formatString(" %d", mp_index); 
                }

                checker_.checkSequence(q.begin(), q.end(), qlabel.c_str());
                QtypeProps qp(qType::Calc, myatoms, coords);
                qp.initializeMoments();
                qp.setQandX(q, coords);
                qp.calcMoments();
                for(auto &mpo : mpoMultiPoles)
                {
                    auto v = qp.getMultipole(mpo);
                    for(auto vv = v.begin(); vv < v.end(); ++vv)
                    {
                        *vv = convertFromGromacs(*vv, mpo_unit2(mpo));
                    }
                    std::string mlabel(mpo_name(mpo));
                    if (molprops.size() > 1)
                    {
                        mlabel += gmx::formatString(" %d", mp_index); 
                    }
                    checker_.checkSequence(v.begin(), v.end(), mlabel.c_str());
                }
                mp_index++;
            }
        }

        static void TearDownTestCase()
        {
        }

};

TEST_F (QtypeTest, ButanolPDB)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::PDB, "1-butanol", 0, qcustom);
}

TEST_F (QtypeTest, TwoMolsSDF)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::SDF, "two_mols", 0, qcustom);
}

TEST_F (QtypeTest, TwoMolsXYZ)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::XYZ, "two_mols", 0, qcustom);
}

TEST_F (QtypeTest, CustomButanolPDB)
{
    // These charges are what is produced by the ButanolPDB test.
    // Therefore they should generate the same moments.
    std::vector<double> qcustom = { -0.37705087402262927,
        0.12475677224623878,
        0.12917658810241425,
        0.12917607547179596,
        -0.18000671022216042,
        0.10253972188098517,
        0.10244050985557694,
        -0.21230084461738286,
        0.12543997208211469,
        0.12543946651948318,
        -0.12516784610966791,
        0.10143385295933868,
        0.10143430364294093,
        -0.34806032890043453,
        0.20074934111138645 };
    testQtype("ACS-g", inputFormat::PDB, "1-butanol", 0, qcustom);
}

TEST_F (QtypeTest, ButanolPDBCoQ)
{
    // Since butanol is dipolar only, the center of charge should have
    // no effect on the dipole but higher moments may be off more.
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::PDB, "1-butanol", 0, qcustom);
}

TEST_F (QtypeTest, AcetatePDB)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::LOG, "acetate", -1, qcustom);
}

TEST_F (QtypeTest, AcetatePDBCoQ)
{
    // For the charged compound the center of charge should matter
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::LOG, "acetate", -1, qcustom);
}

TEST_F (QtypeTest, WaterPDB)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::PDB, "water", 0, qcustom);
}

TEST_F (QtypeTest, WaterPDBCoQ)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::PDB, "water", 0, qcustom);
}

TEST_F (QtypeTest, WaterPDBCustom)
{
    std::vector<double> qcustom = { 0.33, -0.66, 0.33 };
    testQtype("ACS-g", inputFormat::PDB, "water", 0, qcustom);
}

}

}
