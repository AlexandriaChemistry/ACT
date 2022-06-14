/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022
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
#include "act/molprop/multipole_names.h"
#include "act/utility/units.h"
#include "../qtype.h"

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
                       const std::string &molname, bool qSymm,
                       double qtotal,
                       std::vector<double> qcustom,
                       bool useCenterOfCharge)
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
            }
            auto dataName = gmx::test::TestFileManager::getInputFilePath(fileName);
 
            double qtot_babel = qtotal;
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
                EXPECT_TRUE(qtotal == qtot_babel);
            }
            alexandria::MyMol mymol;
            mymol.Merge(&molprop);
            // Generate charges and topology
            t_inputrec      inputrecInstance;
            t_inputrec     *inputrec   = &inputrecInstance;
            fill_inputrec(inputrec);
            mymol.setInputrec(inputrec);

            // Get poldata
            auto pd  = getPoldata(model);
            auto imm = mymol.GenerateTopology(stdout, pd,
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
            auto forceComp = new ForceComputer(pd);
            auto alg = ChargeGenerationAlgorithm::NONE;
            if (!qcustom.empty())
            {
                alg = ChargeGenerationAlgorithm::Custom;
            }
            mymol.symmetrizeCharges(pd, qSymm, nullptr);
            mymol.GenerateCharges(pd, forceComp, mdlog, &cr, alg, qcustom);
                                
            std::vector<double> q;
            auto myatoms = mymol.atomsConst();
            for (int atom = 0; atom < myatoms.nr; atom++)
            {
                q.push_back(myatoms.atom[atom].q);
            }
            checker_.checkSequence(q.begin(), q.end(), "Charge");
            QtypeProps qp(qType::Calc);
            
            qp.setQandX(q, mymol.x());
            if (useCenterOfCharge)
            {
                qp.setCenterOfCharge(mymol.centerOfCharge());
            }
            qp.calcMoments();
            for(auto &mpo : mpoMultiPoles)
            {
                auto v = qp.getMultipole(mpo);
                for(auto vv = v.begin(); vv < v.end(); ++vv)
                {
                    *vv = convertFromGromacs(*vv, mpo_unit2(mpo));
                }
                checker_.checkSequence(v.begin(), v.end(), mpo_name(mpo));
            }
        }

        static void TearDownTestCase()
        {
        }

};

TEST_F (QtypeTest, ButanolPDB)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::PDB, "1-butanol", true, 0, qcustom, false);
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
    testQtype("ACS-g", inputFormat::PDB, "1-butanol", true, 0, qcustom, false);
}

TEST_F (QtypeTest, ButanolPDBCoQ)
{
    // Since butanol is dipolar only, the center of charge should have
    // no effect on the dipole but higher moments may be off more.
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::PDB, "1-butanol", true, 0, qcustom, true);
}

TEST_F (QtypeTest, AcetatePDB)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::LOG, "acetate", false, -1, qcustom, false);
}

TEST_F (QtypeTest, AcetatePDBCoQ)
{
    // For the charged compound the center of charge should matter
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::LOG, "acetate", false, -1, qcustom, true);
}

TEST_F (QtypeTest, WaterPDB)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::PDB, "water", true, 0, qcustom, false);
}

TEST_F (QtypeTest, WaterPDBCoQ)
{
    std::vector<double> qcustom;
    testQtype("ACS-g", inputFormat::PDB, "water", true, 0, qcustom, true);
}

TEST_F (QtypeTest, WaterPDBCustom)
{
    std::vector<double> qcustom = { 0.33, -0.66, 0.33 };
    testQtype("ACS-g", inputFormat::PDB, "water", true, 0, qcustom, false);
}

}

}
