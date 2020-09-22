/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2020
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
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "programs/alexandria/babel_io.h"
#include "programs/alexandria/fill_inputrec.h"
#include "programs/alexandria/mymol.h"
#include "programs/alexandria/poldata.h"
#include "programs/alexandria/poldata_xml.h"
#include "programs/alexandria/qgen_acm.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "poldata_utils.h"

namespace alexandria
{

namespace
{

enum informat{
    einfLOG = 0,
    einfPDB = 1,
    einfNR
};

class AcmTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
        alexandria::MyMol               mp_;
        gmx_atomprop_t                  aps_;

        //init set tolerance
        AcmTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-5);
            checker_.setDefaultTolerance(tolerance);
            aps_ = gmx_atomprop_init();
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        void testAcm(const std::string &model, informat inputformat, bool qSymm)
        {
            int                   maxpot    = 100;
            int                   nsymm     = 0;
            const char           *dihopt[]  = { nullptr, "No", "Single", "All", nullptr };
            const char           *molnm     = (char *)"1-butanol";
            const char           *iupac     = (char *)"1-butanol";
            const char           *conf      = (char *)"minimum";
            const char           *jobtype   = (char *)"Opt";

            std::string           dataName, method, basis;
            alexandria::MolProp   molprop;
            std::vector<MolProp>  vmp;
            
            if (inputformat == einfLOG)
            {
                method.assign("B3LYP");
                basis.assign("GEN");
                dataName = gmx::test::TestFileManager::getInputFilePath("1-butanol-3-oep.log");
            }
            else
            {
                dataName = gmx::test::TestFileManager::getInputFilePath("1-butanol.pdb");
            }

            readBabel(dataName.c_str(),
                      &molprop,
                      molnm,
                      iupac,
                      conf,
                      basis.c_str(),
                      maxpot,
                      nsymm,
                      jobtype,
                      0.0,
                      false);

            vmp.push_back(molprop);
            mp_.Merge(vmp.begin());
            fprintf(stderr, "Read babel for %s\n", dataName.c_str());
            // Generate charges and topology
            eDih            edih       = (eDih) get_option(dihopt);
            t_inputrec      inputrecInstance;
            t_inputrec     *inputrec   = &inputrecInstance;
            fill_inputrec(inputrec);
            mp_.setInputrec(inputrec);

            // Get poldata
            auto pd  = getPoldata(model);
            auto imm = mp_.GenerateTopology(aps_, pd, method, basis, nullptr,
                                            false, false, edih, false, nullptr);
            if (immOK != imm)
            {
                fprintf(stderr, "Error generating topology: %s\n", immsg(imm));
                return;
            }
            fprintf(stderr, "Generated topology for %s\n", dataName.c_str());

            // Needed for GenerateCharges
            t_commrec     *cr       = init_commrec();
            auto           pnc      = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
            gmx::MDLogger  mdlog {};
            auto           hwinfo   = gmx_detect_hardware(mdlog, pnc);
            int            qcycle   = 100;
            real           qtol     = 1e-3;

            mp_.symmetrizeCharges(pd, aps_, qSymm, nullptr);
            mp_.initQgenResp(pd, method, basis, nullptr, 0.0, 100);
            mp_.GenerateCharges(pd, mdlog, cr, nullptr, 
                                hwinfo, qcycle, qtol);
                                
            fprintf(stderr, "Generated charges for %s\n", dataName.c_str());

            std::vector<double> qtotValues;
            for (int atom = 0; atom < mp_.atoms_->nr; atom++)
            {
                qtotValues.push_back(mp_.atoms_->atom[atom].q);
            }
            char buf[256];
            snprintf(buf, sizeof(buf), "qtotValuesEqdAlgorithm_%s", 
                     chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str());
            checker_.checkInteger(static_cast<int>(qtotValues.size()), "qtotSize");
            checker_.checkSequence(qtotValues.begin(), qtotValues.end(), buf);
        }

        static void TearDownTestCase()
        {
        }

};

TEST_F (AcmTest, BultinckLog)
{
    testAcm("Bultinck", einfLOG, true);
}

TEST_F (AcmTest, BultinckPDB)
{
    testAcm("Bultinck", einfPDB, true);
}

TEST_F (AcmTest, VerstraelenLog)
{
    testAcm("Verstraelen", einfLOG, true);
}

TEST_F (AcmTest, VerstraelenPDB)
{
    testAcm("Verstraelen", einfPDB, true);
}

TEST_F (AcmTest, RappeLog)
{
    testAcm("Rappe", einfLOG, true);
}

TEST_F (AcmTest, RappePDB)
{
    testAcm("Rappe", einfPDB, true);
}

TEST_F (AcmTest, YangLog)
{
    testAcm("Yang", einfLOG, true);
}

TEST_F (AcmTest, YangPDB)
{
    testAcm("Yang", einfPDB, true);
}

TEST_F (AcmTest, AXpgLOG)
{
    testAcm("ACM-pg", einfLOG, true);
}

TEST_F (AcmTest, AXpgPDB)
{
    testAcm("ACM-pg", einfPDB, true);
}

TEST_F (AcmTest, AXpgNoSymmLOG)
{
    testAcm("ACM-pg", einfLOG, false);
}

TEST_F (AcmTest, AXpgNoSymmPDB)
{
    testAcm("ACM-pg", einfPDB, false);
}

TEST_F (AcmTest, AXgLOG)
{
    testAcm("ACM-g", einfLOG, true);
}

TEST_F (AcmTest, AXgPDB)
{
    testAcm("ACM-g", einfPDB, true);
}

TEST_F (AcmTest, AXgNoSymmLOG)
{
    testAcm("ACM-g", einfLOG, false);
}

TEST_F (AcmTest, AXgNoSymmPDB)
{
    testAcm("ACM-g", einfPDB, false);
}

}

}
