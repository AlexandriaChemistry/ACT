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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <math.h>

#include <map>

#include <gtest/gtest.h>

#include "act/poldata/poldata.h"
#include "act/poldata/poldata_utils.h"
#include "act/poldata/poldata_xml.h"
#include "act/qgen/qgen_acm.h"
#include "alexandria/atype_mapping.h"
#include "alexandria/babel_io.h"
#include "alexandria/fill_inputrec.h"
#include "alexandria/molhandler.h"
#include "alexandria/mymol.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class MolHandlerTest : public gmx::test::CommandLineTestBase
{
protected:
    void test(const char *molname, const char *forcefield)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        int                             maxpot   = 100;
        int                             nsymm    = 0;
        const char                     *conf     = (char *)"minimum";
        std::string                     method, basis;
        const char                     *jobtype  = (char *)"Opt";
        
        std::string                     dataName;
        alexandria::MolProp             molprop;
        
        dataName = gmx::test::TestFileManager::getInputFilePath(molname);
        double qtot = 0;
        bool readOK = readBabel(dataName.c_str(), &molprop, molname, molname,
                                conf, &method, &basis,
                                maxpot, nsymm, jobtype, &qtot, false);
        EXPECT_TRUE(readOK);
        if (readOK)
        {
            std::map<std::string, std::string> g2a;
            gaffToAlexandria("", &g2a);
            if (!g2a.empty())
            {
                EXPECT_TRUE(renameAtomTypes(&molprop, g2a));
            }
        }
        MyMol mp_;
        mp_.Merge(&molprop);
        // Generate charges and topology
        t_inputrec      inputrecInstance;
        t_inputrec     *inputrec   = &inputrecInstance;
        fill_inputrec(inputrec);
        mp_.setInputrec(inputrec);
        
        // Get poldata
        auto pd  = getPoldata(forcefield);
        auto imm = mp_.GenerateTopology(stdout, pd, method, basis,
                                        missingParameters::Error, true);
        EXPECT_TRUE(immStatus::OK == imm);
        if (immStatus::OK != imm)
        {
            fprintf(stderr, "Could not generate topology because '%s'. Used basis %s and method %s.\n",
                    immsg(imm), basis.c_str(), method.c_str());
            return;
        }
        // Needed for GenerateCharges
        CommunicationRecord cr;
        auto           pnc      = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
        gmx::MDLogger  mdlog {};
        std::string    lot      = gmx::formatString("%s/%s", method.c_str(), basis.c_str());
        auto alg = ChargeGenerationAlgorithm::NONE;
        std::vector<double> qcustom;
        bool qSymm = false;
        mp_.symmetrizeCharges(pd, qSymm, nullptr);
        mp_.GenerateCharges(pd, mdlog, &cr, alg, qcustom, lot);
        
        MolHandler mh;
        
        double rmsd = 0;
        imm = mh.minimizeCoordinates(&mp_, &rmsd);
        EXPECT_TRUE(immStatus::OK == imm);
        if (immStatus::OK != imm)
        {
            fprintf(stderr, "Could not minimize energy because '%s'\n", immsg(imm));
            return;
        }
        checker_.checkReal(rmsd, "Coordinate RMSD after minimizing");
        for(int i = 0; i < F_NRE; i++)
        {
            if (mp_.energyTerms()[i] != 0)
            {
                checker_.checkReal(mp_.energyTerms()[i], 
                                   interaction_function[i].longname);
            }
        }
        std::vector<double> freq, inten;
        mh.nma(&mp_, &freq, &inten, nullptr);
        checker_.checkSequence(freq.begin(), freq.end(), "Frequencies");
        checker_.checkSequence(inten.begin(), inten.end(), "Intensities");
    }
};

#if CMAKE_BUILD_TYPE == CMAKE_BUILD_TYPE_RELEASE
TEST_F (MolHandlerTest, CarbonDioxide)
{
    test("carbon-dioxide.sdf", "ACS-g");
}

TEST_F (MolHandlerTest, HydrogenChloride)
{

    test("hydrogen-chloride.sdf", "ACS-g");
}

TEST_F (MolHandlerTest, Water)
{

    test("water-3-oep.log.pdb", "ACS-g");
}

TEST_F (MolHandlerTest, Acetone)
{
    test("acetone-3-oep.log.pdb", "ACS-g");
}

TEST_F (MolHandlerTest, Uracil)
{
    test("acetone-3-oep.log.pdb", "ACS-g");
}
#endif

} // namespace

} // namespace alexandria
