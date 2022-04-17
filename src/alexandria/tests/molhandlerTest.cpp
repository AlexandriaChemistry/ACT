/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2021
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
        int                             maxpot   = 100;
        int                             nsymm    = 0;
        const char                     *conf     = (char *)"minimum";
        const char                     *basis    = (char *)"";
        const char                     *method   = (char *)"";
        const char                     *jobtype  = (char *)"Opt";
        
        std::string                     dataName;
        alexandria::MolProp             molprop;
        
        dataName = gmx::test::TestFileManager::getInputFilePath(molname);
        double qtot = 0;
        bool readOK = readBabel(dataName.c_str(), &molprop, molname, molname,
                                conf, basis, maxpot, nsymm, jobtype, &qtot, false);
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
        std::string    lot      = gmx::formatString("%s/%s", method, basis);
        auto alg = ChargeGenerationAlgorithm::NONE;
        std::vector<double> qcustom;
        bool qSymm = false;
        mp_.symmetrizeCharges(pd, qSymm, nullptr);
        mp_.GenerateCharges(pd, mdlog, &cr, alg, qcustom, lot);
        
        MolHandler mh;
        
        double rmsd = 0;
        imm = mh.minimizeCoordinates(&mp_, &rmsd);
        EXPECT_TRUE(immStatus::OK == imm);
        checker_.checkReal(rmsd, "Coordinate RMSD after minimizing");
        std::vector<double> freq, inten;
        mh.nma(&mp_, &freq, &inten, nullptr);
        checker_.checkSequence(freq.begin(), freq.end(), "Frequencies");
        checker_.checkSequence(inten.begin(), inten.end(), "Intensities");
    }
};

TEST_F (MolHandlerTest, Acetone)
{
    test("acetone-3-oep.log.pdb", "ACS-g");
}

TEST_F (MolHandlerTest, Uracil)
{
    test("acetone-3-oep.log.pdb", "ACS-g");
}

TEST_F (MolHandlerTest, Ammonia)
{
    test("ammonia.sdf", "ACS-g");
}

} // namespace

} // namespace alexandria
