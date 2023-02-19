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
#include <cstdlib>

#include <gtest/gtest.h>

#include "act/alexandria/atype_mapping.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/actmol.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/qgen/qgen_resp.h"
#include "gromacs/gmxlib/network.h"
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

/*! \brief Class to test the RESP algorithm
 */
class RespTest : public gmx::test::CommandLineTestBase
{
protected:
    //! Checking data structure
    gmx::test::TestReferenceChecker checker_;
    
    //! Init set tolerance
    RespTest () : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
    }
        
    ACTMol readMolecule(ForceField *pd)
    {
        std::vector<alexandria::MolProp> molprops;

        // needed for ReadGauss
        const char *molnm    = (char *)"XXX";
        const char *iupac    = (char *)"";
        const char *conf     = (char *)"minimum";
        std::string basis, method;
        const char *jobtype  = (char *)"Opt";
        int         maxpot   = 100;
        int         nsymm    = 0;

        //Read input file for molprop
        auto dataName = gmx::test::TestFileManager::getInputFilePath("1-butanol-3-oep.log");
        double qtot = 0;
        EXPECT_TRUE(readBabel(pd, dataName.c_str(), &molprops, molnm, iupac, conf, &method, &basis,
                              maxpot, nsymm, jobtype, &qtot, false));
                    
        EXPECT_TRUE(qtot == 0.0);
        ACTMol mp;
        mp.Merge(&molprops[0]);
        return mp;
    }

    //! Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    /*! \brief Actual testing routine
     * \param[in] qdist The charge distribution type
     * \param[in] qSymm Whether or not to use charge symmetry
     */
    void testResp(const std::string &qdist, bool qSymm)
    {
        //Generate charges and topology
        std::string   method("B3LYP");
        std::string   basis("Gen");
        t_inputrec    inputrec;
        fill_inputrec(&inputrec);
        ForceField   *pd = getForceField(qdist);
        auto mp = readMolecule(pd);
        auto imm = mp.GenerateTopology(nullptr, pd, missingParameters::Error, false);
        EXPECT_TRUE(immStatus::OK == imm);
        
        //Needed for GenerateCharges
        CommunicationRecord cr;
        auto           pnc         = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
        gmx::MDLogger  mdlog {};
        auto forceComp = new ForceComputer();
        auto qt = pd->findForcesConst(InteractionType::COULOMB);
        auto ct = name2ChargeType(qt.optionValue("chargetype"));
        
        EXPECT_FALSE(ChargeType::Slater  == ct);

        mp.setInputrec(&inputrec);
        mp.symmetrizeCharges(pd, qSymm, nullptr);
        std::vector<gmx::RVec> coords = mp.xOriginal();
        mp.initQgenResp(pd, coords, 0.0, 100);
        std::vector<double> qcustom;
        std::vector<gmx::RVec> forces(mp.atomsConst().size());
        mp.GenerateCharges(pd, forceComp, mdlog, &cr,
                           ChargeGenerationAlgorithm::ESP,
                           qType::ESP,
                           qcustom, &coords, &forces);
        
        std::vector<double> qtotValues;
        auto atoms = mp.atomsConst();
        for (size_t atom = 0; atom < atoms.size(); atom++)
        {
            qtotValues.push_back(atoms[atom].charge());
        }
        char buf[256];
        snprintf(buf, sizeof(buf), "qtotValuesEqdModel_%s",
                 chargeTypeName(ct).c_str());
        checker_.checkSequence(qtotValues.begin(),
                               qtotValues.end(), buf);
    }
        
    //! Cleanup
    static void TearDownTestCase()
    {
    }
};

TEST_F (RespTest, AXpValues)
{
    testResp("ESP-p", false);
}

TEST_F (RespTest, AXgPolarValues)
{
    testResp("ESP-pg", false);
}

TEST_F (RespTest, AXpSymmetricCharges)
{
    testResp("ESP-p", true);
}

TEST_F (RespTest, AXgSymmetricPolarCharges)
{
    testResp("ESP-pg", true);
}

}

} // namespace alexandria
