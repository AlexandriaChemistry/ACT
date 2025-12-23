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
#include <cstdlib>

#include <gtest/gtest.h>

#include "act/import/atype_mapping.h"
#include "act/import/import.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/actmol.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/qgen/qgen_resp.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/mdmodules.h"
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

    //! \return a specific molecule        
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
        auto dataName   = gmx::test::TestFileManager::getInputFilePath("1-butanol-3-oep.log");
        bool   userqtot = false;
        double qtot     = 0;
        matrix box;
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        importFile(&msghandler, pd, dataName.c_str(), &molprops,
                   molnm, iupac, conf, &method, &basis,
                   maxpot, nsymm, jobtype, userqtot, &qtot,
                   false, box, true);
        EXPECT_TRUE(msghandler.ok());
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
     */
    void testResp(const std::string &qdist)
    {
        //Generate charges and topology
        std::string   method("B3LYP");
        std::string   basis("Gen");
        ForceField   *pd = getForceField(qdist);
        auto mp = readMolecule(pd);
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        mp.GenerateTopology(&msghandler, pd, missingParameters::Ignore);
        EXPECT_TRUE(msghandler.ok());
        if (!msghandler.ok())
        {
            return;
        }
        // Needed for GenerateCharges
        auto forceComp = new ForceComputer();
        auto qt = pd->findForcesConst(InteractionType::ELECTROSTATICS);
        auto ct = potentialToChargeDistributionType(qt.potential());
        
        EXPECT_FALSE(ChargeDistributionType::Slater  == ct);

        std::vector<gmx::RVec> coords = mp.xOriginal();
        
        //auto qprops = mp.qProps();
        //qprops->push_back(ACTQprop(mp.topology()->atoms(), coords));
        std::map<MolPropObservable, iqmType> iqm = {
            { MolPropObservable::POTENTIAL, iqmType::QM }
        };
        mp.getExpProps(&msghandler, pd, iqm, 0, 100);

        std::vector<double> qcustom;
        std::vector<gmx::RVec> forces(mp.atomsConst().size());
        mp.generateCharges(&msghandler, pd, forceComp, ChargeGenerationAlgorithm::ESP,
                           &coords, &forces, true);
        
        std::vector<double> qtotValues;
        auto atoms = mp.atomsConst();
        double qtotal = 0;
        for (size_t atom = 0; atom < atoms.size(); atom++)
        {
            qtotValues.push_back(atoms[atom].charge());
            qtotal += atoms[atom].charge();
        }
        EXPECT_TRUE(std::abs(qtotal) < 1e-3);
        char buf[256];
        snprintf(buf, sizeof(buf), "qtotValuesEqdModel_%s",
                 chargeDistributionTypeName(ct).c_str());
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
    testResp("ESP-p");
}

TEST_F (RespTest, AXgPolarValues)
{
    testResp("ESP-pg");
}

TEST_F (RespTest, AXgSymmetricPolarCharges)
{
    testResp("ESP-pg");
}

}

} // namespace alexandria
