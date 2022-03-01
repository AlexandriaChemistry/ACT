/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2021
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


#include <cmath>
#include <cstdlib>

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
#include "../qgen_resp.h"

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
        //! Alexandria molecular properties class
        alexandria::MyMol               mp_;

        //! Init set tolerance
        RespTest () : checker_(this->rootChecker())
        {
            alexandria::MolProp     molprop;

            // needed for ReadGauss
            const char *molnm    = (char *)"XXX";
            const char *iupac    = (char *)"";
            const char *conf     = (char *)"minimum";
            const char *basis    = (char *)"";
            const char *jobtype  = (char *)"Pop";
            int         maxpot   = 100;
            int         nsymm    = 0;

            //Read input file for molprop
            auto dataName = gmx::test::TestFileManager::getInputFilePath("1-butanol-3-oep.log");
            double qtot;
            if (readBabel(dataName.c_str(), &molprop, molnm, iupac, conf, basis,
                          maxpot, nsymm, jobtype, &qtot, false))
            {
                std::map<std::string, std::string> g2a;
                gaffToAlexandria("", &g2a);
                if (!g2a.empty())
                {
                    EXPECT_TRUE(renameAtomTypes(&molprop, g2a));
                }
                else
                {
                    GMX_THROW(gmx::InternalError("Cannot find atomtype mapping file"));
                }
            }
            else
            {
                fprintf(stderr, "Could not read file %s using OpenBabel\n",
                        dataName.c_str());
                return;
            }
            EXPECT_TRUE(qtot == 0.0);
            molprop.SetTotalCharge(qtot);
            mp_.Merge(&molprop);

            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
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
            mp_.SetForceField("alexandria");
            Poldata      *pd = getPoldata(qdist);
            auto imm = mp_.GenerateTopology(nullptr, pd, method, basis,
                                            missingParameters::Error);
            if (immStatus::OK != imm)
            {
                fprintf(stderr, "Error generating topology: %s\n", immsg(imm));
                return;
            }

            //Needed for GenerateCharges
            CommunicationRecord cr;
            auto           pnc         = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
            gmx::MDLogger  mdlog {};
            auto qt = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
            auto ct = name2ChargeType(qt.optionValue("chargetype"));
            std::string    lot(method);
            lot += "/" + basis;
            
            if (ChargeType::Slater  == ct)
            {
                GMX_THROW(gmx::InternalError("No support for tables anymore"));
            }
            mp_.setInputrec(&inputrec);
            mp_.symmetrizeCharges(pd, qSymm, nullptr);
            mp_.initQgenResp(pd, method, basis, 0.0, 100);
            std::vector<double> qcustom;
            mp_.GenerateCharges(pd, mdlog, &cr,
                                ChargeGenerationAlgorithm::ESP, qcustom, lot);

            std::vector<double> qtotValues;
            auto atoms = mp_.atoms();
            for (int atom = 0; atom < atoms->nr; atom++)
            {
                qtotValues.push_back(atoms->atom[atom].q);
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
