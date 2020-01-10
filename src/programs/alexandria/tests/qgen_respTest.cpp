/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018
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
#include "programs/alexandria/qgen_resp.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "poldata_utils.h"

class RespTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
        alexandria::MyMol               mp_;
        gmx_atomprop_t                  aps_;

        //init set tolecrance
        RespTest () : checker_(this->rootChecker())
        {
            alexandria::MolProp     molprop;
            aps_ = gmx_atomprop_init();

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
            readBabel(dataName.c_str(), &molprop, molnm, iupac, conf, basis,
                      maxpot, nsymm, jobtype, 0.0, false);
            std::vector<MolProp> vmp;
            vmp.push_back(molprop);
            mp_.molProp()->Merge(vmp.begin());

            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        void testResp(ChargeModel qdist)
        {
            //Generate charges and topology
            std::string   method("B3LYP");
            std::string   basis("Gen");
            t_inputrec    inputrec;
            fill_inputrec(&inputrec);
            mp_.SetForceField("gaff");
            auto imm = mp_.GenerateTopology(aps_, getPoldata(qdist), method,
                                            basis, nullptr,
                                            false, false, false,  false,
                                            nullptr);
            if (immOK != imm)
            {
                fprintf(stderr, "Error generating topology: %s\n", immsg(imm));
                return;
            }
            fprintf(stderr, "Generated topology for %s\n",
                    mp_.molProp()->getMolname().c_str());

            //Needed for GenerateCharges
            real           hfac                  = 0;
            real           watoms                = 0;
            char          *symm_string           = (char *)"";
            t_commrec     *cr                    = init_commrec();
            auto           pnc                   = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
            gmx::MDLogger  mdlog {};
            auto           hwinfo      = gmx_detect_hardware(mdlog, pnc);
            int            qcycle      = 1;
            real           qtol        = 1e-3;
            int            maxpot      = 100;
            std::string    tabFile;
            if (getEemtypeSlater(qdist))
            {
                inputrec.coulombtype = eelUSER;
                tabFile              = fileManager().getInputFilePath("table.xvg");
            }
            mp_.setInputrec(&inputrec);
            mp_.GenerateCharges(getPoldata(qdist), mdlog, aps_, watoms,
                                hfac, method, basis, nullptr,
                                false, symm_string, cr,
                                tabFile.empty() ? nullptr : tabFile.c_str(),
                                hwinfo, qcycle, maxpot, qtol, nullptr, nullptr);

            std::vector<double> qtotValues;
            for (int atom = 0; atom < mp_.mtop_->moltype[0].atoms.nr; atom++)
            {
                qtotValues.push_back(mp_.mtop_->moltype[0].atoms.atom[atom].q);
            }
            char buf[256];
            snprintf(buf, sizeof(buf), "qtotValuesEqdModel_%s",
                     getEemtypeName(qdist));
            checker_.checkSequence(qtotValues.begin(),
                                   qtotValues.end(), buf);
        }

        static void TearDownTestCase()
        {
        }

};

TEST_F (RespTest, AXpValues)
{
    testResp(eqdESP_p);
}

TEST_F (RespTest, AXgPolarValues)
{
    testResp(eqdESP_pg);
}

TEST_F (RespTest, AXsPolarValues)
{
    testResp(eqdESP_ps);
}
