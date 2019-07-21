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


#include "gmxpre.h"

#include <math.h>

#include <cstdlib>

#include <gtest/gtest.h>

#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/snprintf.h"
#include "programs/alexandria/molprop.h"
#include "programs/alexandria/molprop_xml.h"
#include "programs/alexandria/poldata_xml.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

class MolpropTest : public gmx::test::CommandLineTestBase
{
    protected:
        std::vector<alexandria::MolProp>  mp_;
        gmx_atomprop_t                    aps_;

        MolpropTest()
        {
            aps_  = gmx_atomprop_init();

            std::string mpFile = fileManager().getInputFilePath("molprop.dat");
            MolPropRead(mpFile.c_str(), mp_);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        void testMolProp ()
        {
            int mol = 1;
            gmx::test::TestReferenceChecker myCheck(this->rootChecker());
            for (alexandria::MolPropIterator mpi = mp_.begin(); (mpi < mp_.end()); ++mpi, ++mol)
            {
                char mbuf[256];
                snprintf(mbuf, sizeof(mbuf), "molecule %d name", mol);
                myCheck.checkString(mpi->getMolname(), mbuf);
                mpi->GenerateFormula(aps_);
                snprintf(mbuf, sizeof(mbuf), "molecule %d formula", mol);
                myCheck.checkString(mpi->formula(), mbuf);
                snprintf(mbuf, sizeof(mbuf), "molecule %d number of bonds", mol);
                myCheck.checkInteger(mpi->NBond(), mbuf);
                int i = 1;
                for (alexandria::BondIterator bi = mpi->BeginBond(); (bi < mpi->EndBond()); ++bi)
                {
                    char buf[256];
                    snprintf(buf, sizeof(buf), "atoms %d %d order %d",
                             bi->getAi(), bi->getAj(), bi->getBondOrder());
                    std::string bond("bond");
                    char        ibuf[256];
                    snprintf(ibuf, sizeof(ibuf), "molecule %d bond %d", mol, i++);
                    myCheck.checkString(buf, ibuf);
                }
            }
        }

        void testExperiments ()
        {
            int mol = 1;
            gmx::test::TestReferenceChecker myCheck(this->rootChecker());
            for (alexandria::MolPropIterator mpi = mp_.begin(); (mpi < mp_.end()); ++mpi, ++mol)
            {
                char mbuf[512];
                int  exp = 1;
                snprintf(mbuf, sizeof(mbuf), "molecule %d number of experiments", mol);
                myCheck.checkInteger(mpi->NExperiment(), mbuf);
                for (alexandria::ExperimentIterator expi = mpi->BeginExperiment(); (expi < mpi->EndExperiment()); ++expi, ++exp)
                {
                    char cbuf[256];
                    snprintf(cbuf, sizeof(cbuf), "molecule %d exper %d", mol, exp);
                    int  nener = 1;
                    for (alexandria::MolecularEnergyIterator ei = expi->BeginEnergy(); (ei < expi->EndEnergy()); ++ei)
                    {
                        char ebuf[256];
                        snprintf(mbuf, sizeof(mbuf), "%s energy %d", cbuf, nener++);
                        snprintf(ebuf, sizeof(ebuf), "%s %g +/- %g %s",
                                 ei->getType().c_str(),
                                 ei->getValue(), ei->getError(),
                                 ei->getUnit().c_str());
                        myCheck.checkString(ebuf, mbuf);
                    }
                }
            }
        }

        void testCalculations ()
        {
            int mol = 1;
            gmx::test::TestReferenceChecker myCheck(this->rootChecker());
            for (alexandria::MolPropIterator mpi = mp_.begin(); (mpi < mp_.end()); ++mpi, ++mol)
            {
                char mbuf[512];
                int  calc = 1;
                snprintf(mbuf, sizeof(mbuf), "molecule %d number of calcs", mol);
                myCheck.checkInteger(mpi->NExperiment(), mbuf);
                for (alexandria::ExperimentIterator ci = mpi->BeginExperiment(); (ci < mpi->EndExperiment()); ++ci, ++calc)
                {
                    char cbuf[256];
                    snprintf(cbuf, sizeof(cbuf), "molecule %d cakc %d", mol, calc);
                    snprintf(mbuf, sizeof(mbuf), "%s program", cbuf);
                    myCheck.checkString(ci->getProgram(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s basisset", cbuf);
                    myCheck.checkString(ci->getBasisset(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s method", cbuf);
                    myCheck.checkString(ci->getMethod(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s number of polar", cbuf);
                    myCheck.checkInteger(ci->NPolar(), mbuf);
                    for (alexandria::MolecularPolarizabilityIterator poli = ci->BeginPolar(); (poli < ci->EndPolar()); ++poli)
                    {
                        snprintf(mbuf, sizeof(mbuf), "%s polar XX", cbuf);
                        myCheck.checkDouble(poli->getXX(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s polar YY", cbuf);
                        myCheck.checkDouble(poli->getYY(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s polar ZZ", cbuf);
                        myCheck.checkDouble(poli->getZZ(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s polar XY", cbuf);
                        myCheck.checkDouble(poli->getXY(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s polar XZ", cbuf);
                        myCheck.checkDouble(poli->getXZ(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s polar YZ", cbuf);
                        myCheck.checkDouble(poli->getYZ(), mbuf);
                    }

                    for (alexandria::MolecularQuadrupoleIterator qi = ci->BeginQuadrupole(); (qi < ci->EndQuadrupole()); ++qi)
                    {
                        snprintf(mbuf, sizeof(mbuf), "%s quadrupole XX", cbuf);
                        myCheck.checkDouble(qi->getXX(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s quadrupole YY", cbuf);
                        myCheck.checkDouble(qi->getYY(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s quadrupole ZZ", cbuf);
                        myCheck.checkDouble(qi->getZZ(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s quadrupole XY", cbuf);
                        myCheck.checkDouble(qi->getXY(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s quadrupole XZ", cbuf);
                        myCheck.checkDouble(qi->getXZ(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s quadrupole YZ", cbuf);
                        myCheck.checkDouble(qi->getYZ(), mbuf);
                    }

                    for (alexandria::MolecularDipoleIterator dip = ci->BeginDipole(); (dip < ci->EndDipole()); ++dip)
                    {
                        snprintf(mbuf, sizeof(mbuf), "%s dipole X", cbuf);
                        myCheck.checkDouble(dip->getX(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s dipole Y", cbuf);
                        myCheck.checkDouble(dip->getY(), mbuf);
                        snprintf(mbuf, sizeof(mbuf), "%s dipole Z", cbuf);
                        myCheck.checkDouble(dip->getZ(), mbuf);
                    }

                    int nener = 1;
                    for (alexandria::MolecularEnergyIterator ei = ci->BeginEnergy(); (ei < ci->EndEnergy()); ++ei)
                    {
                        char ebuf[256];
                        snprintf(mbuf, sizeof(mbuf), "%s energy %d", cbuf, nener++);
                        snprintf(ebuf, sizeof(ebuf), "%s %g +/- %g %s",
                                 ei->getType().c_str(),
                                 ei->getValue(), ei->getError(),
                                 ei->getUnit().c_str());
                        myCheck.checkString(ebuf, mbuf);
                    }
                }
            }
        }

        static void TearDownTestCase()
        {
        }

};

TEST_F (MolpropTest, NameFormulaBonds){
    testMolProp();
}

TEST_F (MolpropTest, Experiments){
    testExperiments();
}
