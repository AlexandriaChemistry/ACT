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


#include "actpre.h"

#include <math.h>

#include <cstdlib>

#include <gtest/gtest.h>

#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/snprintf.h"
#include "alexandria/molprop.h"
#include "alexandria/molpropobservable.h"
#include "alexandria/molprop_xml.h"
#include "alexandria/poldata_xml.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

class MolpropTest : public gmx::test::CommandLineTestBase
{
protected:
    //! Vector of molecule properties
    std::vector<alexandria::MolProp>  mp_;
    //! Structure containing atom properties
    gmx_atomprop_t                    aps_;
    //! Constructor that does initialization and reads an input file
    MolpropTest()
    {
        aps_  = gmx_atomprop_init();
        
        std::string mpFile = fileManager().getInputFilePath("molprop.dat");
        MolPropRead(mpFile.c_str(), &mp_);
    }

    //! Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    //! Do the actual testing
    void testMolProp ()
    {
        int mol = 1;
        gmx::test::TestReferenceChecker myCheck(this->rootChecker());
        for (auto &mpi : mp_)
        {
            char mbuf[256];
            snprintf(mbuf, sizeof(mbuf), "molecule %d name", mol);
            myCheck.checkString(mpi.getMolname(), mbuf);
            mpi.GenerateFormula(aps_);
            snprintf(mbuf, sizeof(mbuf), "molecule %d formula", mol);
            myCheck.checkString(mpi.formula(), mbuf);
            snprintf(mbuf, sizeof(mbuf), "molecule %d number of bonds", mol);
            myCheck.checkInteger(mpi.NBond(), mbuf);
            int i = 1;
            for (auto &bi : mpi.bondsConst())
            {
                char buf[256];
                snprintf(buf, sizeof(buf), "atoms %d %d order %g", bi.getAi(), bi.getAj(), bi.getBondOrder());
                std::string bond("bond");
                char        ibuf[256];
                snprintf(ibuf, sizeof(ibuf), "molecule %d bond %d", mol, i++);
                myCheck.checkString(buf, ibuf);
            }
                mol++;
        }
    }

    //! Test the content of alexandria::Experiment structures containing calculations
    void testCalculations ()
    {
        int mol = 1;
        gmx::test::TestReferenceChecker myCheck(this->rootChecker());
        for (auto &mpi : mp_)
        {
            char mbuf[512];
            int  calc = 1;
            snprintf(mbuf, sizeof(mbuf), "molecule %d number of calcs", mol);
            myCheck.checkInteger(mpi.NExperiment(), mbuf);
            for (auto &ci : mpi.experimentConst())
            {
                char cbuf[256];
                snprintf(cbuf, sizeof(cbuf), "molecule %d cakc %d", mol, calc);
                snprintf(mbuf, sizeof(mbuf), "%s program", cbuf);
                myCheck.checkString(ci.getProgram(), mbuf);
                snprintf(mbuf, sizeof(mbuf), "%s basisset", cbuf);
                myCheck.checkString(ci.getBasisset(), mbuf);
                snprintf(mbuf, sizeof(mbuf), "%s method", cbuf);
                myCheck.checkString(ci.getMethod(), mbuf);
                for (auto &propi : ci.propertyConst())
                {
                    for(auto gp : propi.second)
                    {
                        switch(propi.first)
                        {
                        case alexandria::MolPropObservable::POLARIZABILITY:
                            {
                                auto pp = gp->getTensor();
                                snprintf(mbuf, sizeof(mbuf), "%s polar XX", cbuf);
                                myCheck.checkDouble(pp[XX][XX], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s polar YY", cbuf);
                                myCheck.checkDouble(pp[YY][YY], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s polar ZZ", cbuf);
                                myCheck.checkDouble(pp[ZZ][ZZ], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s polar XY", cbuf);
                                myCheck.checkDouble(pp[XX][YY], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s polar XZ", cbuf);
                                myCheck.checkDouble(pp[XX][ZZ], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s polar YZ", cbuf);
                                myCheck.checkDouble(pp[YY][ZZ], mbuf);
                            }
                            break;
                        case MolPropObservable::QUADRUPOLE:
                            {
                                auto vv = gp->getTensor();
                                snprintf(mbuf, sizeof(mbuf), "%s quadrupole XX", cbuf);
                                myCheck.checkDouble(vv[XX][XX], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s quadrupole YY", cbuf);
                                myCheck.checkDouble(vv[YY][YY], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s quadrupole ZZ", cbuf);
                                myCheck.checkDouble(vv[ZZ][ZZ], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s quadrupole XY", cbuf);
                                myCheck.checkDouble(vv[XX][YY], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s quadrupole XZ", cbuf);
                                myCheck.checkDouble(vv[XX][ZZ], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s quadrupole YZ", cbuf);
                                myCheck.checkDouble(vv[YY][ZZ], mbuf);
                            }
                            break;
                        case MolPropObservable::DIPOLE:
                            {
                                auto mu = gp->getVector();
                                snprintf(mbuf, sizeof(mbuf), "%s dipole X", cbuf);
                                myCheck.checkDouble(mu[XX], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s dipole Y", cbuf);
                                myCheck.checkDouble(mu[YY], mbuf);
                                snprintf(mbuf, sizeof(mbuf), "%s dipole Z", cbuf);
                                myCheck.checkDouble(mu[ZZ], mbuf);
                            }
                            break;
                        default:
                            {
                                int nener = 1;
                                char ebuf[256];
                                snprintf(mbuf, sizeof(mbuf), "%s energy %d", cbuf, nener++);
                                snprintf(ebuf, sizeof(ebuf), "%s %g +/- %g %s",
                                         gp->getType(),
                                         gp->getValue(),
                                         gp->getError(),
                                         gp->getUnit());
                                myCheck.checkString(ebuf, mbuf);
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    //! Clean the test data.
    static void TearDownTestCase()
    {
    }
    
};

TEST_F (MolpropTest, NameFormulaBonds){
    testMolProp();
}

TEST_F (MolpropTest, Calculations){
    testCalculations();
}

}
