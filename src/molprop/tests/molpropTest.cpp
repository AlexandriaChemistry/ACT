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
#include "molprop/molprop.h"
#include "molprop/molpropobservable.h"
#include "molprop/molprop_xml.h"
#include "poldata/poldata_xml.h"

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
        
        std::string mpFile = fileManager().getInputFilePath("molprop.xml");
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
                snprintf(buf, sizeof(buf), "atoms %d %d order %g", 1+bi.aI(), 1+bi.aJ(), bi.bondOrder());
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
            int  calc  = 1;
            int  nener = 1;
            snprintf(mbuf, sizeof(mbuf), "molecule %d, %s, number of calcs",
                     mol++, mpi.getMolname().c_str());
            myCheck.checkInteger(mpi.NExperiment(), mbuf);
            for (auto &ci : mpi.experimentConst())
            {
                char cbuf[256];
                snprintf(cbuf, sizeof(cbuf), "%s calc %d %s",
                         mpi.getMolname().c_str(), calc++,
                         dataSourceName(ci.dataSource()));
                if (ci.dataSource() == dsTheory)
                {
                    snprintf(mbuf, sizeof(mbuf), "%s program", cbuf);
                    myCheck.checkString(ci.getProgram(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s basisset", cbuf);
                    myCheck.checkString(ci.getBasisset(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s method", cbuf);
                    myCheck.checkString(ci.getMethod(), mbuf);
                }
                for (auto &propi : ci.propertyConst())
                {
                    for(auto gp : propi.second)
                    {
                        std::map<const char *,std::pair<int, int> > t_elem = {
                            { "XX", { XX, XX } },
                            { "YY", { YY, YY } },
                            { "ZZ", { ZZ, ZZ } },
                            { "XY", { XX, YY } },
                            { "XZ", { XX, ZZ } },
                            { "YZ", { YY, ZZ } } 
                        }; 
                        const char *mpostr = mpo_name(propi.first);
                        switch(propi.first)
                        {
                        case alexandria::MolPropObservable::POLARIZABILITY:
                            {
                                myCheck.checkDouble(gp->getValue(),
                                                    gmx::formatString("%s %s %s average", mpostr, gp->getType(), cbuf).c_str());
                                myCheck.checkDouble(gp->getError(),
                                                    gmx::formatString("%s %s %s error", mpostr, gp->getType(), cbuf).c_str());
                                auto pp = gp->getTensor();
                                for(auto &t : t_elem)
                                {
                                    myCheck.checkDouble(pp[t.second.first][t.second.second],
                                                        gmx::formatString("%s %s %s %s", mpostr, gp->getType(), cbuf, t.first).c_str());
                                }
                            }
                            break;
                        case MolPropObservable::QUADRUPOLE:
                            {
                                auto vv = gp->getTensor();
                                for(auto &t : t_elem)
                                {
                                    myCheck.checkDouble(vv[t.second.first][t.second.second],
                                                        gmx::formatString("%s %s %s %s", mpostr, gp->getType(), cbuf, t.first).c_str());
                                }
                            }
                            break;
                        case MolPropObservable::DIPOLE:
                            {
                                std::map<const char *,int> v_elem = {
                                    { "X", XX }, { "Y", YY }, { "Z", ZZ }
                                };
                                myCheck.checkDouble(gp->getValue(),
                                                    gmx::formatString("%s %s %s average", mpostr, gp->getType(), cbuf).c_str());
                                myCheck.checkDouble(gp->getError(),
                                                    gmx::formatString("%s %s %s error", mpostr, gp->getType(), cbuf).c_str());
                                auto mu = gp->getVector();
                                for(auto &v : v_elem)
                                {
                                    myCheck.checkDouble(mu[v.second],
                                                        gmx::formatString("%s %s %s %s", mpostr, gp->getType(), cbuf, v.first).c_str());
                                }
                            }
                            break;
                        default:
                            {
                                myCheck.checkDouble(gp->getValue(),
                                                    gmx::formatString("%s %s %s %s %d value",
                                                                      mpostr, gp->getType(), gp->getUnit(),
                                                                      cbuf, nener).c_str());
                                myCheck.checkDouble(gp->getError(),
                                                    gmx::formatString("%s %s %s %s %d error",
                                                                      mpostr, gp->getType(), gp->getUnit(),
                                                                      cbuf, nener).c_str());
                                nener += 1;
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
