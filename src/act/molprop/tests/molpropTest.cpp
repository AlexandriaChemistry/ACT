/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2022
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


#include "actpre.h"

#include <math.h>

#include <cstdlib>

#include <gtest/gtest.h>

#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/snprintf.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molpropobservable.h"
#include "act/molprop/molprop_xml.h"
#include "act/molprop/multipole_names.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/units.h"

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
            myCheck.checkInteger(mpi.experimentConst().size(), mbuf);
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
                for (auto &propi : ci.propertiesConst())
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
                        double fac = convertFromGromacs(1.0, gp->getInputUnit());
                        switch(propi.first)
                        {
                        case alexandria::MolPropObservable::POLARIZABILITY:
                            {
                                auto gpp = static_cast<MolecularPolarizability *>(gp);
                                myCheck.checkDouble(gpp->getValue()*fac,
                                                    gmx::formatString("%s %s %s average", mpostr, gpp->getType(), cbuf).c_str());
                                myCheck.checkDouble(gpp->getError()*fac,
                                                    gmx::formatString("%s %s %s error", mpostr, gpp->getType(), cbuf).c_str());
                                auto pp = gpp->getTensor();
                                for(auto &t : t_elem)
                                {
                                    myCheck.checkDouble(pp[t.second.first][t.second.second]*fac,
                                                        gmx::formatString("%s %s %s %s",
                                                                          mpostr, gpp->getType(), cbuf, t.first).c_str());
                                }
                            }
                            break;
                        case MolPropObservable::DIPOLE:
                        case MolPropObservable::QUADRUPOLE:
                        case MolPropObservable::OCTUPOLE:
                        case MolPropObservable::HEXADECAPOLE:
                            {
                                auto vv = gp->getVector();
                                auto nn = multipoleNames(propi.first);
                                for(size_t i = 0; i < vv.size(); i++)
                                {
                                    myCheck.checkDouble(vv[i]*fac,
                                                        gmx::formatString("%s %s %s %s",
                                                                          mbuf,
                                                                          mpostr, gp->getType(),
                                                                          nn[i].c_str()).c_str());
                                }
                            }
                            break;
                        case MolPropObservable::POTENTIAL:
                            break;
                        default:
                            {
                                myCheck.checkDouble(gp->getValue()*fac,
                                                    gmx::formatString("%s %s %s %s %d value",
                                                                      mpostr, gp->getType(), gp->getUnit(),
                                                                      cbuf, nener).c_str());
                                myCheck.checkDouble(gp->getError()*fac,
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
