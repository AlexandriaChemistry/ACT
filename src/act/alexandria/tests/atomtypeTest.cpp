/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2019-2025
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

#include <map>
#include <tuple>

#include <gtest/gtest.h>

#include "act/alexandria/babel_io.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/actmol.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/qgen/qgen_acm.h"
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

class AtomtypeTest : public gmx::test::CommandLineTestBase
{
private:
    ForceField *pd;

public:   
    AtomtypeTest()
    {
        pd = getForceField("ACS-g");
    }
    
protected:
    void checkAtypes(const std::vector<Bond>  &bonds,
                     const Experiment        &exper)
    {
        // Element name, list of bondorders
        typedef std::tupler<std::string, std::vector<double> > elem_bond;
        std::map<elem_bond, std::string> atype_map = {
            { { "H", 1, 1 }, "h" },
            { { "He", 0, 0 }, "He" },
            { { "Li", 0,0 }, "Li+" },
            { { "Be", 0,0 }, "Be2+" },
            { { "B", 3, 3 }, "b" },
            { { "C", 4, 4 }, "c3" },
            { { "C", 3, 4 }, "c2" },
            { { "C", 2, 4 }, "c1" },
            { { "N", 4, 4 }, "n4" },
            { { "N", 3, 3 }, "n3" },
            { { "N", 2, 3 }, "n2" },
            { { "N", 1, 3 }, "n1" },
            { { "O", 2, 2 }, "o3" },
            { { "O", 1, 2 }, "o2" },
            { { "F", 0, 0 }, "F-" },
            { { "F", 1, 1 }, "f" },
            { { "Ne", 0, 0 }, "Ne" },
            { { "Na", 0, 0 }, "Na+" },
            { { "Mg", 0, 0 }, "Mg2+" },
            { { "Al", 0, 0 }, "Al3+" },
            { { "Si", 4, 4 }, "si3" },
            { { "P", 5, 5 }, "p5" },
            { { "P", 4, 5 }, "p4" },
            { { "P", 3, 3 }, "p3" },
            { { "P", 2, 3 }, "p2" },
            { { "P", 1, 3 }, "p1" },
            { { "S", 6, 6 }, "s6" },
            { { "S", 4, 6 }, "s4" },
            { { "S", 2, 2 }, "s3" },
            { { "S", 1, 2 }, "s2" },
            { { "Cl", 0, 0 }, "Cl-" }, 
            { { "Cl", 1, 1 }, "cl" },
            { { "Ar", 0, 0 }, "Ar" },
            { { "K", 0, 0 }, "K+" },
            { { "Ca", 0, 0 }, "Ca2+" },
            { { "Br", 0, 0 }, "Br-" },
            { { "Br", 1, 1 }, "br" },
            { { "I", 0, 0 }, "I-" },
            { { "I", 1, 1 }, "i" }
        };
        std::vector<elem_bond>   elem_bonds;
        std::vector<std::string> obTypes;
        for (auto &ca : exper.calcAtomConst())
        {
            elem_bonds.push_back( { ca.getName(), 0 } );
            obTypes.push_back(ca.getObtype());
        }
        for (const auto &b: bonds)
        {
            int ai = b.aI();
            int aj = b.aJ();
            elem_bonds[ai].second += 1;
            elem_bonds[aj].second += 1;
            elem_bonds[ai].third += b.bondOrder();
            elem_bonds[aj].third += b.bondOrder();
        }
        std::vector<std::string> atype(elem_bonds.size());
        size_t i = 0;
        for(const auto &eb : elem_bonds)
        {
            auto aptr = atype_map.find(eb);
            if (aptr != atype_map.end())
            {
                printf("Expected atom %zu %s computed %s nbonds %d bondorder %g\n",
                       i, obTypes[i].c_str(), aptr->second.c_str(), eb.second, eb.third);
                //EXPECT_TRUE(aptr->second == obTypes[i]);
            }
            else
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Cannot find atom %s with %g sum of bond oders in atype_map", eb.first.c_str(), eb.second)));
            }
            i += 1;
        }
    }
    void testAtype(const char *molname, bool bondTest=false)
        {
            gmx::test::TestReferenceChecker checker_(this->rootChecker());
            int                             maxpot   = 100;
            int                             nsymm    = 0;
            const char                     *conf     = (char *)"minimum";
            std::string                     basis, method;
            const char                     *jobtype  = (char *)"Opt";

            std::string                     dataName;
            std::vector<alexandria::MolProp> molprops;

            dataName        = gmx::test::TestFileManager::getInputFilePath(molname);
            bool   userqtot = false;
            double qtot     = 0.0;
            matrix box;
            bool readOK = readBabel(nullptr, pd, dataName.c_str(), &molprops,
                                    molname, molname,
                                    conf, &method, &basis, maxpot, nsymm,
                                    jobtype, userqtot, &qtot, false, box, true);
            EXPECT_TRUE(readOK);
            for(auto &molprop: molprops)
            {
                std::vector<std::string> atypes;
                auto                     exper = molprop.experimentConst().begin();
                for (auto &ca : exper->calcAtomConst())
                {
                    atypes.push_back(ca.getObtype());
                }
                // New!
                checkAtypes(molprop.bondsConst(), *exper);
                if (!bondTest)
                {
                    checker_.checkInteger(static_cast<int>(atypes.size()), molname);
                    checker_.checkSequence(atypes.begin(), atypes.end(), "atomtypes");
                }
                else
                {
                    std::vector<std::string> bondorder;
                    for (auto &bond : molprop.bondsConst())
                    {
                        auto ai = bond.aI();
                        auto aj = bond.aJ();
                        bondorder.push_back(gmx::formatString("%s-%d %s-%d: %g",
                                                              atypes[ai].c_str(), ai,
                                                              atypes[aj].c_str(), aj,
                                                              bond.bondOrder()));
                    }
                    checker_.checkInteger(static_cast<int>(bondorder.size()), molname);
                    checker_.checkSequence(bondorder.begin(), bondorder.end(), "bondorder");
                }
            }
        }
};

class BondtypeTest : public AtomtypeTest
{
    protected:
    void testAtype(const char *molname)
    {
        AtomtypeTest::testAtype(molname, true);
    }
};

#include "atom_bond_test_include.cpp"

}

}
