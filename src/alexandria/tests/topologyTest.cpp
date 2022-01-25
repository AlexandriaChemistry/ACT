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

#include <cstdlib>

#include <gtest/gtest.h>

#include "alexandria/topology.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

class TopologyTest : public gmx::test::CommandLineTestBase
{
protected:
    std::vector<Bond> bonds_;
    gmx::HostVector<gmx::RVec> x_, y_;

    //! Constructor that does initialization and reads an input file
    TopologyTest()
    {
        bonds_ = {
            { Bond(1, 2, 1.0) },
            { Bond(2, 3, 2.0) },
            { Bond(2, 4, 1.0) },
            { Bond(4, 5, 1.0) },
            { Bond(5, 6, 1.0) },
            { Bond(6, 1, 1.0) },
        };
        x_ = {
            { 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0 },
            { 1.0, 0.0, 0.0 },
            { 1.6, 0.3, 0.0 },
            { 1.3,-0.6, 0.0 },
            { 1.3, 0.6, 0.0 },
            { 0.6, 0.3, 0.0 },
        };
        y_ = {
            { 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0 },
            { 1.0, 0.0, 0.0 },
            { 2.0, 0.0, 0.0 },
            { 1.6,-0.6, 0.0 },
            { 1.2, 0.6, 0.0 },
            { 0.6, 0.3, 0.0 },
        };
    }

    //! Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    //! Do the actual testing
    void testTopology ()
    {
    }

    //! Test the generation of angles
    void testAngles (const gmx::HostVector<gmx::RVec> &xx)
    {
        gmx::test::TestReferenceChecker myCheck(this->rootChecker());
        Topology top(bonds_);
        top.makeAngles(xx, 175.0);
        auto &angles = top.entry(InteractionType::ANGLES);
        myCheck.checkInteger(angles.size(), "#angles");
        std::vector<std::string> astrings;
        for(auto &a: angles)
        {
            auto ai = a->atomIndices();
            astrings.push_back(gmx::formatString("%d-%d-%d", ai[0], ai[1], ai[2]));
        }
        myCheck.checkSequence(astrings.begin(), astrings.end(), "Angle");
    }

    //! Test the generation of linear angles
    void testLinearAngles (const gmx::HostVector<gmx::RVec> &xx)
    {
        gmx::test::TestReferenceChecker myCheck(this->rootChecker());
        Topology top(bonds_);
        top.makeAngles(xx, 175.0);
        {
            auto &angles = top.entry(InteractionType::ANGLES);
            myCheck.checkInteger(angles.size(), "#angles");
            std::vector<std::string> astrings;
            for(auto &a: angles)
            {
                auto ai = a->atomIndices();
                astrings.push_back(gmx::formatString("%d-%d-%d", ai[0], ai[1], ai[2]));
            }
            myCheck.checkSequence(astrings.begin(), astrings.end(), "Angle");
        }
        {
            auto &linangles = top.entry(InteractionType::LINEAR_ANGLES);
            std::vector<std::string> lastrings;
            for(auto &a: linangles)
            {
                auto ai = a->atomIndices();
                lastrings.push_back(gmx::formatString("%d-%d-%d", ai[0], ai[1], ai[2]));
            }
            myCheck.checkSequence(lastrings.begin(), lastrings.end(), "LinearAngle");
        }
    }

    //! Test the generation of impropers
    void testImpropers (double PlanarAngleMax)
    {
        gmx::test::TestReferenceChecker myCheck(this->rootChecker());
        Topology top(bonds_);
        top.makeImpropers(x_, PlanarAngleMax);
        auto &imps = top.entry(InteractionType::IMPROPER_DIHEDRALS);
        myCheck.checkInteger(imps.size(), "#impropers");
        std::vector<std::string> astrings;
        for(auto &a: imps)
        {
            auto ai = a->atomIndices();
            astrings.push_back(gmx::formatString("%d-%d-%d-%d", ai[0], ai[1], ai[2], ai[3]));
        }
        myCheck.checkSequence(astrings.begin(), astrings.end(), "Improper");
    }

    //! Test the generation of propers
    void testPropers (const gmx::HostVector<gmx::RVec> &xx)
    {
        gmx::test::TestReferenceChecker myCheck(this->rootChecker());
        Topology top(bonds_);
        top.makeAngles(xx, 175.0);
        top.makePropers();
        auto &imps = top.entry(InteractionType::PROPER_DIHEDRALS);
        myCheck.checkInteger(imps.size(), "#propers");
        std::vector<std::string> astrings;
        for(auto &a: imps)
        {
            auto ai = a->atomIndices();
            astrings.push_back(gmx::formatString("%d-%d-%d-%d", ai[0], ai[1], ai[2], ai[3]));
        }
        myCheck.checkSequence(astrings.begin(), astrings.end(), "Proper");
    }

    void testInteraction(InteractionType                     itype,
                         const std::vector<TopologyEntry *> &entry)
    {
        gmx::test::TestReferenceChecker myCheck(this->rootChecker());
        Topology top(bonds_);
        top.addEntry(itype, entry);
        auto &imps = top.entry(itype);
        myCheck.checkInteger(imps.size(), "#custom");
        std::vector<std::string> astrings;
        for(auto &a: imps)
        {
            auto ai = a->atomIndices();
            astrings.push_back(gmx::formatString("%d-%d", ai[0], ai[1]));
        }
        myCheck.checkSequence(astrings.begin(), astrings.end(), interactionTypeToString(itype).c_str());
        
    }
    //! Clean the test data.
    static void TearDownTestCase()
    {
    }
    
};

TEST_F (TopologyTest, OneBond){
    Bond b(2, 3, 2.0);
    EXPECT_TRUE(b.aI() == 2);
    EXPECT_FALSE(b.aJ() == 2);
    EXPECT_TRUE(b.bondOrder() == 2.0);
}

TEST_F (TopologyTest, FindBond){
    Topology top(bonds_);
    for(const auto &b : bonds_)
    {
        auto bb = top.findBond(b.aI(), b.aJ());
        EXPECT_TRUE(bb->bondOrder() == b.bondOrder());
    }
}

TEST_F (TopologyTest, DontFindBond){
    Topology top(bonds_);
    EXPECT_THROW((void) top.findBond(8, 9), gmx::InternalError);
}

TEST_F (TopologyTest, OneAngle) {
    Bond bij(1, 2, 1);
    Bond bjk(2, 3, 2);
    Bond bkl(3, 4, 3);
    Angle a(bij, bjk);
    EXPECT_FALSE(a.isLinear());
    EXPECT_THROW(Angle(bij, bkl), gmx::InvalidInputError);
}

TEST_F (TopologyTest, MakeAngles){
    testLinearAngles(x_);
}

TEST_F (TopologyTest, MakeLinearAngles){
    testLinearAngles(y_);
}

TEST_F (TopologyTest, MakeImpropers){
    testImpropers(5.0);
}

TEST_F (TopologyTest, MakePropers){
    testPropers(x_);
}

TEST_F (TopologyTest, MakePropersLinearAngles){
    testPropers(y_);
}

TEST_F (TopologyTest, AddPolarization){
    std::vector<TopologyEntry *> pols = {
        new Bond(1, 2, 1.0),
        new Bond(3, 4, 1.0)
    };
    testInteraction(InteractionType::POLARIZATION, pols); 
}

} // namespace alexandria
