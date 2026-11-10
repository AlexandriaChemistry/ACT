/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
 * Tests of the compact storage of the non-bonded pairs, see
 * act/topology/nonbondedpairs.h, against the topology entries.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cstring>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/molprop/molprop_xml.h"
#include "act/topology/actmol.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class NonBondedPairsTest : public gmx::test::CommandLineTestBase
{
protected:
    /*! \brief Compare the pairs of a range with the entries of a TopologyEntryRange
     * \param[in] what    What is compared, for the messages
     * \param[in] entries The entries
     * \param[in] pairs   The pairs
     * \param[in] typed   Whether the pairs have parameters and identifiers
     */
    static void compare(const std::string &what, TopologyEntryRange entries, const PairRange &pairs,
                        bool typed)
    {
        SCOPED_TRACE(what);
        ASSERT_EQ(entries.size(), pairs.size());
        size_t k = 0;
        for (const auto &e : entries)
        {
            auto p = pairs[k++];
            ASSERT_EQ(2u, e->atomIndices().size());
            EXPECT_EQ(e->atomIndices()[0], p.atomIndices()[0]) << k;
            EXPECT_EQ(e->atomIndices()[1], p.atomIndices()[1]) << k;
            EXPECT_EQ(e->atomIndex(1), p->atomIndex(1)) << k;
            if (typed)
            {
                // Bitwise, since some parameters are NaN
                const auto &ep = e->params();
                const auto &pp = p.params();
                ASSERT_EQ(ep.size(), pp.size()) << k;
                EXPECT_EQ(0, std::memcmp(ep.data(), pp.data(), ep.size()*sizeof(double))) << k;
                EXPECT_EQ(e->id().id(), p.id().id()) << k;
            }
            else
            {
                EXPECT_TRUE(p.params().empty()) << k;
            }
        }
    }

    /*! \brief Compare the pairs of all non-bonded types of a topology with its entries
     * \param[in] what What topology this is, for the messages
     * \param[in] pd   The force field
     * \param[in] top  The topology
     * \return the number of pairs compared
     */
    static size_t compareTopology(const std::string &what, const ForceField &pd, const Topology &top)
    {
        size_t      ncompared = 0;
        const auto &nb        = top.nonBondedPairs();
        for (auto itype : { InteractionType::VDW, InteractionType::ELECTROSTATICS,
                            InteractionType::VDWCORRECTION, InteractionType::INDUCTIONCORRECTION })
        {
            auto name = what + " " + interactionTypeToString(itype);
            EXPECT_EQ(top.hasEntry(itype), nb.hasType(itype)) << name;
            if (!top.hasEntry(itype) || !nb.hasType(itype))
            {
                continue;
            }
            const auto &entries = top.entry(itype);
            bool        typed   = pd.interactionPresent(itype);
            EXPECT_EQ(typed, nb.pairTypes(itype) != nullptr) << name;
            compare(name + " all", TopologyEntryRange(entries), nb.range(itype), typed);
            compare(name + " intermolecular", top.intermolecular(itype, entries), nb.intermolecular(itype), typed);
            for (bool shells : { false, true })
            {
                compare(name + (shells ? " with shells" : " without shells"),
                        top.shellEntries(itype, entries, shells), nb.shellEntries(itype, shells), typed);
            }
            ncompared += entries.size();
        }
        return ncompared;
    }

    /*! \brief Make the topologies of the molecules of a file and compare their pairs
     * \param[in] ffFile The force field
     * \param[in] mpFile The molecules
     * \param[in] convert Whether to store virtual sites on atoms as second charges
     */
    void compareMolecules(const std::string &ffFile, const std::string &mpFile, bool convert = true)
    {
        ForceField pd;
        readForceField(fileManager().getInputFilePath(ffFile), &pd, nullptr, convert);
        std::vector<MolProp> mps;
        MsgHandler           msghandler;
        MolPropRead(&msghandler, fileManager().getInputFilePath(mpFile).c_str(), &mps);
        ASSERT_FALSE(mps.empty());
        size_t ncompared = 0;
        for (auto &mp : mps)
        {
            ACTMol actmol;
            actmol.Merge(&mp);
            MsgHandler mh;
            actmol.GenerateTopology(&mh, &pd, missingParameters::Ignore);
            if (!mh.ok())
            {
                fprintf(stderr, "%s: no topology for %s\n", ffFile.c_str(), mp.getMolname().c_str());
                continue;
            }
            ncompared += compareTopology(mp.getMolname(), pd, *actmol.topology());
            size_t f = 0;
            for (const auto &top : actmol.fragmentHandler()->topologies())
            {
                ncompared += compareTopology(mp.getMolname() + " fragment " + std::to_string(f++), pd, top);
            }
        }
        EXPECT_GT(ncompared, 0u);
    }
};

// Test: Without shells, with Gaussian charges
TEST_F(NonBondedPairsTest, DimersNoShells)
{
    compareMolecules("ACS-g2.xml", "mpdimer.xml");
}

// Test: Without shells, with virtual sites on the atoms
TEST_F(NonBondedPairsTest, DimersVsites)
{
    compareMolecules("../../constraints/tests/PC+GV-elec.xml", "mpdimer.xml");
}

// Test: With shells on ions and water
TEST_F(NonBondedPairsTest, DimersShells)
{
    compareMolecules("PG-pol-ions.xml", "ion-water.xml");
}

// Test: With virtual sites on the ions and shells
TEST_F(NonBondedPairsTest, DimersVsitesShells)
{
    compareMolecules("PG-pol-ions.xml", "ion-water.xml", false);
}

// Test: The pairs take much less memory than the entries
TEST_F(NonBondedPairsTest, Memory)
{
    ForceField pd;
    readForceField(fileManager().getInputFilePath("ACS-g2.xml"), &pd);
    std::vector<MolProp> mps;
    MsgHandler           msghandler;
    MolPropRead(&msghandler, fileManager().getInputFilePath("mpdimer.xml").c_str(), &mps);
    ASSERT_FALSE(mps.empty());
    ACTMol actmol;
    actmol.Merge(&mps[0]);
    actmol.GenerateTopology(&msghandler, &pd, missingParameters::Ignore);
    const auto &nb = actmol.topology()->nonBondedPairs();
    ASSERT_GT(nb.npair(), 0u);
    // Two particle indices per pair, a position and a type pair per pair
    // and non-bonded type, and the type pairs
    EXPECT_LT(nb.memory(), 32*nb.npair() + 4096);
}

} // namespace

} // namespace alexandria
