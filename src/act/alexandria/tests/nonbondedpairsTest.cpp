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
 * act/topology/nonbondedpairs.h, against their definition.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cstring>
#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/potential.h"
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
    /*! \brief Check the non-bonded pairs of a topology against their definition
     * All pairs of particles that are not excluded, each once, the
     * intermolecular pairs and those with shells last if separated, the
     * same pairs for all non-bonded types of which the force field has
     * the parameters, and the parameters of the pair of particle types.
     * \param[in] what   What topology this is, for the messages
     * \param[in] pd     The force field
     * \param[in] top    The topology
     * \param[in] natom0 The number of particles of the first fragment of a dimer, 0 otherwise
     * \return the number of pairs checked
     */
    static size_t checkTopology(const std::string &what, const ForceField &pd, const Topology &top,
                                size_t natom0)
    {
        SCOPED_TRACE(what);
        const auto &nb    = top.nonBondedPairs();
        const auto &atoms = top.atoms();
        // The pairs that are not excluded
        std::set<std::pair<int, int>> expected;
        const auto &excl = top.exclusions(InteractionType::VDW);
        for (size_t i = 0; i < atoms.size(); i++)
        {
            for (size_t j = i+1; j < atoms.size(); j++)
            {
                if (std::find(excl[i].begin(), excl[i].end(), static_cast<int>(j)) == excl[i].end())
                {
                    expected.insert({ static_cast<int>(i), static_cast<int>(j) });
                }
            }
        }
        EXPECT_EQ(expected.size(), nb.npair());
        if (expected.empty())
        {
            return 0;
        }
        // Coulomb has all of them, each once
        const auto coulomb = nb.range(InteractionType::ELECTROSTATICS);
        std::set<std::pair<int, int>> found;
        for (const auto &p : coulomb)
        {
            EXPECT_TRUE(found.insert({ p.atomIndex(0), p.atomIndex(1) }).second);
            EXPECT_TRUE(p.params().empty());
        }
        EXPECT_EQ(expected, found);
        const auto isShell = [&atoms](const PairView &p)
        {
            return (atoms[p.atomIndex(0)].pType() == ActParticle::Shell ||
                    atoms[p.atomIndex(1)].pType() == ActParticle::Shell);
        };
        const auto isInter = [natom0](const PairView &p)
        {
            return ((static_cast<size_t>(p.atomIndex(0)) < natom0) != (static_cast<size_t>(p.atomIndex(1)) < natom0));
        };
        size_t nchecked = 0;
        for (auto itype : { InteractionType::ELECTROSTATICS, InteractionType::VDW,
                            InteractionType::VDWCORRECTION, InteractionType::INDUCTIONCORRECTION })
        {
            if (!nb.hasType(itype))
            {
                continue;
            }
            SCOPED_TRACE(interactionTypeToString(itype));
            const auto all = nb.range(itype);
            // The separations: without shells first, intermolecular last
            const auto noShells   = nb.shellEntries(itype, false);
            const auto withShells = nb.shellEntries(itype, true);
            EXPECT_EQ(all.size(), noShells.size() + withShells.size());
            if (top.hasShellSplit())
            {
                for (const auto &p : noShells)
                {
                    EXPECT_FALSE(isShell(p));
                }
                for (const auto &p : withShells)
                {
                    EXPECT_TRUE(isShell(p));
                }
            }
            if (top.hasIntermolecular())
            {
                const auto inter = nb.intermolecular(itype);
                size_t ninter = 0;
                for (const auto &p : all)
                {
                    ninter += isInter(p) ? 1 : 0;
                }
                EXPECT_EQ(ninter, inter.size());
                for (const auto &p : inter)
                {
                    EXPECT_TRUE(isInter(p));
                }
            }
            if (InteractionType::ELECTROSTATICS == itype)
            {
                continue;
            }
            // The pairs that the force field has parameters for, in the
            // order of the Coulomb pairs, with those parameters
            EXPECT_NE(nullptr, nb.pairTypes(itype));
            const auto &fs     = pd.findForcesConst(itype);
            const auto &names  = potprops.find(fs.potential())->second.param;
            size_t      k      = 0;
            for (const auto &c : coulomb)
            {
                if (k < all.size() && all[k].atomIndices() == c.atomIndices())
                {
                    const auto &p  = all[k++];
                    const auto  ff = fs.parametersConst().find(p.id());
                    EXPECT_TRUE(fs.parametersConst().end() != ff) << p.id().id();
                    EXPECT_EQ(names.size(), p.params().size());
                    if (fs.parametersConst().end() == ff || names.size() != p.params().size())
                    {
                        continue;
                    }
                    for (size_t i = 0; i < names.size(); i++)
                    {
                        const auto fp    = ff->second.find(names[i]);
                        double value = (ff->second.end() == fp) ? 0 : fp->second.internalValue();
                        // Bitwise, since some parameters are NaN
                        EXPECT_EQ(0, std::memcmp(&value, &p.params()[i], sizeof(double))) << p.id().id() << " " << names[i];
                    }
                }
                // Otherwise the pair was left out, since the force field
                // has no parameters for its type pair
            }
            EXPECT_EQ(all.size(), k);
            nchecked += k;
        }
        return nchecked;
    }

    /*! \brief Make the topologies of the molecules of a file and check their pairs
     * \param[in] ffFile  The force field
     * \param[in] mpFile  The molecules
     * \param[in] convert Whether to store virtual sites on atoms as second charges
     */
    void checkMolecules(const std::string &ffFile, const std::string &mpFile, bool convert = true)
    {
        ForceField pd;
        readForceField(fileManager().getInputFilePath(ffFile), &pd, nullptr, convert);
        std::vector<MolProp> mps;
        MsgHandler           msghandler;
        MolPropRead(&msghandler, fileManager().getInputFilePath(mpFile).c_str(), &mps);
        ASSERT_FALSE(mps.empty());
        size_t nchecked = 0;
        for (auto &mp : mps)
        {
            ACTMol actmol;
            actmol.Merge(&mp);
            MsgHandler mh;
            actmol.GenerateTopology(&mh, &pd, missingParameters::Ignore);
            if (!mh.ok())
            {
                continue;
            }
            const auto &start = actmol.fragmentHandler()->atomStart();
            size_t natom0 = (actmol.isDimer() && start.size() == 2) ? start[1] : 0;
            nchecked += checkTopology(mp.getMolname(), pd, *actmol.topology(), natom0);
            size_t f = 0;
            for (const auto &top : actmol.fragmentHandler()->topologies())
            {
                nchecked += checkTopology(mp.getMolname() + " fragment " + std::to_string(f++), pd, top, 0);
            }
        }
        EXPECT_GT(nchecked, 0u);
    }
};

// Test: Without shells, with Gaussian charges
TEST_F(NonBondedPairsTest, DimersNoShells)
{
    checkMolecules("forcefields/ACS-g2.xml", "molprop/mpdimer.xml");
}

// Test: Without shells, with virtual sites on the atoms
TEST_F(NonBondedPairsTest, DimersVsites)
{
    checkMolecules("forcefields/PC+GV-elec.xml", "molprop/mpdimer.xml");
}

// Test: With shells on ions and water
TEST_F(NonBondedPairsTest, DimersShells)
{
    checkMolecules("forcefields/PG-pol-ions.xml", "molprop/ion-water.xml");
}

// Test: With virtual sites on the ions and shells
TEST_F(NonBondedPairsTest, DimersVsitesShells)
{
    checkMolecules("forcefields/PG-pol-ions.xml", "molprop/ion-water.xml", false);
}

// Test: The pairs take much less memory than the entries
TEST_F(NonBondedPairsTest, Memory)
{
    ForceField pd;
    readForceField(fileManager().getInputFilePath("forcefields/ACS-g2.xml"), &pd);
    std::vector<MolProp> mps;
    MsgHandler           msghandler;
    MolPropRead(&msghandler, fileManager().getInputFilePath("molprop/mpdimer.xml").c_str(), &mps);
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
