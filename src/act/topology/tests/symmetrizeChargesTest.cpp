/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024-2026
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
 * Tests for charge symmetrization functions.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <cstddef>
#include <gtest/gtest.h>
#include <map>
#include <string>
#include <vector>

#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/import/import.h"
#include "act/molprop/molprop.h"
#include "act/topology/actmol.h"
#include "act/topology/symmetrize_charges.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/testfilemanager.h"


namespace alexandria
{

namespace
{

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesNoSymmetry)
{
    // Each atom in its own group (no symmetry)
    std::vector<double> q = { 0.1, 0.2, 0.3 };
    const std::vector<int> sym  = { 0, 1, 2 };
    apply_symmetrized_charges(&q, sym);
    EXPECT_NEAR(q[0], 0.1, 1e-10);
    EXPECT_NEAR(q[1], 0.2, 1e-10);
    EXPECT_NEAR(q[2], 0.3, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesAllSame)
{
    // All atoms in same group
    std::vector<double> q = { 0.1, 0.2, 0.3 };
    const std::vector<int> sym  = { 0, 0, 0 };
    apply_symmetrized_charges(&q, sym);
    double expected = (0.1 + 0.2 + 0.3) / 3.0;
    EXPECT_NEAR(q[0], expected, 1e-10);
    EXPECT_NEAR(q[1], expected, 1e-10);
    EXPECT_NEAR(q[2], expected, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesPartialGroups)
{
    // Water-like: O by itself, two H's symmetrized
    std::vector<double> q = { -0.8, 0.35, 0.45 };
    const std::vector<int> sym  = { 0, 1, 1 };
    apply_symmetrized_charges(&q, sym);
    EXPECT_NEAR(q[0], -0.8, 1e-10);
    double hAvg = (0.35 + 0.45) / 2.0;
    EXPECT_NEAR(q[1], hAvg, 1e-10);
    EXPECT_NEAR(q[2], hAvg, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesMethane)
{
    // Methane-like: C and 4 H's
    std::vector<double> q = { -0.4, 0.08, 0.12, 0.10, 0.10 };
    const std::vector<int> sym  = { 0, 1, 1, 1, 1 };
    apply_symmetrized_charges(&q, sym);
    EXPECT_NEAR(q[0], -0.4, 1e-10);
    double hAvg = (0.08 + 0.12 + 0.10 + 0.10) / 4.0;
    EXPECT_NEAR(q[1], hAvg, 1e-10);
    EXPECT_NEAR(q[2], hAvg, 1e-10);
    EXPECT_NEAR(q[3], hAvg, 1e-10);
    EXPECT_NEAR(q[4], hAvg, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesMultipleGroups)
{
    // Two separate symmetric groups
    std::vector<double> q = { 0.1, 0.2, 0.5, 0.3, 0.4 };
    const std::vector<int> sym  = { 0, 0, 2, 3, 3 };
    apply_symmetrized_charges(&q, sym);
    double avg01 = (0.1 + 0.2) / 2.0;
    EXPECT_NEAR(q[0], avg01, 1e-10);
    EXPECT_NEAR(q[1], avg01, 1e-10);
    EXPECT_NEAR(q[2], 0.5, 1e-10);
    double avg34 = (0.3 + 0.4) / 2.0;
    EXPECT_NEAR(q[3], avg34, 1e-10);
    EXPECT_NEAR(q[4], avg34, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesSizeMismatch)
{
    // Size mismatch should throw
    std::vector<double> q = { 0.1, 0.2, 0.3 };
    const std::vector<int> sym  = { 0, 1 };
    EXPECT_THROW(apply_symmetrized_charges(&q, sym), gmx::InternalError);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesPreservesTotal)
{
    // Total charge should be preserved after symmetrization
    std::vector<double> q = { -0.82, 0.39, 0.43 };
    const std::vector<int> sym  = { 0, 1, 1 };
    double totalBefore = q[0] + q[1] + q[2];
    apply_symmetrized_charges(&q, sym);
    double totalAfter = q[0] + q[1] + q[2];
    EXPECT_NEAR(totalAfter, totalBefore, 1e-10);
}

/*! \brief The symmetry of the charges of 1-butanol
 * \param[in]  methylene Whether the force field has the rule for a carbon
 *                       atom with two hydrogen atoms
 * \param[out] elements  The elements of the atoms
 * \param[out] carbon    For each atom the carbon atom it is bound to, -1
 *                       if there is none or if it is a carbon atom itself
 * \param[out] sym       For each atom the atom of which it has the charge
 */
void butanolSymmetry(bool                      methylene,
                     std::vector<std::string> *elements,
                     std::vector<int>         *carbon,
                     std::vector<int>         *sym)
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACS-g1.xml"), &pd);
    if (methylene)
    {
        pd.addSymcharges("C", "H", 2);
    }
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    std::vector<MolProp> mps;
    double               qtot = 0;
    importFile(&msghandler, &pd, gmx::test::TestFileManager::getInputFilePath("molecules/1-butanol.pdb"),
               &mps, "minimum", JobType::OPT, false, &qtot, true);
    ASSERT_TRUE(msghandler.ok());
    ASSERT_EQ(1u, mps.size());
    ACTMol mol;
    mol.Merge(&mps[0]);
    mol.GenerateTopology(&msghandler, &pd, missingParameters::Ignore);
    ASSERT_TRUE(msghandler.ok());
    const auto &atoms = mol.topology()->atoms();
    for (const auto &atom : atoms)
    {
        elements->push_back(atom.element());
    }
    carbon->assign(atoms.size(), -1);
    ASSERT_TRUE(mol.topology()->hasEntry(InteractionType::BONDS));
    for (const auto &bond : mol.topology()->entry(InteractionType::BONDS))
    {
        const int ai = bond->atomIndex(0);
        const int aj = bond->atomIndex(1);
        if ("C" == (*elements)[ai] && "C" != (*elements)[aj])
        {
            (*carbon)[aj] = ai;
        }
        else if ("C" == (*elements)[aj] && "C" != (*elements)[ai])
        {
            (*carbon)[ai] = aj;
        }
    }
    get_symmetrized_charges(mol.topology(), &pd, nullptr, sym);
    ASSERT_EQ(atoms.size(), sym->size());
}

// Test: With the rule for a carbon atom with two hydrogen atoms, the
// hydrogen atoms of each methylene group of 1-butanol get the same charge,
// as those of the methyl group, which have it without that rule too.
// OpenMM gives the atoms of a residue the charges of its template by the
// bonds, not by the names of the atoms, such that atoms that cannot be told
// apart from the bonds should have the same charge.
TEST(SymmetrizeChargesTest, MethyleneHydrogens)
{
    for (bool methylene : { false, true })
    {
        std::vector<std::string> elements;
        std::vector<int>         carbon, sym;
        butanolSymmetry(methylene, &elements, &carbon, &sym);
        if (::testing::Test::HasFatalFailure())
        {
            return;
        }
        // The hydrogen atoms on each carbon atom
        std::map<int, std::vector<int>> hydrogens;
        for (size_t i = 0; i < elements.size(); i++)
        {
            if ("H" == elements[i] && carbon[i] >= 0)
            {
                hydrogens[carbon[i]].push_back(i);
            }
            else
            {
                // Other atoms keep their own charge
                EXPECT_EQ(static_cast<int>(i), sym[i]) << elements[i] << " " << i;
            }
        }
        int nMethyl = 0, nMethylene = 0;
        for (const auto &ch : hydrogens)
        {
            const auto &hs = ch.second;
            ASSERT_TRUE(hs.size() == 2 || hs.size() == 3);
            const bool same = (3 == hs.size() || methylene);
            (3 == hs.size() ? nMethyl : nMethylene)++;
            for (const auto h : hs)
            {
                EXPECT_EQ(same ? hs[0] : h, sym[h]) << "carbon " << ch.first << " hydrogen " << h
                                                    << " methylene rule " << methylene;
            }
        }
        EXPECT_EQ(1, nMethyl);
        EXPECT_EQ(3, nMethylene);
    }
}

} // namespace

} // namespace alexandria
