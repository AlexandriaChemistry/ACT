/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022-2026
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <algorithm>
#include <fstream>
#include <gtest/gtest.h>
#include <string>
#include <utility>
#include <vector>

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/import/import.h"
#include "act/molprop/experiment.h"
#include "act/molprop/molprop.h"
#include "act/testutils/forcefield_utils.h"
#include "act/topology/actmol.h"
#include "gromacs/utility/stringutil.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/*! \brief The atom types and the bond orders of a molecule, each sorted
 * \param[in] molprops The molecule
 * \return the atom types, and the bond orders as text
 */
std::pair<std::vector<std::string>, std::vector<std::string>>
typesAndBondOrders(const std::vector<MolProp> &molprops)
{
    std::vector<std::string> atypes, orders;
    EXPECT_EQ(1u, molprops.size());
    for (const auto &mp : molprops)
    {
        for (const auto &ca : mp.experimentConst()[0].calcAtomConst())
        {
            atypes.push_back(ca.getObtype());
        }
        for (const auto &bond : mp.bondsConst())
        {
            orders.push_back(gmx::formatString("%g", bond.bondOrder()));
        }
    }
    std::sort(atypes.begin(), atypes.end());
    std::sort(orders.begin(), orders.end());
    return { atypes, orders };
}

/*! \brief The number of times each pair of atoms is in the CONECT records
 * \param[in] file The PDB file
 * \return for each number of times, how many pairs have it
 */
std::vector<int> conectMultiplicities(const std::string &file)
{
    std::ifstream    in(file);
    std::vector<int> count(4, 0);
    for (std::string line; std::getline(in, line); )
    {
        if (line.compare(0, 6, "CONECT") != 0)
        {
            continue;
        }
        // The first partner and how often it is listed
        const auto partners = gmx::splitString(line.substr(11));
        const auto times    = std::count(partners.begin(), partners.end(), partners.front());
        count[std::min<size_t>(times, 3)] += 1;
    }
    return count;
}

class PdbWriterTest : public ::testing::Test
{
protected:
    //! For the files that are written
    gmx::test::TestFileManager fileManager_;
};

/* Test: A molecule written to a PDB file is read back as the same
 * molecule: double and triple bonds are listed two and three times in
 * the CONECT records, and an aromatic ring is written as a Kekule
 * structure. Before, each bond was listed once, and a reader made
 * pyridine piperidine and the oxygen of acetone a hydroxyl.
 */
TEST_F(PdbWriterTest, BondOrdersAreWritten)
{
    auto pd = getForceField("ACS-g1");
    // SMILES, the expected numbers of double and of triple bonds
    const std::vector<std::tuple<std::string, int, int>> molecules = {
        { "c1ccncc1", 3, 0 }, { "CC(=O)C", 1, 0 }, { "CC#N", 0, 1 }, { "c1ccccc1", 3, 0 }
    };
    for (const auto &[smiles, ndouble, ntriple] : molecules)
    {
        SCOPED_TRACE(smiles);
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        std::vector<MolProp> made;
        importSmiles(&msghandler, pd, smiles, &made, false);
        ASSERT_TRUE(msghandler.ok());
        ASSERT_EQ(1u, made.size());
        // Before the molecule takes them over
        const auto expected = typesAndBondOrders(made);
        ACTMol     mol;
        mol.Merge(&made[0]);
        mol.GenerateTopology(&msghandler, pd, missingParameters::Ignore);
        ASSERT_TRUE(msghandler.ok());
        const auto   file = fileManager_.getTemporaryFilePath(gmx::formatString("mol%zu.pdb", smiles.size()));
        const matrix box  = { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } };
        mol.PrintConformation(file.c_str(), mol.xOriginal(), false, box);

        const auto times = conectMultiplicities(file);
        EXPECT_EQ(ndouble, times[2]);
        EXPECT_EQ(ntriple, times[3]);

        std::vector<MolProp> read;
        double               qtot = 0;
        importFile(&msghandler, pd, file, &read, "minimum", JobType::OPT, false, &qtot, false);
        EXPECT_TRUE(msghandler.ok());
        EXPECT_EQ(expected, typesAndBondOrders(read));
    }
}

} // namespace

} // namespace alexandria
