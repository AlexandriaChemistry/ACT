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
/*! \internal \file
 * \brief
 * Tests of making molecules from SMILES strings, see importSmiles.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/import/import.h"
#include "act/molprop/molprop.h"
#include "act/testutils/forcefield_utils.h"
#include "gromacs/math/vectypes.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! What is compared between a molecule from a file and one from SMILES
struct MoleculeSummary
{
    //! The atom types, sorted
    std::vector<std::string> atomTypes;
    //! The bonds as "type-type order" with the types in alphabetical order, sorted
    std::vector<std::string> bonds;
    //! The total charge
    int                      charge = 0;
    //! The coordinates in nm
    std::vector<gmx::RVec>   x;
    //! The bonds as pairs of atom numbers
    std::vector<std::pair<int, int>> pairs;
};

/*! \brief Collect what is compared from a molprop
 * \param[in] mp The molprop
 * \return the summary
 */
MoleculeSummary summarize(MolProp *mp)
{
    MoleculeSummary summary;
    if (mp->experimentConst().empty())
    {
        return summary;
    }
    const auto              &exper = mp->experimentConst()[0];
    std::vector<std::string> types;
    for (const auto &ca : exper.calcAtomConst())
    {
        types.push_back(ca.getObtype());
    }
    summary.x = exper.getCoordinates();
    for (const auto &bond : mp->bondsConst())
    {
        auto ti = types[bond.aI()];
        auto tj = types[bond.aJ()];
        if (tj < ti)
        {
            std::swap(ti, tj);
        }
        summary.bonds.push_back(ti + "-" + tj + " " + std::to_string(bond.bondOrder()));
        summary.pairs.emplace_back(bond.aI(), bond.aJ());
    }
    summary.atomTypes = types;
    std::sort(summary.atomTypes.begin(), summary.atomTypes.end());
    std::sort(summary.bonds.begin(), summary.bonds.end());
    summary.charge = mp->totalCharge();
    return summary;
}

class SmilesTest : public ::testing::Test
{
protected:
    //! The force field with the atom types
    ForceField *pd_ = nullptr;

    //! Get the force field
    void SetUp() override
    {
        pd_ = getForceField("gaff_hydrogen_types");
    }

    /*! \brief Make a molecule from SMILES
     * \param[in]  smiles The SMILES string
     * \param[out] mps    The molecule
     * \return whether there was no error
     */
    bool fromSmiles(const std::string &smiles, std::vector<MolProp> *mps)
    {
        MsgHandler msghandler;
        // No messages on the screen for strings that are wrong on purpose
        msghandler.setPrintLevel(ACTStatus::Fatal);
        importSmiles(&msghandler, pd_, smiles, mps, false);
        return msghandler.ok();
    }

    /*! \brief Compare a molecule from SMILES with the one in a file
     * \param[in] smiles  The SMILES string
     * \param[in] molFile File in src/act/testdata/molecules
     */
    void compareWithFile(const std::string &smiles, const std::string &molFile)
    {
        SCOPED_TRACE(smiles + " and " + molFile);
        std::vector<MolProp> fromString;
        ASSERT_TRUE(fromSmiles(smiles, &fromString));
        ASSERT_EQ(1U, fromString.size());

        std::vector<MolProp> fromFile;
        MsgHandler           msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        double               qtot = 0;
        importFile(&msghandler, pd_,
                   gmx::test::TestFileManager::getInputFilePath("molecules/" + molFile),
                   &fromFile, "minimum", JobType::OPT, false, &qtot, false);
        ASSERT_TRUE(msghandler.ok());
        ASSERT_EQ(1U, fromFile.size());

        const auto s = summarize(&fromString[0]);
        const auto f = summarize(&fromFile[0]);
        EXPECT_EQ(f.atomTypes, s.atomTypes);
        EXPECT_EQ(f.bonds, s.bonds);
        EXPECT_EQ(f.charge, s.charge);
        EXPECT_EQ(fromFile[0].formula(), fromString[0].formula());
        checkGeometry(s);
    }

    /*! \brief Check that the coordinates are those of a real molecule
     * \param[in] s The molecule
     */
    static void checkGeometry(const MoleculeSummary &s)
    {
        ASSERT_EQ(s.atomTypes.size(), s.x.size());
        const auto distance = [&s](int i, int j)
            {
                double d2 = 0;
                for (int m = 0; m < DIM; m++)
                {
                    d2 += (s.x[i][m] - s.x[j][m])*(s.x[i][m] - s.x[j][m]);
                }
                return std::sqrt(d2);
            };
        // Bonds have the length of chemical bonds (nm)
        for (const auto &p : s.pairs)
        {
            EXPECT_GT(distance(p.first, p.second), 0.085);
            EXPECT_LT(distance(p.first, p.second), 0.2);
        }
        // No two atoms are on top of each other
        for (size_t i = 0; i < s.x.size(); i++)
        {
            for (size_t j = i+1; j < s.x.size(); j++)
            {
                EXPECT_GT(distance(i, j), 0.08) << "atoms " << i << " and " << j;
            }
        }
    }
};

TEST_F(SmilesTest, SameMoleculeAsFile)
{
    // The hydrogen atoms are not in the strings, they are added
    compareWithFile("CC(=O)C", "acetone.sdf");
    compareWithFile("O", "water.sdf");
    compareWithFile("CC=O", "acetaldehyde.sdf");
    compareWithFile("CC#N", "acetonitrile.sdf");
    compareWithFile("N", "ammonia.sdf");
    compareWithFile("ClC(Cl)Cl", "chloroform.sdf");
    compareWithFile("ClCCl", "dichloromethane.sdf");
    compareWithFile("NC=O", "formamide.sdf");
    compareWithFile("CN", "methylamine.sdf");
    compareWithFile("Cl", "hydrogen-chloride.sdf");
}

TEST_F(SmilesTest, AromaticRings)
{
    compareWithFile("c1ccoc1", "furan.sdf");
    compareWithFile("c1ccncc1", "pyridine.sdf");
    compareWithFile("c1ccsc1", "thiophene.sdf");
    compareWithFile("O=c1cc[nH]c(=O)[nH]1", "uracil.sdf");
}

TEST_F(SmilesTest, Ion)
{
    compareWithFile("CC(=O)[O-]", "acetate.sdf");
    std::vector<MolProp> mps;
    ASSERT_TRUE(fromSmiles("[NH4+]", &mps));
    ASSERT_EQ(1U, mps.size());
    EXPECT_EQ(1, mps[0].totalCharge());
    EXPECT_EQ(5U, summarize(&mps[0]).atomTypes.size());
}

TEST_F(SmilesTest, ThreeDimensional)
{
    // Ethanol is not flat: the coordinates span all three directions
    std::vector<MolProp> mps;
    ASSERT_TRUE(fromSmiles("CCO", &mps));
    ASSERT_EQ(1U, mps.size());
    const auto s = summarize(&mps[0]);
    ASSERT_EQ(9U, s.x.size());
    for (int m = 0; m < DIM; m++)
    {
        double lo = s.x[0][m];
        double hi = s.x[0][m];
        for (const auto &x : s.x)
        {
            lo = std::min<double>(lo, x[m]);
            hi = std::max<double>(hi, x[m]);
        }
        EXPECT_GT(hi-lo, 0.05) << "direction " << m;
    }
    checkGeometry(s);
}

TEST_F(SmilesTest, SameCoordinatesEveryTime)
{
    // The coordinates are made with a fixed seed, such that results
    // with a SMILES string as input can be repeated.
    std::vector<MolProp> first, second;
    ASSERT_TRUE(fromSmiles("CC(C)Cc1ccccc1", &first));
    ASSERT_TRUE(fromSmiles("CC(C)Cc1ccccc1", &second));
    const auto a = summarize(&first[0]);
    const auto b = summarize(&second[0]);
    ASSERT_EQ(a.x.size(), b.x.size());
    for (size_t i = 0; i < a.x.size(); i++)
    {
        for (int m = 0; m < DIM; m++)
        {
            EXPECT_EQ(a.x[i][m], b.x[i][m]);
        }
    }
}

TEST_F(SmilesTest, WrongStringsGiveErrors)
{
    // No exception and no molecule, but an error for the user
    for (const std::string smiles : { "", "   ", "C(C", "c1ccccc", "XYZ", "C=#C", "not a molecule",
                                      "C1CC", "[Zz]", "C[C@@", "CC(=O)C)" })
    {
        std::vector<MolProp> mps;
        bool                 ok = true;
        EXPECT_NO_THROW(ok = fromSmiles(smiles, &mps)) << "'" << smiles << "'";
        EXPECT_FALSE(ok) << "'" << smiles << "'";
        EXPECT_TRUE(mps.empty()) << "'" << smiles << "'";
    }
}

TEST_F(SmilesTest, VeryLongStringIsRefused)
{
    // A limit protects a server against strings that take minutes
    const std::string smiles(20000, 'C');
    std::vector<MolProp> mps;
    EXPECT_FALSE(fromSmiles(smiles, &mps));
    EXPECT_TRUE(mps.empty());
}

} // namespace

} // namespace alexandria
