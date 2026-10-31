/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2021-2023
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
#include <cmath>

#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/allbondeds.h"
#include "act/alexandria/molselect.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/molprop/molprop_xml.h"
#include "act/topology/actmol.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class AllBondedsTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
 
        AllBondedsTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
            //ForceField *mypd     = getForceField("ACM-g");
        }

        void runOBstats(const std::vector<double> &xx)
        {
            Identifier atom("h");
            OneBonded ob(atom);
            for(const auto &x : xx)
            {
                ob.addPoint(x);
            }
            real   average, sigma;
            size_t N;
            eStats ok = ob.getAverageSigmaN(&average, &sigma, &N);
            if (eStats::OK == ok)
            {
                checker_.checkReal(average, "average");
                checker_.checkReal(sigma, "sigma");
                checker_.checkInteger(static_cast<int>(N), "N");
            }
            else
            {
                checker_.checkString(gmx_stats_message(ok), "Error");
            }
            
        }
        static void TearDownTestCase()
        {
        }

};

TEST_F (AllBondedsTest, OneBondedStats){
    std::vector<double> xx = { 103, 104, 105, 106.5, 107, 109, 106, 104, 107 };
    runOBstats(xx);
}

TEST_F (AllBondedsTest, OneBondedStatsEmpty){
    std::vector<double> xx;
    runOBstats(xx);
}

/* Test: The geometries are extracted as in geometry_ff for molecules with
 * virtual sites. The topologies are made without parameters, such that
 * the positions of the virtual sites cannot be constructed, see
 * https://github.com/AlexandriaChemistry/ACT/issues/506
 */
TEST(AllBondedsVsiteTest, ExtractGeometries)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    gmx::test::TestFileManager fileManager;
    ForceField pd;
    readForceField(fileManager.getInputFilePath("../../forces/tests/ACS-pg-vs3.xml"), &pd);
    // As in geometry_ff, such that there are no bonds to shells
    pd.setPolarizable(false);
    std::vector<MolProp> mps;
    MolPropRead(&msghandler, fileManager.getInputFilePath("../../molprop/tests/molprop.xml").c_str(), &mps);
    auto  selName = fileManager.getTemporaryFilePath("selection.dat");
    FILE *fp      = fopen(selName.c_str(), "w");
    for(const char *mol : { "1-butanol", "acetone", "water" })
    {
        fprintf(fp, "%s|Train\n", mol);
    }
    fclose(fp);
    MolSelect gms;
    gms.read(selName.c_str());

    AllBondeds          bondeds;
    std::vector<ACTMol> actmols;
    EXPECT_NO_THROW(bondeds.extractGeometries(&msghandler, &mps, &actmols, &pd, gms));
    EXPECT_TRUE(msghandler.ok());
    ASSERT_EQ(3u, actmols.size());
    // The molecules have virtual sites without parameters
    size_t nvsite = 0;
    for(auto &mol : actmols)
    {
        size_t molvsite = 0;
        for(const auto &entry : mol.topology()->entries())
        {
            if (isVsite(entry.first))
            {
                for(const auto &vs : entry.second)
                {
                    EXPECT_TRUE(vs->params().empty());
                    molvsite += 1;
                }
            }
        }
        nvsite += molvsite;
        // Their positions cannot be constructed, without them the atoms
        // are at their positions in the experiment.
        auto xexp = mol.experimentConst()[0].getCoordinates();
        if (molvsite > 0)
        {
            EXPECT_THROW(mol.experCoords(xexp), gmx::InternalError) << mol.getMolname();
        }
        std::vector<gmx::RVec> x;
        EXPECT_NO_THROW(x = mol.experCoords(xexp, false)) << mol.getMolname();
        size_t j = 0;
        for(size_t i = 0; i < x.size(); i++)
        {
            if (ActParticle::Atom == mol.topology()->atoms()[i].pType())
            {
                for(int m = 0; m < DIM; m++)
                {
                    EXPECT_EQ(xexp[j][m], x[i][m]);
                }
                j += 1;
            }
        }
        EXPECT_EQ(xexp.size(), j);
    }
    EXPECT_GT(nvsite, 0u);
    // Bonds and angles were found
    gmx::StringOutputStream stream;
    gmx::TextWriter         tw(&stream);
    bondeds.writeSummary(&tw);
    auto summary = stream.toString();
    EXPECT_NE(std::string::npos, summary.find("BONDS")) << summary;
    EXPECT_NE(std::string::npos, summary.find("ANGLES")) << summary;
}

}

} // namespace
