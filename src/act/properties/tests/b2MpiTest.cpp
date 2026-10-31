/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
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
 * Tests for the alexandria b2 tool with a water dimer, using MPI.
 *
 * The tests can be run with any number of MPI ranks and the results
 * should not depend on it. The reference data is generated on a single
 * rank and the test is registered to run on four ranks.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/alex_modules.h"
#include "act/properties/secondvirial.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/stringutil.h"

#include <nlohmann/json.hpp>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! Wait for all MPI ranks
void barrier()
{
#if GMX_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//! Second virial components for one temperature from the B2T output
struct B2Row
{
    //! Temperature
    double T;
    //! Components in cm^3/mol
    double total, classical, force, torque1, torque2;
};

//! Read an xvg file with data columns
std::vector<std::vector<double>> readColumns(const std::string &filename)
{
    std::vector<std::vector<double>> rows;
    std::ifstream in(filename);
    EXPECT_TRUE(in.good()) << "Cannot open " << filename;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#' || line[0] == '@' || line[0] == '&')
        {
            continue;
        }
        std::istringstream  iss(line);
        std::vector<double> row;
        double              v;
        while (iss >> v)
        {
            row.push_back(v);
        }
        if (!row.empty())
        {
            rows.push_back(row);
        }
    }
    return rows;
}

// Test fixture running the b2 tool for a water dimer
class WaterB2Test : public gmx::test::CommandLineTestBase
{
protected:
    /*! \brief Run the b2 tool for the water dimer
     * \param[in]  ffName   Force field file in this directory
     * \param[in]  nDimer   Number of dimer orientations
     * \param[in]  sampling Options for sampling the orientations
     * \param[out] mayer    The Mayer functions, first column distance
     * \return The B2 components for 400 and 500 K
     */
    std::vector<B2Row> runB2(const std::string                                &ffName,
                             int                                               nDimer,
                             const std::vector<std::pair<std::string, int>>   &sampling,
                             std::vector<std::vector<double>>                 *mayer)
    {
        gmx::test::CommandLine cmd;
        cmd.append("b2");
        cmd.addOption("-ff", fileManager().getInputFilePath(ffName));
        cmd.addOption("-f", fileManager().getInputFilePath("../../alexandria/tests/mols/water_dimer.pdb"));
        // A short calculation: few distances and orientations
        cmd.addOption("-mindist", 0.1);
        cmd.addOption("-maxdist", 1.2);
        cmd.addOption("-ndist", 111);
        cmd.addOption("-maxdimer", nDimer);
        // Sobol sampling is deterministic and converges fast, unless a
        // grid is used
        if (sampling.empty())
        {
            cmd.addOption("-randalg", "Sobol");
        }
        else
        {
            cmd.addOption("-randalg", "Gauss");
            for(const auto &s : sampling)
            {
                cmd.addOption(s.first.c_str(), s.second);
            }
        }
        cmd.addOption("-T1", 400.0);
        cmd.addOption("-T2", 500.0);
        cmd.addOption("-dT", 100.0);
        auto label     = gmx::formatString("%s-%d", ffName.c_str(), nDimer);
        auto b2file    = fileManager().getTemporaryFilePath(label + "-B2T.xvg");
        auto mayerfile = fileManager().getTemporaryFilePath(label + "-mayer.xvg");
        cmd.addOption("-b2", b2file);
        cmd.addOption("-eh", mayerfile);
        cmd.addOption("-g", fileManager().getTemporaryFilePath(label + "-b2.log"));
        EXPECT_EQ(0, b2(cmd.argc(), cmd.argv()));
        std::vector<B2Row> rows;
        if (gmx_node_rank() != 0)
        {
            // Only the master writes the output
            return rows;
        }
        for(const auto &r : readColumns(b2file))
        {
            EXPECT_EQ(6u, r.size());
            if (r.size() == 6)
            {
                rows.push_back({ r[0], r[1], r[2], r[3], r[4], r[5] });
            }
        }
        *mayer = readColumns(mayerfile);
        return rows;
    }

    /*! \brief Run the b2 tool and check the results on the master
     * The temporary files are removed at the end of the test by each
     * rank, therefore all ranks wait until the master has checked the
     * results.
     */
    void testWater(const std::string &ffName, int nDimer, bool checkTorqueSymmetry,
                   const std::vector<std::pair<std::string, int>> &sampling = {})
    {
        std::vector<std::vector<double>> mayer;
        auto rows = runB2(ffName, nDimer, sampling, &mayer);
        if (gmx_node_rank() == 0)
        {
            checkWater(rows, mayer, checkTorqueSymmetry);
        }
        barrier();
    }

    /*! \brief Compare the second virial coefficients with reference data
     * A short calculation, with all digits of the results from the JSON
     * output. With a tight tolerance, any change of the energies, forces,
     * torques, virtual sites, shells or sums of the second virial is
     * found, such that optimizations do not change the results unnoticed.
     * \param[in] ffName  Force field file in this directory
     * \param[in] molFile Dimer structure, relative to this directory
     * \param[in] nDimer  Number of dimer orientations
     */
    void testReference(const std::string &ffName, const std::string &molFile, int nDimer)
    {
        gmx::test::CommandLine cmd;
        cmd.append("b2");
        cmd.addOption("-ff", fileManager().getInputFilePath(ffName));
        cmd.addOption("-f", fileManager().getInputFilePath(molFile));
        cmd.addOption("-mindist", 0.15);
        cmd.addOption("-maxdist", 1.2);
        cmd.addOption("-ndist", 106);
        cmd.addOption("-maxdimer", nDimer);
        cmd.addOption("-randalg", "Sobol");
        cmd.addOption("-T1", 300.0);
        cmd.addOption("-T2", 500.0);
        cmd.addOption("-dT", 100.0);
        cmd.addOption("-shelltoler", 1e-8);
        cmd.addOption("-shellMaxIter", 100);
        cmd.append("-json");
        auto label = gmx::formatString("%s-%d", ffName.c_str(), nDimer);
        cmd.addOption("-b2", fileManager().getTemporaryFilePath(label + "-B2T.xvg"));
        cmd.addOption("-eh", fileManager().getTemporaryFilePath(label + "-mayer.xvg"));
        cmd.addOption("-g", fileManager().getTemporaryFilePath(label + "-b2.log"));
        EXPECT_EQ(0, b2(cmd.argc(), cmd.argv()));
        if (gmx_node_rank() == 0)
        {
            // The b2 tool writes b2.json in the working directory
            const char    *jsonFile = "b2.json";
            std::ifstream  in(jsonFile);
            ASSERT_TRUE(in.good()) << "Cannot open " << jsonFile;
            auto json = nlohmann::json::parse(in);
            in.close();
            std::remove(jsonFile);
            auto &b2json = json["SecondVirialCoefficient"]["Second virial coefficients"];
            ASSERT_TRUE(b2json.is_object());
            EXPECT_EQ("cm^3/mol", b2json["Unit"]);
            auto checker = rootChecker();
#if GMX_DOUBLE
            checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-10));
#else
            checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-4));
#endif
            size_t ntemp = 0;
            for(const auto &t : b2json.items())
            {
                if (t.key() == "Unit")
                {
                    continue;
                }
                ntemp += 1;
                for(const auto &b : t.value().items())
                {
                    checker.checkReal(std::stod(b.value().get<std::string>()),
                                      gmx::formatString("%s %s", t.key().c_str(), b.key().c_str()).c_str());
                }
            }
            EXPECT_EQ(3u, ntemp);
        }
        barrier();
    }
    //! Check the results of the b2 tool
    void checkWater(const std::vector<B2Row>               &rows,
                    const std::vector<std::vector<double>> &mayer,
                    bool                                    checkTorqueSymmetry)
    {
        ASSERT_EQ(2u, rows.size());
        EXPECT_EQ(400, rows[0].T);
        EXPECT_EQ(500, rows[1].T);
        for(const auto &r : rows)
        {
            SCOPED_TRACE(gmx::formatString("T = %g", r.T));
            // The total is the classical part plus the quantum corrections,
            // where the torque correction is the average over the two molecules.
            double sum = r.classical + r.force + 0.5*(r.torque1 + r.torque2);
            EXPECT_NEAR(sum, r.total, 1e-4*std::abs(r.total) + 0.01);
            // Water attracts
            EXPECT_LT(r.classical, 0);
            // Quantum corrections are positive and smaller than the classical value
            EXPECT_GT(r.force, 0);
            EXPECT_GT(r.torque1, 0);
            EXPECT_GT(r.torque2, 0);
            EXPECT_LT(r.force + 0.5*(r.torque1 + r.torque2), std::abs(r.classical));
            if (checkTorqueSymmetry)
            {
                // For two identical molecules, the torque corrections of
                // molecule 1 and 2 should converge to the same value.
                EXPECT_NEAR(r.torque1, r.torque2, 0.02*0.5*(r.torque1 + r.torque2));
            }
        }
        // The attraction and the quantum corrections decrease with temperature
        EXPECT_LT(rows[0].classical, rows[1].classical);
        EXPECT_GT(rows[0].force, rows[1].force);
        EXPECT_GT(rows[0].torque1, rows[1].torque1);
        EXPECT_GT(rows[0].torque2, rows[1].torque2);

        // Mayer function: -1 where molecules overlap, zero far away
        ASSERT_FALSE(mayer.empty());
        ASSERT_EQ(3u, mayer.front().size());
        for(int iT = 1; iT <= 2; iT++)
        {
            EXPECT_DOUBLE_EQ(-1, mayer.front()[iT]);
            EXPECT_GE(mayer.back()[iT], -0.05);
            EXPECT_LE(mayer.back()[iT], 0.05);
        }

        // Reference data, the xvg file has six significant digits
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-5));
        for(const auto &r : rows)
        {
            auto tc = checker.checkCompound("Temperature", gmx::formatString("%g", r.T).c_str());
            tc.checkDouble(r.classical, b2TypeToString(b2Type::Classical).c_str());
            tc.checkDouble(r.force, b2TypeToString(b2Type::Force).c_str());
            tc.checkDouble(r.torque1, b2TypeToString(b2Type::Torque1).c_str());
            tc.checkDouble(r.torque2, b2TypeToString(b2Type::Torque2).c_str());
            tc.checkDouble(r.total, b2TypeToString(b2Type::Total).c_str());
        }
    }
};

TEST_F(WaterB2Test, SPCE)
{
    testWater("SPCE.xml", 3200, true);
}

// Test: The number of dimers does not have to be a multiple of the number of ranks
TEST_F(WaterB2Test, SPCEUnevenDimers)
{
    testWater("SPCE.xml", 3203, true);
}

// Test: A grid of orientations as in the program Orient, -maxdimer is not used
TEST_F(WaterB2Test, SPCEGauss)
{
    testWater("SPCE.xml", 1024, true, { { "-nalpha", 9 }, { "-nbeta", 9 }, { "-ngamma", 9 } });
}

TEST_F(WaterB2Test, SWM4NDP)
{
    // For a polarizable model this also checks that the forces of the
    // fully relaxed dimer are used for the torques.
    testWater("SWM4-NDP.xml", 3200, true);
}

// Test: Second virial coefficients compared with reference data, for
// water models with virtual sites, shells, point and Gaussian charges and
// for hydrogen fluoride with a 14-7 potential. These cover the energies,
// forces and torques of the dimers, the virtual sites, the relaxation of
// the shells and the sums over the orientations, also over MPI processes.
TEST_F(WaterB2Test, ReferenceTIP5P)
{
    testReference("TIP5P.xml", "../../alexandria/tests/mols/water_dimer.pdb", 64);
}

TEST_F(WaterB2Test, ReferenceSWM3)
{
    testReference("SWM3.xml", "../../alexandria/tests/mols/water_dimer.pdb", 64);
}

TEST_F(WaterB2Test, ReferenceSWM4NDP)
{
    testReference("SWM4-NDP.xml", "../../alexandria/tests/mols/water_dimer.pdb", 64);
}

TEST_F(WaterB2Test, ReferenceBK3)
{
    testReference("BK3.xml", "../../alexandria/tests/mols/water_dimer.pdb", 64);
}

TEST_F(WaterB2Test, ReferenceSPCE)
{
    testReference("SPCE.xml", "../../alexandria/tests/mols/water_dimer.pdb", 64);
}

TEST_F(WaterB2Test, ReferenceHF14_7)
{
    testReference("HF-LJ14-7.xml", "hydrogen-fluoride#hydrogen-fluoride.pdb", 64);
}

} // namespace

} // namespace alexandria
