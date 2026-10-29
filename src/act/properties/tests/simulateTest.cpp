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
 * Tests for the alexandria simulate tool.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/alex_modules.h"
#include "act/basics/msg_handler.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "gromacs/math/units.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! Contents of an xvg file
struct XvgData
{
    //! Legends of the data columns, not including time
    std::vector<std::string>         legends;
    //! Rows of data, the first column is the time
    std::vector<std::vector<double>> rows;

    //! \return the column with this legend, or -1 if not present
    int column(const std::string &legend) const
    {
        for(size_t i = 0; i < legends.size(); i++)
        {
            if (legends[i] == legend)
            {
                return i+1;
            }
        }
        return -1;
    }

    //! \return all values in a column
    std::vector<double> values(int col) const
    {
        std::vector<double> vals;
        for(const auto &row : rows)
        {
            vals.push_back(row[col]);
        }
        return vals;
    }
};

//! Read an xvg file
XvgData readXvg(const std::string &filename)
{
    XvgData       xvg;
    std::ifstream in(filename);
    EXPECT_TRUE(in.good()) << "Cannot open " << filename;
    std::string   line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        if (line[0] == '@')
        {
            // Legends look like: @ s0 legend "EPOT"
            auto pos = line.find("legend \"");
            if (line.compare(0, 3, "@ s") == 0 && pos != std::string::npos)
            {
                auto start = pos + 8;
                xvg.legends.push_back(line.substr(start, line.rfind('"') - start));
            }
            continue;
        }
        std::istringstream  iss(line);
        std::vector<double> row;
        double              value;
        while (iss >> value)
        {
            row.push_back(value);
        }
        if (!row.empty())
        {
            xvg.rows.push_back(row);
        }
    }
    return xvg;
}

//! Read the atom coordinates (Angstrom) from a pdb file
std::vector<std::vector<double>> readPdbCoordinates(const std::string &filename)
{
    std::vector<std::vector<double>> x;
    std::ifstream in(filename);
    EXPECT_TRUE(in.good()) << "Cannot open " << filename;
    std::string   line;
    while (std::getline(in, line))
    {
        if (line.compare(0, 4, "ATOM") == 0 || line.compare(0, 6, "HETATM") == 0)
        {
            x.push_back({ std::stod(line.substr(30, 8)),
                          std::stod(line.substr(38, 8)),
                          std::stod(line.substr(46, 8)) });
        }
    }
    return x;
}

//! \return whether a file exists and is not empty
bool fileNotEmpty(const std::string &filename)
{
    std::ifstream in(filename);
    return in.good() && in.peek() != std::ifstream::traits_type::eof();
}

//! Output file names for one run
struct SimulateFiles
{
    //! Energy file
    std::string energy;
    //! Trajectory file
    std::string trajectory;
    //! Log file
    std::string log;
};

class SimulateTest : public gmx::test::CommandLineTestBase
{
protected:
    //! \return path to an input file in the alexandria tests
    std::string alexandriaFile(const std::string &name)
    {
        return fileManager().getInputFilePath("../../alexandria/tests/" + name);
    }

    /*! \brief Make a command line for simulate
     * \param[in]  molFile Molecule file in the alexandria/tests/mols directory
     * \param[in]  label   Unique label for the output files of this run
     * \param[out] files   The output file names
     * \return The command line
     */
    gmx::test::CommandLine baseCommand(const std::string &molFile,
                                       const std::string &label,
                                       SimulateFiles     *files)
    {
        gmx::test::CommandLine cmd;
        cmd.append("simulate");
        cmd.addOption("-ff", alexandriaFile("ACS-g.xml"));
        cmd.addOption("-f", alexandriaFile("mols/" + molFile));
        files->energy     = fileManager().getTemporaryFilePath(label + "-energy.xvg");
        files->trajectory = fileManager().getTemporaryFilePath(label + "-traj.pdb");
        files->log        = fileManager().getTemporaryFilePath(label + ".log");
        cmd.addOption("-e", files->energy);
        cmd.addOption("-o", files->trajectory);
        cmd.addOption("-g", files->log);
        return cmd;
    }

    //! Add options for a 100 step MD simulation
    static void addMD(gmx::test::CommandLine *cmd, double temperature, int seed)
    {
        cmd->addOption("-nsteps", 100);
        cmd->addOption("-deltat", 0.0002);
        cmd->addOption("-temp", temperature);
        // A seed of zero would mean a random seed
        cmd->addOption("-seed", seed);
        cmd->addOption("-nstener", 10);
        cmd->addOption("-nstxout", 50);
    }

    //! Run the simulate tool
    static int run(gmx::test::CommandLine *cmd)
    {
        return simulate(cmd->argc(), cmd->argv());
    }

    //! Compare all energy columns with reference data
    void checkEnergies(const XvgData &xvg)
    {
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-6));
        checker.checkSequence(xvg.legends.begin(), xvg.legends.end(), "Legends");
        auto time = xvg.values(0);
        checker.checkSequence(time.begin(), time.end(), "Time");
        for(size_t i = 0; i < xvg.legends.size(); i++)
        {
            auto vals = xvg.values(i+1);
            checker.checkSequence(vals.begin(), vals.end(), xvg.legends[i].c_str());
        }
    }

    //! Check that the total energy is conserved
    static void checkConservation(const XvgData &xvg)
    {
        int cEtot = xvg.column("ETOT");
        int cEkin = xvg.column("EKIN");
        ASSERT_GT(cEtot, 0);
        ASSERT_GT(cEkin, 0);
        double etot0   = xvg.rows[0][cEtot];
        double drift   = 0;
        double ekinAve = 0;
        for(const auto &row : xvg.rows)
        {
            drift    = std::max(drift, std::abs(row[cEtot] - etot0));
            ekinAve += row[cEkin];
        }
        ekinAve /= xvg.rows.size();
        ASSERT_GT(ekinAve, 0);
        // The fluctuation in the total energy should be small compared to
        // the kinetic energy for a short time step.
        EXPECT_LT(drift, 0.02*ekinAve) << "Total energy drift " << drift << " average kinetic energy " << ekinAve;
    }

    //! Run 100 steps of MD and check the energies
    void testMD(const std::string &molFile)
    {
        SimulateFiles files;
        auto cmd = baseCommand(molFile, "md", &files);
        addMD(&cmd, 300, 1993);
        ASSERT_EQ(0, run(&cmd));
        EXPECT_TRUE(fileNotEmpty(files.trajectory));
        auto xvg = readXvg(files.energy);
        // Energies are written every 10 steps
        ASSERT_EQ(10u, xvg.rows.size());
        checkEnergies(xvg);
    }
};

TEST_F(SimulateTest, MDHydrogenChloride)
{
    testMD("hydrogen-chloride.sdf");
}

TEST_F(SimulateTest, MDCarbonDioxide)
{
    testMD("carbon-dioxide.sdf");
}

TEST_F(SimulateTest, MDWater)
{
    testMD("water-3-oep.log.pdb");
}

TEST_F(SimulateTest, MDMethanethiol)
{
    testMD("methanethiol.sdf");
}

TEST_F(SimulateTest, MDAcetone)
{
    testMD("acetone.sdf");
}

TEST_F(SimulateTest, MDUracil)
{
    testMD("uracil.sdf");
}

TEST_F(SimulateTest, MDWaterDimer)
{
    testMD("water_dimer.pdb");
}

// Test: The total energy is conserved
TEST_F(SimulateTest, EnergyConservation)
{
    for(const auto &molFile : { "hydrogen-chloride.sdf", "carbon-dioxide.sdf",
                                "water-3-oep.log.pdb", "methanethiol.sdf",
                                "acetone.sdf", "uracil.sdf", "water_dimer.pdb" })
    {
        SCOPED_TRACE(molFile);
        SimulateFiles files;
        auto cmd = baseCommand(molFile, gmx::formatString("cons-%s", molFile), &files);
        addMD(&cmd, 300, 1993);
        ASSERT_EQ(0, run(&cmd));
        checkConservation(readXvg(files.energy));
    }
}

// Test: The same seed gives the same simulation, a different seed does not
TEST_F(SimulateTest, Reproducibility)
{
    std::vector<XvgData> xvg;
    int                  run_id = 0;
    for(int seed : { 1993, 1993, 1997 })
    {
        SimulateFiles files;
        auto cmd = baseCommand("acetone.sdf", gmx::formatString("seed%d", run_id++), &files);
        addMD(&cmd, 300, seed);
        ASSERT_EQ(0, run(&cmd));
        xvg.push_back(readXvg(files.energy));
    }
    ASSERT_EQ(xvg[0].rows.size(), xvg[1].rows.size());
    for(size_t i = 0; i < xvg[0].rows.size(); i++)
    {
        EXPECT_EQ(xvg[0].rows[i], xvg[1].rows[i]) << "row " << i;
    }
    // The initial kinetic energy is scaled to the temperature, so
    // a different seed shows up only later in the trajectory.
    int cEkin = xvg[0].column("EKIN");
    ASSERT_GT(cEkin, 0);
    ASSERT_EQ(xvg[0].rows.size(), xvg[2].rows.size());
    auto last = xvg[0].rows.size() - 1;
    EXPECT_EQ(xvg[0].rows[0][cEkin], xvg[2].rows[0][cEkin]);
    EXPECT_NE(xvg[0].rows[last][cEkin], xvg[2].rows[last][cEkin]);
}

// Test: The initial temperature is the requested one, computed for
// the degrees of freedom without center of mass motion and constraints
TEST_F(SimulateTest, InitialTemperature)
{
    // Acetone has 10 atoms and 9 bonds
    for(bool constrain : { false, true })
    {
        SimulateFiles files;
        auto cmd = baseCommand("acetone.sdf",
                               constrain ? "tempcons" : "temp", &files);
        addMD(&cmd, 300, 1993);
        if (constrain)
        {
            cmd.append("-constrain");
        }
        ASSERT_EQ(0, run(&cmd));
        auto xvg   = readXvg(files.energy);
        ASSERT_FALSE(xvg.rows.empty());
        int  cEkin = xvg.column("EKIN");
        int  cTemp = xvg.column("TEMPERATURE");
        ASSERT_GT(cEkin, 0);
        ASSERT_GT(cTemp, 0);
        auto ekin  = xvg.rows[0][cEkin];
        auto temp  = xvg.rows[0][cTemp];
        EXPECT_NEAR(300, temp, 1e-4);
        double nDOF = constrain ? 3*9-9 : 3*9;
        EXPECT_NEAR(nDOF, 2*ekin/(BOLTZ*temp), 1e-3);
    }
}

// Test: Without an initial temperature no velocities are generated
TEST_F(SimulateTest, ZeroTemperature)
{
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "zero", &files);
    addMD(&cmd, 0, 1993);
    ASSERT_EQ(0, run(&cmd));
    auto xvg  = readXvg(files.energy);
    ASSERT_FALSE(xvg.rows.empty());
    int cEkin = xvg.column("EKIN");
    ASSERT_GT(cEkin, 0);
    EXPECT_EQ(0, xvg.rows[0][cEkin]);
}

// Test: Minimize before running MD and write the final structure
TEST_F(SimulateTest, MinimizeThenMD)
{
    SimulateFiles files;
    auto cmd     = baseCommand("acetone-3-oep.log.pdb", "minmd", &files);
    auto confout = fileManager().getTemporaryFilePath("minmd-confout.pdb");
    cmd.addOption("-minimize");
    cmd.addOption("-c", confout);
    addMD(&cmd, 300, 1993);
    ASSERT_EQ(0, run(&cmd));
    EXPECT_TRUE(fileNotEmpty(confout));
    EXPECT_TRUE(fileNotEmpty(files.trajectory));
    auto xvg = readXvg(files.energy);
    EXPECT_EQ(10u, xvg.rows.size());

    // Without minimization the simulation starts at a higher energy
    SimulateFiles files2;
    auto cmd2 = baseCommand("acetone-3-oep.log.pdb", "nomin", &files2);
    addMD(&cmd2, 300, 1993);
    ASSERT_EQ(0, run(&cmd2));
    auto xvg2  = readXvg(files2.energy);
    int  cEpot = xvg.column("EPOT");
    ASSERT_GT(cEpot, 0);
    ASSERT_FALSE(xvg.rows.empty());
    ASSERT_FALSE(xvg2.rows.empty());
    EXPECT_LT(xvg.rows[0][cEpot], xvg2.rows[0][cEpot]);
}

// Test: Starting at a minimum without velocities, the molecule stays at rest
TEST_F(SimulateTest, MinimizeThenMDAtRest)
{
    SimulateFiles files;
    auto cmd = baseCommand("acetone-3-oep.log.pdb", "rest", &files);
    cmd.addOption("-minimize");
    addMD(&cmd, 0, 1993);
    ASSERT_EQ(0, run(&cmd));
    auto xvg   = readXvg(files.energy);
    int  cEkin = xvg.column("EKIN");
    int  cEpot = xvg.column("EPOT");
    ASSERT_GT(cEkin, 0);
    ASSERT_GT(cEpot, 0);
    ASSERT_FALSE(xvg.rows.empty());
    for(const auto &row : xvg.rows)
    {
        EXPECT_LT(row[cEkin], 1e-3) << "time " << row[0];
        EXPECT_NEAR(xvg.rows[0][cEpot], row[cEpot], 1e-3) << "time " << row[0];
    }
}

// Test: A single point calculation without MD
TEST_F(SimulateTest, SinglePoint)
{
    for(const auto &molFile : { "acetone.sdf", "water_dimer.pdb" })
    {
        SimulateFiles files;
        auto cmd = baseCommand(molFile, gmx::formatString("sp-%s", molFile), &files);
        cmd.addOption("-sp");
        cmd.addOption("-nsteps", 0);
        EXPECT_EQ(0, run(&cmd)) << molFile;
        EXPECT_FALSE(fileNotEmpty(files.trajectory)) << molFile;
    }
}

// Test: Minimization of a dimer, including the interaction energy
TEST_F(SimulateTest, MinimizeDimer)
{
    SimulateFiles files;
    auto cmd     = baseCommand("water_dimer.pdb", "mindimer", &files);
    auto confout = fileManager().getTemporaryFilePath("mindimer-confout.pdb");
    cmd.addOption("-minimize");
    cmd.addOption("-c", confout);
    EXPECT_EQ(0, run(&cmd));
    EXPECT_TRUE(fileNotEmpty(confout));
}

// Test: Frozen atoms do not move during minimization
TEST_F(SimulateTest, MinimizeFreeze)
{
    std::string   molFile("acetone-3-oep.log.pdb");
    SimulateFiles files;
    auto cmd     = baseCommand(molFile, "freeze", &files);
    auto confout = fileManager().getTemporaryFilePath("freeze-confout.pdb");
    auto ndx     = fileManager().getTemporaryFilePath("freeze.ndx");
    {
        std::ofstream out(ndx);
        out << "[ Frozen ]\n1 2\n";
    }
    cmd.addOption("-minimize");
    cmd.addOption("-c", confout);
    cmd.addOption("-freeze", ndx);
    ASSERT_EQ(0, run(&cmd));
    auto xin  = readPdbCoordinates(alexandriaFile("mols/" + molFile));
    auto xout = readPdbCoordinates(confout);
    ASSERT_EQ(xin.size(), xout.size());
    // Coordinates are written with three decimals
    double toler = 1.5e-3;
    for(size_t i = 0; i < 2; i++)
    {
        for(int m = 0; m < 3; m++)
        {
            EXPECT_NEAR(xin[i][m], xout[i][m], toler) << "atom " << i+1;
        }
    }
    // The other atoms should have moved, otherwise nothing was tested
    double maxMove = 0;
    for(size_t i = 2; i < xin.size(); i++)
    {
        for(int m = 0; m < 3; m++)
        {
            maxMove = std::max(maxMove, std::abs(xin[i][m] - xout[i][m]));
        }
    }
    EXPECT_GT(maxMove, toler);
}

// Test: MD with flat-bottomed position restraints
TEST_F(SimulateTest, PositionRestraints)
{
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "fbpr", &files);
    addMD(&cmd, 300, 1993);
    cmd.addOption("-kfbpr", 1000.0);
    cmd.addOption("-r0fbpr", 0.0);
    ASSERT_EQ(0, run(&cmd));
    auto xvg = readXvg(files.energy);
    bool found = false;
    for(const auto &legend : xvg.legends)
    {
        if (legend.find("RESTRAINT") != std::string::npos)
        {
            found = true;
        }
    }
    EXPECT_TRUE(found);
}

// Test: Writing the molprop file
TEST_F(SimulateTest, WriteMolprop)
{
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "omp", &files);
    auto omp = fileManager().getTemporaryFilePath("omp.xml");
    cmd.addOption("-sp");
    cmd.addOption("-omp", omp);
    ASSERT_EQ(0, run(&cmd));
    MsgHandler           msghandler;
    std::vector<MolProp> mps;
    MolPropRead(&msghandler, omp.c_str(), &mps);
    ASSERT_EQ(1u, mps.size());
    // The atoms are stored in the experiment
    ASSERT_FALSE(mps[0].experimentConst().empty());
    EXPECT_EQ(10, mps[0].experimentConst()[0].NAtom());
    ASSERT_EQ(1u, mps[0].fragments().size());
    EXPECT_EQ("C3H6O", mps[0].fragments()[0].formula());
}

// Test: A molecule file that does not exist
TEST_F(SimulateTest, MissingMolecule)
{
    SimulateFiles files;
    auto cmd = baseCommand("does-not-exist.sdf", "missing", &files);
    cmd.addOption("-sp");
    EXPECT_EQ(1, run(&cmd));
}

} // namespace

} // namespace alexandria
