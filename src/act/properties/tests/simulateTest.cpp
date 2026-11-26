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
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <map>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/basics/msg_handler.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "act/properties/simulate.h"
#include "act/testutils/output_checks.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
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

//! Read all frames of a pdb trajectory, the coordinates in Angstrom
std::vector<std::vector<std::vector<double>>> readPdbFrames(const std::string &filename)
{
    std::vector<std::vector<std::vector<double>>> frames(1);
    std::ifstream in(filename);
    EXPECT_TRUE(in.good()) << "Cannot open " << filename;
    std::string   line;
    while (std::getline(in, line))
    {
        if (line.compare(0, 4, "ATOM") == 0 || line.compare(0, 6, "HETATM") == 0)
        {
            frames.back().push_back({ std::stod(line.substr(30, 8)),
                                      std::stod(line.substr(38, 8)),
                                      std::stod(line.substr(46, 8)) });
        }
        else if (line.compare(0, 6, "ENDMDL") == 0 && !frames.back().empty())
        {
            frames.emplace_back();
        }
    }
    if (frames.back().empty())
    {
        frames.pop_back();
    }
    return frames;
}

//! \return whether a file exists and is not empty
bool fileNotEmpty(const std::string &filename)
{
    std::ifstream in(filename);
    return in.good() && in.peek() != std::ifstream::traits_type::eof();
}

//! \return whether a file has a text
bool fileContains(const std::string &filename, const std::string &text)
{
    std::ifstream in(filename);
    for(std::string line; std::getline(in, line); )
    {
        if (line.find(text) != std::string::npos)
        {
            return true;
        }
    }
    return false;
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
    //! \return path to an input file in src/act/testdata
    std::string alexandriaFile(const std::string &name)
    {
        return fileManager().getInputFilePath(name);
    }

    /*! \brief Make a command line for simulate
     * \param[in]  molFile Molecule file in src/act/testdata/molecules, or an absolute path
     * \param[in]  label   Unique label for the output files of this run
     * \param[out] files   The output file names
     * \param[in]  ffFile  Force field file, forcefields/ACS-g1.xml from src/act/testdata if empty
     * \return The command line
     */
    gmx::test::CommandLine baseCommand(const std::string &molFile,
                                       const std::string &label,
                                       SimulateFiles     *files,
                                       const std::string &ffFile = "")
    {
        gmx::test::CommandLine cmd;
        cmd.append("simulate");
        cmd.addOption("-ff", ffFile.empty() ? alexandriaFile("forcefields/ACS-g1.xml") : ffFile);
        cmd.addOption("-f", molFile[0] == '/' ? molFile : alexandriaFile("molecules/" + molFile));
        files->energy     = fileManager().getTemporaryFilePath(label + "-energy.xvg");
        files->trajectory = fileManager().getTemporaryFilePath(label + "-traj.pdb");
        files->log        = fileManager().getTemporaryFilePath(label + ".log");
        cmd.addOption("-e", files->energy);
        cmd.addOption("-o", files->trajectory);
        cmd.addOption("-g", files->log);
        return cmd;
    }

    //! \return the potential energy of each frame in the log file of a rerun
    static std::vector<double> rerunEpot(const std::string &log)
    {
        std::vector<double> epot;
        std::ifstream       in(log);
        for (std::string line; std::getline(in, line); )
        {
            std::istringstream iss(line);
            std::string        word;
            int                index;
            if (!(iss >> index))
            {
                continue;
            }
            while (iss >> word)
            {
                double value;
                if (word == "EPOT" && iss >> value)
                {
                    epot.push_back(value);
                }
            }
        }
        return epot;
    }

    /*! \brief An energy of a single point calculation from its log file
     * \param[in] log  The log file
     * \param[in] term The energy term
     * \return the energy
     */
    static double singlePointEnergy(const std::string &log, const std::string &term = "EPOT")
    {
        std::ifstream in(log);
        for (std::string line; std::getline(in, line); )
        {
            std::istringstream iss(line);
            std::string        word, valueWord;
            double             value;
            if (iss >> word >> valueWord >> value && word == term + ":")
            {
                return value;
            }
        }
        ADD_FAILURE() << "No energy " << term << " in " << log;
        return 0;
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

    /*! \brief Compare all energy columns with reference data
     * \param[in] xvg       The energies
     * \param[in] shellSlack Difference allowed in addition in single
     *                       precision, for the energies of shells that are
     *                       relaxed to within a tolerance only
     */
    void checkEnergies(const XvgData &xvg, double shellSlack = 0)
    {
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        // The energies are written with five decimals, a tighter
        // tolerance would not make a difference
        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-6));
        checker.checkSequence(xvg.legends.begin(), xvg.legends.end(), "Legends");
        auto time = xvg.values(0);
        checker.checkSequence(time.begin(), time.end(), "Time");
        for(size_t i = 0; i < xvg.legends.size(); i++)
        {
            auto vals = xvg.values(i+1);
            if (!GMX_DOUBLE)
            {
                // In single precision the errors of the forces, relative to
                // the largest terms, accumulate along the trajectory, and
                // the energies differ by at least one in the last decimal
                double vmax = 0;
                for(auto v : vals)
                {
                    vmax = std::max(vmax, std::abs(v));
                }
                double absolute = 2e-5 + 1e-5*vmax + shellSlack;
                checker.setDefaultTolerance(gmx::test::FloatingPointTolerance(absolute, absolute, 1e-6, 1e-6,
                                                                              UINT64_MAX, UINT64_MAX, false));
            }
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
        const auto xvg = readXvg(files.energy);
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
// The energies of simulations with options other than the defaults, as
// reference data, such that a change in how the options reach the work shows

TEST_F(SimulateTest, ResultsOptions)
{
    // Minimized first, with bond constraints and without removing the
    // rotation
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "results-options", &files);
    addMD(&cmd, 300, 1993);
    cmd.addOption("-minimize");
    cmd.addOption("-minalg", "Newton");
    cmd.addOption("-constrain", "bonds");
    cmd.addOption("-comm", "linear");
    cmd.addOption("-nstcomm", 10);
    ASSERT_EQ(0, run(&cmd));
    const auto xvg = readXvg(files.energy);
    ASSERT_EQ(10u, xvg.rows.size());
    checkEnergies(xvg);
}

TEST_F(SimulateTest, ResultsIntervals)
{
    // Another time step, temperature, seed and output intervals
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "results-intervals", &files);
    cmd.addOption("-nsteps", 60);
    cmd.addOption("-deltat", 0.0005);
    cmd.addOption("-temp", 200);
    cmd.addOption("-seed", 7);
    cmd.addOption("-nstener", 15);
    cmd.addOption("-nstxout", 30);
    ASSERT_EQ(0, run(&cmd));
    const auto xvg = readXvg(files.energy);
    ASSERT_EQ(4u, xvg.rows.size());
    checkEnergies(xvg);
}

TEST_F(SimulateTest, ResultsShells)
{
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "results-shells", &files,
                           alexandriaFile("forcefields/ACS-pg1.xml"));
    addMD(&cmd, 300, 1993);
    ASSERT_EQ(0, run(&cmd));
    const auto xvg = readXvg(files.energy);
    ASSERT_EQ(10u, xvg.rows.size());
    // In single precision the induction energies differ by about 1e-3 kJ/mol
    checkEnergies(xvg, 5e-3);
}

TEST_F(SimulateTest, ResultsRestrained)
{
    // Flat-bottom position restraints on a dimer
    SimulateFiles files;
    auto cmd = baseCommand("water_dimer.pdb", "results-restrained", &files);
    addMD(&cmd, 300, 1993);
    cmd.addOption("-kfbpr", 500);
    cmd.addOption("-r0fbpr", 0.05);
    ASSERT_EQ(0, run(&cmd));
    const auto xvg = readXvg(files.energy);
    ASSERT_EQ(10u, xvg.rows.size());
    checkEnergies(xvg);
}

// The simulation without a command line, as a program that uses ACT runs it

TEST_F(SimulateTest, WithoutCommandLineSameAsTheProgram)
{
    // Minimized first, then 100 steps from velocities at 300 K
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "direct-tool", &files);
    addMD(&cmd, 300, 1993);
    cmd.addOption("-minimize");
    ASSERT_EQ(0, run(&cmd));

    SimulateOptions options;
    options.forceFieldFile = alexandriaFile("forcefields/ACS-g1.xml");
    options.compounds.setFileName(alexandriaFile("molecules/acetone.sdf"));
    options.energyFile     = fileManager().getTemporaryFilePath("direct-energy.xvg");
    options.trajectoryFile = fileManager().getTemporaryFilePath("direct-traj.pdb");
    options.simulation.setSimulation(100, 0.0002, 300, 1993, 10, 50);
    options.simulation.setMinimize(true);
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    // Without an output environment and without the processes
    EXPECT_EQ(0, alexandria::runSimulate(&msghandler, options));
    EXPECT_EQ(textWithoutComments(files.energy), textWithoutComments(options.energyFile));
    EXPECT_EQ(textWithoutComments(files.trajectory), textWithoutComments(options.trajectoryFile));

    // Without any file
    SimulateOptions quiet = options;
    quiet.energyFile.clear();
    quiet.trajectoryFile.clear();
    EXPECT_EQ(0, alexandria::runSimulate(&msghandler, quiet));
}

TEST_F(SimulateTest, SettingsAreCheckedWithoutFiles)
{
    // For a molecule in memory, as in Python: no files are needed, and
    // wrong settings are refused with the message of the program
    EXPECT_NO_THROW(SimulateOptions().checkSettings());
    EXPECT_THROW(SimulateOptions().check(), gmx::InvalidInputError);
    {
        SimulateOptions o;
        o.shellTolerance = 0;
        try
        {
            o.checkSettings();
            ADD_FAILURE() << "Accepted a wrong -shelltoler";
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find("-shelltoler")) << ex.what();
        }
    }
    {
        SimulateOptions o;
        o.simulation.setSimulation(10, 0.0, 0, 1, 1, 1);
        try
        {
            o.checkSettings();
            ADD_FAILURE() << "Accepted a wrong -deltat";
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find("-deltat")) << ex.what();
        }
    }
    {
        SimulateOptions o;
        o.simulation.setSimulation(10, 0.001, -1, 1, 1, 1);
        try
        {
            o.checkSettings();
            ADD_FAILURE() << "Accepted a wrong -temp";
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find("-temp")) << ex.what();
        }
    }
    {
        // The generator of velocities takes the seed as unsigned
        SimulateOptions o;
        o.simulation.setSimulation(10, 0.001, 0, -1, 1, 1);
        try
        {
            o.checkSettings();
            ADD_FAILURE() << "Accepted a wrong -seed";
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find("-seed")) << ex.what();
        }
    }
    // Zero steps (no dynamics), zero iterations (until convergence) and a seed of zero (a new one) are fine
    {
        SimulateOptions o;
        o.simulation.setSimulation(0, 0.001, 0, 0, 1, 1);
        o.simulation.setMaxIter(0);
        EXPECT_NO_THROW(o.checkSettings());
    }
}

TEST_F(SimulateTest, WithoutCommandLineWrongOptions)
{
    const auto valid = [this]()
    {
        SimulateOptions options;
        options.forceFieldFile = alexandriaFile("forcefields/ACS-g1.xml");
        options.compounds.setFileName(alexandriaFile("molecules/acetone.sdf"));
        options.simulation.setSimulation(10, 0.0002, 300, 1993, 1, 1);
        return options;
    };
    std::vector<std::pair<SimulateOptions, std::string>> wrong;
    {
        auto o = valid();
        o.forceFieldFile.clear();
        wrong.emplace_back(o, "-ff");
    }
    {
        SimulateOptions o;
        o.forceFieldFile = alexandriaFile("forcefields/ACS-g1.xml");
        wrong.emplace_back(o, "-f");
    }
    {
        auto o = valid();
        o.shellTolerance = 0;
        wrong.emplace_back(o, "-shelltoler");
    }
    {
        auto o = valid();
        o.simulation.setSimulation(10, 0, 300, 1993, 1, 1);
        wrong.emplace_back(o, "-deltat");
    }
    {
        auto o = valid();
        o.simulation.setSimulation(-10, 0.0002, 300, 1993, 1, 1);
        wrong.emplace_back(o, "-nsteps");
    }
    {
        auto o = valid();
        o.simulation.setSimulation(10, 0.0002, -300, 1993, 1, 1);
        wrong.emplace_back(o, "-temp");
    }
    for (const auto &w : wrong)
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        try
        {
            alexandria::runSimulate(&msghandler, w.first);
            ADD_FAILURE() << "Accepted options with a wrong " << w.second;
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find(w.second)) << ex.what();
        }
    }
    // The program says the same and stops, before it writes anything
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "wrongdeltat", &files);
    cmd.addOption("-nsteps", 10);
    cmd.addOption("-deltat", 0);
    EXPECT_EQ(1, run(&cmd));
    EXPECT_FALSE(fileNotEmpty(files.energy));
}

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

// Test: MD of a water dimer with shells and virtual sites. The induction
// energy at the start is that of a single point calculation, and with rigid
// molecules, as the model is meant, the energy is conserved. The flexible
// model needs much shorter time steps.
TEST_F(SimulateTest, MDShells)
{
    const auto ffFile = fileManager().getInputFilePath("forcefields/SWM4-NDP.xml");
    SimulateFiles sp;
    auto cmdsp = baseCommand("water_dimer.pdb", "shells-sp", &sp, ffFile);
    cmdsp.append("-sp");
    ASSERT_EQ(0, run(&cmdsp));

    SimulateFiles rigid;
    auto cmdr = baseCommand("water_dimer.pdb", "shells-rigid", &rigid, ffFile);
    addMD(&cmdr, 300, 1993);
    cmdr.addOption("-constrain", "angles");
    ASSERT_EQ(0, run(&cmdr));
    auto xvgr = readXvg(rigid.energy);
    checkConservation(xvgr);
    int  cetot = xvgr.column("ETOT");
    ASSERT_GT(cetot, 0);
    for (const auto &row : xvgr.rows)
    {
        EXPECT_NEAR(xvgr.rows[0][cetot], row[cetot], 1e-2) << "time " << row[0];
    }

    SimulateFiles md;
    auto cmd = baseCommand("water_dimer.pdb", "shells-md", &md, ffFile);
    addMD(&cmd, 300, 1993);
    ASSERT_EQ(0, run(&cmd));
    auto xvg  = readXvg(md.energy);
    int  cind = xvg.column("INDUCTION");
    ASSERT_GT(cind, 0);
    ASSERT_FALSE(xvg.rows.empty());
    // The single point energy is printed with four decimals
    EXPECT_NEAR(singlePointEnergy(sp.log, "INDUCTION"), xvg.rows[0][cind], 1e-3);
    for (const auto &row : xvg.rows)
    {
        EXPECT_LT(row[cind], 0) << "time " << row[0];
    }
}

// Test: Rerunning the trajectory of a simulation of the water dimer gives
// the potential energies of the simulation, for the frames written
TEST_F(SimulateTest, Rerun)
{
    SimulateFiles md;
    auto cmd = baseCommand("water_dimer.pdb", "md", &md);
    addMD(&cmd, 300, 1993);
    ASSERT_EQ(0, run(&cmd));
    auto xvg   = readXvg(md.energy);
    int  cEpot = xvg.column("EPOT");
    ASSERT_GT(cEpot, 0);

    SimulateFiles rerun;
    auto cmd2 = baseCommand("water_dimer.pdb", "rerun", &rerun);
    cmd2.addOption("-traj", md.trajectory);
    ASSERT_EQ(0, run(&cmd2));
    // Each structure gives a line with its number and the energies
    auto epot = rerunEpot(rerun.log);
    // Frames are written every 50 steps and energies every 10, at steps
    // 0 and 50 of the 100
    ASSERT_EQ(2u, epot.size());
    for (size_t i = 0; i < epot.size(); i++)
    {
        auto row = 5*i;
        ASSERT_LT(row, xvg.rows.size());
        // The trajectory has the coordinates with three decimals in Angstrom,
        // as the starting structure. Later in the simulation that changes
        // the energy of the bonds by up to a few tenths of a kJ/mol.
        EXPECT_NEAR(xvg.rows[row][cEpot], epot[i], i == 0 ? 1e-4 : 0.5) << "frame " << i;
    }
}

// Test: A rerun that a program sets up, without the options of a command
// line, gives the energies of the program
TEST_F(SimulateTest, RerunWithoutCommandLine)
{
    SimulateFiles md;
    auto cmd = baseCommand("water_dimer.pdb", "md-direct", &md);
    addMD(&cmd, 300, 1993);
    ASSERT_EQ(0, run(&cmd));
    SimulateFiles rerun;
    auto cmd2 = baseCommand("water_dimer.pdb", "rerun-tool", &rerun);
    cmd2.addOption("-traj", md.trajectory);
    ASSERT_EQ(0, run(&cmd2));

    const auto log = fileManager().getTemporaryFilePath("rerun-direct.log");
    {
        // A log, which has the energies of the rerun
        MsgHandler            msghandler;
        std::vector<t_pargs>  pa;
        std::vector<t_filenm> fnm;
        msghandler.addOptions(&pa, &fnm, "rerun-direct");
        gmx::test::CommandLine logOnly;
        logOnly.append("test");
        logOnly.addOption("-g", log);
        gmx_output_env_t *oenv = nullptr;
        ASSERT_TRUE(parse_common_args(&logOnly.argc(), logOnly.argv(), 0, fnm.size(), fnm.data(),
                                      pa.size(), pa.data(), 0, nullptr, 0, nullptr, &oenv));
        output_env_done(oenv);
        CommunicationRecord cr;
        cr.init(cr.size());
        msghandler.optionsFinished(fnm, &cr);

        SimulateOptions options;
        options.forceFieldFile = alexandriaFile("forcefields/ACS-g1.xml");
        options.compounds.setFileName(alexandriaFile("molecules/water_dimer.pdb"));
        options.rerun.emplace(false);
        options.rerun->gendimers()->setTrajectory(md.trajectory);
        options.rerun->gendimers()->finishSettings();
        // The options can be copied, with the rerun
        const auto copy = options;
        EXPECT_EQ(0, alexandria::runSimulate(&msghandler, copy, nullptr, &cr));
    }
    const auto tool   = rerunEpot(rerun.log);
    const auto direct = rerunEpot(log);
    ASSERT_EQ(2u, tool.size());
    EXPECT_EQ(tool, direct);
}

// Test: Rerunning a trajectory of the atoms only with a model with shells
// and one with virtual sites gives the energy of a single point calculation
TEST_F(SimulateTest, RerunShellsVsites)
{
    for (const std::string model : { "SWM4-NDP", "TIP5P" })
    {
        const auto ffFile = fileManager().getInputFilePath("forcefields/" + model + ".xml");
        SimulateFiles sp;
        auto cmd = baseCommand("water_dimer.pdb", model + "-sp", &sp, ffFile);
        cmd.append("-sp");
        ASSERT_EQ(0, run(&cmd)) << model;
        SimulateFiles rerun;
        auto cmd2 = baseCommand("water_dimer.pdb", model + "-rerun", &rerun, ffFile);
        cmd2.addOption("-traj", alexandriaFile("molecules/water_dimer.pdb"));
        ASSERT_EQ(0, run(&cmd2)) << model;
        auto epot = rerunEpot(rerun.log);
        ASSERT_EQ(1u, epot.size()) << model;
        // The single point energy is printed with four decimals
        EXPECT_NEAR(singlePointEnergy(sp.log), epot[0], 1e-3) << model;
    }
}

// Test: With -json, part of the output is written in json format to the
// given file, and without it no such file is written
TEST_F(SimulateTest, JsonFile)
{
    const auto ffFile = fileManager().getInputFilePath("forcefields/SWM4-NDP.xml");
    SimulateFiles sp;
    auto cmd  = baseCommand("water_dimer.pdb", "json-sp", &sp, ffFile);
    const auto json = fileManager().getTemporaryFilePath("json-sp.json");
    cmd.append("-sp");
    cmd.addOption("-json", json);
    ASSERT_EQ(0, run(&cmd));
    std::ifstream in(json);
    ASSERT_TRUE(in.good()) << "Cannot open " << json;
    const auto tree = nlohmann::json::parse(in);
    EXPECT_TRUE(tree.contains("simulate")) << tree.dump();

    // Left over from elsewhere, it would hide that it is written by default
    std::remove("simulate.json");
    SimulateFiles noJson;
    auto cmd2 = baseCommand("water_dimer.pdb", "nojson-sp", &noJson, ffFile);
    cmd2.append("-sp");
    ASSERT_EQ(0, run(&cmd2));
    EXPECT_FALSE(std::ifstream("simulate.json").good());
}

// Test: The progress file has one JSON object for each time the energies
// are written, with the same numbers as the energy file, such that a
// program can follow a simulation while it runs
TEST_F(SimulateTest, ProgressFile)
{
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "progress", &files);
    addMD(&cmd, 300, 1993);
    const auto progress = fileManager().getTemporaryFilePath("progress.jsonl");
    cmd.addOption("-progress", progress);
    ASSERT_EQ(0, run(&cmd));

    const auto xvg = readXvg(files.energy);
    ASSERT_EQ(10U, xvg.rows.size());
    std::ifstream               in(progress);
    ASSERT_TRUE(in.good()) << "Cannot open " << progress;
    std::vector<nlohmann::json> lines;
    for (std::string line; std::getline(in, line); )
    {
        const auto js = nlohmann::json::parse(line, nullptr, false);
        ASSERT_TRUE(js.is_object()) << "Not a JSON object: " << line;
        lines.push_back(js);
    }
    // The energies are written every 10 of the 100 steps
    ASSERT_EQ(xvg.rows.size(), lines.size());
    for (size_t i = 0; i < lines.size(); i++)
    {
        const auto &js = lines[i];
        ASSERT_TRUE(js["step"].is_number_integer());
        EXPECT_EQ(static_cast<int>(10*i), js["step"].get<int>());
        EXPECT_EQ(100, js["nsteps"].get<int>());
        ASSERT_TRUE(js["time"]["value"].is_number());
        EXPECT_EQ("ps", js["time"]["unit"].get<std::string>());
        EXPECT_NEAR(xvg.rows[i][0], js["time"]["value"].get<double>(), 1e-9);
        ASSERT_TRUE(js["temperature"]["value"].is_number());
        EXPECT_EQ("K", js["temperature"]["unit"].get<std::string>());
        EXPECT_EQ("kJ/mol", js["energies"]["unit"].get<std::string>());
        // All the columns of the energy file, by their legend
        for (const std::string term : { "EPOT", "EKIN", "ETOT", "BONDS" })
        {
            const int col = xvg.column(term);
            ASSERT_GT(col, 0) << term;
            ASSERT_TRUE(js["energies"].contains(term)) << term;
            ASSERT_TRUE(js["energies"][term].is_number()) << term;
            // The energy file has five decimals
            EXPECT_NEAR(xvg.rows[i][col], js["energies"][term].get<double>(), 1e-4) << term;
        }
        EXPECT_NEAR(xvg.rows[i][xvg.column("TEMPERATURE")],
                    js["temperature"]["value"].get<double>(), 1e-4);
    }
    EXPECT_NEAR(300, lines[0]["temperature"]["value"].get<double>(), 1e-4);
}

// Test: Without the option there is no progress file
TEST_F(SimulateTest, NoProgressFileByDefault)
{
    std::remove("progress.jsonl");
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "noprogress", &files);
    addMD(&cmd, 300, 1993);
    ASSERT_EQ(0, run(&cmd));
    EXPECT_FALSE(std::ifstream("progress.jsonl").good());
}

// Test: The initial temperature is the requested one, computed for
// the degrees of freedom without the translation of and rotation around
// the center of mass and without constraints
TEST_F(SimulateTest, InitialTemperature)
{
    // Acetone has 10 atoms and 9 bonds. Of the 3 angles around the
    // planar carbon two are constrained, and of the 6 angles around
    // each methyl carbon five, leaving the methyl rotations and the
    // out of plane motion of the planar carbon.
    std::map<std::string, double> nDOF = {
        { "no",     3*10-6 },
        { "bonds",  3*10-6-9 },
        { "angles", 3*10-6-9-12 }
    };
    for(const auto &nd : nDOF)
    {
        SimulateFiles files;
        auto cmd = baseCommand("acetone.sdf", "temp" + nd.first, &files);
        addMD(&cmd, 300, 1993);
        cmd.addOption("-constrain", nd.first);
        ASSERT_EQ(0, run(&cmd));
        auto xvg   = readXvg(files.energy);
        ASSERT_FALSE(xvg.rows.empty());
        int  cEkin = xvg.column("EKIN");
        int  cTemp = xvg.column("TEMPERATURE");
        ASSERT_GT(cEkin, 0);
        ASSERT_GT(cTemp, 0);
        auto ekin  = xvg.rows[0][cEkin];
        auto temp  = xvg.rows[0][cTemp];
        EXPECT_NEAR(300, temp, 1e-4) << nd.first;
        EXPECT_NEAR(nd.second, 2*ekin/(BOLTZ*temp), 1e-3) << nd.first;
    }
}

// Test: The number of degrees of freedom depends on the motion of the
// center of mass that is removed, with two rotations for a linear molecule
TEST_F(SimulateTest, CenterOfMassDegreesOfFreedom)
{
    // Number of atoms and whether the molecule is linear
    std::map<std::string, std::pair<int, bool>> mols = {
        { "acetone.sdf",           { 10, false } },
        { "carbon-dioxide.sdf",    { 3, true } },
        { "hydrogen-chloride.sdf", { 2, true } }
    };
    for (const auto &mol : mols)
    {
        int natom = mol.second.first;
        std::map<std::string, int> nDOF = {
            { "none",    3*natom },
            { "linear",  3*natom-3 },
            { "angular", 3*natom-3-(mol.second.second ? 2 : 3) }
        };
        for (const auto &nd : nDOF)
        {
            SCOPED_TRACE(mol.first + " " + nd.first);
            SimulateFiles files;
            auto cmd = baseCommand(mol.first, "comm" + nd.first + mol.first, &files);
            addMD(&cmd, 300, 1993);
            cmd.addOption("-comm", nd.first);
            cmd.addOption("-nstcomm", 10);
            ASSERT_EQ(0, run(&cmd));
            auto xvg   = readXvg(files.energy);
            ASSERT_FALSE(xvg.rows.empty());
            int  cEkin = xvg.column("EKIN");
            int  cTemp = xvg.column("TEMPERATURE");
            EXPECT_NEAR(300, xvg.rows[0][cTemp], 1e-4);
            EXPECT_NEAR(nd.second, 2*xvg.rows[0][cEkin]/(BOLTZ*xvg.rows[0][cTemp]), 1e-3);
            checkConservation(xvg);
        }
    }
}

// Test: Removing the motion of the center of mass often in a cluster
// does not affect the energy conservation, since the momenta are
// conserved anyway
TEST_F(SimulateTest, CenterOfMassRemovalWaterDimer)
{
    std::vector<XvgData> xvg;
    for (int nstcomm : { 0, 1 })
    {
        SimulateFiles files;
        auto cmd = baseCommand("water_dimer.pdb", gmx::formatString("nstcomm%d", nstcomm), &files);
        addMD(&cmd, 300, 1993);
        cmd.addOption("-nstcomm", nstcomm);
        ASSERT_EQ(0, run(&cmd));
        xvg.push_back(readXvg(files.energy));
        checkConservation(xvg.back());
    }
    // The angular momentum is only changed by numerical errors, so the
    // trajectories are nearly the same
    int cEpot = xvg[0].column("EPOT");
    ASSERT_GT(cEpot, 0);
    ASSERT_EQ(xvg[0].rows.size(), xvg[1].rows.size());
    for (size_t i = 0; i < xvg[0].rows.size(); i++)
    {
        EXPECT_NEAR(xvg[0].rows[i][cEpot], xvg[1].rows[i][cEpot], 1e-3) << i;
    }
}

// Test: The total energy is conserved with constraints
TEST_F(SimulateTest, ConstrainedEnergyConservation)
{
    for(const std::string constrain : { "bonds", "angles" })
    {
        SimulateFiles files;
        auto cmd = baseCommand("acetone.sdf", "conserve" + constrain, &files);
        cmd.addOption("-nsteps", 2500);
        cmd.addOption("-deltat", 0.0002);
        cmd.addOption("-temp", 300);
        cmd.addOption("-seed", 1993);
        cmd.addOption("-nstener", 10);
        cmd.addOption("-constrain", constrain);
        ASSERT_EQ(0, run(&cmd));
        auto xvg   = readXvg(files.energy);
        int  cEtot = xvg.column("ETOT");
        ASSERT_GT(cEtot, 0);
        size_t nav = 10;
        ASSERT_GT(xvg.rows.size(), 2*nav);
        double first = 0, last = 0;
        for(size_t i = 0; i < nav; i++)
        {
            first += xvg.rows[i][cEtot];
            last  += xvg.rows[xvg.rows.size()-1-i][cEtot];
        }
        EXPECT_NEAR(first/nav, last/nav, 0.01) << constrain;
    }
}

// Test: The constrained distances stay the same along an MD trajectory,
// the bond lengths with bond constraints and also the distances between
// the outer atoms of the angles with angle constraints. At a planar atom
// with three neighbours two angles are constrained, the third can change
// when the atom moves out of the plane. Without constraints the bonds
// vibrate, which shows that the test can tell.
TEST_F(SimulateTest, ConstraintsKeepDistances)
{
    auto distance = [](const std::vector<double> &a, const std::vector<double> &b)
    {
        double d2 = 0;
        for(int m = 0; m < DIM; m++)
        {
            d2 += (a[m] - b[m])*(a[m] - b[m]);
        }
        return std::sqrt(d2);
    };
    for(const std::string constrain : { "no", "bonds", "angles" })
    {
        SimulateFiles files;
        auto cmd = baseCommand("acetone.sdf", "keep" + constrain, &files);
        cmd.addOption("-nsteps", 1000);
        cmd.addOption("-deltat", 0.0002);
        cmd.addOption("-temp", 300);
        cmd.addOption("-seed", 1993);
        cmd.addOption("-nstxout", 50);
        cmd.addOption("-constrain", constrain);
        ASSERT_EQ(0, run(&cmd)) << constrain;
        auto frames = readPdbFrames(files.trajectory);
        ASSERT_GT(frames.size(), 10u) << constrain;
        const auto &x0    = frames[0];
        size_t      natom = x0.size();
        // The bonds from the first frame, all bonds are shorter than 1.65 A
        // and the other distances longer
        std::vector<std::pair<size_t, size_t>> bonds;
        std::vector<std::vector<size_t>>       neighbours(natom);
        for(size_t i = 0; i < natom; i++)
        {
            for(size_t j = i+1; j < natom; j++)
            {
                if (distance(x0[i], x0[j]) < 1.65)
                {
                    bonds.emplace_back(i, j);
                    neighbours[i].push_back(j);
                    neighbours[j].push_back(i);
                }
            }
        }
        ASSERT_EQ(9u, bonds.size()) << "acetone has nine bonds";
        // The largest change of a distance along the trajectory
        const auto maxChange = [&](const std::vector<std::pair<size_t, size_t>> &pairs)
        {
            double dmax = 0;
            for(const auto &x : frames)
            {
                EXPECT_EQ(natom, x.size());
                for(const auto &p : pairs)
                {
                    dmax = std::max(dmax, std::abs(distance(x[p.first], x[p.second]) -
                                                   distance(x0[p.first], x0[p.second])));
                }
            }
            return dmax;
        };
        // The coordinates are written with three decimals
        double toler = 2e-3;
        if (constrain == "no")
        {
            EXPECT_GT(maxChange(bonds), 0.01);
        }
        else
        {
            EXPECT_LT(maxChange(bonds), toler) << constrain;
            // The number of angles of each atom that stay the same
            size_t nangle = 0, nfixed = 0;
            for(size_t c = 0; c < natom; c++)
            {
                const auto &nb     = neighbours[c];
                size_t      nc     = nb.size()*(nb.size() - 1)/2;
                size_t      nfixc  = 0;
                for(size_t k = 0; k < nb.size(); k++)
                {
                    for(size_t l = k+1; l < nb.size(); l++)
                    {
                        nfixc += maxChange({ { nb[k], nb[l] } }) < toler ? 1 : 0;
                    }
                }
                if (constrain == "angles")
                {
                    EXPECT_GE(nfixc, nb.size() == 3 ? 2 : nc) << "atom " << c;
                }
                nangle += nc;
                nfixed += nfixc;
            }
            if (constrain == "bonds")
            {
                // The angles are free
                EXPECT_LT(nfixed, nangle/2);
            }
        }
    }
}

// Test: Angle constraints work in a ring
TEST_F(SimulateTest, AngleConstraintsUracil)
{
    SimulateFiles files;
    auto cmd = baseCommand("uracil.sdf", "uracilangles", &files);
    cmd.addOption("-nsteps", 2500);
    cmd.addOption("-deltat", 0.0002);
    cmd.addOption("-temp", 300);
    cmd.addOption("-seed", 1993);
    cmd.addOption("-nstener", 10);
    cmd.addOption("-constrain", "angles");
    ASSERT_EQ(0, run(&cmd));
    auto xvg   = readXvg(files.energy);
    ASSERT_FALSE(xvg.rows.empty());
    // Uracil has 12 atoms and 12 bonds. Of the 18 angles, two around
    // each of the six planar ring atoms are constrained, but three of
    // those are redundant since they close the ring.
    int  cEkin = xvg.column("EKIN");
    int  cTemp = xvg.column("TEMPERATURE");
    int  cEtot = xvg.column("ETOT");
    ASSERT_GT(cEtot, 0);
    EXPECT_NEAR(3*12-6-12-9, 2*xvg.rows[0][cEkin]/(BOLTZ*xvg.rows[0][cTemp]), 1e-3);
    EXPECT_NEAR(xvg.rows[0][cEtot], xvg.rows.back()[cEtot], 0.01);
}

// Test: Angle constraints keep the angles at their reference values
TEST_F(SimulateTest, AngleConstraintsWater)
{
    SimulateFiles files;
    auto cmd = baseCommand("water-3-oep.log.pdb", "waterangles", &files);
    addMD(&cmd, 300, 1993);
    cmd.addOption("-constrain", "angles");
    // Keep the rotation, which is all there is to a rigid molecule
    cmd.addOption("-comm", "linear");
    ASSERT_EQ(0, run(&cmd));
    auto xvg = readXvg(files.energy);
    ASSERT_FALSE(xvg.rows.empty());
    // Water is rigid, so only the rotations remain
    int  cEkin = xvg.column("EKIN");
    int  cTemp = xvg.column("TEMPERATURE");
    EXPECT_NEAR(3, 2*xvg.rows[0][cEkin]/(BOLTZ*xvg.rows[0][cTemp]), 1e-3);
    // With a rigid molecule the potential energy does not change
    int  cEpot = xvg.column("EPOT");
    ASSERT_GT(cEpot, 0);
    for(const auto &row : xvg.rows)
    {
        EXPECT_NEAR(xvg.rows[0][cEpot], row[cEpot], 1e-3);
    }
}

// Test: The minimizers find the same structure with constraints
TEST_F(SimulateTest, MinimizeConstrained)
{
    const std::string molFile("methanethiol.sdf");
    const auto distance = [](const std::vector<std::vector<double>> &x, int i, int j)
    {
        double d2 = 0;
        for(int m = 0; m < DIM; m++)
        {
            double dx = x[i][m] - x[j][m];
            d2       += dx*dx;
        }
        return std::sqrt(d2);
    };
    std::map<std::string, std::vector<std::vector<double>>> xmin;
    for(const std::string constrain : { "bonds", "angles" })
    {
        for(const std::string alg : { "LBFGS", "Steep", "Newton" })
        {
            SimulateFiles files;
            const auto label   = "min" + constrain + alg;
            auto cmd     = baseCommand(molFile, label, &files);
            const auto confout = fileManager().getTemporaryFilePath(label + "-confout.pdb");
            cmd.addOption("-minimize");
            cmd.addOption("-minalg", alg);
            cmd.addOption("-maxiter", 1000);
            cmd.addOption("-constrain", constrain);
            cmd.addOption("-c", confout);
            ASSERT_EQ(0, run(&cmd)) << label;
            xmin.insert({ label, readPdbCoordinates(confout) });
        }
        // The structures can be translated or rotated with respect
        // to each other, so compare all distances. Coordinates are
        // written with three decimals.
        double toler = 3e-3;
        auto  &xref  = xmin["min" + constrain + "LBFGS"];
        for(const std::string alg : { "Steep", "Newton" })
        {
            auto &x = xmin["min" + constrain + alg];
            ASSERT_EQ(xref.size(), x.size());
            for(size_t i = 0; i < x.size(); i++)
            {
                for(size_t j = i+1; j < x.size(); j++)
                {
                    EXPECT_NEAR(distance(xref, i, j), distance(x, i, j), toler)
                        << constrain << " " << alg << " atoms " << i+1 << " " << j+1;
                }
            }
        }
    }
    // Both constrain the bonds to the same lengths. The atoms are S H C H H H.
    auto &xb = xmin["minbondsLBFGS"];
    auto &xa = xmin["minanglesLBFGS"];
    for(const auto &bond : std::vector<std::pair<int, int>>{ { 0, 1 }, { 0, 2 }, { 2, 3 }, { 2, 4 }, { 2, 5 } })
    {
        EXPECT_NEAR(distance(xb, bond.first, bond.second),
                    distance(xa, bond.first, bond.second), 3e-3) << "bond " << bond.first+1 << "-" << bond.second+1;
    }
    // The angle constraints change the distances between the hydrogens
    double maxDiff = 0;
    for(size_t i = 0; i < xa.size(); i++)
    {
        for(size_t j = i+1; j < xa.size(); j++)
        {
            maxDiff = std::max(maxDiff, std::abs(distance(xa, i, j) - distance(xb, i, j)));
        }
    }
    EXPECT_LT(0.03, maxDiff);
}

// Test: Freezing atoms in a minimization with constraints
TEST_F(SimulateTest, MinimizeFreezeConstrained)
{
    const std::string molFile("acetone-3-oep.log.pdb");
    SimulateFiles files;
    auto cmd     = baseCommand(molFile, "freezecons", &files);
    const auto confout = fileManager().getTemporaryFilePath("freezecons-confout.pdb");
    const auto ndx     = fileManager().getTemporaryFilePath("freezecons.ndx");
    {
        std::ofstream out(ndx);
        out << "[ Frozen ]\n1 2\n";
    }
    cmd.addOption("-minimize");
    cmd.addOption("-c", confout);
    cmd.addOption("-freeze", ndx);
    cmd.addOption("-constrain", "angles");
    ASSERT_EQ(0, run(&cmd));
    auto xin  = readPdbCoordinates(alexandriaFile("molecules/" + molFile));
    auto xout = readPdbCoordinates(confout);
    ASSERT_EQ(xin.size(), xout.size());
    for(size_t i = 0; i < 2; i++)
    {
        for(int m = 0; m < 3; m++)
        {
            EXPECT_NEAR(xin[i][m], xout[i][m], 1.5e-3) << "atom " << i+1;
        }
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
    const auto confout = fileManager().getTemporaryFilePath("minmd-confout.pdb");
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
    const auto confout = fileManager().getTemporaryFilePath("mindimer-confout.pdb");
    cmd.addOption("-minimize");
    cmd.addOption("-c", confout);
    EXPECT_EQ(0, run(&cmd));
    EXPECT_TRUE(fileNotEmpty(confout));
}

// Test: Frozen atoms do not move during minimization
TEST_F(SimulateTest, MinimizeFreeze)
{
    const std::string molFile("acetone-3-oep.log.pdb");
    SimulateFiles files;
    auto cmd     = baseCommand(molFile, "freeze", &files);
    const auto confout = fileManager().getTemporaryFilePath("freeze-confout.pdb");
    const auto ndx     = fileManager().getTemporaryFilePath("freeze.ndx");
    {
        std::ofstream out(ndx);
        out << "[ Frozen ]\n1 2\n";
    }
    cmd.addOption("-minimize");
    cmd.addOption("-c", confout);
    cmd.addOption("-freeze", ndx);
    ASSERT_EQ(0, run(&cmd));
    auto xin  = readPdbCoordinates(alexandriaFile("molecules/" + molFile));
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

// Test: MD with flat-bottomed position restraints around the center of
// mass of a water dimer far from the origin. Within a sphere that holds the
// dimer they have no energy, the motion of the center of mass is not
// removed, and with a radius of zero they act.
TEST_F(SimulateTest, PositionRestraints)
{
    // The water dimer moved from the origin
    const auto molFile = fileManager().getTemporaryFilePath("shifted.pdb");
    {
        std::ifstream in(alexandriaFile("molecules/water_dimer.pdb"));
        std::ofstream out(molFile);
        for (std::string line; std::getline(in, line); )
        {
            if (line.compare(0, 6, "HETATM") == 0 || line.compare(0, 4, "ATOM") == 0)
            {
                double x[3];
                for (int m = 0; m < 3; m++)
                {
                    x[m] = std::stod(line.substr(30 + 8*m, 8)) + 50 + 10*m;
                }
                line = line.substr(0, 30) + gmx::formatString("%8.3f%8.3f%8.3f", x[0], x[1], x[2]) + line.substr(54);
            }
            out << line << "\n";
        }
    }
    for (double r0 : { 1.0, 0.0 })
    {
        SimulateFiles files;
        auto cmd = baseCommand(molFile, gmx::formatString("fbpr%g", r0), &files);
        addMD(&cmd, 300, 1993);
        cmd.addOption("-kfbpr", 1000.0);
        cmd.addOption("-r0fbpr", r0);
        ASSERT_EQ(0, run(&cmd)) << r0;
        auto xvg = readXvg(files.energy);
        int  cpr = xvg.column("POSITION_RESTRAINT");
        ASSERT_GT(cpr, 0);
        ASSERT_FALSE(xvg.rows.empty());
        double maxPR = 0;
        for (const auto &row : xvg.rows)
        {
            maxPR = std::max(maxPR, row[cpr]);
        }
        if (r0 > 0)
        {
            EXPECT_EQ(0, maxPR);
        }
        else
        {
            EXPECT_GT(xvg.rows[0][cpr], 0);
        }
        const std::ifstream in(files.log);
        std::stringstream log;
        log << in.rdbuf();
        EXPECT_NE(std::string::npos, log.str().find("The motion of the center of mass is not removed")) << r0;
    }
}

// Test: Writing the molprop file
TEST_F(SimulateTest, WriteMolprop)
{
    SimulateFiles files;
    auto cmd = baseCommand("acetone.sdf", "omp", &files);
    const auto omp = fileManager().getTemporaryFilePath("omp.xml");
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

//! \brief Make a command line for the minimize tool
static gmx::test::CommandLine minimizeCommand(const std::string &ff,
                                              const std::string &molFile,
                                              const std::string &confout,
                                              const std::string &log)
{
    gmx::test::CommandLine cmd;
    cmd.append("minimize");
    cmd.addOption("-ff", ff);
    cmd.addOption("-f", molFile);
    cmd.addOption("-c", confout);
    cmd.addOption("-g", log);
    return cmd;
}

// Test: alexandria minimize gives the same structure as simulate -minimize
TEST_F(SimulateTest, MinimizeTool)
{
    for(const std::string constrain : { "no", "bonds" })
    {
        const std::string molFile("water_dimer.pdb");
        SimulateFiles files;
        auto cmd     = baseCommand(molFile, "simmin" + constrain, &files);
        const auto simout  = fileManager().getTemporaryFilePath("simmin" + constrain + "-confout.pdb");
        cmd.addOption("-minimize");
        cmd.addOption("-constrain", constrain);
        cmd.addOption("-c", simout);
        ASSERT_EQ(0, run(&cmd));

        const auto minout  = fileManager().getTemporaryFilePath("minimize" + constrain + "-confout.pdb");
        auto mincmd  = minimizeCommand(alexandriaFile("forcefields/ACS-g1.xml"),
                                       alexandriaFile("molecules/" + molFile), minout,
                                       fileManager().getTemporaryFilePath("minimize" + constrain + ".log"));
        mincmd.addOption("-constrain", constrain);
        ASSERT_EQ(0, minimize(mincmd.argc(), mincmd.argv()));

        auto xsim = readPdbCoordinates(simout);
        auto xmin = readPdbCoordinates(minout);
        ASSERT_EQ(xsim.size(), xmin.size());
        for(size_t i = 0; i < xsim.size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                EXPECT_EQ(xsim[i][m], xmin[i][m]) << constrain << " atom " << i+1;
            }
        }
    }
}

// Test: alexandria minimize that does not converge within the maximum
// number of steps writes the structure with the lowest energy, and the
// energies before and after, which are lower
TEST_F(SimulateTest, MinimizeToolNotConverged)
{
    for(const std::string alg : { "LBFGS", "Steep", "Newton" })
    {
        const auto confout = fileManager().getTemporaryFilePath("notconverged" + alg + "-confout.pdb");
        const auto log     = fileManager().getTemporaryFilePath("notconverged" + alg + ".log");
        auto cmd     = minimizeCommand(alexandriaFile("forcefields/ACS-g1.xml"),
                                       alexandriaFile("molecules/acetone.sdf"), confout, log);
        cmd.addOption("-minalg", alg);
        cmd.addOption("-maxiter", 2);
        EXPECT_EQ(1, minimize(cmd.argc(), cmd.argv())) << alg;
        EXPECT_EQ(10u, readPdbCoordinates(confout).size()) << alg;
        // The potential energy before and after minimization
        std::vector<double> epot;
        std::ifstream       in(log);
        std::string         text;
        for (std::string line; std::getline(in, line); )
        {
            text += line + "\n";
            std::istringstream iss(line);
            std::string        word, valueWord;
            double             value;
            if (iss >> word >> valueWord >> value && word == "EPOT:")
            {
                epot.push_back(value);
            }
        }
        EXPECT_NE(std::string::npos, text.find("Minimization did not converge: Too many steps")) << alg;
        ASSERT_EQ(2u, epot.size()) << alg;
        EXPECT_LT(epot[1], epot[0]) << alg;
    }
}

// Test: alexandria minimize relaxes the shells, the induction energy is
// negative before and after minimization, and so is that of the interaction
// energy of the dimer after minimization
TEST_F(SimulateTest, MinimizeToolShells)
{
    const auto confout = fileManager().getTemporaryFilePath("shells-confout.pdb");
    const auto log     = fileManager().getTemporaryFilePath("minimize-shells.log");
    auto cmd     = minimizeCommand(fileManager().getInputFilePath("forcefields/SWM4-NDP.xml"),
                                   alexandriaFile("molecules/water_dimer.pdb"), confout, log);
    ASSERT_EQ(0, minimize(cmd.argc(), cmd.argv()));
    std::vector<double> induction;
    std::ifstream       in(log);
    for (std::string line; std::getline(in, line); )
    {
        std::istringstream iss(line);
        std::string        word, valueWord;
        double             value;
        if (iss >> word >> valueWord >> value && word == "INDUCTION:")
        {
            induction.push_back(value);
        }
    }
    ASSERT_EQ(3u, induction.size());
    for (auto ind : induction)
    {
        EXPECT_LT(ind, 0);
    }
}

// Test: alexandria minimize fails for a missing molecule
TEST_F(SimulateTest, MinimizeToolMissingMolecule)
{
    auto cmd = minimizeCommand(alexandriaFile("forcefields/ACS-g1.xml"),
                               alexandriaFile("molecules/does-not-exist.sdf"),
                               fileManager().getTemporaryFilePath("missing-confout.pdb"),
                               fileManager().getTemporaryFilePath("minimize-missing.log"));
    EXPECT_EQ(1, minimize(cmd.argc(), cmd.argv()));
}

// Simulations at constant temperature with beads, -nbeads. The random
// numbers of the thermostat differ between standard libraries, so the
// results are compared with independent calculations within their
// statistical errors, and not with reference data.

class SimulateBeadsTest : public SimulateTest
{
protected:
    /*! \brief Make a command line for a pair of noble gas atoms in a wall
     * The wall has a radius of 0.8 nm and a force constant of 1000 kJ/mol nm^2.
     * \param[in]  molFile The pair in src/act/testdata/molecules
     * \param[in]  label   Unique label for the output files
     * \param[out] files   The output file names
     * \param[in]  nbeads  The number of beads
     * \param[in]  T       The temperature
     * \param[in]  nsteps  The number of steps
     * \return The command line
     */
    gmx::test::CommandLine pairCommand(const std::string &molFile,
                                       const std::string &label,
                                       SimulateFiles     *files,
                                       int                nbeads,
                                       double             T,
                                       int                nsteps)
    {
        auto cmd = baseCommand(molFile, label, files,
                               alexandriaFile("forcefields/noble-halgren14-7.xml"));
        cmd.addOption("-nbeads", nbeads);
        cmd.addOption("-temp", T);
        cmd.addOption("-nsteps", nsteps);
        cmd.addOption("-deltat", 0.005);
        cmd.addOption("-tau", 1.0);
        cmd.addOption("-rwall", 0.8);
        cmd.addOption("-kwall", 1000.0);
        cmd.addOption("-nstener", 0);
        cmd.addOption("-nstxout", 0);
        cmd.addOption("-seed", 1993);
        return cmd;
    }

    /*! \brief Run simulate and read the results of the simulation with beads
     * \param[in] cmd   The command line, the JSON file is added
     * \param[in] label Unique label for the JSON file
     * \return what the JSON file has for the simulation with beads
     */
    nlohmann::json runBeads(gmx::test::CommandLine *cmd, const std::string &label)
    {
        const auto json = fileManager().getTemporaryFilePath(label + ".json");
        cmd->addOption("-json", json);
        EXPECT_EQ(0, run(cmd));
        std::ifstream in(json);
        if (!in.good())
        {
            ADD_FAILURE() << "Cannot open " << json;
            return {};
        }
        const auto        tree = nlohmann::json::parse(in);
        const std::string key("Simulation with beads");
        if (!tree.contains("simulate") || !tree["simulate"].contains(key))
        {
            ADD_FAILURE() << "No simulation with beads in " << json;
            return {};
        }
        return tree["simulate"][key];
    }

    //! \return the average of a quantity in the results of a simulation with beads
    static double average(const nlohmann::json &results, const std::string &quantity)
    {
        if (!results.contains(quantity) || !results[quantity].contains("Average"))
        {
            ADD_FAILURE() << "No " << quantity << " in " << results.dump();
            return 0;
        }
        return results[quantity]["Average"].get<double>();
    }

    //! \return the standard error of a quantity in the results of a simulation with beads
    static double error(const nlohmann::json &results, const std::string &quantity)
    {
        if (!results.contains(quantity) || !results[quantity].contains("Standard error"))
        {
            ADD_FAILURE() << "No error of " << quantity << " in " << results.dump();
            return 0;
        }
        return results[quantity]["Standard error"].get<double>();
    }

    /*! \brief Check that simulate refuses a command line
     * \param[in] cmd  The command line
     * \param[in] what What is wrong with it
     */
    void expectRefused(gmx::test::CommandLine *cmd, const std::string &what)
    {
        EXPECT_EQ(1, run(cmd)) << what;
    }
};

TEST_F(SimulateBeadsTest, ArgonPairInWall)
{
    // Two argon atoms with the buffered 14-7 potential at 80 K, kept
    // together by a wall on their distance R. Classically the distribution
    // of R is R^2 times the Boltzmann factor of the potential and the wall,
    // which gives by quadrature <R> = 0.5610 nm, <V> = -0.6962 kJ/mol,
    // 45.98 % of R below 0.5 nm and 6.34 % beyond the wall.
    const double kT = BOLTZ*80;
    {
        SimulateFiles files;
        auto cmd = pairCommand("argon#argon.pdb", "argon-1", &files, 1, 80, 1000000);
        const auto hist = fileManager().getTemporaryFilePath("argon-1-hist.xvg");
        cmd.addOption("-hist", hist);
        const auto res = runBeads(&cmd, "argon-1");
        ASSERT_TRUE(res.contains("all compounds")) << res.dump();
        const auto &all = res["all compounds"];
        EXPECT_EQ(1, all["Beads"].get<int>());
        EXPECT_EQ(900000, all["Steps"].get<int>());
        // The standard error of the distance is about 0.005 nm
        EXPECT_NEAR(0.5610, average(all, "Distance"), 0.03);
        EXPECT_NEAR(-0.6962, average(all, "Potential energy"), 0.08);
        EXPECT_NEAR(0.0634, all["Fraction at the wall"].get<double>(), 0.025);
        // One bead has the classical kinetic energy, 3 kT for two atoms,
        // and no spread
        EXPECT_NEAR(3*kT, average(all, "Kinetic energy"), 1e-9);
        EXPECT_EQ(0.0, average(all, "Spread of the beads"));
        EXPECT_GT(error(all, "Distance"), 0.0);
        EXPECT_LT(error(all, "Distance"), 0.02);
        // The histogram of the distance is a probability density, and the
        // potential of mean force is the potential plus a constant:
        // V(0.501) - V(0.415) = 0.8615 kJ/mol
        const auto xvg = readXvg(hist);
        ASSERT_EQ(2U, xvg.legends.size());
        ASSERT_GT(xvg.rows.size(), 100U);
        double total = 0, below = 0, wmin = 1e9, w415 = 0, w501 = 0;
        for(const auto &row : xvg.rows)
        {
            ASSERT_EQ(3U, row.size());
            total += row[1]*0.002;
            below += (row[0] < 0.5) ? row[1]*0.002 : 0.0;
            wmin   = std::min(wmin, row[2]);
            if (std::abs(row[0] - 0.415) < 1e-6)
            {
                w415 = row[2];
            }
            if (std::abs(row[0] - 0.501) < 1e-6)
            {
                w501 = row[2];
            }
        }
        EXPECT_NEAR(1.0, total, 1e-4);
        EXPECT_NEAR(0.4598, below, 0.06);
        EXPECT_EQ(0.0, wmin);
        EXPECT_NEAR(0.8615, w501 - w415, 0.25);
        // No energy file and no trajectory were asked for
        EXPECT_FALSE(fileNotEmpty(files.trajectory));
    }
    {
        // With 16 beads an independent path integral, in which the ring
        // polymer of the distance vector is drawn from its normal
        // distribution, gives <R> = 0.5618 nm, <V> = -0.68995 kJ/mol and a
        // kinetic energy of 1.99875 kJ/mol, 0.0033 more than 3 kT
        SimulateFiles files;
        auto cmd = pairCommand("argon#argon.pdb", "argon-16", &files, 16, 80, 300000);
        const auto res = runBeads(&cmd, "argon-16");
        ASSERT_TRUE(res.contains("all compounds")) << res.dump();
        const auto &all = res["all compounds"];
        EXPECT_EQ(16, all["Beads"].get<int>());
        // The standard error of the distance is about 0.01 nm
        EXPECT_NEAR(0.5618, average(all, "Distance"), 0.05);
        EXPECT_NEAR(-0.68995, average(all, "Potential energy"), 0.2);
        EXPECT_NEAR(1.99875, average(all, "Kinetic energy"), 0.0012);
        EXPECT_GT(average(all, "Kinetic energy"), 3*kT + 0.002);
        EXPECT_NEAR(average(all, "Potential energy") + average(all, "Kinetic energy"),
                    average(all, "Total energy"), 1e-9);
        // The beads of argon at 80 K are spread over 6 pm
        EXPECT_NEAR(0.0061, average(all, "Spread of the beads"), 0.0005);
    }
}

TEST_F(SimulateBeadsTest, HeliumPairQuantum)
{
    // Two helium atoms at 20 K in the same wall. Classically
    // <R> = 0.5929 nm and <V> = -0.07122 kJ/mol. With 16 beads the
    // independent path integral gives <R> = 0.6080 nm, <V> = -0.05651
    // kJ/mol and a kinetic energy of 0.50714 kJ/mol, where 3 kT is 0.49887:
    // the atoms avoid each other.
    const double kT = BOLTZ*20;
    double       vclass = 0;
    {
        SimulateFiles files;
        auto cmd = pairCommand("helium#helium.pdb", "helium-1", &files, 1, 20, 600000);
        const auto res = runBeads(&cmd, "helium-1");
        ASSERT_TRUE(res.contains("all compounds")) << res.dump();
        const auto &all = res["all compounds"];
        EXPECT_NEAR(0.5929, average(all, "Distance"), 0.03);
        EXPECT_NEAR(-0.07122, average(all, "Potential energy"), 0.015);
        EXPECT_NEAR(3*kT, average(all, "Kinetic energy"), 1e-9);
        vclass = average(all, "Potential energy");
    }
    {
        SimulateFiles files;
        auto cmd = pairCommand("helium#helium.pdb", "helium-16", &files, 16, 20, 400000);
        const auto res = runBeads(&cmd, "helium-16");
        ASSERT_TRUE(res.contains("all compounds")) << res.dump();
        const auto &all = res["all compounds"];
        EXPECT_NEAR(0.6080, average(all, "Distance"), 0.03);
        // The standard errors are about 0.0015 and 0.0004 kJ/mol
        EXPECT_NEAR(-0.05651, average(all, "Potential energy"), 0.008);
        EXPECT_NEAR(0.50714, average(all, "Kinetic energy"), 0.002);
        EXPECT_GT(average(all, "Potential energy"), vclass);
    }
}

TEST_F(SimulateBeadsTest, OneBeadIsClassical)
{
    // Water with a thermostat: the kinetic energy that is written is the
    // classical 9/2 kT in every step, and the momenta have the temperature
    // for the 6 degrees of freedom that are left without the translation
    SimulateFiles files;
    auto cmd = baseCommand("water-3-oep.log.pdb", "one-bead", &files);
    cmd.addOption("-minimize");
    cmd.addOption("-nbeads", 1);
    cmd.addOption("-temp", 300.0);
    cmd.addOption("-nsteps", 60000);
    cmd.addOption("-deltat", 0.00025);
    cmd.addOption("-tau", 0.01);
    cmd.addOption("-equil", 2000);
    cmd.addOption("-nstener", 20);
    cmd.addOption("-nstxout", 0);
    cmd.addOption("-seed", 1993);
    const auto res = runBeads(&cmd, "one-bead");
    ASSERT_TRUE(res.contains("the compound")) << res.dump();
    const auto &one = res["the compound"];
    EXPECT_EQ(58000, one["Steps"].get<int>());
    EXPECT_FALSE(one.contains("Distance"));
    const auto xvg = readXvg(files.energy);
    ASSERT_EQ(3000U, xvg.rows.size());
    const int cEkin = xvg.column("EKIN");
    const int cEtot = xvg.column("ETOT");
    const int cEpot = xvg.column("EPOT");
    const int cTemp = xvg.column("TEMPERATURE");
    const int cRgyr = xvg.column("RGYR");
    ASSERT_GT(cEkin, 0);
    ASSERT_GT(cEtot, 0);
    ASSERT_GT(cEpot, 0);
    ASSERT_GT(cTemp, 0);
    ASSERT_GT(cRgyr, 0);
    // One compound: no wall and no distance
    EXPECT_EQ(-1, xvg.column("WALL"));
    EXPECT_EQ(-1, xvg.column("DISTANCE"));
    double temp = 0, epot = 0;
    for(const auto &row : xvg.rows)
    {
        EXPECT_NEAR(4.5*BOLTZ*300, row[cEkin], 2e-5);
        EXPECT_NEAR(row[cEpot] + row[cEkin], row[cEtot], 2e-5);
        EXPECT_EQ(0.0, row[cRgyr]);
        temp += row[cTemp]/xvg.rows.size();
        epot += row[cEpot]/xvg.rows.size();
    }
    // The temperature of one configuration fluctuates by 170 K, the
    // standard error of the average is about 6 K
    EXPECT_NEAR(300, temp, 30);
    // Three vibrations with kT/2 of potential energy each, above the minimum
    EXPECT_NEAR(1.5*BOLTZ*300, epot - xvg.rows[0][cEpot], 0.3*1.5*BOLTZ*300);
    // The energy terms are averaged over the energies written after the equilibration
    ASSERT_TRUE(one.contains("Energy terms") && one["Energy terms"].contains("EPOT")) << one.dump();
    EXPECT_NEAR(average(one, "Potential energy"), one["Energy terms"]["EPOT"]["Average"].get<double>(),
                0.3*1.5*BOLTZ*300);
}

TEST_F(SimulateBeadsTest, OutputFilesAndBinding)
{
    // A water dimer with shells and four beads in a wall that it does not reach
    SimulateFiles files;
    auto cmd = baseCommand("water_dimer.pdb", "beads-files", &files,
                           alexandriaFile("forcefields/SWM4-NDP.xml"));
    const auto progress = fileManager().getTemporaryFilePath("beads-progress.jsonl");
    const auto hist     = fileManager().getTemporaryFilePath("beads-hist.xvg");
    cmd.addOption("-minimize");
    cmd.addOption("-nbeads", 4);
    cmd.addOption("-temp", 100.0);
    cmd.addOption("-nsteps", 400);
    cmd.addOption("-deltat", 0.00025);
    cmd.addOption("-nstener", 10);
    cmd.addOption("-nstxout", 100);
    cmd.addOption("-seed", 1993);
    cmd.addOption("-rwall", 0.6);
    cmd.addOption("-monomers");
    cmd.addOption("-progress", progress);
    cmd.addOption("-hist", hist);
    const auto res = runBeads(&cmd, "beads-files");
    ASSERT_TRUE(res.contains("all compounds") && res.contains("compound 1") &&
                res.contains("compound 2") && res.contains("Binding")) << res.dump();
    const auto &all = res["all compounds"];
    // A tenth of the steps is equilibration
    EXPECT_EQ(360, all["Steps"].get<int>());
    EXPECT_EQ(0.0, all["Fraction at the wall"].get<double>());
    EXPECT_FALSE(all.contains("Moved apart"));
    EXPECT_NEAR(0.29, average(all, "Distance"), 0.03);
    EXPECT_FALSE(res["compound 1"].contains("Distance"));
    EXPECT_FALSE(res["compound 1"].contains("Fraction at the wall"));
    // The binding energy from the three simulations
    double deltaU = (average(all, "Total energy") - average(res["compound 1"], "Total energy") -
                     average(res["compound 2"], "Total energy"));
    EXPECT_NEAR(deltaU, res["Binding"]["Energy"].get<double>(), 1e-9);
    EXPECT_NEAR(deltaU - BOLTZ*100, res["Binding"]["Enthalpy"].get<double>(), 1e-6);
    EXPECT_GT(res["Binding"]["Standard error"].get<double>(), 0.0);
    // The water molecules bind
    EXPECT_LT(deltaU, -5);
    EXPECT_GT(deltaU, -40);

    // The energy file
    const auto xvg = readXvg(files.energy);
    ASSERT_EQ(40U, xvg.rows.size());
    for(const std::string name : { "BONDS", "INDUCTION", "EPOT", "WALL", "EKIN", "ETOT",
                                   "TEMPERATURE", "RGYR", "DISTANCE" })
    {
        ASSERT_GT(xvg.column(name), 0) << name;
    }
    EXPECT_EQ("DISTANCE", xvg.legends.back());
    double induction = 0;
    for(const auto &row : xvg.rows)
    {
        EXPECT_NEAR(row[xvg.column("EPOT")] + row[xvg.column("EKIN")], row[xvg.column("ETOT")], 2e-5);
        EXPECT_EQ(0.0, row[xvg.column("WALL")]);
        EXPECT_NEAR(0.29, row[xvg.column("DISTANCE")], 0.03);
        induction += row[xvg.column("INDUCTION")]/xvg.rows.size();
    }
    // All beads start at the minimum, and spread
    EXPECT_EQ(0.0, xvg.rows[0][xvg.column("RGYR")]);
    EXPECT_GT(xvg.rows.back()[xvg.column("RGYR")], 0.005);
    EXPECT_LT(induction, 0);

    // The trajectory has the centroids of the atoms, four frames of six atoms
    const auto frames = readPdbFrames(files.trajectory);
    ASSERT_EQ(4U, frames.size());
    for(const auto &frame : frames)
    {
        EXPECT_EQ(6U, frame.size());
    }

    // The progress file, as without beads, and with the beads in addition
    std::ifstream               in(progress);
    ASSERT_TRUE(in.good()) << "Cannot open " << progress;
    std::vector<nlohmann::json> lines;
    for (std::string line; std::getline(in, line); )
    {
        const auto js = nlohmann::json::parse(line, nullptr, false);
        ASSERT_TRUE(js.is_object()) << "Not a JSON object: " << line;
        lines.push_back(js);
    }
    ASSERT_EQ(xvg.rows.size(), lines.size());
    for (size_t i = 0; i < lines.size(); i++)
    {
        const auto &js = lines[i];
        EXPECT_EQ(static_cast<int>(10*i), js["step"].get<int>());
        EXPECT_EQ(400, js["nsteps"].get<int>());
        EXPECT_EQ(4, js["beads"].get<int>());
        EXPECT_EQ("nm", js["gyration"]["unit"].get<std::string>());
        EXPECT_NEAR(xvg.rows[i][xvg.column("RGYR")], js["gyration"]["value"].get<double>(), 1e-6);
        EXPECT_NEAR(xvg.rows[i][xvg.column("DISTANCE")], js["distance"]["value"].get<double>(), 1e-6);
        EXPECT_NEAR(xvg.rows[i][xvg.column("TEMPERATURE")], js["temperature"]["value"].get<double>(), 1e-4);
        for (const std::string term : { "EPOT", "EKIN", "ETOT", "WALL", "INDUCTION" })
        {
            ASSERT_TRUE(js["energies"].contains(term)) << term;
            EXPECT_NEAR(xvg.rows[i][xvg.column(term)], js["energies"][term].get<double>(), 1e-4) << term;
        }
    }

    // The histogram of the distance
    const auto hxvg = readXvg(hist);
    ASSERT_FALSE(hxvg.rows.empty());
    double total = 0;
    for(const auto &row : hxvg.rows)
    {
        total += row[1]*0.002;
    }
    EXPECT_NEAR(1.0, total, 1e-4);
}

TEST_F(SimulateBeadsTest, VirtualSites)
{
    // A water dimer with virtual sites runs, and has no file that was not asked for
    SimulateFiles files;
    auto cmd = baseCommand("water_dimer.pdb", "beads-vsites", &files,
                           alexandriaFile("forcefields/TIP5P.xml"));
    cmd.addOption("-minimize");
    cmd.addOption("-nbeads", 8);
    cmd.addOption("-temp", 100.0);
    cmd.addOption("-nsteps", 400);
    cmd.addOption("-deltat", 0.00025);
    cmd.addOption("-nstener", 20);
    cmd.addOption("-seed", 1993);
    std::remove("distance.xvg");
    const auto res = runBeads(&cmd, "beads-vsites");
    ASSERT_TRUE(res.contains("all compounds")) << res.dump();
    EXPECT_FALSE(res.contains("Binding"));
    EXPECT_FALSE(res["all compounds"].contains("Fraction at the wall"));
    EXPECT_NEAR(0.28, average(res["all compounds"], "Distance"), 0.03);
    EXPECT_FALSE(std::ifstream("distance.xvg").good());
    const auto xvg = readXvg(files.energy);
    ASSERT_EQ(20U, xvg.rows.size());
    // No wall, no column for it
    EXPECT_EQ(-1, xvg.column("WALL"));
    EXPECT_GT(xvg.column("DISTANCE"), 0);
    // The trajectory has the atoms only
    const auto frames = readPdbFrames(files.trajectory);
    ASSERT_EQ(400U, frames.size());
    EXPECT_EQ(6U, frames[0].size());
}

TEST_F(SimulateBeadsTest, Reproducible)
{
    std::vector<SimulateFiles> files(3);
    for(int i = 0; i < 3; i++)
    {
        auto cmd = baseCommand("water_dimer.pdb", "beads-repro" + std::to_string(i), &files[i]);
        cmd.addOption("-nbeads", 4);
        cmd.addOption("-temp", 100.0);
        cmd.addOption("-nsteps", 100);
        cmd.addOption("-nstener", 10);
        cmd.addOption("-nstxout", 50);
        cmd.addOption("-seed", i < 2 ? 7 : 8);
        ASSERT_EQ(0, run(&cmd));
    }
    EXPECT_EQ(textWithoutComments(files[0].energy), textWithoutComments(files[1].energy));
    EXPECT_EQ(textWithoutComments(files[0].trajectory), textWithoutComments(files[1].trajectory));
    EXPECT_NE(textWithoutComments(files[0].energy), textWithoutComments(files[2].energy));
}

TEST_F(SimulateBeadsTest, WithoutCommandLineSameAsTheProgram)
{
    SimulateFiles files;
    auto cmd = baseCommand("water_dimer.pdb", "beads-tool", &files);
    cmd.addOption("-nbeads", 4);
    cmd.addOption("-temp", 100.0);
    cmd.addOption("-nsteps", 100);
    cmd.addOption("-deltat", 0.0002);
    cmd.addOption("-nstener", 10);
    cmd.addOption("-nstxout", 50);
    cmd.addOption("-seed", 1993);
    cmd.addOption("-tau", 0.1);
    cmd.addOption("-equil", 20);
    cmd.addOption("-rwall", 0.25);
    cmd.addOption("-kwall", 500.0);
    ASSERT_EQ(0, run(&cmd));

    SimulateOptions options;
    options.forceFieldFile = alexandriaFile("forcefields/ACS-g1.xml");
    options.compounds.setFileName(alexandriaFile("molecules/water_dimer.pdb"));
    options.energyFile     = fileManager().getTemporaryFilePath("beads-direct-energy.xvg");
    options.trajectoryFile = fileManager().getTemporaryFilePath("beads-direct-traj.pdb");
    options.simulation.setSimulation(100, 0.0002, 100, 1993, 10, 50);
    options.simulation.setBeads(4, 0.1, 20);
    options.simulation.setWall(0.25, 500);
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    EXPECT_EQ(0, alexandria::runSimulate(&msghandler, options));
    EXPECT_EQ(textWithoutComments(files.energy), textWithoutComments(options.energyFile));
    EXPECT_EQ(textWithoutComments(files.trajectory), textWithoutComments(options.trajectoryFile));
    // The wall is closer than the molecules are, and pushes them together
    const auto xvg = readXvg(files.energy);
    ASSERT_GT(xvg.column("WALL"), 0);
    EXPECT_GT(xvg.rows[0][xvg.column("WALL")], 0.0);
}

TEST_F(SimulateBeadsTest, MovedApart)
{
    // Two argon atoms at 300 K without a wall do not stay together
    SimulateFiles files;
    auto cmd = baseCommand("argon#argon.pdb", "apart", &files,
                           alexandriaFile("forcefields/noble-halgren14-7.xml"));
    cmd.addOption("-nbeads", 1);
    cmd.addOption("-temp", 300.0);
    cmd.addOption("-nsteps", 20000);
    cmd.addOption("-deltat", 0.005);
    cmd.addOption("-tau", 1.0);
    cmd.addOption("-nstener", 0);
    cmd.addOption("-nstxout", 0);
    cmd.addOption("-seed", 1993);
    cmd.addOption("-monomers");
    const auto res = runBeads(&cmd, "apart");
    ASSERT_TRUE(res.contains("all compounds")) << res.dump();
    ASSERT_TRUE(res["all compounds"].contains("Moved apart")) << res.dump();
    EXPECT_TRUE(res["all compounds"]["Moved apart"].get<bool>());
    EXPECT_GT(average(res["all compounds"], "Distance"), 1.5);
    // The log file says so, and that the binding energy has no meaning
    std::ifstream in(files.log);
    const std::string log((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(std::string::npos, log.find("The compounds moved apart"));
    EXPECT_NE(std::string::npos, log.find("not bound all the time"));
}

TEST_F(SimulateBeadsTest, ThermochemistryHydrogenFluoride)
{
    // Hydrogen fluoride with the anharmonic bond of hh-intra.xml at
    // 298.15 K. The levels of the vibration and the rotation in this
    // potential, from the radial Schrodinger equation on a grid for each
    // angular momentum (NumPy, not part of ACT), give exactly, relative to
    // the minimum: E = 30.778 kJ/mol, cV = 20.82 J/mol K and
    // S = 173.83 J/mol K at 1 bar. The rigid rotor with a harmonic
    // oscillator has 31.118, 20.786 and 173.639. With 16 beads harmonic
    // oscillators lack 3.8 kJ/mol and have 20 J/mol K too much, which is
    // corrected for.
    SimulateFiles files;
    auto cmd = baseCommand("hydrogen-fluoride.sdf", "thermo-hf", &files,
                           alexandriaFile("forcefields/hh-intra.xml"));
    const auto progress = fileManager().getTemporaryFilePath("thermo-hf.jsonl");
    cmd.addOption("-nbeads", 16);
    cmd.addOption("-temp", 298.15);
    cmd.addOption("-nsteps", 80000);
    cmd.addOption("-tau", 0.01);
    cmd.addOption("-nstener", 2000);
    cmd.addOption("-nstxout", 0);
    cmd.addOption("-seed", 1993);
    cmd.addOption("-thermo");
    cmd.addOption("-progress", progress);
    const auto res = runBeads(&cmd, "thermo-hf");
    ASSERT_TRUE(res.contains("Thermochemistry")) << res.dump();
    const auto &tc = res["Thermochemistry"];
    EXPECT_EQ(16, tc["Beads"].get<int>());
    EXPECT_EQ(1, tc["Symmetry number"].get<int>());
    // The temperature is a real number
    EXPECT_NEAR(298.15, tc["Temperature"]["value"].get<double>(), 1e-4);
    EXPECT_NEAR(-595.587, tc["Minimum energy"]["value"].get<double>(), 0.01);
    // The normal mode analysis
    EXPECT_NEAR(24.921, tc["Zero point energy"]["Harmonic"].get<double>(), 0.01);
    EXPECT_FALSE(tc["Zero point energy"].contains("Average"));
    EXPECT_NEAR(31.118, tc["Internal energy"]["Harmonic"].get<double>(), 0.01);
    EXPECT_NEAR(20.786, tc["Heat capacity cV"]["Harmonic"].get<double>(), 0.01);
    EXPECT_NEAR(173.639, tc["Standard entropy"]["Harmonic"].get<double>(), 0.02);
    EXPECT_NEAR(3.82, tc["Bead correction"]["Energy"]["value"].get<double>(), 0.02);
    EXPECT_NEAR(-20.03, tc["Bead correction"]["Heat capacity"]["value"].get<double>(), 0.05);
    // The simulations, within about six standard errors
    EXPECT_NEAR(30.778, average(tc, "Internal energy"), 0.5);
    EXPECT_NEAR(average(tc, "Internal energy") + BOLTZ*298.15, average(tc, "Enthalpy"), 1e-4);
    EXPECT_NEAR(20.82, average(tc, "Heat capacity cV"), 8.0);
    EXPECT_NEAR(average(tc, "Heat capacity cV") + RGAS, average(tc, "Heat capacity cP"), 1e-6);
    EXPECT_NEAR(173.83, average(tc, "Standard entropy"), 2.0);
    EXPECT_NEAR(average(tc, "Enthalpy") - 298.15*average(tc, "Standard entropy")/KILO,
                average(tc, "Gibbs energy"), 1e-4);
    for(const std::string name : { "Internal energy", "Heat capacity cV", "Standard entropy", "Gibbs energy" })
    {
        EXPECT_GT(error(tc, name), 0.0) << name;
        EXPECT_EQ("J/mol K" == tc[name]["Unit"].get<std::string>(), name.find("energy") == std::string::npos) << name;
    }
    EXPECT_LT(error(tc, "Internal energy"), 0.2);
    EXPECT_LT(error(tc, "Standard entropy"), 1.0);
    // The internal energy is the total energy with the correction, relative to the minimum
    EXPECT_NEAR(average(res["the compound"], "Total energy") - tc["Minimum energy"]["value"].get<double>() +
                tc["Bead correction"]["Energy"]["value"].get<double>(),
                average(tc, "Internal energy"), 1e-6);
    // The quantum effects are close to those of the harmonic oscillator,
    // the anharmonicity of the classical molecule is small
    ASSERT_TRUE(tc.contains("Free energy")) << tc.dump();
    const auto &fe = tc["Free energy"];
    EXPECT_NEAR(17.481, fe["Quantum effects"]["Harmonic"].get<double>(), 0.01);
    EXPECT_NEAR(17.15, average(fe, "Quantum effects"), 0.3);
    EXPECT_NEAR(-0.05, average(fe, "Anharmonicity"), 0.4);
    // The progress: the simulation, then a line for each of the eight
    // simulations of the entropy, of which the one with one bead is not done
    std::ifstream            in(progress);
    std::vector<std::string> entropy;
    size_t                   nenergy = 0;
    for(std::string line; std::getline(in, line); )
    {
        const auto js = nlohmann::json::parse(line);
        if (js.contains("entropy"))
        {
            entropy.push_back(js["entropy"].dump());
        }
        else
        {
            // All lines of the simulation come first
            EXPECT_TRUE(entropy.empty());
            nenergy++;
        }
    }
    EXPECT_EQ(40U, nenergy);
    ASSERT_EQ(8U, entropy.size());
    EXPECT_EQ(R"({"done":1,"total":8})", entropy.front());
    EXPECT_EQ(R"({"done":8,"total":8})", entropy.back());
    // The log file has the table
    EXPECT_TRUE(fileContains(files.log, "Thermochemistry at 298.15 K and 1 bar from simulations with 16 beads"));
    EXPECT_TRUE(fileContains(files.log, "Standard entropy (J/mol K)"));
}

TEST_F(SimulateBeadsTest, ThermochemistryClassicalAndWithoutEntropy)
{
    // With one bead the thermochemistry is classical: water has 6 kT of
    // energy and a heat capacity of 6 k, and only the anharmonicity is
    // computed. The harmonic values next to it are quantum mechanical.
    const double T = 300;
    {
        SimulateFiles files;
        auto cmd = baseCommand("water-3-oep.log.pdb", "thermo-classical", &files);
        cmd.addOption("-nbeads", 1);
        cmd.addOption("-temp", T);
        cmd.addOption("-nsteps", 100000);
        cmd.addOption("-tau", 0.01);
        cmd.addOption("-nstener", 0);
        cmd.addOption("-nstxout", 0);
        cmd.addOption("-seed", 1993);
        cmd.addOption("-thermo");
        cmd.addOption("-nquad", 3);
        const auto res = runBeads(&cmd, "thermo-classical");
        ASSERT_TRUE(res.contains("Thermochemistry")) << res.dump();
        const auto &tc = res["Thermochemistry"];
        EXPECT_EQ(1, tc["Beads"].get<int>());
        EXPECT_EQ(2, tc["Symmetry number"].get<int>());
        EXPECT_FALSE(tc.contains("Bead correction"));
        EXPECT_NEAR(6*BOLTZ*T, average(tc, "Internal energy"), 0.5);
        EXPECT_NEAR(6*RGAS, average(tc, "Heat capacity cV"), 8.0);
        // The quantum oscillators have a zero point energy and hardly a heat capacity
        EXPECT_GT(tc["Internal energy"]["Harmonic"].get<double>(), 40.0);
        EXPECT_GT(tc["Heat capacity cV"]["Harmonic"].get<double>(), 3*RGAS);
        EXPECT_LT(tc["Heat capacity cV"]["Harmonic"].get<double>(), 3*RGAS + 1.0);
        ASSERT_TRUE(tc.contains("Free energy")) << tc.dump();
        EXPECT_FALSE(tc["Free energy"].contains("Quantum effects"));
        EXPECT_NEAR(0.0, average(tc["Free energy"], "Anharmonicity"), 0.5);
        // and hardly an entropy, where that of a classical oscillator,
        // R(1 - ln(h nu/kT)), is negative for these frequencies: -40 J/mol K
        // for the bend and the two stretches together
        const double sharmonic = tc["Standard entropy"]["Harmonic"].get<double>();
        EXPECT_NEAR(sharmonic - 40, average(tc, "Standard entropy"), 5.0);
        EXPECT_TRUE(fileContains(files.log, "With one bead the simulation is classical"));
    }
    {
        // Without the simulations of the integrals there is no entropy.
        // The structure is minimized also without -minimize.
        SimulateFiles files;
        auto cmd = baseCommand("water-3-oep.log.pdb", "thermo-nquad", &files);
        cmd.addOption("-nbeads", 8);
        cmd.addOption("-temp", T);
        cmd.addOption("-nsteps", 4000);
        cmd.addOption("-nstener", 0);
        cmd.addOption("-nstxout", 0);
        cmd.addOption("-seed", 1993);
        cmd.addOption("-thermo");
        cmd.addOption("-nquad", 0);
        const auto res = runBeads(&cmd, "thermo-nquad");
        ASSERT_TRUE(res.contains("Thermochemistry")) << res.dump();
        const auto &tc = res["Thermochemistry"];
        EXPECT_FALSE(tc.contains("Free energy"));
        EXPECT_TRUE(tc["Standard entropy"].contains("Harmonic"));
        EXPECT_FALSE(tc["Standard entropy"].contains("Average"));
        EXPECT_TRUE(tc["Internal energy"].contains("Average"));
        EXPECT_TRUE(tc.contains("Bead correction"));
        EXPECT_FALSE(fileContains(files.log, "Entropy:"));
    }
    {
        // The thermochemistry is that of one compound
        SimulateFiles files;
        auto cmd = baseCommand("water_dimer.pdb", "thermo-dimer", &files);
        cmd.addOption("-nbeads", 4);
        cmd.addOption("-temp", 100.0);
        cmd.addOption("-nsteps", 200);
        cmd.addOption("-nstener", 0);
        cmd.addOption("-nstxout", 0);
        cmd.addOption("-seed", 1993);
        cmd.addOption("-thermo");
        const auto res = runBeads(&cmd, "thermo-dimer");
        EXPECT_TRUE(res.contains("all compounds"));
        EXPECT_FALSE(res.contains("Thermochemistry"));
        EXPECT_TRUE(fileContains(files.log, "the thermochemistry (-thermo) is that of one compound"));
    }
}

TEST_F(SimulateBeadsTest, WrongOptions)
{
    const auto beads = [this](const std::string &label, SimulateFiles *files)
    {
        auto cmd = baseCommand("water_dimer.pdb", label, files);
        cmd.addOption("-nsteps", 20);
        cmd.addOption("-temp", 100.0);
        return cmd;
    };
    SimulateFiles files;
    {
        auto cmd = beads("wrong-nbeads", &files);
        cmd.addOption("-nbeads", -1);
        expectRefused(&cmd, "a negative number of beads");
    }
    {
        // No temperature
        auto cmd = baseCommand("water_dimer.pdb", "wrong-temp", &files);
        cmd.addOption("-nsteps", 20);
        cmd.addOption("-nbeads", 4);
        expectRefused(&cmd, "no temperature");
    }
    {
        auto cmd = beads("wrong-constrain", &files);
        cmd.addOption("-nbeads", 4);
        cmd.addOption("-constrain", "bonds");
        expectRefused(&cmd, "constraints");
    }
    {
        auto cmd = beads("wrong-nstvout", &files);
        cmd.addOption("-nbeads", 4);
        cmd.addOption("-nstvout", 10);
        expectRefused(&cmd, "velocities");
    }
    {
        auto cmd = beads("wrong-equil", &files);
        cmd.addOption("-nbeads", 4);
        cmd.addOption("-equil", 20);
        expectRefused(&cmd, "equilibration as long as the simulation");
    }
    {
        auto cmd = beads("wrong-tau", &files);
        cmd.addOption("-nbeads", 4);
        cmd.addOption("-tau", 0.0);
        expectRefused(&cmd, "no time constant");
    }
    {
        auto cmd = beads("wrong-rwall", &files);
        cmd.addOption("-nbeads", 4);
        cmd.addOption("-rwall", -1.0);
        expectRefused(&cmd, "a negative radius");
    }
    {
        // A wall and the binding energy need beads
        auto cmd = beads("wrong-wall", &files);
        cmd.addOption("-rwall", 1.0);
        expectRefused(&cmd, "a wall without beads");
        auto cmd2 = beads("wrong-monomers", &files);
        cmd2.addOption("-monomers");
        expectRefused(&cmd2, "monomers without beads");
    }
    {
        // The thermochemistry needs beads, and one compound without a wall
        auto cmd = beads("wrong-thermo", &files);
        cmd.addOption("-thermo");
        expectRefused(&cmd, "thermochemistry without beads");
        auto cmd2 = beads("wrong-thermo-wall", &files);
        cmd2.addOption("-nbeads", 4);
        cmd2.addOption("-thermo");
        cmd2.addOption("-rwall", 1.0);
        expectRefused(&cmd2, "thermochemistry with a wall");
        auto cmd3 = beads("wrong-thermo-monomers", &files);
        cmd3.addOption("-nbeads", 4);
        cmd3.addOption("-thermo");
        cmd3.addOption("-monomers");
        expectRefused(&cmd3, "thermochemistry with the binding energy");
        auto cmd4 = beads("wrong-nquad", &files);
        cmd4.addOption("-nbeads", 4);
        cmd4.addOption("-thermo");
        cmd4.addOption("-nquad", -1);
        expectRefused(&cmd4, "a negative number of simulations");
        auto cmd5 = beads("wrong-nquad-large", &files);
        cmd5.addOption("-nbeads", 4);
        cmd5.addOption("-thermo");
        cmd5.addOption("-nquad", 17);
        expectRefused(&cmd5, "too many simulations");
    }
    {
        // The springs of 64 beads at 300 K need a time step below 0.4 fs
        auto cmd = baseCommand("water_dimer.pdb", "wrong-deltat", &files);
        cmd.addOption("-nsteps", 20);
        cmd.addOption("-temp", 300.0);
        cmd.addOption("-nbeads", 64);
        cmd.addOption("-deltat", 0.002);
        expectRefused(&cmd, "a time step too large for the beads");
        EXPECT_FALSE(fileNotEmpty(files.energy));
    }
    // The settings say which option is wrong, without a command line as well
    const auto message = [](const SimulationConfigHandler &sch)
    {
        try
        {
            sch.validate();
        }
        catch (const gmx::InvalidInputError &ex)
        {
            return std::string(ex.what());
        }
        return std::string();
    };
    {
        SimulationConfigHandler sch;
        sch.setSimulation(20, 0.0002, 0, 1, 1, 1);
        sch.setBeads(4, 0.05, -1);
        EXPECT_NE(std::string::npos, message(sch).find("-temp")) << message(sch);
        sch.setSimulation(20, 0.0002, 100, 1, 1, 1);
        EXPECT_EQ("", message(sch));
        sch.setBeads(4, 0.05, 20);
        EXPECT_NE(std::string::npos, message(sch).find("-equil")) << message(sch);
        sch.setBeads(4, 0.05, -2);
        EXPECT_NE(std::string::npos, message(sch).find("-equil")) << message(sch);
        sch.setBeads(4, -0.05, -1);
        EXPECT_NE(std::string::npos, message(sch).find("-tau")) << message(sch);
        sch.setBeads(0, 0.05, -1);
        sch.setWall(1.0, 1000);
        EXPECT_NE(std::string::npos, message(sch).find("-rwall")) << message(sch);
        sch.setWall(0, 1000);
        sch.setMonomers(true);
        EXPECT_NE(std::string::npos, message(sch).find("-monomers")) << message(sch);
        sch.setMonomers(false);
        sch.setThermochemistry(true, 4);
        EXPECT_NE(std::string::npos, message(sch).find("-thermo")) << message(sch);
        sch.setBeads(4, 0.05, -1);
        sch.setSimulation(20, 0.0002, 100, 1, 1, 1);
        EXPECT_EQ("", message(sch));
        EXPECT_TRUE(sch.thermochemistry());
        EXPECT_EQ(4, sch.nquad());
        sch.setWall(1.0, 1000);
        EXPECT_NE(std::string::npos, message(sch).find("-thermo")) << message(sch);
        sch.setWall(0, 1000);
        sch.setThermochemistry(true, 17);
        EXPECT_NE(std::string::npos, message(sch).find("-nquad")) << message(sch);
        // Without steps there is no simulation that needs a temperature
        SimulationConfigHandler none;
        none.setBeads(4, 0.05, -1);
        EXPECT_EQ("", message(none));
    }
}

} // namespace

} // namespace alexandria
