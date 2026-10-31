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
#include <map>
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
        // The energies are written with five decimals, a tighter
        // tolerance would not make a difference
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
                    bonds.push_back({ i, j });
                    neighbours[i].push_back(j);
                    neighbours[j].push_back(i);
                }
            }
        }
        ASSERT_EQ(9u, bonds.size()) << "acetone has nine bonds";
        // The largest change of a distance along the trajectory
        auto maxChange = [&](const std::vector<std::pair<size_t, size_t>> &pairs)
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
    std::string molFile("methanethiol.sdf");
    auto distance = [](const std::vector<std::vector<double>> &x, int i, int j)
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
            auto label   = "min" + constrain + alg;
            auto cmd     = baseCommand(molFile, label, &files);
            auto confout = fileManager().getTemporaryFilePath(label + "-confout.pdb");
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
    std::string   molFile("acetone-3-oep.log.pdb");
    SimulateFiles files;
    auto cmd     = baseCommand(molFile, "freezecons", &files);
    auto confout = fileManager().getTemporaryFilePath("freezecons-confout.pdb");
    auto ndx     = fileManager().getTemporaryFilePath("freezecons.ndx");
    {
        std::ofstream out(ndx);
        out << "[ Frozen ]\n1 2\n";
    }
    cmd.addOption("-minimize");
    cmd.addOption("-c", confout);
    cmd.addOption("-freeze", ndx);
    cmd.addOption("-constrain", "angles");
    ASSERT_EQ(0, run(&cmd));
    auto xin  = readPdbCoordinates(alexandriaFile("mols/" + molFile));
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
        std::string   molFile("water_dimer.pdb");
        SimulateFiles files;
        auto cmd     = baseCommand(molFile, "simmin" + constrain, &files);
        auto simout  = fileManager().getTemporaryFilePath("simmin" + constrain + "-confout.pdb");
        cmd.addOption("-minimize");
        cmd.addOption("-constrain", constrain);
        cmd.addOption("-c", simout);
        ASSERT_EQ(0, run(&cmd));

        auto minout  = fileManager().getTemporaryFilePath("minimize" + constrain + "-confout.pdb");
        auto mincmd  = minimizeCommand(alexandriaFile("ACS-g.xml"),
                                       alexandriaFile("mols/" + molFile), minout,
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

// Test: alexandria minimize fails for a missing molecule
TEST_F(SimulateTest, MinimizeToolMissingMolecule)
{
    auto cmd = minimizeCommand(alexandriaFile("ACS-g.xml"),
                               alexandriaFile("mols/does-not-exist.sdf"),
                               fileManager().getTemporaryFilePath("missing-confout.pdb"),
                               fileManager().getTemporaryFilePath("minimize-missing.log"));
    EXPECT_EQ(1, minimize(cmd.argc(), cmd.argv()));
}

} // namespace

} // namespace alexandria
