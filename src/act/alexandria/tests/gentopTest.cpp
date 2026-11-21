/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2023,2026
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
 * Tests of the gentop module, in particular the GROMACS topology it writes.
 */

#include "actpre.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <map>
#include <nlohmann/json.hpp>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/gentop.h"
#include "act/basics/chargemodel.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/molprop/experiment.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "act/testutils/output_checks.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! A line of a section of a topology, without the comment
struct TopLine
{
    //! The fields
    std::vector<std::string> fields;
    //! The comment after the semicolon, without spaces around it
    std::string              comment;
};

//! The sections of a topology, in order, sections may occur more than once
using Topology = std::vector<std::pair<std::string, std::vector<TopLine>>>;

//! \return the topology in a file
Topology readTopology(const std::string &fileName)
{
    Topology      top;
    std::ifstream in(fileName);
    std::string   line;
    while (std::getline(in, line))
    {
        auto        semicolon = line.find(';');
        const std::string content   = line.substr(0, semicolon);
        const std::string comment   = semicolon == std::string::npos ? "" : line.substr(semicolon + 1);
        auto        open      = content.find('[');
        if (open != std::string::npos)
        {
            std::istringstream name(content.substr(open + 1, content.find(']') - open - 1));
            std::string        section;
            name >> section;
            top.push_back({ section, {} });
            continue;
        }
        std::istringstream words(content);
        TopLine            tl;
        for (std::string w; words >> w; )
        {
            tl.fields.push_back(w);
        }
        if (tl.fields.empty() || top.empty())
        {
            continue;
        }
        std::istringstream c(comment);
        c >> tl.comment;
        top.back().second.push_back(tl);
    }
    return top;
}

//! \return all lines of the sections with a name
std::vector<TopLine> section(const Topology &top, const std::string &name)
{
    std::vector<TopLine> lines;
    for (const auto &s : top)
    {
        if (s.first == name)
        {
            lines.insert(lines.end(), s.second.begin(), s.second.end());
        }
    }
    return lines;
}

//! \return the excluded pairs of particles, each with the lowest index first
std::set<std::pair<int, int>> exclusionPairs(const Topology &top)
{
    std::set<std::pair<int, int>> pairs;
    for (const auto &e : section(top, "exclusions"))
    {
        int i = std::stoi(e.fields[0]);
        for (size_t k = 1; k < e.fields.size(); k++)
        {
            int j = std::stoi(e.fields[k]);
            pairs.insert({ std::min(i, j), std::max(i, j) });
        }
    }
    return pairs;
}

/*! \brief The pairs that GROMACS excludes for nrexcl, extended to the shells
 * \return the pairs of particles at most nrexcl bonds apart, where a
 *         shell has the place of its core, each with the lowest index first
 */
std::set<std::pair<int, int>> expectedExclusionPairs(const Topology &top)
{
    int  nrexcl = std::stoi(section(top, "moleculetype")[0].fields[1]);
    auto natoms = section(top, "atoms").size();
    // The core of each particle, starting at 1 as in the topology
    std::vector<int> core(natoms + 1);
    for (size_t i = 0; i <= natoms; i++)
    {
        core[i] = i;
    }
    for (const auto &p : section(top, "polarization"))
    {
        core[std::stoi(p.fields[1])] = std::stoi(p.fields[0]);
    }
    std::vector<std::vector<int>> neighbors(natoms + 1);
    for (const auto &b : section(top, "bonds"))
    {
        int i = std::stoi(b.fields[0]);
        int j = std::stoi(b.fields[1]);
        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
    }
    std::set<std::pair<int, int>> pairs;
    for (size_t i = 1; i <= natoms; i++)
    {
        // The particles at most nrexcl bonds from the core of i
        std::vector<int> distance(natoms + 1, -1);
        std::vector<int> front = { core[i] };
        distance[core[i]] = 0;
        for (int n = 1; n <= nrexcl; n++)
        {
            std::vector<int> next;
            for (const auto j : front)
            {
                for (const auto k : neighbors[j])
                {
                    if (distance[k] < 0)
                    {
                        distance[k] = n;
                        next.push_back(k);
                    }
                }
            }
            front = next;
        }
        for (size_t j = i + 1; j <= natoms; j++)
        {
            if (distance[core[j]] >= 0)
            {
                pairs.insert({ i, j });
            }
        }
    }
    return pairs;
}

//! \return the parameters of an identifier, given as text
const ForceFieldParameterMap &parameters(const ForceField &pd, InteractionType itype, const std::string &id)
{
    for (const auto &p : pd.findForcesConst(itype).parametersConst())
    {
        if (p.first.id() == id)
        {
            return p.second;
        }
    }
    GMX_THROW(gmx::InvalidInputError("No parameters for " + id));
}

//! A cube file, the lengths are in Bohr
struct Cube
{
    //! Origin of the grid
    double                           origin[3];
    //! Number of grid points in each direction
    int                              n[3];
    //! Spacing in each direction
    double                           space[3];
    //! Atomic number, charge and coordinates of each atom
    std::vector<std::array<double, 5>> atoms;
    //! The values on the grid, with z running fastest
    std::vector<double>              values;
};

//! \return the cube in a file
Cube readCube(const std::string &fileName)
{
    Cube          cube;
    std::ifstream in(fileName);
    std::string   line;
    // Title and type
    std::getline(in, line);
    std::getline(in, line);
    int natoms;
    in >> natoms >> cube.origin[0] >> cube.origin[1] >> cube.origin[2];
    for (int m = 0; m < 3; m++)
    {
        double zero;
        in >> cube.n[m];
        for (int k = 0; k < 3; k++)
        {
            in >> (k == m ? cube.space[m] : zero);
        }
    }
    cube.atoms.resize(natoms);
    for (auto &a : cube.atoms)
    {
        for (auto &f : a)
        {
            in >> f;
        }
    }
    for (double v; in >> v; )
    {
        cube.values.push_back(v);
    }
    return cube;
}

class GentopTest : public gmx::test::CommandLineTestBase
{
protected:
    /*! \brief Run gentop
     * \param[in] ffName   Force field in this directory
     * \param[in] molecule Structure in mols
     * \param[in] options  Output files and other options, the file names are
     *                     made temporary
     */
    void run(const std::string &ffName, const std::string &molecule,
             std::map<std::string, std::string> *files)
    {
        run(ffName, { { "-f", fileManager().getInputFilePath("molecules/" + molecule) } }, files);
    }

    /*! \brief Run gentop
     * \param[in] ffName  Force field in this directory
     * \param[in] input   Options that are passed as they are
     * \param[in] options Output files, the file names are made temporary
     */
    void run(const std::string &ffName, const std::vector<std::pair<std::string, std::string>> &input,
             std::map<std::string, std::string> *files)
    {
        ASSERT_EQ(0, status(ffName, input, files));
    }

    /*! \brief Run gentop
     * \param[in] ffName  Force field in this directory
     * \param[in] input   Options that are passed as they are
     * \param[in] options Output files, the file names are made temporary
     * \return the exit status of gentop
     */
    int status(const std::string &ffName, const std::vector<std::pair<std::string, std::string>> &input,
               std::map<std::string, std::string> *files)
    {
        gmx::test::CommandLine cmd;
        cmd.append("gentop");
        cmd.addOption("-ff", fileManager().getInputFilePath(ffName));
        for (const auto &i : input)
        {
            cmd.addOption(i.first.c_str(), i.second);
        }
        cmd.addOption("-g", fileManager().getTemporaryFilePath("gentop.log"));
        for (auto &f : *files)
        {
            f.second = fileManager().getTemporaryFilePath(f.second);
            cmd.addOption(f.first.c_str(), f.second);
        }
        return gentop(cmd.argc(), cmd.argv());
    }

    /*! \brief Check the bonded sections of a topology
     * Each line has the atoms, the function type and the number of
     * parameters of that GROMACS function, and none is unsupported.
     * \param[in] top The topology
     */
    static void checkBondeds(const Topology &top)
    {
        // Number of atoms of each section and parameters of each function
        std::map<std::string, size_t> natoms = { { "bonds", 2 }, { "angles", 3 }, { "dihedrals", 4 } };
        std::map<std::pair<std::string, int>, size_t> nparams = {
            { { "bonds", 1 }, 2 }, { { "bonds", 3 }, 3 },
            { { "angles", 1 }, 2 }, { { "angles", 5 }, 4 }, { { "angles", 9 }, 2 },
            { { "dihedrals", 1 }, 3 }, { { "dihedrals", 2 }, 2 }, { { "dihedrals", 3 }, 6 }
        };
        for (const auto &s : top)
        {
            if (natoms.count(s.first) == 0)
            {
                continue;
            }
            size_t na = natoms[s.first];
            for (const auto &l : s.second)
            {
                ASSERT_GT(l.fields.size(), na) << s.first;
                int  type = std::stoi(l.fields[na]);
                const auto np   = nparams.find({ s.first, type });
                ASSERT_TRUE(np != nparams.end()) << s.first << " function type " << type;
                EXPECT_EQ(na + 1 + np->second, l.fields.size()) << s.first << " " << l.comment;
            }
        }
    }
};

// Test: The topology of acetone with harmonic bonds and angles and proper
// dihedrals, the parameters are those of the force field in the units of
// GROMACS, as are the van der Waals parameters of the pairs of atom types.
// With GROMACS this gives the energies of ACT, apart from the constant
// bond energies.
TEST_F(GentopTest, GromacsTopology)
{
    std::map<std::string, std::string> files = { { "-p", "acetone.top" }, { "-c", "acetone.pdb" } };
    run("forcefields/ACS-g1.xml", "acetone.sdf", &files);
    const auto top = readTopology(files["-p"]);
    checkBondeds(top);
    for (const char *s : { "defaults", "atomtypes", "moleculetype", "atoms", "bonds", "angles", "dihedrals" })
    {
        EXPECT_FALSE(section(top, s).empty()) << s;
    }
    auto atoms = section(top, "atoms");
    ASSERT_EQ(10u, atoms.size());
    double qtot = 0;
    for (const auto &a : atoms)
    {
        qtot += std::stod(a.fields[6]);
    }
    // The charges have six significant digits
    EXPECT_NEAR(0, qtot, 1e-5);
    EXPECT_EQ(9u, section(top, "bonds").size());
    // The same exclusions as GROMACS generates from the bonds
    EXPECT_FALSE(exclusionPairs(top).empty());
    EXPECT_EQ(expectedExclusionPairs(top), exclusionPairs(top));
    ForceField pd;
    readForceField(fileManager().getInputFilePath("forcefields/ACS-g1.xml"), &pd);
    const auto angles = section(top, "angles");
    EXPECT_EQ(15u, angles.size());
    for (const auto &a : angles)
    {
        auto &p = parameters(pd, InteractionType::ANGLES, a.comment);
        EXPECT_NEAR(p.find("angle")->second.value(), std::stod(a.fields[4]), 1e-4) << a.comment;
        EXPECT_NEAR(p.find("kt")->second.value(), std::stod(a.fields[5]), 1e-4) << a.comment;
    }
    for (const auto &d : section(top, "dihedrals"))
    {
        if (std::stoi(d.fields[4]) != 1)
        {
            continue;
        }
        auto &p = parameters(pd, InteractionType::PROPER_DIHEDRALS, d.comment);
        EXPECT_NEAR(p.find("phi0")->second.value(), std::stod(d.fields[5]), 1e-4) << d.comment;
        EXPECT_NEAR(p.find("kp")->second.value(), std::stod(d.fields[6]), 1e-4) << d.comment;
        EXPECT_EQ(std::lround(p.find("mult")->second.value()), std::stol(d.fields[7])) << d.comment;
    }
    // Lennard-Jones with sigma and epsilon of the atom types, and those of
    // the pairs as computed by ACT
    auto defaults = section(top, "defaults");
    ASSERT_EQ(1u, defaults.size());
    const std::vector<std::string> expected = { "1", "3", "no", "1", "1" };
    EXPECT_EQ(expected, defaults[0].fields);
    for (const auto &a : section(top, "atomtypes"))
    {
        ASSERT_EQ(7u, a.fields.size()) << a.fields[0];
        EXPECT_TRUE(a.fields[4] == "A" || a.fields[4] == "S" || a.fields[4] == "V") << a.fields[0];
    }
    auto &vdw   = pd.findForcesConst(InteractionType::VDW);
    const auto pairs = section(top, "nonbond_params");
    EXPECT_FALSE(pairs.empty());
    for (const auto &np : pairs)
    {
        ASSERT_EQ(5u, np.fields.size());
        const auto vi = pd.findParticleType(np.fields[0])->interactionTypeToIdentifier(InteractionType::VDW);
        const auto vj = pd.findParticleType(np.fields[1])->interactionTypeToIdentifier(InteractionType::VDW);
        auto &p = vdw.findParameterMapConst(Identifier({ vi.id(), vj.id() }, { 1 }, CanSwap::Yes));
        EXPECT_NEAR(p.find("sigma")->second.internalValue(), std::stod(np.fields[3]), 1e-5) << np.fields[0] << " " << np.fields[1];
        EXPECT_NEAR(p.find("epsilon")->second.internalValue(), std::stod(np.fields[4]), 1e-5) << np.fields[0] << " " << np.fields[1];
    }
    // The coordinates of all atoms, with the names in the topology
    std::ifstream            pdb(files["-c"]);
    std::vector<std::string> names;
    for (std::string line; std::getline(pdb, line); )
    {
        if (line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0)
        {
            std::istringstream name(line.substr(12, 4));
            names.emplace_back("");
            name >> names.back();
        }
    }
    ASSERT_EQ(10u, names.size());
    for (size_t i = 0; i < atoms.size(); i++)
    {
        EXPECT_EQ(names[i], atoms[i].fields[4]);
    }
}

// Test: The shells are excluded from their core and have the same
// exclusions as their core
// Test: The particles with their charges and coordinates in JSON format,
// for programs that show the molecule and its charge distribution
TEST_F(GentopTest, Json)
{
    std::map<std::string, std::string> files = { { "-json", "acetone.json" }, { "-c", "acetone.pdb" } };
    run("forcefields/ACS-pg1.xml", "acetone.sdf", &files);
    std::ifstream in(files["-json"]);
    ASSERT_TRUE(in.good()) << "Cannot open " << files["-json"];
    const auto doc = nlohmann::json::parse(in, nullptr, false);
    ASSERT_TRUE(doc.is_object());
    ASSERT_TRUE(doc.contains("gentop"));
    const auto &mols = doc["gentop"]["molecules"];
    ASSERT_TRUE(mols.is_array());
    ASSERT_EQ(1U, mols.size());
    const auto &mol = mols[0];
    EXPECT_EQ("C3H6O", mol["molecule"].value("formula", ""));
    const auto &particles = mol["molecule"]["particles"];
    // Ten atoms, each with a shell
    ASSERT_EQ(20U, particles.size());
    EXPECT_EQ(9U, mol["molecule"]["bonds"].size());
    EXPECT_EQ("nm", mol["coordinates"].value("unit", ""));
    const auto &x = mol["coordinates"]["positions"];
    ASSERT_TRUE(x.is_array());
    ASSERT_EQ(particles.size(), x.size());
    double qtot = 0;
    for (size_t i = 0; i < particles.size(); i++)
    {
        ASSERT_TRUE(particles[i]["charge"].is_number());
        qtot += particles[i]["charge"].get<double>();
        ASSERT_EQ(3U, x[i].size());
        if (particles[i].value("particletype", "") == "Shell")
        {
            // A shell is close to its atom
            const int parent = particles[i]["parents"][0].get<int>();
            double    d2     = 0;
            for (int m = 0; m < 3; m++)
            {
                d2 += std::pow(x[i][m].get<double>() - x[parent][m].get<double>(), 2);
            }
            EXPECT_LT(std::sqrt(d2), 0.05) << "shell " << i;
        }
    }
    EXPECT_NEAR(0.0, qtot, 1e-4);
    // The charges are those of the topology
    std::map<std::string, std::string> topFiles = { { "-p", "acetone.top" } };
    run("forcefields/ACS-pg1.xml", "acetone.sdf", &topFiles);
    const auto atoms = section(readTopology(topFiles["-p"]), "atoms");
    ASSERT_EQ(particles.size(), atoms.size());
    for (size_t i = 0; i < atoms.size(); i++)
    {
        // The charge is the seventh field of a line in the atoms section
        EXPECT_NEAR(std::stod(atoms[i].fields[6]), particles[i]["charge"].get<double>(), 1e-5) << "particle " << i;
    }
}

// Test: Without the option no JSON file is written
TEST_F(GentopTest, NoJsonByDefault)
{
    std::remove("gentop.json");
    std::map<std::string, std::string> files = { { "-c", "acetone.pdb" } };
    run("forcefields/ACS-g1.xml", "acetone.sdf", &files);
    EXPECT_FALSE(std::ifstream("gentop.json").good());
}

TEST_F(GentopTest, GromacsTopologyShells)
{
    std::map<std::string, std::string> files = { { "-p", "acetone.top" } };
    run("forcefields/ACS-pg1.xml", "acetone.sdf", &files);
    const auto top = readTopology(files["-p"]);
    ASSERT_EQ(20u, section(top, "atoms").size());
    ASSERT_EQ(10u, section(top, "polarization").size());
    EXPECT_EQ(expectedExclusionPairs(top), exclusionPairs(top));
}

//! \return the total charge of each atom with its shells, in the order of the atoms
static std::vector<double> atomCharges(const Topology &top)
{
    const auto atoms = section(top, "atoms");
    std::vector<int> core(atoms.size() + 1);
    for (size_t i = 0; i <= atoms.size(); i++)
    {
        core[i] = i;
    }
    for (const auto &p : section(top, "polarization"))
    {
        core[std::stoi(p.fields[1])] = std::stoi(p.fields[0]);
    }
    std::map<int, double> q;
    for (const auto &a : atoms)
    {
        q[core[std::stoi(a.fields[0])]] += std::stod(a.fields[6]);
    }
    std::vector<double> charges;
    for (const auto &qq : q)
    {
        charges.push_back(qq.second);
    }
    return charges;
}

// Test: The charges of quantum chemistry from the input, with and without
// shells, which take their charge from their core
TEST_F(GentopTest, ChargesRead)
{
    MsgHandler           msghandler;
    std::vector<MolProp> mps;
    auto molpropFile = fileManager().getInputFilePath("molprop/molprop.xml");
    MolPropRead(&msghandler, molpropFile.c_str(), &mps);
    std::vector<double> mulliken;
    for (const auto &mp : mps)
    {
        if (mp.getMolname() == "water")
        {
            for (const auto &ca : mp.findExperimentConst(JobType::OPT)->calcAtomConst())
            {
                EXPECT_TRUE(ca.hasCharge("qMulliken"));
                mulliken.push_back(ca.charge("qMulliken"));
            }
        }
    }
    ASSERT_EQ(3u, mulliken.size());
    for (const char *ffName : { "forcefields/ACS-g1.xml", "forcefields/ACS-pg1.xml" })
    {
        std::map<std::string, std::string> files = { { "-p", "water.top" } };
        run(ffName, { { "-charges", molpropFile }, { "-db", "water" }, { "-qqm", "qMulliken" } }, &files);
        auto charges = atomCharges(readTopology(files["-p"]));
        ASSERT_EQ(mulliken.size(), charges.size()) << ffName;
        for (size_t i = 0; i < charges.size(); i++)
        {
            EXPECT_NEAR(mulliken[i], charges[i], 1e-5) << ffName << " atom " << i;
        }
    }
}

// Test: Custom charges, of which there must be one for each particle.
// The space keeps the command line parser from taking the charges for an option.
TEST_F(GentopTest, ChargesCustom)
{
    auto molpropFile = fileManager().getInputFilePath("molprop/molprop.xml");
    std::map<std::string, std::string> files = { { "-p", "water.top" } };
    run("forcefields/ACS-g1.xml", { { "-charges", molpropFile }, { "-db", "water" }, { "-qcustom", " -0.8 0.4 0.4" } }, &files);
    std::vector<double> expected = { -0.8, 0.4, 0.4 };
    auto charges = atomCharges(readTopology(files["-p"]));
    ASSERT_EQ(expected.size(), charges.size());
    for (size_t i = 0; i < charges.size(); i++)
    {
        EXPECT_NEAR(expected[i], charges[i], 1e-6) << "atom " << i;
    }
    files = { { "-p", "water2.top" } };
    EXPECT_NE(0, status("forcefields/ACS-g1.xml", { { "-charges", molpropFile }, { "-db", "water" }, { "-qcustom", " -0.8 0.8" } }, &files));
}

// Test: The charges of the particle types in the force field
TEST_F(GentopTest, ChargesNone)
{
    std::map<std::string, std::string> files = { { "-p", "water.top" } };
    run("forcefields/SPCE.xml",
        { { "-f", fileManager().getInputFilePath("molecules/water_dimer.pdb") }, { "-qalg", "None" } }, &files);
    std::vector<double> expected = { -0.8476, 0.4238, 0.4238, -0.8476, 0.4238, 0.4238 };
    auto charges = atomCharges(readTopology(files["-p"]));
    ASSERT_EQ(expected.size(), charges.size());
    for (size_t i = 0; i < charges.size(); i++)
    {
        EXPECT_NEAR(expected[i], charges[i], 1e-6) << "atom " << i;
    }
}

// Test: Morse bonds, Urey-Bradley angles and Fourier dihedrals, which are
// the Ryckaert-Bellemans dihedrals of GROMACS
TEST_F(GentopTest, OtherFunctions)
{
    std::map<std::string, std::string> files = { { "-p", "acetone.top" } };
    run("forcefields/ACS-g-dih.xml", "acetone.sdf", &files);
    const auto top = readTopology(files["-p"]);
    checkBondeds(top);
    for (const auto &b : section(top, "bonds"))
    {
        EXPECT_EQ("3", b.fields[2]);
    }
    for (const auto &a : section(top, "angles"))
    {
        EXPECT_EQ("5", a.fields[3]);
    }
    ForceField pd;
    readForceField(fileManager().getInputFilePath("forcefields/ACS-g-dih.xml"), &pd);
    int nfourier = 0;
    for (const auto &d : section(top, "dihedrals"))
    {
        if (std::stoi(d.fields[4]) != 3)
        {
            continue;
        }
        nfourier++;
        auto &p = parameters(pd, InteractionType::PROPER_DIHEDRALS, d.comment);
        for (int k = 0; k < 6; k++)
        {
            const auto c = p.find(gmx::formatString("c%d", k));
            ASSERT_TRUE(c != p.end());
            EXPECT_NEAR(c->second.value(), std::stod(d.fields[5+k]), 1e-4) << d.comment << " c" << k;
        }
    }
    EXPECT_GT(nfourier, 0);
}

// Test: An include topology has only the molecule
TEST_F(GentopTest, IncludeTopology)
{
    std::map<std::string, std::string> files = { { "-oi", "acetone.itp" } };
    run("forcefields/ACS-g1.xml", "acetone.sdf", &files);
    const auto top = readTopology(files["-oi"]);
    EXPECT_TRUE(section(top, "defaults").empty());
    EXPECT_TRUE(section(top, "atomtypes").empty());
    EXPECT_EQ(10u, section(top, "atoms").size());
    checkBondeds(top);
}

// Test: The potential and the density of the charges of water on a grid,
// the potential far from the atoms is that of point charges and the
// density integrates to the total charge. Reading the potential as the
// reference gives no difference.
TEST_F(GentopTest, Cubes)
{
    std::vector<std::pair<std::string, std::string>> water = {
        { "-charges", fileManager().getInputFilePath("molprop/molprop.xml") },
        { "-db", "water" }, { "-spacing", "0.01" }
    };
    std::map<std::string, std::string> files = {
        { "-c", "water.pdb" }, { "-pot", "pot.cub" }, { "-rho", "rho.cub" },
        { "-his", "his.xvg" }, { "-pc", "pc.xvg" }, { "-pdbdiff", "diff.pdb" }
    };
    run("forcefields/ACS-g1.xml", water, &files);
    auto pot = readCube(files["-pot"]);
    ASSERT_EQ(3u, pot.atoms.size());
    double qtot = 0;
    for (size_t i = 0; i < 3; i++)
    {
        EXPECT_EQ(i == 1 ? 8 : 1, pot.atoms[i][0]);
        qtot += pot.atoms[i][1];
    }
    EXPECT_NEAR(0, qtot, 1e-5);
    ASSERT_EQ(static_cast<size_t>(pot.n[0]*pot.n[1]*pot.n[2]), pot.values.size());
    // The first point is the origin, at least 0.2 nm from the atoms
    double V = 0;
    for (const auto &a : pot.atoms)
    {
        V += a[1]/std::hypot(a[2] - pot.origin[0], a[3] - pot.origin[1], a[4] - pot.origin[2]);
    }
    EXPECT_NEAR(V, pot.values[0], 0.01*std::abs(V));

    const auto rho = readCube(files["-rho"]);
    ASSERT_EQ(pot.values.size(), rho.values.size());
    double sum = 0, sumAbs = 0;
    for (auto r : rho.values)
    {
        sum    += r;
        sumAbs += std::abs(r);
    }
    double dV = rho.space[0]*rho.space[1]*rho.space[2];
    EXPECT_GT(sumAbs*dV, 0.1);
    EXPECT_NEAR(0, sum*dV, 1e-2*sumAbs*dV);

    for (const auto &f : { "-his", "-pc" })
    {
        EXPECT_GT(std::filesystem::file_size(files[f]), 0u) << f;
    }
    auto pdbdiff = std::filesystem::path(files["-pdbdiff"]);
    pdbdiff.replace_filename("0_" + pdbdiff.filename().string());
    EXPECT_GT(std::filesystem::file_size(pdbdiff), 0u);
    std::filesystem::remove(pdbdiff);

    water.emplace_back("-ref", files["-pot"]);
    std::map<std::string, std::string> diffFiles = {
        { "-c", "water2.pdb" }, { "-diff", "diff.cub" }, { "-diffhist", "diffhist.xvg" }
    };
    run("forcefields/ACS-g1.xml", water, &diffFiles);
    auto diff = readCube(diffFiles["-diff"]);
    ASSERT_EQ(pot.values.size(), diff.values.size());
    // The differences are due to the six digits of the values in the file
    // and of the origin and spacing of the grid
    double maxPot = 0, maxDiff = 0;
    for (size_t i = 0; i < diff.values.size(); i++)
    {
        maxPot  = std::max(maxPot, std::abs(pot.values[i]));
        maxDiff = std::max(maxDiff, std::abs(diff.values[i]));
    }
    EXPECT_LT(maxDiff, 1e-5*maxPot);
    const auto corr = std::filesystem::path(diffFiles["-diffhist"]).replace_extension().string() + "-corr.xvg";
    EXPECT_GT(std::filesystem::file_size(corr), 0u);
    std::filesystem::remove(corr);
}

//! Runs gentop with a force field from which parameters are removed
/*! \brief Run gentop and compare what it writes with the reference data
 *
 * These tests fix the results of the tool with the given options, such
 * that a change in how the options reach the work shows.
 */
class GentopResultsTest : public GentopTest
{
protected:
    /*! \brief Run gentop and check the files it writes
     * \param[in] ffName Force field file, relative to the test data
     * \param[in] input  Options that are passed as they are
     * \param[in] files  Output files by option, among -json, -p, -oi, -c and
     *                   -openmm; the names are made temporary
     */
    void check(const std::string                                      &ffName,
               const std::vector<std::pair<std::string, std::string>> &input,
               std::map<std::string, std::string>                      files)
    {
        run(ffName, input, &files);
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        if (files.count("-json") > 0)
        {
            const auto tree = readJsonOutput(files["-json"]);
            ASSERT_TRUE(tree.is_object() && tree.contains("gentop")) << files["-json"];
            // The charges are generated, in single precision they differ in
            // the sixth decimal
            checker.setDefaultTolerance(buildTolerance(2e-5, 1e-8, 1e-5));
            checkJson(&checker, tree["gentop"], "gentop");
        }
        // Topologies with charges in five decimals, structures in Angstrom
        // with three, and OpenMM files with all digits
        const std::map<std::string, gmx::test::FloatingPointTolerance> tolerances = {
            { "-p", buildTolerance(3e-5, 2e-5, 1e-5) },
            { "-oi", buildTolerance(3e-5, 2e-5, 1e-5) },
            { "-c", buildTolerance(2e-3, 2e-3, 0) },
            { "-openmm", buildTolerance(1e-5, 1e-8, 1e-5) }
        };
        for (const auto &t : tolerances)
        {
            if (files.count(t.first) > 0)
            {
                checker.setDefaultTolerance(t.second);
                checkTextFile(&checker, files[t.first], {}, t.first.c_str());
            }
        }
    }
};

// The results of the tool for several sets of options, as reference data

TEST_F(GentopResultsTest, Defaults)
{
    check("forcefields/ACS-g1.xml",
          { { "-f", fileManager().getInputFilePath("molecules/acetone.sdf") } },
          { { "-json", "defaults.json" }, { "-p", "defaults.top" }, { "-c", "defaults.pdb" } });
}

TEST_F(GentopResultsTest, Polarizable)
{
    check("forcefields/ACS-pg1.xml",
          { { "-f", fileManager().getInputFilePath("molecules/acetone.sdf") } },
          { { "-json", "pol.json" }, { "-p", "pol.top" } });
}

TEST_F(GentopResultsTest, OpenMM)
{
    // With -openmm no GROMACS topology is written
    check("forcefields/ACS-pg1.xml",
          { { "-f", fileManager().getInputFilePath("molecules/acetone.sdf") } },
          { { "-openmm", "openmm.xml" }, { "-c", "openmm.pdb" } });
}

TEST_F(GentopResultsTest, OpenMMOptions)
{
    // The structure with the shells, a heavier shell, all parameters also
    // those that were not trained, and the types of the force field
    check("forcefields/ACS-pg1.xml",
          { { "-f", fileManager().getInputFilePath("molecules/acetone.sdf") },
            { "-mDrude", "0.4" }, { "-ntrain", "0" }, { "-numberAtypes", "no" }, { "-ws", "yes" } },
          { { "-openmm", "options-openmm.xml" }, { "-c", "options.pdb" } });
}

TEST_F(GentopResultsTest, Dimer)
{
    check("forcefields/ACS-g1.xml",
          { { "-f", fileManager().getInputFilePath("molecules/methanol#water.pdb") } },
          { { "-json", "dimer.json" }, { "-p", "dimer.top" }, { "-c", "dimer.pdb" } });
}

TEST_F(GentopResultsTest, DimerOpenMM)
{
    // A residue for each compound, in both files
    check("forcefields/ACS-g1.xml",
          { { "-f", fileManager().getInputFilePath("molecules/methanol#water.pdb") } },
          { { "-openmm", "dimer-openmm.xml" }, { "-c", "dimer-openmm.pdb" } });
}

TEST_F(GentopResultsTest, IncludeTopology)
{
    check("forcefields/ACS-g1.xml",
          { { "-f", fileManager().getInputFilePath("molecules/acetone.sdf") } },
          { { "-oi", "include.itp" } });
}

TEST_F(GentopResultsTest, ChargesOfTypes)
{
    check("forcefields/ACS-g1.xml",
          { { "-f", fileManager().getInputFilePath("molecules/acetone.sdf") }, { "-qalg", "None" } },
          { { "-json", "qnone.json" } });
}

TEST_F(GentopResultsTest, ChargesRead)
{
    check("forcefields/ACS-g1.xml",
          { { "-charges", fileManager().getInputFilePath("molprop/molprop.xml") }, { "-db", "water" },
            { "-qqm", "qMulliken" } },
          { { "-json", "qread.json" } });
}

TEST_F(GentopResultsTest, ChargesCustom)
{
    check("forcefields/ACS-g1.xml",
          // A leading space, such that the value is not taken for an option
          { { "-f", fileManager().getInputFilePath("molecules/water.pdb") }, { "-qcustom", " -0.8 0.4 0.4" } },
          { { "-json", "qcustom.json" } });
}

TEST_F(GentopResultsTest, TotalCharge)
{
    check("forcefields/ACS-g1.xml",
          { { "-f", fileManager().getInputFilePath("molecules/sodium#water.pdb") }, { "-qtot", "1" } },
          { { "-json", "qtot.json" } });
}

//! \return the contents of a file
static std::string fileContents(const std::string &fileName)
{
    const std::ifstream in(fileName);
    std::stringstream text;
    text << in.rdbuf();
    return text.str();
}

// Test: A topology is written only when it is asked for, also with an OpenMM file
TEST_F(GentopTest, TopologyOnlyWhenAskedFor)
{
    const auto acetone = fileManager().getInputFilePath("molecules/acetone.sdf");
    const auto dir     = std::filesystem::path(fileManager().getTemporaryFilePath("x")).parent_path();
    const auto old     = std::filesystem::current_path();
    // In a directory of its own, where the default name would be used
    const auto work    = dir / "topology-only-when-asked";
    std::filesystem::create_directories(work);
    std::filesystem::current_path(work);
    {
        std::map<std::string, std::string> files = { { "-c", "structure.pdb" } };
        run("forcefields/ACS-g1.xml", { { "-f", acetone } }, &files);
        EXPECT_TRUE(std::filesystem::exists(files["-c"]));
        EXPECT_TRUE(std::filesystem::is_empty(work));
    }
    {
        std::map<std::string, std::string> files = { { "-openmm", "both.xml" }, { "-p", "both.top" } };
        run("forcefields/ACS-g1.xml", { { "-f", acetone } }, &files);
        EXPECT_TRUE(std::filesystem::exists(files["-openmm"]));
        EXPECT_TRUE(std::filesystem::exists(files["-p"]));
        EXPECT_FALSE(section(readTopology(files["-p"]), "atoms").empty());
    }
    std::filesystem::current_path(old);
}

// The topologies without a command line, as a program that uses ACT makes them

TEST_F(GentopTest, WithoutCommandLineSameAsTheProgram)
{
    const auto acetone = fileManager().getInputFilePath("molecules/acetone.sdf");
    // A topology, a structure and the compounds in JSON
    {
        std::map<std::string, std::string> files = { { "-p", "tool.top" }, { "-c", "tool.pdb" },
                                                     { "-json", "tool.json" } };
        run("forcefields/ACS-g1.xml", { { "-f", acetone } }, &files);
        GentopOptions options;
        options.forceFieldFile = fileManager().getInputFilePath("forcefields/ACS-g1.xml");
        options.compounds.setFileName(acetone);
        options.topologyFile  = fileManager().getTemporaryFilePath("direct.top");
        options.structureFile = fileManager().getTemporaryFilePath("direct.pdb");
        options.jsonFile      = fileManager().getTemporaryFilePath("direct.json");
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        EXPECT_EQ(0, alexandria::runGentop(&msghandler, options));
        EXPECT_EQ(fileContents(files["-p"]), fileContents(options.topologyFile));
        EXPECT_EQ(fileContents(files["-c"]), fileContents(options.structureFile));
        EXPECT_EQ(fileContents(files["-json"]), fileContents(options.jsonFile));
    }
    // OpenMM with other options, and a structure with the shells in another box
    {
        std::map<std::string, std::string> files = { { "-openmm", "tool-openmm.xml" }, { "-c", "tool-openmm.pdb" } };
        run("forcefields/ACS-pg1.xml", { { "-f", acetone }, { "-mDrude", "0.4" }, { "-ntrain", "0" },
                                         { "-numberAtypes", "no" }, { "-ws", "yes" } }, &files);
        GentopOptions options;
        options.forceFieldFile = fileManager().getInputFilePath("forcefields/ACS-pg1.xml");
        options.compounds.setFileName(acetone);
        options.openmmFile    = fileManager().getTemporaryFilePath("direct-openmm.xml");
        options.structureFile = fileManager().getTemporaryFilePath("direct-openmm.pdb");
        // As the program reads it, a real
        options.drudeMass     = static_cast<real>(0.4);
        options.minTraining   = 0;
        options.numberTypes   = false;
        options.writeShells   = true;
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        EXPECT_EQ(0, alexandria::runGentop(&msghandler, options));
        EXPECT_EQ(fileContents(files["-openmm"]), fileContents(options.openmmFile));
        EXPECT_EQ(fileContents(files["-c"]), fileContents(options.structureFile));
        // The box of the structure file
        options.box = std::array<double, 3>{ 3, 4, 5 };
        EXPECT_EQ(0, alexandria::runGentop(&msghandler, options));
        EXPECT_NE(std::string::npos, fileContents(options.structureFile).find("CRYST1   30.000   40.000   50.000"));
    }
}

TEST_F(GentopTest, WithoutCommandLineWrongOptions)
{
    const auto valid = [this]()
    {
        GentopOptions options;
        options.forceFieldFile = fileManager().getInputFilePath("forcefields/ACS-g1.xml");
        options.compounds.setFileName(fileManager().getInputFilePath("molecules/acetone.sdf"));
        return options;
    };
    std::vector<std::pair<GentopOptions, std::string>> wrong;
    {
        auto o = valid();
        o.forceFieldFile.clear();
        wrong.emplace_back(o, "-ff");
    }
    {
        GentopOptions o;
        o.forceFieldFile = fileManager().getInputFilePath("forcefields/ACS-g1.xml");
        wrong.emplace_back(o, "-f");
    }
    {
        auto o = valid();
        o.drudeMass = -0.1;
        wrong.emplace_back(o, "-mDrude");
    }
    {
        auto o = valid();
        o.minTraining = -1;
        wrong.emplace_back(o, "-ntrain");
    }
    {
        auto o = valid();
        o.spacing = 0;
        wrong.emplace_back(o, "-spacing");
    }
    {
        auto o = valid();
        o.box = std::array<double, 3>{ 3, -1, 3 };
        wrong.emplace_back(o, "-box");
    }
    for (const auto &w : wrong)
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        try
        {
            alexandria::runGentop(&msghandler, w.first);
            ADD_FAILURE() << "Accepted options with a wrong " << w.second;
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find(w.second)) << ex.what();
        }
    }
    // The program says the same and stops
    std::map<std::string, std::string> files = { { "-p", "wrong.top" } };
    EXPECT_EQ(1, status("forcefields/ACS-g1.xml",
                        { { "-f", fileManager().getInputFilePath("molecules/acetone.sdf") }, { "-ntrain", "-2" } },
                        &files));
    EXPECT_FALSE(std::filesystem::exists(files["-p"]));
}

class GentopMissingTest : public GentopTest
{
protected:
    /*! \brief Run gentop for acetone with a force field from which a
     * torsion or one of its parameters is removed
     * \param[in]  id        Identifier of the torsion
     * \param[in]  parameter The parameter to remove, all if empty
     * \param[out] status    The exit status of gentop
     * \return the number of messages about missing parameters of each
     *         torsion, and of the other ones as "other"
     */
    std::map<std::string, int> missing(const std::string &id, const std::string &parameter,
                                       int *status)
    {
        ForceField pd;
        readForceField(fileManager().getInputFilePath("forcefields/ACS-pg1.xml"), &pd);
        auto propers = pd.findForces(InteractionType::PROPER_DIHEDRALS)->parameters();
        for (auto p = propers->begin(); p != propers->end(); ++p)
        {
            if (p->first.id() == id)
            {
                if (parameter.empty())
                {
                    propers->erase(p);
                }
                else
                {
                    p->second.erase(parameter);
                }
                break;
            }
        }
        const auto ffName = fileManager().getTemporaryFilePath("missing.xml");
        writeForceField(ffName, &pd, false);

        const auto logName = fileManager().getTemporaryFilePath("missing.log");
        gmx::test::CommandLine cmd;
        cmd.append("gentop");
        cmd.addOption("-ff", ffName);
        cmd.addOption("-f", fileManager().getInputFilePath("molecules/acetone.sdf"));
        cmd.addOption("-c", fileManager().getTemporaryFilePath("missing.pdb"));
        cmd.addOption("-g", logName);
        *status = gentop(cmd.argc(), cmd.argv());
        std::ifstream              in(logName);
        std::map<std::string, int> torsions;
        for (std::string line; std::getline(in, line); )
        {
            const std::string key("Missing parameter in force field: FOURIER_DIHEDRALS ");
            auto        pos = line.find(key);
            if (pos != std::string::npos)
            {
                std::istringstream iss(line.substr(pos + key.size()));
                std::string        name;
                iss >> name;
                torsions[name] += 1;
            }
            else if (line.find("Missing parameter") != std::string::npos)
            {
                torsions["other"] += 1;
            }
        }
        return torsions;
    }
};

// Test: Each of the six torsions without parameters is reported, and not
// the torsions after them, which have all their parameters, gentop fails
TEST_F(GentopMissingTest, Torsion)
{
    int status;
    EXPECT_EQ((std::map<std::string, int>{ { "h_b~c3_b~c2_b=o2_b", 6 } }),
              missing("h_b~c3_b~c2_b=o2_b", "", &status));
    EXPECT_NE(0, status);
}

// Test: A parameter that is missing from a torsion is zero, with a
// warning about that torsion only
TEST_F(GentopMissingTest, Parameter)
{
    int  status;
    auto torsions = missing("h_b~c3_b~c2_b=o2_b", "c5", &status);
    ASSERT_EQ(1u, torsions.size());
    EXPECT_EQ("h_b~c3_b~c2_b=o2_b", torsions.begin()->first);
    EXPECT_EQ(0, status);
}


/* Test: Without bond charge corrections for all bonds of water, SQE gives
 * no charges, and gentop, simulate and b2 leave water out and fail.
 */
TEST_F(GentopMissingTest, BondCorrection)
{
    ForceField pd;
    readForceField(fileManager().getInputFilePath("forcefields/ACS-g1.xml"), &pd);
    ASSERT_EQ(ChargeGenerationAlgorithm::SQE, pd.chargeGenerationAlgorithm());
    auto   *bcc     = pd.findForces(InteractionType::BONDCORRECTIONS)->parameters();
    size_t  nerased = 0;
    for (auto b = bcc->begin(); b != bcc->end(); )
    {
        if (b->first.id() == "h_z~o3_z" || b->first.id() == "o3_z~h_z")
        {
            b = bcc->erase(b);
            nerased++;
        }
        else
        {
            ++b;
        }
    }
    ASSERT_GT(nerased, 0u);
    const auto ffName = fileManager().getTemporaryFilePath("nobcc.xml");
    writeForceField(ffName, &pd, false);

    std::vector<std::pair<std::string, int (*)(int, char *[])>> programs = {
        { "gentop", gentop }, { "simulate", simulate }, { "b2", b2 }
    };
    for (const auto &program : programs)
    {
        const auto logName = fileManager().getTemporaryFilePath(program.first + ".log");
        gmx::test::CommandLine cmd;
        cmd.append(program.first.c_str());
        cmd.addOption("-ff", ffName);
        cmd.addOption("-charges", fileManager().getInputFilePath("molprop/molprop.xml"));
        cmd.addOption("-db", "water");
        cmd.addOption("-g", logName);
        if (program.first == "gentop")
        {
            cmd.addOption("-p", fileManager().getTemporaryFilePath("water.top"));
        }
        else if (program.first == "simulate")
        {
            cmd.append("-sp");
            cmd.addOption("-c", fileManager().getTemporaryFilePath("water.pdb"));
        }
        EXPECT_NE(0, program.second(cmd.argc(), cmd.argv())) << program.first;
        const std::ifstream in(logName);
        std::stringstream log;
        log << in.rdbuf();
        EXPECT_NE(std::string::npos, log.str().find("Parameters for charge generation missing: water")) << program.first;
    }
}
} // namespace

} // namespace alexandria
