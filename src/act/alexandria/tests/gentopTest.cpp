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

#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/alex_modules.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"

#include "testutils/cmdlinetest.h"
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
        std::string content   = line.substr(0, semicolon);
        std::string comment   = semicolon == std::string::npos ? "" : line.substr(semicolon + 1);
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
        gmx::test::CommandLine cmd;
        cmd.append("gentop");
        cmd.addOption("-ff", fileManager().getInputFilePath(ffName));
        cmd.addOption("-f", fileManager().getInputFilePath("mols/" + molecule));
        cmd.addOption("-g", fileManager().getTemporaryFilePath("gentop.log"));
        for (auto &f : *files)
        {
            f.second = fileManager().getTemporaryFilePath(f.second);
            cmd.addOption(f.first.c_str(), f.second);
        }
        ASSERT_EQ(0, gentop(cmd.argc(), cmd.argv()));
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
                auto np   = nparams.find({ s.first, type });
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
    run("ACS-g.xml", "acetone.sdf", &files);
    auto top = readTopology(files["-p"]);
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
    ForceField pd;
    readForceField(fileManager().getInputFilePath("ACS-g.xml"), &pd);
    auto angles = section(top, "angles");
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
    std::vector<std::string> expected = { "1", "3", "no", "1", "1" };
    EXPECT_EQ(expected, defaults[0].fields);
    for (const auto &a : section(top, "atomtypes"))
    {
        ASSERT_EQ(7u, a.fields.size()) << a.fields[0];
        EXPECT_TRUE(a.fields[4] == "A" || a.fields[4] == "S" || a.fields[4] == "V") << a.fields[0];
    }
    auto &vdw   = pd.findForcesConst(InteractionType::VDW);
    auto  pairs = section(top, "nonbond_params");
    EXPECT_FALSE(pairs.empty());
    for (const auto &np : pairs)
    {
        ASSERT_EQ(5u, np.fields.size());
        auto vi = pd.findParticleType(np.fields[0])->interactionTypeToIdentifier(InteractionType::VDW);
        auto vj = pd.findParticleType(np.fields[1])->interactionTypeToIdentifier(InteractionType::VDW);
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
            names.push_back("");
            name >> names.back();
        }
    }
    ASSERT_EQ(10u, names.size());
    for (size_t i = 0; i < atoms.size(); i++)
    {
        EXPECT_EQ(names[i], atoms[i].fields[4]);
    }
}

// Test: Morse bonds, Urey-Bradley angles and Fourier dihedrals, which are
// the Ryckaert-Bellemans dihedrals of GROMACS
TEST_F(GentopTest, OtherFunctions)
{
    std::map<std::string, std::string> files = { { "-p", "acetone.top" } };
    run("ACS-g-dih.xml", "acetone.sdf", &files);
    auto top = readTopology(files["-p"]);
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
    readForceField(fileManager().getInputFilePath("ACS-g-dih.xml"), &pd);
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
            auto c = p.find(gmx::formatString("c%d", k));
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
    run("ACS-g.xml", "acetone.sdf", &files);
    auto top = readTopology(files["-oi"]);
    EXPECT_TRUE(section(top, "defaults").empty());
    EXPECT_TRUE(section(top, "atomtypes").empty());
    EXPECT_EQ(10u, section(top, "atoms").size());
    checkBondeds(top);
}

} // namespace

} // namespace alexandria
