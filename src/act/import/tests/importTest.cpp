/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2019-2026
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cmath>

#include <fstream>
#include <map>
#include <sstream>

#include <gtest/gtest.h>

#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/import/import.h"
#include "act/import/import_utils.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class ImportTest : public ::testing::TestWithParam<std::tuple<std::tuple<std::string, double>, bool, bool > >
{
private:
    ForceField *pd;
    std::string molname;
    double      qtot;
    bool        bondTest;
    bool        oneH;
    gmx::test::TestReferenceData    refData_;
    gmx::test::TestReferenceChecker checker_;

public:   
    ImportTest()
    {
        pd       = getForceField("gaff_hydrogen_types");
        molname  = std::get<0>(std::get<0>(GetParam()));
        qtot     = std::get<1>(std::get<0>(GetParam()));
        oneH     = std::get<1>(GetParam());
        bondTest = std::get<2>(GetParam());
    }
    
 protected:
    void runTest()
    {
        gmx::test::TestReferenceChecker  checker_(refData_.rootChecker());
        const char                      *conf     = (char *)"minimum";
        std::string                      dataName;
        std::vector<alexandria::MolProp> molprops;

        std::string dirmolname("molecules/");
        dirmolname     += molname;
        dataName        = gmx::test::TestFileManager::getInputFilePath(dirmolname);
        bool       userqtot = false;
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        importFile(&msghandler, pd, dataName, &molprops,
                   conf, JobType::OPT, userqtot, &qtot, oneH);
        EXPECT_TRUE(msghandler.ok());
        for(auto &molprop: molprops)
        {
            checker_.checkString(molname, "molname");
            checker_.checkDouble(qtot, "qtot");
            EXPECT_TRUE(!molprop.experimentConst().empty());
            if (molprop.experimentConst().empty())
            {
                continue;
            }
            const auto               exper = molprop.experimentConst().begin();
            std::vector<std::string> atypes;
            for (auto &ca : exper->calcAtomConst())
            {
                atypes.push_back(ca.getObtype());
            }
            if (!bondTest)
            {
                checker_.checkInteger(static_cast<int>(atypes.size()), molname.c_str());
                checker_.checkSequence(atypes.begin(), atypes.end(), "atomtypes");
                // Check for fragment charges in case of dimers
                const auto& fragments = molprop.fragments();
                if (fragments.size() == 2)
                {
                    for(int i = 0; i < 2; ++i)
                    {
                        const std::string qi = gmx::formatString("fragment_charge_%d", i);
                        checker_.checkInteger(fragments[i].charge(), qi.c_str());
                    }
                }
            }
            else
            {
                std::vector<std::string> bondorder;
                for (auto &bond : molprop.bondsConst())
                {
                    auto ai = bond.aI();
                    auto aj = bond.aJ();
                    bondorder.push_back(gmx::formatString("%s-%d %s-%d: %g",
                                                          atypes[ai].c_str(), ai,
                                                          atypes[aj].c_str(), aj,
                                                          bond.bondOrder()));
                }
                checker_.checkInteger(static_cast<int>(bondorder.size()), molname.c_str());
                checker_.checkSequence(bondorder.begin(), bondorder.end(), "bondorder");
            }
        }
    }
};

std::vector<std::tuple<std::string, double>> get_files()
{
    return {
        { "1--ethoxyethylphosphoryl-oxyethane-3-oep.log.pdb", 0 },
        { "1-1-dimethylguanidinium.sdf", 1 },
        { "1-amino-1-hydroxyguanidine-3-oep.log.pdb", 0 },
        { "1-buten-3-yne.sdf", 0 },
        { "1-methylsulfinylethene.sdf", 0 },
        { "1-nitroethyne.pdb", 0 },
        { "12-dimethyl-imidazole.sdf", 0 },
        { "123-trimethyl-imidazolium.sdf", 0 },
        { "125-thiadiazole.sdf", 0 },
        { "3H-diphosphole.sdf", 0 },
        { "N--ethenyl-N-hydroxyethanimidamide.sdf", 0 },
        { "N-methylmethanimine.sdf", 0 },
        { "acetaldehyde.sdf", 0 },
        { "acetate-3-oep.log.pdb", -1 },
        { "acetic-acid-3-oep.log.pdb", 0 },
        { "acetone-3-oep.log.pdb", 0 },
        { "acetone-nonplanar.pdb", 0 },
        { "acetone.sdf", 0 },
        { "acetonitrile.sdf", 0 },
        { "ammonia.sdf", 0 },
        { "aniline.sdf", 0 },
        { "carbon-dioxide.sdf", 0 },
        { "cyclobutene.sdf", 0 },
        { "cyclopropene.sdf", 0 },
        { "diethyl-sulfate-3-oep.log.pdb", 0 },
        { "dihydrogen-phosphate.sdf", -1 },
        { "dihydrogen-phosphate.xyz", -1 },
        { "dimethyl-carbonate.pdb", 0 },
        { "dimethyl-carbonate.xyz", 0 },
        { "dimethyl-sulfide.sdf", 0 },
        { "dimethylether.sdf", 0 },
        { "disulfur-monoxide-3-oep.log.pdb", 0 },
        { "ethane-12-diamine.sdf", 0 },
        { "ethyl-sulfate-3-oep.log.pdb", -1 },
        { "ethylsulfonylformaldehyde.sdf", 0 },
        { "fluorane.sdf", 0 },
        { "fluoranthene.sdf", 0 },
        { "formamide.sdf", 0 },
        { "formylphosphonic-acid.sdf", 0 },
        { "furan.sdf", 0 },
        { "glutamate-3-oep.log.pdb", -1 },
        { "glutamic-acid-3-oep.log.pdb", 0 },
        { "guanidine.pdb", 0 },
        { "histidine-hdhe-3-oep.log.pdb", 1 },
        { "hydrogen-bromide#hydrogen-fluoride.pdb", 0 },
        { "hydrogen-chloride.sdf", 0 },
        { "hydrogen-sulfide.sdf", 0 },
        { "isoxazole.sdf", 0 },
        { "isoxazole.xyz", 0 },
        { "merged.pdb", 0 },
        { "methanethiol.sdf", 0 },
        { "methyl-acetate.sdf", 0 },
        { "nitromethane-3-oep.log.pdb", 0 },
        { "nitromethane.xyz", 0 },
        { "oxophosphane.sdf", 0 },
        { "pentachlorophosphorane.sdf", 0 },
        { "pentafluorophosphorane.sdf", 0 },
        { "phosphine.sdf", 0 },
        { "phosphinine.sdf", 0 },
        { "phosphoethanoamine.sdf", 0 },
        { "phosphorus-nitride.sdf", 0 },
        { "propa-1-2-dienylidenephosphane.sdf", 0 },
        { "pyridine.sdf", 0 },
        { "sulfate.sdf", -2 },
        { "sulfate.xyz", -2 },
        { "sulfur-dioxide-3-oep.log.pdb", 0 },
        { "sulfur-monoxide.sdf", 0 },
        { "thiazirene.sdf", 0 },
        { "thiazyl-fluoride-3-oep.log.pdb", 0 },
        { "thiophene.sdf", 0 },
        { "uracil.sdf", 0 },
        { "water-3-oep.log.pdb", 0 },
        { "water_dimer.pdb", 0 },
        { "chloroform.sdf", 0 },
        { "dichloromethane.sdf", 0 },
        { "1,3,5-triazine.sdf", 0 },
        { "trimethylorthoformate.sdf", 0 },
        { "chloroform.xyz", 0 },
        { "dichloromethane.xyz", 0 },
        { "1,3,5-triazine.xyz", 0 },
        { "trimethylorthoformate.xyz", 0 },
        { "KBr.sdf", 0 },
        { "NaCl.sdf", 0 },
        { "ammonium#formate.pdb", 0 },
        { "1-methyl-adenine-conect.pdb", 0 },
        { "1-methyl-adenine.pdb", 0 },
        { "phosphor-hexafluoride.sdf", -1 },
        { "boron-tetrafluoride.sdf", -1 },
        { "helium-4.sdf", 0 },
        { "argon-krypton.sdf", 0 },
        { "phosphor-hexafluoride.pdb", -1 },
        { "boron-tetrafluoride.pdb", -1 },
        { "3-bromophenol.sdf", 0 },
        { "3-bromophenol-sodium.sdf", 1 }
    };
}

/*! \return booleans controlling tests */
static std::vector<bool> get_ab()
{
    return { false, true };
}

/*! \return booleans controlling tests */
static std::vector<bool> get_ba()
{
    return { true, false };
}

TEST_P (ImportTest, All)
{
    runTest();
}

//! Initiate file structure
auto my_files = get_files();
//! Initiate booleans for testing
auto my_btest = get_ab();
//! Initiate booleans for testing
auto my_oneH  = get_ba();

INSTANTIATE_TEST_CASE_P(IT, ImportTest, ::testing::Combine(::testing::ValuesIn(my_files), ::testing::ValuesIn(my_oneH),  ::testing::ValuesIn(my_btest)));

//! \return the atom types of the first experiment of each compound in a file
static std::vector<std::string> importTypes(const ForceField *pd, const std::string &fileName)
{
    std::vector<MolProp> molprops;
    MsgHandler           msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    double qtot = 0;
    importFile(&msghandler, pd, fileName, &molprops, "minimum",
               JobType::OPT, false, &qtot, false);
    EXPECT_TRUE(msghandler.ok());
    std::vector<std::string> atypes;
    for (const auto &mp : molprops)
    {
        for (const auto &ca : mp.experimentConst()[0].calcAtomConst())
        {
            atypes.push_back(ca.getObtype());
        }
    }
    return atypes;
}

// Test: In a box with many copies of a compound, a pattern can match more
// than 1000 times, the default maximum of RDKit. Each copy gets the atom
// types of the single compound.
TEST(ImportBoxTest, ManyMatches)
{
    auto pd      = getForceField("gaff_hydrogen_types");
    const auto single  = gmx::test::TestFileManager::getInputFilePath("molecules/acetone-3-oep.log.pdb");
    auto ref     = importTypes(pd, single);
    ASSERT_EQ(10u, ref.size());
    // Six hydrogens of the methyl groups in each copy
    const int ncopy = 180;
    std::vector<std::string> atoms, conect;
    {
        std::ifstream in(single);
        for (std::string line; std::getline(in, line); )
        {
            if (line.compare(0, 6, "HETATM") == 0)
            {
                atoms.push_back(line);
            }
            else if (line.compare(0, 6, "CONECT") == 0)
            {
                conect.push_back(line);
            }
        }
    }
    gmx::test::TestFileManager fileManager;
    const auto boxName = fileManager.getTemporaryFilePath("acetone-box.pdb");
    {
        std::ofstream out(boxName);
        int natom = atoms.size();
        for (int c = 0; c < ncopy; c++)
        {
            // On a grid with 8 Angstrom between the copies
            const double shift[3] = { 8.0*(c % 6), 8.0*((c/6) % 6), 8.0*(c/36) };
            for (int a = 0; a < natom; a++)
            {
                const auto &line = atoms[a];
                std::string coords;
                for (int m = 0; m < 3; m++)
                {
                    coords += gmx::formatString("%8.3f", std::stod(line.substr(30 + 8*m, 8)) + shift[m]);
                }
                out << gmx::formatString("HETATM%5d", c*natom + a + 1) << line.substr(11, 11)
                    << gmx::formatString("%4d", c + 1) << line.substr(26, 4) << coords << line.substr(54) << "\n";
            }
        }
        for (int c = 0; c < ncopy; c++)
        {
            for (const auto &line : conect)
            {
                std::istringstream iss(line.substr(6));
                std::string        rec = "CONECT";
                for (int i; iss >> i; )
                {
                    rec += gmx::formatString("%5d", c*natom + i);
                }
                out << rec << "\n";
            }
        }
        out << "END\n";
    }
    auto types = importTypes(pd, boxName);
    ASSERT_EQ(ncopy*ref.size(), types.size());
    for (size_t i = 0; i < types.size(); i++)
    {
        EXPECT_EQ(ref[i % ref.size()], types[i]) << "atom " << i+1;
    }
}
}

}
