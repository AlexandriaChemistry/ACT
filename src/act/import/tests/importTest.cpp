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

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "act/basics/identifier.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/particletype.h"
#include "act/import/import.h"
#include "act/import/import_utils.h"
#include "act/molprop/experiment.h"
#include "act/molprop/molprop.h"
#include "act/testutils/forcefield_utils.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/stringutil.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
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

/*! \brief The atom types and the bond orders of a molecule, each sorted
 * \param[in] molprops The molecule
 * \return the atom types, and the bond orders as text
 */
static std::pair<std::vector<std::string>, std::vector<std::string>>
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

// Test: A PDB file as ACT writes it, with each bond once in the CONECT
// records and so without bond orders, gives the same molecule as its
// SMILES string: the bond orders are inferred. Before, RDKit gave the
// atoms implicit hydrogens, which made pyridine piperidine and the oxygen
// of acetone a hydroxyl.
TEST(ImportPdbTest, BondOrdersOfConectWithoutThem)
{
    auto pd = getForceField("gaff_hydrogen_types");
    for (const auto &[file, smiles] : std::vector<std::pair<std::string, std::string>>{
             { "pyridine-act.pdb", "c1ccncc1" }, { "acetone-act.pdb", "CC(=O)C" } })
    {
        SCOPED_TRACE(file);
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        std::vector<MolProp> fromFile, fromSmiles;
        double               qtot = 0;
        importFile(&msghandler, pd, gmx::test::TestFileManager::getInputFilePath("molecules/" + file),
                   &fromFile, "minimum", JobType::OPT, false, &qtot, false);
        importSmiles(&msghandler, pd, smiles, &fromSmiles, false);
        EXPECT_TRUE(msghandler.ok());
        const auto read = typesAndBondOrders(fromFile);
        const auto made = typesAndBondOrders(fromSmiles);
        EXPECT_EQ(made.first, read.first);
        EXPECT_EQ(made.second, read.second);
        // No hydrogens were added
        ASSERT_EQ(1u, fromFile.size());
        EXPECT_EQ(fromSmiles[0].experimentConst()[0].calcAtomConst().size(),
                  fromFile[0].experimentConst()[0].calcAtomConst().size());
    }
}

// Test: Helium is read from pdb and xyz files, in which RDKit determines
// the bonds. RDKit gives a helium atom without bonds the formal charge that
// fills its octet, minus six, and then found no bond orders for any file
// with helium and another atom.
TEST(ImportHeliumTest, FilesWithoutBonds)
{
    auto       pd     = getForceField("gaff_hydrogen_types");
    const auto helium = importTypes(pd, gmx::test::TestFileManager::getInputFilePath("molecules/helium-4.sdf"));
    ASSERT_EQ(1u, helium.size());
    auto       water  = importTypes(pd, gmx::test::TestFileManager::getInputFilePath("molecules/water-3-oep.log.pdb"));
    ASSERT_EQ(3u, water.size());
    std::sort(water.begin(), water.end());
    // The file, and the numbers of helium atoms and of bonds
    for (const auto &[file, nhelium, nbond] : std::vector<std::tuple<std::string, size_t, size_t>>{
             { "helium#helium.pdb", 2, 0 }, { "helium#helium.xyz", 2, 0 }, { "helium#water.xyz", 1, 2 } })
    {
        SCOPED_TRACE(file);
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        std::vector<MolProp> molprops;
        double               qtot = 0;
        importFile(&msghandler, pd, gmx::test::TestFileManager::getInputFilePath("molecules/" + file),
                   &molprops, "minimum", JobType::OPT, false, &qtot, false);
        EXPECT_TRUE(msghandler.ok());
        ASSERT_EQ(1u, molprops.size());
        EXPECT_EQ(0, qtot);
        std::vector<std::string> atypes;
        for (const auto &ca : molprops[0].experimentConst()[0].calcAtomConst())
        {
            atypes.push_back(ca.getObtype());
        }
        // The helium atoms come first in these files
        ASSERT_LE(nhelium, atypes.size());
        for (size_t i = 0; i < nhelium; i++)
        {
            EXPECT_EQ(helium[0], atypes[i]);
        }
        std::vector<std::string> other(atypes.begin() + nhelium, atypes.end());
        std::sort(other.begin(), other.end());
        EXPECT_EQ(nhelium == 1 ? water : std::vector<std::string>{}, other);
        EXPECT_EQ(nbond, molprops[0].bondsConst().size());
        // Two compounds without charge
        const auto &fragments = molprops[0].fragments();
        ASSERT_EQ(2u, fragments.size());
        for (const auto &fragment : fragments)
        {
            EXPECT_EQ(0, fragment.charge());
        }
    }
}

// Test: The bond orders of formate are written with whole bonds, also with
// a helium atom next to it.
TEST(ImportHeliumTest, WholeBondOrders)
{
    // Carbon, hydrogen and two oxygens, in nm
    std::vector<int>                       elements = { 6, 1, 8, 8 };
    std::vector<gmx::RVec>                 x        = { { 0, 0, 0 }, { 0, -0.112, 0 },
                                                        { 0.110, 0.062, 0 }, { -0.110, 0.062, 0 } };
    const std::vector<std::pair<int, int>> bonds    = { { 0, 1 }, { 0, 2 }, { 0, 3 } };
    const std::vector<double>              orders   = { 1, 1.5, 1.5 };
    for (bool withHelium : { false, true })
    {
        SCOPED_TRACE(withHelium ? "with helium" : "formate");
        if (withHelium)
        {
            elements.push_back(2);
            x.emplace_back(0, 0, 0.35);
        }
        const auto whole = wholeBondOrders(elements, x, bonds, orders, -1);
        ASSERT_EQ(3u, whole.size());
        EXPECT_EQ(1, whole[0]);
        EXPECT_EQ(3, whole[1] + whole[2]);
    }
}

//! What the tests of ions look at in a file
struct ImportedIons
{
    //! Whether the file was read without error
    bool                     ok = false;
    //! The atom types
    std::vector<std::string> atypes;
    //! The number of bonds
    size_t                   nbond = 0;
    //! The charges of the fragments
    std::vector<int>         charges;
};

/*! \brief Read a file with ions
 * \param[in] file       The file in the molecules directory of the test data
 * \param[in] userqtot   Whether the total charge is given
 * \param[in] qtot       The total charge
 * \param[in] pd         The force field
 * \return what was read
 */
static ImportedIons importIons(const std::string &file, bool userqtot, double qtot,
                               const ForceField *pd)
{
    MsgHandler   msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    std::vector<MolProp> molprops;
    importFile(&msghandler, pd, gmx::test::TestFileManager::getInputFilePath("molecules/" + file),
               &molprops, "minimum", JobType::OPT, userqtot, &qtot, false);
    ImportedIons ions;
    ions.ok = msghandler.ok() && molprops.size() == 1;
    if (ions.ok)
    {
        for (const auto &ca : molprops[0].experimentConst()[0].calcAtomConst())
        {
            ions.atypes.push_back(ca.getObtype());
        }
        ions.nbond = molprops[0].bondsConst().size();
        for (const auto &fragment : molprops[0].fragments())
        {
            ions.charges.push_back(fragment.charge());
        }
    }
    return ions;
}

/*! \brief Read a file with ions
 * \param[in] file       The file in the molecules directory of the test data
 * \param[in] userqtot   Whether the total charge is given
 * \param[in] qtot       The total charge
 * \param[in] forceField The name of the force field
 * \return what was read
 */
static ImportedIons importIons(const std::string &file, bool userqtot, double qtot,
                               const std::string &forceField = "gaff_hydrogen_types")
{
    return importIons(file, userqtot, qtot, getForceField(forceField));
}

// Test: An alkali halide from an xyz file is the same pair of ions as from
// an sdf file with the charges. RDKit made a bond between the atoms from
// their distance and then found no bond orders.
TEST(ImportIonTest, AlkaliHalideFromXyz)
{
    for (const std::string salt : { "NaCl", "KBr" })
    {
        SCOPED_TRACE(salt);
        const auto sdf = importIons(salt + ".sdf", false, 0);
        const auto xyz = importIons(salt + ".xyz", false, 0);
        ASSERT_TRUE(sdf.ok);
        ASSERT_TRUE(xyz.ok);
        EXPECT_EQ(sdf.atypes, xyz.atypes);
        EXPECT_EQ(0u, xyz.nbond);
        EXPECT_EQ((std::vector<int>{ 1, -1 }), xyz.charges);
        EXPECT_EQ(sdf.charges, xyz.charges);
    }
}

// Test: A sodium ion next to a water molecule, without bonds in the file,
// is the same as in the file with the bonds of water. The ion has the
// charge one and no bond to the oxygen.
TEST(ImportIonTest, SodiumAndWater)
{
    const auto bonds = importIons("sodium#water.pdb", true, 1);
    ASSERT_TRUE(bonds.ok);
    ASSERT_EQ(4u, bonds.atypes.size());
    EXPECT_EQ("Na+", bonds.atypes[0]);
    // The total charge is given for the xyz file, and is in the pdb file
    for (const auto &[file, userqtot] : std::vector<std::pair<std::string, bool>>{
             { "sodium#water.xyz", true }, { "sodium#water-nobonds.pdb", false } })
    {
        SCOPED_TRACE(file);
        const auto ions = importIons(file, userqtot, userqtot ? 1 : 0);
        ASSERT_TRUE(ions.ok);
        EXPECT_EQ(bonds.atypes, ions.atypes);
        EXPECT_EQ(2u, ions.nbond);
        EXPECT_EQ((std::vector<int>{ 1, 0 }), ions.charges);
        EXPECT_EQ(bonds.charges, ions.charges);
    }
    // Without the charge of the ion in the total charge there is no structure
    EXPECT_FALSE(importIons("sodium#water.xyz", false, 0).ok);
}

// Test: A single atom has the total charge, and an alkaline earth metal
// without bonds the charge two.
TEST(ImportIonTest, SingleIonAndMagnesium)
{
    const auto sodium = importIons("sodium.xyz", true, 1);
    ASSERT_TRUE(sodium.ok);
    EXPECT_EQ((std::vector<std::string>{ "Na+" }), sodium.atypes);
    EXPECT_EQ((std::vector<int>{ 1 }), sodium.charges);
    const auto chloride = importIons("chloride.xyz", true, -1);
    ASSERT_TRUE(chloride.ok);
    EXPECT_EQ((std::vector<std::string>{ "Cl-" }), chloride.atypes);
    EXPECT_EQ((std::vector<int>{ -1 }), chloride.charges);
    // Without a charge it is an atom
    const auto chlorine = importIons("chloride.xyz", false, 0);
    ASSERT_TRUE(chlorine.ok);
    EXPECT_EQ((std::vector<std::string>{ "cl" }), chlorine.atypes);
    EXPECT_EQ((std::vector<int>{ 0 }), chlorine.charges);
    // This force field has the alkaline earth metals
    const auto magnesium = importIons("magnesium#water.xyz", true, 2, "ACS-g1");
    ASSERT_TRUE(magnesium.ok);
    ASSERT_EQ(4u, magnesium.atypes.size());
    EXPECT_EQ("Mg2+", magnesium.atypes[0]);
    EXPECT_EQ(2u, magnesium.nbond);
    EXPECT_EQ((std::vector<int>{ 2, 0 }), magnesium.charges);
}

// Test: Aluminium is an ion with the charge three in a file without bonds,
// also with four chlorides around it, for which RDKit found bonds.
TEST(ImportIonTest, Aluminium)
{
    // The force fields of the tests have no aluminium. This one gets the
    // ion, with the parameters of the magnesium ion.
    ForceField  pd = *getForceField("ACS-g1");
    const auto *mg = pd.findParticleType("Mg2+");
    ASSERT_NE(nullptr, mg);
    ParticleType al(Identifier("Al3+"), "aluminium ion", mg->apType());
    for (const auto &[key, value] : mg->optionsConst())
    {
        if (key != "element" && key != "atomnumber")
        {
            al.setOption(key, value);
        }
    }
    al.setOption("element", "Al");
    al.setOption("atomnumber", "13");
    for (const auto &[type, param] : mg->parametersConst())
    {
        al.addForceFieldParameter(type, param);
    }
    pd.addParticleType(al);

    const auto water = importIons("aluminium#water.xyz", true, 3, &pd);
    ASSERT_TRUE(water.ok);
    ASSERT_EQ(4u, water.atypes.size());
    EXPECT_EQ("Al3+", water.atypes[0]);
    EXPECT_EQ(2u, water.nbond);
    EXPECT_EQ((std::vector<int>{ 3, 0 }), water.charges);
    // Without the charge of the ion in the total charge there is no structure
    EXPECT_FALSE(importIons("aluminium#water.xyz", false, 0, &pd).ok);

    const auto chloride = importIons("tetrachloroaluminate.xyz", true, -1, &pd);
    ASSERT_TRUE(chloride.ok);
    EXPECT_EQ((std::vector<std::string>{ "Al3+", "Cl-", "Cl-", "Cl-", "Cl-" }), chloride.atypes);
    EXPECT_EQ(0u, chloride.nbond);
    EXPECT_EQ((std::vector<int>{ 3, -1, -1, -1, -1 }), chloride.charges);

    // Formate with an aluminium ion next to it, see the test with helium
    const std::vector<int>                 elements = { 6, 1, 8, 8, 13 };
    const std::vector<gmx::RVec>           x        = { { 0, 0, 0 }, { 0, -0.112, 0 }, { 0.110, 0.062, 0 },
                                                        { -0.110, 0.062, 0 }, { 0, 0, 0.35 } };
    const std::vector<std::pair<int, int>> bonds    = { { 0, 1 }, { 0, 2 }, { 0, 3 } };
    const auto whole = wholeBondOrders(elements, x, bonds, { 1, 1.5, 1.5 }, 2);
    ASSERT_EQ(3u, whole.size());
    EXPECT_EQ(1, whole[0]);
    EXPECT_EQ(3, whole[1] + whole[2]);
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
