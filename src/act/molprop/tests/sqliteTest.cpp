/*
 * This source file is part of the Alexandria program.
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

/*! \internal \brief
 * Tests of reading molecular properties from an SQLite3 database.
 */

#include "actpre.h"

#include "config.h"

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#if HAVE_LIBSQLITE3
#include <sqlite3.h>
#endif

#include "act/basics/msg_handler.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_sqlite3.h"
#include "act/utility/units.h"

#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

#if HAVE_LIBSQLITE3

//! The database, with the tables and columns that ReadSqlite3 uses
const char *databaseSql = R"(
CREATE TABLE molecules(molid INTEGER, iupac TEXT, cas TEXT, csid TEXT);
CREATE TABLE synonyms(molid INTEGER, name TEXT);
CREATE TABLE classification(classid INTEGER, class TEXT);
CREATE TABLE link_mol_class(molid INTEGER, classid INTEGER);
CREATE TABLE proptypes(propid INTEGER, prop TEXT, unit_text TEXT);
CREATE TABLE datasource(srcid INTEGER, theory INTEGER);
CREATE TABLE phasetype(phaseid INTEGER, phase TEXT);
CREATE TABLE reference(refid INTEGER, ref TEXT);
CREATE TABLE molproperty(molpropid INTEGER, molid INTEGER, propid INTEGER, phaseid INTEGER,
                         srcid INTEGER, temperature REAL, value REAL, error REAL, preferred INTEGER);
CREATE TABLE link_molprop_ref(molpropid INTEGER, refid INTEGER);
INSERT INTO molecules VALUES (1, 'water', '7732-18-5', '937'), (2, 'acetone', '67-64-1', '175');
INSERT INTO synonyms VALUES (1, 'water'), (2, 'propanone');
INSERT INTO classification VALUES (1, 'inorganic'), (2, 'solvent'), (3, 'ketone');
INSERT INTO link_mol_class VALUES (1, 1), (1, 2), (2, 2), (2, 3);
INSERT INTO proptypes VALUES (1, 'dipole', 'D'), (2, 'DeltaHform', 'kJ/mol'),
                             (3, 'polarizability', 'Angstrom3'), (4, 'boiling point', 'K');
INSERT INTO datasource VALUES (1, 0), (2, 1);
INSERT INTO phasetype VALUES (1, 'gas'), (2, 'liquid');
INSERT INTO reference VALUES (1, 'Lide2004a'), (2, 'Spoel2022a');
-- water: preferred experiment, not preferred experiment, other temperature,
-- liquid, theory, unknown property
INSERT INTO molproperty VALUES (1, 1, 1, 1, 1, 298.15, 1.85, 0.01, 1),
                               (2, 1, 1, 1, 1, 298.15, 1.80, 0.02, 0),
                               (3, 1, 1, 1, 1, 0,      1.86, 0.01, 1),
                               (4, 1, 1, 2, 1, 298.15, 2.90, 0.10, 1),
                               (5, 1, 2, 1, 1, 298.15, -241.8, 0.1, 1),
                               (6, 1, 3, 1, 2, 298.15, 1.47, 0, 0),
                               (7, 1, 4, 1, 1, 298.15, 373.15, 0, 1),
                               (8, 2, 1, 1, 1, 298.15, 2.88, 0.03, 1);
INSERT INTO link_molprop_ref VALUES (1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 2), (7, 1), (8, 1);
)";

class Sqlite3Test : public ::testing::Test
{
protected:
    //! Removes the database at the end of the test
    gmx::test::TestFileManager fileManager_;
    //! The database file
    std::string                database_;

    Sqlite3Test() : database_(fileManager_.getTemporaryFilePath("molprops.db"))
    {
        std::remove(database_.c_str());
        sqlite3 *db = nullptr;
        EXPECT_EQ(SQLITE_OK, sqlite3_open(database_.c_str(), &db));
        char *error = nullptr;
        EXPECT_EQ(SQLITE_OK, sqlite3_exec(db, databaseSql, nullptr, nullptr, &error))
            << (error ? error : "");
        sqlite3_free(error);
        sqlite3_close(db);
    }

    //! \return molecules with these names, read from the database at T
    std::vector<MolProp> read(const std::vector<std::string> &names, double T)
    {
        std::vector<MolProp> mps(names.size());
        for (size_t i = 0; i < names.size(); i++)
        {
            mps[i].SetMolname(names[i]);
        }
        MsgHandler msghandler;
        ReadSqlite3(&msghandler, database_.c_str(), &mps, T);
        return mps;
    }
};

//! \return the properties mpo of all experiments of a molecule
std::vector<const GenericProperty *> properties(const MolProp &mp, MolPropObservable mpo)
{
    std::vector<const GenericProperty *> props;
    for (const auto &ei : mp.experimentConst())
    {
        const auto p = ei.propertiesConst().find(mpo);
        if (p != ei.propertiesConst().end())
        {
            for (const auto &gp : p->second)
            {
                props.push_back(gp.get());
            }
        }
    }
    return props;
}

// Test: Only the preferred experiments and the calculations of the gas
// phase at the temperature are read, with the reference, categories and
// identifiers of the molecule
TEST_F(Sqlite3Test, Water)
{
    auto mps   = read({ "water" }, 298.15);
    auto &water = mps[0];
    // Dipole, enthalpy of formation and polarizability
    EXPECT_EQ(3u, water.experimentConst().size());
    auto dipole = properties(water, MolPropObservable::DIPOLE);
    ASSERT_EQ(1u, dipole.size());
    EXPECT_NEAR(convertToGromacs(1.85, "D"), dipole[0]->getValue(), 1e-12);
    EXPECT_DOUBLE_EQ(298.15, dipole[0]->getTemperature());
    auto dhf = properties(water, MolPropObservable::DHFORM);
    ASSERT_EQ(1u, dhf.size());
    EXPECT_NEAR(convertToGromacs(-241.8, "kJ/mol"), dhf[0]->getValue(), 1e-12);
    auto alpha = properties(water, MolPropObservable::POLARIZABILITY);
    ASSERT_EQ(1u, alpha.size());
    EXPECT_NEAR(convertToGromacs(1.47, "Angstrom3"), alpha[0]->getValue(), 1e-12);
    for (const auto &ei : water.experimentConst())
    {
        EXPECT_EQ(ei.hasProperty(MolPropObservable::POLARIZABILITY) ? "Spoel2022a" : "Lide2004a",
                  ei.getReference());
        EXPECT_EQ("minimum", ei.getConformation());
    }
    EXPECT_EQ((std::vector<std::string>{ "inorganic", "solvent" }), water.categoryConst());
    EXPECT_EQ("7732-18-5", water.getCas());
    EXPECT_EQ("937", water.getCid());
}

// Test: A molecule is found by its synonym, and one without a synonym is
// left as it is
TEST_F(Sqlite3Test, SynonymsAndMissing)
{
    auto mps = read({ "propanone", "methane" }, 298.15);
    auto dipole = properties(mps[0], MolPropObservable::DIPOLE);
    ASSERT_EQ(1u, dipole.size());
    EXPECT_NEAR(convertToGromacs(2.88, "D"), dipole[0]->getValue(), 1e-12);
    EXPECT_EQ((std::vector<std::string>{ "solvent", "ketone" }), mps[0].categoryConst());
    EXPECT_EQ("67-64-1", mps[0].getCas());
    EXPECT_TRUE(mps[1].experimentConst().empty());
    EXPECT_TRUE(mps[1].categoryConst().empty());
    EXPECT_TRUE(mps[1].getCas().empty());
}

// Test: Only the properties at the temperature are read
TEST_F(Sqlite3Test, Temperature)
{
    auto mps    = read({ "water" }, 0);
    auto dipole = properties(mps[0], MolPropObservable::DIPOLE);
    ASSERT_EQ(1u, dipole.size());
    EXPECT_NEAR(convertToGromacs(1.86, "D"), dipole[0]->getValue(), 1e-12);
    EXPECT_EQ(1u, mps[0].experimentConst().size());
}

// Test: Without a database nothing is read
TEST_F(Sqlite3Test, NoDatabase)
{
    std::vector<MolProp> mps(1);
    mps[0].SetMolname("water");
    MsgHandler msghandler;
    ReadSqlite3(&msghandler, nullptr, &mps, 298.15);
    EXPECT_TRUE(mps[0].experimentConst().empty());
}

#endif

} // namespace

} // namespace alexandria
