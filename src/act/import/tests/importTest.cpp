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

#include <map>

#include <gtest/gtest.h>

//#include "act/alexandria/actmol.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/import/import.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class ImportTest : public ::testing::TestWithParam<std::tuple<bool, std::tuple<std::string, double>, bool > >
//: public gmx::test::CommandLineTestBase
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
        pd       = getForceField("../../alexandria/tests/ACS-g");
        oneH     = std::get<0>(GetParam());
        bondTest = std::get<2>(GetParam());
        molname  = std::get<0>(std::get<1>(GetParam()));
        qtot     = std::get<1>(std::get<1>(GetParam()));
    }
    
 protected:
    void runTest()
    {
        gmx::test::TestReferenceChecker  checker_(refData_.rootChecker());
        const char                      *conf     = (char *)"minimum";
        std::string                      dataName;
        std::vector<alexandria::MolProp> molprops;

        std::string dirmolname("../../alexandria/tests/mols/");
        dirmolname     += molname;
        dataName        = gmx::test::TestFileManager::getInputFilePath(dirmolname);
        bool       userqtot = false;
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        importFile(&msghandler, pd, dataName.c_str(), &molprops,
                   conf, JobType::OPT, userqtot, &qtot, oneH);
        EXPECT_TRUE(msghandler.ok());
        for(auto &molprop: molprops)
        {
            checker_.checkString(molname, "molname");
            checker_.checkDouble(qtot, "qtot");
            EXPECT_TRUE(molprop.experimentConst().size() > 0);
            if (molprop.experimentConst().size() == 0)
            {
                continue;
            }
            auto                     exper = molprop.experimentConst().begin();
            std::vector<std::string> atypes;
            for (auto &ca : exper->calcAtomConst())
            {
                atypes.push_back(ca.getObtype());
            }
            if (!bondTest)
            {
                checker_.checkInteger(static_cast<int>(atypes.size()), molname.c_str());
                checker_.checkSequence(atypes.begin(), atypes.end(), "atomtypes");
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
        { "water_dimer.pdb", 0 }
    };
}

static std::vector<bool> get_ab()
{
    return { false, true };
}

static std::vector<bool> get_ba()
{
    return { true, false };
}

TEST_P (ImportTest, All)
{
    runTest();
}

auto my_files = get_files();
auto my_btest = get_ab();
auto my_oneH  = get_ba();

INSTANTIATE_TEST_CASE_P(IT, ImportTest, ::testing::Combine(::testing::ValuesIn(my_oneH),  ::testing::ValuesIn(my_files), ::testing::ValuesIn(my_btest)));

}

}
