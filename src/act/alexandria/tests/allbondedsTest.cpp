/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2021-2023,2026
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cmath>

#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/allbondeds.h"
#include "act/alexandria/molselect.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/import/import.h"
#include "act/molprop/molprop_xml.h"
#include "act/topology/actmol.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class AllBondedsTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
 
        AllBondedsTest () : checker_(this->rootChecker())
        {
            const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, GMX_DOUBLE ? 1e-10 : 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
            //ForceField *mypd     = getForceField("ACM-g");
        }

        void runOBstats(const std::vector<double> &xx)
        {
            const Identifier atom("h");
            OneBonded ob(atom);
            for(const auto &x : xx)
            {
                ob.addPoint(x);
            }
            real   average, sigma;
            size_t N;
            eStats ok = ob.getAverageSigmaN(&average, &sigma, &N);
            if (eStats::OK == ok)
            {
                checker_.checkReal(average, "average");
                checker_.checkReal(sigma, "sigma");
                checker_.checkInteger(static_cast<int>(N), "N");
            }
            else
            {
                checker_.checkString(gmx_stats_message(ok), "Error");
            }
            
        }
        static void TearDownTestCase()
        {
        }

};

TEST_F (AllBondedsTest, OneBondedStats){
    const std::vector<double> xx = { 103, 104, 105, 106.5, 107, 109, 106, 104, 107 };
    runOBstats(xx);
}

TEST_F (AllBondedsTest, OneBondedStatsEmpty){
    const std::vector<double> xx;
    runOBstats(xx);
}

/* Test: The geometries are extracted as in geometry_ff for molecules with
 * virtual sites. The topologies are made without parameters, such that
 * the positions of the virtual sites cannot be constructed, see
 * https://github.com/AlexandriaChemistry/ACT/issues/506
 */
TEST(AllBondedsVsiteTest, ExtractGeometries)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    gmx::test::TestFileManager fileManager;
    ForceField pd;
    readForceField(fileManager.getInputFilePath("forcefields/ACS-pg-vs3.xml"), &pd);
    // As in geometry_ff, such that there are no bonds to shells
    pd.setPolarizable(false);
    std::vector<MolProp> mps;
    MolPropRead(&msghandler, fileManager.getInputFilePath("molprop/molprop.xml").c_str(), &mps);
    const auto selName = fileManager.getTemporaryFilePath("selection.dat");
    FILE *fp      = fopen(selName.c_str(), "w");
    for(const char *mol : { "1-butanol", "acetone", "water" })
    {
        fprintf(fp, "%s|Train\n", mol);
    }
    fclose(fp);
    MolSelect gms;
    gms.read(selName.c_str());

    AllBondeds          bondeds;
    std::vector<ACTMol> actmols;
    EXPECT_NO_THROW(bondeds.extractGeometries(&msghandler, &mps, &actmols, &pd, gms));
    EXPECT_TRUE(msghandler.ok());
    ASSERT_EQ(3u, actmols.size());
    // The molecules have virtual sites without parameters
    size_t nvsite = 0;
    for(auto &mol : actmols)
    {
        size_t molvsite = 0;
        for(const auto &entry : mol.topology()->entries())
        {
            if (isVsite(entry.first))
            {
                for(const auto &vs : entry.second)
                {
                    EXPECT_TRUE(vs->params().empty());
                    molvsite += 1;
                }
            }
        }
        nvsite += molvsite;
        // Their positions cannot be constructed, without them the atoms
        // are at their positions in the experiment.
        auto xexp = mol.experimentConst()[0].getCoordinates();
        if (molvsite > 0)
        {
            EXPECT_THROW(mol.experCoords(xexp), gmx::InternalError) << mol.getMolname();
        }
        std::vector<gmx::RVec> x;
        EXPECT_NO_THROW(x = mol.experCoords(xexp, false)) << mol.getMolname();
        size_t j = 0;
        for(size_t i = 0; i < x.size(); i++)
        {
            if (ActParticle::Atom == mol.topology()->atoms()[i].pType())
            {
                for(int m = 0; m < DIM; m++)
                {
                    EXPECT_EQ(xexp[j][m], x[i][m]);
                }
                j += 1;
            }
        }
        EXPECT_EQ(xexp.size(), j);
    }
    EXPECT_GT(nvsite, 0u);
    // Bonds and angles were found
    gmx::StringOutputStream stream;
    gmx::TextWriter         tw(&stream);
    bondeds.writeSummary(&tw);
    const auto summary = stream.toString();
    EXPECT_NE(std::string::npos, summary.find("BONDS")) << summary;
    EXPECT_NE(std::string::npos, summary.find("ANGLES")) << summary;
}

/*! \brief Make new bonded parameters from the compounds of molprop.xml
 * \param[in]  ffName       Force field in this directory
 * \param[in]  bondFunction Bond function to use instead of that of the force field
 * \param[in]  structures   Structure files in mols with more compounds
 * \param[out] pd           The force field with the new parameters
 */
void updateBondeds(const std::string              &ffName,
                   const char                     *bondFunction,
                   const std::vector<std::string> &structures,
                   ForceField                     *pd)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    gmx::test::TestFileManager fileManager;
    readForceField(fileManager.getInputFilePath(ffName), pd);
    if (bondFunction)
    {
        pd->findForces(InteractionType::BONDS)->setFunction(bondFunction);
    }
    std::vector<MolProp> mps;
    MolPropRead(&msghandler, fileManager.getInputFilePath("molprop/molprop.xml").c_str(), &mps);
    for (const auto &s : structures)
    {
        double qtot = 0;
        importFile(&msghandler, pd, fileManager.getInputFilePath("molecules/" + s), &mps, "",
                   JobType::OPT, false, &qtot, false);
        ASSERT_TRUE(msghandler.ok()) << s;
        // The compound is named after its file
        mps.back().SetMolname(s.substr(0, s.find('.')));
    }
    const auto selName = fileManager.getTemporaryFilePath("selection.dat");
    FILE *fp      = fopen(selName.c_str(), "w");
    for (const auto &mp : mps)
    {
        fprintf(fp, "%s|Train\n", mp.getMolname().c_str());
    }
    fclose(fp);
    MolSelect gms;
    gms.read(selName.c_str());

    AllBondeds          bondeds;
    std::vector<ACTMol> actmols;
    bondeds.extractGeometries(&msghandler, &mps, &actmols, pd, gms);
    ASSERT_TRUE(msghandler.ok());
    ASSERT_EQ(mps.size(), actmols.size());
    gmx::StringOutputStream stream;
    gmx::TextWriter         tw(&stream);
    bondeds.updateForceField(&tw, pd);
}

//! \return The value of a parameter
double value(const ForceFieldParameterMap &params, const std::string &name)
{
    const auto p = params.find(name);
    if (p == params.end())
    {
        ADD_FAILURE() << "No parameter " << name;
        return 0;
    }
    return p->second.value();
}

/* Test: Bond parameters for all bond functions. The lengths are the averages
 * of the compounds, the other parameters are the defaults. Interactions that
 * are not in the compounds are removed.
 */
TEST(AllBondedsUpdateTest, BondFunctions)
{
    std::map<std::string, std::map<std::string, double>> defaults = {
        { "MORSE_BONDS",    { { "De", 400 }, { "beta", 20 }, { "D0", 0 } } },
        { "HUA_BONDS",      { { "De", 400 }, { "b", 20 }, { "c", 0 } } },
        { "HARMONIC_BONDS", { { "kb", 300000 }, { "bondenergy", -200 } } },
        { "CUBIC_BONDS",    { { "kb", 300000 }, { "De", 400 } } }
    };
    for (const auto &fn : defaults)
    {
        ForceField pd;
        updateBondeds("forcefields/ACS-g1.xml", fn.first.c_str(), {}, &pd);
        auto &bonds = pd.findForcesConst(InteractionType::BONDS);
        EXPECT_EQ(fn.first, potentialToString(bonds.potential()));
        ASSERT_FALSE(bonds.parametersConst().empty()) << fn.first;
        for (const auto &b : bonds.parametersConst())
        {
            auto &length = b.second.find("bondlength")->second;
            EXPECT_GT(length.ntrain(), 0u) << b.first.id();
            EXPECT_LT(80, length.value()) << b.first.id();
            EXPECT_GT(160, length.value()) << b.first.id();
            for (const auto &d : fn.second)
            {
                EXPECT_EQ(d.second, value(b.second, d.first)) << fn.first << " " << b.first.id() << " " << d.first;
            }
            if (fn.first == "CUBIC_BONDS")
            {
                // The parameters are stored with the precision of real
                EXPECT_NEAR(3*length.value(), value(b.second, "rmax"),
                            gmx::test::precisionTolerance(1e-6, 1e-7*3*length.value())) << b.first.id();
            }
            EXPECT_EQ(potentialToParameterName(bonds.potential()).size(), b.second.size()) << fn.first;
        }
        // None of the compounds is linear, the linear angles of the force field are removed
        EXPECT_TRUE(pd.findForcesConst(InteractionType::LINEAR_ANGLES).parametersConst().empty()) << fn.first;
    }
}

/* Test: The Urey-Bradley distance is the 1-3 distance from the new bond
 * lengths and the angle, and the Fourier dihedral parameters are one.
 */
TEST(AllBondedsUpdateTest, UreyBradleyFourier)
{
    ForceField pd;
    updateBondeds("forcefields/ACS-g-dih.xml", nullptr, {}, &pd);
    auto &bonds  = pd.findForcesConst(InteractionType::BONDS);
    auto &angles = pd.findForcesConst(InteractionType::ANGLES);
    ASSERT_EQ(Potential::UREY_BRADLEY_ANGLES, angles.potential());
    ASSERT_FALSE(angles.parametersConst().empty());
    for (const auto &a : angles.parametersConst())
    {
        auto atoms = a.first.atoms();
        auto bo    = a.first.bondOrders();
        const Identifier bij({ atoms[0], atoms[1] }, { bo[0] }, CanSwap::Yes);
        const Identifier bkj({ atoms[2], atoms[1] }, { bo[1] }, CanSwap::Yes);
        double r12   = bonds.findParameterTypeConst(bij, "bondlength").value();
        double r23   = bonds.findParameterTypeConst(bkj, "bondlength").value();
        double theta = value(a.second, "ub_angle")*M_PI/180;
        double r13   = std::sqrt(r12*r12 + r23*r23 - 2*r12*r23*std::cos(theta));
        EXPECT_NEAR(r13, value(a.second, "r13"), 1e-3*r13) << a.first.id();
        EXPECT_EQ(30000, value(a.second, "kub")) << a.first.id();
        EXPECT_EQ(300, value(a.second, "kt")) << a.first.id();
    }
    auto &dihs = pd.findForcesConst(InteractionType::PROPER_DIHEDRALS);
    ASSERT_EQ(Potential::FOURIER_DIHEDRALS, dihs.potential());
    ASSERT_FALSE(dihs.parametersConst().empty());
    for (const auto &d : dihs.parametersConst())
    {
        EXPECT_EQ(potentialToParameterName(dihs.potential()).size(), d.second.size()) << d.first.id();
        for (const auto &p : d.second)
        {
            EXPECT_EQ(1, p.second.value()) << d.first.id() << " " << p.first;
        }
    }
}


/* Test: The linear angle of acetonitrile, of which the parameter a puts
 * the middle atom between the outer atoms according to the bond lengths.
 */
TEST(AllBondedsUpdateTest, LinearAngles)
{
    ForceField pd;
    updateBondeds("forcefields/ACS-g1.xml", nullptr, { "acetonitrile.sdf" }, &pd);
    auto &bonds = pd.findForcesConst(InteractionType::BONDS);
    auto &lin   = pd.findForcesConst(InteractionType::LINEAR_ANGLES);
    ASSERT_EQ(1u, lin.parametersConst().size());
    for (const auto &a : lin.parametersConst())
    {
        auto atoms = a.first.atoms();
        auto bo    = a.first.bondOrders();
        const Identifier bij({ atoms[0], atoms[1] }, { bo[0] }, CanSwap::Yes);
        const Identifier bjk({ atoms[1], atoms[2] }, { bo[1] }, CanSwap::Yes);
        double rij = bonds.findParameterTypeConst(bij, "bondlength").value();
        double rjk = bonds.findParameterTypeConst(bjk, "bondlength").value();
        EXPECT_NEAR(rjk/(rij+rjk), value(a.second, "a"), 1e-6) << a.first.id();
        EXPECT_EQ(20000, value(a.second, "klin")) << a.first.id();
        EXPECT_EQ(2u, a.second.size()) << a.first.id();
        auto &alin = a.second.find("a")->second;
        EXPECT_EQ(1u, alin.ntrain()) << a.first.id();
        EXPECT_LE(alin.minimum(), alin.value()) << a.first.id();
        EXPECT_GE(alin.maximum(), alin.value()) << a.first.id();
    }
}

}

} // namespace
