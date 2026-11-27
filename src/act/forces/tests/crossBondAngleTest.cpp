/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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
 * Tests of the bond-angle cross term, the potential CROSS_BOND_ANGLES.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/gromacs_top.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/basics/mutability.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/potential.h"
#include "act/forces/forcecomputer.h"
#include "act/import/import.h"
#include "act/molprop/molprop.h"
#include "act/topology/actmol.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/* The bond-angle cross term of water in the flexible TIP4P model of
 * Lawrence and Skinner, Chem. Phys. Lett. 372 (2003) 842, as in the
 * GROMACS topology tip4pflex.itp: the angle H O H of type 4 with
 * r1e, r2e, r3e (nm) and krtheta (kJ/mol nm^2).
 */
//! Equilibrium length of the first bond (nm)
const double cR1e = 0.09347;
//! Equilibrium length of the second bond (nm)
const double cR2e = 0.09347;
//! Equilibrium distance between the hydrogen atoms (nm)
const double cR3e = 0.1526;
//! Force constant (kJ/mol nm^2)
const double cKrt = -81901.5;

/* A water molecule that is slightly off equilibrium (nm): the bonds are
 * 0.09522 and 0.09728 nm long and the hydrogen atoms 0.15209 nm apart.
 */
//! Position of the oxygen atom
const double cXO[3]  = { 0.001, -0.002, 0.0005 };
//! Position of the first hydrogen atom
const double cXH1[3] = { 0.096, 0.004, -0.002 };
//! Position of the second hydrogen atom
const double cXH2[3] = { -0.029, 0.0905, 0.003 };

/* Energy and forces for these parameters and positions from Eqns. 194 and
 * 195 of the GROMACS reference manual (bonded interactions, 2026.4),
 * evaluated with 40 decimal digits in a Python script.
 */
//! Energy (kJ/mol)
const double cEnergy     = 0.230815456837317;
//! Force on the first hydrogen atom (kJ/mol nm)
const double cForceH1[3] = { 332.650776798099, -261.484635044086, -13.8730423173901 };
//! Force on the oxygen atom (kJ/mol nm)
const double cForceO[3]  = { 28.6267557605897, 42.1104093206863, -0.023016123940194 };
//! Force on the second hydrogen atom (kJ/mol nm)
const double cForceH2[3] = { -361.277532558689, 219.374225723399, 13.8960584413303 };

/* The same from GROMACS 2026.1 (mixed precision): a single point with
 * tip4pflex.itp for these positions, the energy term BA-Cross, and the
 * forces from a topology with this angle only, as gmx dump prints them.
 */
//! Energy (kJ/mol)
const double cGromacsEnergy     = 0.230822;
//! Force on the first hydrogen atom (kJ/mol nm)
const double cGromacsForceH1[3] = { 332.649, -261.484, -13.8730 };
//! Force on the oxygen atom (kJ/mol nm)
const double cGromacsForceO[3]  = { 28.6276, 42.1117, -0.0230188 };
//! Force on the second hydrogen atom (kJ/mol nm)
const double cGromacsForceH2[3] = { -361.277, 219.373, 13.8960 };

//! A water molecule with its force field
struct Water
{
    //! The force field
    ForceField             pd;
    //! The molecule
    ACTMol                 mol;
    //! The indices of the first hydrogen, the oxygen and the second hydrogen atom
    std::vector<int>       hoh;
    //! The positions
    std::vector<gmx::RVec> x;
};

/*! \brief Make a water molecule
 * \param[in]  pd  The force field
 * \param[out] mol The molecule with its topology
 */
void makeMolecule(ForceField *pd, ACTMol *mol)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    std::vector<MolProp> mps;
    double               qtot = 0;
    importFile(&msghandler, pd, gmx::test::TestFileManager::getInputFilePath("molecules/water.sdf"),
               &mps, "minimum", JobType::OPT, false, &qtot, true);
    ASSERT_TRUE(msghandler.ok());
    ASSERT_EQ(1u, mps.size());
    mol->Merge(&mps[0]);
    mol->GenerateTopology(&msghandler, pd, missingParameters::Ignore);
    ASSERT_TRUE(msghandler.ok());
}

/*! \brief Make a water molecule with the bond-angle cross term for its angle
 *
 * The force field ACS-g2 gets CROSS_BOND_ANGLES as the potential for the
 * angles, with one angle, that of water.
 * \param[in]  krt   The force constant (kJ/mol nm^2)
 * \param[out] water The molecule, off equilibrium
 */
void makeWater(double krt, Water *water)
{
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACS-g2.xml"), &water->pd);
    Identifier angleId;
    {
        // The identifier of the angle, from the force field as it is
        ACTMol harmonic;
        makeMolecule(&water->pd, &harmonic);
        ASSERT_TRUE(harmonic.topology()->hasEntry(InteractionType::ANGLES));
        const auto &angles = harmonic.topology()->entry(InteractionType::ANGLES);
        ASSERT_EQ(1u, angles.size());
        angleId = angles[0]->id();
    }
    auto *fs = water->pd.findForces(InteractionType::ANGLES);
    fs->eraseParameter();
    fs->setFunction("CROSS_BOND_ANGLES");
    // Lengths are in pm in a force field
    const std::vector<std::pair<const char *, double>> lengths = {
        { "r1e", 1000*cR1e }, { "r2e", 1000*cR2e }, { "r3e", 1000*cR3e }
    };
    for (const auto &len : lengths)
    {
        fs->addParameter(angleId, len.first,
                         ForceFieldParameter("pm", len.second, 0, 1, 0.9*len.second, 1.1*len.second,
                                             Mutability::Bounded, true, true));
    }
    fs->addParameter(angleId, "krt",
                     ForceFieldParameter("kJ/mol/nm2", krt, 0, 1, -200000, 200000,
                                         Mutability::Bounded, true, false));
    makeMolecule(&water->pd, &water->mol);
    ASSERT_TRUE(water->mol.topology()->hasEntry(InteractionType::ANGLES));
    const auto &angles = water->mol.topology()->entry(InteractionType::ANGLES);
    ASSERT_EQ(1u, angles.size());
    water->hoh = angles[0]->atomIndices();
    ASSERT_EQ(3u, water->hoh.size());
    ASSERT_EQ(3u, water->mol.topology()->atoms().size());
    ASSERT_EQ("O", water->mol.topology()->atoms()[water->hoh[1]].element());
    water->x.resize(3);
    for (int m = 0; m < DIM; m++)
    {
        water->x[water->hoh[0]][m] = cXH1[m];
        water->x[water->hoh[1]][m] = cXO[m];
        water->x[water->hoh[2]][m] = cXH2[m];
    }
}

/*! \brief Energies and forces of a water molecule
 * \param[in]  water    The molecule
 * \param[out] forces   The forces
 * \param[out] energies The energies
 */
void compute(const Water                       &water,
             std::vector<gmx::RVec>            *forces,
             std::map<InteractionType, double> *energies)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    const ForceComputer forceComp;
    auto                x = water.x;
    forces->assign(x.size(), { 0, 0, 0 });
    forceComp.compute(&msghandler, &water.pd, water.mol.topology(), &x, forces, energies);
    ASSERT_TRUE(msghandler.ok());
}

// Test: The names of the potential and of its parameters, in the order of GROMACS
TEST(CrossBondAngleTest, Names)
{
    Potential p = Potential::NONE;
    ASSERT_TRUE(stringToPotential("CROSS_BOND_ANGLES", &p));
    EXPECT_EQ(Potential::CROSS_BOND_ANGLES, p);
    EXPECT_EQ("CROSS_BOND_ANGLES", potentialToString(Potential::CROSS_BOND_ANGLES));
    // The name in GROMACS
    p = Potential::NONE;
    ASSERT_TRUE(stringToPotential("CROSS_BOND_ANGLE", &p));
    EXPECT_EQ(Potential::CROSS_BOND_ANGLES, p);
    EXPECT_EQ(F_CROSS_BOND_ANGLES, potentialToGromacsType(Potential::CROSS_BOND_ANGLES));
    const auto names = potentialToParameterName(Potential::CROSS_BOND_ANGLES);
    ASSERT_EQ(static_cast<size_t>(cbaNR), names.size());
    EXPECT_STREQ("r1e", names[cbaR1E]);
    EXPECT_STREQ("r2e", names[cbaR2E]);
    EXPECT_STREQ("r3e", names[cbaR3E]);
    EXPECT_STREQ("krt", names[cbaKRT]);
}

// Test: Energy and forces of a water molecule that is off equilibrium, with
// the parameters of flexible TIP4P, against Eqns. 194 and 195 of the GROMACS
// manual and against GROMACS
TEST(CrossBondAngleTest, FlexibleTip4pWater)
{
    Water water;
    makeWater(cKrt, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    // The parameters arrive in the topology in nm and kJ/mol nm^2
    const auto &params = water.mol.topology()->entry(InteractionType::ANGLES)[0]->params();
    ASSERT_EQ(static_cast<size_t>(cbaNR), params.size());
    EXPECT_NEAR(cR1e, params[cbaR1E], 1e-12);
    EXPECT_NEAR(cR2e, params[cbaR2E], 1e-12);
    EXPECT_NEAR(cR3e, params[cbaR3E], 1e-12);
    EXPECT_NEAR(cKrt, params[cbaKRT], 1e-6);

    std::vector<gmx::RVec>            forces;
    std::map<InteractionType, double> energies;
    compute(water, &forces, &energies);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    // The positions are rounded to single precision in a single precision
    // build, which changes the small deviations from the equilibrium lengths
    const double etol = gmx::test::precisionTolerance(1e-10, 1e-4);
    const double ftol = gmx::test::precisionTolerance(1e-7, 5e-2);
    ASSERT_EQ(1u, energies.count(InteractionType::ANGLES));
    EXPECT_NEAR(cEnergy, energies[InteractionType::ANGLES], etol);
    // GROMACS has single precision positions
    EXPECT_NEAR(cGromacsEnergy, energies[InteractionType::ANGLES], gmx::test::precisionTolerance(2e-5, 1e-4));

    // The same from the equation, written out here
    double rij[3], rkj[3], rik[3], nij2 = 0, nkj2 = 0, nik2 = 0;
    for (int m = 0; m < DIM; m++)
    {
        // The positions as the force computer got them
        double xi = water.x[water.hoh[0]][m];
        double xj = water.x[water.hoh[1]][m];
        double xk = water.x[water.hoh[2]][m];
        rij[m]    = xi - xj;
        rkj[m]    = xk - xj;
        rik[m]    = xi - xk;
        nij2     += rij[m]*rij[m];
        nkj2     += rkj[m]*rkj[m];
        nik2     += rik[m]*rik[m];
    }
    double nij = std::sqrt(nij2), nkj = std::sqrt(nkj2), nik = std::sqrt(nik2);
    EXPECT_NEAR(cKrt*(nik - cR3e)*(nij - cR1e + nkj - cR2e), energies[InteractionType::ANGLES],
                gmx::test::precisionTolerance(1e-10, 1e-7));

    // The forces of the cross term alone: the difference with the forces
    // without it
    Water without;
    makeWater(0, &without);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    std::vector<gmx::RVec>            forces0;
    std::map<InteractionType, double> energies0;
    compute(without, &forces0, &energies0);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    EXPECT_EQ(0, energies0[InteractionType::ANGLES]);
    const std::vector<const double *> reference = { cForceH1, cForceO, cForceH2 };
    const std::vector<const double *> gromacs   = { cGromacsForceH1, cGromacsForceO, cGromacsForceH2 };
    double fsum[3] = { 0, 0, 0 };
    for (size_t a = 0; a < reference.size(); a++)
    {
        for (int m = 0; m < DIM; m++)
        {
            double f = forces[water.hoh[a]][m] - forces0[water.hoh[a]][m];
            EXPECT_NEAR(reference[a][m], f, ftol) << "atom " << a << " dimension " << m;
            EXPECT_NEAR(gromacs[a][m], f, gmx::test::precisionTolerance(3e-3, 5e-2))
                << "GROMACS, atom " << a << " dimension " << m;
            fsum[m] += f;
        }
    }
    // Eqn. 195 for the first hydrogen atom
    for (int m = 0; m < DIM; m++)
    {
        double f195 = -cKrt*((nik - cR3e)*rij[m]/nij + (nij - cR1e + nkj - cR2e)*rik[m]/nik);
        EXPECT_NEAR(f195, forces[water.hoh[0]][m] - forces0[water.hoh[0]][m],
                    gmx::test::precisionTolerance(1e-7, 2e-3)) << "dimension " << m;
        // No net force
        EXPECT_NEAR(0, fsum[m], gmx::test::precisionTolerance(1e-9, 2e-3)) << "dimension " << m;
    }
}

// Test: The energy is zero with the hydrogen atoms at their equilibrium
// distance, and it changes sign with the sum of the bond lengths
TEST(CrossBondAngleTest, SignOfTheEnergy)
{
    Water water;
    makeWater(cKrt, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    const double theta = 2*std::asin(0.5*cR3e/cR1e);
    // Factors for the bond length and for the angle
    const auto   energy = [&](double scale, double dtheta)
    {
        double r  = scale*cR1e;
        double th = theta + dtheta;
        water.x[water.hoh[1]] = { 0, 0, 0 };
        water.x[water.hoh[0]] = { static_cast<real>(r), 0, 0 };
        water.x[water.hoh[2]] = { static_cast<real>(r*std::cos(th)), static_cast<real>(r*std::sin(th)), 0 };
        std::vector<gmx::RVec>            forces;
        std::map<InteractionType, double> energies;
        compute(water, &forces, &energies);
        return energies[InteractionType::ANGLES];
    };
    // At the equilibrium lengths of the bonds any angle has no energy
    for (double dtheta : { -0.2, 0.0, 0.2 })
    {
        EXPECT_NEAR(0, energy(1.0, dtheta), gmx::test::precisionTolerance(1e-9, 1e-3)) << dtheta;
    }
    // Longer bonds at the equilibrium angle: both factors are positive
    // and the force constant is negative
    EXPECT_LT(energy(1.02, 0), -0.1);
    // Longer bonds with a distance between the hydrogen atoms that is
    // shorter than at equilibrium
    EXPECT_GT(energy(1.02, -0.2), 0.1);
}

// Test: The angle in a GROMACS topology has type 4 and the parameters in
// the order of GROMACS
TEST(CrossBondAngleTest, GromacsTopology)
{
    Water water;
    makeWater(cKrt, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    gmx::test::TestFileManager fileManager;
    const auto topName = fileManager.getTemporaryFilePath("water.top");
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        FILE *fp = fopen(topName.c_str(), "w");
        ASSERT_NE(nullptr, fp);
        write_top(&msghandler, fp, "water", water.mol.topology(), &water.pd);
        fclose(fp);
    }
    std::ifstream in(topName);
    ASSERT_TRUE(in.good());
    std::string line;
    bool        inAngles = false;
    int         nangle   = 0;
    while (std::getline(in, line))
    {
        if (line.find('[') != std::string::npos)
        {
            inAngles = line.find("angles") != std::string::npos;
            continue;
        }
        if (!inAngles || line.empty() || line[0] == ';')
        {
            continue;
        }
        std::istringstream iss(line);
        int                ai, aj, ak, type;
        double             r1e, r2e, r3e, krt;
        ASSERT_TRUE(static_cast<bool>(iss >> ai >> aj >> ak >> type >> r1e >> r2e >> r3e >> krt)) << line;
        EXPECT_EQ(water.hoh[0] + 1, ai);
        EXPECT_EQ(water.hoh[1] + 1, aj);
        EXPECT_EQ(water.hoh[2] + 1, ak);
        EXPECT_EQ(4, type);
        EXPECT_NEAR(cR1e, r1e, 1e-9);
        EXPECT_NEAR(cR2e, r2e, 1e-9);
        EXPECT_NEAR(cR3e, r3e, 1e-9);
        EXPECT_NEAR(cKrt, krt, 1e-3);
        nangle += 1;
    }
    EXPECT_EQ(1, nangle);
}

} // namespace

} // namespace alexandria
