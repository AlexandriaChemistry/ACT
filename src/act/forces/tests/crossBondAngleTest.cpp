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

#include <algorithm>
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
#include "act/forcefield/particletype.h"
#include "act/forcefield/potential.h"
#include "act/forces/forcecomputer.h"
#include "act/import/import.h"
#include "act/molprop/molprop.h"
#include "act/molprop/topologyentry.h"
#include "act/topology/actmol.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/* The terms of the flexible TIP4P model of water of Lawrence and Skinner,
 * Chem. Phys. Lett. 372 (2003) 842, that involve the angle, as in the
 * GROMACS topology tip4pflex.itp: the bond-angle cross term (angle of type
 * 4), the bond-bond cross term (angle of type 3) and the harmonic bond
 * between the hydrogen atoms.
 */
//! Equilibrium length of the first bond (nm)
const double cR1e = 0.09347;
//! Equilibrium length of the second bond (nm)
const double cR2e = 0.09347;
//! Equilibrium distance between the hydrogen atoms (nm)
const double cR3e = 0.1526;
//! Force constant of the bond-angle cross term (kJ/mol nm^2)
const double cKrt = -81901.5;
//! Force constant of the bond-bond cross term (kJ/mol nm^2)
const double cKrr = 52995.1;
//! Force constant of the harmonic potential between the hydrogen atoms (kJ/mol nm^2)
const double cKub = 103098.98;

/* A water molecule that is slightly off equilibrium (nm): the bonds are
 * 0.09522 and 0.09728 nm long and the hydrogen atoms 0.15209 nm apart.
 */
//! Position of the oxygen atom
const double cXO[3]  = { 0.001, -0.002, 0.0005 };
//! Position of the first hydrogen atom
const double cXH1[3] = { 0.096, 0.004, -0.002 };
//! Position of the second hydrogen atom
const double cXH2[3] = { -0.029, 0.0905, 0.003 };

/* Energies and forces for these parameters and positions from Eqns. 192 to
 * 195 of the GROMACS reference manual (bonded interactions, 2026.4) and
 * the harmonic potential, evaluated with 40 decimal digits in a Python
 * script.
 */
//! Energy of the bond-angle cross term (kJ/mol)
const double cEnergyBondAngle = 0.230815456837317;
//! Energy of the bond-bond cross term (kJ/mol)
const double cEnergyBondBond  = 0.353341846388098;
//! Energy of the harmonic potential between the hydrogen atoms (kJ/mol)
const double cEnergyHarmonic  = 0.0132560120902048;
//! Force of the bond-angle cross term on the first hydrogen atom (kJ/mol nm)
const double cForceBondAngleH1[3] = { 332.650776798099, -261.484635044086, -13.8730423173901 };
//! Force of the three terms on the first hydrogen atom (kJ/mol nm)
const double cForceH1[3] = { 174.423007331907, -303.925933187013, -10.2971429828126 };
//! Force of the three terms on the oxygen atom (kJ/mol nm)
const double cForceO[3]  = { 201.186796271703, 143.112428634474, -2.931305149574 };
//! Force of the three terms on the second hydrogen atom (kJ/mol nm)
const double cForceH2[3] = { -375.60980360361, 160.81350455254, 13.2284481323866 };

/* The same from GROMACS 2026.1 (mixed precision) for a topology with the
 * three terms as angles of type 4, 3 and 5 (Urey-Bradley without angle
 * part): the energy terms BA-Cross, Bond-Cross and U-B, and the forces as
 * gmx dump prints them. With tip4pflex.itp itself the energy terms are the
 * same, with Bond for the harmonic term.
 */
//! Energy term BA-Cross (kJ/mol)
const double cGromacsBondAngle  = 0.230822;
//! Energy term Bond-Cross (kJ/mol)
const double cGromacsBondBond   = 0.353341;
//! Energy term U-B (kJ/mol)
const double cGromacsHarmonic   = 0.0132568;
//! Force on the first hydrogen atom (kJ/mol nm)
const double cGromacsForceH1[3] = { 174.423, -303.927, -10.2971 };
//! Force on the oxygen atom (kJ/mol nm)
const double cGromacsForceO[3]  = { 201.187, 143.114, -2.93131 };
//! Force on the second hydrogen atom (kJ/mol nm)
const double cGromacsForceH2[3] = { -375.610, 160.813, 13.2285 };

/* Flexible TIP4P as a whole in GROMACS 2026.1, with tip4pflex.itp and the
 * virtual site at its place: the potential energy, with Morse bonds that
 * are zero at their minimum, and the forces.
 */
//! Potential energy (kJ/mol)
const double cGromacsPotential       = 5.11491;
//! Potential energy with 40 digits (kJ/mol)
const double cPotential              = 5.11489709058226;
//! Total force on the first hydrogen atom (kJ/mol nm)
const double cGromacsTotalForceH1[3] = { -743.802, -361.920, 13.8667 };
//! Total force on the oxygen atom (kJ/mol nm)
const double cGromacsTotalForceO[3]  = { 549.617, 1957.98, 20.3878 };
//! Total force on the second hydrogen atom (kJ/mol nm)
const double cGromacsTotalForceH2[3] = { 194.185, -1596.06, -34.2545 };

//! A water molecule with its force field
struct Water
{
    //! The force field
    ForceField             pd;
    //! The molecule
    ACTMol                 mol;
    //! The indices of the first hydrogen, the oxygen and the second hydrogen atom
    std::vector<int>       hoh;
    //! The positions, off equilibrium
    std::vector<gmx::RVec> x;
};

/*! \brief Make a water molecule
 * \param[in]  pd       The force field
 * \param[out] mol      The molecule with its topology
 * \param[in]  fileName The file in the test data, for another molecule
 */
void makeMolecule(ForceField *pd, ACTMol *mol, const std::string &fileName = "water.sdf")
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    std::vector<MolProp> mps;
    double               qtot = 0;
    importFile(&msghandler, pd, gmx::test::TestFileManager::getInputFilePath("molecules/" + fileName),
               &mps, "minimum", JobType::OPT, false, &qtot, true);
    ASSERT_TRUE(msghandler.ok());
    ASSERT_EQ(1u, mps.size());
    mol->Merge(&mps[0]);
    mol->GenerateTopology(&msghandler, pd, missingParameters::Ignore);
    ASSERT_TRUE(msghandler.ok());
}

/*! \brief Find the angle of a water molecule and put the atoms off equilibrium
 *
 * Other particles, such as a virtual site, are put on the oxygen atom: the
 * force computer constructs their positions.
 * \param[in,out] water The molecule with its topology
 */
void placeAtoms(Water *water)
{
    ASSERT_TRUE(water->mol.topology()->hasEntry(InteractionType::ANGLES));
    const auto &angles = water->mol.topology()->entry(InteractionType::ANGLES);
    ASSERT_EQ(1u, angles.size());
    water->hoh = angles[0]->atomIndices();
    ASSERT_EQ(3u, water->hoh.size());
    const auto &atoms = water->mol.topology()->atoms();
    ASSERT_EQ("O", atoms[water->hoh[1]].element());
    water->x.assign(atoms.size(), { static_cast<real>(cXO[XX]), static_cast<real>(cXO[YY]),
                                    static_cast<real>(cXO[ZZ]) });
    for (int m = 0; m < DIM; m++)
    {
        water->x[water->hoh[0]][m] = cXH1[m];
        water->x[water->hoh[2]][m] = cXH2[m];
    }
}

/*! \brief Make a water molecule with CROSS_BOND_ANGLES for its angle
 *
 * The force field ACS-g2 gets CROSS_BOND_ANGLES as the potential for the
 * angles, with one angle, that of water.
 * \param[in]  krt   The force constant of the bond-angle cross term (kJ/mol nm^2)
 * \param[in]  krr   The force constant of the bond-bond cross term (kJ/mol nm^2)
 * \param[in]  kub   The force constant of the harmonic potential (kJ/mol nm^2)
 * \param[out] water The molecule, off equilibrium
 */
void makeWater(double krt, double krr, double kub, Water *water)
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
    const std::vector<std::pair<const char *, double>> forceConstants = {
        { "krt", krt }, { "krr", krr }, { "kub", kub }
    };
    for (const auto &fc : forceConstants)
    {
        fs->addParameter(angleId, fc.first,
                         ForceFieldParameter("kJ/mol/nm2", fc.second, 0, 1, -200000, 200000,
                                             Mutability::Bounded, true, false));
    }
    makeMolecule(&water->pd, &water->mol);
    ASSERT_EQ(3u, water->mol.topology()->atoms().size());
    placeAtoms(water);
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

/*! \brief Energy of the angle for a set of force constants
 * \param[in] krt The force constant of the bond-angle cross term (kJ/mol nm^2)
 * \param[in] krr The force constant of the bond-bond cross term (kJ/mol nm^2)
 * \param[in] kub The force constant of the harmonic potential (kJ/mol nm^2)
 * \return the energy (kJ/mol)
 */
double angleEnergy(double krt, double krr, double kub)
{
    Water water;
    makeWater(krt, krr, kub, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return 0;
    }
    std::vector<gmx::RVec>            forces;
    std::map<InteractionType, double> energies;
    compute(water, &forces, &energies);
    return energies[InteractionType::ANGLES];
}

// Test: The names of the potential and of its parameters: those of the
// bond-angle cross term in the order of GROMACS, then the other two
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
    EXPECT_STREQ("krr", names[cbaKRR]);
    EXPECT_STREQ("kub", names[cbaKUB]);
}

// Test: Each of the three terms by itself, against the equations and
// against the energy terms of GROMACS
TEST(CrossBondAngleTest, Terms)
{
    // The positions are rounded to single precision in a single precision
    // build, which changes the small deviations from the equilibrium
    // lengths, as in GROMACS.
    const double etol = gmx::test::precisionTolerance(1e-10, 1e-4);
    const double gtol = gmx::test::precisionTolerance(2e-5, 1e-4);
    double       eba  = angleEnergy(cKrt, 0, 0);
    EXPECT_NEAR(cEnergyBondAngle, eba, etol);
    EXPECT_NEAR(cGromacsBondAngle, eba, gtol);
    double       ebb  = angleEnergy(0, cKrr, 0);
    EXPECT_NEAR(cEnergyBondBond, ebb, etol);
    EXPECT_NEAR(cGromacsBondBond, ebb, gtol);
    double       eub  = angleEnergy(0, 0, cKub);
    EXPECT_NEAR(cEnergyHarmonic, eub, etol);
    EXPECT_NEAR(cGromacsHarmonic, eub, gtol);
    EXPECT_EQ(0, angleEnergy(0, 0, 0));
}

// Test: Forces of the bond-angle cross term by itself, against Eqn. 195
// of the GROMACS manual
TEST(CrossBondAngleTest, BondAngleForces)
{
    Water water, without;
    makeWater(cKrt, 0, 0, &water);
    makeWater(0, 0, 0, &without);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    std::vector<gmx::RVec>            forces, forces0;
    std::map<InteractionType, double> energies, energies0;
    compute(water, &forces, &energies);
    compute(without, &forces0, &energies0);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    // The equation, written out here with the positions as the force
    // computer got them
    double rij[3], rik[3], nij2 = 0, nkj2 = 0, nik2 = 0;
    for (int m = 0; m < DIM; m++)
    {
        double xi = water.x[water.hoh[0]][m];
        double xj = water.x[water.hoh[1]][m];
        double xk = water.x[water.hoh[2]][m];
        rij[m]    = xi - xj;
        rik[m]    = xi - xk;
        nij2     += rij[m]*rij[m];
        nkj2     += (xk - xj)*(xk - xj);
        nik2     += rik[m]*rik[m];
    }
    double nij = std::sqrt(nij2), nkj = std::sqrt(nkj2), nik = std::sqrt(nik2);
    EXPECT_NEAR(cKrt*(nik - cR3e)*(nij - cR1e + nkj - cR2e), energies[InteractionType::ANGLES],
                gmx::test::precisionTolerance(1e-10, 1e-7));
    for (int m = 0; m < DIM; m++)
    {
        // The other forces of the molecule are the same in both
        double f    = forces[water.hoh[0]][m] - forces0[water.hoh[0]][m];
        double f195 = -cKrt*((nik - cR3e)*rij[m]/nij + (nij - cR1e + nkj - cR2e)*rik[m]/nik);
        EXPECT_NEAR(f195, f, gmx::test::precisionTolerance(1e-7, 2e-3)) << "dimension " << m;
        EXPECT_NEAR(cForceBondAngleH1[m], f, gmx::test::precisionTolerance(1e-7, 5e-2)) << "dimension " << m;
    }
}

// Test: Energy and forces of the three terms together for a water molecule
// that is off equilibrium, with the parameters of flexible TIP4P, against
// the equations and against GROMACS
TEST(CrossBondAngleTest, FlexibleTip4pWater)
{
    Water water;
    makeWater(cKrt, cKrr, cKub, &water);
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
    EXPECT_NEAR(cKrr, params[cbaKRR], 1e-6);
    EXPECT_NEAR(cKub, params[cbaKUB], 1e-6);

    std::vector<gmx::RVec>            forces;
    std::map<InteractionType, double> energies;
    compute(water, &forces, &energies);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    ASSERT_EQ(1u, energies.count(InteractionType::ANGLES));
    EXPECT_NEAR(cEnergyBondAngle + cEnergyBondBond + cEnergyHarmonic, energies[InteractionType::ANGLES],
                gmx::test::precisionTolerance(1e-10, 2e-4));
    // GROMACS has single precision positions
    EXPECT_NEAR(cGromacsBondAngle + cGromacsBondBond + cGromacsHarmonic, energies[InteractionType::ANGLES],
                gmx::test::precisionTolerance(2e-5, 2e-4));

    // The forces of the angle alone: the difference with the forces
    // without it
    Water without;
    makeWater(0, 0, 0, &without);
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
            EXPECT_NEAR(reference[a][m], f, gmx::test::precisionTolerance(1e-7, 5e-2))
                << "atom " << a << " dimension " << m;
            EXPECT_NEAR(gromacs[a][m], f, gmx::test::precisionTolerance(3e-3, 5e-2))
                << "GROMACS, atom " << a << " dimension " << m;
            fsum[m] += f;
        }
    }
    for (int m = 0; m < DIM; m++)
    {
        // No net force
        EXPECT_NEAR(0, fsum[m], gmx::test::precisionTolerance(1e-9, 2e-3)) << "dimension " << m;
    }
}

// Test: The force field file of flexible TIP4P as a whole, with Morse bonds
// and a virtual site, against GROMACS with tip4pflex.itp
TEST(CrossBondAngleTest, FlexibleTip4pForceField)
{
    Water water;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/TIP4P-Flex.xml"), &water.pd);
    EXPECT_EQ(Potential::MORSE_BONDS, water.pd.findForcesConst(InteractionType::BONDS).potential());
    EXPECT_EQ(Potential::CROSS_BOND_ANGLES, water.pd.findForcesConst(InteractionType::ANGLES).potential());
    makeMolecule(&water.pd, &water.mol);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    // Three atoms and the virtual site
    ASSERT_EQ(4u, water.mol.topology()->atoms().size());
    placeAtoms(&water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    std::vector<gmx::RVec>            forces;
    std::map<InteractionType, double> energies;
    compute(water, &forces, &energies);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    // The stiff Morse bonds make the energy sensitive to the rounding of
    // the positions in single precision
    EXPECT_NEAR(cEnergyBondAngle + cEnergyBondBond + cEnergyHarmonic, energies[InteractionType::ANGLES],
                gmx::test::precisionTolerance(1e-10, 2e-4));
    EXPECT_NEAR(cPotential, energies[InteractionType::EPOT], gmx::test::precisionTolerance(1e-9, 1e-3));
    EXPECT_NEAR(cGromacsPotential, energies[InteractionType::EPOT], gmx::test::precisionTolerance(3e-5, 1e-3));
    EXPECT_NEAR(cPotential - cEnergyBondAngle - cEnergyBondBond - cEnergyHarmonic,
                energies[InteractionType::BONDS], gmx::test::precisionTolerance(1e-9, 1e-3));
    // There are no interactions between the charges of one molecule
    if (energies.count(InteractionType::ELECTROSTATICS) > 0)
    {
        EXPECT_EQ(0, energies[InteractionType::ELECTROSTATICS]);
    }
    const std::vector<const double *> gromacs = { cGromacsTotalForceH1, cGromacsTotalForceO, cGromacsTotalForceH2 };
    for (size_t a = 0; a < gromacs.size(); a++)
    {
        for (int m = 0; m < DIM; m++)
        {
            // GROMACS prints six digits
            EXPECT_NEAR(gromacs[a][m], forces[water.hoh[a]][m], 6e-6*std::abs(gromacs[a][m]) + 5e-2)
                << "GROMACS, atom " << a << " dimension " << m;
        }
    }
}

// Test: The bond-angle cross term is zero with the bonds at their
// equilibrium lengths, and it changes sign with the distance between the
// hydrogen atoms
TEST(CrossBondAngleTest, SignOfTheBondAngleTerm)
{
    Water water;
    makeWater(cKrt, 0, 0, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    const double theta = 2*std::asin(0.5*cR3e/cR1e);
    // Factor for the bond length and change of the angle
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

/*! \brief Read the angles of a GROMACS topology
 * \param[in]  water The molecule
 * \param[out] lines For each angle the type and its parameters
 */
void readTopologyAngles(const Water &water, std::vector<std::pair<int, std::vector<double>>> *lines)
{
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
        std::istringstream iss(line.substr(0, line.find(';')));
        int                ai, aj, ak, type;
        ASSERT_TRUE(static_cast<bool>(iss >> ai >> aj >> ak >> type)) << line;
        EXPECT_EQ(water.hoh[0] + 1, ai);
        EXPECT_EQ(water.hoh[1] + 1, aj);
        EXPECT_EQ(water.hoh[2] + 1, ak);
        std::vector<double> params;
        double              value;
        while (iss >> value)
        {
            params.push_back(value);
        }
        lines->emplace_back(type, params);
    }
}

// Test: The angle becomes three angles in a GROMACS topology: the bond-angle
// cross term (type 4), the bond-bond cross term (type 3) and a Urey-Bradley
// angle (type 5) without angle part
TEST(CrossBondAngleTest, GromacsTopology)
{
    Water water;
    makeWater(cKrt, cKrr, cKub, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    std::vector<std::pair<int, std::vector<double>>> lines;
    readTopologyAngles(water, &lines);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    const std::vector<std::pair<int, std::vector<double>>> expected = {
        { 4, { cR1e, cR2e, cR3e, cKrt } },
        { 3, { cR1e, cR2e, cKrr } },
        { 5, { 0, 0, cR3e, cKub } }
    };
    ASSERT_EQ(expected.size(), lines.size());
    for (size_t i = 0; i < expected.size(); i++)
    {
        EXPECT_EQ(expected[i].first, lines[i].first);
        ASSERT_EQ(expected[i].second.size(), lines[i].second.size()) << "type " << expected[i].first;
        for (size_t k = 0; k < expected[i].second.size(); k++)
        {
            // Nine digits in the file, six would not do for the last
            // force constant
            EXPECT_NEAR(expected[i].second[k], lines[i].second[k], 1e-8*std::abs(expected[i].second[k]))
                << "type " << expected[i].first << " parameter " << k;
        }
    }
}

// Test: Terms without force constant are left out of a GROMACS topology
TEST(CrossBondAngleTest, GromacsTopologyBondAngleOnly)
{
    Water water;
    makeWater(cKrt, 0, 0, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    std::vector<std::pair<int, std::vector<double>>> lines;
    readTopologyAngles(water, &lines);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    ASSERT_EQ(1u, lines.size());
    EXPECT_EQ(4, lines[0].first);
    EXPECT_EQ(4u, lines[0].second.size());
}

/* The lengths r1e and r2e belong to the bonds in the order of the
 * identifier of the angle. For an angle such as C-O-H the atoms can be in
 * the reverse order in the molecule.
 */

//! The force constants of the angles of 1-butanol (kJ/mol nm^2)
const double cButanolKrt = -80000, cButanolKrr = 50000, cButanolKub = 30000;

/*! \brief Make 1-butanol with CROSS_BOND_ANGLES for its angles
 *
 * The lengths r1e and r2e of each angle are the bond lengths of the force
 * field ACS-g2, such that they differ in an angle with different bonds.
 * \param[out] pd  The force field
 * \param[out] mol The molecule
 */
void makeButanol(ForceField *pd, ACTMol *mol)
{
    const std::string fileName("1-butanol.pdb");
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACS-g2.xml"), pd);
    std::vector<Identifier> angleIds;
    {
        ACTMol harmonic;
        makeMolecule(pd, &harmonic, fileName);
        ASSERT_TRUE(harmonic.topology()->hasEntry(InteractionType::ANGLES));
        for (const auto &angle : harmonic.topology()->entry(InteractionType::ANGLES))
        {
            if (std::find(angleIds.begin(), angleIds.end(), angle->id()) == angleIds.end())
            {
                angleIds.push_back(angle->id());
            }
        }
    }
    const auto &bonds = pd->findForcesConst(InteractionType::BONDS);
    auto       *fs    = pd->findForces(InteractionType::ANGLES);
    fs->eraseParameter();
    fs->setFunction("CROSS_BOND_ANGLES");
    for (const auto &id : angleIds)
    {
        const auto atoms = id.atoms();
        const auto bo    = id.bondOrders();
        ASSERT_EQ(3u, atoms.size());
        const Identifier bij({ atoms[0], atoms[1] }, { bo[0] }, CanSwap::Yes);
        const Identifier bkj({ atoms[2], atoms[1] }, { bo[1] }, CanSwap::Yes);
        ASSERT_TRUE(bonds.parameterExists(bij)) << bij.id();
        ASSERT_TRUE(bonds.parameterExists(bkj)) << bkj.id();
        const auto &lij = bonds.findParameterTypeConst(bij, "bondlength");
        const auto &lkj = bonds.findParameterTypeConst(bkj, "bondlength");
        ASSERT_EQ("pm", lij.unit());
        const std::vector<std::pair<const char *, double>> lengths = {
            { "r1e", lij.value() }, { "r2e", lkj.value() }, { "r3e", 0.8*(lij.value() + lkj.value()) }
        };
        for (const auto &len : lengths)
        {
            fs->addParameter(id, len.first,
                             ForceFieldParameter("pm", len.second, 0, 1, 0.9*len.second, 1.1*len.second,
                                                 Mutability::Bounded, true, true));
        }
        const std::vector<std::pair<const char *, double>> forceConstants = {
            { "krt", cButanolKrt }, { "krr", cButanolKrr }, { "kub", cButanolKub }
        };
        for (const auto &fc : forceConstants)
        {
            fs->addParameter(id, fc.first,
                             ForceFieldParameter("kJ/mol/nm2", fc.second, 0, 1, -200000, 200000,
                                                 Mutability::Bounded, true, false));
        }
    }
    makeMolecule(pd, mol, fileName);
}

/*! \brief Tell from the atoms whether an angle is in the reverse order of
 * its identifier
 * \param[in]  pd      The force field
 * \param[in]  mol     The molecule
 * \param[in]  angle   The angle
 * \param[out] reverse Whether the first atom of the angle is the last one
 *                     of the identifier
 * \return whether the outer atoms differ, such that the order can be told
 */
bool angleOrder(const ForceField    &pd,
                const ACTMol        &mol,
                const TopologyEntry &angle,
                bool                *reverse)
{
    const auto &atoms = mol.topology()->atoms();
    const auto  bondType = [&](int index)
    {
        return pd.findParticleType(atoms[index].ffType())->interactionTypeToIdentifier(InteractionType::BONDS).id();
    };
    const auto idAtoms = angle.id().atoms();
    *reverse = false;
    if (idAtoms[0] == idAtoms[2])
    {
        return false;
    }
    const auto first = bondType(angle.atomIndex(0));
    const auto last  = bondType(angle.atomIndex(2));
    EXPECT_TRUE((first == idAtoms[0] && last == idAtoms[2]) ||
                (first == idAtoms[2] && last == idAtoms[0])) << angle.id().id();
    *reverse = (first == idAtoms[2]);
    return true;
}

// Test: The energy of angles with different bonds, of which the atoms are
// in the order of the identifier or in the reverse order, is that of the
// equations with r1e for the bond of the first atom of the identifier
TEST(CrossBondAngleTest, OrderOfTheAtoms)
{
    ForceField pd;
    ACTMol     mol;
    makeButanol(&pd, &mol);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    auto         x        = mol.xOriginal();
    const auto  &angles   = mol.topology()->entry(InteractionType::ANGLES);
    double       expected = 0;
    int          nForward = 0, nReverse = 0;
    const auto   distance = [&](int i, int j)
    {
        double d2 = 0;
        for (int m = 0; m < DIM; m++)
        {
            d2 += (x[i][m] - x[j][m])*(x[i][m] - x[j][m]);
        }
        return std::sqrt(d2);
    };
    for (const auto &angle : angles)
    {
        bool reverse = false;
        if (angleOrder(pd, mol, *angle->self(), &reverse))
        {
            (reverse ? nReverse : nForward)++;
        }
        EXPECT_EQ(reverse, angle->swapped()) << angle->id().id();
        const auto &params = angle->params();
        ASSERT_EQ(static_cast<size_t>(cbaNR), params.size());
        // The first atom of the identifier
        const int    ai = angle->atomIndex(reverse ? 2 : 0);
        const int    aj = angle->atomIndex(1);
        const int    ak = angle->atomIndex(reverse ? 0 : 2);
        const double d1 = distance(ai, aj) - params[cbaR1E];
        const double d2 = distance(ak, aj) - params[cbaR2E];
        const double d3 = distance(ai, ak) - params[cbaR3E];
        expected += params[cbaKRT]*d3*(d1 + d2) + params[cbaKRR]*d1*d2 + 0.5*params[cbaKUB]*d3*d3;
    }
    // Both orders are there, with different lengths of the bonds
    EXPECT_GT(nForward, 0);
    EXPECT_GT(nReverse, 0);

    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    const ForceComputer               forceComp;
    std::vector<gmx::RVec>            forces(x.size(), { 0, 0, 0 });
    std::map<InteractionType, double> energies;
    forceComp.compute(&msghandler, &pd, mol.topology(), &x, &forces, &energies);
    ASSERT_TRUE(msghandler.ok());
    ASSERT_EQ(1u, energies.count(InteractionType::ANGLES));
    EXPECT_NEAR(expected, energies[InteractionType::ANGLES],
                gmx::test::precisionTolerance(1e-9, 2e-3)*std::max(1.0, std::abs(expected)));
}

// Test: In a GROMACS topology the lengths of the bonds are in the order of
// the atoms of each angle
TEST(CrossBondAngleTest, GromacsTopologyOrderOfTheAtoms)
{
    ForceField pd;
    ACTMol     mol;
    makeButanol(&pd, &mol);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    gmx::test::TestFileManager fileManager;
    const auto topName = fileManager.getTemporaryFilePath("butanol.top");
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        FILE *fp = fopen(topName.c_str(), "w");
        ASSERT_NE(nullptr, fp);
        write_top(&msghandler, fp, "butanol", mol.topology(), &pd);
        fclose(fp);
    }
    const auto   &angles = mol.topology()->entry(InteractionType::ANGLES);
    std::ifstream in(topName);
    ASSERT_TRUE(in.good());
    std::string   line;
    bool          inAngles = false;
    int           nForward = 0, nReverse = 0;
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
        std::istringstream iss(line.substr(0, line.find(';')));
        int                ai, aj, ak, type;
        double             r1, r2;
        ASSERT_TRUE(static_cast<bool>(iss >> ai >> aj >> ak >> type >> r1 >> r2)) << line;
        if (4 != type && 3 != type)
        {
            // The Urey-Bradley angle has no lengths of the bonds
            continue;
        }
        bool found = false;
        for (const auto &angle : angles)
        {
            if (angle->atomIndex(0) == ai - 1 && angle->atomIndex(1) == aj - 1 && angle->atomIndex(2) == ak - 1)
            {
                bool reverse = false;
                if (angleOrder(pd, mol, *angle->self(), &reverse) && 4 == type)
                {
                    (reverse ? nReverse : nForward)++;
                }
                const auto &params = angle->params();
                EXPECT_NEAR(params[reverse ? cbaR2E : cbaR1E], r1, 1e-8) << line;
                EXPECT_NEAR(params[reverse ? cbaR1E : cbaR2E], r2, 1e-8) << line;
                found = true;
            }
        }
        EXPECT_TRUE(found) << line;
    }
    EXPECT_GT(nForward, 0);
    EXPECT_GT(nReverse, 0);
}

} // namespace

} // namespace alexandria
