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
 * Tests of the quartic angle and the quartic bond: polynomials of the
 * fourth order in the deviation from the reference angle or bond length.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/gromacs_top.h"
#include "act/alexandria/openmm_xml.h"
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

/* The quartic angle is Eqn. 196 of the GROMACS reference manual (bonded
 * interactions, 2026.4), the sum over n from 0 to 4 of C_n (theta-theta0)^n,
 * and the quartic bond is the same polynomial in r - bondlength.
 */
//! Reference angle (degrees)
const double cTheta0    = 95;
//! Coefficients of the quartic angle (kJ/mol/rad^n)
const double cAngleC[5] = { 1.5, -3, 200, -150, 400 };
//! Reference bond length (nm)
const double cBond0     = 0.09;
//! Coefficients of the quartic bond (kJ/mol/nm^n)
const double cBondC[5]  = { -250, 120, 150000, -3000000, 35000000 };
//! No coefficients
const double cZero[5]   = { 0, 0, 0, 0, 0 };

/* A water molecule that is off equilibrium (nm): the bonds are 0.09522
 * and 0.09728 nm long and the angle is 104.385 degrees.
 */
//! Position of the oxygen atom
const double cXO[3]  = { 0.001, -0.002, 0.0005 };
//! Position of the first hydrogen atom
const double cXH1[3] = { 0.096, 0.004, -0.002 };
//! Position of the second hydrogen atom
const double cXH2[3] = { -0.029, 0.0905, 0.003 };

/* Energies and forces for these parameters and positions from the
 * equation, evaluated with 40 decimal digits in a Python script. The
 * forces were checked there against numerical derivatives of the energy.
 */
//! Energy of the quartic angle (kJ/mol)
const double cAngleEnergy       = 6.00363582384278;
//! The five terms of the energy of the quartic angle (kJ/mol)
const double cAngleTerms[5]     = { 1.5, -0.491411625010449, 5.36634189323135,
                                    -0.659270697528617, 0.287976253150498 };
//! Force of the quartic angle on the first hydrogen atom (kJ/mol nm)
const double cAngleForceH1[3]   = { -37.7278613479851, 602.337374849742, 11.9509684159492 };
//! Force of the quartic angle on the oxygen atom (kJ/mol nm)
const double cAngleForceO[3]    = { -524.13012548484, -784.889174742313, 0.169785615279163 };
//! Force of the quartic angle on the second hydrogen atom (kJ/mol nm)
const double cAngleForceH2[3]   = { 561.857986832825, 182.551799892571, -12.1207540312284 };
//! Energy of the two quartic bonds (kJ/mol)
const double cBondEnergy        = -487.928482583733;
//! The five terms of the energy of the two quartic bonds (kJ/mol)
const double cBondTerms[5]      = { -500, 1.49969892828478, 12.0302411359242,
                                    -1.58251117266449, 0.124088524722013 };
//! Force of the quartic bonds on the first hydrogen atom (kJ/mol nm)
const double cBondForceH1[3]    = { -1457.72797721759, -92.0670301400583, 38.3612625583576 };
//! Force of the quartic bonds on the oxygen atom (kJ/mol nm)
const double cBondForceO[3]     = { 877.885306897708, 1879.91526362636, 9.95895996829924 };
//! Force of the quartic bonds on the second hydrogen atom (kJ/mol nm)
const double cBondForceH2[3]    = { 579.842670319883, -1787.8482334863, -48.3202225266569 };

/* The same from GROMACS 2026.1 (mixed precision). For the angle a topology
 * with an angle of type 6 with the parameters above: the energy term
 * Quartic Angles and the forces as gmx dump prints them. GROMACS has no
 * quartic bond: its cubic bond (bond of type 4), kb (r-b0)^2 +
 * kb kcub (r-b0)^3, is the quartic bond with the second and third
 * coefficient only. With b0 = 0.09 nm, kb = 150000 kJ/mol nm^2 and
 * kcub = -20 /nm for both bonds: the energy term Cubic Bonds and the forces.
 */
//! Energy term Quartic Angles (kJ/mol)
const double cGromacsAngleEnergy     = 6.00363;
//! Force on the first hydrogen atom (kJ/mol nm)
const double cGromacsAngleForceH1[3] = { -37.7279, 602.337, 11.9509 };
//! Force on the oxygen atom (kJ/mol nm)
const double cGromacsAngleForceO[3]  = { -524.130, -784.889, 0.169847 };
//! Force on the second hydrogen atom (kJ/mol nm)
const double cGromacsAngleForceH2[3] = { 561.858, 182.552, -12.1208 };
//! Coefficients of the quartic bond that is the cubic bond of GROMACS (kJ/mol/nm^n)
const double cCubicC[5]              = { 0, 0, 150000, -3000000, 0 };
//! Energy term Cubic Bonds (kJ/mol)
const double cGromacsCubicEnergy     = 10.4477;
//! The same with 40 digits (kJ/mol)
const double cCubicEnergy            = 10.4477299632597;
//! Force of the cubic bonds on the first hydrogen atom (kJ/mol nm)
const double cGromacsCubicForceH1[3] = { -1318.11, -83.2494, 34.6873 };
//! Force of the cubic bonds on the oxygen atom (kJ/mol nm)
const double cGromacsCubicForceO[3]  = { 791.907, 1705.72, 9.16328 };
//! Force of the cubic bonds on the second hydrogen atom (kJ/mol nm)
const double cGromacsCubicForceH2[3] = { 526.208, -1622.47, -43.8506 };

//! A parameter of an interaction: its name and unit in a force field and its value
struct Parameter
{
    //! The name
    const char *name;
    //! The unit
    const char *unit;
    //! The value
    double      value;
};

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

/*! \brief Give the bonds or the angles of the force field another function
 *
 * All parameters of the interaction type are removed, and those of the
 * bonds or the angle of water are added.
 * \param[in,out] pd       The force field
 * \param[in]     itype    Bonds or angles
 * \param[in]     function The name of the potential
 * \param[in]     params   Its parameters, the same for all interactions
 */
void setFunction(ForceField                   *pd,
                 InteractionType               itype,
                 const char                   *function,
                 const std::vector<Parameter> &params)
{
    std::vector<Identifier> ids;
    {
        // The identifiers of the interactions, from the force field as it is
        ACTMol before;
        makeMolecule(pd, &before);
        ASSERT_TRUE(before.topology()->hasEntry(itype));
        for (const auto &entry : before.topology()->entry(itype))
        {
            bool found = false;
            for (const auto &id : ids)
            {
                found = found || id.id() == entry->id().id();
            }
            if (!found)
            {
                ids.push_back(entry->id());
            }
        }
    }
    ASSERT_FALSE(ids.empty());
    auto *fs = pd->findForces(itype);
    fs->eraseParameter();
    fs->setFunction(function);
    for (const auto &id : ids)
    {
        for (const auto &p : params)
        {
            const double range = std::abs(p.value) + 1;
            fs->addParameter(id, p.name,
                             ForceFieldParameter(p.unit, p.value, 0, 1, p.value - range, p.value + range,
                                                 Mutability::Bounded, true, false));
        }
    }
}

/*! \brief The parameters of a quartic bond as they are in a force field
 * \param[in] length The reference bond length (nm)
 * \param[in] c      The five coefficients (kJ/mol/nm^n)
 * \return the parameters, with the length in pm
 */
std::vector<Parameter> quarticBond(double length, const double c[5])
{
    return { { "bondlength", "pm", 1000*length },
             { "cb0", "kJ/mol", c[0] },
             { "cb1", "kJ/mol/nm", c[1] },
             { "cb2", "kJ/mol/nm2", c[2] },
             { "cb3", "kJ/mol/nm3", c[3] },
             { "cb4", "kJ/mol/nm4", c[4] } };
}

/*! \brief The parameters of a quartic angle as they are in a force field
 * \param[in] angle The reference angle (degrees)
 * \param[in] c     The five coefficients (kJ/mol/rad^n)
 * \return the parameters
 */
std::vector<Parameter> quarticAngle(double angle, const double c[5])
{
    return { { "angle", "degree", angle },
             { "ct0", "kJ/mol", c[0] },
             { "ct1", "kJ/mol/rad", c[1] },
             { "ct2", "kJ/mol/rad2", c[2] },
             { "ct3", "kJ/mol/rad3", c[3] },
             { "ct4", "kJ/mol/rad4", c[4] } };
}

/*! \brief Make a water molecule that is off equilibrium
 *
 * The force field is ACS-g2, of which the bonds or the angles or both get
 * another function.
 * \param[in]  bondFunction  The potential for the bonds, nullptr for that of the force field
 * \param[in]  bondParams    Its parameters
 * \param[in]  angleFunction The potential for the angles, nullptr for that of the force field
 * \param[in]  angleParams   Its parameters
 * \param[out] water         The molecule
 */
void makeWater(const char                   *bondFunction,
               const std::vector<Parameter> &bondParams,
               const char                   *angleFunction,
               const std::vector<Parameter> &angleParams,
               Water                        *water)
{
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACS-g2.xml"), &water->pd);
    if (bondFunction)
    {
        setFunction(&water->pd, InteractionType::BONDS, bondFunction, bondParams);
        ASSERT_FALSE(::testing::Test::HasFatalFailure());
    }
    if (angleFunction)
    {
        setFunction(&water->pd, InteractionType::ANGLES, angleFunction, angleParams);
        ASSERT_FALSE(::testing::Test::HasFatalFailure());
    }
    makeMolecule(&water->pd, &water->mol);
    ASSERT_FALSE(::testing::Test::HasFatalFailure());
    ASSERT_EQ(3u, water->mol.topology()->atoms().size());
    ASSERT_TRUE(water->mol.topology()->hasEntry(InteractionType::ANGLES));
    const auto &angles = water->mol.topology()->entry(InteractionType::ANGLES);
    ASSERT_EQ(1u, angles.size());
    water->hoh = angles[0]->atomIndices();
    ASSERT_EQ(3u, water->hoh.size());
    ASSERT_EQ("O", water->mol.topology()->atoms()[water->hoh[1]].element());
    water->x.resize(3);
    for (int m = 0; m < DIM; m++)
    {
        water->x[water->hoh[0]][m] = cXH1[m];
        water->x[water->hoh[1]][m] = cXO[m];
        water->x[water->hoh[2]][m] = cXH2[m];
    }
}

//! Make a water molecule with quartic bonds, see makeWater
void makeQuarticBonds(double length, const double c[5], Water *water)
{
    makeWater("QUARTIC_BONDS", quarticBond(length, c), nullptr, {}, water);
}

//! Make a water molecule with a quartic angle, see makeWater
void makeQuarticAngle(double angle, const double c[5], Water *water)
{
    makeWater(nullptr, {}, "QUARTIC_ANGLES", quarticAngle(angle, c), water);
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

/*! \brief Energy of one interaction type of a water molecule
 * \param[in] water The molecule
 * \param[in] itype Bonds or angles
 * \return the energy (kJ/mol), zero after a failure
 */
double energyOf(const Water &water, InteractionType itype)
{
    if (::testing::Test::HasFatalFailure())
    {
        return 0;
    }
    std::vector<gmx::RVec>            forces;
    std::map<InteractionType, double> energies;
    compute(water, &forces, &energies);
    return energies[itype];
}

/*! \brief Compare the forces of one interaction type with reference values
 *
 * The forces are the difference between those of two molecules that
 * differ in that interaction type only.
 * \param[in] with      The molecule with the interactions
 * \param[in] without   The molecule in which they have no energy
 * \param[in] reference The forces on the first hydrogen, the oxygen and the second hydrogen atom
 * \param[in] tolerance The absolute tolerance
 * \param[in] what      Text for the messages
 */
void checkForces(const Water                       &with,
                 const Water                       &without,
                 const std::vector<const double *> &reference,
                 double                             tolerance,
                 const char                        *what)
{
    std::vector<gmx::RVec>            forces, forces0;
    std::map<InteractionType, double> energies, energies0;
    compute(with, &forces, &energies);
    compute(without, &forces0, &energies0);
    ASSERT_FALSE(::testing::Test::HasFatalFailure());
    ASSERT_EQ(3u, reference.size());
    for (size_t a = 0; a < reference.size(); a++)
    {
        for (int m = 0; m < DIM; m++)
        {
            double f = forces[with.hoh[a]][m] - forces0[without.hoh[a]][m];
            EXPECT_NEAR(reference[a][m], f, tolerance) << what << ", atom " << a << " dimension " << m;
        }
    }
}

// Test: The names of the potentials and of their parameters: the reference
// value first, then the coefficients, as the quartic angle in GROMACS
TEST(QuarticTest, Names)
{
    Potential p = Potential::NONE;
    ASSERT_TRUE(stringToPotential("QUARTIC_ANGLES", &p));
    EXPECT_EQ(Potential::QUARTIC_ANGLES, p);
    EXPECT_EQ("QUARTIC_ANGLES", potentialToString(Potential::QUARTIC_ANGLES));
    ASSERT_TRUE(stringToPotential("QUARTIC_BONDS", &p));
    EXPECT_EQ(Potential::QUARTIC_BONDS, p);
    EXPECT_EQ("QUARTIC_BONDS", potentialToString(Potential::QUARTIC_BONDS));
    // The name in GROMACS
    p = Potential::NONE;
    ASSERT_TRUE(stringToPotential("QANGLES", &p));
    EXPECT_EQ(Potential::QUARTIC_ANGLES, p);
    EXPECT_EQ(F_QUARTIC_ANGLES, potentialToGromacsType(Potential::QUARTIC_ANGLES));
    EXPECT_EQ(F_QUARTIC_BONDS, potentialToGromacsType(Potential::QUARTIC_BONDS));
    EXPECT_EQ(2, interaction_function[F_QUARTIC_BONDS].nratoms);
    EXPECT_EQ(3, interaction_function[F_QUARTIC_ANGLES].nratoms);

    const auto angle = potentialToParameterName(Potential::QUARTIC_ANGLES);
    ASSERT_EQ(static_cast<size_t>(qaNR), angle.size());
    EXPECT_STREQ("angle", angle[qaANGLE]);
    EXPECT_STREQ("ct0", angle[qaC0]);
    EXPECT_STREQ("ct1", angle[qaC1]);
    EXPECT_STREQ("ct2", angle[qaC2]);
    EXPECT_STREQ("ct3", angle[qaC3]);
    EXPECT_STREQ("ct4", angle[qaC4]);
    const auto bond = potentialToParameterName(Potential::QUARTIC_BONDS);
    ASSERT_EQ(static_cast<size_t>(qbNR), bond.size());
    // Virtual sites and linear angles look up the bond length by this name
    EXPECT_STREQ("bondlength", bond[qbLENGTH]);
    EXPECT_STREQ("cb0", bond[qbC0]);
    EXPECT_STREQ("cb1", bond[qbC1]);
    EXPECT_STREQ("cb2", bond[qbC2]);
    EXPECT_STREQ("cb3", bond[qbC3]);
    EXPECT_STREQ("cb4", bond[qbC4]);
    // The coefficients are in the order of the powers
    EXPECT_EQ(qaC0 + 4, qaC4);
    EXPECT_EQ(qbC0 + 4, qbC4);
}

// Test: Each of the five terms of the quartic angle by itself and their
// sum, against the equation and against GROMACS
TEST(QuarticTest, AngleEnergy)
{
    // The positions are rounded to single precision in a single precision build
    const double etol = gmx::test::precisionTolerance(1e-10, 1e-4);
    for (int n = 0; n < 5; n++)
    {
        double c[5] = { 0, 0, 0, 0, 0 };
        c[n]        = cAngleC[n];
        Water water;
        makeQuarticAngle(cTheta0, c, &water);
        EXPECT_NEAR(cAngleTerms[n], energyOf(water, InteractionType::ANGLES), etol) << "term " << n;
    }
    Water water;
    makeQuarticAngle(cTheta0, cAngleC, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    // The parameters arrive in the topology in radians and kJ/mol/rad^n
    const auto &params = water.mol.topology()->entry(InteractionType::ANGLES)[0]->params();
    ASSERT_EQ(static_cast<size_t>(qaNR), params.size());
    EXPECT_NEAR(cTheta0*M_PI/180, params[qaANGLE], 1e-12);
    for (int n = 0; n < 5; n++)
    {
        EXPECT_NEAR(cAngleC[n], params[qaC0 + n], 1e-9) << "coefficient " << n;
    }
    const double energy = energyOf(water, InteractionType::ANGLES);
    EXPECT_NEAR(cAngleEnergy, energy, etol);
    // GROMACS prints six digits
    EXPECT_NEAR(cGromacsAngleEnergy, energy, gmx::test::precisionTolerance(1e-5, 1e-4));
    // At the reference angle only the constant is left
    Water atReference;
    makeQuarticAngle(104.385270705588, cAngleC, &atReference);
    EXPECT_NEAR(cAngleC[0], energyOf(atReference, InteractionType::ANGLES), gmx::test::precisionTolerance(1e-10, 1e-5));
}

// Test: The forces of the quartic angle, against the equation and against
// GROMACS
TEST(QuarticTest, AngleForces)
{
    Water water, without;
    makeQuarticAngle(cTheta0, cAngleC, &water);
    makeQuarticAngle(cTheta0, cZero, &without);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    EXPECT_EQ(0, energyOf(without, InteractionType::ANGLES));
    checkForces(water, without, { cAngleForceH1, cAngleForceO, cAngleForceH2 },
                gmx::test::precisionTolerance(1e-8, 5e-2), "equation");
    // GROMACS has single precision positions and prints six digits
    checkForces(water, without, { cGromacsAngleForceH1, cGromacsAngleForceO, cGromacsAngleForceH2 },
                gmx::test::precisionTolerance(5e-3, 5e-2), "GROMACS");
}

// Test: With the quadratic term only the quartic angle is the harmonic
// angle with a force constant of twice the coefficient
TEST(QuarticTest, AngleIsHarmonic)
{
    const double c[5] = { 0, 0, cAngleC[2], 0, 0 };
    Water        quartic, harmonic, without;
    makeQuarticAngle(cTheta0, c, &quartic);
    makeQuarticAngle(cTheta0, cZero, &without);
    makeWater(nullptr, {}, "HARMONIC_ANGLES",
              { { "kt", "kJ/mol/rad2", 2*cAngleC[2] }, { "angle", "degree", cTheta0 } }, &harmonic);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    EXPECT_EQ(Potential::HARMONIC_ANGLES, harmonic.pd.findForcesConst(InteractionType::ANGLES).potential());
    const double tol = gmx::test::precisionTolerance(1e-10, 1e-5);
    EXPECT_NEAR(energyOf(harmonic, InteractionType::ANGLES), energyOf(quartic, InteractionType::ANGLES), tol);
    EXPECT_GT(energyOf(quartic, InteractionType::ANGLES), 5);
    std::vector<gmx::RVec>            fq, fh;
    std::map<InteractionType, double> eq, eh;
    compute(quartic, &fq, &eq);
    compute(harmonic, &fh, &eh);
    for (size_t a = 0; a < fq.size(); a++)
    {
        for (int m = 0; m < DIM; m++)
        {
            EXPECT_NEAR(fh[a][m], fq[a][m], gmx::test::precisionTolerance(1e-8, 2e-2)) << "atom " << a << " dimension " << m;
        }
    }
}

// Test: Each of the five terms of the quartic bond by itself and their
// sum for the two bonds of water, against the equation
TEST(QuarticTest, BondEnergy)
{
    // The bonds are stiff: the energy is sensitive to the rounding of the
    // positions in single precision
    const double etol = gmx::test::precisionTolerance(1e-9, 2e-3);
    for (int n = 0; n < 5; n++)
    {
        double c[5] = { 0, 0, 0, 0, 0 };
        c[n]        = cBondC[n];
        Water water;
        makeQuarticBonds(cBond0, c, &water);
        EXPECT_NEAR(cBondTerms[n], energyOf(water, InteractionType::BONDS), etol) << "term " << n;
    }
    Water water;
    makeQuarticBonds(cBond0, cBondC, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    // The parameters arrive in the topology in nm and kJ/mol/nm^n
    const auto &bonds = water.mol.topology()->entry(InteractionType::BONDS);
    ASSERT_EQ(2u, bonds.size());
    for (const auto &b : bonds)
    {
        const auto &params = b->params();
        ASSERT_EQ(static_cast<size_t>(qbNR), params.size());
        EXPECT_NEAR(cBond0, params[qbLENGTH], 1e-12);
        for (int n = 0; n < 5; n++)
        {
            EXPECT_NEAR(cBondC[n], params[qbC0 + n], 1e-9*(1 + std::abs(cBondC[n]))) << "coefficient " << n;
        }
    }
    EXPECT_NEAR(cBondEnergy, energyOf(water, InteractionType::BONDS), etol);
}

// Test: The forces of the quartic bonds, against the equation
TEST(QuarticTest, BondForces)
{
    Water water, without;
    makeQuarticBonds(cBond0, cBondC, &water);
    makeQuarticBonds(cBond0, cZero, &without);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    EXPECT_EQ(0, energyOf(without, InteractionType::BONDS));
    checkForces(water, without, { cBondForceH1, cBondForceO, cBondForceH2 },
                gmx::test::precisionTolerance(1e-7, 1), "equation");
}

// Test: With the quadratic and cubic terms only the quartic bond is the
// cubic bond of GROMACS
TEST(QuarticTest, BondIsGromacsCubic)
{
    Water water, without;
    makeQuarticBonds(cBond0, cCubicC, &water);
    makeQuarticBonds(cBond0, cZero, &without);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    const double energy = energyOf(water, InteractionType::BONDS);
    EXPECT_NEAR(cCubicEnergy, energy, gmx::test::precisionTolerance(1e-9, 2e-3));
    // GROMACS has single precision positions and prints six digits
    EXPECT_NEAR(cGromacsCubicEnergy, energy, gmx::test::precisionTolerance(1e-4, 2e-3));
    checkForces(water, without, { cGromacsCubicForceH1, cGromacsCubicForceO, cGromacsCubicForceH2 },
                gmx::test::precisionTolerance(2e-2, 1), "GROMACS");
}

// Test: With the constant and the quadratic term only the quartic bond is
// the harmonic bond with a force constant of twice the coefficient
TEST(QuarticTest, BondIsHarmonic)
{
    const double c[5] = { cBondC[0], 0, cBondC[2], 0, 0 };
    Water        quartic, harmonic;
    makeQuarticBonds(cBond0, c, &quartic);
    makeWater("HARMONIC_BONDS",
              { { "kb", "kJ/mol/nm2", 2*cBondC[2] }, { "bondlength", "pm", 1000*cBond0 },
                { "bondenergy", "kJ/mol", cBondC[0] } }, nullptr, {}, &harmonic);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    EXPECT_EQ(Potential::HARMONIC_BONDS, harmonic.pd.findForcesConst(InteractionType::BONDS).potential());
    EXPECT_NEAR(energyOf(harmonic, InteractionType::BONDS), energyOf(quartic, InteractionType::BONDS),
                gmx::test::precisionTolerance(1e-9, 1e-4));
    std::vector<gmx::RVec>            fq, fh;
    std::map<InteractionType, double> eq, eh;
    compute(quartic, &fq, &eq);
    compute(harmonic, &fh, &eh);
    for (size_t a = 0; a < fq.size(); a++)
    {
        for (int m = 0; m < DIM; m++)
        {
            EXPECT_NEAR(fh[a][m], fq[a][m], gmx::test::precisionTolerance(1e-7, 1e-1)) << "atom " << a << " dimension " << m;
        }
    }
}

// Test: The quartic bond with the coefficients that geometry_ff gives it,
// the expansion of a Morse potential, follows that Morse potential close
// to the minimum much better than the harmonic bond does
TEST(QuarticTest, BondFollowsMorse)
{
    const double De   = 400;
    const double beta = 20;
    // (1 - exp(-x))^2 = x^2 - x^3 + 7 x^4/12 - x^5/4 + ..., with x = beta (r - bondlength)
    const double c2          = De*beta*beta;
    const double expansion[5] = { 0, 0, c2, -c2*beta, c2*beta*beta*7.0/12.0 };
    const double parabola[5]  = { 0, 0, c2, 0, 0 };
    Water        quartic, harmonic, morse;
    makeQuarticBonds(cBond0, expansion, &quartic);
    makeQuarticBonds(cBond0, parabola, &harmonic);
    // A Morse potential that is zero at its minimum
    makeWater("MORSE_BONDS",
              { { "beta", "1/nm", beta }, { "De", "kJ/mol", De }, { "D0", "kJ/mol", De },
                { "bondlength", "pm", 1000*cBond0 } }, nullptr, {}, &morse);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    const double em = energyOf(morse, InteractionType::BONDS);
    const double eq = energyOf(quartic, InteractionType::BONDS);
    const double eh = energyOf(harmonic, InteractionType::BONDS);
    // The bonds are 0.0052 and 0.0073 nm longer than the reference: x is
    // 0.10 and 0.15, and the first term that is left out is De x^5/4
    EXPECT_GT(em, 10);
    EXPECT_NEAR(em, eq, 0.01);
    EXPECT_GT(std::abs(em - eh), 50*std::abs(em - eq));
}

/*! \brief Write a GROMACS topology and read the lines of one directive
 * \param[in]  water     The molecule
 * \param[in]  directive The name of the directive, bonds or angles
 * \param[in]  natoms    The number of atoms of an interaction
 * \param[out] lines     For each interaction the atoms, the function type and the parameters
 * \param[out] text      The topology as a whole, with its header
 * \param[out] status    What the writer reported
 */
void readTopology(const Water                                      &water,
                  const char                                       *directive,
                  int                                               natoms,
                  std::vector<std::pair<std::vector<int>, std::vector<double>>> *lines,
                  std::string                                      *text,
                  ACTStatus                                        *status)
{
    gmx::test::TestFileManager fileManager;
    const auto topName = fileManager.getTemporaryFilePath("water.top");
    {
        // Warnings are counted but not printed
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Error);
        FILE *fp = fopen(topName.c_str(), "w");
        ASSERT_NE(nullptr, fp);
        print_top_header(&msghandler, fp, &water.pd, water.mol.topology(), {}, false);
        write_top(&msghandler, fp, "water", water.mol.topology(), &water.pd);
        fclose(fp);
        *status = msghandler.status();
    }
    std::ifstream in(topName);
    ASSERT_TRUE(in.good());
    std::string line;
    bool        inDirective = false;
    while (std::getline(in, line))
    {
        *text += line + "\n";
        if (line.find('[') != std::string::npos)
        {
            inDirective = line.find(directive) != std::string::npos;
            continue;
        }
        if (!inDirective || line.empty() || line[0] == ';')
        {
            continue;
        }
        std::istringstream iss(line.substr(0, line.find(';')));
        std::vector<int>   atoms(natoms + 1);
        for (auto &a : atoms)
        {
            ASSERT_TRUE(static_cast<bool>(iss >> a)) << line;
        }
        std::vector<double> params;
        double              value;
        while (iss >> value)
        {
            params.push_back(value);
        }
        lines->emplace_back(atoms, params);
    }
}

// Test: The quartic angle is an angle of type 6 in a GROMACS topology,
// with the angle in degrees and the five coefficients
TEST(QuarticTest, AngleInGromacsTopology)
{
    Water water;
    makeQuarticAngle(cTheta0, cAngleC, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    std::vector<std::pair<std::vector<int>, std::vector<double>>> lines;
    std::string                                                   text;
    ACTStatus                                                     status = ACTStatus::Fatal;
    readTopology(water, "angles", 3, &lines, &text, &status);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    // GROMACS has the potential: it is not in the list of what is missing
    EXPECT_EQ(std::string::npos, text.find("QUARTIC_ANGLES")) << text;
    EXPECT_EQ(std::string::npos, text.find("unsupported function")) << text;
    ASSERT_EQ(1u, lines.size());
    const auto &atoms = lines[0].first;
    EXPECT_EQ(water.hoh[0] + 1, atoms[0]);
    EXPECT_EQ(water.hoh[1] + 1, atoms[1]);
    EXPECT_EQ(water.hoh[2] + 1, atoms[2]);
    EXPECT_EQ(6, atoms[3]);
    const auto &params = lines[0].second;
    ASSERT_EQ(6u, params.size());
    EXPECT_NEAR(cTheta0, params[0], 1e-7);
    for (int n = 0; n < 5; n++)
    {
        EXPECT_NEAR(cAngleC[n], params[1 + n], 1e-8*std::abs(cAngleC[n])) << "coefficient " << n;
    }
}

// Test: GROMACS has no quartic bond: the bonds are in the topology without
// parameters, such that GROMACS does not accept it, the header of the
// topology says so and there is a warning
TEST(QuarticTest, BondInGromacsTopology)
{
    Water water;
    makeQuarticBonds(cBond0, cBondC, &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    std::vector<std::pair<std::vector<int>, std::vector<double>>> lines;
    std::string                                                   text;
    ACTStatus                                                     status = ACTStatus::Fatal;
    readTopology(water, "bonds", 2, &lines, &text, &status);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    ASSERT_EQ(2u, lines.size());
    for (const auto &line : lines)
    {
        EXPECT_TRUE(line.second.empty());
    }
    EXPECT_NE(std::string::npos, text.find("QUARTIC_BONDS, GROMACS will not accept this topology")) << text;
    EXPECT_NE(std::string::npos, text.find("unsupported function QBONDS")) << text;
    EXPECT_EQ(ACTStatus::Warning, status);
}

// Test: In a force field file for OpenMM the bonds and the angle are
// custom forces with the polynomial as energy, and the reference values
// and coefficients are parameters of each bond and angle, in nm, radians
// and kJ/mol
TEST(QuarticTest, OpenMM)
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath("forcefields/ACS-g2.xml"), &pd);
    setFunction(&pd, InteractionType::BONDS, "QUARTIC_BONDS", quarticBond(cBond0, cBondC));
    setFunction(&pd, InteractionType::ANGLES, "QUARTIC_ANGLES", quarticAngle(cTheta0, cAngleC));
    std::vector<ACTMol> mols(1);
    makeMolecule(&pd, &mols[0]);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    gmx::test::TestFileManager fileManager;
    const auto xmlName = fileManager.getTemporaryFilePath("openmm.xml");
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        writeOpenMM(&msghandler, xmlName, &pd, mols, 0.1, 1, true);
        EXPECT_TRUE(msghandler.ok());
    }
    std::ifstream in(xmlName);
    ASSERT_TRUE(in.good());
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    for (const char *expected : {
             "<CustomBondForce energy=\"(cb0+dr*(cb1+dr*(cb2+dr*(cb3+dr*cb4))));dr=r-bondlength\">",
             "<GlobalParameter name=\"pot-QUARTIC-BONDS\" defaultValue=\"1\"/>",
             "<PerBondParameter name=\"bondlength\"/>",
             "<PerBondParameter name=\"cb0\"/>",
             "<PerBondParameter name=\"cb4\"/>",
             "bondlength=\"0.09\" cb0=\"-250\" cb1=\"120\" cb2=\"150000\" cb3=\"-3000000\" cb4=\"35000000\"/>",
             "<CustomAngleForce energy=\"(ct0+dt*(ct1+dt*(ct2+dt*(ct3+dt*ct4))));dt=theta-angle\">",
             "<PerAngleParameter name=\"angle\"/>",
             "<PerAngleParameter name=\"ct0\"/>",
             "<PerAngleParameter name=\"ct4\"/>",
             "angle=\"1.658062789395\" ct0=\"1.5\" ct1=\"-3\" ct2=\"200\" ct3=\"-150\" ct4=\"400\"/>" })
    {
        EXPECT_NE(std::string::npos, text.find(expected)) << expected << " is not in\n" << text;
    }
}

// Test: A force field with the two potentials is the same after writing it
// to a file and reading it again
TEST(QuarticTest, ForceFieldFile)
{
    Water water;
    makeWater("QUARTIC_BONDS", quarticBond(cBond0, cBondC),
              "QUARTIC_ANGLES", quarticAngle(cTheta0, cAngleC), &water);
    if (::testing::Test::HasFatalFailure())
    {
        return;
    }
    const double energy = energyOf(water, InteractionType::BONDS) + energyOf(water, InteractionType::ANGLES);
    EXPECT_NEAR(cBondEnergy + cAngleEnergy, energy, gmx::test::precisionTolerance(1e-9, 2e-3));

    gmx::test::TestFileManager fileManager;
    const auto ffName = fileManager.getTemporaryFilePath("quartic.xml");
    writeForceField(ffName, &water.pd, false);
    ForceField again;
    readForceField(ffName, &again);
    for (auto itype : { InteractionType::BONDS, InteractionType::ANGLES })
    {
        const auto &before = water.pd.findForcesConst(itype);
        const auto &after  = again.findForcesConst(itype);
        EXPECT_EQ(before.potential(), after.potential());
        ASSERT_EQ(1u, before.parametersConst().size());
        ASSERT_EQ(before.parametersConst().size(), after.parametersConst().size());
        const auto &pb = before.parametersConst().begin()->second;
        const auto &pa = after.parametersConst().begin()->second;
        ASSERT_EQ(6u, pb.size());
        ASSERT_EQ(pb.size(), pa.size());
        for (const auto &p : pb)
        {
            const auto found = pa.find(p.first);
            ASSERT_TRUE(found != pa.end()) << p.first;
            EXPECT_EQ(p.second.unit(), found->second.unit()) << p.first;
            EXPECT_NEAR(p.second.value(), found->second.value(), 1e-6*(1 + std::abs(p.second.value()))) << p.first;
            EXPECT_NEAR(p.second.internalValue(), found->second.internalValue(),
                        1e-6*(1 + std::abs(p.second.internalValue()))) << p.first;
        }
    }
    EXPECT_EQ(Potential::QUARTIC_BONDS, again.findForcesConst(InteractionType::BONDS).potential());
    EXPECT_EQ(Potential::QUARTIC_ANGLES, again.findForcesConst(InteractionType::ANGLES).potential());
}

} // namespace

} // namespace alexandria
