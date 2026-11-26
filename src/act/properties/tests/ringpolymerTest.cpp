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
 * Tests of the path integral molecular dynamics of a ring polymer.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/basics/act_particle.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forces/forcecomputer.h"
#include "act/import/compound_reader.h"
#include "act/properties/ringpolymer.h"
#include "act/topology/actmol.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! Planck's constant divided by 2 pi in kJ/mol ps
const double hbar = PLANCK/(2*M_PI);

//! A molecule with its force field, ready for a ring polymer simulation
struct TestMolecule
{
    //! The force field
    ForceField                     pd;
    //! The force computer
    std::unique_ptr<ForceComputer> forceComp;
    //! The molecule
    std::vector<ACTMol>            mols;
    //! The minimized coordinates of all particles
    std::vector<gmx::RVec>         x;
    //! \return the topology of the molecule
    const Topology *top() const { return mols[0].topology(); }
};

class RingPolymerTest : public gmx::test::CommandLineTestBase
{
protected:
    /*! \brief Read a molecule and minimize its energy
     * \param[in]  ffFile  Force field in src/act/testdata/forcefields
     * \param[in]  molFile Molecule in src/act/testdata/molecules
     * \param[out] tm      The molecule
     */
    void load(const std::string &ffFile,
              const std::string &molFile,
              TestMolecule      *tm)
    {
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        readForceField(fileManager().getInputFilePath("forcefields/" + ffFile), &tm->pd, &msghandler);
        tm->forceComp = std::make_unique<ForceComputer>(1e-8, 100, tm->pd.epsilonr());
        CompoundReader compR;
        compR.setFileName(fileManager().getInputFilePath("molecules/" + molFile));
        compR.read(&msghandler, tm->pd, tm->forceComp.get(), &tm->mols);
        ASSERT_EQ(1U, tm->mols.size());
        tm->x = tm->mols[0].xOriginal();
        const MolHandler              molhandler;
        const SimulationConfigHandler sch;
        std::map<InteractionType, double> energies;
        const auto status = molhandler.minimizeCoordinates(&msghandler, &tm->pd, tm->top(),
                                                           tm->forceComp.get(), sch,
                                                           &tm->x, &energies, {});
        ASSERT_EQ(eMinimizeStatus::OK, status);
    }

    /*! \brief Largest change of the Hamiltonian of a ring polymer without thermostat
     * \param[in] tm       The molecule
     * \param[in] settings The settings, the thermostat is turned off
     * \param[in] nsteps   The number of steps
     * \return the largest difference with the Hamiltonian at the start in kJ/mol
     */
    double hamiltonianChange(const TestMolecule  &tm,
                             RingPolymerSettings  settings,
                             int                  nsteps)
    {
        settings.thermostat = false;
        RingPolymerSampler sampler(&tm.pd, tm.top(), tm.forceComp.get(), settings);
        MsgHandler         msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        sampler.start(&msghandler, tm.x);
        double h0     = sampler.hamiltonian();
        double change = 0;
        for(int i = 0; i < nsteps; i++)
        {
            sampler.run(&msghandler, 1);
            change = std::max(change, std::abs(sampler.hamiltonian() - h0));
        }
        return change;
    }

    /*! \brief The vibrational frequencies of a molecule
     * \param[in] tm The molecule
     * \return the angular frequencies in 1/ps
     */
    std::vector<double> frequencies(const TestMolecule &tm)
    {
        const MolHandler    molhandler;
        std::vector<double> freq, inten, omega;
        auto                x = tm.x;
        molhandler.nma(&tm.pd, &tm.mols[0], tm.forceComp.get(), &x, &freq, &inten);
        for(auto f : freq)
        {
            // The vibrations only, in 1/ps
            omega.push_back(2*M_PI*f);
        }
        return omega;
    }
};

/*! \brief The potential energy of harmonic oscillators as ring polymers
 * The exact average for P beads, which becomes the quantum mechanical
 * one, hbar omega/4 coth(beta hbar omega/2), for many beads and the
 * classical one, kT/2, for one.
 * \param[in] omega  The angular frequencies in 1/ps
 * \param[in] T      The temperature
 * \param[in] nbeads The number of beads
 * \return the potential energy in kJ/mol
 */
double harmonicPotentialEnergy(const std::vector<double> &omega,
                               double                     T,
                               int                        nbeads)
{
    double kT     = BOLTZ*T;
    double omegaP = nbeads*kT/hbar;
    double epot   = 0;
    for(auto w : omega)
    {
        for(int k = 0; k < nbeads; k++)
        {
            double wk = 2*omegaP*std::sin(M_PI*k/nbeads);
            epot += 0.5*kT*w*w/(w*w + wk*wk);
        }
    }
    return epot;
}

TEST_F(RingPolymerTest, SettingsAreChecked)
{
    TestMolecule tm;
    load("SPCE.xml", "water.pdb", &tm);
    {
        RingPolymerSettings s;
        s.nbeads = 0;
        EXPECT_THROW(RingPolymerSampler(&tm.pd, tm.top(), tm.forceComp.get(), s), gmx::InvalidInputError);
    }
    {
        RingPolymerSettings s;
        s.temperature = 0;
        EXPECT_THROW(RingPolymerSampler(&tm.pd, tm.top(), tm.forceComp.get(), s), gmx::InvalidInputError);
    }
    {
        RingPolymerSettings s;
        s.timeStep = -1;
        EXPECT_THROW(RingPolymerSampler(&tm.pd, tm.top(), tm.forceComp.get(), s), gmx::InvalidInputError);
    }
    {
        RingPolymerSettings s;
        s.centroidTau = 0;
        EXPECT_THROW(RingPolymerSampler(&tm.pd, tm.top(), tm.forceComp.get(), s), gmx::InvalidInputError);
    }
    {
        // The springs of many beads need a small time step
        RingPolymerSettings s;
        s.nbeads   = 64;
        s.timeStep = 0.002;
        EXPECT_THROW(RingPolymerSampler(&tm.pd, tm.top(), tm.forceComp.get(), s), gmx::InvalidInputError);
        s.nbeads   = 1;
        EXPECT_NO_THROW(RingPolymerSampler(&tm.pd, tm.top(), tm.forceComp.get(), s));
    }
    {
        const RingPolymerSettings s;
        RingPolymerSampler  sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
        // No coordinates yet
        EXPECT_THROW(sampler.run(nullptr, 1), gmx::InternalError);
        // The wrong number of them
        const std::vector<gmx::RVec> x(2);
        EXPECT_THROW(sampler.start(nullptr, x), gmx::InvalidInputError);
    }
}

TEST_F(RingPolymerTest, HamiltonianConservedWater)
{
    TestMolecule tm;
    load("SPCE.xml", "water.pdb", &tm);
    for(int nbeads : { 1, 2, 5, 8 })
    {
        RingPolymerSettings s;
        s.nbeads      = nbeads;
        s.temperature = 300;
        s.timeStep    = 0.0001;
        // The kinetic energy is about 9/2 kT per bead of the polymer at P T
        double scale  = 4.5*nbeads*nbeads*BOLTZ*s.temperature;
        EXPECT_LT(hamiltonianChange(tm, s, 2000), 2e-3*scale) << nbeads << " beads";
    }
}

TEST_F(RingPolymerTest, HamiltonianConservedVirtualSites)
{
    TestMolecule tm;
    load("TIP5P.xml", "water.pdb", &tm);
    RingPolymerSettings s;
    s.nbeads      = 4;
    s.temperature = 300;
    s.timeStep    = 0.0001;
    double scale  = 4.5*s.nbeads*s.nbeads*BOLTZ*s.temperature;
    EXPECT_LT(hamiltonianChange(tm, s, 2000), 2e-3*scale);
}

TEST_F(RingPolymerTest, HamiltonianConservedShells)
{
    TestMolecule tm;
    load("SWM4-NDP.xml", "water.pdb", &tm);
    RingPolymerSettings s;
    s.nbeads      = 4;
    s.temperature = 300;
    s.timeStep    = 0.0001;
    double scale  = 4.5*s.nbeads*s.nbeads*BOLTZ*s.temperature;
    EXPECT_LT(hamiltonianChange(tm, s, 2000), 2e-3*scale);
    // All shells were relaxed
    EXPECT_EQ(0, tm.forceComp->numShellConvergenceFailed());
}

TEST_F(RingPolymerTest, ClassicalEnergyWater)
{
    TestMolecule tm;
    load("SPCE.xml", "water.pdb", &tm);
    RingPolymerSettings s;
    s.nbeads      = 1;
    s.temperature = 300;
    s.timeStep    = 0.00025;
    s.centroidTau = 0.01;
    s.seed        = 1234;
    RingPolymerSampler       sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
    std::vector<RingPolymer> pool;
    sampler.start(nullptr, tm.x);
    const auto aver = sampler.sample(nullptr, 2000, 20, 20000, &pool);
    ASSERT_EQ(20000U, pool.size());
    // Three vibrations with kT/2 each, relative to the minimum
    std::map<InteractionType, double> emin;
    std::vector<gmx::RVec>            forces;
    auto                              x = tm.x;
    tm.forceComp->compute(nullptr, &tm.pd, tm.top(), &x, &forces, &emin);
    double kT = BOLTZ*s.temperature;
    EXPECT_NEAR(1.5*kT, aver.potentialEnergy() - emin[InteractionType::EPOT], 0.06*1.5*kT);
    // One bead has no spread, and the classical kinetic energy
    EXPECT_EQ(0.0, aver.gyrationSquared());
    EXPECT_NEAR(4.5*kT, aver.kineticEnergy(), 1e-10);
}

TEST_F(RingPolymerTest, QuantumEnergyWater)
{
    TestMolecule tm;
    load("SPCE.xml", "water.pdb", &tm);
    const auto omega = frequencies(tm);
    ASSERT_EQ(3U, omega.size());
    std::map<InteractionType, double> emin;
    std::vector<gmx::RVec>            forces;
    auto                              x = tm.x;
    tm.forceComp->compute(nullptr, &tm.pd, tm.top(), &x, &forces, &emin);
    for(int nbeads : { 4, 16 })
    {
        RingPolymerSettings s;
        s.nbeads      = nbeads;
        s.temperature = 300;
        s.timeStep    = 0.00025;
        s.centroidTau = 0.01;
        s.seed        = 4321;
        RingPolymerSampler       sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
        std::vector<RingPolymer> pool;
        sampler.start(nullptr, tm.x);
        const auto aver = sampler.sample(nullptr, 2000, 10, 8000, &pool);
        double     kT   = BOLTZ*s.temperature;
        // The exact energy of harmonic vibrations with this number of beads.
        // The molecule is not harmonic, and it rotates.
        double     eref = harmonicPotentialEnergy(omega, s.temperature, nbeads);
        EXPECT_GT(eref, 2*1.5*kT);
        EXPECT_NEAR(eref, aver.potentialEnergy() - emin[InteractionType::EPOT], 0.05*eref) << nbeads << " beads";
        // The kinetic energy of the vibrations is their potential energy,
        // the translation and rotation have nearly kT/2 each
        EXPECT_NEAR(eref + 3*kT, aver.kineticEnergy(), 0.05*(eref + 3*kT)) << nbeads << " beads";
        // The beads of the hydrogen atoms are spread over a few pm
        EXPECT_GT(aver.gyrationSquared(), 0.0);
        EXPECT_LT(aver.gyrationSquared(), 0.02*0.02);
    }
}

TEST_F(RingPolymerTest, PolymersAreCentered)
{
    TestMolecule tm;
    load("SWM4-NDP.xml", "water.pdb", &tm);
    auto &atoms = tm.top()->atoms();
    RingPolymerSettings s;
    s.nbeads      = 4;
    s.temperature = 300;
    RingPolymerSampler       sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
    std::vector<RingPolymer> pool;
    sampler.start(nullptr, tm.x);
    (void) sampler.sample(nullptr, 200, 50, 5, &pool);
    ASSERT_EQ(5U, pool.size());
    double shellMax = shellMaxDistance(&tm.pd);
    for(const auto &rp : pool)
    {
        ASSERT_EQ(4U, rp.x.size());
        ASSERT_EQ(4U, rp.epot.size());
        gmx::DVec com   = { 0, 0, 0 };
        double    mtot  = 0;
        for(size_t j = 0; j < rp.x.size(); j++)
        {
            ASSERT_EQ(atoms.size(), rp.x[j].size());
            EXPECT_TRUE(std::isfinite(rp.epot[j]));
            for(size_t a = 0; a < atoms.size(); a++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    com[m] += atoms[a].mass()*rp.x[j][a][m];
                }
                mtot += atoms[a].mass();
                if (ActParticle::Shell == atoms[a].pType())
                {
                    // The shell is close to its parent. In this model it
                    // only moves in the field of another molecule.
                    gmx::RVec dx;
                    rvec_sub(rp.x[j][a], rp.x[j][atoms[a].parents()[0]], dx);
                    if (shellMax > 0)
                    {
                        EXPECT_LE(norm(dx), shellMax*1.0001);
                    }
                }
            }
        }
        for(int m = 0; m < DIM; m++)
        {
            EXPECT_NEAR(0.0, com[m]/mtot, GMX_DOUBLE ? 1e-12 : 1e-6);
        }
    }
    // The configurations differ
    EXPECT_NE(pool[0].epot[0], pool[1].epot[0]);
}

TEST_F(RingPolymerTest, NoStrideGivesTheSamePolymers)
{
    TestMolecule tm;
    load("SPCE.xml", "water.pdb", &tm);
    RingPolymerSettings s;
    s.nbeads = 3;
    RingPolymerSampler       sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
    std::vector<RingPolymer> pool;
    sampler.start(nullptr, tm.x);
    (void) sampler.sample(nullptr, 0, 0, 3, &pool);
    ASSERT_EQ(3U, pool.size());
    for(const auto &rp : pool)
    {
        for(size_t j = 0; j < rp.x.size(); j++)
        {
            EXPECT_EQ(pool[0].epot[0], rp.epot[j]);
            for(size_t a = 0; a < rp.x[j].size(); a++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    EXPECT_EQ(pool[0].x[0][a][m], rp.x[j][a][m]);
                }
            }
        }
    }
}

TEST_F(RingPolymerTest, SameSeedSamePolymers)
{
    TestMolecule tm;
    load("SPCE.xml", "water.pdb", &tm);
    std::vector<RingPolymer> pools[3];
    for(int i = 0; i < 3; i++)
    {
        RingPolymerSettings s;
        s.nbeads = 4;
        s.seed   = (i < 2) ? 7 : 8;
        RingPolymerSampler sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
        sampler.start(nullptr, tm.x);
        (void) sampler.sample(nullptr, 100, 10, 2, &pools[i]);
    }
    EXPECT_EQ(pools[0][1].epot[2], pools[1][1].epot[2]);
    EXPECT_NE(pools[0][1].epot[2], pools[2][1].epot[2]);
}

} // namespace

} // namespace alexandria
