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
#include "act/properties/fragmentwall.h"
#include "act/properties/pathintegralmd.h"
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
     * \param[in]  tm         The molecule
     * \param[in]  settings   The settings, the thermostat is turned off
     * \param[in]  nsteps     The number of steps
     * \param[in]  potential  A potential of the centroids, may be empty
     * \param[out] maxCentroidEnergy The largest energy of that potential, may be nullptr
     * \return the largest difference with the Hamiltonian at the start in kJ/mol
     */
    double hamiltonianChange(const TestMolecule      &tm,
                             RingPolymerSettings      settings,
                             int                      nsteps,
                             const CentroidPotential &potential = {},
                             double                  *maxCentroidEnergy = nullptr)
    {
        settings.thermostat = false;
        RingPolymerSampler sampler(&tm.pd, tm.top(), tm.forceComp.get(), settings);
        sampler.setCentroidPotential(potential);
        MsgHandler         msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        sampler.start(&msghandler, tm.x);
        double h0     = sampler.hamiltonian();
        double change = 0;
        double emax   = sampler.centroidEnergy();
        for(int i = 0; i < nsteps; i++)
        {
            sampler.run(&msghandler, 1);
            change = std::max(change, std::abs(sampler.hamiltonian() - h0));
            emax   = std::max(emax, sampler.centroidEnergy());
        }
        if (maxCentroidEnergy)
        {
            *maxCentroidEnergy = emax;
        }
        return change;
    }

    /*! \brief The first particle of each compound and the number of particles
     * \param[in] tm The molecule
     * \return the indices
     */
    static std::vector<size_t> fragmentStart(const TestMolecule &tm)
    {
        std::vector<size_t> start = { 0 };
        for(const auto &top : tm.mols[0].fragmentHandler()->topologies())
        {
            start.push_back(start.back() + top.atoms().size());
        }
        return start;
    }

    /*! \brief The masses of all particles
     * \param[in] tm The molecule
     * \return the masses
     */
    static std::vector<double> masses(const TestMolecule &tm)
    {
        std::vector<double> mass;
        for(const auto &atom : tm.top()->atoms())
        {
            mass.push_back(atom.mass());
        }
        return mass;
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

TEST_F(RingPolymerTest, HamiltonianConservedWithWall)
{
    // A wall on the distance between the centers of mass of the centroids
    // of two water molecules, at a radius below their distance, such that
    // it pushes them together all the time: its energy times the number
    // of beads is part of the Hamiltonian
    for(const std::string ff : { "SPCE.xml", "SWM4-NDP.xml" })
    {
        TestMolecule tm;
        load(ff, "water_dimer.pdb", &tm);
        const FragmentWall wall(fragmentStart(tm), masses(tm), 0.2, 2000);
        ASSERT_TRUE(wall.active());
        RingPolymerSettings s;
        s.nbeads      = 4;
        s.temperature = 100;
        s.timeStep    = 0.0001;
        double emax   = 0;
        double change = hamiltonianChange(tm, s, 2000,
                                          [&wall](const std::vector<gmx::DVec> &x, std::vector<gmx::DVec> *f)
                                          { return wall.compute(x, f); }, &emax);
        // The kinetic energy is about 9 kT per bead of the polymer at P T
        double scale  = 9.0*s.nbeads*s.nbeads*BOLTZ*s.temperature;
        EXPECT_LT(change, 2e-3*scale) << ff;
        // The wall takes part: its energy is much more than the change
        EXPECT_GT(s.nbeads*emax, 20*change) << ff;
    }
}

TEST_F(RingPolymerTest, EnergyTermsAreAveragesOverTheBeads)
{
    for(const std::string ff : { "SPCE.xml", "TIP5P.xml", "SWM4-NDP.xml" })
    {
        TestMolecule tm;
        load(ff, "water_dimer.pdb", &tm);
        RingPolymerSettings s;
        s.nbeads      = 3;
        s.temperature = 100;
        RingPolymerSampler sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
        // Not asked for, not there
        sampler.start(nullptr, tm.x);
        EXPECT_TRUE(sampler.energyTerms().empty()) << ff;
        // All beads start as the molecule, with its energy terms
        sampler.start(nullptr, tm.x, true);
        std::map<InteractionType, double> eref;
        std::vector<gmx::RVec>            forces;
        auto                              x = tm.x;
        tm.forceComp->compute(nullptr, &tm.pd, tm.top(), &x, &forces, &eref);
        const auto &terms = sampler.energyTerms();
        ASSERT_FALSE(terms.empty()) << ff;
        double toler = GMX_DOUBLE ? 1e-7 : 2e-3;
        for(const auto &e : eref)
        {
            ASSERT_EQ(1U, terms.count(e.first)) << ff << " " << interactionTypeToString(e.first);
            EXPECT_NEAR(e.second, terms.at(e.first), toler) << ff << " " << interactionTypeToString(e.first);
        }
        // After some steps the beads differ, and the terms are those of
        // the last step only
        sampler.run(nullptr, 20, true);
        const auto &later = sampler.energyTerms();
        ASSERT_EQ(1U, later.count(InteractionType::EPOT)) << ff;
        EXPECT_NEAR(sampler.potentialEnergy(), later.at(InteractionType::EPOT), toler) << ff;
        EXPECT_NE(eref[InteractionType::EPOT], later.at(InteractionType::EPOT)) << ff;
        if (tm.pd.polarizable())
        {
            // The shells lower the energy
            ASSERT_EQ(1U, later.count(InteractionType::INDUCTION)) << ff;
            EXPECT_LT(later.at(InteractionType::INDUCTION), 0) << ff;
        }
        sampler.run(nullptr, 1);
        EXPECT_TRUE(sampler.energyTerms().empty()) << ff;
    }
}

TEST_F(RingPolymerTest, EnergyTermsDoNotChangeTheTrajectory)
{
    // With shells the terms take one more energy calculation, which
    // should leave the relaxed shells where they are
    for(const std::string ff : { "SPCE.xml", "SWM4-NDP.xml" })
    {
        TestMolecule tm;
        load(ff, "water_dimer.pdb", &tm);
        double epot[2];
        for(int i = 0; i < 2; i++)
        {
            RingPolymerSettings s;
            s.nbeads      = 2;
            s.temperature = 100;
            s.seed        = 5;
            RingPolymerSampler sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
            sampler.start(nullptr, tm.x, i == 1);
            for(int k = 0; k < 20; k++)
            {
                sampler.run(nullptr, 1, i == 1);
            }
            epot[i] = sampler.potentialEnergy();
        }
        // The shells are relaxed to within the tolerance from another start
        EXPECT_NEAR(epot[0], epot[1], tm.pd.polarizable() ? (GMX_DOUBLE ? 1e-4 : 5e-3) : 0.0) << ff;
    }
}

TEST_F(RingPolymerTest, Centroids)
{
    TestMolecule tm;
    load("SWM4-NDP.xml", "water_dimer.pdb", &tm);
    auto &atoms = tm.top()->atoms();
    RingPolymerSettings s;
    s.nbeads      = 4;
    s.temperature = 100;
    RingPolymerSampler     sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
    std::vector<gmx::DVec> xc;
    sampler.start(nullptr, tm.x);
    sampler.centroids(&xc);
    ASSERT_EQ(atoms.size(), xc.size());
    for(size_t a = 0; a < atoms.size(); a++)
    {
        if (atoms[a].mass() > 0)
        {
            for(int m = 0; m < DIM; m++)
            {
                EXPECT_NEAR(tm.x[a][m], xc[a][m], GMX_DOUBLE ? 1e-12 : 1e-6);
            }
        }
    }
    // The beads of the sampler are relative to the center of mass of the
    // centroids, which is the same shift for all particles
    sampler.run(nullptr, 50);
    sampler.centroids(&xc);
    const auto rp = sampler.polymer();
    gmx::DVec  shift0 = { 0, 0, 0 };
    for(size_t a = 0; a < atoms.size(); a++)
    {
        for(int m = 0; m < DIM; m++)
        {
            double mean = 0;
            for(size_t j = 0; j < rp.x.size(); j++)
            {
                mean += rp.x[j][a][m]/rp.x.size();
            }
            if (0 == a)
            {
                shift0[m] = xc[a][m] - mean;
            }
            EXPECT_NEAR(shift0[m], xc[a][m] - mean, GMX_DOUBLE ? 1e-10 : 2e-6) << "particle " << a;
        }
    }
    // The centroids moved
    EXPECT_GT(std::abs(xc[0][XX] - tm.x[0][XX]) + std::abs(xc[0][YY] - tm.x[0][YY]), 0.0);
}

//! Simulations with ring polymers, see pathIntegralMD
class PathIntegralMDTest : public RingPolymerTest
{
protected:
    /*! \brief The settings of a simulation with beads
     * \param[in] nbeads The number of beads
     * \param[in] T      The temperature
     * \param[in] nsteps The number of steps
     * \param[in] deltat The time step
     * \param[in] tau    The time constant of the thermostat
     * \param[in] seed   The seed
     * \return the settings, without output
     */
    static SimulationConfigHandler config(int nbeads, double T, int nsteps, double deltat, double tau, int seed)
    {
        SimulationConfigHandler sch;
        sch.setSimulation(nsteps, deltat, T, seed, 0, 0);
        sch.setBeads(nbeads, tau, -1);
        sch.validate();
        return sch;
    }
};

TEST_F(PathIntegralMDTest, QuantumEnergyWaterDimer)
{
    // The energies of a cold water dimer with point charges only, with
    // virtual sites and with shells: twelve vibrations, of which six
    // between the molecules, as harmonic ring polymers
    for(const std::string ff : { "SPCE.xml", "TIP5P.xml", "SWM4-NDP.xml" })
    {
        TestMolecule tm;
        load(ff, "water_dimer.pdb", &tm);
        const auto omega = frequencies(tm);
        ASSERT_EQ(12U, omega.size()) << ff;
        std::map<InteractionType, double> emin;
        std::vector<gmx::RVec>            forces;
        auto                              x = tm.x;
        tm.forceComp->compute(nullptr, &tm.pd, tm.top(), &x, &forces, &emin);
        const double T      = 50;
        const int    nbeads = 16;
        const auto   sch    = config(nbeads, T, 24000, 0.00025, 0.02, 2468);
        MsgHandler   msghandler;
        msghandler.setPrintLevel(ACTStatus::Fatal);
        const auto   res = pathIntegralMD(&msghandler, &tm.pd, tm.top(), fragmentStart(tm), tm.x,
                                          tm.forceComp.get(), sch, sch.seed(), PathIntegralFiles());
        EXPECT_TRUE(res.finished);
        EXPECT_EQ(24000, res.nsteps);
        EXPECT_EQ(2400, res.nequil);
        EXPECT_EQ(21600U, res.epot.n);
        EXPECT_EQ(2U, res.nfragments);
        EXPECT_FALSE(res.movedApart) << ff;
        double kT   = BOLTZ*T;
        double eref = harmonicPotentialEnergy(omega, T, nbeads);
        EXPECT_NEAR(eref, res.epot.average - emin[InteractionType::EPOT], 0.06*eref) << ff;
        // The kinetic energy of harmonic vibrations is their potential
        // energy, the rotation and the translation have kT/2 each. The
        // motions of the molecules with respect to each other are not
        // harmonic, which lowers the kinetic energy by 3 to 9 %. That it
        // is correct is tested in KineticEnergyEstimators.
        EXPECT_NEAR(eref + 3*kT, res.ekin.average, 0.13*(eref + 3*kT)) << ff;
        EXPECT_LT(res.ekin.average, eref + 3*kT) << ff;
        EXPECT_NEAR(res.epot.average + res.ekin.average, res.etot.average, 1e-9) << ff;
        EXPECT_GT(res.epot.error, 0) << ff;
        EXPECT_LT(res.epot.error, 0.03*eref) << ff;
        // The molecules stay where they were
        EXPECT_GT(res.distance.average, 0.25) << ff;
        EXPECT_LT(res.distance.average, 0.35) << ff;
        EXPECT_GT(res.gyration.average, 0.005) << ff;
        EXPECT_TRUE(res.terms.empty()) << ff;
    }
}

TEST_F(PathIntegralMDTest, ObserverGetsEnergiesAndCentroids)
{
    // A program that uses the simulation as a library gets the energies
    // and the centroids without files, and can end the simulation
    TestMolecule tm;
    load("SWM4-NDP.xml", "water_dimer.pdb", &tm);
    SimulationConfigHandler sch;
    sch.setSimulation(200, 0.00025, 100, 4321, 10, 20);
    sch.setBeads(4, 0.05, 0);
    sch.validate();
    std::vector<int>    energySteps, coordinateSteps;
    std::vector<double> epot, temperature;
    size_t              ncoord = 0;
    SimulationObserver  observer;
    observer.energies = [&](int step, double time, double temp, double ekin,
                            const std::map<InteractionType, double> &energies)
    {
        energySteps.push_back(step);
        // The time step is a real number
        EXPECT_NEAR(step*0.00025, time, GMX_DOUBLE ? 1e-12 : 1e-8);
        EXPECT_GT(ekin, 0.0);
        temperature.push_back(temp);
        epot.push_back(energies.at(InteractionType::EPOT));
        EXPECT_EQ(1U, energies.count(InteractionType::INDUCTION));
        return true;
    };
    observer.coordinates = [&](int step, double, const std::vector<gmx::RVec> &x)
    {
        coordinateSteps.push_back(step);
        ncoord = x.size();
        // Enough after the frame of step 40
        return step < 40;
    };
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    const auto res = pathIntegralMD(&msghandler, &tm.pd, tm.top(), fragmentStart(tm), tm.x,
                                    tm.forceComp.get(), sch, sch.seed(), PathIntegralFiles(), &observer);
    EXPECT_FALSE(res.finished);
    EXPECT_EQ(40, res.nsteps);
    // What was due for the last step was given
    EXPECT_EQ(std::vector<int>({ 0, 10, 20, 30, 40 }), energySteps);
    EXPECT_EQ(std::vector<int>({ 0, 20, 40 }), coordinateSteps);
    // All particles, with the shells
    EXPECT_EQ(tm.top()->atoms().size(), ncoord);
    ASSERT_EQ(5U, epot.size());
    EXPECT_NE(epot[0], epot[4]);
    // The averages are over the steps that were done, without standard
    // errors: of the eight blocks of 25 steps only one is complete
    EXPECT_EQ(41U, res.epot.n);
    EXPECT_EQ(0.0, res.epot.error);
    ASSERT_EQ(1U, res.terms.count(InteractionType::EPOT));
    EXPECT_EQ(5U, res.terms.at(InteractionType::EPOT).n);
    double mean = 0;
    for(auto e : epot)
    {
        mean += e/epot.size();
    }
    EXPECT_NEAR(mean, res.terms.at(InteractionType::EPOT).average, 1e-9);
}

TEST_F(RingPolymerTest, KineticEnergyEstimators)
{
    // The centroid virial estimator of the kinetic energy has the same
    // average as the primitive estimator, 3NP kT/2 minus the energy of
    // the springs over P, for any potential. A water dimer is not harmonic.
    TestMolecule tm;
    load("SPCE.xml", "water_dimer.pdb", &tm);
    RingPolymerSettings s;
    s.nbeads      = 16;
    s.temperature = 50;
    s.timeStep    = 0.00025;
    s.centroidTau = 0.02;
    s.seed        = 99;
    RingPolymerSampler sampler(&tm.pd, tm.top(), tm.forceComp.get(), s);
    sampler.start(nullptr, tm.x);
    sampler.run(nullptr, 4000);
    const int    nsteps = 40000;
    const double natom  = 6;
    // The center of mass of the centroids does not move
    const double ndof   = 3*natom*s.nbeads - 3;
    double       virial = 0, primitive = 0;
    for(int i = 0; i < nsteps; i++)
    {
        sampler.run(nullptr, 1);
        double spring = (sampler.hamiltonian() - 0.5*ndof*BOLTZ*sampler.beadTemperature() -
                         s.nbeads*sampler.potentialEnergy());
        virial    += sampler.kineticEnergy()/nsteps;
        primitive += (1.5*natom*s.nbeads*BOLTZ*s.temperature - spring/s.nbeads)/nsteps;
    }
    // Much more than the classical 9 kT
    EXPECT_GT(virial, 5*9*BOLTZ*s.temperature);
    // The primitive estimator has the larger fluctuations, about 1 %
    EXPECT_NEAR(primitive, virial, 0.04*virial);
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

