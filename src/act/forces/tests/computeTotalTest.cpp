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
 * Tests comparing ForceComputer::computeTotal and computeInteraction,
 * which use the SIMD routines, with compute, computeOnce and
 * ACTMol::calculateInteractionEnergy.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "../forcecomputer.h"

#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/import/import.h"
#include "act/topology/actmol.h"
#include "act/topology/fragmenthandler.h"
#include "gromacs/math/vectypes.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/*! \brief Read a molecule and make its topology
 * \param[in]  msghandler For errors
 * \param[in]  pd         The force field
 * \param[in]  molName    File name of the molecule without .sdf
 * \param[out] mol        The molecule
 * \return Whether it worked
 */
bool readMolecule(MsgHandler        *msghandler,
                  ForceField        *pd,
                  const std::string &molName,
                  ACTMol            *mol)
{
    std::vector<MolProp> molprops;
    bool   userqtot   = false;
    double qtot_babel = 0;
    auto   dataName   = gmx::test::TestFileManager::getInputFilePath(molName + ".sdf");
    importFile(msghandler, pd, dataName.c_str(), &molprops,
               "minimum", JobType::OPT, userqtot, &qtot_babel, true);
    if (!msghandler->ok() || molprops.size() != 1)
    {
        return false;
    }
    mol->Merge(&molprops[0]);
    mol->GenerateTopology(msghandler, pd, missingParameters::Ignore);
    return msghandler->ok();
}

/*! \brief Compare two energy maps
 * \param[in] ref   The reference energies
 * \param[in] test  The energies to test
 * \param[in] toler Relative tolerance, also used as absolute tolerance
 * \param[in] label For messages
 */
void compareEnergies(const std::map<InteractionType, double> &ref,
                     const std::map<InteractionType, double> &test,
                     double                                   toler,
                     const std::string                       &label)
{
    for(const auto &e : ref)
    {
        auto t = test.find(e.first);
        if (test.end() == t)
        {
            // Terms that are zero may be left out, as VDWCORRECTION
            EXPECT_EQ(0, e.second) << label << " missing " << interactionTypeToString(e.first);
            continue;
        }
        EXPECT_NEAR(e.second, t->second, toler*std::max(1.0, std::abs(e.second)))
            << label << " " << interactionTypeToString(e.first);
    }
    for(const auto &t : test)
    {
        EXPECT_TRUE(ref.end() != ref.find(t.first))
            << label << " extra " << interactionTypeToString(t.first);
    }
}

/*! \brief Compare two sets of forces
 * \param[in] ref   The reference forces
 * \param[in] test  The forces to test
 * \param[in] toler Relative tolerance, also used as absolute tolerance
 * \param[in] label For messages
 */
void compareForces(const std::vector<gmx::RVec> &ref,
                   const std::vector<gmx::RVec> &test,
                   double                        toler,
                   const std::string            &label)
{
    ASSERT_EQ(ref.size(), test.size()) << label;
    for(size_t i = 0; i < ref.size(); i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            EXPECT_NEAR(ref[i][m], test[i][m], toler*std::max(1.0, std::abs(ref[i][m])))
                << label << " particle " << i << " dim " << m;
        }
    }
}

/*! \brief Coordinates displaced a little, differently for each conformation
 * \param[in] x The coordinates
 * \param[in] c The conformation
 * \return The displaced coordinates
 */
std::vector<gmx::RVec> displaced(const std::vector<gmx::RVec> &x, size_t c)
{
    auto y = x;
    for(size_t i = 0; i < y.size(); i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            y[i][m] += 0.002*std::sin(1.0 + c + 3*i + m);
        }
    }
    return y;
}

class ComputeTotalTest : public ::testing::Test
{
protected:
    MsgHandler    msghandler_;
    //! Tight convergence, such that both methods find the same minimum
    ForceComputer forceComp_ = ForceComputer(1e-12, 200, 0.04, 1);

    ComputeTotalTest()
    {
        msghandler_.setPrintLevel(ACTStatus::Warning);
    }

    /*! \brief Compare computeTotal with compute for one or more conformations
     * \param[in] ffName  The force field
     * \param[in] molName The molecule
     * \param[in] nconf   The number of conformations
     * \param[in] options The shell options
     */
    void testTotal(const std::string  &ffName,
                   const std::string  &molName,
                   size_t              nconf,
                   const ShellOptions &options)
    {
        auto pd = getForceField(ffName);
        ASSERT_NE(nullptr, pd);
        ACTMol mol;
        ASSERT_TRUE(readMolecule(&msghandler_, pd, molName, &mol));
        auto top = mol.topology();
        std::vector<std::vector<gmx::RVec>> x;
        for(size_t c = 0; c < nconf; c++)
        {
            x.push_back(displaced(mol.xOriginal(), c));
        }
        ConformationSet set;
        set.top = top;
        set.setCoordinates(x);
        std::vector<std::map<InteractionType, double>> energies;
        forceComp_.computeTotal(&msghandler_, pd, &set, options, &energies);
        ASSERT_EQ(nconf, energies.size());
        for(size_t c = 0; c < nconf; c++)
        {
            auto label = gmx::formatString("%s %s conformation %zu", ffName.c_str(), molName.c_str(), c);
            auto xref  = x[c];
            std::vector<gmx::RVec>            fref(xref.size());
            std::map<InteractionType, double> eref;
            forceComp_.compute(&msghandler_, pd, top, &xref, &fref, &eref,
                               options.field, options.resetShells, options.relaxShells);
            compareEnergies(eref, energies[c], 1e-8, label);
            compareForces(fref, set.confs[c].fcalc, 1e-6, label);
        }
    }

    /*! \brief Compare computeInteraction with calculateInteractionEnergy
     * \param[in] ffName     The force field
     * \param[in] molName    The dimer
     * \param[in] nconf      The number of conformations
     * \param[in] separate   Whether to treat the induction correction separately
     * \param[in] correction Parameters for the induction correction potential, if not empty
     */
    void testInteraction(const std::string         &ffName,
                         const std::string         &molName,
                         size_t                     nconf,
                         bool                       separate   = true,
                         const std::vector<double> &correction = {})
    {
        auto pd = getForceField(ffName);
        ASSERT_NE(nullptr, pd);
        ACTMol mol;
        ASSERT_TRUE(readMolecule(&msghandler_, pd, molName, &mol));
        ASSERT_TRUE(mol.isDimer());
        if (!correction.empty())
        {
            // In the dimer and in the monomers
            std::vector<std::map<InteractionType, TopologyEntryVector> *> allEntries = { mol.topologyPtr()->entriesPtr() };
            for(auto &top : mol.fragmentHandler()->topologiesPtr())
            {
                allEntries.push_back(top.entriesPtr());
            }
            for(auto entries : allEntries)
            {
                // The monomers may have none, if all pairs are excluded
                auto ee = entries->find(InteractionType::INDUCTIONCORRECTION);
                if (entries->end() == ee)
                {
                    continue;
                }
                for(auto &e : ee->second)
                {
                    e->setParams(correction);
                }
            }
        }
        std::vector<std::vector<gmx::RVec>> x;
        for(size_t c = 0; c < nconf; c++)
        {
            x.push_back(displaced(mol.xOriginal(), c));
        }
        auto xinter = x;
        std::vector<std::map<InteractionType, double>> einter;
        std::vector<double>                            epot;
        std::vector<std::vector<gmx::RVec>>            finter;
        forceComp_.computeInteraction(&msghandler_, pd, &mol, &xinter, separate,
                                      &einter, &epot, &finter);
        ASSERT_EQ(nconf, einter.size());
        for(size_t c = 0; c < nconf; c++)
        {
            auto label = gmx::formatString("%s %s conformation %zu", ffName.c_str(), molName.c_str(), c);
            auto xref  = x[c];
            std::map<InteractionType, double> eref;
            std::vector<gmx::RVec>            fref;
            double epotref = mol.calculateInteractionEnergy(&msghandler_, pd, &forceComp_,
                                                            &eref, &fref, &xref, separate);
            // The split of the electrostatics between COULOMB and the
            // INDUCTIONCORRECTION differs, since calculateInteractionEnergy
            // takes the electrostatics before relaxation with the shells
            // in the relaxed monomers and computeRelaxedEnergies with the
            // shells on their cores. The sums and all other terms agree.
            auto icoul = InteractionType::ELECTROSTATICS;
            auto icorr = InteractionType::INDUCTIONCORRECTION;
            auto iexci = InteractionType::EXCHIND;
            auto icoulCorr = [&](const std::map<InteractionType, double> &e)
            {
                return e.at(icoul) + (e.count(icorr) ? e.at(icorr) : 0.0);
            };
            EXPECT_NEAR(icoulCorr(eref), icoulCorr(einter[c]),
                        1e-7*std::max(1.0, std::abs(icoulCorr(eref)))) << label;
            auto eref2 = eref, einter2 = einter[c];
            for(auto t : { icoul, icorr, iexci })
            {
                eref2.erase(t);
                einter2.erase(t);
            }
            // calculateInteractionEnergy moves the polarization to the electrostatics
            EXPECT_EQ(0, einter2[InteractionType::POLARIZATION]) << label;
            einter2.erase(InteractionType::POLARIZATION);
            compareEnergies(eref2, einter2, 1e-7, label);
            EXPECT_NEAR(epotref, epot[c], 1e-7*std::max(1.0, std::abs(epotref))) << label;
            compareForces(fref, finter[c], 1e-5, label);
        }
    }
};

TEST_F(ComputeTotalTest, PolarizableVsites)
{
    testTotal("ACS-pg-vs2", "hydrogen-fluoride", 1, ShellOptions());
}

TEST_F(ComputeTotalTest, PolarizableVsitesConformations)
{
    testTotal("ACS-pg-vs2", "hydrogen-fluoride", 3, ShellOptions());
}

TEST_F(ComputeTotalTest, PolarizableWater)
{
    testTotal("ACS-pg-vs3", "water", 3, ShellOptions());
}

TEST_F(ComputeTotalTest, HyperpolarizableDimer)
{
    testTotal("water_hyperpol", "water_dimer", 3, ShellOptions());
}

TEST_F(ComputeTotalTest, PolarizableDimer)
{
    testTotal("ACS-pg-vs2", "hbr-hf", 3, ShellOptions());
}

TEST_F(ComputeTotalTest, InductionCorrectionPotential)
{
    // MACDANIEL_SCHMIDT for the INDUCTIONCORRECTION. Its parameters are
    // zero in the force field, therefore they are set in the topology.
    auto pd = getForceField("ACS-pol");
    ASSERT_NE(nullptr, pd);
    ACTMol mol;
    ASSERT_TRUE(readMolecule(&msghandler_, pd, "water_dimer", &mol));
    auto itype   = InteractionType::INDUCTIONCORRECTION;
    auto entries = mol.topologyPtr()->entriesPtr();
    ASSERT_TRUE(entries->end() != entries->find(itype));
    for(auto &e : entries->find(itype)->second)
    {
        e->setParams({ 2000.0, 1000.0, 20.0 });
    }
    auto top = mol.topology();
    std::vector<std::vector<gmx::RVec>> x = { displaced(mol.xOriginal(), 0),
                                              displaced(mol.xOriginal(), 1),
                                              displaced(mol.xOriginal(), 2) };
    ConformationSet set;
    set.top = top;
    set.setCoordinates(x);
    std::vector<std::map<InteractionType, double>> energies;
    forceComp_.computeTotal(&msghandler_, pd, &set, ShellOptions(), &energies);
    for(size_t c = 0; c < x.size(); c++)
    {
        auto label = gmx::formatString("conformation %zu", c);
        auto xref  = x[c];
        std::vector<gmx::RVec>            fref(xref.size());
        std::map<InteractionType, double> eref;
        forceComp_.compute(&msghandler_, pd, top, &xref, &fref, &eref);
        // The correction is booked as such, not as exchange
        ASSERT_TRUE(eref.end() != eref.find(itype));
        EXPECT_NE(0, eref[itype]) << label;
        compareEnergies(eref, energies[c], 1e-8, label);
        compareForces(fref, set.confs[c].fcalc, 1e-6, label);
    }
}

TEST_F(ComputeTotalTest, ElectricField)
{
    ShellOptions options;
    options.field = { 0.5, -1.0, 2.0 };
    testTotal("ACS-pg-vs2", "hydrogen-fluoride", 3, options);
}

TEST_F(ComputeTotalTest, NoReset)
{
    ShellOptions options;
    options.resetShells = false;
    testTotal("ACS-pg-vs3", "water", 2, options);
}

TEST_F(ComputeTotalTest, NoRelaxation)
{
    // Compare with computeOnce, which does not relax the shells
    auto pd = getForceField("ACS-pg-vs2");
    ASSERT_NE(nullptr, pd);
    ACTMol mol;
    ASSERT_TRUE(readMolecule(&msghandler_, pd, "hbr-hf", &mol));
    auto top = mol.topology();
    std::vector<std::vector<gmx::RVec>> x = { displaced(mol.xOriginal(), 0),
                                              displaced(mol.xOriginal(), 1) };
    // The shells away from their cores, as during a minimization
    for(size_t c = 0; c < x.size(); c++)
    {
        forceComp_.constructVsiteCoordinates(top, &x[c]);
        for(size_t i = 0; i < x[c].size(); i++)
        {
            if (top->atoms()[i].pType() == ActParticle::Shell)
            {
                x[c][i][XX] += 0.001*(c+1);
            }
        }
    }
    ShellOptions options;
    options.relax       = false;
    options.resetShells = false;
    ConformationSet set;
    set.top = top;
    set.setCoordinates(x);
    std::vector<std::map<InteractionType, double>> energies;
    forceComp_.computeTotal(&msghandler_, pd, &set, options, &energies);
    for(size_t c = 0; c < x.size(); c++)
    {
        auto label = gmx::formatString("conformation %zu", c);
        auto xref  = x[c];
        std::vector<gmx::RVec>            fref(xref.size());
        std::map<InteractionType, double> eref;
        forceComp_.computeOnce(&msghandler_, pd, top, &xref, &fref, &eref, { 0, 0, 0 });
        forceComp_.spreadVsiteForces(top, &xref, &fref);
        // computeOnce gives the potential terms only
        for(const auto &e : eref)
        {
            auto t = energies[c].find(e.first);
            ASSERT_TRUE(energies[c].end() != t) << label << " " << interactionTypeToString(e.first);
            EXPECT_NEAR(e.second, t->second, 1e-10*std::max(1.0, std::abs(e.second)))
                << label << " " << interactionTypeToString(e.first);
        }
        compareForces(fref, set.confs[c].fcalc, 1e-10, label);
        // The shells did not move
        for(size_t i = 0; i < x[c].size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                EXPECT_EQ(x[c][i][m], set.confs[c].x[i][m]) << label;
            }
        }
    }
}

TEST_F(ComputeTotalTest, RelaxSomeShells)
{
    // Relax the shell of the first monomer only
    auto pd = getForceField("ACS-pg-vs2");
    ASSERT_NE(nullptr, pd);
    ACTMol mol;
    ASSERT_TRUE(readMolecule(&msghandler_, pd, "hbr-hf", &mol));
    auto  top    = mol.topology();
    auto &atoms  = top->atoms();
    auto  natom1 = mol.fragmentHandler()->atomStart()[1];
    ShellOptions options;
    for(size_t i = 0; i < natom1; i++)
    {
        if (atoms[i].pType() == ActParticle::Shell)
        {
            options.relaxShells.insert(i);
        }
    }
    ASSERT_FALSE(options.relaxShells.empty());
    std::vector<std::vector<gmx::RVec>> x = { displaced(mol.xOriginal(), 0) };
    ConformationSet set;
    set.top = top;
    set.setCoordinates(x);
    std::vector<std::map<InteractionType, double>> energies;
    forceComp_.computeTotal(&msghandler_, pd, &set, options, &energies);
    // The other shells stay on their cores
    for(size_t i = natom1; i < atoms.size(); i++)
    {
        if (atoms[i].pType() == ActParticle::Shell)
        {
            auto core = atoms[i].parents()[0];
            for(int m = 0; m < DIM; m++)
            {
                EXPECT_EQ(set.confs[0].x[core][m], set.confs[0].x[i][m]);
            }
        }
    }
    // The relaxed shells give the same energies as compute
    auto xref = x[0];
    std::vector<gmx::RVec>            fref(xref.size());
    std::map<InteractionType, double> eref;
    forceComp_.compute(&msghandler_, pd, top, &xref, &fref, &eref,
                       options.field, options.resetShells, options.relaxShells);
    compareEnergies(eref, energies[0], 1e-8, "relax some");
}

TEST_F(ComputeTotalTest, Statistics)
{
    auto pd = getForceField("ACS-pg-vs3");
    ASSERT_NE(nullptr, pd);
    ACTMol mol;
    ASSERT_TRUE(readMolecule(&msghandler_, pd, "water", &mol));
    ConformationSet set;
    set.top = mol.topology();
    set.setCoordinates({ mol.xOriginal(), displaced(mol.xOriginal(), 1) });
    std::vector<std::map<InteractionType, double>> energies;
    forceComp_.resetStatistics();
    forceComp_.computeTotal(&msghandler_, pd, &set, ShellOptions(), &energies);
    EXPECT_EQ(2u, forceComp_.numEvaluations());
    EXPECT_GE(forceComp_.numShellIterations(), 2u);
    EXPECT_EQ(0u, forceComp_.numShellConvergenceFailed());
}

TEST_F(ComputeTotalTest, InteractionPolarizable)
{
    testInteraction("ACS-pg-vs2", "hbr-hf", 3);
}

TEST_F(ComputeTotalTest, InteractionHyperpolarizable)
{
    testInteraction("water_hyperpol", "water_dimer", 3);
}

TEST_F(ComputeTotalTest, InteractionNotSeparated)
{
    testInteraction("water_hyperpol", "water_dimer", 3, false);
}

TEST_F(ComputeTotalTest, InteractionCorrectionPotential)
{
    // MACDANIEL_SCHMIDT for the INDUCTIONCORRECTION, with parameters
    // since they are zero in the force field
    testInteraction("ACS-pol", "water_dimer", 3, true, { 2000.0, 1000.0, 20.0 });
}

TEST_F(ComputeTotalTest, InteractionCorrectionPotentialNotSeparated)
{
    testInteraction("ACS-pol", "water_dimer", 3, false, { 2000.0, 1000.0, 20.0 });
}

} // namespace

} // namespace alexandria
