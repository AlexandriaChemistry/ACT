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
 * which use the SIMD routines, with compute and computeOnce.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "../forcecomputer.h"

#include <algorithm>
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
#include "testutils/testasserts.h"
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
    const auto dataName   = gmx::test::TestFileManager::getInputFilePath("molecules/" + molName + ".sdf");
    importFile(msghandler, pd, dataName, &molprops,
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
 * In single precision the energies are compared relative to the largest
 * one, since some terms, like the induction, are differences of others.
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
    double emax = 1;
    for(const auto &e : ref)
    {
        emax = std::max(emax, std::abs(e.second));
    }
    for(const auto &e : ref)
    {
        const auto t = test.find(e.first);
        if (test.end() == t)
        {
            // Terms that are zero may be left out, as VDWCORRECTION
            EXPECT_EQ(0, e.second) << label << " missing " << interactionTypeToString(e.first);
            continue;
        }
        EXPECT_NEAR(e.second, t->second, toler*(GMX_DOUBLE ? std::max(1.0, std::abs(e.second)) : emax))
            << label << " " << interactionTypeToString(e.first);
    }
    for(const auto &t : test)
    {
        EXPECT_TRUE(ref.end() != ref.find(t.first))
            << label << " extra " << interactionTypeToString(t.first);
    }
}

/*! \brief Compare two sets of forces
 * In single precision the forces are compared relative to the largest
 * one, since the terms of the forces cancel.
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
    double fmax = 1;
    for(const auto &f : ref)
    {
        for(int m = 0; m < DIM; m++)
        {
            fmax = std::max<double>(fmax, std::abs(f[m]));
        }
    }
    for(size_t i = 0; i < ref.size(); i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            EXPECT_NEAR(ref[i][m], test[i][m], toler*(GMX_DOUBLE ? std::max<double>(1.0, std::abs(ref[i][m])) : fmax))
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
    /*! Tight convergence, such that both methods find the same minimum, in
     * single precision as far as the precision of the forces allows
     */
    ForceComputer forceComp_ = ForceComputer(gmx::test::precisionTolerance(1e-12, 1e-6), 200, 1);

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
            const auto label = gmx::formatString("%s %s conformation %zu", ffName.c_str(), molName.c_str(), c);
            auto xref  = x[c];
            std::vector<gmx::RVec>            fref(xref.size());
            std::map<InteractionType, double> eref;
            forceComp_.compute(&msghandler_, pd, top, &xref, &fref, &eref,
                               options.field, options.resetShells, options.relaxShells);
            compareEnergies(eref, energies[c], gmx::test::precisionTolerance(1e-8, 1e-6), label);
            compareForces(fref, set.data.f(c).toRVec(), gmx::test::precisionTolerance(1e-6, 1e-5), label);
        }
    }

    /*! \brief Compare computeInteraction for many dimers with subsets of them
     * Each dimer is computed independently, also when the shells of some
     * converge earlier than those of others, such that computing the
     * dimers together, in two subsets or in another order gives the same
     * energies, forces and shell positions. The second monomer is moved by
     * a different distance in each dimer, such that the shells converge
     * after different numbers of iterations.
     * \param[in] ffName  The force field
     * \param[in] molName The dimer
     * \param[in] nconf   The number of dimers
     */
    void testInteractionSubsets(const std::string &ffName,
                                const std::string &molName,
                                size_t             nconf)
    {
        auto pd = getForceField(ffName);
        ASSERT_NE(nullptr, pd);
        ACTMol mol;
        ASSERT_TRUE(readMolecule(&msghandler_, pd, molName, &mol));
        ASSERT_TRUE(mol.isDimer());
        size_t natom1 = mol.fragmentHandler()->atomStart()[1];
        std::vector<std::vector<gmx::RVec>> x;
        for(size_t c = 0; c < nconf; c++)
        {
            auto y = displaced(mol.xOriginal(), c);
            for(size_t i = natom1; i < y.size(); i++)
            {
                y[i][XX] += 0.03*c;
            }
            x.push_back(y);
        }
        struct Result
        {
            std::vector<std::map<InteractionType, double>> einter;
            std::vector<double>                            epot;
            std::vector<std::vector<gmx::RVec>>            finter;
            std::vector<std::vector<gmx::RVec>>            x;
        };
        const auto run = [&](const std::vector<size_t> &confs)
        {
            Result r;
            for(auto c : confs)
            {
                r.x.push_back(x[c]);
            }
            forceComp_.computeInteraction(&msghandler_, pd, &mol, &r.x, true,
                                          &r.einter, &r.epot, &r.finter);
            EXPECT_TRUE(msghandler_.ok());
            return r;
        };
        std::vector<size_t> all(nconf), even, oddReversed;
        for(size_t c = 0; c < nconf; c++)
        {
            all[c] = c;
            if (c % 2 == 0)
            {
                even.push_back(c);
            }
            else
            {
                oddReversed.insert(oddReversed.begin(), c);
            }
        }
        auto together = run(all);
        ASSERT_EQ(nconf, together.einter.size());
        for(const auto &subset : { even, oddReversed })
        {
            auto part = run(subset);
            ASSERT_EQ(subset.size(), part.einter.size());
            for(size_t k = 0; k < subset.size(); k++)
            {
                size_t c = subset[k];
                SCOPED_TRACE(gmx::formatString("%s %s dimer %zu", ffName.c_str(), molName.c_str(), c));
                EXPECT_EQ(together.einter[c], part.einter[k]);
                EXPECT_EQ(together.epot[c], part.epot[k]);
                ASSERT_EQ(together.finter[c].size(), part.finter[k].size());
                for(size_t i = 0; i < part.finter[k].size(); i++)
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        EXPECT_EQ(together.finter[c][i][m], part.finter[k][i][m]) << "force on " << i;
                        EXPECT_EQ(together.x[c][i][m], part.x[k][i][m]) << "position of " << i;
                    }
                }
            }
        }
    }

    /*! \brief Compare computeInteraction with the energies of the dimer and the monomers
     * The interaction energy is the energy of the relaxed dimer minus those
     * of the relaxed monomers, each computed with compute, and the
     * interaction forces are the differences of the forces.
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
            std::vector<Topology *> tops = { mol.topologyPtr() };
            for(auto &top : mol.fragmentHandler()->topologiesPtr())
            {
                tops.push_back(&top);
            }
            for(auto top : tops)
            {
                // The monomers may have none, if all pairs are excluded
                if (top->nonBondedPairs().hasType(InteractionType::INDUCTIONCORRECTION))
                {
                    top->setNonBondedParameters(InteractionType::INDUCTIONCORRECTION, correction);
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
        const auto &tops   = mol.fragmentHandler()->topologies();
        const auto &astart = mol.fragmentHandler()->atomStart();
        for(size_t c = 0; c < nconf; c++)
        {
            const auto label = gmx::formatString("%s %s conformation %zu", ffName.c_str(), molName.c_str(), c);
            // The dimer
            auto                              xdimer = x[c];
            std::vector<gmx::RVec>            fref(xdimer.size());
            std::map<InteractionType, double> edimer;
            forceComp_.compute(&msghandler_, pd, mol.topology(), &xdimer, &fref, &edimer);
            double eref = edimer[InteractionType::EPOT];
            EXPECT_NEAR(eref, epot[c], gmx::test::precisionTolerance(1e-7, 1e-5)*std::max(1.0, std::abs(eref))) << label;
            // Minus the monomers
            for(size_t ff = 0; ff < tops.size(); ff++)
            {
                size_t                            natom = tops[ff].atoms().size();
                std::vector<gmx::RVec>            xmon(x[c].begin()+astart[ff], x[c].begin()+astart[ff]+natom);
                std::vector<gmx::RVec>            fmon(natom);
                std::map<InteractionType, double> emon;
                forceComp_.compute(&msghandler_, pd, &tops[ff], &xmon, &fmon, &emon);
                eref -= emon[InteractionType::EPOT];
                for(size_t i = 0; i < natom; i++)
                {
                    rvec_dec(fref[astart[ff]+i], fmon[i]);
                }
            }
            double einterPot = einter[c][InteractionType::EPOT];
            EXPECT_NEAR(eref, einterPot, gmx::test::precisionTolerance(1e-7, 1e-5)*std::max(1.0, std::abs(eref))) << label;
            compareForces(fref, finter[c], 1e-5, label);
            // The terms add up to the interaction energy, except those that
            // are sums of other terms.
            const std::set<InteractionType> sums = { InteractionType::ALLELEC, InteractionType::EXCHIND,
                                               InteractionType::EPOT };
            double esum = 0;
            for(const auto &e : einter[c])
            {
                if (sums.count(e.first) == 0)
                {
                    esum += e.second;
                }
            }
            EXPECT_NEAR(einterPot, esum, gmx::test::precisionTolerance(1e-7, 1e-5)*std::max(1.0, std::abs(einterPot))) << label;
            // The polarization is part of the induction
            EXPECT_TRUE(einter[c].end() == einter[c].find(InteractionType::POLARIZATION)) << label;
            if (!separate)
            {
                EXPECT_TRUE(einter[c].end() == einter[c].find(InteractionType::INDUCTIONCORRECTION)) << label;
            }
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
    testTotal("water_hyperpol2", "water_dimer", 3, ShellOptions());
}

// Test: The maximum distance between shells and cores is an option of the
// polarization in the force field, 0.02 nm if it is not there.
TEST_F(ComputeTotalTest, ShellMaxDistanceOption)
{
    auto pd0 = getForceField("water_hyperpol2");
    ASSERT_NE(nullptr, pd0);
    ForceField pd = *pd0;
    auto       fs = pd.findForces(InteractionType::POLARIZATION);
    ASSERT_FALSE(fs->optionExists(shellMaxDistanceName));
    EXPECT_DOUBLE_EQ(0.02, shellMaxDistance(&pd));
    fs->addOption(shellMaxDistanceName, "0.035");
    EXPECT_DOUBLE_EQ(0.035, shellMaxDistance(&pd));
}

// Test: The shells of a hyperpolarizable dimer converge in a few steps.
// The quartic term of the polarization of oxygen makes the stiffness of the
// spring along the core-shell vector three times larger than perpendicular
// to it, which is taken into account in the step.
TEST_F(ComputeTotalTest, HyperpolarizableConvergence)
{
    auto pd0 = getForceField("water_hyperpol2");
    ASSERT_NE(nullptr, pd0);
    // Without a maximum shell distance
    ForceField pd = *pd0;
    pd.findForces(InteractionType::POLARIZATION)->addOption(shellMaxDistanceName, "0");
    EXPECT_EQ(0, shellMaxDistance(&pd));
    ACTMol mol;
    ASSERT_TRUE(readMolecule(&msghandler_, &pd, "water_dimer", &mol));
    const ForceComputer forceComp(1e-6, 100, 1);
    ConformationSet set;
    set.top = mol.topology();
    set.setCoordinates({ mol.xOriginal() });
    std::vector<std::map<InteractionType, double>> energies;
    forceComp.computeTotal(&msghandler_, &pd, &set, ShellOptions(), &energies);
    // In single precision the forces on the shells, with a charge of
    // -4.9 e, are differences of Coulomb forces of about 1e5 kJ/mol/nm,
    // which are too large for the tolerance
    if (GMX_DOUBLE)
    {
        EXPECT_EQ(0u, forceComp.numShellConvergenceFailed());
        EXPECT_LE(forceComp.numShellIterations(), 8u);
    }
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
    ASSERT_TRUE(mol.topology()->nonBondedPairs().hasType(itype));
    mol.topologyPtr()->setNonBondedParameters(itype, { 2000.0, 1000.0, 20.0 });
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
        const auto label = gmx::formatString("conformation %zu", c);
        auto xref  = x[c];
        std::vector<gmx::RVec>            fref(xref.size());
        std::map<InteractionType, double> eref;
        forceComp_.compute(&msghandler_, pd, top, &xref, &fref, &eref);
        // The correction is booked as such, not as exchange
        ASSERT_TRUE(eref.end() != eref.find(itype));
        EXPECT_NE(0, eref[itype]) << label;
        compareEnergies(eref, energies[c], gmx::test::precisionTolerance(1e-8, 1e-6), label);
        compareForces(fref, set.data.f(c).toRVec(), gmx::test::precisionTolerance(1e-6, 1e-5), label);
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

// Test: Starting the relaxation from relaxed shells gives the same energies,
// with the induction from the shells on their parents, in fewer iterations
TEST_F(ComputeTotalTest, ShellsFromInput)
{
    for (const std::string molName : { "water", "hbr-hf" })
    {
        auto pd = getForceField("ACS-pg-vs3");
        ASSERT_NE(nullptr, pd);
        ACTMol mol;
        ASSERT_TRUE(readMolecule(&msghandler_, pd, molName, &mol));
        auto top = mol.topology();
        auto x   = mol.xOriginal();
        std::vector<gmx::RVec>            fref(x.size()), f(x.size());
        std::map<InteractionType, double> eref, e;
        forceComp_.resetStatistics();
        forceComp_.compute(&msghandler_, pd, top, &x, &fref, &eref);
        auto fromParents = forceComp_.numShellIterations();
        ASSERT_GT(fromParents, 2u) << molName;
        // The shells in x are relaxed now
        forceComp_.resetStatistics();
        forceComp_.compute(&msghandler_, pd, top, &x, &f, &e, { 0, 0, 0 }, true, {}, true);
        EXPECT_TRUE(msghandler_.ok());
        EXPECT_LE(forceComp_.numShellIterations(), 2u) << molName;
        ASSERT_NE(eref.end(), eref.find(InteractionType::INDUCTION));
        compareEnergies(eref, e, gmx::test::precisionTolerance(1e-8, 1e-6), molName + " shells from input");
        // In single precision the relaxed shells are returned in x as
        // floats, which rounds the integer coordinates, see coordinate.h.
        // The relaxation from these shells then stops elsewhere within the
        // shell tolerance, which changes the forces by about 3e-5 of the
        // largest one.
        compareForces(fref, f, gmx::test::precisionTolerance(1e-5, 1e-4), molName + " shells from input");
    }
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
        const auto label = gmx::formatString("conformation %zu", c);
        auto xref  = x[c];
        std::vector<gmx::RVec>            fref(xref.size());
        std::map<InteractionType, double> eref;
        forceComp_.computeOnce(&msghandler_, pd, top, &xref, &fref, &eref, { 0, 0, 0 });
        forceComp_.spreadVsiteForces(top, &xref, &fref);
        // computeOnce gives the potential terms only
        for(const auto &e : eref)
        {
            const auto t = energies[c].find(e.first);
            ASSERT_TRUE(energies[c].end() != t) << label << " " << interactionTypeToString(e.first);
            EXPECT_NEAR(e.second, t->second, gmx::test::precisionTolerance(1e-10, 1e-6)*std::max(1.0, std::abs(e.second)))
                << label << " " << interactionTypeToString(e.first);
        }
        compareForces(fref, set.data.f(c).toRVec(), gmx::test::precisionTolerance(1e-10, 1e-5), label);
        // The shells did not move, the coordinates are rounded to the
        // unit of the integer coordinates, see coordinate.h. The virtual
        // sites are constructed again from the rounded coordinates, in
        // single precision from those relative to their first parent.
        for(size_t i = 0; i < x[c].size(); i++)
        {
            bool vsite = top->atoms()[i].pType() == ActParticle::Vsite;
            for(int m = 0; m < DIM; m++)
            {
                double xr = fromCoordinate(toCoordinate(x[c][i][m]));
                EXPECT_NEAR(xr, static_cast<double>(set.data.x(c)[i][m]),
                            vsite ? gmx::test::precisionTolerance(1e-15, 1e-6)*std::max(1.0, std::abs(xr)) : 0)
                    << label << " particle " << i;
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
                EXPECT_EQ(set.data.x(0)[core][m], set.data.x(0)[i][m]);
            }
        }
    }
    // The relaxed shells give the same energies as compute
    auto xref = x[0];
    std::vector<gmx::RVec>            fref(xref.size());
    std::map<InteractionType, double> eref;
    forceComp_.compute(&msghandler_, pd, top, &xref, &fref, &eref,
                       options.field, options.resetShells, options.relaxShells);
    compareEnergies(eref, energies[0], gmx::test::precisionTolerance(1e-8, 1e-6), "relax some");
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
    // Within the precision of the forces, which in single precision without
    // explicit SIMD are computed from single precision coordinates
    const ForceComputer forceComp(gmx::test::precisionTolerance(1e-12, 1e-4), 200, 1);
    forceComp.computeTotal(&msghandler_, pd, &set, ShellOptions(), &energies);
    EXPECT_EQ(2u, forceComp.numEvaluations());
    EXPECT_GE(forceComp.numShellIterations(), 2u);
    EXPECT_EQ(0u, forceComp.numShellConvergenceFailed());
}

TEST_F(ComputeTotalTest, InteractionPolarizable)
{
    testInteraction("ACS-pg-vs2", "hbr-hf", 3);
}

TEST_F(ComputeTotalTest, InteractionHyperpolarizable)
{
    testInteraction("water_hyperpol2", "water_dimer", 3);
}

TEST_F(ComputeTotalTest, InteractionNotSeparated)
{
    testInteraction("water_hyperpol2", "water_dimer", 3, false);
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

// Test: Dimers computed together or in subsets give the same results, with
// virtual sites and without shells
TEST_F(ComputeTotalTest, InteractionSubsetsVsites)
{
    testInteractionSubsets("ACS-v3fds", "water_dimer", 13);
}

// Test: The same with shells and virtual sites
TEST_F(ComputeTotalTest, InteractionSubsetsPolarizable)
{
    testInteractionSubsets("ACS-pg-vs2", "hbr-hf", 13);
}

// Test: The same with hyperpolarizable shells
TEST_F(ComputeTotalTest, InteractionSubsetsHyperpolarizable)
{
    testInteractionSubsets("water_hyperpol2", "water_dimer", 13);
}

// Test: Computing with integer coordinates gives the same as with
// coordinates in nm, which are rounded to the integers, and returns the
// relaxed shells as integers
TEST_F(ComputeTotalTest, IntegerCoordinates)
{
    auto pd = getForceField("ACS-pg-vs2");
    ASSERT_NE(nullptr, pd);
    ACTMol mol;
    ASSERT_TRUE(readMolecule(&msghandler_, pd, "hbr-hf", &mol));
    auto top  = mol.topology();
    auto xr   = mol.xOriginal();
    Coordinates xi(xr);
    std::vector<gmx::RVec>            fr(xr.size()), fi(xr.size());
    std::map<InteractionType, double> er, ei;
    forceComp_.compute(&msghandler_, pd, top, &xr, &fr, &er);
    forceComp_.compute(&msghandler_, pd, top, &xi, &fi, &ei);
    compareEnergies(er, ei, 0, "integer coordinates");
    compareForces(fr, fi, 0, "integer coordinates");
    ASSERT_EQ(xr.size(), xi.size());
    for(size_t i = 0; i < xr.size(); i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            EXPECT_EQ(xr[i][m], static_cast<real>(xi[i][m])) << "particle " << i;
        }
    }
    // Without relaxation, as computeOnce
    std::map<InteractionType, double> eo;
    forceComp_.computeOnce(&msghandler_, pd, top, &xi, &fi, &eo, { 0, 0, 0 });
    EXPECT_TRUE(eo.end() != eo.find(InteractionType::EPOT));
}

// Test: Conformations computed in double precision, see
// ConformationData::doublePrecision, have the same energies and forces
// within the precision of real, with shells and virtual sites and without.
// In a double precision build they are the same arrays.
TEST_F(ComputeTotalTest, DoublePrecisionArrays)
{
    for (const auto &ffMol : { std::pair<std::string, std::string>{ "ACS-pg-vs2", "hbr-hf" },
                               std::pair<std::string, std::string>{ "ACS-v3fds", "water_dimer" } })
    {
        auto pd = getForceField(ffMol.first);
        ASSERT_NE(nullptr, pd);
        ACTMol mol;
        ASSERT_TRUE(readMolecule(&msghandler_, pd, ffMol.second, &mol));
        std::vector<std::vector<gmx::RVec>> x;
        for(size_t c = 0; c < 9; c++)
        {
            x.push_back(displaced(mol.xOriginal(), c));
        }
        ConformationSet set[2];
        std::vector<std::map<InteractionType, double>> energies[2];
        for(int k = 0; k < 2; k++)
        {
            set[k].top = mol.topology();
            set[k].setCoordinates(x);
            set[k].data.doublePrecision = (k == 1);
            forceComp_.computeTotal(&msghandler_, pd, &set[k], ShellOptions(), &energies[k]);
            ASSERT_EQ(x.size(), energies[k].size());
        }
        // The arrays hold the conformations of the last shell iteration
        EXPECT_EQ(!GMX_DOUBLE, set[1].data.arraysDouble.nconf > 0);
        for(size_t c = 0; c < x.size(); c++)
        {
            const auto label = gmx::formatString("%s %s conformation %zu", ffMol.first.c_str(),
                                                 ffMol.second.c_str(), c);
            compareEnergies(energies[0][c], energies[1][c], gmx::test::precisionTolerance(1e-14, 1e-6), label);
            // As in the comparison of the explicit SIMD in single and double precision
            compareForces(set[0].data.f(c).toRVec(), set[1].data.f(c).toRVec(),
                          gmx::test::precisionTolerance(1e-14, 1e-4), label);
        }
    }
}

} // namespace

} // namespace alexandria
