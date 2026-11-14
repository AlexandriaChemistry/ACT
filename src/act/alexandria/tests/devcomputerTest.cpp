/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
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
 * Unit tests for DevComputer classes.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include "act/alexandria/devcomputer.h"
#include "act/alexandria/molgen.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/potential.h"
#include "act/forcefield/generate_dependent.h"
#include "act/forces/conformationdata.h"
#include "act/gpu/gpu.h"
#include "act/qgen/qgen_resp.h"
#include "act/qgen/qtype.h"
#include "act/utility/timings.h"
#include "act/forces/simd_information.h"
#include "act/molprop/molprop_xml.h"
#include "act/topology/actmol.h"

#include "testutils/cmdlinetest.h"
#include "testutils/gputests.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

//! Total weight and chi squared, keyed by molecule name and target
typedef std::map<std::string, std::pair<double, double>> DevResults;

//! Options for running the ForceEnergyDevComputer
struct DevOptions
{
    //! Weight for each fitting target
    std::map<eRMS, double> weights;
    //! Boltzmann temperature for each target
    std::map<eRMS, double> boltzmann;
    //! The loss function to use
    LossFunction           lossFunction                = LossFunction::MSE;
    //! Whether to split off induction correction from induction
    bool                   separateInductionCorrection = true;
    //! Relative tolerance for comparing to reference data
    double                 tolerance                   = 1e-7;
    //! Also check the van der Waals energies and forces of each conformation
    bool                   checkConformations          = false;
    //! Only use these molecules, all if empty
    std::set<std::string>  molecules;
};

//! \return pointers to the molecules, as the DevComputers take them
static std::vector<ACTMol *> pointers(std::vector<ACTMol> *mols)
{
    std::vector<ACTMol *> ptrs;
    for(auto &mol : *mols)
    {
        ptrs.push_back(&mol);
    }
    return ptrs;
}

//! \return key for the DevResults map
static std::string resultKey(const std::string &molname, eRMS erms)
{
    return gmx::formatString("%s-%s", molname.c_str(), rmsName(erms));
}

//! Check that two sets of results are identical
static void expectSameResults(const DevResults &a, const DevResults &b)
{
    EXPECT_EQ(a.size(), b.size());
    for(const auto &ra : a)
    {
        const auto rb = b.find(ra.first);
        ASSERT_TRUE(b.end() != rb) << ra.first;
        EXPECT_DOUBLE_EQ(ra.second.first, rb->second.first) << ra.first << " weight";
        EXPECT_DOUBLE_EQ(ra.second.second, rb->second.second) << ra.first << " chi2";
    }
}

//! \return a weight of one for each target
static std::map<eRMS, double> unitWeights(const std::vector<eRMS> &erms)
{
    std::map<eRMS, double> weights;
    for(const auto e : erms)
    {
        weights.insert({e, 1.0});
    }
    return weights;
}

//! \return the same temperature for each target
static std::map<eRMS, double> sameTemperature(const std::vector<eRMS> &erms, double T)
{
    std::map<eRMS, double> boltzmann;
    for(const auto e : erms)
    {
        boltzmann.insert({e, T});
    }
    return boltzmann;
}

//! Targets for monomers
static const std::vector<eRMS> monomerTargets = { eRMS::EPOT, eRMS::Force2 };

//! Targets for dimers
static const std::vector<eRMS> dimerTargets = { eRMS::Interaction, eRMS::Electrostatics,
                                                eRMS::Dispersion, eRMS::Exchange };

//! Targets for dimers including the induction terms
static const std::vector<eRMS> inductionTargets = { eRMS::Interaction, eRMS::Electrostatics,
                                                    eRMS::Dispersion, eRMS::Exchange,
                                                    eRMS::Induction, eRMS::DeltaHF };

// Test fixture for ForceEnergyDevComputer
class ForceEnergyDevComputerTest : public gmx::test::CommandLineTestBase
{
protected:
    ForceEnergyDevComputerTest() = default;

    //! \return the monomer test file
    std::string monomerFile() { return fileManager().getInputFilePath("molprop/molprop.xml"); }

    //! \return the dimer test file
    std::string dimerFile() { return fileManager().getInputFilePath("molprop/mpdimer.xml"); }

    //! \return the dimer test file with induction correction
    std::string deltaHFFile() { return fileManager().getInputFilePath("molprop/methanol-ethane-deltaHF.xml"); }

    //! Run the test with each Boltzmann key as a target with weight 1
    void testCalcDev(const std::string      &mpFile,
                     const std::string      &ffFile,
                     std::map<eRMS, double>  boltzmann)
    {
        DevOptions opts;
        opts.boltzmann = boltzmann;
        for(const auto &erms : boltzmann)
        {
            opts.weights.insert({erms.first, 1.0});
        }
        testCalcDev(mpFile, ffFile, opts);
    }

    //! Run the test and compare to reference data
    void testCalcDev(const std::string &mpFile,
                     const std::string &ffFile,
                     const DevOptions  &opts)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, opts.tolerance);
        checker_.setDefaultTolerance(tolerance);
        computeDev(mpFile, ffFile, opts, &checker_);
    }

    /*! \brief Run the ForceEnergyDevComputer on all molecules in a file
     * Targets are shared between molecules, so the results accumulate.
     * \param[in] mpFile  Molprop file
     * \param[in] ffFile  Force field name
     * \param[in] opts    The options
     * \param[in] checker Reference data checker, may be nullptr
     * \return The weight and chi squared after each molecule
     */
    DevResults computeDev(const std::string               &mpFile,
                          const std::string               &ffFile,
                          const DevOptions                &opts,
                          gmx::test::TestReferenceChecker *checker)
    {
        const ForceComputer     forceComputer;
        MsgHandler              msghandler;
        // Use next statement for debugging only
        // msghandler.setPrintLevel(ACTStatus::Debug);
        //! FF
        auto pd = getForceField(ffFile);

        //! Vector of molecule properties
        std::vector<alexandria::MolProp>  mps;
        MolPropRead(&msghandler, mpFile.c_str(), &mps);

        ForceEnergyDevComputer  computer(opts.boltzmann);
        computer.setSeparateInductionCorrection(opts.separateInductionCorrection);
        std::map<eRMS, FittingTarget> targets;
        for(const auto &erms : opts.weights)
        {
            targets.insert({erms.first,
                    FittingTarget(erms.first, iMolSelect::Train)});
            targets.find(erms.first)->second.setWeight(erms.second);
        };

        DevResults results;
        for(auto &mp : mps)
        {
            if (!opts.molecules.empty() && opts.molecules.count(mp.getMolname()) == 0)
            {
                continue;
            }
            std::vector<ACTMol> actmols(1);
            actmols[0].Merge(&mp);
            actmols[0].GenerateTopology(&msghandler, pd, missingParameters::Ignore);
            EXPECT_TRUE(msghandler.ok());
            if (!msghandler.ok())
            {
                continue;
            }
            auto &atoms = actmols[0].atomsConst();
            std::vector<gmx::RVec> forces(actmols[0].atomsConst().size());
            std::vector<gmx::RVec> coords = actmols[0].xOriginal();
            actmols[0].generateCharges(&msghandler, pd, &forceComputer,
                                       pd->chargeGenerationAlgorithm(), &coords, &forces);
            std::vector<double> qq;
            for(size_t j = 0; j < atoms.size(); j++)
            {
                qq.push_back(atoms[j].charge());
            }
            if (checker)
            {
                std::string label(actmols[0].getMolname());
                label += "-charge";
                setSingleTolerance(checker, 1e-6, 1e-6);
                checker->checkSequence(qq.begin(), qq.end(), label.c_str());
            }
            computer.calcDeviation(&msghandler,
                                   &forceComputer,
                                   pointers(&actmols),
                                   &targets,
                                   pd,
                                   opts.lossFunction);
            if (checker && opts.checkConformations)
            {
                checkConformations(actmols[0], checker);
            }
            for(const auto &t : targets)
            {
                results[resultKey(actmols[0].getMolname(), t.first)] =
                    { t.second.totalWeight(), t.second.chiSquared() };
                if (checker)
                {
                    const std::string weight = gmx::formatString("%s-%s-weight",
                                                                 actmols[0].getMolname().c_str(), rmsName(t.first));
                    checker->checkDouble(t.second.totalWeight(), weight.c_str());
                    const std::string chi2 = gmx::formatString("%s-%s-chi2",
                                                               actmols[0].getMolname().c_str(), rmsName(t.first));
                    // The deviations are differences of energies or forces
                    // that are much larger, with errors relative to those.
                    // The interaction energies of dimers have errors of up
                    // to about 1e-3 kJ/mol, which change a chi2 by twice the
                    // deviations times that, also when the chi2 is small.
                    setSingleTolerance(checker, 1e-4, 1e-3);
                    checker->checkDouble(t.second.chiSquared(), chi2.c_str());
                }
            }
        }
        return results;
    }

    /*! \brief Check that conformations kept their coordinates and reference data
     * \param[in] before The conformations before
     * \param[in] after  The conformations after
     */
    static void checkConformationsKept(const ConformationData &before,
                                       const ConformationData &after)
    {
        ASSERT_EQ(before.nconf(), after.nconf());
        for(size_t c = 0; c < before.nconf(); c++)
        {
            const auto &cb = before.confs[c];
            const auto &ca = after.confs[c];
            EXPECT_EQ(cb.dataFile, ca.dataFile) << "conformation " << c;
            EXPECT_EQ(cb.eMin, ca.eMin) << "conformation " << c;
            const auto xb = before.x(c);
            const auto xa = after.x(c);
            ASSERT_EQ(xb.size(), xa.size()) << "conformation " << c;
            for(size_t i = 0; i < xb.size(); i++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    EXPECT_EQ(xb[i][m], xa[i][m]) << "conformation " << c;
                }
            }
            ASSERT_EQ(cb.eMap.size(), ca.eMap.size()) << "conformation " << c;
            for(const auto &e : cb.eMap)
            {
                const auto a = ca.eMap.find(e.first);
                ASSERT_TRUE(a != ca.eMap.end()) << "conformation " << c;
                EXPECT_EQ(e.second.haveQM(), a->second.haveQM()) << "conformation " << c;
                if (e.second.haveQM())
                {
                    EXPECT_EQ(e.second.eqm(), a->second.eqm()) << "conformation " << c;
                }
            }
        }
    }

    /*! \brief Check the QM reference for induction against the molprop data
     * \param[in] mpFile  Molprop file with dimers
     * \param[in] sepIndCorr Whether to separate the induction correction
     * \return The number of energies with a non-zero induction correction
     */
    int checkInductionQM(const std::string &mpFile, bool sepIndCorr)
    {
        const ForceComputer forceComputer;
        MsgHandler    msghandler;
        auto          pd = getForceField("ACS-pol");
        std::vector<alexandria::MolProp> mps;
        MolPropRead(&msghandler, mpFile.c_str(), &mps);
        int nonZeroDHF = 0;
        for(auto &mp : mps)
        {
            ACTMol actmol;
            actmol.Merge(&mp);
            actmol.GenerateTopology(&msghandler, pd, missingParameters::Ignore);
            EXPECT_TRUE(msghandler.ok());
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            std::vector<gmx::RVec> coords = actmol.xOriginal();
            actmol.generateCharges(&msghandler, pd, &forceComputer,
                                   pd->chargeGenerationAlgorithm(), &coords, &forces);
            // QM induction and induction correction for each experiment
            std::map<int, std::pair<double, double>> qmInduction;
            // Experiments that have an induction correction
            std::set<int>                            haveDHF;
            for(const auto &exper : actmol.experimentConst())
            {
                if (exper.hasProperty(MolPropObservable::INDUCTION))
                {
                    double dhf = 0;
                    if (exper.hasProperty(MolPropObservable::INDUCTIONCORRECTION))
                    {
                        dhf = exper.propertyConst(MolPropObservable::INDUCTIONCORRECTION)[0]->getValue();
                        haveDHF.insert(exper.id());
                    }
                    qmInduction.insert({exper.id(),
                            { exper.propertyConst(MolPropObservable::INDUCTION)[0]->getValue(), dhf }});
                }
            }
            std::vector<std::vector<std::pair<double, double>>>               forceMap;
            std::vector<ACTEnergy>                                            energyMap;
            ACTEnergyMapVector                                                interactionEnergyMap;
            std::vector<std::pair<double, std::map<InteractionType, double>>> energyComponentMap;
            // The conformations hold the reference data for training
            const auto before = *actmol.conformationData(2);
            EXPECT_GT(before.nconf(), 0u);
            actmol.forceEnergyMaps(&msghandler, pd, &forceComputer, &forceMap, &energyMap,
                                   &interactionEnergyMap, &energyComponentMap, sepIndCorr);
            checkConformationsKept(before, *actmol.conformationData(2));
            EXPECT_EQ(qmInduction.size(), interactionEnergyMap.size());
            for(const auto &iem : interactionEnergyMap)
            {
                const auto ind = iem.find(InteractionType::INDUCTION);
                if (iem.end() == ind || !ind->second.haveQM())
                {
                    ADD_FAILURE() << "No QM induction for " << actmol.getMolname();
                    continue;
                }
                const auto qm = qmInduction.find(ind->second.id());
                if (qmInduction.end() == qm)
                {
                    ADD_FAILURE() << "No experiment " << ind->second.id() << " for " << actmol.getMolname();
                    continue;
                }
                double qmInd = qm->second.first;
                double qmDHF = qm->second.second;
                if (qmDHF != 0)
                {
                    nonZeroDHF += 1;
                }
                const auto icorr = iem.find(InteractionType::INDUCTIONCORRECTION);
                if (sepIndCorr)
                {
                    // Induction correction is a separate term
                    EXPECT_DOUBLE_EQ(qmInd, ind->second.eqm());
                    if (haveDHF.count(ind->second.id()) > 0)
                    {
                        EXPECT_TRUE(iem.end() != icorr);
                        if (iem.end() != icorr)
                        {
                            EXPECT_DOUBLE_EQ(qmDHF, icorr->second.eqm());
                        }
                    }
                    else if (iem.end() != icorr)
                    {
                        // Without data there should not be a QM reference
                        EXPECT_FALSE(icorr->second.haveQM());
                    }
                }
                else
                {
                    // Induction correction is added to the induction
                    EXPECT_DOUBLE_EQ(qmInd + qmDHF, ind->second.eqm());
                    EXPECT_TRUE(iem.end() == icorr);
                }
            }
        }
        return nonZeroDHF;
    }

    /*! \brief Check the van der Waals energies and forces of each conformation
     * \param[in] actmol  The molecule after calcDeviation
     * \param[in] checker Reference data checker
     */
    /*! \brief Set a tolerance in single precision only
     * In double precision the default tolerance of the checker is kept.
     * \param[inout] checker  The checker
     * \param[in]    absolute Absolute tolerance
     * \param[in]    relative Relative tolerance
     */
    static void setSingleTolerance(gmx::test::TestReferenceChecker *checker,
                                   double                           absolute,
                                   double                           relative)
    {
        if (!GMX_DOUBLE)
        {
            checker->setDefaultTolerance(gmx::test::FloatingPointTolerance(absolute, absolute, relative, relative,
                                                                           UINT64_MAX, UINT64_MAX, false));
        }
    }

    static void checkConformations(ACTMol                          &actmol,
                                   gmx::test::TestReferenceChecker *checker)
    {
        const auto *cdata = actmol.conformationData(actmol.isDimer() ? 2 : 0);
        EXPECT_GT(cdata->nconf(), 0u) << actmol.getMolname();
        for(size_t i = 0; i < cdata->nconf(); i++)
        {
            auto        ener  = cdata->energies.energies(i);
            const auto  label = gmx::formatString("%s-conf%zu", actmol.getMolname().c_str(), i);
            // In single precision relative to the largest energy and force
            double emax = 0;
            for(const auto itype : { InteractionType::EXCHANGE,
                                     InteractionType::DISPERSION,
                                     InteractionType::EPOT })
            {
                const auto ee = ener.find(itype);
                if (ener.end() != ee)
                {
                    emax = std::max(emax, std::abs(ee->second));
                }
            }
            setSingleTolerance(checker, 1e-5*emax, 1e-6);
            for(const auto itype : { InteractionType::EXCHANGE,
                                     InteractionType::DISPERSION,
                                     InteractionType::EPOT })
            {
                const auto ee = ener.find(itype);
                ASSERT_TRUE(ener.end() != ee)
                    << label << " " << interactionTypeToString(itype);
                checker->checkDouble(ee->second,
                                     gmx::formatString("%s-%s", label.c_str(),
                                                       interactionTypeToString(itype).c_str()).c_str());
            }
            // Repulsion and attraction
            EXPECT_GT(ener[InteractionType::EXCHANGE], 0) << label;
            EXPECT_LT(ener[InteractionType::DISPERSION], 0) << label;
            std::vector<double> f;
            for(const auto &fi : cdata->f(i).toRVec())
            {
                for(int m = 0; m < DIM; m++)
                {
                    EXPECT_TRUE(std::isfinite(fi[m])) << label;
                    f.push_back(fi[m]);
                }
            }
            // The coordinates are rounded, which changes the forces by up
            // to the stiffness of the bonds times about 3e-8 nm
            double fmax = 0;
            for(auto fm : f)
            {
                fmax = std::max(fmax, std::abs(fm));
            }
            setSingleTolerance(checker, 5e-5*fmax, 1e-6);
            checker->checkSequence(f.begin(), f.end(), (label + "-forces").c_str());
        }
    }

    //! Run monomer and dimer tests with a loss function
    void testLossFunction(LossFunction lf)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-11);
        checker_.setDefaultTolerance(tolerance);
        DevOptions mono;
        mono.weights      = unitWeights(monomerTargets);
        mono.lossFunction = lf;
        computeDev(monomerFile(), "ACS-g2", mono, &checker_);
        DevOptions dim;
        dim.weights      = unitWeights(dimerTargets);
        dim.lossFunction = lf;
        computeDev(dimerFile(), "ACS-g2", dim, &checker_);
    }
};

// Test: Constructor initializes with correct name
TEST_F(ForceEnergyDevComputerTest, ConstructorSetsCorrectName)
{
    std::map<eRMS, double> boltzmann;
    const ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Constructor with empty Boltzmann map
TEST_F(ForceEnergyDevComputerTest, ConstructorWithEmptyBoltzmannMap)
{
    std::map<eRMS, double> boltzmann;
    const ForceEnergyDevComputer computer(boltzmann);
    EXPECT_TRUE(computer.separateInductionCorrection());
}

// Test: Constructor with non-empty Boltzmann map
TEST_F(ForceEnergyDevComputerTest, ConstructorWithNonEmptyBoltzmannMap)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 298.15;
    boltzmann[eRMS::Force2] = 0.0;
    
    const ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
    EXPECT_TRUE(computer.separateInductionCorrection());
}

// Test: setSeparateInductionCorrection and separateInductionCorrection methods
TEST_F(ForceEnergyDevComputerTest, SetAndGetSeparateInductionCorrection)
{
    std::map<eRMS, double> boltzmann;
    ForceEnergyDevComputer computer(boltzmann);
    
    // Default should be true (from header initialization)
    EXPECT_TRUE(computer.separateInductionCorrection());
    
    // Set to false
    computer.setSeparateInductionCorrection(false);
    EXPECT_FALSE(computer.separateInductionCorrection());
    
    // Set back to true
    computer.setSeparateInductionCorrection(true);
    EXPECT_TRUE(computer.separateInductionCorrection());
}

// Test: Multiple flag toggles
TEST_F(ForceEnergyDevComputerTest, MultipleSeparateInductionCorrectionToggles)
{
    std::map<eRMS, double> boltzmann;
    ForceEnergyDevComputer computer(boltzmann);
    
    for (int i = 0; i < 5; ++i)
    {
        computer.setSeparateInductionCorrection(false);
        EXPECT_FALSE(computer.separateInductionCorrection());
        
        computer.setSeparateInductionCorrection(true);
        EXPECT_TRUE(computer.separateInductionCorrection());
    }
}

// Test: Boltzmann temperature with EPOT
TEST_F(ForceEnergyDevComputerTest, BoltzmannTemperatureEPOT)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 298.15;
    
    const ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Boltzmann temperature with Force2
TEST_F(ForceEnergyDevComputerTest, BoltzmannTemperatureForce2)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::Force2] = 300.0;
    
    const ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Boltzmann temperature with multiple eRMS terms
TEST_F(ForceEnergyDevComputerTest, BoltzmannTemperatureMultipleTerms)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 298.15;
    boltzmann[eRMS::Force2] = 0.0;
    boltzmann[eRMS::Interaction] = 500.0;
    
    const ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
    EXPECT_TRUE(computer.separateInductionCorrection());
}

// Test: Inheritance from DevComputer
TEST_F(ForceEnergyDevComputerTest, InheritsFromDevComputer)
{
    std::map<eRMS, double> boltzmann;
    const ForceEnergyDevComputer computer(boltzmann);
    
    // Verify that DevComputer interface is accessible
    const std::string& name = computer.name();
    EXPECT_EQ("ForceEnergy", name);
}

// Test: Zero Boltzmann temperature
TEST_F(ForceEnergyDevComputerTest, ZeroBoltzmannTemperature)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 0.0;
    
    const ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Large Boltzmann temperature values
TEST_F(ForceEnergyDevComputerTest, LargeBoltzmannTemperatureValues)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 1e6;
    boltzmann[eRMS::Force2] = 1e10;
    
    const ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Negative Boltzmann temperature (edge case)
TEST_F(ForceEnergyDevComputerTest, NegativeBoltzmannTemperature)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = -100.0;
    
    const ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: CalcDeviation for monomer without Boltzmann weighting
TEST_F(ForceEnergyDevComputerTest, CalcDevMonomerNoBoltz)
{
    std::map<eRMS, double> boltzmann;
    for(const auto erms : { eRMS::EPOT,
                            eRMS::Force2 })
    {
        boltzmann.insert({erms, 0.0});
    }
    const std::string mpFile = fileManager().getInputFilePath("molprop/molprop.xml");
    testCalcDev(mpFile, "ACS-g2", boltzmann);
}

// Test: Buffered 14-7 van der Waals energies and forces (force field ACS-v3fds)
TEST_F(ForceEnergyDevComputerTest, CalcDevMonomerLJ14_7)
{
    DevOptions opts;
    // Forces are only computed if they are fitted
    opts.weights            = unitWeights({ eRMS::EPOT, eRMS::Force2 });
    opts.checkConformations = true;
    // The force field has all parameters for these molecules, water has a
    // VSITE3FD virtual site.
    opts.molecules          = { "1-butanol", "water" };
    testCalcDev(monomerFile(), "ACS-v3fds", opts);
}

// Test: Forces of a molecule with a virtual site are compared atom by atom
TEST_F(ForceEnergyDevComputerTest, CalcDevMonomerForcesVsite)
{
    DevOptions opts;
    opts.weights = unitWeights({ eRMS::EPOT, eRMS::Force2 });
    // Water with QM forces for four conformations and a VSITE3FD virtual site
    gmx::test::TestReferenceChecker checker(this->rootChecker());
    checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, opts.tolerance));
    auto results = computeDev(fileManager().getInputFilePath("molprop/water-forces.xml"),
                              "ACS-v3fds", opts, &checker);
    // Three atoms with three force components for each conformation
    EXPECT_EQ(4*3*DIM, results[resultKey("water", eRMS::Force2)].first);
}

// Test: CalcDeviation for dimer without Boltzmann weighting
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerNoBoltz)
{
    std::map<eRMS, double> boltzmann;
    for(const auto erms : { eRMS::Interaction,
                            eRMS::Electrostatics,
                            eRMS::Dispersion,
                            eRMS::Exchange })
    {
        boltzmann.insert({erms, 0.0});
    }
    const std::string mpFile = fileManager().getInputFilePath("molprop/mpdimer.xml");
    testCalcDev(mpFile, "ACS-g2", boltzmann);
}

// Test: CalcDeviation for monomer with Boltzmann weighting
TEST_F(ForceEnergyDevComputerTest, CalcDevMonomerBoltz)
{
    std::map<eRMS, double> boltzmann;
    for(const auto erms : { eRMS::EPOT,
                            eRMS::Force2 })
    {
        boltzmann.insert({erms, 500.0});
    }
    const std::string mpFile = fileManager().getInputFilePath("molprop/molprop.xml");
    testCalcDev(mpFile, "ACS-g2", boltzmann);
}

// Test: CalcDeviation for dimer with Boltzmann weighting
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerBoltz)
{
    std::map<eRMS, double> boltzmann;
    for(const auto erms : { eRMS::Interaction,
                            eRMS::Electrostatics, 
                            eRMS::Dispersion,
                            eRMS::Exchange })
    {
        boltzmann.insert({erms, 500.0});
    }
    const std::string mpFile = fileManager().getInputFilePath("molprop/mpdimer.xml");
    testCalcDev(mpFile, "ACS-g2", boltzmann);
}

// Test: CalcDeviation for polarizable dimer with Boltzmann weighting
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerBoltzPol)
{
    std::map<eRMS, double> boltzmann;
    for(const auto erms : { eRMS::Interaction,
                            eRMS::Electrostatics, 
                            eRMS::Dispersion,
                            eRMS::Exchange })
    {
        boltzmann.insert({erms, 500.0});
    }
    const std::string mpFile = fileManager().getInputFilePath("molprop/mpdimer.xml");
    testCalcDev(mpFile, "ACS-pol", boltzmann);
}

// Test: CalcDeviation for polarizable dimer without Boltzmann weighting
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerNoBoltzPol)
{
    std::map<eRMS, double> boltzmann;
    for(const auto erms : { eRMS::Interaction,
                            eRMS::Electrostatics, 
                            eRMS::Dispersion,
                            eRMS::Exchange })
    {
        boltzmann.insert({erms, 0.0});
    }
    const std::string mpFile = fileManager().getInputFilePath("molprop/mpdimer.xml");
    testCalcDev(mpFile, "ACS-pol", boltzmann);
}

// Test: CalcDeviation for dimer with Induction Correction weighting
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerBoltzInductionCorrection)
{
    std::map<eRMS, double> boltzmann;
    for(const auto erms : { eRMS::Interaction,
                            eRMS::Electrostatics, 
                            eRMS::Dispersion,
                            eRMS::Exchange,
                            eRMS::Induction,
                            eRMS::DeltaHF } )
    {
        boltzmann.insert({erms, 500.0});
    }
    const std::string mpFile = fileManager().getInputFilePath("molprop/methanol-ethane-deltaHF.xml");
    testCalcDev(mpFile, "ACS-pol", boltzmann);
}

// Test: CalcDeviation for dimer with Induction Correction weighting and No Boltzmann
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerNoBoltzInductionCorrection)
{
    std::map<eRMS, double> boltzmann;
    for(const auto erms : { eRMS::Interaction,
                            eRMS::Electrostatics, 
                            eRMS::Dispersion,
                            eRMS::Exchange,
                            eRMS::Induction,
                            eRMS::DeltaHF } )
    {
        boltzmann.insert({erms, 0.0});
    }
    const std::string mpFile = fileManager().getInputFilePath("molprop/methanol-ethane-deltaHF.xml");
    testCalcDev(mpFile, "ACS-pol", boltzmann);
}

// Test: Induction correction is summed into induction rather than separated
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerInductionCorrectionNotSeparated)
{
    DevOptions opts;
    opts.weights                     = unitWeights(inductionTargets);
    opts.separateInductionCorrection = false;
    testCalcDev(deltaHFFile(), "ACS-pol", opts);
}

// Test: Induction target without DeltaHF, as set up by ACMFitnessComputer
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerInductionWithoutDeltaHF)
{
    DevOptions opts;
    opts.weights                     = unitWeights({ eRMS::Interaction,
                                                     eRMS::Electrostatics,
                                                     eRMS::Induction });
    opts.separateInductionCorrection = false;
    testCalcDev(deltaHFFile(), "ACS-pol", opts);
}

// Test: The QM induction includes the induction correction if it is not separated
TEST_F(ForceEnergyDevComputerTest, QMInductionIncludesCorrectionNotSeparated)
{
    // Make sure the test data actually has non-zero induction corrections
    EXPECT_GT(checkInductionQM(deltaHFFile(), false), 0);
}

// Test: The QM induction and induction correction are kept apart if separated
TEST_F(ForceEnergyDevComputerTest, QMInductionSeparateCorrection)
{
    EXPECT_GT(checkInductionQM(deltaHFFile(), true), 0);
}

// Test: Without induction correction data the QM induction is unchanged
TEST_F(ForceEnergyDevComputerTest, QMInductionWithoutCorrectionData)
{
    EXPECT_EQ(0, checkInductionQM(dimerFile(), false));
}

// Test: Separating the induction correction without induction correction data
TEST_F(ForceEnergyDevComputerTest, QMInductionSeparateWithoutCorrectionData)
{
    EXPECT_EQ(0, checkInductionQM(dimerFile(), true));
}

// Test: ExchInd and AllElec targets for a polarizable dimer
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerExchIndAllElecPol)
{
    DevOptions opts;
    opts.weights = unitWeights({ eRMS::Interaction, eRMS::ExchInd, eRMS::AllElec });
    testCalcDev(dimerFile(), "ACS-pol", opts);
}

// Test: ExchInd and AllElec without induction correction data. These
// combine three QM terms, and one missing term is only allowed if the
// induction correction is not separated.
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerExchIndAllElecPolNotSeparated)
{
    DevOptions opts;
    opts.weights                     = unitWeights({ eRMS::Interaction, eRMS::ExchInd, eRMS::AllElec });
    opts.separateInductionCorrection = false;
    // ACS-pol has no exchange for these dimers, so ExchInd is just the induction,
    // which is about 1e-6 kJ/mol here and therefore depends on the convergence
    // of the shell positions.
    opts.tolerance                   = 1e-6;
    testCalcDev(dimerFile(), "ACS-pol", opts);
}

// Test: ExchInd and AllElec targets with induction correction data
TEST_F(ForceEnergyDevComputerTest, CalcDevDimerExchIndAllElecInductionCorrection)
{
    DevOptions opts;
    opts.weights = unitWeights({ eRMS::Interaction, eRMS::ExchInd, eRMS::AllElec,
                                 eRMS::Induction, eRMS::DeltaHF });
    testCalcDev(deltaHFFile(), "ACS-pol", opts);
}

// Test: Mean absolute error loss function
TEST_F(ForceEnergyDevComputerTest, CalcDevLossMAE)
{
    testLossFunction(LossFunction::MAE);
}

// Test: Huber loss function
TEST_F(ForceEnergyDevComputerTest, CalcDevLossHuber)
{
    testLossFunction(LossFunction::Huber);
}

// Test: Asinh loss function
TEST_F(ForceEnergyDevComputerTest, CalcDevLossAsinh)
{
    testLossFunction(LossFunction::Asinh);
}

// Test: A monomer target with weight zero is ignored
TEST_F(ForceEnergyDevComputerTest, ZeroWeightTargetIgnoredMonomer)
{
    DevOptions both;
    both.weights = { { eRMS::EPOT, 1.0 }, { eRMS::Force2, 0.0 } };
    auto rboth   = computeDev(monomerFile(), "ACS-g2", both, nullptr);
    DevOptions epot;
    epot.weights = { { eRMS::EPOT, 1.0 } };
    const auto repot   = computeDev(monomerFile(), "ACS-g2", epot, nullptr);
    DevResults rforce;
    for(const auto &r : rboth)
    {
        if (r.first.find(rmsName(eRMS::Force2)) != std::string::npos)
        {
            rforce.insert(r);
        }
    }
    ASSERT_FALSE(rforce.empty());
    for(const auto &r : rforce)
    {
        EXPECT_EQ(0, r.second.first) << r.first;
        EXPECT_EQ(0, r.second.second) << r.first;
        rboth.erase(r.first);
    }
    expectSameResults(repot, rboth);
}

// Test: A dimer component target with weight zero is ignored
TEST_F(ForceEnergyDevComputerTest, ZeroWeightTargetIgnoredDimer)
{
    DevOptions both;
    both.weights = { { eRMS::Interaction, 1.0 }, { eRMS::Electrostatics, 0.0 } };
    auto rboth   = computeDev(dimerFile(), "ACS-g2", both, nullptr);
    DevOptions inter;
    inter.weights = { { eRMS::Interaction, 1.0 } };
    const auto rinter   = computeDev(dimerFile(), "ACS-g2", inter, nullptr);
    DevResults relec;
    for(const auto &r : rboth)
    {
        if (r.first.find(rmsName(eRMS::Electrostatics)) != std::string::npos)
        {
            relec.insert(r);
        }
    }
    ASSERT_FALSE(relec.empty());
    for(const auto &r : relec)
    {
        EXPECT_EQ(0, r.second.first) << r.first;
        EXPECT_EQ(0, r.second.second) << r.first;
        rboth.erase(r.first);
    }
    expectSameResults(rinter, rboth);
}

// Test: Nothing is computed if there are no force or energy targets
TEST_F(ForceEnergyDevComputerTest, NoForceEnergyTargets)
{
    DevOptions opts;
    opts.weights = unitWeights({ eRMS::ESP, eRMS::MU });
    for(const auto &mpFile : { monomerFile(), dimerFile() })
    {
        const auto results = computeDev(mpFile, "ACS-g2", opts, nullptr);
        EXPECT_FALSE(results.empty());
        for(const auto &r : results)
        {
            EXPECT_EQ(0, r.second.first) << r.first;
            EXPECT_EQ(0, r.second.second) << r.first;
        }
    }
}

// Test: A negative temperature means no Boltzmann weighting
TEST_F(ForceEnergyDevComputerTest, NegativeTemperatureEqualsNoBoltzmann)
{
    for(const auto &setup : { std::make_pair(monomerFile(), monomerTargets),
                              std::make_pair(dimerFile(), dimerTargets) })
    {
        DevOptions noBoltz;
        noBoltz.weights   = unitWeights(setup.second);
        noBoltz.boltzmann = sameTemperature(setup.second, 0.0);
        DevOptions negative;
        negative.weights   = noBoltz.weights;
        negative.boltzmann = sameTemperature(setup.second, -300.0);
        expectSameResults(computeDev(setup.first, "ACS-g2", noBoltz, nullptr),
                          computeDev(setup.first, "ACS-g2", negative, nullptr));
    }
}

// Test: Only the Interaction temperature is used for weighting dimer components
TEST_F(ForceEnergyDevComputerTest, ComponentTemperatureIgnored)
{
    DevOptions noBoltz;
    noBoltz.weights = unitWeights(dimerTargets);
    DevOptions components;
    components.weights   = noBoltz.weights;
    components.boltzmann = sameTemperature({ eRMS::Electrostatics,
                                             eRMS::Dispersion,
                                             eRMS::Exchange }, 500.0);
    expectSameResults(computeDev(dimerFile(), "ACS-g2", noBoltz, nullptr),
                      computeDev(dimerFile(), "ACS-g2", components, nullptr));
}

// Test: At very low temperature weights underflow but results stay finite
TEST_F(ForceEnergyDevComputerTest, VeryLowTemperatureMonomer)
{
    DevOptions noBoltz;
    noBoltz.weights = unitWeights({ eRMS::EPOT });
    DevOptions cold;
    cold.weights   = noBoltz.weights;
    cold.boltzmann = sameTemperature({ eRMS::EPOT }, 1.0);
    auto rref  = computeDev(monomerFile(), "ACS-g2", noBoltz, nullptr);
    const auto rcold = computeDev(monomerFile(), "ACS-g2", cold, nullptr);
    ASSERT_EQ(rref.size(), rcold.size());
    for(const auto &r : rcold)
    {
        EXPECT_TRUE(std::isfinite(r.second.first)) << r.first;
        EXPECT_TRUE(std::isfinite(r.second.second)) << r.first;
        // The lowest energy structure always has weight 1
        EXPECT_GE(r.second.first, 1) << r.first;
        EXPECT_LE(r.second.first, rref[r.first].first) << r.first;
    }
}

// Test: At very low temperature the dimer weights hit the lower limit
TEST_F(ForceEnergyDevComputerTest, VeryLowTemperatureDimer)
{
    DevOptions noBoltz;
    noBoltz.weights = unitWeights(dimerTargets);
    DevOptions cold;
    cold.weights   = noBoltz.weights;
    cold.boltzmann = sameTemperature({ eRMS::Interaction }, 1.0);
    auto rref  = computeDev(dimerFile(), "ACS-g2", noBoltz, nullptr);
    const auto rcold = computeDev(dimerFile(), "ACS-g2", cold, nullptr);
    ASSERT_EQ(rref.size(), rcold.size());
    for(const auto &r : rcold)
    {
        EXPECT_TRUE(std::isfinite(r.second.first)) << r.first;
        EXPECT_TRUE(std::isfinite(r.second.second)) << r.first;
        EXPECT_GT(r.second.first, 0) << r.first;
        EXPECT_LT(r.second.first, rref[r.first].first) << r.first;
    }
}

// Test: Monomers have no interaction energies, the interaction targets get no data
TEST_F(ForceEnergyDevComputerTest, InteractionTargetsOnMonomers)
{
    DevOptions opts;
    opts.weights = unitWeights(dimerTargets);
    const auto results = computeDev(monomerFile(), "ACS-g2", opts, nullptr);
    EXPECT_FALSE(results.empty());
    for(const auto &r : results)
    {
        EXPECT_EQ(0, r.second.first) << r.first;
        EXPECT_EQ(0, r.second.second) << r.first;
    }
}

/*! \brief Check that a very high Boltzmann temperature gives the same
 * results as no Boltzmann weighting, since all weights are then one.
 * \param[in] rref The results without Boltzmann weighting
 * \param[in] rhot The results at a very high temperature
 */
static void compareVeryHighTemperature(const DevResults &rref,
                                       const DevResults &rhot)
{
    ASSERT_EQ(rref.size(), rhot.size());
    for(const auto &r : rhot)
    {
        const auto ref = rref.find(r.first);
        ASSERT_TRUE(rref.end() != ref) << r.first;
        EXPECT_NEAR(ref->second.first, r.second.first, 1e-6*std::abs(ref->second.first)) << r.first;
        EXPECT_NEAR(ref->second.second, r.second.second, 1e-6*std::abs(ref->second.second)) << r.first;
    }
}

// Test: At very high temperature all monomer conformations have weight one
TEST_F(ForceEnergyDevComputerTest, VeryHighTemperatureMonomer)
{
    DevOptions noBoltz;
    noBoltz.weights = unitWeights({ eRMS::EPOT });
    DevOptions hot;
    hot.weights   = noBoltz.weights;
    hot.boltzmann = sameTemperature({ eRMS::EPOT }, 1e12);
    compareVeryHighTemperature(computeDev(monomerFile(), "ACS-g2", noBoltz, nullptr),
                               computeDev(monomerFile(), "ACS-g2", hot, nullptr));
}

// Test: At very high temperature all dimer conformations have weight one
TEST_F(ForceEnergyDevComputerTest, VeryHighTemperatureDimer)
{
    DevOptions noBoltz;
    noBoltz.weights = unitWeights(dimerTargets);
    DevOptions hot;
    hot.weights   = noBoltz.weights;
    hot.boltzmann = sameTemperature({ eRMS::Interaction }, 1e12);
    compareVeryHighTemperature(computeDev(dimerFile(), "ACS-g2", noBoltz, nullptr),
                               computeDev(dimerFile(), "ACS-g2", hot, nullptr));
}

//! Molprop observables to read when preparing molecules
static const std::map<MolPropObservable, iqmType> propertyIqm = {
    { MolPropObservable::POTENTIAL,      iqmType::QM },
    { MolPropObservable::DIPOLE,         iqmType::QM },
    { MolPropObservable::QUADRUPOLE,     iqmType::QM },
    { MolPropObservable::OCTUPOLE,       iqmType::QM },
    { MolPropObservable::HEXADECAPOLE,   iqmType::QM },
    { MolPropObservable::POLARIZABILITY, iqmType::QM }
};

// Test fixture for the DevComputers of molecular properties
//! \return the energies of all terms and conformations of the molecules
static std::vector<double> allEnergies(std::vector<ACTMol> *mols)
{
    std::vector<double> e;
    for(auto &mol : *mols)
    {
        for(size_t what = 0; what < 3; what++)
        {
            auto &ener = mol.conformationData(what)->energies;
            for(auto t : ener.types())
            {
                const auto term = ener.find(t);
                for(size_t c = 0; c < ener.nconf(); c++)
                {
                    e.push_back(term.has(c) ? term.act[c] : 0.0);
                }
            }
        }
    }
    return e;
}

/*! \brief Compute the energies of all molecules
 * \param[in]    forceComputer To compute the energies
 * \param[in]    pd            The force field
 * \param[inout] mols          The molecules
 * \param[in]    fresh         Whether to forget the energies of earlier calls first
 * \return the energies, see allEnergies
 */
static std::vector<double> computeAll(const ForceComputer *forceComputer,
                                      const ForceField    *pd,
                                      std::vector<ACTMol> *mols,
                                      bool                 fresh)
{
    if (fresh)
    {
        for(auto &mol : *mols)
        {
            for(size_t what = 0; what < 3; what++)
            {
                mol.conformationData(what)->cache.clear();
            }
        }
    }
    MsgHandler msghandler;
    forceComputer->computeRelaxedEnergies(&msghandler, pd, pointers(mols), true, false, false);
    EXPECT_TRUE(msghandler.ok());
    return allEnergies(mols);
}

/*! \brief Fill the parameters of the topologies of the molecules
 * \param[in]    pd   The force field
 * \param[inout] mols The molecules
 */
static void fillAll(const ForceField *pd, std::vector<ACTMol> *mols)
{
    MsgHandler msghandler;
    for(auto &mol : *mols)
    {
        mol.topologyPtr()->fillParameters(&msghandler, pd, missingParameters::Ignore);
        for(auto &ft : mol.fragmentHandler()->topologiesPtr())
        {
            ft.fillParameters(&msghandler, pd, missingParameters::Ignore);
        }
    }
}

/*! \brief Change the first parameter of entries of an interaction type
 * \param[inout] pd          The force field
 * \param[in]    itype       The interaction type
 * \param[in]    factor      To multiply the parameters with
 * \param[in]    singlesOnly Only change the entries of single particle
 *                           types, since pairs are made from them
 */
static void changeParameters(ForceField      *pd,
                             InteractionType  itype,
                             double           factor,
                             bool             singlesOnly)
{
    auto params = pd->findForces(itype)->parametersNoRemoval();
    std::vector<ForceFieldParameter *> changed;
    for(auto &entry : *params)
    {
        if (!singlesOnly || entry.first.atoms().size() <= 1)
        {
            auto &p = entry.second.begin()->second;
            p.forceSetValue(factor*p.value());
            p.setUpdated(true);
            changed.push_back(&p);
        }
    }
    EXPECT_FALSE(changed.empty());
    MsgHandler msghandler;
    generateDependentParameter(&msghandler, pd);
    for(auto p : changed)
    {
        p->setUpdated(false);
    }
}

class PropertyDevComputerTest : public gmx::test::CommandLineTestBase
{
protected:
    PropertyDevComputerTest() = default;

    //! \return the monomer test file
    std::string monomerFile() { return fileManager().getInputFilePath("molprop/molprop.xml"); }

    /*! \brief Read molecules, generate topologies and charges like MolGen does
     * \param[in] msghandler    For messages
     * \param[in] pd            The force field
     * \param[in] forceComputer The force computer
     * \param[in] molecules     The molecules to use, all if empty
     * \param[in] mpFile        The molecule file, the monomer file if empty
     * \param[in] modify        If not empty, changes the data of each molecule
     *                          before its topology is made
     * \return The molecules
     */
    std::vector<ACTMol> prepareMolecules(MsgHandler                           *msghandler,
                                         ForceField                           *pd,
                                         const ForceComputer                  *forceComputer,
                                         const std::set<std::string>          &molecules = {},
                                         const std::string                    &mpFile    = "",
                                         const std::function<void(MolProp *)> &modify    = nullptr)
    {
        std::vector<alexandria::MolProp> mps;
        MolPropRead(msghandler, mpFile.empty() ? monomerFile().c_str() : mpFile.c_str(), &mps);
        std::vector<ACTMol> mols;
        for(auto &mp : mps)
        {
            if (!molecules.empty() && molecules.count(mp.getMolname()) == 0)
            {
                continue;
            }
            if (modify)
            {
                modify(&mp);
            }
            ACTMol actmol;
            actmol.Merge(&mp);
            actmol.GenerateTopology(msghandler, pd, missingParameters::Ignore);
            EXPECT_TRUE(msghandler->ok());
            if (!msghandler->ok())
            {
                continue;
            }
            actmol.getExpProps(msghandler, pd, propertyIqm, 0, 100);
            EXPECT_TRUE(msghandler->ok());
            std::vector<gmx::RVec> forces(actmol.atomsConst().size());
            std::vector<gmx::RVec> coords = actmol.xOriginal();
            actmol.generateCharges(msghandler, pd, forceComputer,
                                   pd->chargeGenerationAlgorithm(), &coords, &forces, true);
            EXPECT_TRUE(msghandler->ok());
            mols.push_back(std::move(actmol));
        }
        EXPECT_FALSE(mols.empty());
        return mols;
    }

    /*! \brief Compare virtual sites on top of atoms to the second charges that replace them
     * The force field is read with and without replacing the virtual sites,
     * see ForceField::convertVsite1, and the charges, the energies of all
     * terms and conformations and the forces of the atoms have to be the
     * same, apart from rounding.
     * \param[in] ffFile      The force field
     * \param[in] mpFile      The molecules
     * \param[in] vsiteType   A virtual site type that is replaced
     * \param[in] atomType    Its atom type
     * \param[in] polarizable Whether the force field has shells
     * \param[in] tolerance   Relative tolerance for charges and energies,
     *                        the one for forces is 100 times larger
     * \param[in] alphas      Polarizabilities of shells to change, in the
     *                        unit of the force field
     */
    void compareSecondChargesToVsites(const std::string                   &ffFile,
                                      const std::string                   &mpFile,
                                      const std::string                   &vsiteType,
                                      const std::string                   &atomType,
                                      bool                                 polarizable,
                                      double                               tolerance,
                                      const std::map<std::string, double> &alphas = {})
    {
        ForceField pdVsite, pdCharge;
        readForceField(fileManager().getInputFilePath(ffFile), &pdVsite, nullptr, false);
        readForceField(fileManager().getInputFilePath(ffFile), &pdCharge);
        for(auto pd : { &pdVsite, &pdCharge })
        {
            for(const auto &a : alphas)
            {
                auto params = pd->findForces(InteractionType::POLARIZATION)->parametersNoRemoval();
                const auto pp     = params->find(Identifier(a.first));
                ASSERT_TRUE(params->end() != pp) << a.first;
                pp->second.find("alpha")->second.forceSetValue(a.second);
            }
        }
        ASSERT_TRUE(pdVsite.hasParticleType(vsiteType));
        ASSERT_FALSE(pdCharge.hasParticleType(vsiteType));
        ASSERT_TRUE(pdCharge.findParticleType(atomType)->hasChargeVS());
        ASSERT_EQ(polarizable, pdCharge.polarizable());
        const ForceComputer forceComputer;
        MsgHandler    msghandler;
        auto molsV  = prepareMolecules(&msghandler, &pdVsite, &forceComputer, {},
                                       fileManager().getInputFilePath(mpFile));
        auto molsC  = prepareMolecules(&msghandler, &pdCharge, &forceComputer, {},
                                       fileManager().getInputFilePath(mpFile));
        ASSERT_EQ(molsV.size(), molsC.size());
        ASSERT_FALSE(molsV.empty());
        if (polarizable)
        {
            // The shells have to be there, and stay excluded from their cores
            size_t nshell = 0;
            for(const auto &mol : molsC)
            {
                for(const auto &atom : mol.atomsConst())
                {
                    nshell += atom.pType() == ActParticle::Shell ? 1 : 0;
                }
            }
            EXPECT_LT(0u, nshell);
        }
        for(auto pmols : { &molsV, &molsC })
        {
            MsgHandler mh;
            forceComputer.computeRelaxedEnergies(&mh, pmols == &molsV ? &pdVsite : &pdCharge,
                                                 pointers(pmols), true, false, true);
            EXPECT_TRUE(mh.ok());
        }
        // Whether a particle is a virtual site that was replaced by a second charge
        auto replaced = [&pdCharge](const ActAtom &atom)
        {
            return atom.pType() == ActParticle::Vsite && !pdCharge.hasParticleType(atom.ffType());
        };
        // Compare forces of the other particles only, the forces on virtual
        // sites have been spread over the atoms.
        const auto compareForces = [&replaced, tolerance](const std::vector<ActAtom> &atomsV,
                                      const ConformationData     &coordsV,
                                      const ConformationData     &coordsC)
        {
            ASSERT_EQ(coordsV.nconf(), coordsC.nconf());
            for(size_t c = 0; c < coordsV.nconf(); c++)
            {
                const auto fV = coordsV.f(c);
                const auto fC = coordsC.f(c);
                // In single precision the errors are relative to the
                // largest force, since the terms of the forces cancel
                double fmax = 0;
                for(const auto &f : fV.toRVec())
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        fmax = std::max(fmax, static_cast<double>(std::abs(f[m])));
                    }
                }
                size_t k = 0;
                for(size_t i = 0; i < atomsV.size(); i++)
                {
                    if (replaced(atomsV[i]))
                    {
                        continue;
                    }
                    ASSERT_LT(k, fC.size());
                    for(int m = 0; m < DIM; m++)
                    {
                        double fmag = GMX_DOUBLE ? 1 + std::abs(fV[i][m]) : fmax;
                        EXPECT_NEAR(fV[i][m], fC[k][m], 100*tolerance*fmag);
                    }
                    k++;
                }
                EXPECT_EQ(k, fC.size());
            }
        };
        for(size_t m = 0; m < molsV.size(); m++)
        {
            auto &molV = molsV[m];
            auto &molC = molsC[m];
            SCOPED_TRACE(molV.getMolname());
            // Charges, the second charge is the charge of the virtual site
            const auto &atomsV = molV.atomsConst();
            const auto &atomsC = molC.atomsConst();
            size_t k = 0;
            for(size_t i = 0; i < atomsV.size(); i++)
            {
                if (replaced(atomsV[i]))
                {
                    continue;
                }
                ASSERT_LT(k, atomsC.size());
                double qvs = 0;
                for(auto v : atomsV[i].vsites())
                {
                    if (replaced(atomsV[v]))
                    {
                        qvs += atomsV[v].charge();
                    }
                }
                EXPECT_NEAR(atomsV[i].charge(), atomsC[k].charge(), tolerance);
                EXPECT_NEAR(qvs, atomsC[k].chargeVS(), tolerance);
                k++;
            }
            EXPECT_EQ(k, atomsC.size());
            EXPECT_LT(atomsC.size(), atomsV.size());
            // Energies of the dimer and the monomers
            for(size_t what = 0; what < 3; what++)
            {
                auto &eV = molV.conformationData(what)->energies;
                auto &eC = molC.conformationData(what)->energies;
                ASSERT_EQ(eV.types(), eC.types());
                ASSERT_EQ(eV.nconf(), eC.nconf());
                for(auto t : eV.types())
                {
                    const auto tV = eV.find(t);
                    const auto tC = eC.find(t);
                    for(size_t c = 0; c < eV.nconf(); c++)
                    {
                        ASSERT_EQ(tV.has(c), tC.has(c));
                        if (tV.has(c))
                        {
                            EXPECT_NEAR(tV.act[c], tC.act[c], tolerance*(1 + std::abs(tV.act[c])))
                                << interactionTypeToString(t) << " conformation " << c << " part " << what;
                        }
                    }
                }
            }
            // Forces
            compareForces(atomsV, *molV.conformationData(2), *molC.conformationData(2));
            auto &fragV = molV.fragmentHandler()->topologies();
            for(size_t f = 0; f < fragV.size(); f++)
            {
                compareForces(fragV[f].atoms(), *molV.conformationData(f), *molC.conformationData(f));
            }
            // The coordinates, with the virtual sites on their atoms
            const auto         &xC = molC.xOriginal();
            std::vector<gmx::RVec> xV(atomsV.size(), { 0, 0, 0 });
            k = 0;
            for(size_t i = 0; i < atomsV.size(); i++)
            {
                if (!replaced(atomsV[i]))
                {
                    xV[i] = xC[k++];
                    for(auto v : atomsV[i].vsites())
                    {
                        if (replaced(atomsV[v]))
                        {
                            xV[v] = xC[k-1];
                        }
                    }
                }
            }
            compareMultipoles(atomsV, xV, atomsC, xC, tolerance);
            compareEsp(&pdVsite, molV, xV, &pdCharge, molC, xC, molV.totalCharge(), tolerance);
        }
        compareExpanded(mpFile, pdCharge, &molsV, tolerance);
    }

    /*! \brief Check that second charges turned back into virtual sites,
     * as for GROMACS and OpenMM, give the molecules with virtual sites.
     * \param[in] mpFile    The molecules
     * \param[in] pdCharge  The force field with second charges
     * \param[in] molsV     The molecules with virtual sites, with energies
     * \param[in] tolerance The tolerance
     */
    void compareExpanded(const std::string         &mpFile,
                         const ForceField          &pdCharge,
                         std::vector<ACTMol>       *molsV,
                         double                     tolerance)
    {
        ForceField pdExpand = pdCharge;
        MsgHandler msghandler;
        pdExpand.expandSecondCharges(&msghandler);
        EXPECT_TRUE(msghandler.ok());
        for(const auto &pt : pdExpand.particleTypesConst())
        {
            EXPECT_FALSE(pt.second.hasChargeVS()) << pt.first.id();
        }
        ASSERT_TRUE(pdExpand.interactionPresent(InteractionType::VSITE1));
        const ForceComputer forceComputer;
        auto molsE = prepareMolecules(&msghandler, &pdExpand, &forceComputer, {},
                                      fileManager().getInputFilePath(mpFile));
        ASSERT_EQ(molsV->size(), molsE.size());
        forceComputer.computeRelaxedEnergies(&msghandler, &pdExpand, pointers(&molsE), true, false, true);
        EXPECT_TRUE(msghandler.ok());
        for(size_t m = 0; m < molsV->size(); m++)
        {
            auto &molV = (*molsV)[m];
            SCOPED_TRACE(molV.getMolname());
            const auto &atomsV = molV.atomsConst();
            const auto &atomsE = molsE[m].atomsConst();
            ASSERT_EQ(atomsV.size(), atomsE.size());
            for(size_t i = 0; i < atomsV.size(); i++)
            {
                EXPECT_EQ(atomsV[i].pType(), atomsE[i].pType()) << i;
                EXPECT_NEAR(atomsV[i].charge(), atomsE[i].charge(), tolerance) << i;
                EXPECT_EQ(0, atomsE[i].chargeVS()) << i;
            }
            for(size_t what = 0; what < 3; what++)
            {
                auto &eV = molV.conformationData(what)->energies;
                auto &eE = molsE[m].conformationData(what)->energies;
                ASSERT_EQ(eV.types(), eE.types());
                for(auto t : eV.types())
                {
                    const auto tV = eV.find(t);
                    const auto tE = eE.find(t);
                    for(size_t c = 0; c < eV.nconf(); c++)
                    {
                        ASSERT_EQ(tV.has(c), tE.has(c));
                        if (tV.has(c))
                        {
                            EXPECT_NEAR(tV.act[c], tE.act[c], tolerance*(1 + std::abs(tV.act[c])))
                                << interactionTypeToString(t) << " conformation " << c << " part " << what;
                        }
                    }
                }
            }
        }
    }

    /*! \brief Compare the multipoles of two representations of a molecule
     * \param[in] atomsA    The atoms of the first
     * \param[in] xA        Their coordinates
     * \param[in] atomsB    The atoms of the second
     * \param[in] xB        Their coordinates
     * \param[in] tolerance The tolerance
     */
    void compareMultipoles(const std::vector<ActAtom>   &atomsA,
                           const std::vector<gmx::RVec> &xA,
                           const std::vector<ActAtom>   &atomsB,
                           const std::vector<gmx::RVec> &xB,
                           double                        tolerance)
    {
        QtypeProps qA(qPropertyType::ACM, atomsA, xA);
        QtypeProps qB(qPropertyType::ACM, atomsB, xB);
        MsgHandler msghandler;
        for(auto q : { &qA, &qB })
        {
            q->initializeMoments();
            q->calcMoments(&msghandler);
        }
        for(auto mpo : { MolPropObservable::DIPOLE, MolPropObservable::QUADRUPOLE })
        {
            auto mA = qA.getMultipole(&msghandler, mpo);
            auto mB = qB.getMultipole(&msghandler, mpo);
            ASSERT_EQ(mA.size(), mB.size());
            double norm = 0;
            for(size_t m = 0; m < mA.size(); m++)
            {
                EXPECT_NEAR(mA[m], mB[m], tolerance*(1 + std::abs(mA[m])))
                    << mpo_name(mpo) << " " << m;
                norm += std::abs(mA[m]);
            }
            if (MolPropObservable::DIPOLE == mpo)
            {
                EXPECT_LT(0, norm);
            }
        }
    }

    /*! \brief Compare the electrostatic potential and the RESP charges
     * of two representations of a molecule
     * The potential of the first is used as the reference for fitting.
     * \param[in] pdA       The force field of the first
     * \param[in] molA      The first molecule
     * \param[in] xA        Its coordinates
     * \param[in] pdB       The force field of the second
     * \param[in] molB      The second molecule
     * \param[in] xB        Its coordinates
     * \param[in] qtot      The total charge
     * \param[in] tolerance The tolerance
     */
    void compareEsp(const ForceField             *pdA,
                    const ACTMol                 &molA,
                    const std::vector<gmx::RVec> &xA,
                    const ForceField             *pdB,
                    const ACTMol                 &molB,
                    const std::vector<gmx::RVec> &xB,
                    int                           qtot,
                    double                        tolerance)
    {
        MsgHandler msghandler;
        auto pot = pdA->potential(InteractionType::ELECTROSTATICS);
        std::vector<QgenResp> resp(2);
        std::vector<const ForceField *>             pds   = { pdA, pdB };
        std::vector<const ACTMol *>                 mols  = { &molA, &molB };
        std::vector<const std::vector<gmx::RVec> *> xs    = { &xA, &xB };
        for(size_t r = 0; r < resp.size(); r++)
        {
            resp[r].setChargeDistributionType(potentialToChargeDistributionType(pot));
            resp[r].setAtomInfo(&msghandler, mols[r]->atomsConst(), pds[r], qtot);
            resp[r].updateAtomCharges(&msghandler, mols[r]->atomsConst());
            resp[r].updateAtomCoords(*xs[r]);
            resp[r].setAtomSymmetry({});
        }
        resp[0].makeGrid(0.1, 0.3, xA);
        ASSERT_LT(0u, resp[0].nEsp());
        resp[0].calcPot(&msghandler, 1);
        // The potential of the first is the reference for the second
        for(auto &ep : resp[0].espPoints())
        {
            auto x = ep.esp();
            resp[1].addEspPoint(x[XX], x[YY], x[ZZ], ep.vCalc());
        }
        resp[1].calcPot(&msghandler, 1);
        ASSERT_EQ(resp[0].nEsp(), resp[1].nEsp());
        for(size_t j = 0; j < resp[0].nEsp(); j++)
        {
            double v = resp[0].espPoint(j).vCalc();
            EXPECT_NEAR(v, resp[1].espPoint(j).vCalc(), tolerance*(1 + std::abs(v))) << j;
        }
        // Fitting the charges of the second to the potential of the first
        // gives the charges back, the second charges are fixed.
        resp[1].optimizeCharges(&msghandler, 1);
        EXPECT_TRUE(msghandler.ok());
        const auto &atomsB = molB.atomsConst();
        for(size_t i = 0; i < atomsB.size(); i++)
        {
            EXPECT_NEAR(atomsB[i].charge(), resp[1].getCharge(i), 1e-6) << i;
        }
    }

    /*! \brief Check that energies used again are those computed from scratch
     * The energies of an interaction type are used again if its parameters,
     * the coordinates and, if used, the charges did not change, see
     * CachedEnergies. After changing each of them, the energies have to
     * change and to be the same as without the energies of earlier calls.
     * \param[in] ffName    The force field
     * \param[in] mpFile    The molecules, the monomer file if empty
     * \param[in] molecules The molecules to use, all if empty
     * \param[in] vsite     Whether to change a virtual site parameter too
     */
    void testEnergiesUsedAgain(const std::string             &ffName,
                               const std::string             &mpFile,
                               const std::set<std::string>   &molecules,
                               bool                           vsite)
    {
        const ForceComputer forceComputer;
        MsgHandler    msghandler;
        // A copy, since the parameters are changed
        ForceField    pd   = *getForceField(ffName);
        ASSERT_FALSE(pd.polarizable());
        auto          mols = prepareMolecules(&msghandler, &pd, &forceComputer, molecules, mpFile);
        ASSERT_FALSE(mols.empty());
        auto &timings = Timings::instance();
        auto  reused  = timings.index("Energies used again");
        const auto e0      = computeAll(&forceComputer, &pd, &mols, false);
        // Nothing changed, all energies are used again
        double ncalls = timings.calls()[reused];
        EXPECT_EQ(e0, computeAll(&forceComputer, &pd, &mols, false));
        EXPECT_LT(ncalls, timings.calls()[reused]);
        // Van der Waals parameters
        changeParameters(&pd, InteractionType::VDW, 1.1, true);
        fillAll(&pd, &mols);
        const auto e1 = computeAll(&forceComputer, &pd, &mols, false);
        EXPECT_NE(e0, e1);
        EXPECT_EQ(e1, computeAll(&forceComputer, &pd, &mols, true));
        // Charges, virtual sites get a charge such that their position matters
        const auto newCharge = [](ActAtom *atom)
        {
            atom->setCharge(0.9*atom->charge() + (atom->pType() == ActParticle::Vsite ? 0.1 : 0.0));
        };
        for(auto &mol : mols)
        {
            for(auto &atom : *mol.atoms())
            {
                newCharge(&atom);
            }
            for(auto &ft : mol.fragmentHandler()->topologiesPtr())
            {
                for(auto &atom : *ft.atomsPtr())
                {
                    newCharge(&atom);
                }
            }
        }
        const auto e2 = computeAll(&forceComputer, &pd, &mols, false);
        EXPECT_NE(e1, e2);
        EXPECT_EQ(e2, computeAll(&forceComputer, &pd, &mols, true));
        if (vsite)
        {
            // Virtual sites move, and with them the coordinates
            changeParameters(&pd, InteractionType::VSITE3FD, 1.2, false);
            fillAll(&pd, &mols);
            const auto e3 = computeAll(&forceComputer, &pd, &mols, false);
            EXPECT_NE(e2, e3);
            EXPECT_EQ(e3, computeAll(&forceComputer, &pd, &mols, true));
        }
    }

    //! Set the bounds of all ACM charges in a force field
    static void setAcmChargeBounds(ForceField *pd, double qmin, double qmax)
    {
        for(auto &pt : *pd->particleTypes())
        {
            if (pt.second.hasParameter("charge"))
            {
                auto p = pt.second.parameter("charge");
                if (p->mutability() == Mutability::ACM)
                {
                    p->setMaximum(qmax);
                    EXPECT_TRUE(p->setMinimum(qmin));
                }
            }
        }
    }

    /*! \brief Compute the charge bounds deviation for each molecule
     * \param[in] ffFile Force field name
     * \param[in] change Whether to change the bounds of the ACM charges
     * \param[in] qmin   New minimum charge
     * \param[in] qmax   New maximum charge
     * \return The chi squared for each molecule
     */
    std::vector<double> testBounds(const std::string &ffFile,
                                   bool               change,
                                   double             qmin,
                                   double             qmax)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-11);
        checker_.setDefaultTolerance(tolerance);
        const ForceComputer forceComputer;
        MsgHandler    msghandler;
        auto          pd   = getForceField(ffFile);
        auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
        // Change the bounds in a copy only, the original is shared between tests
        ForceField    pdBounds(*pd);
        if (change)
        {
            setAcmChargeBounds(&pdBounds, qmin, qmax);
        }
        BoundsDevComputer   computer(nullptr);
        EXPECT_EQ("Bounds", computer.name());
        const auto newTargets = []()
        {
            std::map<eRMS, FittingTarget> targets;
            targets.insert({eRMS::BOUNDS, FittingTarget(eRMS::BOUNDS, iMolSelect::Train)});
            targets.find(eRMS::BOUNDS)->second.setWeight(1);
            return targets;
        };
        // One molecule at a time
        std::vector<double> chi2;
        for(auto &mol : mols)
        {
            const std::vector<ACTMol *> one = { &mol };
            auto targets = newTargets();
            computer.calcDeviation(&msghandler, &forceComputer, one,
                                   &targets, &pdBounds, LossFunction::MSE);
            auto &t = targets.find(eRMS::BOUNDS)->second;
            EXPECT_EQ(1, t.totalWeight());
            checker_.checkDouble(t.chiSquared(),
                                 gmx::formatString("%s-Bounds-chi2", mol.getMolname().c_str()).c_str());
            chi2.push_back(t.chiSquared());
        }
        // All molecules at once give the sum
        auto targets = newTargets();
        computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols),
                               &targets, &pdBounds, LossFunction::MSE);
        auto &t = targets.find(eRMS::BOUNDS)->second;
        EXPECT_EQ(mols.size(), t.totalWeight());
        double sum = std::accumulate(chi2.begin(), chi2.end(), 0.0);
        EXPECT_NEAR(sum, t.chiSquared(), 1e-10*(std::abs(sum) + 1));
        return chi2;
    }

    //! Compute the multipole deviation for each molecule
    void testMultipole(MolPropObservable mpo, eRMS erms, const std::string &ffFile)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-8);
        checker_.setDefaultTolerance(tolerance);
        const ForceComputer forceComputer;
        MsgHandler    msghandler;
        auto          pd   = getForceField(ffFile);
        auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
        MultiPoleDevComputer computer(mpo);
        EXPECT_EQ(mpo_name(mpo), computer.name());
        const auto newTargets = [erms]()
        {
            std::map<eRMS, FittingTarget> targets;
            targets.insert({erms, FittingTarget(erms, iMolSelect::Train)});
            targets.find(erms)->second.setWeight(1);
            return targets;
        };
        // Compute the moments first, as ACMFitnessComputer::computeMultipoles does
        for(auto &mol : mols)
        {
            auto qProps = mol.qProps();
            for(auto qp = qProps->begin(); qp < qProps->end(); ++qp)
            {
                auto qcalc = qp->qPact();
                qcalc->setQ(mol.atomsConst());
                qcalc->calcMoments();
            }
        }
        // One molecule at a time
        double weightSum = 0, chi2Sum = 0;
        for(auto &mol : mols)
        {
            const std::vector<ACTMol *> one = { &mol };
            auto targets = newTargets();
            computer.calcDeviation(&msghandler, &forceComputer, one,
                                   &targets, pd, LossFunction::MSE);
            auto &t = targets.find(erms)->second;
            EXPECT_GT(t.totalWeight(), 0) << mol.getMolname();
            EXPECT_TRUE(std::isfinite(t.chiSquared())) << mol.getMolname();
            checker_.checkDouble(t.totalWeight(),
                                 gmx::formatString("%s-%s-weight", mol.getMolname().c_str(), rmsName(erms)).c_str());
            checker_.checkDouble(t.chiSquared(),
                                 gmx::formatString("%s-%s-chi2", mol.getMolname().c_str(), rmsName(erms)).c_str());
            weightSum += t.totalWeight();
            chi2Sum   += t.chiSquared();
        }
        // All molecules at once give the sum
        auto targets = newTargets();
        computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols),
                               &targets, pd, LossFunction::MSE);
        auto &t = targets.find(erms)->second;
        EXPECT_DOUBLE_EQ(weightSum, t.totalWeight());
        EXPECT_NEAR(chi2Sum, t.chiSquared(), 1e-10*(chi2Sum + 1));
    }

    /*! \brief Check the force deviation of dimers
     * The reference forces differ by one from the computed ones in each
     * component of each atom, the forces of shells and virtual sites are
     * not used.
     * \param[in] ffFile The force field
     */
    void testDimerForces(const std::string &ffFile)
    {
        const ForceComputer forceComputer;
        MsgHandler    msghandler;
        auto          pd   = getForceField(ffFile);
        auto          mols = prepareMolecules(&msghandler, pd, &forceComputer, {},
                                              fileManager().getInputFilePath("molprop/mpdimer.xml"));
        ASSERT_FALSE(mols.empty());
        forceComputer.computeRelaxedEnergies(&msghandler, pd, pointers(&mols), true, true, true);
        double ncomponent = 0;
        size_t nother     = 0;
        for(auto &mol : mols)
        {
            ASSERT_TRUE(mol.isDimer());
            const auto &atoms = mol.topology()->atoms();
            auto *cdata = mol.conformationData(2);
            for(size_t c = 0; c < cdata->nconf(); c++)
            {
                auto &dim = cdata->confs[c];
                dim.fobs.clear();
                for(size_t i = 0; i < atoms.size(); i++)
                {
                    if (atoms[i].pType() == ActParticle::Atom)
                    {
                        gmx::RVec f = cdata->f(c)[i].rvec();
                        f[XX] += 1;
                        f[YY] -= 1;
                        f[ZZ] += 1;
                        dim.fobs.push_back(f);
                        ncomponent += DIM;
                    }
                    else
                    {
                        nother += 1;
                    }
                }
            }
        }
        if (pd->polarizable())
        {
            EXPECT_LT(0u, nother);
        }
        ForceEnergyDevComputer computer({});
        std::map<eRMS, FittingTarget> targets;
        targets.insert({eRMS::Force2, FittingTarget(eRMS::Force2, iMolSelect::Train)});
        targets.find(eRMS::Force2)->second.setWeight(1);
        computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols),
                               &targets, pd, LossFunction::MSE);
        EXPECT_TRUE(msghandler.ok());
        auto &t = targets.find(eRMS::Force2)->second;
        EXPECT_DOUBLE_EQ(ncomponent, t.totalWeight());
        EXPECT_NEAR(ncomponent, t.chiSquared(), 1e-6*ncomponent);
    }

    /*! \brief Check the ESP deviation
     * It is the sum of the squared differences between the QM and the
     * computed potential over all ESP points, its weight their number.
     * \param[in] ffFile The force field
     */
    void testEsp(const std::string &ffFile)
    {
        const ForceComputer forceComputer;
        MsgHandler    msghandler;
        auto          pd   = getForceField(ffFile);
        auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
        EspDevComputer computer;
        double weightSum = 0, chi2Sum = 0;
        for(auto &mol : mols)
        {
            const std::vector<ACTMol *> one = { &mol };
            std::map<eRMS, FittingTarget> targets;
            targets.insert({eRMS::ESP, FittingTarget(eRMS::ESP, iMolSelect::Train)});
            targets.find(eRMS::ESP)->second.setWeight(1);
            computer.calcDeviation(&msghandler, &forceComputer, one, &targets, pd, LossFunction::MSE);
            auto  &t    = targets.find(eRMS::ESP)->second;
            size_t nesp = 0;
            double sum2 = 0;
            auto   qProps = mol.qProps();
            for(auto qp = qProps->begin(); qp < qProps->end(); ++qp)
            {
                for(const auto &ep : qp->qgenResp()->espPoints())
                {
                    sum2 += gmx::square(ep.v() - ep.vCalc());
                    nesp += 1;
                }
            }
            EXPECT_LT(0u, nesp) << mol.getMolname();
            EXPECT_LT(0, sum2) << mol.getMolname();
            EXPECT_DOUBLE_EQ(nesp, t.totalWeight()) << mol.getMolname();
            EXPECT_NEAR(sum2, t.chiSquared(), 1e-8*sum2) << mol.getMolname();
            weightSum += t.totalWeight();
            chi2Sum   += t.chiSquared();
        }
        // All molecules at once give the sum
        std::map<eRMS, FittingTarget> targets;
        targets.insert({eRMS::ESP, FittingTarget(eRMS::ESP, iMolSelect::Train)});
        targets.find(eRMS::ESP)->second.setWeight(1);
        computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols), &targets, pd, LossFunction::MSE);
        auto &t = targets.find(eRMS::ESP)->second;
        EXPECT_DOUBLE_EQ(weightSum, t.totalWeight());
        EXPECT_NEAR(chi2Sum, t.chiSquared(), 1e-8*chi2Sum);
    }

    /*! \brief Check the deviation of frequencies or intensities
     * The reference values are those computed for water plus a shift,
     * which then gives the deviation. If the number of reference values is
     * wrong the molecule is not used.
     * \param[in] mpo  The observable
     * \param[in] erms The fitting target
     */
    void testHarmonics(MolPropObservable mpo, eRMS erms)
    {
        const ForceComputer forceComputer;
        MsgHandler    msghandler;
        auto          pd   = getForceField("ACS-g2");
        auto          mols = prepareMolecules(&msghandler, pd, &forceComputer, { "water" });
        ASSERT_EQ(1u, mols.size());
        ASSERT_EQ(JobType::OPT, mols[0].jobType());
        // The frequencies and intensities as HarmonicsDevComputer computes them
        const MolHandler        handler;
        const SimulationConfigHandler simConfig;
        auto coords = mols[0].xOriginal();
        ASSERT_EQ(eMinimizeStatus::OK,
                  handler.minimizeCoordinates(&msghandler, pd, mols[0].topology(), &forceComputer,
                                              simConfig, &coords, nullptr, {}));
        std::vector<double> frequencies, intensities;
        handler.nma(pd, &mols[0], &forceComputer, &coords, &frequencies, &intensities);
        const auto &values = MolPropObservable::FREQUENCY == mpo ? frequencies : intensities;
        ASSERT_FALSE(values.empty());
        const char *unit  = mpo_unit2(mpo);
        double      shift = 2;
        for(size_t nref : { values.size(), values.size() - 1 })
        {
            const auto addHarmonics = [&](MolProp *mp)
            {
                auto harm = std::make_unique<Harmonics>(unit, 0, mpo);
                for(size_t k = 0; k < nref; k++)
                {
                    harm->addValue(convertFromGromacs(values[k], unit) + shift);
                }
                (*mp->experiment())[0].addProperty(mpo, std::move(harm));
            };
            auto withRef = prepareMolecules(&msghandler, pd, &forceComputer, { "water" }, "", addHarmonics);
            ASSERT_EQ(1u, withRef.size());
            HarmonicsDevComputer computer(mpo);
            std::map<eRMS, FittingTarget> targets;
            targets.insert({erms, FittingTarget(erms, iMolSelect::Train)});
            targets.find(erms)->second.setWeight(1);
            computer.calcDeviation(&msghandler, &forceComputer, pointers(&withRef),
                                   &targets, pd, LossFunction::MSE);
            auto &t = targets.find(erms)->second;
            if (nref == values.size())
            {
                double expected = values.size()*gmx::square(convertToGromacs(shift, unit));
                EXPECT_EQ(1, t.totalWeight());
                EXPECT_NEAR(expected, t.chiSquared(), 1e-6*expected);
            }
            else
            {
                EXPECT_EQ(0, t.totalWeight());
            }
        }
    }
};

/*! \brief Tolerance for comparing second charges with virtual sites
 * In single precision the energy of each interaction is computed in single
 * precision, the virtual sites and the second charges round them
 * differently. That changes the Coulomb energies by about 1e-6 relative,
 * and the interaction energies of dimers, differences of sums of hundreds
 * of terms, by up to 1e-3 kJ/mol.
 * \return the tolerance relative to one plus the energy
 */
static double secondChargesTolerance()
{
    return gmx::test::precisionTolerance(1e-10, 1e-3);
}


// Test: Virtual sites on top of atoms give the same charges, energies and
// forces as the second charges that replace them, see ForceField::convertVsite1
TEST_F(PropertyDevComputerTest, SecondChargesAsVsites)
{
    compareSecondChargesToVsites("forcefields/PC+GV-elec.xml", "molprop/mpdimer.xml",
                                 "v1ow", "ow", false, secondChargesTolerance());
}

// Test: The same for a polarizable force field, with shells on the ions,
// which are excluded from their own second charge as from the virtual site.
TEST_F(PropertyDevComputerTest, SecondChargesAsVsitesPolarizable)
{
    compareSecondChargesToVsites("forcefields/PG-pol-ions.xml", "molprop/ion-water.xml",
                                 "v1Na+", "Na+", true, secondChargesTolerance());
}

// Test: The same with polarizable ions, the ion shells have a negligible
// polarizability in the force field file, such that they hardly move.
TEST_F(PropertyDevComputerTest, SecondChargesAsVsitesPolarizableIons)
{
    compareSecondChargesToVsites("forcefields/PG-pol-ions.xml", "molprop/ion-water.xml",
                                 "v1Na+", "Na+", true, secondChargesTolerance(),
                                 { { "Na+_s", 0.18 }, { "Cl-_s", 3.4 } });
}

// Test: After training polarizabilities the shells are relaxed with the new
// force constants, giving the same energies as molecules made with them.
TEST_F(PropertyDevComputerTest, ShellsRelaxedAfterChangingAlpha)
{
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    ForceField    pd;
    readForceField(fileManager().getInputFilePath("forcefields/PG-pol-ions.xml"), &pd);
    ASSERT_TRUE(pd.polarizable());
    const auto mpFile  = fileManager().getInputFilePath("molprop/ion-water.xml");
    auto mols    = prepareMolecules(&msghandler, &pd, &forceComputer, {}, mpFile);
    ASSERT_FALSE(mols.empty());
    const auto before  = computeAll(&forceComputer, &pd, &mols, false);
    // Polarizable ions, as in training
    for(const auto &a : std::map<std::string, double>{ { "Na+_s", 0.18 }, { "Cl-_s", 3.4 } })
    {
        auto params = pd.findForces(InteractionType::POLARIZATION)->parametersNoRemoval();
        const auto pp     = params->find(Identifier(a.first));
        ASSERT_TRUE(params->end() != pp) << a.first;
        pp->second.find("alpha")->second.forceSetValue(a.second);
    }
    fillAll(&pd, &mols);
    auto after   = computeAll(&forceComputer, &pd, &mols, false);
    EXPECT_NE(before, after);
    // Molecules made with the new polarizabilities
    auto newMols = prepareMolecules(&msghandler, &pd, &forceComputer, {}, mpFile);
    auto fresh   = computeAll(&forceComputer, &pd, &newMols, false);
    ASSERT_EQ(fresh.size(), after.size());
    for(size_t i = 0; i < fresh.size(); i++)
    {
        EXPECT_NEAR(fresh[i], after[i], 1e-10*(1 + std::abs(fresh[i]))) << i;
    }
}

//! Computes on the processor while it exists, and then as before
class ProcessorOnly
{
public:
    ProcessorOnly() : before_(gpu::setUseGpu(false)) {}
    ~ProcessorOnly() { (void) gpu::setUseGpu(before_); }
private:
    //! Whether the GPU was used before
    bool before_;
};

/*! \brief Check which shells were kept, see ForceComputer::setKeepShells
 * \param[in] mols The molecules
 * \param[in] keep Whether the shells of the conformations that converged
 *                 were to be kept, otherwise none
 * \return the number of conformations with kept shells
 */
static size_t checkKeptShells(std::vector<ACTMol> *mols, bool keep)
{
    size_t nkept = 0;
    for(auto &mol : *mols)
    {
        auto &tops = mol.fragmentHandler()->topologies();
        for(size_t what = 0; what < 3; what++)
        {
            if ((what < 2 && what >= tops.size()) || (what == 2 && !mol.isDimer()))
            {
                continue;
            }
            const Topology *top    = (what < 2) ? &tops[what] : mol.topology();
            size_t          nshell = 0;
            for(const auto &a : top->atoms())
            {
                nshell += (a.pType() == ActParticle::Shell) ? 1 : 0;
            }
            auto *cdata = mol.conformationData(what);
            for(size_t c = 0; c < cdata->nconf(); c++)
            {
                SCOPED_TRACE(gmx::formatString("%s what %zu conformation %zu", mol.getMolname().c_str(), what, c));
                const auto &conf = cdata->confs[c];
                // Shells that did not converge start from their cores next time
                EXPECT_EQ(keep && conf.done ? DIM*nshell : 0, conf.xkept.size());
                nkept += conf.xkept.empty() ? 0 : 1;
            }
        }
    }
    return nkept;
}

/*! \brief Compare starting the relaxation of the shells from where they were
 * relaxed the previous time (train_ff -shell_keep) to starting from the cores
 * The energies have to be the same within the tolerance on the forces, and
 * fewer iterations are needed after a small change of the parameters, as in
 * an MCMC step. Each force computer has its own molecules, since the
 * conformations remember whether their shells converged.
 * \param[in] pd        The force field, its polarizabilities are changed
 * \param[in] makeMols  Makes the molecules
 * \param[in] resident  Whether to compute on the GPU, see ForceComputer::setGpuResident
 * \param[in] tolerance Relative tolerance for the energies after the change
 * \return false if the GPU cannot compute the molecules, true otherwise
 */
static bool compareShellKeep(ForceField                                                   *pd,
                             const std::function<std::vector<ACTMol>(const ForceComputer *)> &makeMols,
                             bool                                                          resident,
                             double                                                        tolerance)
{
    // A tight tolerance, since the kept shells converge to slightly other positions
    ForceComputer fromCores(1e-8, 100, pd->epsilonr());
    ForceComputer kept(1e-8, 100, pd->epsilonr());
    kept.setKeepShells(true);
    EXPECT_TRUE(kept.keepShells());
    EXPECT_FALSE(fromCores.keepShells());
    for(auto fc : { &fromCores, &kept })
    {
        fc->setGpuResident(resident);
    }
    auto molsCores = makeMols(&fromCores);
    auto molsKept  = makeMols(&kept);
    EXPECT_FALSE(molsCores.empty());
    if (resident && (!fromCores.prepareGpuResident(pd, pointers(&molsCores)) ||
                     !kept.prepareGpuResident(pd, pointers(&molsKept))))
    {
        return false;
    }
    // The first time both start from the cores
    const auto first     = computeAll(&fromCores, pd, &molsCores, false);
    const auto firstKept = computeAll(&kept, pd, &molsKept, false);
    EXPECT_EQ(first, firstKept);
    EXPECT_EQ(fromCores.numShellIterations(), kept.numShellIterations());
    if (!resident)
    {
        // The GPU keeps the shells on the GPU
        EXPECT_EQ(0u, checkKeptShells(&molsCores, false));
        EXPECT_GT(checkKeptShells(&molsKept, true), 0u);
    }
    // Slightly larger polarizabilities
    changeParameters(pd, InteractionType::POLARIZATION, 1.02, true);
    fillAll(pd, &molsCores);
    fillAll(pd, &molsKept);
    kept.resetStatistics();
    fromCores.resetStatistics();
    auto reference = computeAll(&fromCores, pd, &molsCores, false);
    auto after     = computeAll(&kept, pd, &molsKept, false);
    EXPECT_NE(first, reference);
    EXPECT_EQ(reference.size(), after.size());
    for(size_t i = 0; i < reference.size() && i < after.size(); i++)
    {
        EXPECT_NEAR(reference[i], after[i], tolerance*(1 + std::abs(reference[i]))) << i;
    }
    EXPECT_EQ(0u, fromCores.numShellConvergenceFailed());
    EXPECT_EQ(0u, kept.numShellConvergenceFailed());
    EXPECT_GT(fromCores.numShellIterations(), 0u);
    EXPECT_LT(kept.numShellIterations(), fromCores.numShellIterations());
    return true;
}

// Test: The relaxation stops when the shells no longer move, since their
// displacements are below the resolution of the coordinates. With the very
// stiff spring of the shell of chloride (about 1e13 kJ/mol nm^2) this
// happens at a force of about 1e-4 kJ/mol nm, above a tolerance of 1e-8 on
// the mean square force. In double precision, the relaxation of one
// conformation then went on until the maximum number of iterations. On the
// processor.
TEST_F(PropertyDevComputerTest, StiffShellsStopAtCoordinateResolution)
{
    const ProcessorOnly cpu;
    MsgHandler    msghandler;
    ForceField    pd;
    readForceField(fileManager().getInputFilePath("forcefields/PG-pol-ions.xml"), &pd);
    const int     maxiter = 100;
    const ForceComputer forceComputer(1e-8, maxiter, pd.epsilonr());
    auto mols = prepareMolecules(&msghandler, &pd, &forceComputer, { "water#chloride" },
                                 fileManager().getInputFilePath("molprop/ion-water.xml"));
    ASSERT_EQ(1u, mols.size());
    (void) computeAll(&forceComputer, &pd, &mols, false);
    EXPECT_EQ(0u, forceComputer.numShellConvergenceFailed());
    // Slightly larger polarizabilities, after which one conformation stalled
    changeParameters(&pd, InteractionType::POLARIZATION, 1.02, true);
    fillAll(&pd, &mols);
    forceComputer.resetStatistics();
    (void) computeAll(&forceComputer, &pd, &mols, false);
    EXPECT_EQ(0u, forceComputer.numShellConvergenceFailed());
    EXPECT_LT(forceComputer.numShellIterations(), static_cast<size_t>(maxiter));
    for(size_t what = 0; what < 3; what++)
    {
        for(const auto &conf : mols[0].conformationData(what)->confs)
        {
            EXPECT_TRUE(conf.done);
        }
    }
}

// Test: -shell_keep on the processor gives the energies of starting from the
// cores, with fewer iterations, and keeps the shells that converged
TEST_F(PropertyDevComputerTest, ShellKeepSameEnergies)
{
    const ProcessorOnly cpu;
    ForceField    pd;
    readForceField(fileManager().getInputFilePath("forcefields/PG-pol-ions.xml"), &pd);
    ASSERT_TRUE(pd.polarizable());
    auto mpFile = fileManager().getInputFilePath("molprop/ion-water.xml");
    const auto make   = [this, &pd, &mpFile](const ForceComputer *fc)
    {
        MsgHandler msghandler;
        return prepareMolecules(&msghandler, &pd, fc, {}, mpFile);
    };
    EXPECT_TRUE(compareShellKeep(&pd, make, false, 1e-6));
}

// Test: -shell_keep on the GPU gives the energies of starting from the cores,
// within single precision, with fewer iterations
TEST_F(PropertyDevComputerTest, ShellKeepSameEnergiesGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && !GMX_DOUBLE, "no GPU or double precision"))
    {
        return;
    }
    ForceField pd;
    readForceField(fileManager().getInputFilePath("forcefields/PG-pol-ions.xml"), &pd);
    ASSERT_TRUE(pd.polarizable());
    auto mpFile = fileManager().getInputFilePath("molprop/ion-water.xml");
    const auto make   = [this, &pd, &mpFile](const ForceComputer *fc)
    {
        MsgHandler msghandler;
        return prepareMolecules(&msghandler, &pd, fc, {}, mpFile);
    };
    EXPECT_TRUE(compareShellKeep(&pd, make, true, 1e-5)) << "The GPU cannot compute these molecules";
}

// Test: With -gpu_resident, molecules that the GPU cannot take are computed
// on the processor, with the same energies as without -gpu_resident, and
// from then on the processor computes everything. The GPU is turned off to
// make it refuse the molecules, as without a GPU.
TEST_F(PropertyDevComputerTest, GpuResidentFallsBackToProcessor)
{
    MsgHandler msghandler;
    ForceField pd;
    readForceField(fileManager().getInputFilePath("forcefields/PG-pol-ions.xml"), &pd);
    ASSERT_TRUE(pd.polarizable());
    ForceComputer resident;
    resident.setGpuResident(true);
    resident.setKeepShells(false);
    const ForceComputer processor;
    const auto mpFile = fileManager().getInputFilePath("molprop/ion-water.xml");
    auto mols   = prepareMolecules(&msghandler, &pd, &processor, {}, mpFile);
    ASSERT_FALSE(mols.empty());
    std::vector<double> energies, reference;
    {
        const ProcessorOnly cpu;
        EXPECT_FALSE(resident.gpuCharges(&pd));
        EXPECT_FALSE(resident.prepareGpuResident(&pd, pointers(&mols)));
        energies  = computeAll(&resident, &pd, &mols, true);
        reference = computeAll(&processor, &pd, &mols, true);
    }
    EXPECT_EQ(reference, energies);
    // Also when the GPU could be used again, the processor computes
    EXPECT_FALSE(resident.prepareGpuResident(&pd, pointers(&mols)));
    EXPECT_FALSE(resident.gpuCharges(&pd));
    EXPECT_EQ(reference, computeAll(&resident, &pd, &mols, true));
    // Without polarization the GPU is not used either
    ForceField nonPolarizable;
    readForceField(fileManager().getInputFilePath("forcefields/ACS-g1.xml"), &nonPolarizable);
    ASSERT_FALSE(nonPolarizable.polarizable());
    ForceComputer residentNonPol;
    residentNonPol.setGpuResident(true);
    auto monomers = prepareMolecules(&msghandler, &nonPolarizable, &residentNonPol);
    ASSERT_FALSE(monomers.empty());
    EXPECT_FALSE(residentNonPol.prepareGpuResident(&nonPolarizable, pointers(&monomers)));
    EXPECT_FALSE(residentNonPol.gpuCharges(&nonPolarizable));
}

// Test: Several sets of parameters computed at once on the GPU, as replicas
// of the molecules, give the same energies and charges as computing each set
// by itself. The sets differ in the polarizabilities and in the
// electronegativities, of which the GPU generates the charges with SQE.
TEST_F(PropertyDevComputerTest, ReplicasAsSeparateEvaluationsGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && !GMX_DOUBLE, "no GPU or double precision"))
    {
        return;
    }
    ForceComputer forceComputer;
    forceComputer.setGpuResident(true);
    forceComputer.setKeepShells(false);
    // As many replicas as sets of parameters below
    forceComputer.setGpuReplicas(3);
    MsgHandler    msghandler;
    ForceField    pd;
    readForceField(fileManager().getInputFilePath("forcefields/PG-pol-ions.xml"), &pd);
    ASSERT_TRUE(pd.polarizable());
    const auto mpFile = fileManager().getInputFilePath("molprop/ion-water.xml");
    auto mols   = prepareMolecules(&msghandler, &pd, &forceComputer, {}, mpFile);
    ASSERT_FALSE(mols.empty());
    ASSERT_TRUE(forceComputer.prepareGpuResident(&pd, pointers(&mols))) << "The GPU cannot compute these molecules";
    std::vector<std::pair<ForceFieldParameter *, double>> original;
    for(auto itype : { InteractionType::POLARIZATION, InteractionType::ELECTRONEGATIVITYEQUALIZATION })
    {
        for(auto &param : *pd.findForces(itype)->parametersNoRemoval())
        {
            for(auto &pp : param.second)
            {
                if (pp.first == "alpha" || pp.first == "chi")
                {
                    original.emplace_back(&pp.second, pp.second.value());
                }
            }
        }
    }
    ASSERT_FALSE(original.empty());
    const std::vector<double> factors = { 1.0, 1.1, 0.95 };
    std::function<void(size_t)> useSet = [&original, &factors, &pd, &mols](size_t s)
    {
        for(auto &o : original)
        {
            o.first->forceSetValue(o.second*factors[s]);
        }
        fillAll(&pd, &mols);
    };
    const auto allCharges = [&mols]()
    {
        std::vector<double> q;
        for(const auto &mol : mols)
        {
            for(const auto &a : mol.atomsConst())
            {
                q.push_back(a.charge());
            }
        }
        return q;
    };
    // Each set by itself, all replicas with the same parameters
    std::vector<std::vector<double>> energies, charges;
    for(size_t s = 0; s < factors.size(); s++)
    {
        useSet(s);
        energies.push_back(computeAll(&forceComputer, &pd, &mols, false));
        charges.push_back(allCharges());
    }
    EXPECT_NE(energies[0], energies[1]);
    EXPECT_NE(energies[0], energies[2]);
    // All at once
    ASSERT_EQ(factors.size(), forceComputer.gpuReplicas());
    ASSERT_TRUE(forceComputer.computeRelaxedEnergiesReplicas(&msghandler, &pd, pointers(&mols),
                                                             useSet, false));
    EXPECT_TRUE(msghandler.ok());
    for(size_t s = 0; s < factors.size(); s++)
    {
        SCOPED_TRACE(gmx::formatString("set %zu", s));
        forceComputer.selectGpuReplica(s);
        auto e = allEnergies(&mols);
        auto q = allCharges();
        ASSERT_EQ(energies[s].size(), e.size());
        ASSERT_EQ(charges[s].size(), q.size());
        for(size_t i = 0; i < e.size(); i++)
        {
            EXPECT_NEAR(energies[s][i], e[i], 1e-6*(1 + std::abs(e[i]))) << "energy " << i;
        }
        for(size_t i = 0; i < q.size(); i++)
        {
            EXPECT_NEAR(charges[s][i], q[i], 1e-7) << "charge " << i;
        }
    }
}

// Test: The deviations of the energies computed from the energies as the GPU
// gives them, without the conformation data, on the processor and on the
// GPU, are those of calcDeviation, exactly, for all energy targets, loss
// functions and Boltzmann weighting.
TEST_F(PropertyDevComputerTest, CompactDeviationsAsCalcDeviationGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && !GMX_DOUBLE, "no GPU or double precision"))
    {
        return;
    }
    ForceComputer forceComputer;
    forceComputer.setGpuResident(true);
    forceComputer.setKeepShells(false);
    forceComputer.setKeepGpuEnergies(true);
    MsgHandler    msghandler;
    ForceField    pd;
    readForceField(fileManager().getInputFilePath("forcefields/PG-pol-ions.xml"), &pd);
    const auto mpFile = fileManager().getInputFilePath("molprop/ion-water.xml");
    auto mols   = prepareMolecules(&msghandler, &pd, &forceComputer, {}, mpFile);
    ASSERT_FALSE(mols.empty());
    ASSERT_TRUE(forceComputer.prepareGpuResident(&pd, pointers(&mols))) << "The GPU cannot compute these molecules";
    const std::vector<eRMS> all = { eRMS::EPOT, eRMS::Interaction, eRMS::Electrostatics,
                                    eRMS::Dispersion, eRMS::Exchange, eRMS::Induction,
                                    eRMS::DeltaHF, eRMS::ExchInd, eRMS::AllElec };
    for(auto lf : { LossFunction::MSE, LossFunction::MAE, LossFunction::Huber, LossFunction::Asinh })
    {
        for(double T : { 0.0, 300.0 })
        {
            SCOPED_TRACE(gmx::formatString("loss function %d temperature %g", static_cast<int>(lf), T));
            ForceEnergyDevComputer computer(sameTemperature({ eRMS::EPOT, eRMS::Interaction }, T));
            computer.setSeparateInductionCorrection(false);
            std::map<eRMS, FittingTarget> targets;
            for(auto erms : all)
            {
                targets.insert({ erms, FittingTarget(erms, iMolSelect::Train) });
                targets.find(erms)->second.setWeight(1);
            }
            auto compact = targets;
            auto onGpu   = targets;
            // The GPU computes the deviations with the energies
            ASSERT_TRUE(computer.setupGpuDeviations(&forceComputer, pointers(&mols), targets, lf));
            MsgHandler mh;
            computer.calcDeviation(&mh, &forceComputer, pointers(&mols), &targets, &pd, lf);
            EXPECT_TRUE(mh.ok());
            ASSERT_TRUE(computer.calcDeviationCompact(&forceComputer, pointers(&mols), &compact, lf));
            ASSERT_TRUE(computer.calcDeviationGpu(&forceComputer, &onGpu));
            for(auto erms : all)
            {
                const auto &a = targets.find(erms)->second;
                const auto &b = compact.find(erms)->second;
                const auto &g = onGpu.find(erms)->second;
                EXPECT_EQ(a.chiSquared(), b.chiSquared()) << rmsName(erms);
                EXPECT_EQ(a.totalWeight(), b.totalWeight()) << rmsName(erms);
                EXPECT_EQ(a.totalWeight(), g.totalWeight()) << rmsName(erms) << " on the GPU";
                if (lf == LossFunction::MSE || lf == LossFunction::MAE)
                {
                    EXPECT_EQ(a.chiSquared(), g.chiSquared()) << rmsName(erms) << " on the GPU";
                }
                else
                {
                    // The square root and logarithm of the GPU may differ in the last bit
                    EXPECT_NEAR(a.chiSquared(), g.chiSquared(), 1e-12*std::abs(a.chiSquared()))
                        << rmsName(erms) << " on the GPU";
                }
            }
            EXPECT_LT(0, targets.find(eRMS::AllElec)->second.totalWeight());
            // Only the deviations, with the energies on the GPU
            auto only = compact;
            for(auto &t : only)
            {
                t.second.reset();
            }
            ASSERT_TRUE(forceComputer.computeGpuDeviations(&mh, &pd, pointers(&mols)));
            ASSERT_TRUE(computer.calcDeviationGpu(&forceComputer, &only));
            for(auto erms : all)
            {
                EXPECT_EQ(onGpu.find(erms)->second.chiSquared(), only.find(erms)->second.chiSquared())
                    << rmsName(erms) << " deviations only";
                EXPECT_EQ(onGpu.find(erms)->second.totalWeight(), only.find(erms)->second.totalWeight())
                    << rmsName(erms) << " deviations only";
            }
        }
    }
}

// Test: Without polarization and forces, the interaction energies of dimers
// are computed from the intermolecular interactions only. They are the same
// as the dimer energies minus those of the monomers, computed with forces.
static void compareIntermolecular(const ForceField *pd, std::vector<ACTMol> *mols)
{
    ASSERT_FALSE(pd->polarizable());
    const ForceComputer forceComputer;
    std::vector<std::vector<std::map<InteractionType, std::vector<double>>>> energies;
    for(bool forces : { true, false })
    {
        MsgHandler mh;
        forceComputer.computeRelaxedEnergies(&mh, pd, pointers(mols), true, false, forces);
        EXPECT_TRUE(mh.ok());
        std::vector<std::map<InteractionType, std::vector<double>>> e;
        for(auto &mol : *mols)
        {
            std::map<InteractionType, std::vector<double>> emol;
            if (mol.isDimer())
            {
                EXPECT_TRUE(mol.topology()->hasIntermolecular());
                auto &ener = mol.conformationData(2)->energies;
                for(auto t : ener.types())
                {
                    const auto term = ener.find(t);
                    for(size_t c = 0; c < ener.nconf(); c++)
                    {
                        emol[t].push_back(term.has(c) ? term.act[c] : 0.0);
                    }
                }
            }
            e.push_back(emol);
        }
        energies.push_back(e);
    }
    size_t ndimer = 0;
    for(size_t m = 0; m < mols->size(); m++)
    {
        SCOPED_TRACE((*mols)[m].getMolname());
        auto &full  = energies[0][m];
        auto &inter = energies[1][m];
        ndimer += full.empty() ? 0 : 1;
        for(const auto &f : full)
        {
            const auto i = inter.find(f.first);
            ASSERT_TRUE(inter.end() != i) << interactionTypeToString(f.first);
            ASSERT_EQ(f.second.size(), i->second.size());
            for(size_t c = 0; c < f.second.size(); c++)
            {
                // In single precision the energy of each interaction is
                // rounded differently with and without forces, see
                // calcCoulombKernel, which changes the sums of hundreds of
                // Coulomb interactions by up to about 1e-3 kJ/mol
                double tolerance = (GMX_DOUBLE ? 1e-8*(1 + std::abs(f.second[c])) :
                                    5e-3 + 1e-5*std::abs(f.second[c]));
                EXPECT_NEAR(f.second[c], i->second[c], tolerance)
                    << interactionTypeToString(f.first) << " conformation " << c;
            }
        }
    }
    EXPECT_LT(0u, ndimer);
}

TEST_F(PropertyDevComputerTest, IntermolecularEnergies)
{
    ForceField    pd;
    readForceField(fileManager().getInputFilePath("forcefields/PC+GV-elec.xml"), &pd);
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto mols = prepareMolecules(&msghandler, &pd, &forceComputer, {},
                                 fileManager().getInputFilePath("molprop/mpdimer.xml"));
    compareIntermolecular(&pd, &mols);
}

TEST_F(PropertyDevComputerTest, IntermolecularEnergiesACSg2)
{
    ForceField    pd = *getForceField("ACS-g2");
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto mols = prepareMolecules(&msghandler, &pd, &forceComputer, {},
                                 fileManager().getInputFilePath("molprop/mpdimer.xml"));
    compareIntermolecular(&pd, &mols);
}

// Test: The same with a virtual site, that only interacts in the dimer
TEST_F(PropertyDevComputerTest, IntermolecularEnergiesVsite)
{
    ForceField    pd = *getForceField("ACS-v3fds");
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto mols = prepareMolecules(&msghandler, &pd, &forceComputer, { "water#1-butanol" },
                                 fileManager().getInputFilePath("molprop/water-butanol.xml"));
    compareIntermolecular(&pd, &mols);
}

// Test: Energies used again are the same as computed from scratch
TEST_F(PropertyDevComputerTest, EnergiesUsedAgain)
{
    testEnergiesUsedAgain("ACS-g2", "", {}, false);
}

// Test: The same with a dimer with a virtual site, of which the parameter changes.
// Within water the virtual site has no interactions, in the dimer it has.
TEST_F(PropertyDevComputerTest, EnergiesUsedAgainVsite)
{
    testEnergiesUsedAgain("ACS-v3fds", fileManager().getInputFilePath("molprop/water-butanol.xml"),
                          { "water#1-butanol" }, true);
}

// Test: Without a BOUNDS target nothing is computed
TEST_F(PropertyDevComputerTest, BoundsNoTarget)
{
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto          pd   = getForceField("ACS-g2");
    auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
    BoundsDevComputer computer(nullptr);

    std::map<eRMS, FittingTarget> targets;
    targets.insert({eRMS::EPOT, FittingTarget(eRMS::EPOT, iMolSelect::Train)});
    targets.find(eRMS::EPOT)->second.setWeight(1);
    computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols),
                           &targets, pd, LossFunction::MSE);
    EXPECT_EQ(1u, targets.size());
    EXPECT_EQ(0, targets.find(eRMS::EPOT)->second.totalWeight());
    EXPECT_EQ(0, targets.find(eRMS::EPOT)->second.chiSquared());
}

// Test: ACM charges without a range (minimum == maximum) are not checked
TEST_F(PropertyDevComputerTest, BoundsNoRange)
{
    for(auto chi2 : testBounds("ACS-g2", false, 0, 0))
    {
        EXPECT_EQ(0, chi2);
    }
}

// Test: Charges inside wide bounds give zero deviation
TEST_F(PropertyDevComputerTest, BoundsWide)
{
    for(auto chi2 : testBounds("ACS-g2", true, -100, 100))
    {
        EXPECT_EQ(0, chi2);
    }
}

// Test: Charges below the minimum
TEST_F(PropertyDevComputerTest, BoundsBelowMinimum)
{
    for(auto chi2 : testBounds("ACS-g2", true, 50, 100))
    {
        EXPECT_GT(chi2, 0);
    }
}

// Test: Charges above the maximum
TEST_F(PropertyDevComputerTest, BoundsAboveMaximum)
{
    for(auto chi2 : testBounds("ACS-g2", true, -100, -50))
    {
        EXPECT_GT(chi2, 0);
    }
}

// Test: Shell charges are added to the core charge
TEST_F(PropertyDevComputerTest, BoundsPolarizable)
{
    for(auto chi2 : testBounds("ACS-pg1", true, 50, 100))
    {
        EXPECT_GT(chi2, 0);
    }
}

// Test: Dipole, non-polarizable
TEST_F(PropertyDevComputerTest, MultipoleDipole)
{
    testMultipole(MolPropObservable::DIPOLE, eRMS::MU, "ACS-g2");
}

// Test: Quadrupole, non-polarizable
TEST_F(PropertyDevComputerTest, MultipoleQuadrupole)
{
    testMultipole(MolPropObservable::QUADRUPOLE, eRMS::QUAD, "ACS-g2");
}

// Test: Octupole, non-polarizable
TEST_F(PropertyDevComputerTest, MultipoleOctupole)
{
    testMultipole(MolPropObservable::OCTUPOLE, eRMS::OCT, "ACS-g2");
}

// Test: Hexadecapole, non-polarizable
TEST_F(PropertyDevComputerTest, MultipoleHexadecapole)
{
    testMultipole(MolPropObservable::HEXADECAPOLE, eRMS::HEXADEC, "ACS-g2");
}

// Test: Dipole, polarizable, shells are relaxed
TEST_F(PropertyDevComputerTest, MultipoleDipolePol)
{
    testMultipole(MolPropObservable::DIPOLE, eRMS::MU, "ACS-pg1");
}

// Test: Quadrupole, polarizable, shells are relaxed
TEST_F(PropertyDevComputerTest, MultipoleQuadrupolePol)
{
    testMultipole(MolPropObservable::QUADRUPOLE, eRMS::QUAD, "ACS-pg1");
}

// Test: Octupole, polarizable, shells are relaxed
TEST_F(PropertyDevComputerTest, MultipoleOctupolePol)
{
    testMultipole(MolPropObservable::OCTUPOLE, eRMS::OCT, "ACS-pg1");
}

// Test: Hexadecapole, polarizable, shells are relaxed
TEST_F(PropertyDevComputerTest, MultipoleHexadecapolePol)
{
    testMultipole(MolPropObservable::HEXADECAPOLE, eRMS::HEXADEC, "ACS-pg1");
}

// Test: An observable that is not a multipole throws
TEST_F(PropertyDevComputerTest, MultipoleUnsupportedThrows)
{
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto          pd   = getForceField("ACS-g2");
    auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
    ASSERT_FALSE(mols.empty());
    MultiPoleDevComputer computer(MolPropObservable::POLARIZABILITY);
    std::map<eRMS, FittingTarget> targets;
    targets.insert({eRMS::Polar, FittingTarget(eRMS::Polar, iMolSelect::Train)});
    const std::vector<gmx::RVec> coords = mols[0].xOriginal();
    EXPECT_THROW(computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols),
                                        &targets, pd, LossFunction::MSE),
                 gmx::InvalidInputError);
}

// Test: The force deviation of dimers uses the atoms only
TEST_F(PropertyDevComputerTest, ForceDeviationDimer)
{
    testDimerForces("ACS-g2");
}

// Test: The same with shells, which are left out
TEST_F(PropertyDevComputerTest, ForceDeviationDimerPol)
{
    testDimerForces("ACS-pol");
}

// Test: Calling the deviation computer again with the same targets adds
// the same contribution, also when the reference data is used again
TEST_F(PropertyDevComputerTest, ForceEnergyAccumulates)
{
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto          pd   = getForceField("ACS-pol");
    auto          mols = prepareMolecules(&msghandler, pd, &forceComputer, {},
                                          fileManager().getInputFilePath("molprop/mpdimer.xml"));
    ForceEnergyDevComputer computer({});
    std::map<eRMS, FittingTarget> targets;
    for(auto erms : dimerTargets)
    {
        targets.insert({erms, FittingTarget(erms, iMolSelect::Train)});
        targets.find(erms)->second.setWeight(1);
    }
    computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols), &targets, pd, LossFunction::MSE);
    std::map<eRMS, std::pair<double, double>> once;
    for(const auto &t : targets)
    {
        once.insert({t.first, { t.second.totalWeight(), t.second.chiSquared() }});
        EXPECT_LT(0, t.second.totalWeight()) << rmsName(t.first);
    }
    computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols), &targets, pd, LossFunction::MSE);
    for(const auto &t : targets)
    {
        const auto o = once[t.first];
        EXPECT_NEAR(2*o.first, t.second.totalWeight(), 1e-12*o.first) << rmsName(t.first);
        EXPECT_NEAR(2*o.second, t.second.chiSquared(), 1e-10*std::abs(o.second)) << rmsName(t.first);
    }
}

// Test: Energies that are not finite are left out, but only for the terms
// that are not finite
TEST_F(PropertyDevComputerTest, NonFiniteEnergiesIgnored)
{
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto          pd   = getForceField("ACS-g2");
    auto          mols = prepareMolecules(&msghandler, pd, &forceComputer, {},
                                          fileManager().getInputFilePath("molprop/mpdimer.xml"));
    ASSERT_FALSE(mols.empty());
    // Van der Waals parameters that give no finite energies
    auto itype = InteractionType::VDW;
    auto nparam = potprops.find(pd->potential(itype))->second.param.size();
    const std::vector<double> nan(nparam, std::numeric_limits<double>::quiet_NaN());
    for(auto &mol : mols)
    {
        std::vector<Topology *> tops = { mol.topologyPtr() };
        for(auto &top : mol.fragmentHandler()->topologiesPtr())
        {
            tops.push_back(&top);
        }
        for(auto top : tops)
        {
            if (top->nonBondedPairs().hasType(itype))
            {
                top->setNonBondedParameters(itype, nan);
            }
        }
    }
    ForceEnergyDevComputer computer({});
    std::map<eRMS, FittingTarget> targets;
    for(auto erms : { eRMS::Interaction, eRMS::Electrostatics })
    {
        targets.insert({erms, FittingTarget(erms, iMolSelect::Train)});
        targets.find(erms)->second.setWeight(1);
    }
    computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols), &targets, pd, LossFunction::MSE);
    auto &tinter = targets.find(eRMS::Interaction)->second;
    EXPECT_EQ(0, tinter.totalWeight());
    EXPECT_EQ(0, tinter.chiSquared());
    auto &telec = targets.find(eRMS::Electrostatics)->second;
    EXPECT_LT(0, telec.totalWeight());
    EXPECT_TRUE(std::isfinite(telec.chiSquared()));
}

// Test: The bounds of an ACM charge apply to the total charge of the
// particle, its second charge and the charges of its shells included, see
// BoundsDevComputer::moleculeBound. With -lockstep, the second charges were
// once those of another individual.
TEST_F(PropertyDevComputerTest, BoundsShellsAndSecondCharge)
{
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    ForceField    pd;
    readForceField(fileManager().getInputFilePath("forcefields/ACS-pol.xml"), &pd);
    // The ACM charge of the carbons, which have shells, gets a range
    auto qparam = pd.findParticleType("c3")->parameter("charge");
    ASSERT_EQ(Mutability::ACM, qparam->mutability());
    qparam->setMaximum(qparam->minimum() + 1);
    double qmin = qparam->minimum();
    double qmax = qparam->maximum();
    auto   mols = prepareMolecules(&msghandler, &pd, &forceComputer, { "water#1-butanol" },
                                   fileManager().getInputFilePath("molprop/mpdimer.xml"));
    ASSERT_EQ(1u, mols.size());
    auto  &mol   = mols[0];
    auto  *atoms = mol.atoms();
    size_t core  = atoms->size();
    for (size_t i = 0; i < atoms->size() && core == atoms->size(); i++)
    {
        if ((*atoms)[i].ffType() == "c3" && !(*atoms)[i].shells().empty())
        {
            core = i;
        }
    }
    ASSERT_LT(core, atoms->size());
    // All charges within their bounds, the ACM charges in the middle
    for (auto &a : *atoms)
    {
        a.setChargeVS(0, a.zetaVS(), a.rowVS());
        a.setCharge(0);
        const auto &p = pd.findParticleType(a.ffType())->parameterConst("charge");
        if (p.mutability() == Mutability::ACM)
        {
            a.setCharge(0.5*(p.minimum() + p.maximum()));
        }
    }
    const auto bound = [&mol, &pd](LossFunction lf)
    {
        return BoundsDevComputer::moleculeBound(&mol, &pd, lf);
    };
    EXPECT_EQ(0, bound(LossFunction::MSE));
    // The second charge makes the total too large
    auto  &c     = (*atoms)[core];
    auto  &shell = (*atoms)[c.shells()[0]];
    double d     = 0.25;
    c.setChargeVS(qmax - c.charge() + d, c.zetaVS(), c.rowVS());
    EXPECT_DOUBLE_EQ(d*d, bound(LossFunction::MSE));
    EXPECT_DOUBLE_EQ(d, bound(LossFunction::MAE));
    // The shell brings it back to the maximum
    shell.setCharge(-d);
    EXPECT_EQ(0, bound(LossFunction::MSE));
    // and below the minimum
    shell.setCharge(qmin - qmax - 2*d);
    EXPECT_DOUBLE_EQ(d*d, bound(LossFunction::MSE));
    // The deviation adds the bound of each molecule
    std::map<eRMS, FittingTarget> targets;
    targets.insert({ eRMS::BOUNDS, FittingTarget(eRMS::BOUNDS, iMolSelect::Train) });
    targets.find(eRMS::BOUNDS)->second.setWeight(1);
    BoundsDevComputer computer(nullptr);
    computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols), &targets, &pd, LossFunction::MSE);
    EXPECT_DOUBLE_EQ(d*d, targets.find(eRMS::BOUNDS)->second.chiSquared());
}

// Test: ESP deviation, non-polarizable
TEST_F(PropertyDevComputerTest, EspNonPolarizable)
{
    testEsp("ACS-g2");
}

// Changing a cloud width on an existing molecule must give the same ESP
// as a newly initialized RESP object. No parameter-selection flag is needed.
// Exercise both retained virtual sites and converted second charges.
TEST_F(PropertyDevComputerTest, EspCloudWidthsStayCurrent)
{
    for(auto distribution : { ChargeDistributionType::Gaussian, ChargeDistributionType::Slater })
    {
        for(int row : { 1, 2 })
        {
            for(bool converted : { false, true })
            {
                SCOPED_TRACE(::testing::Message() << "distribution " << static_cast<int>(distribution)
                                                   << " row " << row << " converted " << converted);
                ForceField pd;
                readForceField(fileManager().getInputFilePath("forcefields/PC+GV-elec.xml"),
                                 &pd, nullptr, converted);
                pd.setChargeDistribution(distribution);
                auto cloud = pd.findParticleType(converted ? "ow" : "v1ow");
                cloud->removeOption(converted ? "rowVS" : "row");
                cloud->setOption(converted ? "rowVS" : "row", std::to_string(row));
                ASSERT_EQ(row, converted ? cloud->rowVS() : cloud->row());
                const char *widthName = converted ? "zetaVS" : "zeta";
                auto width = cloud->parameter(widthName);
                width->forceSetValue(10);
                const ForceComputer forceComputer;
                MsgHandler msghandler;
                const auto waterTypes = [](MolProp *mp)
                {
                    for(auto &experiment : *mp->experiment())
                    {
                        for(auto &atom : *experiment.calcAtom())
                        {
                            atom.setObtype(atom.getObtype() == "o3" ? "ow" : "hw");
                        }
                    }
                };
                auto mols = prepareMolecules(&msghandler, &pd, &forceComputer, { "water" }, "", waterTypes);
                ASSERT_EQ(1u, mols.size());
                auto &mol = mols.front();
                ACTQprop grid(mol.atomsConst(), mol.xOriginal());
                auto resp = grid.qgenResp();
                resp->setChargeDistributionType(distribution);
                resp->setAtomInfo(&msghandler, mol.atomsConst(), &pd, mol.totalCharge());
                resp->updateAtomCoords(mol.xOriginal());
                resp->makeGrid(0.15, 0.25, mol.xOriginal());
                ASSERT_GT(resp->nEsp(), 0u);
                mol.qProps()->clear();
                mol.qProps()->push_back(std::move(grid));
                EspDevComputer computer;
                std::vector<double> initial;
                const auto initialAtoms = mol.atomsConst();
                for(double value : { 10.0, 20.0, 10.0 })
                {
                    width->forceSetValue(value);
                    auto coords = mol.xOriginal();
                    std::vector<gmx::RVec> forces(mol.atomsConst().size());
                    mol.generateCharges(&msghandler, &pd, &forceComputer,
                                          pd.chargeGenerationAlgorithm(), &coords, &forces);
                    std::map<eRMS, FittingTarget> targets;
                    targets.insert({eRMS::ESP, FittingTarget(eRMS::ESP, iMolSelect::Train)});
                    targets.find(eRMS::ESP)->second.setWeight(1);
                    computer.calcDeviation(&msghandler, &forceComputer, { &mol }, &targets, &pd,
                                             LossFunction::MSE);
                    auto cached = mol.qProps()->front().qgenResp();
                    QgenResp fresh;
                    fresh.setChargeDistributionType(distribution);
                    fresh.setAtomInfo(&msghandler, mol.atomsConst(), &pd, mol.totalCharge());
                    fresh.updateAtomCoords(cached->coords());
                    fresh.updateAtomCharges(&msghandler, mol.atomsConst());
                    for(const auto &point : cached->espPoints())
                    {
                        auto x = point.esp();
                        fresh.addEspPoint(x[XX], x[YY], x[ZZ], point.v());
                    }
                    fresh.calcPot(&msghandler, pd.epsilonr());
                    double difference = 0;
                    for(size_t i = 0; i < cached->nEsp(); i++)
                    {
                        double potential = cached->espPoint(i).vCalc();
                        EXPECT_NEAR(potential, fresh.espPoint(i).vCalc(),
                                    gmx::test::precisionTolerance(1e-10, 1e-6)*(1 + std::abs(potential))) << i;
                        if (initial.size() == cached->nEsp())
                        {
                            difference += std::abs(potential - initial[i]);
                        }
                    }
                    if (initial.empty())
                    {
                        for(const auto &point : cached->espPoints())
                        {
                            initial.push_back(point.vCalc());
                        }
                    }
                    else if (value == 20)
                    {
                        EXPECT_GT(difference, 1e-5);
                        // Reintroduce a stale RESP width deliberately. Cube
                        // output must refresh it even without an ESP objective.
                        width->forceSetValue(10);
                        cached->updateZeta(initialAtoms, &pd);
                        width->forceSetValue(value);
                        const auto cubeFile = fileManager().getTemporaryFilePath(
                                  "width-" + std::to_string(static_cast<int>(distribution)) + "-" +
                                  std::to_string(row) + "-" + std::to_string(converted) + ".cube");
                        mol.GenerateCube(&msghandler, &pd, mol.xOriginal(), &forceComputer,
                                           0.15, 0.25, nullptr, nullptr, nullptr, cubeFile.c_str(),
                                           nullptr, nullptr, nullptr, nullptr, nullptr);
                        QgenResp cube;
                        cube.setChargeDistributionType(distribution);
                        cube.setAtomInfo(&msghandler, mol.atomsConst(), &pd, mol.totalCharge());
                        cube.updateAtomCoords(mol.xOriginal());
                        cube.readCube(cubeFile, true);
                        cube.updateAtomCharges(&msghandler, mol.atomsConst());
                        cube.calcPot(&msghandler, pd.epsilonr());
                        ASSERT_GT(cube.nEsp(), 0u);
                        for(const auto &point : cube.espPoints())
                        {
                            // Cube coordinates and potentials are rounded.
                            EXPECT_NEAR(point.v(), point.vCalc(), 1e-4*(1 + std::abs(point.vCalc())));
                        }
                    }
                    else
                    {
                        EXPECT_NEAR(0, difference, gmx::test::precisionTolerance(1e-8, 1e-4));
                    }
                    EXPECT_TRUE(msghandler.ok());
                }
            }
        }
    }
}

// Test: ESP deviation, polarizable, shells are relaxed first
TEST_F(PropertyDevComputerTest, EspPolarizable)
{
    testEsp("ACS-pg1");
}

// Test: The polarizability deviation is the sum of the squared differences
// of the tensor elements in the unit of the polarizability, its weight the
// number of reference tensors
TEST_F(PropertyDevComputerTest, Polarizability)
{
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto          pd   = getForceField("ACS-pg1");
    ASSERT_TRUE(pd->polarizable());
    auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
    PolarDevComputer computer;
    double convert = convertFromGromacs(1.0, mpo_unit2(MolPropObservable::POLARIZABILITY));
    EXPECT_NE(1, convert);
    for(auto &mol : mols)
    {
        const std::vector<ACTMol *> one = { &mol };
        std::map<eRMS, FittingTarget> targets;
        targets.insert({eRMS::Polar, FittingTarget(eRMS::Polar, iMolSelect::Train)});
        targets.find(eRMS::Polar)->second.setWeight(1);
        computer.calcDeviation(&msghandler, &forceComputer, one, &targets, pd, LossFunction::MSE);
        auto  &t    = targets.find(eRMS::Polar)->second;
        size_t nref = 0;
        double sum2 = 0;
        auto   qProps = mol.qProps();
        for(auto qp = qProps->begin(); qp < qProps->end(); ++qp)
        {
            auto &qref = qp->qPqmConst();
            if (!qref.hasPolarizability())
            {
                continue;
            }
            nref += 1;
            auto aref  = qref.polarizabilityTensor();
            auto acalc = qp->qPact()->polarizabilityTensor();
            for(int i = 0; i < DIM; i++)
            {
                // The molecules are polarizable in every direction
                EXPECT_LT(0, acalc[i][i]) << mol.getMolname();
                for(int j = 0; j < DIM; j++)
                {
                    sum2 += gmx::square(convert*(aref[i][j] - acalc[i][j]));
                }
            }
        }
        EXPECT_LT(0u, nref) << mol.getMolname();
        EXPECT_DOUBLE_EQ(nref, t.totalWeight()) << mol.getMolname();
        EXPECT_NEAR(sum2, t.chiSquared(), 1e-10*sum2) << mol.getMolname();
    }
}

// Test: Deviation of the frequencies
TEST_F(PropertyDevComputerTest, HarmonicsFrequencies)
{
    testHarmonics(MolPropObservable::FREQUENCY, eRMS::FREQUENCY);
}

// Test: Deviation of the intensities
TEST_F(PropertyDevComputerTest, HarmonicsIntensities)
{
    testHarmonics(MolPropObservable::INTENSITY, eRMS::INTENSITY);
}

// Test: With Boltzmann weighting the potential energy deviation is the sum
// of the weights times the loss, as for the interaction energies. The
// scale of the Huber loss is not changed by the weight.
TEST_F(PropertyDevComputerTest, EpotBoltzmannHuber)
{
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto          pd   = getForceField("ACS-g2");
    auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
    double        T    = 500;
    ForceEnergyDevComputer computer({ { eRMS::EPOT, T } });
    std::map<eRMS, FittingTarget> targets;
    targets.insert({eRMS::EPOT, FittingTarget(eRMS::EPOT, iMolSelect::Train)});
    targets.find(eRMS::EPOT)->second.setWeight(1);
    computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols), &targets, pd, LossFunction::Huber);
    auto  &t    = targets.find(eRMS::EPOT)->second;
    double beta = 1.0/(BOLTZ*T);
    double wsum = 0, lsum = 0;
    for(auto &mol : mols)
    {
        const auto &coords = *mol.monomerCoords(0);
        const auto  epot   = mol.conformationData(0)->energies.find(InteractionType::EPOT);
        for(size_t i = 0; i < coords.size(); i++)
        {
            const auto ee = coords[i].eMap.find(InteractionType::EPOT);
            if (coords[i].eMap.end() == ee || !ee->second.haveQM() || !epot.has(i))
            {
                continue;
            }
            double eqm = ee->second.eqm();
            double w   = std::exp(-beta*(eqm - coords[i].eMin));
            wsum += w;
            lsum += w*loss(LossFunction::Huber, 1, eqm - epot.act[i]);
        }
    }
    EXPECT_LT(0, wsum);
    // Some weights are well below one
    EXPECT_LT(wsum + 0.5, t.totalWeight() + mols.size() + 1);
    EXPECT_NEAR(wsum, t.totalWeight(), 1e-10*wsum);
    EXPECT_NEAR(lsum, t.chiSquared(), 1e-10*lsum);
}

// Test: A molecule of which the minimization fails is left out, the
// frequencies of the following molecules are used
TEST_F(PropertyDevComputerTest, HarmonicsAfterFailedMinimization)
{
    const ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto          pd = getForceField("ACS-g2");
    // Reference frequencies, as many as each molecule has
    const auto addHarmonics = [](MolProp *mp)
    {
        size_t natom = (*mp->experiment())[0].NAtom();
        auto   mpo   = MolPropObservable::FREQUENCY;
        auto   harm  = std::make_unique<Harmonics>(mpo_unit2(mpo), 0, mpo);
        for(size_t k = 0; k < 3*natom - 6; k++)
        {
            harm->addValue(1000);
        }
        (*mp->experiment())[0].addProperty(mpo, std::move(harm));
    };
    auto mols = prepareMolecules(&msghandler, pd, &forceComputer, { "acetone", "water" }, "", addHarmonics);
    ASSERT_EQ(2u, mols.size());
    // Acetone first, its bonds can not be minimized
    if (mols[0].getMolname() != "acetone")
    {
        std::swap(mols[0], mols[1]);
    }
    auto entries = mols[0].topologyPtr()->entriesPtr();
    const auto bonds   = entries->find(InteractionType::BONDS);
    ASSERT_TRUE(entries->end() != bonds);
    for(auto &e : bonds->second)
    {
        e->setParams(std::vector<double>(e->params().size(), std::numeric_limits<double>::quiet_NaN()));
    }
    HarmonicsDevComputer computer(MolPropObservable::FREQUENCY);
    std::map<eRMS, FittingTarget> targets;
    targets.insert({eRMS::FREQUENCY, FittingTarget(eRMS::FREQUENCY, iMolSelect::Train)});
    targets.find(eRMS::FREQUENCY)->second.setWeight(1);
    computer.calcDeviation(&msghandler, &forceComputer, pointers(&mols), &targets, pd, LossFunction::MSE);
    auto &t = targets.find(eRMS::FREQUENCY)->second;
    // Only water
    EXPECT_EQ(1, t.totalWeight());
    EXPECT_TRUE(std::isfinite(t.chiSquared()));
}

} // namespace alexandria
