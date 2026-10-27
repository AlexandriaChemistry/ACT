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
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include "act/alexandria/devcomputer.h"
#include "act/alexandria/molgen.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/molprop/molprop_xml.h"
#include "act/topology/actmol.h"

#include "testutils/cmdlinetest.h"
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
        auto rb = b.find(ra.first);
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
    std::string monomerFile() { return fileManager().getInputFilePath("../../molprop/tests/molprop.xml"); }

    //! \return the dimer test file
    std::string dimerFile() { return fileManager().getInputFilePath("mpdimer.xml"); }

    //! \return the dimer test file with induction correction
    std::string deltaHFFile() { return fileManager().getInputFilePath("methanol-ethane-deltaHF.xml"); }

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
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, opts.tolerance);
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
        ForceComputer           forceComputer;
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
                checker->checkSequence(qq.begin(), qq.end(), label.c_str());
            }
            computer.calcDeviation(&msghandler,
                                   &forceComputer,
                                   &actmols,
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
                    std::string weight = gmx::formatString("%s-%s-weight",
                                                           actmols[0].getMolname().c_str(), rmsName(t.first));
                    checker->checkDouble(t.second.totalWeight(), weight.c_str());
                    std::string chi2 = gmx::formatString("%s-%s-chi2",
                                                         actmols[0].getMolname().c_str(), rmsName(t.first));
                    checker->checkDouble(t.second.chiSquared(), chi2.c_str());
                }
            }
        }
        return results;
    }

    /*! \brief Check the QM reference for induction against the molprop data
     * \param[in] mpFile  Molprop file with dimers
     * \param[in] sepIndCorr Whether to separate the induction correction
     * \return The number of energies with a non-zero induction correction
     */
    int checkInductionQM(const std::string &mpFile, bool sepIndCorr)
    {
        ForceComputer forceComputer;
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
            actmol.forceEnergyMaps(&msghandler, pd, &forceComputer, &forceMap, &energyMap,
                                   &interactionEnergyMap, &energyComponentMap, sepIndCorr);
            EXPECT_EQ(qmInduction.size(), interactionEnergyMap.size());
            for(const auto &iem : interactionEnergyMap)
            {
                auto ind = iem.find(InteractionType::INDUCTION);
                if (iem.end() == ind || !ind->second.haveQM())
                {
                    ADD_FAILURE() << "No QM induction for " << actmol.getMolname();
                    continue;
                }
                auto qm = qmInduction.find(ind->second.id());
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
                auto icorr = iem.find(InteractionType::INDUCTIONCORRECTION);
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
    static void checkConformations(ACTMol                          &actmol,
                                   gmx::test::TestReferenceChecker *checker)
    {
        const auto *coords = actmol.isDimer() ? actmol.dimerCoords() : actmol.monomerCoords(0);
        EXPECT_FALSE(coords->empty()) << actmol.getMolname();
        for(size_t i = 0; i < coords->size(); i++)
        {
            const auto &ce    = (*coords)[i];
            auto        label = gmx::formatString("%s-conf%zu", actmol.getMolname().c_str(), i);
            for(const auto itype : { InteractionType::EXCHANGE,
                                     InteractionType::DISPERSION,
                                     InteractionType::EPOT })
            {
                auto ee = ce.eMap.find(itype);
                ASSERT_TRUE(ce.eMap.end() != ee && ee->second.haveACT())
                    << label << " " << interactionTypeToString(itype);
                checker->checkDouble(ee->second.eact(),
                                     gmx::formatString("%s-%s", label.c_str(),
                                                       interactionTypeToString(itype).c_str()).c_str());
            }
            // Repulsion and attraction
            EXPECT_GT(ce.eMap.find(InteractionType::EXCHANGE)->second.eact(), 0) << label;
            EXPECT_LT(ce.eMap.find(InteractionType::DISPERSION)->second.eact(), 0) << label;
            std::vector<double> f;
            for(const auto &fi : ce.fcalc)
            {
                for(int m = 0; m < DIM; m++)
                {
                    EXPECT_TRUE(std::isfinite(fi[m])) << label;
                    f.push_back(fi[m]);
                }
            }
            checker->checkSequence(f.begin(), f.end(), (label + "-forces").c_str());
        }
    }

    //! Run monomer and dimer tests with a loss function
    void testLossFunction(LossFunction lf)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-11);
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
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Constructor with empty Boltzmann map
TEST_F(ForceEnergyDevComputerTest, ConstructorWithEmptyBoltzmannMap)
{
    std::map<eRMS, double> boltzmann;
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_TRUE(computer.separateInductionCorrection());
}

// Test: Constructor with non-empty Boltzmann map
TEST_F(ForceEnergyDevComputerTest, ConstructorWithNonEmptyBoltzmannMap)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 298.15;
    boltzmann[eRMS::Force2] = 0.0;
    
    ForceEnergyDevComputer computer(boltzmann);
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
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Boltzmann temperature with Force2
TEST_F(ForceEnergyDevComputerTest, BoltzmannTemperatureForce2)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::Force2] = 300.0;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Boltzmann temperature with multiple eRMS terms
TEST_F(ForceEnergyDevComputerTest, BoltzmannTemperatureMultipleTerms)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 298.15;
    boltzmann[eRMS::Force2] = 0.0;
    boltzmann[eRMS::Interaction] = 500.0;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
    EXPECT_TRUE(computer.separateInductionCorrection());
}

// Test: Inheritance from DevComputer
TEST_F(ForceEnergyDevComputerTest, InheritsFromDevComputer)
{
    std::map<eRMS, double> boltzmann;
    ForceEnergyDevComputer computer(boltzmann);
    
    // Verify that DevComputer interface is accessible
    const std::string& name = computer.name();
    EXPECT_EQ("ForceEnergy", name);
}

// Test: Zero Boltzmann temperature
TEST_F(ForceEnergyDevComputerTest, ZeroBoltzmannTemperature)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 0.0;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Large Boltzmann temperature values
TEST_F(ForceEnergyDevComputerTest, LargeBoltzmannTemperatureValues)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 1e6;
    boltzmann[eRMS::Force2] = 1e10;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Negative Boltzmann temperature (edge case)
TEST_F(ForceEnergyDevComputerTest, NegativeBoltzmannTemperature)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = -100.0;
    
    ForceEnergyDevComputer computer(boltzmann);
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
    std::string mpFile = fileManager().getInputFilePath("../../molprop/tests/molprop.xml");
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
    testCalcDev(monomerFile(), "../../forces/tests/ACS-v3fds", opts);
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
    std::string mpFile = fileManager().getInputFilePath("mpdimer.xml");
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
    std::string mpFile = fileManager().getInputFilePath("../../molprop/tests/molprop.xml");
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
    std::string mpFile = fileManager().getInputFilePath("mpdimer.xml");
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
    std::string mpFile = fileManager().getInputFilePath("mpdimer.xml");
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
    std::string mpFile = fileManager().getInputFilePath("mpdimer.xml");
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
    std::string mpFile = fileManager().getInputFilePath("methanol-ethane-deltaHF.xml");
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
    std::string mpFile = fileManager().getInputFilePath("methanol-ethane-deltaHF.xml");
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
    auto repot   = computeDev(monomerFile(), "ACS-g2", epot, nullptr);
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
    auto rinter   = computeDev(dimerFile(), "ACS-g2", inter, nullptr);
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
        auto results = computeDev(mpFile, "ACS-g2", opts, nullptr);
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
    auto rcold = computeDev(monomerFile(), "ACS-g2", cold, nullptr);
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
    auto rcold = computeDev(dimerFile(), "ACS-g2", cold, nullptr);
    ASSERT_EQ(rref.size(), rcold.size());
    for(const auto &r : rcold)
    {
        EXPECT_TRUE(std::isfinite(r.second.first)) << r.first;
        EXPECT_TRUE(std::isfinite(r.second.second)) << r.first;
        EXPECT_GT(r.second.first, 0) << r.first;
        EXPECT_LT(r.second.first, rref[r.first].first) << r.first;
    }
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
class PropertyDevComputerTest : public gmx::test::CommandLineTestBase
{
protected:
    PropertyDevComputerTest() = default;

    //! \return the monomer test file
    std::string monomerFile() { return fileManager().getInputFilePath("../../molprop/tests/molprop.xml"); }

    //! Read molecules, generate topologies and charges like MolGen does
    std::vector<ACTMol> prepareMolecules(MsgHandler          *msghandler,
                                         ForceField          *pd,
                                         const ForceComputer *forceComputer)
    {
        std::vector<alexandria::MolProp> mps;
        MolPropRead(msghandler, monomerFile().c_str(), &mps);
        std::vector<ACTMol> mols;
        for(auto &mp : mps)
        {
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
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-11);
        checker_.setDefaultTolerance(tolerance);
        ForceComputer forceComputer;
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
        auto newTargets = []()
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
            std::vector<ACTMol> one;
            one.push_back(std::move(mol));
            auto targets = newTargets();
            computer.calcDeviation(&msghandler, &forceComputer, &one,
                                   &targets, &pdBounds, LossFunction::MSE);
            mol = std::move(one[0]);
            auto &t = targets.find(eRMS::BOUNDS)->second;
            EXPECT_EQ(1, t.totalWeight());
            checker_.checkDouble(t.chiSquared(),
                                 gmx::formatString("%s-Bounds-chi2", mol.getMolname().c_str()).c_str());
            chi2.push_back(t.chiSquared());
        }
        // All molecules at once give the sum
        auto targets = newTargets();
        computer.calcDeviation(&msghandler, &forceComputer, &mols,
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
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-8);
        checker_.setDefaultTolerance(tolerance);
        ForceComputer forceComputer;
        MsgHandler    msghandler;
        auto          pd   = getForceField(ffFile);
        auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
        MultiPoleDevComputer computer(mpo);
        EXPECT_EQ(mpo_name(mpo), computer.name());
        auto newTargets = [erms]()
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
            std::vector<ACTMol> one;
            one.push_back(std::move(mol));
            auto targets = newTargets();
            computer.calcDeviation(&msghandler, &forceComputer, &one,
                                   &targets, pd, LossFunction::MSE);
            mol = std::move(one[0]);
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
        computer.calcDeviation(&msghandler, &forceComputer, &mols,
                               &targets, pd, LossFunction::MSE);
        auto &t = targets.find(erms)->second;
        EXPECT_DOUBLE_EQ(weightSum, t.totalWeight());
        EXPECT_NEAR(chi2Sum, t.chiSquared(), 1e-10*(chi2Sum + 1));
    }
};

// Test: Without a BOUNDS target nothing is computed
TEST_F(PropertyDevComputerTest, BoundsNoTarget)
{
    ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto          pd   = getForceField("ACS-g2");
    auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
    BoundsDevComputer computer(nullptr);

    std::map<eRMS, FittingTarget> targets;
    targets.insert({eRMS::EPOT, FittingTarget(eRMS::EPOT, iMolSelect::Train)});
    targets.find(eRMS::EPOT)->second.setWeight(1);
    computer.calcDeviation(&msghandler, &forceComputer, &mols,
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
    for(auto chi2 : testBounds("ACS-pg", true, 50, 100))
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
    testMultipole(MolPropObservable::DIPOLE, eRMS::MU, "ACS-pg");
}

// Test: Quadrupole, polarizable, shells are relaxed
TEST_F(PropertyDevComputerTest, MultipoleQuadrupolePol)
{
    testMultipole(MolPropObservable::QUADRUPOLE, eRMS::QUAD, "ACS-pg");
}

// Test: Octupole, polarizable, shells are relaxed
TEST_F(PropertyDevComputerTest, MultipoleOctupolePol)
{
    testMultipole(MolPropObservable::OCTUPOLE, eRMS::OCT, "ACS-pg");
}

// Test: Hexadecapole, polarizable, shells are relaxed
TEST_F(PropertyDevComputerTest, MultipoleHexadecapolePol)
{
    testMultipole(MolPropObservable::HEXADECAPOLE, eRMS::HEXADEC, "ACS-pg");
}

// Test: An observable that is not a multipole throws
TEST_F(PropertyDevComputerTest, MultipoleUnsupportedThrows)
{
    ForceComputer forceComputer;
    MsgHandler    msghandler;
    auto          pd   = getForceField("ACS-g2");
    auto          mols = prepareMolecules(&msghandler, pd, &forceComputer);
    ASSERT_FALSE(mols.empty());
    MultiPoleDevComputer computer(MolPropObservable::POLARIZABILITY);
    std::map<eRMS, FittingTarget> targets;
    targets.insert({eRMS::Polar, FittingTarget(eRMS::Polar, iMolSelect::Train)});
    std::vector<gmx::RVec> coords = mols[0].xOriginal();
    EXPECT_THROW(computer.calcDeviation(&msghandler, &forceComputer, &mols,
                                        &targets, pd, LossFunction::MSE),
                 gmx::InvalidInputError);
}

} // namespace alexandria
