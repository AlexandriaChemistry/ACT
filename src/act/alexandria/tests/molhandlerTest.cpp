/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022-2026
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

#include "actpre.h"

#include <cmath>

#include <map>

#include <gtest/gtest.h>

#include "act/import/import.h"
#include "act/alexandria/molhandler.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/properties/thermochemistry.h"
#include "act/qgen/qgen_acm.h"
#include "act/topology/actmol.h"
#include "act/utility/units.h"
#include "gromacs/utility/fatalerror.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! \brief Add energies to a checker
static void add_energies(const ForceField                        *pd,
                         gmx::test::TestReferenceChecker         *checker,
                         const std::map<InteractionType, double> &energies,
                         const char                              *label)
{
       for(auto &imf : energies)
    {
        const std::string fname   = interactionTypeToString(imf.first);
        const std::string mylabel = gmx::formatString("%s %s", fname.c_str(), label);

        checker->checkReal(imf.second, mylabel.c_str());
    }
}

class MolHandlerTest : public gmx::test::CommandLineTestBase
{
protected:
    void test(const char *molname, const char *forcefield, bool nma, double ftoler=1e-15)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        const auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        const char                      *conf     = (char *)"minimum";
        std::string                      dataName;
        std::vector<alexandria::MolProp> molprops;
        
        // Get forcefield
        auto pd         = getForceField(forcefield);
        std::string moldir("molecules/");
        moldir += molname;
        dataName        = gmx::test::TestFileManager::getInputFilePath(moldir);
        double qtot     = 0;
        bool   userqtot = false;
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        importFile(&msghandler, pd, dataName, &molprops,
                   conf, JobType::OPT, userqtot, &qtot, true);
        EXPECT_TRUE(msghandler.ok());
        std::vector<ACTMol> mps;
        // Needed for GenerateCharges
        auto alg = pd->chargeGenerationAlgorithm();
        double shellTolerance = ftoler;
        int    shellMaxIter   = 100;
        double epsilonr = 1;
        if (!ffOption(*pd, InteractionType::ELECTROSTATICS, 
                      "epsilonr", &epsilonr))
        {
            epsilonr = 1;
        }
        auto forceComp = new ForceComputer(shellTolerance, shellMaxIter, epsilonr);
        const std::vector<double> qcustom;
        for(auto &molprop: molprops)
        {
            ACTMol mm;
            mm.Merge(&molprop);
            // Generate charges and topology
            mm.GenerateTopology(&msghandler, pd,
                                missingParameters::Ignore);
            EXPECT_TRUE(msghandler.ok());
            if (!msghandler.ok())
            {
                return;
            }
            std::vector<gmx::RVec> forces(mm.atomsConst().size());
            std::vector<gmx::RVec> coords = mm.xOriginal();
            mm.generateCharges(&msghandler, pd, forceComp, alg, &coords, &forces);
            if (msghandler.ok())
            {
                mps.push_back(std::move(mm));
            }
        }
        
        for(auto &mp : mps)
        {
            std::vector<gmx::RVec> forces(mp.atomsConst().size());
            std::vector<gmx::RVec> coords = mp.xOriginal();
            std::map<InteractionType, double> eBefore;
            forceComp->compute(&msghandler, pd, mp.topology(), &coords, &forces, &eBefore);
            add_energies(pd, &checker_, eBefore, "before");
            
            const MolHandler mh;
            
            std::vector<gmx::RVec> xmin = coords;
            double rmsd = mh.coordinateRmsd(mp.atomsConst(), coords, &xmin);
            checker_.checkReal(rmsd, "Coordinate RMSD before minimizing");
            // Reset coordinates, this should not be necessary though.
            xmin = coords;
            // Infinite number of shell iterations, i.e. until convergence.
            std::map<InteractionType, double> eAfter;
            SimulationConfigHandler simConfig;
            simConfig.setForceTolerance(ftoler);
            // Change nullptr to stdout for debugging.
            auto eMin = mh.minimizeCoordinates(&msghandler, pd, mp.topology(),
                                               forceComp,
                                               simConfig, &xmin, &eAfter, {});
            EXPECT_TRUE(eMinimizeStatus::OK == eMin);
            if (eMinimizeStatus::OK != eMin)
            {
                continue;
            }
            // Let's see which algorithm we ended up using.
            checker_.checkString(eMinimizeAlgorithmToString(simConfig.minAlg()), "algorithm");
            rmsd = mh.coordinateRmsd(mp.atomsConst(), coords, &xmin);
            checker_.checkReal(rmsd, "Coordinate RMSD after minimizing");
            add_energies(pd, &checker_, eAfter, "after");
            
            // Verify that the energy has gone down, not up.
            EXPECT_TRUE(eAfter[InteractionType::EPOT] <= eBefore[InteractionType::EPOT]);
            
            if (nma && eMinimizeStatus::OK == eMin)
            {
                // First, test calculation of the Hessian
                std::vector<double> forceZero;
                std::map<InteractionType, double> energyZero;
                std::vector<int> atomIndex;
                auto &atoms = mp.atomsConst();
                for(size_t atom = 0; atom < atoms.size(); atom++)
                {
                    if (atoms[atom].pType() == ActParticle::Atom)
                    {
                        atomIndex.push_back(atom);
                    }
                }
                const int     matrixSide = DIM*atomIndex.size();
                {
                    MatrixWrapper hessian(matrixSide, matrixSide);
                    mh.computeHessian(pd, mp.topology(), forceComp, &xmin, atomIndex,
                                      &hessian, &forceZero, &energyZero);
                    double msf2 = 0;
                    for(const auto f : forceZero)
                    {
                        msf2 += f*f;
                    }
                    checker_.checkReal(std::sqrt(msf2/forceZero.size()), "RMS force");
                    // In single precision the forces at the minimum are the
                    // rounding errors of the forces, after relaxing the
                    // shells to within them, of about 0.1 kJ/mol/nm
                    double feq = gmx::test::precisionTolerance(5e-2, 0.2);
                    checker_.setDefaultTolerance(gmx::test::FloatingPointTolerance(feq, feq, 5e-2, 5e-2,
                                                                                    UINT64_MAX, UINT64_MAX, false));
                    checker_.checkSequence(forceZero.begin(), forceZero.end(), "Equilibrium force");
                    checker_.setDefaultTolerance(tolerance);
                    // Now test the solver used in minimization
                    std::vector<double> deltaX(DIM*atomIndex.size(), 0.0);
                    int result = hessian.solve(forceZero, &deltaX);
                    EXPECT_TRUE(0 == result);
                    // checker_.checkSequence(deltaX.begin(), deltaX.end(), "DeltaX");
                    
                }
                std::vector<double> freq, freq_extern, inten, inten_extern;
                mh.nma(pd, &mp, forceComp, &xmin, &freq, &inten, nullptr);
                auto mpo = MolPropObservable::FREQUENCY;
                const char *unit = mpo_unit2(mpo);
                for(auto f = freq.begin(); f < freq.end(); ++f)
                {
                    freq_extern.push_back(convertFromGromacs(*f, unit));
                }
                auto mpoi = MolPropObservable::INTENSITY;
                const char *uniti = mpo_unit2(mpoi);
                for(auto f = inten.begin(); f < inten.end(); ++f)
                {
                    inten_extern.push_back(convertFromGromacs(*f, uniti));
                }
                checker_.checkSequence(freq_extern.begin(), freq_extern.end(), "Frequencies");
                checker_.checkSequence(inten_extern.begin(), inten_extern.end(), "Intensities");
                
                double scale_factor = 1;
                AtomizationEnergy atomenergy;
                atomenergy.read();
                const ThermoChemistry tc(nullptr, &mp, coords, atomenergy, freq, eAfter[InteractionType::EPOT],
                                         298.15, 1, scale_factor);
                checker_.checkReal(tc.ZPE(),  "Zero point energy (kJ/mol)");
                checker_.checkReal(tc.DHform(), "Delta H form (kJ/mol)");
                for(const auto &tcc : tccmap())
                {
                    checker_.checkReal(tc.S0(tcc.first), gmx::formatString("Standard entropy - %11s  (J/mol K)",
                                                                           tcc.second.c_str()).c_str());
                }
                for(const auto &tcc : tccmap())
                {
                    checker_.checkReal(tc.cv(tcc.first), gmx::formatString("Heat capacity cV - %11s (J/mol K)", 
                                                                           tcc.second.c_str()).c_str());
                }
                for(const auto &tcc : tccmap())
                {
                    checker_.checkReal(tc.Einternal(tcc.first), gmx::formatString("Internal energy  - %11s (kJ/mol)",
                                                                                  tcc.second.c_str()).c_str());
                }
            }
        }
    }
};

/*! \brief A molecule with its topology and charges
 * \param[in]  pd       The force field
 * \param[in]  file     The structure in mols
 * \param[in]  fc       The force computer for the charges
 * \param[out] mol      The molecule
 * \param[out] x        Its coordinates
 */
static void prepareMolecule(ForceField *pd, const char *file, const ForceComputer &fc,
                            ACTMol *mol, std::vector<gmx::RVec> *x)
{
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    std::vector<alexandria::MolProp> molprops;
    double qtot = 0;
    importFile(&msghandler, pd, gmx::test::TestFileManager::getInputFilePath(std::string("molecules/") + file),
               &molprops, "minimum", JobType::OPT, false, &qtot, true);
    ASSERT_TRUE(msghandler.ok());
    ASSERT_EQ(1u, molprops.size());
    mol->Merge(&molprops[0]);
    mol->GenerateTopology(&msghandler, pd, missingParameters::Ignore);
    ASSERT_TRUE(msghandler.ok());
    *x = mol->xOriginal();
    std::vector<gmx::RVec> forces(x->size());
    mol->generateCharges(&msghandler, pd, &fc, pd->chargeGenerationAlgorithm(), x, &forces);
    ASSERT_TRUE(msghandler.ok());
}

// Test: All minimizers reach the same minimum of acetone with a tolerance on
// the mean square force that they can reach
TEST_F (MolHandlerTest, MinimizersAgree)
{
    auto pd = getForceField("ACS-g1");
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    const ForceComputer    forceComp(1e-10, 100, 1);
    ACTMol                 mol;
    std::vector<gmx::RVec> x0;
    prepareMolecule(pd, "acetone.sdf", forceComp, &mol, &x0);
    const MolHandler    mh;
    std::vector<double> epot;
    for (auto alg : { eMinimizeAlgorithm::LBFGS, eMinimizeAlgorithm::Newton, eMinimizeAlgorithm::Steep })
    {
        SCOPED_TRACE(eMinimizeAlgorithmToString(alg));
        SimulationConfigHandler simConfig;
        simConfig.setMinimizeAlgorithm(alg);
        simConfig.setForceTolerance(1e-6);
        simConfig.setMaxIter(20000);
        std::vector<gmx::RVec>            x = x0;
        std::map<InteractionType, double> energies;
        EXPECT_EQ(eMinimizeStatus::OK, mh.minimizeCoordinates(&msghandler, pd, mol.topology(), &forceComp,
                                                              simConfig, &x, &energies, {}));
        epot.push_back(energies[InteractionType::EPOT]);
        EXPECT_LT(epot.back(), epot.front() + 1e-3);
    }
    EXPECT_NEAR(epot[0], epot[1], 1e-3);
    EXPECT_NEAR(epot[1], epot[2], 1e-3);
}

/* Test: Newton and steepest descent converge for a water dimer with shells
 * to a tolerance on the mean square force that single precision allows.
 * This needs a step of the finite difference Hessian that is not too small
 * for the precision, and energy changes from the forces close to the minimum.
 */
TEST_F (MolHandlerTest, MinimizersConvergeWithShells)
{
    auto pd = getForceField("ACS-pg1");
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    // The shells converge within the precision of the forces
    const ForceComputer    forceComp(gmx::test::precisionTolerance(1e-10, 1e-6), 100, 1);
    ACTMol                 mol;
    std::vector<gmx::RVec> x0;
    prepareMolecule(pd, "water_dimer.pdb", forceComp, &mol, &x0);
    const MolHandler mh;
    for (auto alg : { eMinimizeAlgorithm::Newton, eMinimizeAlgorithm::Steep })
    {
        SCOPED_TRACE(eMinimizeAlgorithmToString(alg));
        SimulationConfigHandler simConfig;
        simConfig.setMinimizeAlgorithm(alg);
        simConfig.setForceTolerance(1e-4);
        simConfig.setMaxIter(20000);
        std::vector<gmx::RVec>            x = x0;
        std::map<InteractionType, double> energies;
        EXPECT_EQ(eMinimizeStatus::OK, mh.minimizeCoordinates(&msghandler, pd, mol.topology(), &forceComp,
                                                              simConfig, &x, &energies, {}));
    }
}

/*! \brief The vibrations at a minimum are 3N - 6, or 3N - 5 for a linear
 * molecule, without the translations and rotations, independent of the
 * orientation and position of the molecule
 */
static void testNormalModes(const char *file, size_t nvibration)
{
    SCOPED_TRACE(file);
    auto pd = getForceField("ACS-g1");
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    const ForceComputer    forceComp(1e-12, 100, 1);
    ACTMol                 mol;
    std::vector<gmx::RVec> x;
    prepareMolecule(pd, file, forceComp, &mol, &x);
    const MolHandler        mh;
    SimulationConfigHandler simConfig;
    simConfig.setForceTolerance(1e-10);
    std::map<InteractionType, double> energies;
    ASSERT_EQ(eMinimizeStatus::OK, mh.minimizeCoordinates(&msghandler, pd, mol.topology(), &forceComp,
                                                          simConfig, &x, &energies, {}));
    std::vector<double> freq, intens;
    mh.nma(pd, &mol, &forceComp, &x, &freq, &intens, nullptr, false);
    ASSERT_EQ(nvibration, freq.size());
    for (auto f : freq)
    {
        EXPECT_GT(f, 1) << "a translation or rotation";
    }
    // Rotated and moved
    double c = std::cos(0.7), sn = std::sin(0.7);
    const auto r   = [](double v) { return static_cast<real>(v); };
    matrix rot = { { r(c), r(-sn), 0 }, { r(sn*0.6), r(c*0.6), r(0.8) }, { r(-sn*0.8), r(-c*0.8), r(0.6) } };
    std::vector<gmx::RVec> xr(x.size());
    for (size_t i = 0; i < x.size(); i++)
    {
        mvmul(rot, x[i], xr[i]);
        rvec shift = { 0.3, -1.2, 2.5 };
        rvec_inc(xr[i], shift);
    }
    std::vector<double> freqr, intensr;
    mh.nma(pd, &mol, &forceComp, &xr, &freqr, &intensr, nullptr, false);
    ASSERT_EQ(freq.size(), freqr.size());
    for (size_t k = 0; k < freq.size(); k++)
    {
        EXPECT_NEAR(freq[k], freqr[k], 1e-3*freq[k]) << "mode " << k;
    }
}

// Test: The normal modes of acetone and of linear carbon dioxide
TEST_F (MolHandlerTest, NormalModesWithoutRigidBody)
{
    testNormalModes("acetone.sdf", 24);
    testNormalModes("carbon-dioxide.sdf", 4);
}

TEST_F (MolHandlerTest, MethaneThiolNoFreq)
{
    test("methanethiol.sdf", "ACS-g1", false, 1e-13);
}

TEST_F (MolHandlerTest, CarbonDioxideNoFreq)
{
    test("carbon-dioxide.sdf", "ACS-g1", false, 1e-6);
}

TEST_F (MolHandlerTest, HydrogenChlorideNoFreq)
{
    test("hydrogen-chloride.sdf", "ACS-g1", false);
}

TEST_F (MolHandlerTest, WaterNoFreq)
{
    test("water-3-oep.log.pdb", "ACS-g1", false);
}

TEST_F (MolHandlerTest, AcetoneNoFreq)
{
    test("acetone-3-oep.log.pdb", "ACS-g1", false);
}

TEST_F (MolHandlerTest, UracilNoFreq)
{
    test("uracil.sdf", "ACS-g1", false, 1e-14);
}

TEST_F (MolHandlerTest, CarbonDioxideNoFreqPol)
{
    test("carbon-dioxide.sdf", "ACS-pg1", false, 1e-8);
}

TEST_F (MolHandlerTest, HydrogenChlorideNoFreqPol)
{

    test("hydrogen-chloride.sdf", "ACS-pg1", false, 1e-8);
}

TEST_F (MolHandlerTest, WaterNoFreqPol)
{

    test("water-3-oep.log.pdb", "ACS-pg1", false, 1e-8);
}

TEST_F (MolHandlerTest, AcetoneNoFreqPol)
{
    test("acetone-3-oep.log.pdb", "ACS-pg1", false, 1e-12);
}

TEST_F (MolHandlerTest, UracilNoFreqPol)
{
    // In single precision the minimization of polarizable uracil reaches
    // forces within the rounding of the energies on a plateau, from where
    // it either finds the minimum or another one, depending on that
    // rounding. The result is not reproducible, see the integer
    // coordinates in coordinate.h.
    if (!GMX_DOUBLE)
    {
        return;
    }
    test("uracil.sdf", "ACS-pg1", false, 1e-12);
}

TEST_F (MolHandlerTest, CarbonDioxide)
{
    test("carbon-dioxide.sdf", "ACS-g1", true, 1e-8);
}

TEST_F (MolHandlerTest, HydrogenChloride)
{

    test("hydrogen-chloride.sdf", "ACS-g1", true);
}

TEST_F (MolHandlerTest, Water)
{

    test("water-3-oep.log.pdb", "ACS-g1", true);
}

TEST_F (MolHandlerTest, Acetone)
{
    test("acetone-3-oep.log.pdb", "ACS-g1", true);
}

TEST_F (MolHandlerTest, Uracil)
{
    test("uracil.sdf", "ACS-g1", true, 1e-14);
}

TEST_F (MolHandlerTest, CarbonDioxidePol)
{
    test("carbon-dioxide.sdf", "ACS-pg1", true, 1e-10);
}

TEST_F (MolHandlerTest, HydrogenChloridePol)
{
    test("hydrogen-chloride.sdf", "ACS-pg1", true, 1e-12);
}

TEST_F (MolHandlerTest, WaterPol)
{
    test("water-3-oep.log.pdb", "ACS-pg1", true, 1e-12);
}

TEST_F (MolHandlerTest, AcetonePol)
{
    test("acetone-3-oep.log.pdb", "ACS-pg1", true, 1e-12);
}

TEST_F (MolHandlerTest, UracilPol)
{
    // In single precision the minimization of polarizable uracil reaches
    // forces within the rounding of the energies on a plateau, from where
    // it either finds the minimum or another one, depending on that
    // rounding. The result is not reproducible, see the integer
    // coordinates in coordinate.h.
    if (!GMX_DOUBLE)
    {
        return;
    }
    test("uracil.sdf", "ACS-pg1", true, 1e-12);
}

// Test: Without a lower minimum the structure is usable and the message
// gives its energy, without convergence in the steps it is not usable
TEST(MinimizeStatusTest, NoLowerMinimum)
{
    EXPECT_TRUE(minimumUsable(eMinimizeStatus::OK));
    EXPECT_TRUE(minimumUsable(eMinimizeStatus::NoMinimum));
    EXPECT_FALSE(minimumUsable(eMinimizeStatus::TooManySteps));
    EXPECT_FALSE(minimumUsable(eMinimizeStatus::Solver));
    EXPECT_EQ("Could not find a lower minimum than -1982.11 kJ/mol",
              minimizeStatusMessage(eMinimizeStatus::NoMinimum, -1982.11));
    EXPECT_EQ("Too many steps", minimizeStatusMessage(eMinimizeStatus::TooManySteps, -1982.11));
}

} // namespace

} // namespace alexandria
