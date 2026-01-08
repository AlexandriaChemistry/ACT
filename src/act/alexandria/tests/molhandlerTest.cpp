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

#include "act/import/atype_mapping.h"
#include "act/import/import.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/actmol.h"
#include "act/basics/msg_handler.h"
#include "act/alexandria/thermochemistry.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/qgen/qgen_acm.h"
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
    auto fsc = pd->forcesConst();
    for(auto &imf : energies)
    {
        std::string fname   = interactionTypeToString(imf.first);
        std::string mylabel = gmx::formatString("%s %s", fname.c_str(), label);

        checker->checkReal(imf.second, mylabel.c_str());
    }
}

class MolHandlerTest : public gmx::test::CommandLineTestBase
{
protected:
    void test(const char *molname, const char *forcefield, bool nma, double ftoler=1e-15)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        int                              maxpot   = 100;
        int                              nsymm    = 0;
        const char                      *conf     = (char *)"minimum";
        std::string                      method, basis;
        const char                      *jobtype  = (char *)"Opt";
        
        std::string                      dataName;
        std::vector<alexandria::MolProp> molprops;
        
        // Get forcefield
        auto pd         = getForceField(forcefield);
        std::string moldir("mols/");
        moldir += molname;
        dataName        = gmx::test::TestFileManager::getInputFilePath(moldir);
        double qtot     = 0;
        bool   userqtot = false;
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);
        matrix box;
        importFile(&msghandler, pd, dataName.c_str(), &molprops,
                   molname, molname,
                   conf, &method, &basis,
                   maxpot, nsymm, jobtype, userqtot, &qtot, false, box, true);
        EXPECT_TRUE(msghandler.ok());
        std::vector<ACTMol> mps;
        // Needed for GenerateCharges
        auto alg = pd->chargeGenerationAlgorithm();
        double shellTolerance = ftoler;
        int    shellMaxIter   = 100;
        auto forceComp = new ForceComputer(shellTolerance, shellMaxIter);
        std::vector<double>    qcustom;
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
            
            MolHandler mh;
            
            std::vector<gmx::RVec> xmin = coords;
            double rmsd = mh.coordinateRmsd(&mp, coords, &xmin);
            checker_.checkReal(rmsd, "Coordinate RMSD before minimizing");
            // Infinite number of shell iterations, i.e. until convergence.
            std::map<InteractionType, double> eAfter;
            SimulationConfigHandler simConfig;
            simConfig.setForceTolerance(ftoler);
            // Change nullptr to stdout for debugging.
            auto eMin = mh.minimizeCoordinates(&msghandler, pd, &mp, forceComp,
                                               simConfig, &xmin, &eAfter, {});
            EXPECT_TRUE(eMinimizeStatus::OK == eMin);
            if (eMinimizeStatus::OK != eMin)
            {
                continue;
            }
            // Let's see which algorithm we ended up using.
            checker_.checkString(eMinimizeAlgorithmToString(simConfig.minAlg()), "algorithm");
            rmsd = mh.coordinateRmsd(&mp, coords, &xmin);
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
                    mh.computeHessian(pd, &mp, forceComp, &xmin, atomIndex,
                                      &hessian, &forceZero, &energyZero);
                    double msf2 = 0;
                    for(const auto f : forceZero)
                    {
                        msf2 += f*f;
                    }
                    checker_.checkReal(std::sqrt(msf2/forceZero.size()), "RMS force");
                    checker_.checkSequence(forceZero.begin(), forceZero.end(), "Equilibrium force");
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
                ThermoChemistry tc(&mp, coords, atomenergy, freq, eAfter[InteractionType::EPOT],
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

TEST_F (MolHandlerTest, MethaneThiolNoFreq)
{
    test("methanethiol.sdf", "ACS-g", false, 1e-13);
}

TEST_F (MolHandlerTest, CarbonDioxideNoFreq)
{
    test("carbon-dioxide.sdf", "ACS-g", false, 1e-6);
}

TEST_F (MolHandlerTest, HydrogenChlorideNoFreq)
{
    test("hydrogen-chloride.sdf", "ACS-g", false);
}

TEST_F (MolHandlerTest, WaterNoFreq)
{
    test("water-3-oep.log.pdb", "ACS-g", false);
}

TEST_F (MolHandlerTest, AcetoneNoFreq)
{
    test("acetone-3-oep.log.pdb", "ACS-g", false);
}

TEST_F (MolHandlerTest, UracilNoFreq)
{
    test("uracil.sdf", "ACS-g", false, 1e-14);
}

TEST_F (MolHandlerTest, CarbonDioxideNoFreqPol)
{
    test("carbon-dioxide.sdf", "ACS-pg", false, 1e-8);
}

TEST_F (MolHandlerTest, HydrogenChlorideNoFreqPol)
{

    test("hydrogen-chloride.sdf", "ACS-pg", false, 1e-8);
}

TEST_F (MolHandlerTest, WaterNoFreqPol)
{

    test("water-3-oep.log.pdb", "ACS-pg", false, 1e-8);
}

TEST_F (MolHandlerTest, AcetoneNoFreqPol)
{
    test("acetone-3-oep.log.pdb", "ACS-pg", false, 1e-12);
}

TEST_F (MolHandlerTest, UracilNoFreqPol)
{
    test("uracil.sdf", "ACS-pg", false, 1e-12);
}

TEST_F (MolHandlerTest, CarbonDioxide)
{
    test("carbon-dioxide.sdf", "ACS-g", true, 1e-8);
}

TEST_F (MolHandlerTest, HydrogenChloride)
{

    test("hydrogen-chloride.sdf", "ACS-g", true);
}

TEST_F (MolHandlerTest, Water)
{

    test("water-3-oep.log.pdb", "ACS-g", true);
}

TEST_F (MolHandlerTest, Acetone)
{
    test("acetone-3-oep.log.pdb", "ACS-g", true);
}

TEST_F (MolHandlerTest, Uracil)
{
    test("uracil.sdf", "ACS-g", true, 1e-14);
}

TEST_F (MolHandlerTest, CarbonDioxidePol)
{
    test("carbon-dioxide.sdf", "ACS-pg", true, 1e-10);
}

TEST_F (MolHandlerTest, HydrogenChloridePol)
{
    test("hydrogen-chloride.sdf", "ACS-pg", true, 1e-12);
}

TEST_F (MolHandlerTest, WaterPol)
{
    test("water-3-oep.log.pdb", "ACS-pg", true, 1e-12);
}

TEST_F (MolHandlerTest, AcetonePol)
{
    test("acetone-3-oep.log.pdb", "ACS-pg", true, 1e-12);
}

TEST_F (MolHandlerTest, UracilPol)
{
    test("uracil.sdf", "ACS-pg", true, 1e-12);
}

} // namespace

} // namespace alexandria
