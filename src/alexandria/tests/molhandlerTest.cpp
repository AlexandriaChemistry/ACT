/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022
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

#include <math.h>

#include <map>

#include <gtest/gtest.h>

#include "act/poldata/poldata.h"
#include "act/poldata/poldata_utils.h"
#include "act/poldata/poldata_xml.h"
#include "act/qgen/qgen_acm.h"
#include "act/utility/units.h"
#include "alexandria/atype_mapping.h"
#include "alexandria/babel_io.h"
#include "alexandria/fill_inputrec.h"
#include "alexandria/molhandler.h"
#include "alexandria/mymol.h"
#include "alexandria/thermochemistry.h"
#include "gromacs/linearalgebra/eigensolver.h"
#include "gromacs/utility/fatalerror.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

static void add_energies(const Poldata                           *pd,
                         gmx::test::TestReferenceChecker         *checker,
                         const std::map<InteractionType, double> &energies,
                         const char                              *label)
{
    std::map<InteractionType, int> i2f = {
        { InteractionType::EPOT,       F_EPOT },
        { InteractionType::REPULSION,  F_REPULSION },
        { InteractionType::DISPERSION, F_DISPERSION }
    };
    auto fsc = pd->forcesConst();
    for(auto &imf : energies)
    {
        int ftype;
        if (i2f.find(imf.first) != i2f.end())
        {
            ftype = i2f.find(imf.first)->second;
        }
        else
        {
            auto i = fsc.find(imf.first);
            EXPECT_TRUE(fsc.end() != i);
            ftype = i->second.fType();
        }
        std::string mylabel = gmx::formatString("%s %s",
                                                interaction_function[ftype].longname, label);
        
        checker->checkReal(imf.second, mylabel.c_str());
    }
}

class MolHandlerTest : public gmx::test::CommandLineTestBase
{
protected:
    void test(const char *molname, const char *forcefield, bool nma)
    {
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        int                             maxpot   = 100;
        int                             nsymm    = 0;
        const char                     *conf     = (char *)"minimum";
        std::string                     method, basis;
        const char                     *jobtype  = (char *)"Opt";
        
        std::string                     dataName;
        auto molprop = new alexandria::MolProp;
        
        dataName = gmx::test::TestFileManager::getInputFilePath(molname);
        double qtot = 0;
        bool readOK = readBabel(dataName.c_str(), molprop, molname, molname,
                                conf, &method, &basis,
                                maxpot, nsymm, jobtype, &qtot, false);
        EXPECT_TRUE(readOK);
        if (readOK)
        {
            std::map<std::string, std::string> g2a;
            gaffToAlexandria("", &g2a);
            if (!g2a.empty())
            {
                EXPECT_TRUE(renameAtomTypes(molprop, g2a));
            }
        }
        MyMol mp_;
        mp_.Merge(molprop);
        // Generate charges and topology
        t_inputrec      inputrecInstance;
        t_inputrec     *inputrec   = &inputrecInstance;
        fill_inputrec(inputrec);
        mp_.setInputrec(inputrec);
        
        // Get poldata
        auto pd  = getPoldata(forcefield);
        auto imm = mp_.GenerateTopology(stdout, pd,
                                        missingParameters::Error, false);
        EXPECT_TRUE(immStatus::OK == imm);
        if (immStatus::OK != imm)
        {
            fprintf(stderr, "Could not generate topology because '%s'. Used basis %s and method %s.\n",
                    immsg(imm), basis.c_str(), method.c_str());
            return;
        }
        // Needed for GenerateCharges
        CommunicationRecord cr;
        gmx::MDLogger  mdlog {};
        auto alg = ChargeGenerationAlgorithm::NONE;
        double shellTolerance = 1e-12;
        int    shellMaxIter   = 100;
        auto forceComp = new ForceComputer(shellTolerance, shellMaxIter);
        std::vector<double>    qcustom;
        bool                   qSymm = false;
        std::vector<gmx::RVec> forces(mp_.atomsConst().size());
        std::vector<gmx::RVec> coords = mp_.x();
        mp_.symmetrizeCharges(pd, qSymm, nullptr);
        mp_.GenerateCharges(pd, forceComp, mdlog, &cr, alg, qcustom, &coords, &forces);
        
        std::map<InteractionType, double> eBefore;
        (void) forceComp->compute(pd, mp_.topology(), &coords, &forces, &eBefore);
        add_energies(pd, &checker_, eBefore, "before");
        
        MolHandler mh;
        
        std::vector<gmx::RVec> xmin = coords;
        double rmsd = mh.coordinateRmsd(&mp_, coords, &xmin);
        checker_.checkReal(rmsd, "Coordinate RMSD before minimizing");
        // Infinite number of shell iterations, i.e. until convergence.
        std::map<InteractionType, double> eAfter;
        SimulationConfigHandler simConfig;
        simConfig.setForceTolerance(1e-4);
        auto eMin = mh.minimizeCoordinates(pd, &mp_, forceComp, simConfig,
                                           &xmin, &eAfter, nullptr);
        if (eMinimizeStatus::OK != eMin)
        {
            // New try using steepest descents
            simConfig.setMinimizeAlgorithm(eMinimizeAlgorithm::Steep);
            xmin    = coords;
            simConfig.setMaxIter(5000);
            eMin    = mh.minimizeCoordinates(pd, &mp_, forceComp, simConfig,
                                             &xmin, &eAfter, nullptr);
        }
        EXPECT_TRUE(eMinimizeStatus::OK == eMin);
        // Let's see which algorithm we ended up using.
        checker_.checkString(eMinimizeAlgorithmToString(simConfig.minAlg()), "algorithm");
        rmsd = mh.coordinateRmsd(&mp_, coords, &xmin);
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
            auto &atoms = mp_.atomsConst();
            for(size_t atom = 0; atom < atoms.size(); atom++)
            {
                if (atoms[atom].pType() == eptAtom)
                {
                    atomIndex.push_back(atom);
                }
            }
            const int     matrixSide = DIM*atomIndex.size();
            {
                MatrixWrapper hessian(matrixSide, matrixSide);
                mh.computeHessian(pd, &mp_, forceComp, &xmin, atomIndex,
                                  &hessian, &forceZero, &energyZero);
                checker_.checkSequence(forceZero.begin(), forceZero.end(), "Equilibrium force");
                // Now test the solver used in minimization
                std::vector<double> deltaX(DIM*atomIndex.size(), 0.0);
                int result = hessian.solve(forceZero, &deltaX);
                EXPECT_TRUE(0 == result);
                checker_.checkSequence(deltaX.begin(), deltaX.end(), "DeltaX");
                
            }
            std::vector<double> freq, freq_extern, inten, inten_extern;
            mh.nma(pd, &mp_, forceComp, &xmin, &freq, &inten, nullptr);
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
            ThermoChemistry tc(&mp_, atomenergy, freq, 298.15, 1, scale_factor);
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
};

TEST_F (MolHandlerTest, MethaneThiolNoFreq)
{
    test("methanethiol.sdf", "ACS-g", false);
}

TEST_F (MolHandlerTest, CarbonDioxideNoFreq)
{
    test("carbon-dioxide.sdf", "ACS-g", false);
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
    test("uracil.sdf", "ACS-g", false);
}

TEST_F (MolHandlerTest, CarbonDioxideNoFreqPol)
{
    test("carbon-dioxide.sdf", "ACS-pg", false);
}

TEST_F (MolHandlerTest, HydrogenChlorideNoFreqPol)
{

    test("hydrogen-chloride.sdf", "ACS-pg", false);
}

TEST_F (MolHandlerTest, WaterNoFreqPol)
{

    test("water-3-oep.log.pdb", "ACS-pg", false);
}

TEST_F (MolHandlerTest, AcetoneNoFreqPol)
{
    test("acetone-3-oep.log.pdb", "ACS-pg", false);
}

TEST_F (MolHandlerTest, UracilNoFreqPol)
{
    test("uracil.sdf", "ACS-pg", false);
}

// We cannot run these tests in debug mode because the LAPACK library
// performs a 1/0 calculation to test the exception handling.
// #if CMAKE_BUILD_TYPE != CMAKE_BUILD_TYPE_DEBUG
TEST_F (MolHandlerTest, CarbonDioxide)
{
    test("carbon-dioxide.sdf", "ACS-g", true);
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
    test("uracil.sdf", "ACS-g", true);
}

TEST_F (MolHandlerTest, CarbonDioxidePol)
{
    test("carbon-dioxide.sdf", "ACS-pg", true);
}

TEST_F (MolHandlerTest, HydrogenChloridePol)
{

    test("hydrogen-chloride.sdf", "ACS-pg", true);
}

TEST_F (MolHandlerTest, WaterPol)
{

    test("water-3-oep.log.pdb", "ACS-pg", true);
}

TEST_F (MolHandlerTest, AcetonePol)
{
    test("acetone-3-oep.log.pdb", "ACS-pg", true);
}

TEST_F (MolHandlerTest, UracilPol)
{
    test("uracil.sdf", "ACS-pg", true);
}
// #endif

} // namespace

} // namespace alexandria
