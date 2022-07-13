/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

static void add_energies(gmx::test::TestReferenceChecker *checker,
                         const real                       ener[F_NRE],
                         const char                      *label)
{
    for(int i = 0; i < F_NRE; i++)
    {
        real ee = ener[i];
        if (ee != 0)
        {
            std::string mylabel = gmx::formatString("%s %s",
                                                    interaction_function[i].longname, label);
            checker->checkReal(ee, mylabel.c_str());
        }
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
        auto           pnc      = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
        gmx::MDLogger  mdlog {};
        auto alg = ChargeGenerationAlgorithm::NONE;
        auto forceComp = new ForceComputer(pd);
        std::vector<double> qcustom;
        bool qSymm = false;
        mp_.symmetrizeCharges(pd, qSymm, nullptr);
        mp_.GenerateCharges(pd, forceComp, mdlog, &cr, alg, qcustom);
        
        // real shellForceRMS;
        // (void) mp_.calculateEnergy(cr.commrec(), &shellForceRMS);
        mp_.calculateEnergy(forceComp);
        add_energies(&checker_, mp_.energyTerms(), "before");

        MolHandler mh;
        
        std::map<coordSet, std::vector<gmx::RVec> > xrmsd; 
        double rmsd = mh.coordinateRmsd(&mp_, &xrmsd);
        checker_.checkReal(rmsd, "Coordinate RMSD before minimizing");
        double overRelax  = 1;
        // MS force tolerance
        double forceToler = 1e-8;
        // Infinite number of shell iterations, i.e. until convergence.
        int    maxIter    = 0;
        (void) mh.minimizeCoordinates(&mp_, forceComp, nullptr, maxIter, overRelax, forceToler);

        rmsd = mh.coordinateRmsd(&mp_, &xrmsd);
        checker_.checkReal(rmsd, "Coordinate RMSD after minimizing");
        add_energies(&checker_, mp_.energyTerms(), "after");

        if (nma)
        {
            std::vector<double> freq, freq_extern, inten;
            mh.nma(&mp_, forceComp, &freq, &inten, nullptr);
            auto mpo = MolPropObservable::FREQUENCY;
            const char *unit = mpo_unit2(mpo);
            for(auto f = freq.begin(); f < freq.end(); ++f)
            {
                freq_extern.push_back(convertFromGromacs(*f, unit));
            }
            checker_.checkSequence(freq_extern.begin(), freq_extern.end(), "Frequencies");
            checker_.checkSequence(inten.begin(), inten.end(), "Intensities");
            
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

// We cannot run these tests in debug mode because the LAPACK library
// performs a 1/0 calculation to test the exception handling.
#if CMAKE_BUILD_TYPE != CMAKE_BUILD_TYPE_DEBUG
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
#endif

} // namespace

} // namespace alexandria
