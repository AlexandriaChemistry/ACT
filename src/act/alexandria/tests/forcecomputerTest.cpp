/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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
/*! \internal \file
 * \brief
 * Implements test of bonded force routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_act_tests
 */
#include "actpre.h"

#include "act/forces/forcecomputer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>

#include <memory>

#include <gtest/gtest.h>

#include "act/import/import.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_util.h"
#include "gromacs/utility/fatalerror.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class ForceComputerTest : public ::testing::Test
{
protected:
    gmx::test::TestReferenceData    refData_;
    gmx::test::TestReferenceChecker checker_;
 
    ForceComputerTest( ) :
        checker_(refData_.rootChecker())
    {
        const gmx::test::FloatingPointTolerance tolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-6));
        checker_.setDefaultTolerance(tolerance);
    }

    /*! \brief Set an absolute tolerance in single precision
     * The relative tolerance is the same as in double precision.
     * \param[in] absolute The absolute tolerance
     */
    void setSingleTolerance(double absolute)
    {
        if (!GMX_DOUBLE)
        {
            checker_.setDefaultTolerance(gmx::test::FloatingPointTolerance(absolute, absolute, 1e-6, 1e-6,
                                                                           UINT64_MAX, UINT64_MAX, false));
        }
    }

    //! The largest absolute energy term
    static double largestEnergy(const std::map<InteractionType, double> &energies)
    {
        double emax = 0;
        for (const auto &e : energies)
        {
            emax = std::max(emax, std::abs(e.second));
        }
        return emax;
    }

    //! The largest absolute force component on an atom that is not a shell
    static double largestForce(const std::vector<gmx::RVec> &forces,
                               const std::vector<ActAtom>   &atoms)
    {
        double fmax = 0;
        for (size_t i = 0; i < forces.size(); i++)
        {
            if (atoms[i].pType() != ActParticle::Shell)
            {
                for (int m = 0; m < DIM; m++)
                {
                    fmax = std::max(fmax, static_cast<double>(std::abs(forces[i][m])));
                }
            }
        }
        return fmax;
    }

    /*! \brief Estimate the errors due to single precision
     * The coordinates are rounded to single precision, by up to half a
     * unit in the last place. The resulting errors of the energies and
     * forces are estimated from displacements of each coordinate of the
     * atoms by one unit in the last place. The other errors are relative
     * to the largest energy term or force, since the terms cancel, like in
     * the induction or the forces on atoms with shells. In double
     * precision the estimates are zero.
     * \param[in]  msghandler  For messages
     * \param[in]  pd          The force field
     * \param[in]  fcomp       The force computer
     * \param[in]  mp          The molecule
     * \param[in]  coordinates The coordinates, with relaxed shells
     * \param[in]  energies    The energies at the coordinates
     * \param[in]  forces      The forces at the coordinates
     * \param[out] eerror      Estimated error of each energy term
     * \param[out] ferror      Estimated error of each force component
     */
    static void singleErrors(MsgHandler                              *msghandler,
                             const ForceField                        *pd,
                             ForceComputer                           *fcomp,
                             const ACTMol                            &mp,
                             const std::vector<gmx::RVec>            &coordinates,
                             const std::map<InteractionType, double> &energies,
                             const std::vector<gmx::RVec>            &forces,
                             std::map<InteractionType, double>       *eerror,
                             std::vector<gmx::DVec>                  *ferror)
    {
        ferror->assign(forces.size(), { 0, 0, 0 });
        eerror->clear();
        for (const auto &e : energies)
        {
            (*eerror)[e.first] = 0;
        }
        if (GMX_DOUBLE)
        {
            return;
        }
        const auto &atoms = mp.atomsConst();
        // The differences of the coordinates, from the integer coordinates,
        // are rounded to float, see coordinate.h, by up to the rounding of
        // the largest difference. The change of the energies and forces
        // when a coordinate moves by that much is their rounding error.
        double extent = 0;
        for (const auto &xi : coordinates)
        {
            for (int m = 0; m < DIM; m++)
            {
                extent = std::max(extent, 2*std::abs(static_cast<double>(xi[m])));
            }
        }
        double step = 0.5*GMX_FLOAT_EPS*extent;
        for (size_t i = 0; i < coordinates.size(); i++)
        {
            if (atoms[i].pType() == ActParticle::Shell)
            {
                continue;
            }
            for (int m = 0; m < DIM; m++)
            {
                auto xd = coordinates;
                xd[i][m] = std::nextafter(static_cast<real>(xd[i][m] + step), std::numeric_limits<real>::max());
                std::vector<gmx::RVec>            fd(forces.size(), { 0, 0, 0 });
                std::map<InteractionType, double> ed;
                fcomp->compute(msghandler, pd, mp.topology(), &xd, &fd, &ed);
                for (auto &e : *eerror)
                {
                    e.second += 0.5*std::abs(ed[e.first] - energies.find(e.first)->second);
                }
                for (size_t k = 0; k < forces.size(); k++)
                {
                    for (int n = 0; n < DIM; n++)
                    {
                        (*ferror)[k][n] += 0.5*std::abs(fd[k][n] - forces[k][n]);
                    }
                }
            }
        }
        double emax = largestEnergy(energies);
        for (auto &e : *eerror)
        {
            e.second += 1e-5*emax;
        }
        double fmax = largestForce(forces, atoms);
        for (auto &f : *ferror)
        {
            for (int n = 0; n < DIM; n++)
            {
                f[n] += 1e-5*fmax;
            }
        }
    }

    void test(const char *molname, const char *forcefield, 
              bool testPolarizability, double stretch = 1, int nexcl = 0)
    {
        // Get forcefield
        auto qd  = getForceField(forcefield);
        // Copy FF to prevent that changes mess up opther tests
        ForceField qqd = *qd;
        ForceField *pd = &qqd;
        if (nexcl > 2)
        {
            pd->setNexcl(nexcl);
        }
        
        double rmsToler = 1e-12;
        auto fcomp = new ForceComputer(rmsToler, 25, 1.0);
        
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);

        // The molecule
        std::vector<ACTMol> mps;
        std::string moldir("molecules/");
        moldir += molname;
        const auto dataName = gmx::test::TestFileManager::getInputFilePath(moldir);
        initACTMol(&msghandler, dataName, pd, fcomp, &mps);
   
        for(auto &mp : mps)
        {
            std::vector<gmx::RVec> coordinates, forces;
            auto xo = mp.xOriginal();
            for(size_t i = 0; i < xo.size(); i++)
            {
                rvec xxx;
                for(int m = 0; m < DIM; m++)
                {
                    xxx[m] = xo[i][m]*stretch;
                }
                coordinates.emplace_back(xxx);
                forces.emplace_back(0, 0, 0);
            }
            if (testPolarizability)
            {
                for(auto qp = mp.qProps()->begin(); qp < mp.qProps()->end(); ++qp)
                {
                    auto qCalc = qp->qPact();
                    qCalc->initializeMoments();
                    qCalc->calcPolarizability(&msghandler, pd, mp.topology(), fcomp);
                    auto alpha = qCalc->polarizabilityTensor();
                    const char *xyz[DIM] = { "X", "Y", "Z" };
                    // In single precision the errors are relative to the
                    // largest element, due to the resolution of the shell
                    // coordinates compared to their displacements, and of
                    // the forces in single precision on them.
                    double amax = 0;
                    for(int m = 0; m < DIM; m++)
                    {
                        for(int n = 0; n < DIM; n++)
                        {
                            amax = std::max(amax, std::abs(convertFromGromacs(alpha[m][n], "A^3")));
                        }
                    }
                    setSingleTolerance(2e-5*amax);
                
                    for(int m = 0; m < DIM; m++)
                    {
                        for(int n = 0; n < DIM; n++)
                        {
                            const std::string label = gmx::formatString("alpha[%s][%s]", xyz[m], xyz[n]);
                            checker_.checkReal(convertFromGromacs(alpha[m][n], "A^3"),
                                               label.c_str());
                        }
                    }
                }
            }
            else
            {
                std::map<InteractionType, double> actEnergies;
                fcomp->compute(&msghandler, pd, mp.topology(), &coordinates, &forces, &actEnergies);
                std::map<InteractionType, double> eerror;
                std::vector<gmx::DVec>            ferror;
                singleErrors(&msghandler, pd, fcomp, mp, coordinates, actEnergies, forces, &eerror, &ferror);
                for(auto &ifm : actEnergies)
                {
                    std::string label = interactionTypeToString(ifm.first);
                    if (stretch != 1)
                    {
                        label += gmx::formatString("%g", stretch);
                    }
                    const auto actlabel = label+"_act";
                    double actEner  = actEnergies[ifm.first];
                    setSingleTolerance(eerror[ifm.first]);
                    checker_.checkReal(actEner, actlabel.c_str());
                }
                const auto &atoms = mp.atomsConst();
                const char *xyz[DIM] = { "X", "Y", "Z" };
                for(size_t i = 0; i < forces.size(); i++)
                {
                    bool shell = atoms[i].pType() == ActParticle::Shell;
                    for(int m = 0; m < DIM; m++)
                    {
                        if (!shell)
                        {
                            std::string stretchName;
                            if (stretch != 1)
                            {
                                stretchName += gmx::formatString("%g", stretch);
                            }
                            setSingleTolerance(ferror[i][m]);
                            checker_.checkReal(forces[i][m], gmx::formatString("%s-%zu_act%s f%s", 
                                                                               atoms[i].ffType().c_str(),
                                                                               i+1, stretchName.c_str(),
                                                                               xyz[m]).c_str());
                        }
                    }
                }
            }
        }
    }
    
    void testSimple(const char *molname, const char *forcefield)
    {
        // Get forcefield
        auto pd  = getForceField(forcefield);
        
        double rmsToler = 0.0000001;
        auto fcomp = new ForceComputer(rmsToler, 25, 1.0);
        
        MsgHandler msghandler;
        msghandler.setPrintLevel(ACTStatus::Warning);

        // The molecule
        std::vector<ACTMol> mps;
        std::string moldir("molecules/");
        moldir += molname;
        const auto dataName = gmx::test::TestFileManager::getInputFilePath(moldir);
        initACTMol(&msghandler, dataName, pd, fcomp, &mps);
   
        for(auto &mp : mps)
        {
            std::vector<gmx::RVec> coordinates, forces;
            auto xo = mp.xOriginal();
            for(size_t i = 0; i < xo.size(); i++)
            {
                coordinates.push_back(xo[i]);
                forces.emplace_back(0, 0, 0);
            }
            std::map<InteractionType, double> actEnergies;
            const auto fsc = pd->forcesConst();
            fcomp->compute(&msghandler, pd, mp.topology(), &coordinates, &forces, &actEnergies);
            std::map<InteractionType, double> eerror;
            std::vector<gmx::DVec>            ferror;
            singleErrors(&msghandler, pd, fcomp, mp, coordinates, actEnergies, forces, &eerror, &ferror);
            for(auto &ifm : actEnergies)
            {
                auto actEner  = actEnergies[ifm.first];
                setSingleTolerance(eerror[ifm.first]);
                checker_.checkReal(actEner, interactionTypeToString(ifm.first).c_str());
                const auto &atoms = mp.atomsConst();
                const char *xyz[DIM] = { "X", "Y", "Z" };
                for(size_t i = 0; i < forces.size(); i++)
                {
                    bool shell = atoms[i].pType() == ActParticle::Shell;
                    for(int m = 0; m < DIM; m++)
                    {
                        if (!shell)
                        {
                            setSingleTolerance(ferror[i][m]);
                            checker_.checkReal(forces[i][m], gmx::formatString("%s-%zu f%s", 
                                                                               atoms[i].ffType().c_str(),
                                                                               i+1, xyz[m]).c_str());
                        }
                    }
                }
            }
        }
    }
};

TEST_F (ForceComputerTest, MethaneThiol)
{
    test("methanethiol.sdf", "ACS-g1", false);
}

TEST_F (ForceComputerTest, CarbonDioxide)
{
    test("carbon-dioxide.sdf", "ACS-g1", false);
}

TEST_F (ForceComputerTest, HydrogenChloride)
{

    test("hydrogen-chloride.sdf", "ACS-g1", false);
}

TEST_F (ForceComputerTest, HydrogenChlorideStretch)
{
    test("hydrogen-chloride.sdf", "ACS-g1", false, 0.98);
    test("hydrogen-chloride.sdf", "ACS-g1", false, 1);
    test("hydrogen-chloride.sdf", "ACS-g1", false, 1.02);
}

TEST_F (ForceComputerTest, Water)
{
    test("water-3-oep.log.pdb", "ACS-g1", false);
}

TEST_F (ForceComputerTest, WaterStretch)
{
    test("water-3-oep.log.pdb", "ACS-g1", false, 0.98);
    test("water-3-oep.log.pdb", "ACS-g1", false, 1.02);
}

TEST_F (ForceComputerTest, Acetone)
{
    test("acetone-3-oep.log.pdb", "ACS-g1", false);
}

TEST_F (ForceComputerTest, AcetoneSdf)
{
    test("acetone.sdf", "ACS-g1", false);
}

TEST_F (ForceComputerTest, AcetoneStretch)
{
    test("acetone-3-oep.log.pdb", "ACS-g1", false, 0.98);
    test("acetone-3-oep.log.pdb", "ACS-g1", false, 1.07);
}

TEST_F (ForceComputerTest, AcetoneNonPlanar)
{
    test("acetone-nonplanar.pdb", "ACS-g1", false);
}

TEST_F (ForceComputerTest, Uracil)
{
    test("uracil.sdf", "ACS-g1", false);
}

TEST_F (ForceComputerTest, UracilStretch)
{
    test("uracil.sdf", "ACS-g1", false, 0.9);
    test("uracil.sdf", "ACS-g1", false, 1.1);
}

TEST_F (ForceComputerTest, AcetoneDih)
{
    test("acetone-3-oep.log.pdb", "ACS-g-dih", false);
}

TEST_F (ForceComputerTest, UracilDih)
{
    test("uracil.sdf", "ACS-g-dih", false);
}

TEST_F (ForceComputerTest, WaterUB)
{

    test("water-3-oep.log.pdb", "ACS-g-dih", false);
}

TEST_F (ForceComputerTest, HydrogenChloridePolStretch)
{
    test("hydrogen-chloride.sdf", "ACS-pg1", false, 0.98);
    test("hydrogen-chloride.sdf", "ACS-pg1", false, 1.05);
}

TEST_F (ForceComputerTest, WaterPol)
{
    test("water-3-oep.log.pdb", "ACS-pg1", false);
}

TEST_F (ForceComputerTest, WaterPolStretch)
{
    test("water-3-oep.log.pdb", "ACS-pg1", false, 0.98);
    test("water-3-oep.log.pdb", "ACS-pg1", false, 1.04);
}

TEST_F (ForceComputerTest, AcetonePol)
{
    test("acetone-3-oep.log.pdb", "ACS-pg1", false);
}

TEST_F (ForceComputerTest, AcetonePolStretch)
{
    test("acetone-3-oep.log.pdb", "ACS-pg1", false, 0.97);
    test("acetone-3-oep.log.pdb", "ACS-pg1", false, 1.01);
}

TEST_F (ForceComputerTest, AcetoneNonPlanarPol)
{
    test("acetone-nonplanar.pdb", "ACS-pg1", false);
}

TEST_F (ForceComputerTest, UracilPol)
{
    test("uracil.sdf", "ACS-pg1", false);
}

TEST_F (ForceComputerTest, UracilPolStretch)
{
    test("uracil.sdf", "ACS-pg1", false, 0.99);
    test("uracil.sdf", "ACS-pg1", false, 1.06);
}

TEST_F (ForceComputerTest, CarbonDioxidePolarizability)
{
  test("carbon-dioxide.sdf", "ACS-pg1", true);
}

TEST_F (ForceComputerTest, HydrogenChloridePolarizability)
{
    test("hydrogen-chloride.sdf", "ACS-pg1", true);
}

TEST_F (ForceComputerTest, WaterPolarizability)
{
    test("water-3-oep.log.pdb", "ACS-pg1", true);
}

TEST_F (ForceComputerTest, AcetonePolarizability)
{
    test("acetone-3-oep.log.pdb", "ACS-pg1", true);
}

TEST_F (ForceComputerTest, AcetonePolGbham)
{
    testSimple("acetone-3-oep.log.pdb", "ACS-pg-gbham");
}

TEST_F (ForceComputerTest, AcetonePolarizabilityGbham)
{
    test("acetone-3-oep.log.pdb", "ACS-pg-gbham", true);
}

TEST_F (ForceComputerTest, UracilPolarizability)
{
    test("uracil.sdf", "ACS-pg1", true);
}

// These have to be last since they modify the force field that is shared,
TEST_F (ForceComputerTest, MethaneThiolNexcl3)
{
    test("methanethiol.sdf", "ACS-g1", false, 1, 3);
}

TEST_F (ForceComputerTest, UracilPolNexcl3)
{
    test("uracil.sdf", "ACS-pg1", false, 1, 3);
}

//! \return the energies of the uracil with nexcl 3, optionally without proper dihedrals
static std::map<InteractionType, double> uracilEnergies(bool propers)
{
    ForceField pd = *getForceField("ACS-pg1");
    for (auto &fs : *pd.forces())
    {
        if (fs.second.optionExists("nexcl"))
        {
            fs.second.removeOption("nexcl");
            fs.second.addOption("nexcl", "3");
        }
    }
    if (!propers)
    {
        pd.findForces(InteractionType::PROPER_DIHEDRALS)->parameters()->clear();
    }
    ForceComputer fcomp(1e-12, 25, 1.0);
    MsgHandler    msghandler;
    msghandler.setPrintLevel(ACTStatus::Warning);
    std::vector<ACTMol> mps;
    initACTMol(&msghandler, gmx::test::TestFileManager::getInputFilePath("molecules/uracil.sdf"),
               &pd, &fcomp, &mps);
    EXPECT_EQ(1u, mps.size());
    auto                              xo = mps[0].xOriginal();
    std::vector<gmx::RVec>            coordinates(xo.begin(), xo.end());
    std::vector<gmx::RVec>            forces(xo.size(), { 0, 0, 0 });
    std::map<InteractionType, double> energies;
    fcomp.compute(&msghandler, &pd, mps[0].topology(), &coordinates, &forces, &energies);
    return energies;
}

// Test: The exclusions follow from the bonds, not from the proper
// dihedrals in the force field
TEST (ForceComputerExclusionTest, Nexcl3WithoutPropers)
{
    auto with    = uracilEnergies(true);
    auto without = uracilEnergies(false);
    for (auto itype : { InteractionType::DISPERSION, InteractionType::EXCHANGE,
                        InteractionType::ELECTROSTATICS, InteractionType::INDUCTION })
    {
        ASSERT_TRUE(with.find(itype) != with.end());
        ASSERT_TRUE(without.find(itype) != without.end());
        EXPECT_NEAR(with[itype], without[itype], 1e-8*std::max(1.0, std::abs(with[itype])))
            << interactionTypeToString(itype);
    }
}

TEST_F (ForceComputerTest, AcetonePolarizabilityGbhamNexcl3)
{
    test("acetone-3-oep.log.pdb", "ACS-pg-gbham", true, 1, 3);
}


}  // namespace

}  // namespace gmx
