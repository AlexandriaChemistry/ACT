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
 * Physical checks of complete molecules that any correct implementation
 * passes: the forces are minus the derivatives of the energy, with virtual
 * sites and relaxed shells, the total force and torque are zero, and the
 * energies, forces and charges do not change when the molecule is rotated,
 * translated or its atoms are numbered in another order. Unlike reference
 * data, these do not have to be updated when results change on purpose.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <cmath>
#include <fstream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/basics/chargemodel.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/forces/conformationdata.h"
#include "act/forces/forcecomputer.h"
#include "act/forces/forces_simd.h"
#include "act/import/import.h"
#include "act/topology/actmol.h"
#include "act/utility/timings.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/*! \brief Write an sdf file with the atoms in reverse order
 * The bonds and the formal charges refer to the new numbers. Only the
 * V2000 format of the test molecules is supported.
 * \param[in] input  The sdf file
 * \param[in] output The sdf file to write
 * \return Whether it worked
 */
bool writeReversedSdf(const std::string &input, const std::string &output)
{
    std::ifstream            in(input);
    std::vector<std::string> lines;
    std::string              line;
    while (std::getline(in, line))
    {
        lines.push_back(line);
    }
    if (lines.size() < 4 || lines[3].find("V2000") == std::string::npos)
    {
        return false;
    }
    int natom = std::stoi(lines[3].substr(0, 3));
    int nbond = std::stoi(lines[3].substr(3, 3));
    if (static_cast<int>(lines.size()) < 4 + natom + nbond)
    {
        return false;
    }
    std::ofstream out(output);
    for (int i = 0; i < 4; i++)
    {
        out << lines[i] << "\n";
    }
    for (int i = natom - 1; i >= 0; i--)
    {
        out << lines[4 + i] << "\n";
    }
    auto renumber = [natom](int a) { return natom + 1 - a; };
    for (int b = 0; b < nbond; b++)
    {
        const auto &bl = lines[4 + natom + b];
        out << gmx::formatString("%3d%3d", renumber(std::stoi(bl.substr(0, 3))),
                                 renumber(std::stoi(bl.substr(3, 3))))
            << bl.substr(6) << "\n";
    }
    for (size_t i = 4 + natom + nbond; i < lines.size(); i++)
    {
        if (lines[i].compare(0, 6, "M  CHG") == 0)
        {
            std::istringstream iss(lines[i].substr(6));
            int                n;
            iss >> n;
            out << gmx::formatString("M  CHG%3d", n);
            for (int k = 0; k < n; k++)
            {
                int a, q;
                iss >> a >> q;
                out << gmx::formatString(" %3d %3d", renumber(a), q);
            }
            out << "\n";
        }
        else
        {
            out << lines[i] << "\n";
        }
    }
    return out.good();
}

/*! \brief Make a rotation matrix
 * \param[in]  axis  The rotation axis, not necessarily normalized
 * \param[in]  angle The angle in radians
 * \param[out] R     The matrix
 */
void rotationMatrix(gmx::RVec axis, double angle, matrix R)
{
    unitv(axis, axis);
    double c = std::cos(angle);
    double s = std::sin(angle);
    double t = 1 - c;
    double x = axis[XX], y = axis[YY], z = axis[ZZ];
    R[XX][XX] = t*x*x + c;   R[XX][YY] = t*x*y - s*z; R[XX][ZZ] = t*x*z + s*y;
    R[YY][XX] = t*x*y + s*z; R[YY][YY] = t*y*y + c;   R[YY][ZZ] = t*y*z - s*x;
    R[ZZ][XX] = t*x*z - s*y; R[ZZ][YY] = t*y*z + s*x; R[ZZ][ZZ] = t*z*z + c;
}

/*! \brief Rotate and translate coordinates
 * \param[in] R     The rotation matrix
 * \param[in] shift The translation, after the rotation
 * \param[in] x     The coordinates
 * \return the new coordinates
 */
std::vector<gmx::RVec> transformed(const matrix R, const gmx::RVec &shift,
                                   const std::vector<gmx::RVec> &x)
{
    std::vector<gmx::RVec> y(x.size());
    for (size_t i = 0; i < x.size(); i++)
    {
        mvmul(R, x[i], y[i]);
        rvec_inc(y[i], shift);
    }
    return y;
}

//! \return the indices of the atoms, not shells or virtual sites, in the order of the input
std::vector<size_t> atomIndices(const ACTMol &mol)
{
    std::vector<size_t> atoms;
    const auto         &particles = mol.atomsConst();
    for (size_t i = 0; i < particles.size(); i++)
    {
        if (particles[i].pType() == ActParticle::Atom)
        {
            atoms.push_back(i);
        }
    }
    return atoms;
}

/*! \brief The charge of each atom with those of its shells
 * Shells that belong to one atom are added to it, such that the charges
 * of molecules with their atoms in another order can be compared.
 * \param[in] mol The molecule
 * \return the charges, in the order of atomIndices
 */
std::vector<double> atomCharges(const ACTMol &mol)
{
    const auto         &particles = mol.atomsConst();
    std::vector<double> q(particles.size(), 0.0);
    for (size_t i = 0; i < particles.size(); i++)
    {
        size_t owner = i;
        if (particles[i].pType() == ActParticle::Shell && particles[i].parents().size() == 1)
        {
            owner = particles[i].parents()[0];
        }
        q[owner] += particles[i].charge() + particles[i].chargeVS();
    }
    std::vector<double> qatom;
    for (auto a : atomIndices(mol))
    {
        qatom.push_back(q[a]);
    }
    return qatom;
}

//! \return the total charge of all particles
double totalCharge(const ACTMol &mol)
{
    double qtot = 0;
    for (const auto &p : mol.atomsConst())
    {
        qtot += p.charge() + p.chargeVS();
    }
    return qtot;
}

class PhysicsTest : public ::testing::Test
{
protected:
    MsgHandler                  msghandler_;
    //! For temporary files, which are removed after the test
    gmx::test::TestFileManager  fileManager_;
    //! Tight convergence, such that the shells are at their minimum
    ForceComputer forceComp_ = ForceComputer(1e-12, 500, 1);
    //! Modified copies of force fields, see forceField
    std::list<ForceField>       copies_;

    /*! \brief A force field, possibly modified
     * \param[in] ffName The force field
     * \param[in] nexcl  If not negative, the number of excluded bonds for
     *                   the non-bonded interactions
     * \param[in] alg    If not NONE, the charge generation algorithm
     * \return the force field, a copy if it is modified
     */
    ForceField *forceField(const std::string        &ffName,
                           int                       nexcl = -1,
                           ChargeGenerationAlgorithm alg   = ChargeGenerationAlgorithm::NONE)
    {
        auto pd = getForceField(ffName);
        if (nullptr == pd || (nexcl < 0 && ChargeGenerationAlgorithm::NONE == alg))
        {
            return pd;
        }
        copies_.push_back(*pd);
        auto copy = &copies_.back();
        if (nexcl >= 0)
        {
            std::string exclstr("nexcl");
            for (auto &ffpl : *copy->forces())
            {
                if (ffpl.second.optionExists(exclstr))
                {
                    ffpl.second.removeOption(exclstr);
                    ffpl.second.addOption(exclstr, gmx::formatString("%d", nexcl));
                }
            }
        }
        if (ChargeGenerationAlgorithm::NONE != alg)
        {
            copy->setChargeGenerationAlgorithm(alg);
        }
        return copy;
    }

    /*! \brief Check that the topology has the interaction types
     * Such that a test of them is not empty without notice.
     * \param[in] mol    The molecule
     * \param[in] itypes The interaction types
     */
    static void expectEntries(const ACTMol &mol, const std::vector<InteractionType> &itypes)
    {
        for (auto itype : itypes)
        {
            EXPECT_TRUE(mol.topology()->hasEntry(itype)) << "no " << interactionTypeToString(itype);
        }
    }

    PhysicsTest()
    {
        msghandler_.setPrintLevel(ACTStatus::Warning);
    }

    /*! \brief Read a molecule, make its topology and generate the charges
     * \param[in]  pd      The force field
     * \param[in]  sdfFile The structure, a file name relative to this directory or absolute
     * \param[in]  qfrag   The charge of each compound, if empty all are neutral
     * \param[in]  x       The coordinates for the charges, if empty those of the file
     * \param[out] mol     The molecule
     * \return Whether it worked
     */
    bool readMolecule(ForceField                   *pd,
                      const std::string            &sdfFile,
                      const std::vector<double>    &qfrag,
                      const std::vector<gmx::RVec> &x,
                      ACTMol                       *mol)
    {
        std::vector<MolProp> molprops;
        double qtot = 0;
        for (auto q : qfrag)
        {
            qtot += q;
        }
        auto dataName = sdfFile[0] == '/' ? sdfFile : gmx::test::TestFileManager::getInputFilePath(sdfFile);
        importFile(&msghandler_, pd, dataName.c_str(), &molprops,
                   "minimum", JobType::OPT, !qfrag.empty(), &qtot, true);
        if (!msghandler_.ok() || molprops.size() != 1)
        {
            return false;
        }
        if (!qfrag.empty())
        {
            auto fptr = molprops[0].fragmentPtr();
            if (fptr->size() != qfrag.size())
            {
                return false;
            }
            for (size_t i = 0; i < qfrag.size(); i++)
            {
                (*fptr)[i].setCharge(qfrag[i]);
            }
        }
        mol->Merge(&molprops[0]);
        mol->GenerateTopology(&msghandler_, pd, missingParameters::Ignore);
        if (!msghandler_.ok())
        {
            return false;
        }
        std::vector<gmx::RVec> coords = x.empty() ? mol->xOriginal() : x;
        std::vector<gmx::RVec> forces(coords.size());
        mol->generateCharges(&msghandler_, pd, &forceComp_, pd->chargeGenerationAlgorithm(),
                             &coords, &forces);
        return msghandler_.ok();
    }

    /*! \brief Compute the energy and forces, with relaxed shells
     * \param[in]    pd  The force field
     * \param[in]    mol The molecule
     * \param[inout] x   The coordinates, on return with the virtual sites and relaxed shells
     * \param[out]   f   The forces, if not nullptr
     * \return the potential energy
     */
    double energy(ForceField             *pd,
                  const ACTMol           &mol,
                  std::vector<gmx::RVec> *x,
                  std::vector<gmx::RVec> *f = nullptr)
    {
        std::vector<gmx::RVec>            forces(x->size());
        std::map<InteractionType, double> e;
        forceComp_.compute(&msghandler_, pd, mol.topology(), x, &forces, &e);
        EXPECT_TRUE(msghandler_.ok());
        if (f)
        {
            *f = forces;
        }
        return e[InteractionType::EPOT];
    }

    //! \return the largest force component
    static double maxForce(const std::vector<gmx::RVec> &f)
    {
        double fmax = 0;
        for (const auto &fi : f)
        {
            for (int m = 0; m < DIM; m++)
            {
                fmax = std::max(fmax, std::abs(fi[m]));
            }
        }
        return fmax;
    }

    /*! \brief The forces on the atoms are minus the derivatives of the energy
     * With virtual sites, whose forces are spread over the atoms, and
     * shells, which are relaxed for each energy, such that the derivative
     * with respect to an atom is the force on it. The total force and
     * torque on the molecule are zero.
     * \param[in] ffName  The force field
     * \param[in] sdfFile The molecule
     * \param[in] itypes  Interaction types the topology should have
     */
    void testForces(const std::string                  &ffName,
                    const std::string                  &sdfFile,
                    const std::vector<InteractionType> &itypes = {})
    {
        auto pd = getForceField(ffName);
        ASSERT_NE(nullptr, pd);
        ACTMol mol;
        ASSERT_TRUE(readMolecule(pd, sdfFile, {}, {}, &mol));
        expectEntries(mol, itypes);
        auto                   x0 = mol.xOriginal();
        auto                   x  = x0;
        std::vector<gmx::RVec> f;
        energy(pd, mol, &x, &f);
        double fmax = maxForce(f);
        ASSERT_GT(fmax, 0);
        // Central differences, the error of which is of order h^2
        double h = 1e-5;
        for (auto a : atomIndices(mol))
        {
            for (int m = 0; m < DIM; m++)
            {
                auto xplus  = x0;
                auto xminus = x0;
                xplus[a][m]  += h;
                xminus[a][m] -= h;
                double fnum = -(energy(pd, mol, &xplus) - energy(pd, mol, &xminus))/(2*h);
                EXPECT_NEAR(fnum, f[a][m], 1e-5*std::max(1.0, fmax))
                    << ffName << " " << sdfFile << " atom " << a << " dim " << m;
            }
        }
        // Newton's third law, also for the virtual sites and shells, of
        // which the force is zero at the minimum
        gmx::RVec fsum = { 0, 0, 0 }, tsum = { 0, 0, 0 };
        for (size_t i = 0; i < f.size(); i++)
        {
            gmx::RVec t;
            cprod(x[i], f[i], t);
            rvec_inc(fsum, f[i]);
            rvec_inc(tsum, t);
        }
        for (int m = 0; m < DIM; m++)
        {
            EXPECT_NEAR(0, fsum[m], 1e-8*std::max(1.0, fmax)) << ffName << " total force, dim " << m;
            EXPECT_NEAR(0, tsum[m], 1e-8*std::max(1.0, fmax)) << ffName << " total torque, dim " << m;
        }
    }

    /*! \brief Rotating and translating the molecule changes nothing
     * The energy is the same and the forces rotate with the molecule. A
     * mirror image has the same energy as well, which also checks the
     * virtual sites out of a plane, which change places.
     * \param[in] ffName  The force field
     * \param[in] sdfFile The molecule
     * \param[in] reflect Whether to make a mirror image as well
     */
    void testRotation(const std::string &ffName, const std::string &sdfFile, bool reflect = false)
    {
        auto pd = getForceField(ffName);
        ASSERT_NE(nullptr, pd);
        ACTMol mol;
        ASSERT_TRUE(readMolecule(pd, sdfFile, {}, {}, &mol));
        auto                   x = mol.xOriginal();
        std::vector<gmx::RVec> f;
        double                 e    = energy(pd, mol, &x, &f);
        double                 fmax = maxForce(f);
        matrix                 R;
        rotationMatrix({ 0.3, -0.5, 0.8 }, 1.1, R);
        if (reflect)
        {
            // A mirror image in the plane x = 0, then the rotation
            for (int m = 0; m < DIM; m++)
            {
                R[m][XX] = -R[m][XX];
            }
        }
        gmx::RVec              shift = { 1.7, -0.6, 2.3 };
        auto                   y     = transformed(R, shift, mol.xOriginal());
        std::vector<gmx::RVec> g;
        double                 er = energy(pd, mol, &y, &g);
        EXPECT_NEAR(e, er, 1e-10*std::max(1.0, std::abs(e))) << ffName << " " << sdfFile;
        // The forces on the virtual sites are zero after spreading, also
        // when a mirror image makes them change places
        for (size_t i = 0; i < f.size(); i++)
        {
            gmx::RVec fr;
            mvmul(R, f[i], fr);
            for (int m = 0; m < DIM; m++)
            {
                EXPECT_NEAR(fr[m], g[i][m], 1e-8*std::max(1.0, fmax))
                    << ffName << " " << sdfFile << " particle " << i << " dim " << m;
            }
        }
    }

    /*! \brief Numbering the atoms in reverse order changes nothing
     * The energy, the forces on the atoms and their charges are the same.
     * The topologies are made separately, which also checks that the
     * interactions do not depend on the order of the atoms.
     * \param[in] ffName  The force field
     * \param[in] sdfFile The molecule
     * \param[in] qfrag   The charge of each compound, if empty all are neutral
     * \param[in] itypes  Interaction types the topology should have
     * \param[in] nexcl   If not negative, the number of excluded bonds
     */
    void testRenumbered(const std::string                  &ffName,
                        const std::string                  &sdfFile,
                        const std::vector<double>          &qfrag  = {},
                        const std::vector<InteractionType> &itypes = {},
                        int                                 nexcl  = -1)
    {
        auto pd = forceField(ffName, nexcl);
        ASSERT_NE(nullptr, pd);
        ACTMol mol;
        ASSERT_TRUE(readMolecule(pd, sdfFile, qfrag, {}, &mol));
        expectEntries(mol, itypes);
        auto label    = sdfFile.substr(sdfFile.find_last_of('/') + 1);
        auto reversed = fileManager_.getTemporaryFilePath(label);
        ASSERT_TRUE(writeReversedSdf(gmx::test::TestFileManager::getInputFilePath(sdfFile), reversed));
        // The compounds are in reverse order as well
        std::vector<double> qrev(qfrag.rbegin(), qfrag.rend());
        ACTMol rmol;
        ASSERT_TRUE(readMolecule(pd, reversed, qrev, {}, &rmol));
        expectEntries(rmol, itypes);
        // The same interactions
        for (const auto &entry : mol.topology()->entries())
        {
            ASSERT_TRUE(rmol.topology()->hasEntry(entry.first)) << interactionTypeToString(entry.first);
            EXPECT_EQ(entry.second.size(), rmol.topology()->entry(entry.first).size())
                << interactionTypeToString(entry.first);
        }
        auto atoms  = atomIndices(mol);
        auto ratoms = atomIndices(rmol);
        ASSERT_EQ(atoms.size(), ratoms.size());
        size_t natom = atoms.size();
        auto   q     = atomCharges(mol);
        auto   rq    = atomCharges(rmol);
        for (size_t k = 0; k < natom; k++)
        {
            EXPECT_NEAR(q[k], rq[natom-1-k], 1e-10) << ffName << " " << label << " charge of atom " << k;
        }
        auto                   x = mol.xOriginal();
        auto                   rx = rmol.xOriginal();
        std::vector<gmx::RVec> f, rf;
        double                 e    = energy(pd, mol, &x, &f);
        double                 re   = energy(pd, rmol, &rx, &rf);
        double                 fmax = maxForce(f);
        EXPECT_NEAR(e, re, 1e-10*std::max(1.0, std::abs(e))) << ffName << " " << label;
        for (size_t k = 0; k < natom; k++)
        {
            for (int m = 0; m < DIM; m++)
            {
                EXPECT_NEAR(f[atoms[k]][m], rf[ratoms[natom-1-k]][m], 1e-8*std::max(1.0, fmax))
                    << ffName << " " << label << " atom " << k << " dim " << m;
            }
        }
    }

    /*! \brief The generated charges add up to the total charge and do not
     * change when the molecule is rotated, translated or renumbered.
     * \param[in] ffName  The force field
     * \param[in] sdfFile The molecule
     * \param[in] qfrag   The charge of each compound
     * \param[in] alg     If not NONE, the charge generation algorithm
     */
    void testCharges(const std::string         &ffName,
                     const std::string         &sdfFile,
                     const std::vector<double> &qfrag,
                     ChargeGenerationAlgorithm  alg = ChargeGenerationAlgorithm::NONE)
    {
        auto pd = forceField(ffName, -1, alg);
        ASSERT_NE(nullptr, pd);
        double qtot = 0;
        for (auto qf : qfrag)
        {
            qtot += qf;
        }
        ACTMol mol;
        ASSERT_TRUE(readMolecule(pd, sdfFile, qfrag, {}, &mol));
        EXPECT_NEAR(qtot, totalCharge(mol), 1e-12) << ffName << " " << sdfFile;
        auto q = atomCharges(mol);
        // Rotated and translated
        matrix R;
        rotationMatrix({ -0.2, 0.9, 0.4 }, 2.3, R);
        ACTMol tmol;
        ASSERT_TRUE(readMolecule(pd, sdfFile, qfrag, transformed(R, { -3.1, 0.4, 1.2 }, mol.xOriginal()), &tmol));
        EXPECT_NEAR(qtot, totalCharge(tmol), 1e-12) << ffName << " " << sdfFile << " rotated";
        auto tq = atomCharges(tmol);
        ASSERT_EQ(q.size(), tq.size());
        for (size_t k = 0; k < q.size(); k++)
        {
            EXPECT_NEAR(q[k], tq[k], 1e-10) << ffName << " " << sdfFile << " rotated, charge of atom " << k;
        }
        // A mirror image
        matrix M = { { -1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        ACTMol mmol;
        ASSERT_TRUE(readMolecule(pd, sdfFile, qfrag, transformed(M, { 0, 0, 0 }, mol.xOriginal()), &mmol));
        auto mq = atomCharges(mmol);
        for (size_t k = 0; k < q.size(); k++)
        {
            EXPECT_NEAR(q[k], mq[k], 1e-10) << ffName << " " << sdfFile << " mirror image, charge of atom " << k;
        }
    }
};

// The molecules of the force and charge tests
const std::string waterDimer = "water_dimer.sdf";
const std::string hbrHf      = "hbr-hf.sdf";
const std::string tip5p      = "../../properties/tests/TIP5P";

// Test: Forces against finite differences of the energy, without shells
// and with the virtual sites VSITE2 and VSITE3FD
TEST_F(PhysicsTest, ForcesVsites)
{
    testForces("ACS-v3fds", waterDimer);
}

// Test: The same with the virtual sites of TIP5P, VSITE3OUT
TEST_F(PhysicsTest, ForcesVsite3OUT)
{
    testForces(tip5p, waterDimer);
}

// Test: The same with shells and virtual sites
TEST_F(PhysicsTest, ForcesShellsVsites)
{
    testForces("ACS-pg-vs2", hbrHf);
}

// Test: The same with shells
TEST_F(PhysicsTest, ForcesShells)
{
    testForces("ACS-pol", waterDimer);
}

// Test: The same with hyperpolarizable shells
TEST_F(PhysicsTest, ForcesHyperpolarizable)
{
    testForces("water_hyperpol", waterDimer);
}

// Test: Rotation and translation with virtual sites
TEST_F(PhysicsTest, RotationVsites)
{
    testRotation("ACS-v3fds", waterDimer);
}

// Test: Rotation and translation with the virtual sites of TIP5P
TEST_F(PhysicsTest, RotationVsite3OUT)
{
    testRotation(tip5p, waterDimer);
}

// Test: Rotation and translation with shells and virtual sites
TEST_F(PhysicsTest, RotationShellsVsites)
{
    testRotation("ACS-pg-vs2", hbrHf);
}

// Test: Rotation and translation with hyperpolarizable shells
TEST_F(PhysicsTest, RotationHyperpolarizable)
{
    testRotation("water_hyperpol", waterDimer);
}

// Test: Reverse numbering with virtual sites
TEST_F(PhysicsTest, RenumberedVsites)
{
    testRenumbered("ACS-v3fds", waterDimer);
}

// Test: Reverse numbering with the virtual sites of TIP5P
TEST_F(PhysicsTest, RenumberedVsite3OUT)
{
    testRenumbered(tip5p, waterDimer);
}

// Test: Reverse numbering with shells and virtual sites
TEST_F(PhysicsTest, RenumberedShellsVsites)
{
    testRenumbered("ACS-pg-vs2", hbrHf);
}

// Test: Reverse numbering with shells
TEST_F(PhysicsTest, RenumberedShells)
{
    testRenumbered("ACS-pol", waterDimer);
}

// Test: Charges of a neutral dimer with point charges, SQE
TEST_F(PhysicsTest, ChargesMethanolWater)
{
    testCharges("../../qgen/tests/ACS-g", "../../qgen/tests/methanol-water.sdf", { 0, 0 });
}

// Test: Charges of a charged dimer with Gaussian charges and shells
TEST_F(PhysicsTest, ChargesAcetateWaterPolarizable)
{
    testCharges("../../qgen/tests/ACS-pg", "../../qgen/tests/acetate-water.sdf", { -1, 0 });
}

// Test: Charges of a charged molecule
TEST_F(PhysicsTest, ChargesGuanidinium)
{
    testCharges("../../qgen/tests/ACS-g", "../../qgen/tests/guanidinium.sdf", { 1 });
}

// Test: Reverse numbering of a charged dimer, energies, forces and charges
TEST_F(PhysicsTest, RenumberedAcetateWater)
{
    testRenumbered("../../qgen/tests/ACS-pg", "../../qgen/tests/acetate-water.sdf", { -1, 0 });
}

// Test: Reverse numbering of a neutral dimer with point charges
TEST_F(PhysicsTest, RenumberedMethanolWater)
{
    testRenumbered("../../qgen/tests/ACS-g", "../../qgen/tests/methanol-water.sdf", { 0, 0 });
}

// The molecules of the topology tests
const std::string diamine = "../../alexandria/tests/mols/ethane-12-diamine.sdf";
const std::string water   = "water.sdf";
const std::string acsG    = "../../alexandria/tests/ACS-g";

// Test: Forces with the virtual sites VSITE3 and VSITE3OUT of water
TEST_F(PhysicsTest, ForcesVsite3Water)
{
    testForces("ACS-pg-vs3", water, { InteractionType::VSITE3, InteractionType::VSITE3OUT });
}

// Test: Mirror images, with the lone pairs of TIP5P changing places
TEST_F(PhysicsTest, ReflectionVsite3OUT)
{
    testRotation(tip5p, waterDimer, true);
}

// Test: Mirror images with shells and virtual sites
TEST_F(PhysicsTest, ReflectionShellsVsites)
{
    testRotation("ACS-pg-vs2", hbrHf, true);
}

// Test: Reverse numbering with VSITE3 and VSITE3OUT
TEST_F(PhysicsTest, RenumberedVsite3Water)
{
    testRenumbered("ACS-pg-vs3", water, {}, { InteractionType::VSITE3, InteractionType::VSITE3OUT });
}

// Test: Reverse numbering with all exclusions settings, nexcl 0 to 3
TEST_F(PhysicsTest, RenumberedExclusions)
{
    for (int nexcl = 0; nexcl <= 3; nexcl++)
    {
        SCOPED_TRACE(gmx::formatString("nexcl %d", nexcl));
        testRenumbered(acsG, diamine, {},
                       { InteractionType::BONDS, InteractionType::ANGLES, InteractionType::PROPER_DIHEDRALS },
                       nexcl);
    }
}

// Test: Charges with EEM rather than SQE
TEST_F(PhysicsTest, ChargesEEM)
{
    testCharges("../../qgen/tests/ACS-g", "../../qgen/tests/methanol-water.sdf", { 0, 0 },
                ChargeGenerationAlgorithm::EEM);
}

// Test: Charges of a charged dimer with EEM, Gaussian charges and shells
TEST_F(PhysicsTest, ChargesEEMPolarizable)
{
    testCharges("../../qgen/tests/ACS-pg", "../../qgen/tests/acetate-water.sdf", { -1, 0 },
                ChargeGenerationAlgorithm::EEM);
}

// Test: Charges of the ACM force fields
TEST_F(PhysicsTest, ChargesACM)
{
    testCharges("../../qgen/tests/ACM-g", "../../qgen/tests/acetate-water.sdf", { -1, 0 });
    testCharges("../../qgen/tests/ACM-pg", "../../qgen/tests/methanol-water.sdf", { 0, 0 });
}

// Test: Energies are used again when nothing changed since the previous
// calculation, and computed again, the same as without them, after the
// coordinates, the charges or the parameters change
TEST_F(PhysicsTest, EnergyCache)
{
    // Energies are only used again without polarization
    auto pd0 = getForceField("ACS-v3fds");
    ASSERT_NE(nullptr, pd0);
    ASSERT_FALSE(pd0->polarizable());
    // A copy, of which a parameter is changed
    ForceField pd = *pd0;
    ACTMol     mol;
    ASSERT_TRUE(readMolecule(&pd, waterDimer, {}, {}, &mol));
    auto                      x0 = mol.xOriginal();
    std::vector<CoordsEnergy> coords(5);
    for (size_t c = 0; c < coords.size(); c++)
    {
        coords[c].x = x0;
        for (size_t i = 0; i < x0.size(); i++)
        {
            for (int m = 0; m < DIM; m++)
            {
                coords[c].x[i][m] += 0.005*std::sin(1.0 + c + 3*i + m);
            }
        }
    }
    // The potential energies without forces, as in training
    auto epot = [&](ConformationData *cdata)
    {
        computeConformations(&msghandler_, &pd, mol.topology(), &coords, cdata, false, false);
        EXPECT_TRUE(msghandler_.ok());
        auto                e = cdata->energies.find(InteractionType::EPOT);
        std::vector<double> ep;
        for (size_t c = 0; c < coords.size(); c++)
        {
            ep.push_back(e.has(c) ? e.act[c] : 0.0);
        }
        return ep;
    };
    auto fresh = [&]()
    {
        ConformationData cdata;
        return epot(&cdata);
    };
    auto  &timings  = Timings::instance();
    size_t reusedId = timings.index("Energies used again");
    auto   reused   = [&timings, reusedId]() { return timings.calls()[reusedId]; };

    ConformationData cdata;
    auto             e1 = epot(&cdata);
    double           n0 = reused();
    auto             e2 = epot(&cdata);
    EXPECT_GT(reused(), n0) << "energies used again";
    EXPECT_EQ(e1, e2);
    EXPECT_EQ(fresh(), e2);

    // Coordinates, the change is made known to the conformation data
    for (auto &ce : coords)
    {
        ce.x[0][XX] += 0.01;
    }
    cdata.invalidate();
    n0      = reused();
    auto e3 = epot(&cdata);
    EXPECT_EQ(n0, reused()) << "nothing used again after the coordinates changed";
    EXPECT_EQ(fresh(), e3);
    EXPECT_NE(e1, e3);

    // Charges, keeping the total charge
    auto  &atoms = *mol.topologyPtr()->atomsPtr();
    size_t qi = atoms.size(), qj = atoms.size();
    for (size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].charge() != 0)
        {
            (qi == atoms.size() ? qi : qj) = i;
            if (qj < atoms.size())
            {
                break;
            }
        }
    }
    ASSERT_LT(qj, atoms.size());
    atoms[qi].setCharge(atoms[qi].charge() + 0.05);
    atoms[qj].setCharge(atoms[qj].charge() - 0.05);
    auto e4 = epot(&cdata);
    EXPECT_EQ(fresh(), e4);
    EXPECT_NE(e3, e4);

    // A bond length, filled in the topology
    auto  itype  = InteractionType::BONDS;
    ASSERT_TRUE(mol.topology()->hasEntry(itype));
    auto  id     = mol.topology()->entry(itype)[0]->id();
    auto *params = pd.findForces(itype)->findParameters(id);
    ASSERT_NE(nullptr, params);
    ASSERT_FALSE(params->empty());
    auto &param  = params->begin()->second;
    param.forceSetValue(1.1*param.value());
    mol.topologyPtr()->fillParameters(&msghandler_, &pd, missingParameters::Ignore);
    auto e5 = epot(&cdata);
    EXPECT_EQ(fresh(), e5);
    EXPECT_NE(e4, e5);
}

} // namespace

} // namespace alexandria
