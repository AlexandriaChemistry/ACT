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
 * Tests of the removal of the motion of the center of mass.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <vector>

#include <gtest/gtest.h>

#include "act/properties/velocityhandler.h"
#include "act/topology/actatom.h"
#include "gromacs/utility/exceptions.h"

namespace alexandria
{

namespace
{

//! Make atoms with the masses given, a mass of zero is a shell
std::vector<ActAtom> makeAtoms(const std::vector<double> &mass)
{
    std::vector<ActAtom> atoms;
    for (auto m : mass)
    {
        atoms.push_back(ActAtom("A", "C", "c", m > 0 ? ActParticle::Atom : ActParticle::Shell,
                                6, m, 0, 2));
    }
    return atoms;
}

//! Deterministic velocities, different for each particle and dimension
std::vector<gmx::RVec> makeVelocities(size_t n)
{
    std::vector<gmx::RVec> v;
    for (size_t i = 0; i < n; i++)
    {
        v.push_back({ static_cast<real>(0.3 + 0.7*i), static_cast<real>(-0.5 + 0.2*i*i),
                      static_cast<real>(1.1 - 0.4*i) });
    }
    return v;
}

//! The linear momentum and the angular momentum around the center of mass
void momentum(const std::vector<ActAtom>   &atoms,
              const std::vector<gmx::RVec> &x,
              const std::vector<gmx::RVec> &v,
              dvec                          p,
              dvec                          L)
{
    double tm  = 0;
    dvec   xcm = { 0, 0, 0 };
    clear_dvec(p);
    clear_dvec(L);
    for (size_t i = 0; i < atoms.size(); i++)
    {
        double m = atoms[i].mass();
        tm += m;
        for (int d = 0; d < DIM; d++)
        {
            xcm[d] += m*x[i][d];
            p[d]   += m*v[i][d];
        }
    }
    for (size_t i = 0; i < atoms.size(); i++)
    {
        double m = atoms[i].mass();
        dvec   dx, dv, a;
        for (int d = 0; d < DIM; d++)
        {
            dx[d] = x[i][d] - xcm[d]/tm;
            dv[d] = v[i][d] - p[d]/tm;
        }
        dcprod(dx, dv, a);
        for (int d = 0; d < DIM; d++)
        {
            L[d] += m*a[d];
        }
    }
}

//! Expect two lists of velocities to be identical
void expectEqual(const std::vector<gmx::RVec> &a, const std::vector<gmx::RVec> &b)
{
    ASSERT_EQ(a.size(), b.size());
    for (size_t i = 0; i < a.size(); i++)
    {
        for (int d = 0; d < DIM; d++)
        {
            EXPECT_EQ(a[i][d], b[i][d]) << i << " " << d;
        }
    }
}

//! The kinetic energy
double kineticEnergy(const std::vector<ActAtom> &atoms, const std::vector<gmx::RVec> &v)
{
    double ekin = 0;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        ekin += 0.5*atoms[i].mass()*iprod(v[i], v[i]);
    }
    return ekin;
}

//! A non-linear molecule with a shell
class COMMRemoverTest : public ::testing::Test
{
protected:
    std::vector<ActAtom>   atoms_ = makeAtoms({ 12, 1, 16, 0, 14 });
    std::vector<gmx::RVec> x_     = { { 0.1, 0.2, 0.0 }, { 0.2, 0.1, 0.05 },
                                      { -0.1, 0.0, 0.12 }, { -0.1, 0.01, 0.12 },
                                      { 0.0, -0.15, -0.1 } };
};

// Test: Angular removal leaves no linear or angular momentum, and does
// not change the velocities of shells
TEST_F(COMMRemoverTest, Angular)
{
    auto v = makeVelocities(x_.size());
    auto v0 = v;
    COMMRemover remover(COMMRemoval::Angular, atoms_, x_);
    EXPECT_EQ(6, remover.nDOF());
    remover.remove(x_, &v);
    dvec p, L;
    momentum(atoms_, x_, v, p, L);
    for (int d = 0; d < DIM; d++)
    {
        EXPECT_NEAR(0, p[d], 1e-5) << d;
        EXPECT_NEAR(0, L[d], 1e-5) << d;
    }
    expectEqual({ v0[3] }, { v[3] });
    // Removing again changes nothing
    auto v2 = v;
    remover.remove(x_, &v2);
    for (size_t i = 0; i < v.size(); i++)
    {
        for (int d = 0; d < DIM; d++)
        {
            EXPECT_NEAR(v[i][d], v2[i][d], 1e-5) << i << " " << d;
        }
    }
}

// Test: Linear removal leaves no linear momentum but keeps the angular
// momentum, and no removal changes nothing
TEST_F(COMMRemoverTest, LinearAndNone)
{
    auto v0 = makeVelocities(x_.size());
    dvec p0, L0;
    momentum(atoms_, x_, v0, p0, L0);
    auto v  = v0;
    COMMRemover linear(COMMRemoval::Linear, atoms_, x_);
    EXPECT_EQ(3, linear.nDOF());
    linear.remove(x_, &v);
    dvec p, L;
    momentum(atoms_, x_, v, p, L);
    for (int d = 0; d < DIM; d++)
    {
        EXPECT_NEAR(0, p[d], 1e-5) << d;
        EXPECT_NEAR(L0[d], L[d], 1e-5) << d;
    }
    COMMRemover none(COMMRemoval::None, atoms_, x_);
    EXPECT_EQ(0, none.nDOF());
    auto vn = v0;
    none.remove(x_, &vn);
    expectEqual(v0, vn);
}

// Test: The rotation that is removed is orthogonal to the remaining
// velocities, so the kinetic energy decreases by that of the rotation
TEST_F(COMMRemoverTest, KineticEnergy)
{
    auto v = makeVelocities(x_.size());
    COMMRemover linear(COMMRemoval::Linear, atoms_, x_);
    COMMRemover angular(COMMRemoval::Angular, atoms_, x_);
    auto vl = v, va = v;
    linear.remove(x_, &vl);
    angular.remove(x_, &va);
    // The difference between linear and angular removal is the rigid rotation
    std::vector<gmx::RVec> rot(v.size());
    for (size_t i = 0; i < v.size(); i++)
    {
        rvec_sub(vl[i], va[i], rot[i]);
    }
    EXPECT_NEAR(kineticEnergy(atoms_, vl),
                kineticEnergy(atoms_, va) + kineticEnergy(atoms_, rot), 1e-5);
    EXPECT_LT(kineticEnergy(atoms_, va), kineticEnergy(atoms_, vl));
}

// Test: For a linear molecule the rotation around its axis is not removed
TEST(COMMRemoverLinearTest, LinearMolecule)
{
    auto atoms = makeAtoms({ 16, 12, 16 });
    std::vector<gmx::RVec> x = { { -0.116, 0, 0 }, { 0, 0, 0 }, { 0.116, 0, 0 } };
    COMMRemover remover(COMMRemoval::Angular, atoms, x);
    EXPECT_EQ(5, remover.nDOF());
    auto v = makeVelocities(x.size());
    remover.remove(x, &v);
    dvec p, L;
    momentum(atoms, x, v, p, L);
    for (int d = 0; d < DIM; d++)
    {
        EXPECT_NEAR(0, p[d], 1e-5) << d;
        EXPECT_NEAR(0, L[d], 1e-5) << d;
    }
    // A bent geometry has three rotations, a single atom none
    x[2][YY] = 0.05;
    EXPECT_EQ(6, COMMRemover(COMMRemoval::Angular, atoms, x).nDOF());
    EXPECT_EQ(3, COMMRemover(COMMRemoval::Angular, makeAtoms({ 12 }), { { 0, 0, 0 } }).nDOF());
}

// Test: The names of the removal types
TEST(COMMRemovalTest, Strings)
{
    for (auto cr : { COMMRemoval::None, COMMRemoval::Linear, COMMRemoval::Angular })
    {
        EXPECT_EQ(cr, stringToCOMMRemoval(commRemovalToString(cr)));
    }
    EXPECT_THROW(stringToCOMMRemoval("rotational"), gmx::InvalidInputError);
}

} // namespace

} // namespace alexandria
