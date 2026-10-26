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
 * Tests for the rotational symmetry number.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "act/properties/symmetrynumber.h"

namespace alexandria
{

namespace
{

//! A test structure
struct Structure
{
    //! Atomic numbers
    std::vector<int>       Z;
    //! Coordinates (nm)
    std::vector<gmx::RVec> x;

    //! Add an atom
    void add(int z, double x0, double y0, double z0)
    {
        Z.push_back(z);
        x.push_back({ static_cast<real>(x0), static_cast<real>(y0), static_cast<real>(z0) });
    }
};

//! \return a mass for each element (the value is not important)
std::vector<double> massesOf(const Structure &s)
{
    std::vector<double> m;
    for(auto z : s.Z)
    {
        m.push_back(z == 1 ? 1.008 : 2.0*z);
    }
    return m;
}

//! \return the symmetry number of a structure
int sigma(const Structure &s, double tolerance = 0.01)
{
    return rotationalSymmetryNumber(s.Z, massesOf(s), s.x, tolerance);
}

//! Tetrahedral directions
const double tet[4][3] = { { 1, 1, 1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 } };

//! \return a tetrahedral XY4 molecule
Structure tetrahedral(int center, const std::vector<int> &ligands, double r = 0.11)
{
    Structure s;
    s.add(center, 0, 0, 0);
    double f = r/std::sqrt(3.0);
    for(int i = 0; i < 4; i++)
    {
        s.add(ligands[i], f*tet[i][0], f*tet[i][1], f*tet[i][2]);
    }
    return s;
}

//! \return a regular polygon of atoms in the xy plane with hydrogens further out
Structure ring(int n, int z, double r = 0.14, double rH = 0.25)
{
    Structure s;
    for(int i = 0; i < n; i++)
    {
        double phi = 2*M_PI*i/n;
        s.add(z, r*std::cos(phi), r*std::sin(phi), 0);
        s.add(1, rH*std::cos(phi), rH*std::sin(phi), 0);
    }
    return s;
}

TEST(SymmetryNumberTest, SingleAtom)
{
    Structure s;
    s.add(18, 0, 0, 0);
    EXPECT_EQ(1, sigma(s));
}

TEST(SymmetryNumberTest, HeteronuclearDiatomic)
{
    Structure s;
    s.add(1, 0, 0, 0);
    s.add(17, 0, 0, 0.127);
    EXPECT_EQ(1, sigma(s));
}

TEST(SymmetryNumberTest, HomonuclearDiatomic)
{
    Structure s;
    s.add(7, 0, 0, -0.055);
    s.add(7, 0, 0, 0.055);
    EXPECT_EQ(2, sigma(s));
}

TEST(SymmetryNumberTest, CarbonDioxide)
{
    Structure s;
    s.add(6, 0, 0, 0);
    s.add(8, 0, 0, 0.116);
    s.add(8, 0, 0, -0.116);
    EXPECT_EQ(2, sigma(s));
}

TEST(SymmetryNumberTest, HydrogenCyanide)
{
    Structure s;
    s.add(1, 0, 0, -0.107);
    s.add(6, 0, 0, 0);
    s.add(7, 0, 0, 0.116);
    EXPECT_EQ(1, sigma(s));
}

TEST(SymmetryNumberTest, Water)
{
    Structure s;
    s.add(8, 0, 0, 0);
    s.add(1, 0.0757, 0, 0.0586);
    s.add(1, -0.0757, 0, 0.0586);
    EXPECT_EQ(2, sigma(s));
}

TEST(SymmetryNumberTest, Ammonia)
{
    Structure s;
    s.add(7, 0, 0, 0);
    for(int i = 0; i < 3; i++)
    {
        double phi = 2*M_PI*i/3;
        s.add(1, 0.094*std::cos(phi), 0.094*std::sin(phi), -0.038);
    }
    EXPECT_EQ(3, sigma(s));
}

TEST(SymmetryNumberTest, Methane)
{
    EXPECT_EQ(12, sigma(tetrahedral(6, { 1, 1, 1, 1 })));
}

TEST(SymmetryNumberTest, Chloromethane)
{
    // C3v
    EXPECT_EQ(3, sigma(tetrahedral(6, { 17, 1, 1, 1 })));
}

TEST(SymmetryNumberTest, Dichloromethane)
{
    // C2v
    EXPECT_EQ(2, sigma(tetrahedral(6, { 17, 17, 1, 1 })));
}

TEST(SymmetryNumberTest, MirrorPlaneOnly)
{
    // CH2FCl has a mirror plane (Cs) but no rotation axis, the
    // reflection should not be counted.
    EXPECT_EQ(1, sigma(tetrahedral(6, { 9, 17, 1, 1 })));
}

TEST(SymmetryNumberTest, Chiral)
{
    // CHFClBr has no symmetry
    EXPECT_EQ(1, sigma(tetrahedral(6, { 1, 9, 17, 35 })));
}

TEST(SymmetryNumberTest, SulfurHexafluoride)
{
    Structure s;
    s.add(16, 0, 0, 0);
    for(int m = 0; m < 3; m++)
    {
        for(double sign : { 1.0, -1.0 })
        {
            double v[3] = { 0, 0, 0 };
            v[m] = sign*0.156;
            s.add(9, v[0], v[1], v[2]);
        }
    }
    EXPECT_EQ(24, sigma(s));
}

TEST(SymmetryNumberTest, Benzene)
{
    EXPECT_EQ(12, sigma(ring(6, 6)));
}

TEST(SymmetryNumberTest, Cyclopentadienyl)
{
    // D5h
    EXPECT_EQ(10, sigma(ring(5, 6)));
}

TEST(SymmetryNumberTest, StaggeredEthane)
{
    // D3d
    Structure s;
    s.add(6, 0, 0, 0.0765);
    s.add(6, 0, 0, -0.0765);
    for(int i = 0; i < 3; i++)
    {
        double phi1 = 2*M_PI*i/3;
        double phi2 = phi1 + M_PI/3;
        s.add(1, 0.102*std::cos(phi1), 0.102*std::sin(phi1), 0.115);
        s.add(1, 0.102*std::cos(phi2), 0.102*std::sin(phi2), -0.115);
    }
    EXPECT_EQ(6, sigma(s));
}

TEST(SymmetryNumberTest, Allene)
{
    // D2d: the CH2 groups are perpendicular
    Structure s;
    s.add(6, 0, 0, 0);
    s.add(6, 0, 0, 0.131);
    s.add(6, 0, 0, -0.131);
    s.add(1, 0.093, 0, 0.185);
    s.add(1, -0.093, 0, 0.185);
    s.add(1, 0, 0.093, -0.185);
    s.add(1, 0, -0.093, -0.185);
    EXPECT_EQ(4, sigma(s));
}

TEST(SymmetryNumberTest, Translated)
{
    // The result does not depend on the position of the molecule
    auto s = tetrahedral(6, { 1, 1, 1, 1 });
    for(auto &x : s.x)
    {
        x[XX] += 1.3;
        x[YY] -= 0.7;
        x[ZZ] += 2.1;
    }
    EXPECT_EQ(12, sigma(s));
}

TEST(SymmetryNumberTest, Tolerance)
{
    // A small distortion is within the tolerance
    auto s = tetrahedral(6, { 1, 1, 1, 1 });
    s.x[1][XX] += 0.002;
    EXPECT_EQ(12, sigma(s));
    // A large one is not
    s.x[1][XX] += 0.05;
    EXPECT_EQ(1, sigma(s));
    // Stretching one bond along its direction leaves a C3 axis
    auto c3v = tetrahedral(6, { 1, 1, 1, 1 });
    for(int m = 0; m < 3; m++)
    {
        c3v.x[1][m] += 0.05*tet[0][m]/std::sqrt(3.0);
    }
    EXPECT_EQ(3, sigma(c3v));
}

TEST(SymmetryNumberTest, SearchLimit)
{
    // With a too small search limit the result is unknown
    auto s = ring(6, 6);
    EXPECT_EQ(0, rotationalSymmetryNumber(s.Z, massesOf(s), s.x, 0.01, 3));
}

} // namespace

} // namespace alexandria
