/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2025
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
 * Unit tests for qPropertyType free functions and QtypeProps class methods.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cmath>

#include <gtest/gtest.h>

#include "act/qgen/qtype.h"
#include "act/molprop/multipole_names.h"

#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

// ============================================================
// Tests for qPropertyType free functions
// ============================================================

TEST(QPropertyTypeTest, NameForACM)
{
    EXPECT_EQ("ACM", qPropertyTypeName(qPropertyType::ACM));
}

TEST(QPropertyTypeTest, NameForESP)
{
    EXPECT_EQ("ESP", qPropertyTypeName(qPropertyType::ESP));
}

TEST(QPropertyTypeTest, NameForElec)
{
    EXPECT_EQ("Electronic", qPropertyTypeName(qPropertyType::Elec));
}

TEST(QPropertyTypeTest, StringToQtypeACM)
{
    EXPECT_EQ(qPropertyType::ACM, stringToQtype("ACM"));
}

TEST(QPropertyTypeTest, StringToQtypeESP)
{
    EXPECT_EQ(qPropertyType::ESP, stringToQtype("ESP"));
}

TEST(QPropertyTypeTest, StringToQtypeElectronic)
{
    EXPECT_EQ(qPropertyType::Elec, stringToQtype("Electronic"));
}

TEST(QPropertyTypeTest, StringToQtypeInvalidThrows)
{
    EXPECT_THROW(stringToQtype("InvalidType"), gmx::InvalidInputError);
}

TEST(QPropertyTypeTest, StringToQtypeEmptyThrows)
{
    EXPECT_THROW(stringToQtype(""), gmx::InvalidInputError);
}

TEST(QPropertyTypeTest, QPropertyTypesMapComplete)
{
    auto &types = qPropertyTypes();
    EXPECT_EQ(3u, types.size());
    EXPECT_TRUE(types.find(qPropertyType::ACM)  != types.end());
    EXPECT_TRUE(types.find(qPropertyType::ESP)  != types.end());
    EXPECT_TRUE(types.find(qPropertyType::Elec) != types.end());
}

TEST(QPropertyTypeTest, RoundTripConversion)
{
    for (const auto &entry : qPropertyTypes())
    {
        EXPECT_EQ(entry.first, stringToQtype(entry.second));
        EXPECT_EQ(entry.second, qPropertyTypeName(entry.first));
    }
}

// ============================================================
// Tests for QtypeProps - empty constructor and basic operations
// ============================================================

TEST(QtypePropsTest, EmptyConstructorDefaults)
{
    QtypeProps qp;
    EXPECT_EQ(qPropertyType::ACM, qp.qtype());
    EXPECT_TRUE(qp.charge().empty());
    EXPECT_TRUE(qp.x().empty());
    EXPECT_FALSE(qp.hasPolarizability());
    EXPECT_DOUBLE_EQ(0.0, qp.isotropicPolarizability());
    EXPECT_DOUBLE_EQ(0.0, qp.anisotropicPolarizability());
}

TEST(QtypePropsTest, SetQtype)
{
    QtypeProps qp;
    qp.setQtype(qPropertyType::ESP);
    EXPECT_EQ(qPropertyType::ESP, qp.qtype());
    qp.setQtype(qPropertyType::Elec);
    EXPECT_EQ(qPropertyType::Elec, qp.qtype());
}

TEST(QtypePropsTest, HasMultipoleReturnsFalseBeforeInit)
{
    QtypeProps qp;
    EXPECT_FALSE(qp.hasMultipole(MolPropObservable::DIPOLE));
    EXPECT_FALSE(qp.hasMultipole(MolPropObservable::QUADRUPOLE));
    EXPECT_FALSE(qp.hasMultipole(MolPropObservable::OCTUPOLE));
    EXPECT_FALSE(qp.hasMultipole(MolPropObservable::HEXADECAPOLE));
}

TEST(QtypePropsTest, HasMultipoleReturnsTrueAfterInit)
{
    QtypeProps qp;
    qp.initializeMoments();
    EXPECT_TRUE(qp.hasMultipole(MolPropObservable::DIPOLE));
    EXPECT_TRUE(qp.hasMultipole(MolPropObservable::QUADRUPOLE));
    EXPECT_TRUE(qp.hasMultipole(MolPropObservable::OCTUPOLE));
    EXPECT_TRUE(qp.hasMultipole(MolPropObservable::HEXADECAPOLE));
}

TEST(QtypePropsTest, GetMultipoleThrowsWhenNotPresent)
{
    QtypeProps qp;
    EXPECT_THROW(qp.getMultipole(MolPropObservable::DIPOLE), gmx::InternalError);
}

TEST(QtypePropsTest, DipoleThrowsWhenNotPresent)
{
    QtypeProps qp;
    EXPECT_THROW(qp.dipole(), gmx::InternalError);
}

TEST(QtypePropsTest, SetMultipole)
{
    QtypeProps qp;
    std::vector<double> dipValues = { 1.0, 2.0, 3.0 };
    qp.setMultipole(MolPropObservable::DIPOLE, dipValues);
    EXPECT_TRUE(qp.hasMultipole(MolPropObservable::DIPOLE));
    auto retrieved = qp.getMultipole(MolPropObservable::DIPOLE);
    EXPECT_EQ(3u, retrieved.size());
    EXPECT_DOUBLE_EQ(1.0, retrieved[0]);
    EXPECT_DOUBLE_EQ(2.0, retrieved[1]);
    EXPECT_DOUBLE_EQ(3.0, retrieved[2]);
}

TEST(QtypePropsTest, SetMultipoleTwiceThrows)
{
    QtypeProps qp;
    std::vector<double> dipValues = { 1.0, 2.0, 3.0 };
    qp.setMultipole(MolPropObservable::DIPOLE, dipValues);
    EXPECT_THROW(qp.setMultipole(MolPropObservable::DIPOLE, dipValues), gmx::InternalError);
}

TEST(QtypePropsTest, DipoleComputation)
{
    QtypeProps qp;
    std::vector<double> dipValues = { 3.0, 4.0, 0.0 };
    qp.setMultipole(MolPropObservable::DIPOLE, dipValues);
    EXPECT_DOUBLE_EQ(5.0, qp.dipole());
}

// ============================================================
// Tests for QtypeProps - polarizability methods
// ============================================================

TEST(QtypePropsTest, SetPolarizabilityTensor)
{
    QtypeProps qp;
    EXPECT_FALSE(qp.hasPolarizability());

    tensor alpha = { { 0 } };
    alpha[XX][XX] = 10.0;
    alpha[YY][YY] = 20.0;
    alpha[ZZ][ZZ] = 30.0;
    alpha[XX][YY] = 0.0;
    alpha[XX][ZZ] = 0.0;
    alpha[YY][ZZ] = 0.0;

    qp.setPolarizabilityTensor(alpha);
    EXPECT_TRUE(qp.hasPolarizability());
}

TEST(QtypePropsTest, IsotropicPolarizability)
{
    QtypeProps qp;
    tensor alpha = { { 0 } };
    alpha[XX][XX] = 10.0;
    alpha[YY][YY] = 20.0;
    alpha[ZZ][ZZ] = 30.0;
    qp.setPolarizabilityTensor(alpha);
    // Isotropic = (xx + yy + zz) / 3 = (10 + 20 + 30) / 3 = 20
    EXPECT_DOUBLE_EQ(20.0, qp.isotropicPolarizability());
}

TEST(QtypePropsTest, IsotropicPolarizabilityUniform)
{
    QtypeProps qp;
    tensor alpha = { { 0 } };
    alpha[XX][XX] = 5.0;
    alpha[YY][YY] = 5.0;
    alpha[ZZ][ZZ] = 5.0;
    qp.setPolarizabilityTensor(alpha);
    EXPECT_DOUBLE_EQ(5.0, qp.isotropicPolarizability());
    // Anisotropy should be zero for uniform tensor
    EXPECT_NEAR(0.0, qp.anisotropicPolarizability(), 1e-10);
}

TEST(QtypePropsTest, AnisotropicPolarizability)
{
    QtypeProps qp;
    tensor alpha = { { 0 } };
    alpha[XX][XX] = 10.0;
    alpha[YY][YY] = 20.0;
    alpha[ZZ][ZZ] = 30.0;
    qp.setPolarizabilityTensor(alpha);

    // anisotropy = sqrt(0.5) * sqrt((xx-yy)^2 + (xx-zz)^2 + (zz-yy)^2 + 6*(xy^2+xz^2+zy^2))
    // = sqrt(0.5) * sqrt(100 + 400 + 100 + 0) = sqrt(0.5) * sqrt(600)
    double expected = std::sqrt(0.5) * std::sqrt(600.0);
    EXPECT_NEAR(expected, qp.anisotropicPolarizability(), 1e-10);
}

TEST(QtypePropsTest, PolarizabilityTensorRoundTrip)
{
    QtypeProps qp;
    tensor alpha = { { 0 } };
    alpha[XX][XX] = 1.0;
    alpha[YY][YY] = 2.0;
    alpha[ZZ][ZZ] = 3.0;
    alpha[XX][YY] = 0.5;
    alpha[YY][XX] = 0.5;
    alpha[XX][ZZ] = 0.1;
    alpha[ZZ][XX] = 0.1;
    alpha[YY][ZZ] = 0.2;
    alpha[ZZ][YY] = 0.2;

    qp.setPolarizabilityTensor(alpha);
    const tensor &retrieved = qp.polarizabilityTensor();
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            EXPECT_DOUBLE_EQ(alpha[i][j], retrieved[i][j]);
        }
    }
}

TEST(QtypePropsTest, AnisotropicPolarizabilityWithOffDiagonal)
{
    QtypeProps qp;
    tensor alpha = { { 0 } };
    alpha[XX][XX] = 10.0;
    alpha[YY][YY] = 10.0;
    alpha[ZZ][ZZ] = 10.0;
    alpha[XX][YY] = 1.0;
    alpha[XX][ZZ] = 2.0;
    alpha[ZZ][YY] = 3.0;

    qp.setPolarizabilityTensor(alpha);
    // Diagonal terms equal -> only off-diagonal contribution:
    // anisotropy = sqrt(0.5) * sqrt(0 + 0 + 0 + 6*(1 + 4 + 9))
    //            = sqrt(0.5) * sqrt(84)
    double expected = std::sqrt(0.5) * std::sqrt(84.0);
    EXPECT_NEAR(expected, qp.anisotropicPolarizability(), 1e-10);
}

// ============================================================
// Tests for QtypeProps - setCenterOfCharge
// ============================================================

TEST(QtypePropsTest, SetCenterOfCharge)
{
    QtypeProps qp;
    rvec coc = { 1.0, 2.0, 3.0 };
    // Should not throw
    qp.setCenterOfCharge(coc);
}

// ============================================================
// Tests for QtypeProps - setQ and setX
// ============================================================

TEST(QtypePropsTest, SetQVector)
{
    QtypeProps qp;
    std::vector<double> q = { 0.5, -0.25, -0.25 };
    qp.setQ(q);
    auto charges = qp.charge();
    EXPECT_EQ(3u, charges.size());
    EXPECT_DOUBLE_EQ(0.5,   charges[0]);
    EXPECT_DOUBLE_EQ(-0.25, charges[1]);
    EXPECT_DOUBLE_EQ(-0.25, charges[2]);
}

TEST(QtypePropsTest, SetQandXVectorRequiresAtomNumbers)
{
    // setQandX(q, x) calls computeCoC() which requires atomNumber_
    // to be set. Calling it on a fresh QtypeProps without atomNumber_
    // raises an error because atomNumber_.size() != x_.size().
    QtypeProps qp;
    std::vector<double> q = { 1.0, -0.5, -0.5 };
    std::vector<gmx::RVec> x = {
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 }
    };
    EXPECT_THROW(qp.setQandX(q, x), gmx::InternalError);
}

// ============================================================
// Tests for QtypeProps - moments with simple data
// ============================================================

TEST(QtypePropsTest, InitializeMomentsCreatesAll)
{
    QtypeProps qp;
    qp.initializeMoments();
    // Check that all 4 multipole types are initialized
    for (const auto &mpo : mpoMultiPoles)
    {
        EXPECT_TRUE(qp.hasMultipole(mpo));
        auto values = qp.getMultipole(mpo);
        // All values should be zero after initialization
        for (const auto &v : values)
        {
            EXPECT_DOUBLE_EQ(0.0, v);
        }
    }
}

}

}
