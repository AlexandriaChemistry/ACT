/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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
 * Implements test of Coulomb routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include <config.h>

#include <math.h>

#include <map>

#include <gtest/gtest.h>

#include "act/basics/chargemodel.h"
#include "act/coulombintegrals/gaussian_integrals.h"
#include "act/coulombintegrals/slater_integrals.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace
{

void check_force(const std::vector<double> &potential,
                 const std::vector<double> &force,
                 double                     dx)
{
    // Check whether force is derivative of potential
    double toler = 0.15;
    // but only if the real force is large enough which is not always
    // the case for Slaters close to a distance of zero.
    double minforce = 0.2;
    for(size_t i = 1; i < force.size()-1; i++)
    {
        double deriv = -(potential[i+1]-potential[i-1])/(2*dx);
        double relerror = std::abs(deriv-force[i])/(std::abs(deriv)+std::abs(force[i]));
        if (force[i] > minforce)
        {
            if (relerror >= toler)
            {
                printf("Analytical force[%lu] %g Numerical force %g relerror %g\n"
                       "V %.8f %.8f %.8f F %.8f %.8f %.8f\n",
                       i, force[i], deriv, relerror,
                       potential[i-1], potential[i], potential[i+1],
                       force[i-1], force[i], force[i+1]);
            }
            EXPECT_TRUE(relerror < toler);
        }
    }
}

/*! \brief Utility to do the real testing
 *
 * \param[in] cd    Charge distribution type
 * \param[in] irow  Row number for atom i
 * \param[in] jrow  Row number for atom j
 * \param[in] izeta Distribution width atom i (may be 0)
 * \param[in] jzeta Distribution width atom j (may be 0)
 * \param[in] checker The checker data structure
 */
void testCoulomb(alexandria::ChargeType           cd,
                 int                              irow,
                 int                              jrow,
                 double                           izeta,
                 double                           jzeta,
                 gmx::test::TestReferenceChecker *checker)
{
    std::vector<double> coulomb;
    std::vector<double> force;
    std::vector<double> ncoulomb;
    std::vector<double> nforce;
    double dx = 0.002;
    for(int i = 0; i <= 50; i++)
    {
        double r = dx*i;
        
        switch (cd)
        {
        case alexandria::ChargeType::Gaussian:
            coulomb.push_back(Coulomb_GG(r,  izeta, jzeta));
            force.push_back(-DCoulomb_GG(r,  izeta, jzeta));
            ncoulomb.push_back(Nuclear_GG(r, izeta));
            nforce.push_back(-DNuclear_GG(r, izeta));
            break;
        case alexandria::ChargeType::Slater:
            coulomb.push_back(Coulomb_SS(r,  irow, jrow, izeta, jzeta));
            force.push_back(-DCoulomb_SS(r,  irow, jrow, izeta, jzeta));
            ncoulomb.push_back(Nuclear_SS(r, irow, izeta));
            nforce.push_back(-DNuclear_SS(r, irow, izeta));
            break;
        default:
            break;
        }
    }
    auto name = alexandria::chargeTypeName(cd).c_str();
    char buf[256];
    if (cd == alexandria::ChargeType::Slater)
    {
        checker->checkInt64(irow, "irow");
        checker->checkInt64(jrow, "jrow");
    }
    checker->checkDouble(izeta, "izeta");
    checker->checkDouble(jzeta, "jzeta");
    snprintf(buf, sizeof(buf), "Potential-%s", name);
    checker->checkSequence(coulomb.begin(), coulomb.end(), buf);
    snprintf(buf, sizeof(buf), "Force-%s", name);
    checker->checkSequence(force.begin(), force.end(), buf);
    snprintf(buf, sizeof(buf), "NuclearPotential-%s", name);
    checker->checkSequence(ncoulomb.begin(), ncoulomb.end(), buf);
    snprintf(buf, sizeof(buf), "NuclearForce-%s", name);
    checker->checkSequence(nforce.begin(), nforce.end(), buf);
    check_force(coulomb, force, dx);
    check_force(ncoulomb, nforce, dx);
}

class SlaterTest : public ::testing::TestWithParam<std::tuple<std::tuple<int, int>, std::tuple<double, double> > >
{
    protected:
        int                        irow_;
        int                        jrow_;
        double                     xi_;
        double                     xj_;
        test::TestReferenceData    refData_;
        test::TestReferenceChecker checker_;
        
        SlaterTest () : checker_(refData_.rootChecker())
        {
            double toler     = 2e-2;
            auto   tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, toler);
            checker_.setDefaultTolerance(tolerance);
            irow_ = std::get<0>(std::get<0>(GetParam()));
            jrow_ = std::get<1>(std::get<0>(GetParam()));
            xi_   = std::get<0>(std::get<1>(GetParam()));
            xj_   = std::get<1>(std::get<1>(GetParam()));
        }
        void runTest()
        {
            testCoulomb(alexandria::ChargeType::Slater, irow_, jrow_, xi_, xj_, &checker_);
        }
};

class GaussianTest : public ::testing::TestWithParam<std::tuple<double, double> >
{
    protected:
        double                     xi_;
        double                     xj_;
        test::TestReferenceData    refData_;
        test::TestReferenceChecker checker_;
        
        GaussianTest () : checker_(refData_.rootChecker())
        {
            double toler     = 1e-5;
            auto   tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, toler);
            checker_.setDefaultTolerance(tolerance);
            xi_   = std::get<0>(GetParam());
            xj_   = std::get<1>(GetParam());
        }
        void runTest()
        {
            testCoulomb(alexandria::ChargeType::Gaussian, 0, 0, xi_, xj_, &checker_);
        }

};

TEST_P (GaussianTest, All)
{
    runTest();
}

TEST_P (SlaterTest, All)
{
    runTest();
}

//! Rows for Slater tests
const std::vector<std::tuple<int, int> > make_rows()
{
#if !HAVE_LIBCLN
    return {
        { 1, 1 }
    };
#else
    return {
        { 1, 1 }, { 1, 2 }, { 2, 2 },
        { 1, 3 }, { 2, 3 }, { 3, 3 }
    };
#endif
};

//! Instantiate the combination of Slater row-row interactions to test
static const std::vector<std::tuple<int, int> > c_rows = make_rows();

//! xi and xj for tests
const std::vector<std::tuple<double, double> > make_xi()
{
    return {
        {  5.6,   5.7  },
        {  5.7,   5.71 },
        {  5.91,  5.9  },
        { 15.8,  16.0  },
        {  6.1,   6.6  },
        { 22.3,  22.4  },
        { 34.6,  34.5  }
    };
};

//! Instantiate the combination of Slater widths to test
static auto c_xi = make_xi();

INSTANTIATE_TEST_CASE_P(Xi, SlaterTest, ::testing::Combine(::testing::ValuesIn(c_rows), ::testing::ValuesIn(c_xi)));

INSTANTIATE_TEST_CASE_P(Xi, GaussianTest, ::testing::ValuesIn(c_xi));

//! integer xi and xj for tests
const std::vector<std::tuple<double, double> > make_xiInteger()
{
    return {
        {  3.0,  4.0 },
        { 17.0, 18.0 },
        { 25.0, 26.0 },
        { 29.0, 28.0 },
        { 30.0, 29.0 },
        { 31.0, 33.0 },
        { 37.0, 38.0 }
    };
};

//! Instantiate xi and xj pair values
static auto c_xiInteger = make_xiInteger();

INSTANTIATE_TEST_CASE_P(IntegerXi, SlaterTest, ::testing::Combine(::testing::ValuesIn(c_rows), ::testing::ValuesIn(c_xiInteger)));

//! integer xi and xj for tests
const std::vector<std::tuple<double, double> > make_xiIdentical()
{
    return {
        {  3.0,  3.0 },
        { 17.0, 17.0 },
        { 26.1, 26.1 },
        {  9.4,  9.4 },
        { 10.0, 10.0 }
    };
};

static auto c_xiIdentical = make_xiIdentical();

INSTANTIATE_TEST_CASE_P(IdenticalXi, SlaterTest, ::testing::Combine(::testing::ValuesIn(c_rows), ::testing::ValuesIn(c_xiIdentical)));

} // namespace

} // namespace gmx

