/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2020
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
#include <gtest/gtest.h>

#include "alexandria/units.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
//#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class UnitsTest : public gmx::test::CommandLineTestBase
{
    protected:
    gmx::test::TestReferenceChecker checker_;
    std::vector<std::string> units_ = { "Angstrom",  "nm", "pm", "Bohr",      
        "kcal/mol", "kJ/mol", "J/mol K", "cal/mol K", "Hartree",   
        "Hartree/e", "Angstrom3", "Coulomb", "Debye", "Electron", "Buckingham",
        "1/nm", "degree", "kJ/mol/rad2", "kJ/mol/nm2" };
    UnitsTest () : checker_(this->rootChecker())
    {
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    static void TearDownTestCase()
    {
    }
    
    void toTest(double value)
    {
        for (const auto &unit : units_)
        {
            double x = convertToGromacs(value, unit);
            checker_.checkDouble(x, unit.c_str());
        }
    }
    void fromTest(double value)
    {
        for (const auto &unit : units_)
        {
            double y = convertFromGromacs(value, unit);
            checker_.checkDouble(y, unit.c_str());
        }
    }
};

TEST_F(UnitsTest, ToGromacs) 
{
    toTest(3.0);
}

TEST_F(UnitsTest, FromGromacs) 
{
    fromTest(1.23e-5);
}

#if CMAKE_BUILD_TYPE == CMAKE_BUILD_TYPE_DEBUG
TEST(UnitsTestSimple, WrongUnitThrows)
{
    double x, y;
    EXPECT_THROW(x = convertToGromacs(1.0, "Foo"), gmx::InternalError);
    EXPECT_THROW(y = convertFromGromacs(2.2, "Bar"), gmx::InternalError);
}
#endif
}

} // namespace
