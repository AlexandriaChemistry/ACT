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
 * Unit tests for DevComputer classes.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <gtest/gtest.h>
#include <map>

#include "act/alexandria/devcomputer.h"

namespace alexandria
{
namespace
{

// Test fixture for ForceEnergyDevComputer
class ForceEnergyDevComputerTest : public ::testing::Test
{
protected:
    ForceEnergyDevComputerTest() = default;
};

// Test: Constructor initializes with correct name
TEST_F(ForceEnergyDevComputerTest, ConstructorSetsCorrectName)
{
    std::map<eRMS, double> boltzmann;
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Constructor with empty Boltzmann map
TEST_F(ForceEnergyDevComputerTest, ConstructorWithEmptyBoltzmannMap)
{
    std::map<eRMS, double> boltzmann;
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_FALSE(computer.separateInductionCorrection());
}

// Test: Constructor with non-empty Boltzmann map
TEST_F(ForceEnergyDevComputerTest, ConstructorWithNonEmptyBoltzmannMap)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 298.15;
    boltzmann[eRMS::Force2] = 0.0;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
    EXPECT_FALSE(computer.separateInductionCorrection());
}

// Test: setSeparateInductionCorrection and separateInductionCorrection methods
TEST_F(ForceEnergyDevComputerTest, SetAndGetSeparateInductionCorrection)
{
    std::map<eRMS, double> boltzmann;
    ForceEnergyDevComputer computer(boltzmann);
    
    // Default should be true (from header initialization)
    EXPECT_TRUE(computer.separateInductionCorrection());
    
    // Set to false
    computer.setSeparateInductionCorrection(false);
    EXPECT_FALSE(computer.separateInductionCorrection());
    
    // Set back to true
    computer.setSeparateInductionCorrection(true);
    EXPECT_TRUE(computer.separateInductionCorrection());
}

// Test: Multiple flag toggles
TEST_F(ForceEnergyDevComputerTest, MultipleSeparateInductionCorrectionToggles)
{
    std::map<eRMS, double> boltzmann;
    ForceEnergyDevComputer computer(boltzmann);
    
    for (int i = 0; i < 5; ++i)
    {
        computer.setSeparateInductionCorrection(false);
        EXPECT_FALSE(computer.separateInductionCorrection());
        
        computer.setSeparateInductionCorrection(true);
        EXPECT_TRUE(computer.separateInductionCorrection());
    }
}

// Test: Boltzmann temperature with EPOT
TEST_F(ForceEnergyDevComputerTest, BoltzmannTemperatureEPOT)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 298.15;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Boltzmann temperature with Force2
TEST_F(ForceEnergyDevComputerTest, BoltzmannTemperatureForce2)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::Force2] = 300.0;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Boltzmann temperature with multiple eRMS terms
TEST_F(ForceEnergyDevComputerTest, BoltzmannTemperatureMultipleTerms)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 298.15;
    boltzmann[eRMS::Force2] = 0.0;
    boltzmann[eRMS::Interaction] = 500.0;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
    EXPECT_TRUE(computer.separateInductionCorrection());
}

// Test: Inheritance from DevComputer
TEST_F(ForceEnergyDevComputerTest, InheritsFromDevComputer)
{
    std::map<eRMS, double> boltzmann;
    ForceEnergyDevComputer computer(boltzmann);
    
    // Verify that DevComputer interface is accessible
    const std::string& name = computer.name();
    EXPECT_EQ("ForceEnergy", name);
}

// Test: Zero Boltzmann temperature
TEST_F(ForceEnergyDevComputerTest, ZeroBoltzmannTemperature)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 0.0;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Large Boltzmann temperature values
TEST_F(ForceEnergyDevComputerTest, LargeBoltzmannTemperatureValues)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = 1e6;
    boltzmann[eRMS::Force2] = 1e10;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

// Test: Negative Boltzmann temperature (edge case)
TEST_F(ForceEnergyDevComputerTest, NegativeBoltzmannTemperature)
{
    std::map<eRMS, double> boltzmann;
    boltzmann[eRMS::EPOT] = -100.0;
    
    ForceEnergyDevComputer computer(boltzmann);
    EXPECT_EQ("ForceEnergy", computer.name());
}

} // namespace
} // namespace alexandria
