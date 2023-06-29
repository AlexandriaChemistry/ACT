/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023
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
 * \ingroup module_listed-forces
 */
#include "actpre.h"

#include "act/forces/forcecomputer.h"

#include <cmath>

#include <memory>

#include <gtest/gtest.h>

#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield_utils.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

TEST(Vsite2, HF)
{
    std::string forcefield("ACS-pg-vs2");
    // Get forcefield
    auto pd  = getForceField(forcefield);
    
    auto fh = Identifier({ "f_b", "h_b", "v2" }, { 1, 9 }, CanSwap::Vsite2);            
    auto hf = Identifier({ "h_b", "f_b", "v2" }, { 1, 9 }, CanSwap::Vsite2);        
    // Compare identifiers
    EXPECT_TRUE(fh == hf);
    auto itype = InteractionType::VSITE2;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(fh));
    auto param_fh = fs.findParameterType(fh, vsite2_name[vsite2A]);
    EXPECT_TRUE(param_fh->internalValue() == 1.05);
    if (fs.parameterExists(hf))
    {
        auto param_hf = fs.findParameterType(hf, vsite2_name[vsite2A]);
        EXPECT_TRUE(param_fh->internalValue() == param_hf->internalValue());
    }
}

TEST(Vsite2, HFCanSwapNo)
{
    std::string forcefield("ACS-pg-vs2");
    // Get forcefield
    auto pd  = getForceField(forcefield);
    
    auto fh = Identifier({ "f_b", "h_b", "v2" }, { 1, 9 }, CanSwap::No); 
    auto hf = Identifier({ "h_b", "f_b", "v2" }, { 1, 9 }, CanSwap::No);
    // Compare identifiers
    EXPECT_FALSE(fh == hf);
    auto itype = InteractionType::VSITE2;
    EXPECT_TRUE(pd->interactionPresent(itype));
    auto fs = pd->findForcesConst(itype);
    // Compare values
    EXPECT_TRUE(fs.parameterExists(fh));
    auto param_fh = fs.findParameterType(fh, vsite2_name[vsite2A]);
    EXPECT_TRUE(param_fh->internalValue() == 1.05);
    EXPECT_FALSE(fs.parameterExists(hf));
}

}  // namespace

}  // namespace gmx
