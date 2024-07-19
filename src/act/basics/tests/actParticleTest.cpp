/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2024
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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include "../act_particle.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

TEST(ActParticleTest, stringToActParticle) {
    bool        b;
    ActParticle apType;
    
    b = stringToActParticle("Atom", &apType);
    EXPECT_TRUE(b && apType == ActParticle::Atom);
    b = stringToActParticle("VSite", &apType);
    EXPECT_TRUE(b && apType == ActParticle::Vsite);
    b = stringToActParticle("Shell", &apType);
    EXPECT_TRUE(b && apType == ActParticle::Shell);
    b = stringToActParticle("SigmaHole", &apType);
    EXPECT_TRUE(b && apType == ActParticle::SigmaHole);

    EXPECT_FALSE(stringToActParticle("Atoms", &apType));
    EXPECT_FALSE(stringToActParticle("Vsite", &apType));
    EXPECT_FALSE(stringToActParticle("Sunflower", &apType));
}

}

}
