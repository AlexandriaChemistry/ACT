/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2025
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

#ifndef FORCEFIELDPARAMETERNAME_H
#define FORCEFIELDPARAMETERNAME_H

namespace alexandria
{

    enum { lj12_6SIGMA = 0, lj12_6EPSILON = 1, lj12_6NR = 2 };

    enum { lj8_6SIGMA = 0, lj8_6EPSILON = 1, lj8_6NR = 2 };

    enum { wbhSIGMA = 0, wbhEPSILON = 1, wbhGAMMA = 2, wbhNR = 3 };

    enum { bhA = 0, bhB = 1, bhC6 = 2, bhNR = 3 };

    enum { ttA, ttB, ttC6, ttC8, ttC10, ttNR };

    enum { tt2bA, tt2bBexch, tt2bBdisp, tt2bC6, tt2bC8, tt2bC10, tt2bNR };

    enum { lj14_7SIGMA = 0, lj14_7EPSILON = 1, lj14_7GAMMA = 2, lj14_7DELTA = 3,
           lj14_7NR = 4 };

    enum { expA = 0, expB = 1, expNR = 2 };

    enum { dexpA1 = 0, dexpA2 = 1, dexpB = 2, dexpNR = 3 };

    enum { gbhRMIN = 0, gbhEPSILON = 1, gbhGAMMA = 2, gbhDELTA = 3, gbhNR = 4 };

    enum { coulZETA, coulZETA2, coulNR };

    enum { bondKB = 0, bondLENGTH = 1, bondENERGY = 2, bondNR = 3 };

    enum { cubicLENGTH = 0, cubicRMAX = 1, cubicKB = 2, cubicDE = 3, cubicNR = 4 };

    enum { angleKT = 0, angleANGLE = 1, angleNR = 2 };

    enum { ubKT = 0, ubANGLE = 1, ubR13 = 2, ubKUB = 3, ubNR = 4 };

    enum { psANGLE = 0, psRIJ0 = 1, psRJK0 = 2, psNR = 3 };

    enum { polALPHA = 0, polRHYPER = 1, polFCHYPER = 2, polNR = 3 };

    enum { morseBETA = 0, morseDE = 1, morseD0 = 2, morseLENGTH = 3, morseNR = 4 };

    enum { huaLENGTH = 0, huaDE = 1, huaB = 2, huaC = 3, huaNR = 4 };

    enum { linangA = 0, linangKLIN = 1, linangNR = 2 };

    enum { idihKPHI = 0, idihNR = 1 };

    enum { fdihC0 = 0, fdihC1 = 1, fdihC2 = 2, fdihC3 = 3, fdihC4 = 4, fdihC5 = 5, fdihNR = 6 };

    enum { pdihANGLE = 0, pdihKP = 1, pdihMULT = 2, pdihNR = 3 };

    enum { vsite1A = 0, vsite1NR = 1 };

    enum { vsite2A = 0, vsite2NR = 1 };

    enum { vsite2fdA = 0, vsite2fdNR = 1 };

    enum { vsite3A = 0,vsite3B = 1, vsite3NR = 2 };

    enum { vsite3sA = 0, vsite3sNR = 1 };

    enum { vsite3fdA = 0, vsite3fdB = 1, vsite3fdNR = 2 };

    enum { vsite3fadA = 0, vsite3fadB = 1, vsite3fadNR = 2 };

    enum { vsite3outA = 0, vsite3outB = 1, vsite3outC = 2,  vsite3outNR = 3 };

    enum { vsite3outsA = 0, vsite3outsC = 1,  vsite3outsNR = 2 };

} // namespace alexandria

#endif
