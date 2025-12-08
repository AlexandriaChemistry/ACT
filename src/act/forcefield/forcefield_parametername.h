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

    enum { lj12_6SIGMA, lj12_6EPSILON, lj12_6NR };
    
    enum { lj8_6SIGMA, lj8_6EPSILON, lj8_6NR };
    
    enum { lj12_6_4SIGMA, lj12_6_4EPSILON, lj12_6_4GAMMA, lj12_6_4NR };

    enum { wbhSIGMA, wbhEPSILON, wbhGAMMA, wbhNR };

    enum { bhA, bhB, bhC6, bhNR };

    enum { ttA, ttB, ttC6, ttC8, ttC10, ttNR };

    enum { SIttA, SIttBexch, SIttBdisp, SIttC6, SIttC8, SIttC10, SIttNR };

    enum { tt2bA, tt2bBexch, tt2bBdisp, tt2bC6, tt2bC8, tt2bC10, tt2bNR };

    enum { lj14_7SIGMA, lj14_7EPSILON, lj14_7GAMMA, lj14_7DELTA,
           lj14_7NR };

    enum { expA, expB, expVSite, expNR };

    enum { dexpA1, dexpA2, dexpB, dexpNR };

    enum { gbhRMIN, gbhEPSILON, gbhGAMMA, gbhDELTA, gbhNR };

    enum { coulZETA, coulZETA2, coulNR };

    enum { fbprK, fbprR0, fbprNR };

    enum { bondKB, bondLENGTH, bondENERGY, bondNR };

    enum { cubicLENGTH, cubicRMAX, cubicKB, cubicDE, cubicNR };

    enum { angleKT, angleANGLE, angleNR };

    enum { ubKT, ubANGLE, ubR13, ubKUB, ubNR };

    enum { psANGLE, psRIJ0, psRJK0, psNR };

    enum { polALPHA, polRHYPER, polFCHYPER, polNR };

    enum { morseBETA, morseDE, morseD0, morseLENGTH, morseNR };

    enum { huaLENGTH, huaDE, huaB, huaC, huaNR };

    enum { linangA, linangKLIN, linangNR };

    enum { idihKPHI, idihNR };

    enum { fdihC0, fdihC1, fdihC2, fdihC3, fdihC4, fdihC5, fdihNR };

    enum { pdihANGLE, pdihKP, pdihMULT, pdihNR };

    enum { vsite1A, vsite1NR };

    enum { vsite2A, vsite2NR };

    enum { vsite2fdA, vsite2fdNR };

    enum { vsite3A,vsite3B, vsite3NR };

    enum { vsite3sA, vsite3sNR };

    enum { vsite3fdA, vsite3fdB, vsite3fdNR };

    enum { vsite3fadA, vsite3fadB, vsite3fadNR };

    enum { vsite3outA, vsite3outB, vsite3outC,  vsite3outNR };

    enum { vsite3outsA, vsite3outsC,  vsite3outsNR };

    enum { vsite4A, vsite4B, vsite4C, vsite4NR };

    enum { vsite4sA, vsite4sB, vsite4sNR };

    enum { vsite4s3A, vsite4s3NR };

} // namespace alexandria

#endif
