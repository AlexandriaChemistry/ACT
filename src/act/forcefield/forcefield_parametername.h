/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024
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

    enum { lj12_6SIGMA = 0, lj12_6EPSILON = 1, lj12_6SIGMA_IJ = 2, lj12_6EPSILON_IJ = 3, lj12_6NR = 4 };

    extern const char *lj12_6_name[lj12_6NR];

    enum { lj8_6SIGMA = 0, lj8_6EPSILON = 1, lj8_6SIGMA_IJ = 2, lj8_6EPSILON_IJ = 3, lj8_6NR = 4 };

    extern const char *lj8_6_name[lj8_6NR];

    enum { wbhSIGMA = 0, wbhEPSILON = 1, wbhGAMMA = 2, wbhSIGMA_IJ = 3, wbhEPSILON_IJ = 4, wbhGAMMA_IJ = 5, wbhNR = 6 };

    extern const char *wbh_name[wbhNR];

    enum { bhA = 0, bhB = 1, bhC6 = 2, bhA_IJ = 3, bhB_IJ = 4, bhC6_IJ = 5, bhNR = 6 };

    extern const char *bh_name[bhNR];

    enum { ttA = 0, ttB = 1, ttC6 = 2, ttC8 = 3, ttC10 = 4, ttA_IJ = 5, ttB_IJ = 6, ttC6_IJ = 7, ttC8_IJ = 8, ttC10_IJ = 9, ttNR = 10 };

    extern const char *tt_name[ttNR];

    enum { lj14_7SIGMA = 0, lj14_7EPSILON = 1, lj14_7GAMMA = 2, lj14_7DELTA = 3,
           lj14_7SIGMA_IJ = 4, lj14_7EPSILON_IJ = 5, lj14_7GAMMA_IJ = 6, lj14_7DELTA_IJ = 7, lj14_7NR = 8 };

    extern const char *lj14_7_name[lj14_7NR];

    enum { expA = 0, expB = 1, expA_IJ = 2, expB_IJ = 3, expNR = 4 };

    extern const char *exp_name[expNR];

    enum { dexpA1 = 0, dexpA2 = 1, dexpB = 2, dexpA1_IJ = 3, dexpA2_IJ = 4, dexpB_IJ = 5, dexpNR = 6 };

    extern const char *dexp_name[dexpNR];

    enum { gbhRMIN = 0, gbhEPSILON = 1, gbhGAMMA = 2, gbhDELTA = 3, gbhRMIN_IJ = 4, gbhEPSILON_IJ = 5, gbhGAMMA_IJ = 6, gbhDELTA_IJ = 7, gbhNR = 8 };

    extern const char *gbh_name[gbhNR];

    enum { coulZETA = 0, coulZETAI = 1, coulZETAJ = 2, coulNR = 3 };

    extern const char *coul_name[coulNR];

    enum { bondKB = 0, bondLENGTH = 1, bondENERGY = 2, bondNR = 3 };

    extern const char *bond_name[bondNR];

    enum { cubicLENGTH = 0, cubicRMAX = 1, cubicKB = 2, cubicDE = 3, cubicNR = 4 };

    extern const char *cubic_name[cubicNR];

    enum { angleKT = 0, angleANGLE = 1, angleNR = 2 };

    extern const char *angle_name[angleNR];

    enum { ubKT = 0, ubANGLE = 1, ubR13 = 2, ubKUB = 3, ubNR = 4 };

    extern const char *ub_name[ubNR];

    enum { psANGLE = 0, psRIJ0 = 1, psRJK0 = 2, psNR = 3 };

    extern const char *ps_names[psNR];

    enum { polALPHA = 0, polRHYPER = 1, polFCHYPER = 2, polNR = 3 };

    extern const char *pol_name[polNR];

    enum { morseBETA = 0, morseDE = 1, morseD0 = 2, morseLENGTH = 3, morseNR = 4 };

    extern const char *morse_name[morseNR];

    enum { huaLENGTH = 0, huaDE = 1, huaB = 2, huaC = 3, huaNR = 4 };

    extern const char *hua_name[huaNR];

    enum { linangA = 0, linangKLIN = 1, linangNR = 2 };

    extern const char *linang_name[linangNR];

    enum { idihKPHI = 0, idihNR = 1 };

    extern const char *idih_name[idihNR];

    enum { fdihC0 = 0, fdihC1 = 1, fdihC2 = 2, fdihC3 = 3, fdihC4 = 4, fdihC5 = 5, fdihNR = 6 };

    extern const char *fdih_name[fdihNR];

    enum { pdihANGLE = 0, pdihKP = 1, pdihMULT = 2, pdihNR = 3 };

    extern const char *pdih_name[pdihNR];

    enum { vsite1A = 0, vsite1NR = 1 };

    extern const char *vsite1_name[vsite1NR];

    enum { vsite2A = 0, vsite2NR = 1 };

    extern const char *vsite2_name[vsite2NR];

    enum { vsite2fdA = 0, vsite2fdNR = 1 };

    extern const char *vsite2fd_name[vsite2fdNR];

    enum { vsite3A = 0,vsite3B = 1, vsite3NR = 2 };

    extern const char *vsite3_name[vsite3NR];

    enum { vsite3sA = 0, vsite3sNR = 1 };

    extern const char *vsite3s_name[vsite3sNR];

    enum { vsite3fdA = 0, vsite3fdB = 1, vsite3fdNR = 2 };

    extern const char *vsite3fd_name[vsite3fdNR];

    enum { vsite3fadA = 0, vsite3fadB = 1, vsite3fadNR = 2 };

    extern const char *vsite3fad_name[vsite3fadNR];

    enum { vsite3outA = 0, vsite3outB = 1, vsite3outC = 2,  vsite3outNR = 3 };

    extern const char *vsite3out_name[vsite3outNR];

    enum { vsite3outsA = 0, vsite3outsC = 1,  vsite3outsNR = 2 };

    extern const char *vsite3outs_name[vsite3outsNR];

} // namespace alexandria

#endif
