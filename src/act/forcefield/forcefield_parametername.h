/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

    enum { ljSIGMA = 0, ljEPSILON = 1, ljC6_IJ = 2, ljC12_IJ = 3, ljNR = 4 };
    
    extern const char *lj_name[ljNR];

    enum { lj_86SIGMA = 0, lj_86EPSILON = 1, lj_C6_IJ = 2, lj_C8_IJ = 3, lj_NR_86 = 4 };

    extern const char *lj_86_name[lj_NR_86];
    
    enum { wbhSIGMA = 0, wbhEPSILON = 1, wbhGAMMA = 2, wbhSIGMA_IJ = 3, wbhEPSILON_IJ = 4, wbhGAMMA_IJ = 5, wbhNR = 6 };
    
    extern const char *wbh_name[wbhNR];

    enum { lj_147SIGMA = 0, lj_147EPSILON = 1, lj_147GAMMA = 2, lj_147DELTA = 3, lj_147SIGMA_IJ = 4, lj_147EPSILON_IJ = 5, lj_147GAMMA_IJ = 6, lj_147DELTA_IJ = 7, lj_NR_147 = 8 };

    extern const char *lj_147_name[lj_NR_147];
    
    enum { gbhRMIN = 0, gbhEPSILON = 1, gbhGAMMA = 2, gbhDELTA = 3, gbhRMIN_IJ = 4, gbhEPSILON_IJ = 5, gbhGAMMA_IJ = 6, gbhDELTA_IJ = 7, gbhNR = 8 };
    
    extern const char *gbh_name[gbhNR];
    
    enum { coulZETAI = 0, coulZETAJ = 1, coulNR = 2 };
    
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
    
    enum { polALPHA = 0, polKSH = 1, polNR = 2 };
    
    extern const char *pol_name[polNR];
    
    enum { morseBETA = 0, morseDE = 1, morseD0 = 2, morseLENGTH = 3, morseNR = 4 };
    
    extern const char *morse_name[morseNR];
    
    enum { linangA = 0, linangKLIN = 1, linangNR = 2 };
    
    extern const char *linang_name[linangNR];
    
    enum { idihKPHI = 0, idihNR = 1 };
    
    extern const char *idih_name[idihNR];
    
    enum { fdihC0 = 0, fdihC1 = 1, fdihC2 = 2, fdihC3 = 3, fdihC4 = 4, fdihC5 = 5, fdihC6 = 6, fdihNR = 7 };
    
    extern const char *fdih_name[fdihNR];
    
    enum { pdihANGLE = 0, pdihKP = 1, pdihMULT = 2, pdihNR = 3 };
    
    extern const char *pdih_name[pdihNR];
    
} // namespace alexandria

#endif
