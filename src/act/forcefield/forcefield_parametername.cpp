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

#include "forcefield_parametername.h"

namespace alexandria
{

const char *lj12_6_name[lj12_6NR] = { "sigma", "epsilon", "sigma_ij", "epsilon_ij" };

const char *lj8_6_name[lj8_6NR] = { "sigma", "epsilon", "sigma_ij", "epsilon_ij" };

const char *lj14_7_name[lj14_7NR] = { "sigma", "epsilon", "gamma", "delta",
                                      "sigma_ij", "epsilon_ij", "gamma_ij", "delta_ij" };

const char *qt_name[qtNR] = { "aqt", "bqt", "aqt_ij", "bqt_ij" };

const char *wbh_name[wbhNR] = { "sigma", "epsilon", "gamma", "sigma_ij", "epsilon_ij", "gamma_ij" };

const char *gbh_name[gbhNR] = { "rmin", "epsilon", "gamma", "delta", "rmin_ij", "epsilon_ij", "gamma_ij", "delta_ij" };

const char *coul_name[coulNR] = { "zeta", "zeta_i", "zeta_j" };

const char *bond_name[bondNR] = { "kb", "bondlength", "bondenergy" };

const char *cubic_name[cubicNR] = { "bondlength", "rmax", "kb", "De" };

const char *angle_name[angleNR] = { "kt", "angle" };

const char *ub_name[ubNR] = { "kt", "ub_angle", "r13", "kub" };

const char *ps_names[psNR] = { "ps_angle", "rij0", "rjk0" };

const char *pol_name[polNR] = { "alpha", "kshell" };

const char *morse_name[morseNR] = { "beta", "De", "D0", "bondlength" };

const char *linang_name[linangNR] = { "a", "klin" };

const char *idih_name[idihNR] = { "kimp" };

const char *fdih_name[fdihNR] = { "c0", "c1", "c2", "c3", "c4", "c5", "c6" };

const char *pdih_name[pdihNR] = { "phi0", "kp", "mult" };

const char *vsite1_name[vsite1NR] = { "vs1a" };

const char *vsite2_name[vsite2NR] = { "vs2a" };

const char *vsite2fd_name[vsite2fdNR] = { "vs2fd_a" };

const char *vsite3_name[vsite3NR] = { "vs3a", "vs3b" };

const char *vsite3fd_name[vsite3fdNR] = { "vs3fd_a", "vs3fd_b" };

const char *vsite3fad_name[vsite3fadNR] = { "vs3fad_a", "vs3fad_b" };

const char *vsite3out_name[vsite3outNR] = { "vs3out_a", "vs3out_b", "vs3out_c" };

const char *vsite3outs_name[vsite3outsNR] = { "vs3outs_a", "vs3outs_c" };

} // namespace alexandria
