/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2008-2022
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

#ifndef _SLATER_LOW_H
#define _SLATER_LOW_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_LIBCLN
// CLN uses some code that is not allowed by modern compilers
// with the -Werror flag
#pragma clang diagnostic ignored "-Wcast-function-type-mismatch"
#pragma clang diagnostic ignored "-Wmismatched-tags"
#include <cln/cln.h>
#pragma clang diagnostic pop
#pragma clang diagnostic pop

using namespace cln;
#define PRECISION 80

static cl_R           ZERO      = "0.0_80";
static cl_R           ONE       = "1.0_80";
static cl_R           TWO       = "2.0_80";
static cl_R           THREE     = "3.0_80";
static cl_R           FOUR      = "4.0_80";
static cl_R           FIVE      = "5.0_80";
static cl_R           SIX       = "6.0_80";
static cl_R           SEVEN     = "7.0_80";
static cl_R           EIGHT     = "8.0_80";
static cl_R           NINE      = "9.0_80";
static float_format_t precision = float_format(80);

extern cl_R Power(cl_R a, int b);

extern cl_R Slater_1S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Nuclear_1S(cl_R r, cl_R xi);

extern cl_R Nuclear_2S(cl_R r, cl_R xi);

extern cl_R Nuclear_3S(cl_R r, cl_R xi);

extern cl_R Nuclear_4S(cl_R r, cl_R xi);

extern cl_R Nuclear_5S(cl_R r, cl_R xi);

extern cl_R Nuclear_6S(cl_R r, cl_R xi);

extern cl_R DNuclear_1S(cl_R r, cl_R xi);

extern cl_R DNuclear_2S(cl_R r, cl_R xi);

extern cl_R DNuclear_3S(cl_R r, cl_R xi);

extern cl_R DNuclear_4S(cl_R r, cl_R xi);

extern cl_R DNuclear_5S(cl_R r, cl_R xi);

extern cl_R DNuclear_6S(cl_R r, cl_R xi);

#else
#include <cmath>

extern double Slater_1S_1S(double r, double xi, double xj);

extern double Slater_1S_2S(double r, double xi, double xj);

extern double Slater_2S_1S(double r, double xi, double xj);

extern double Slater_1S_3S(double r, double xi, double xj);

extern double Slater_3S_1S(double r, double xi, double xj);

extern double Slater_2S_2S(double r, double xi, double xj);

extern double Slater_2S_3S(double r, double xi, double xj);

extern double Slater_3S_2S(double r, double xi, double xj);

extern double Slater_3S_3S(double r, double xi, double xj);

extern double DSlater_1S_1S(double r, double xi, double xj);

extern double DSlater_1S_2S(double r, double xi, double xj);

extern double DSlater_2S_1S(double r, double xi, double xj);

extern double DSlater_1S_3S(double r, double xi, double xj);

extern double DSlater_3S_1S(double r, double xi, double xj);

extern double DSlater_2S_2S(double r, double xi, double xj);

extern double DSlater_2S_3S(double r, double xi, double xj);

extern double DSlater_3S_2S(double r, double xi, double xj);

extern double DSlater_3S_3S(double r, double xi, double xj);

extern double Nuclear_1S(double r, double xi);

extern double Nuclear_2S(double r, double xi);

extern double Nuclear_3S(double r, double xi);

extern double DNuclear_1S(double r, double xi);

extern double DNuclear_2S(double r, double xi);

extern double DNuclear_3S(double r, double xi);

#endif
#endif
