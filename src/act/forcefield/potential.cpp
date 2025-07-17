/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024,2025
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
#include "potential.h"

#include <cstring>
#include <map>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

std::map<Potential, PotentialProperties> potprops = {
    { Potential::NONE, { "", -1, { }, "0" } },
    { Potential::LJ12_6, { "LJ12_6" , F_LJ, { "sigma", "epsilon" }, "0" } },
    { Potential::LJ8_6, {  "LJ8_6", -1, { "sigma", "epsilon" }, "0" } },
    { Potential::LJ14_7, { "LJ14_7", -1, { "sigma", "epsilon", "gamma", "delta" }, "0" } },
    { Potential::EXPONENTIAL,
      { "EXPONENTIAL", -1, { "aexp", "bexp", "vsite" },
        // Implement xor function for Kronecker combination rule
        "-A*exp(-b*r); b=0.5*(bexp1+bexp2); A=select(vsite1*vsite2+(1-vsite1)*(1-vsite2),0,sqrt(aexp1*aexp2))" } },
    { Potential::MACDANIEL_SCHMIDT,
      { "MACDANIEL_SCHMIDT", -1, { "a1dexp", "a2dexp", "bdexp" },
        "(Ab-Aa)*exp(-b*r); b=0.5*(bdexp1+bdexp2); Aa=sqrt(a1dexp1*a1dexp2); Ab=sqrt(a2dexp1*a2dexp2)" } },
    { Potential::WANG_BUCKINGHAM, { "WANG_BUCKINGHAM", F_WBHAM, { "sigma", "epsilon", "gamma" }, "0" } },
    { Potential::BUCKINGHAM, { "BUCKINGHAM", -1, { "Abh", "bbh", "c6bh" }, "0" } },
    { Potential::TANG_TOENNIES, { "TANG_TOENNIES", -1, { "Att", "btt", "c6tt", "c8tt", "c10tt" }, "0" } },
    { Potential::TT2b, { "TT2b", -1, { "Att2b", "bExchtt2b", "bDisptt2b", "c6tt2b", "c8tt2b", "c10tt2b" }, "0" } },
    { Potential::GENERALIZED_BUCKINGHAM, { "GENERALIZED_BUCKINGHAM", F_GBHAM, { "rmin", "epsilon", "gamma", "delta" }, "0" } },
    { Potential::COULOMB_GAUSSIAN, { "COULOMB_GAUSSIAN", -1, { "zeta", "zeta2" }, "0" } },
    { Potential::COULOMB_POINT, { "COULOMB_POINT", F_COUL_SR, { "zeta", "zeta2" }, "0" } },
    { Potential::COULOMB_SLATER, { "COULOMB_SLATER", -1, { "zeta", "zeta2" }, "0" } },
    { Potential::HARMONIC_BONDS, { "HARMONIC_BONDS", F_BONDS, { "kb", "bondlength", "bondenergy" }, "0" } },
    { Potential::CUBIC_BONDS,
      { "CUBIC_BONDS", F_CUBICBONDS, { "bondlength", "rmax", "kb", "De" },
        "(kb*(r-bondlength)^2 * (rmax-r) - D_e)*step(2*rmax+bondlength-3*r) + step(3*r-2*rmax-bondlength)*(-D_e - (4*kb/27)*(bondlength-rmax)^2);"} },
    { Potential::HARMONIC_ANGLES, { "HARMONIC_ANGLES", F_ANGLES, { "kt", "angle" }, "0" } },
    { Potential::MORSE_BONDS,
      { "MORSE_BONDS", F_MORSE, { "beta", "De", "D0", "bondlength" },
        "(De*((1 - exp(-beta*(r-bondlength)))^2-1)+D0);" } },
    { Potential::HUA_BONDS,
      { "HUA_BONDS", F_HUA, { "bondlength", "De", "b", "c" },
        "(De*(((1-myexp)/(1-c*myexp))^2 -1));myexp=exp(-b*(r-bondlength))" } },
    { Potential::LINEAR_ANGLES, { "LINEAR_ANGLES", F_LINEAR_ANGLES, { "a", "klin" }, "0" } },
    { Potential::UREY_BRADLEY_ANGLES, { "UREY_BRADLEY_ANGLES", F_UREY_BRADLEY, { "kt", "ub_angle", "r13", "kub" }, "0" } },
    { Potential::HARMONIC_DIHEDRALS, { "HARMONIC_DIHEDRALS", F_IDIHS, { "kimp" }, "0" } },
    { Potential::FOURIER_DIHEDRALS, { "FOURIER_DIHEDRALS", F_FOURDIHS, { "c0", "c1", "c2", "c3", "c4", "c5" }, "0" } },
    { Potential::PROPER_DIHEDRALS, { "PROPER_DIHEDRALS", F_PDIHS, { "phi0", "kp", "mult" }, "0" } },
    { Potential::POLARIZATION, { "POLARIZATION", F_POLARIZATION, { "alpha", "rhyper", "fchyper" }, "0" } },
    { Potential::VSITE1, { "VSITE1", F_VSITEN, { "vs1a" }, "0" } },
    { Potential::VSITE2, { "VSITE2", F_VSITE2, { "vs2a" }, "0" } },
    { Potential::VSITE2FD, { "VSITE2FD", F_VSITE2FD, { "vs2fd_a" }, "0" } },
    { Potential::VSITE3, { "VSITE3", F_VSITE3, { "vs3a", "vs3b" }, "0" } },
    { Potential::VSITE3S, { "VSITE3S", F_VSITE3, { "vs3sa" }, "0" } },
    { Potential::VSITE3FD, { "VSITE3FD", F_VSITE3FD, { "vs3fd_a", "vs3fd_b" }, "0" } },
    { Potential::VSITE3FAD, { "VSITE3FAD", F_VSITE3FAD, { "vs3fad_a", "vs3fad_b" }, "0" } },
    { Potential::VSITE3OUT, { "VSITE3OUT", F_VSITE3OUT, { "vs3out_a", "vs3out_b", "vs3out_c" }, "0" } },
    { Potential::VSITE3OUTS, { "VSITE3OUTS", F_VSITE3OUTS, { "vs3outs_a", "vs3outs_c" }, "0" } }
};

const std::string &potentialToString(Potential p)
{
    auto pp = potprops.find(p);
    if (potprops.end() == pp)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No support for potential %d", static_cast<int>(p))));
    }
    return pp->second.name;
}

const std::vector<const char *> potentialToParameterName(Potential p)
{
    auto pp = potprops.find(p);
    if (potprops.end() == pp)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No support for potential %d", static_cast<int>(p))));
    }
    return pp->second.param;
}

const std::string &potentialToEnergy(Potential p)
{
    auto pp = potprops.find(p);
    if (potprops.end() == pp)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No support for potential %d", static_cast<int>(p))));
    }
    return pp->second.energy;
}

bool stringToPotential(const std::string &pname, Potential *p)
{
    // Search for native ACT names first.
    for (const auto &ps : potprops)
    {
        if (pname == ps.second.name)
        {
            *p = ps.first;
            return true;
        }
    }
    // Now search for old-style GROMACS names
    for (const auto &ps : potprops)
    {
        if (ps.second.ftype >= 0 && ps.second.ftype < F_NRE)
        {
            if (pname.compare(interaction_function[ps.second.ftype].name) == 0 ||
                pname.compare(interaction_function[ps.second.ftype].longname) == 0)
            {
                *p = ps.first;
                return true;
            }
        }
    }
    return false;
}

int potentialToGromacsType(Potential p)
{
    auto ppp = potprops.find(p);
    if (potprops.end() != ppp)
    {
        return ppp->second.ftype;
    }
    return -1;
}

const char *potentialToGromacsString(Potential p)
{
    int ftype = potentialToGromacsType(p);
    if (-1 == ftype)
    {
        return nullptr;
    }
    else
    {
        return interaction_function[ftype].name;
    }
}

static const std::map<ChargeType, Potential> cp = {
    { ChargeType::Point, Potential::COULOMB_POINT },
    { ChargeType::Gaussian, Potential::COULOMB_GAUSSIAN },
    { ChargeType::Slater, Potential::COULOMB_SLATER }
};

Potential chargeTypeToPotential(ChargeType c)
{
    return cp.find(c)->second;
}

ChargeType potentialToChargeType(Potential p)
{
    ChargeType ct = ChargeType::Point;
    for (const auto &cccp : cp)
    {
        if (cccp.second == p)
        {
            ct = cccp.first;
        }
    }
    return ct;
}

} // namespace alexandria
