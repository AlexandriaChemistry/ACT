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
    { Potential::LJ12_6,
      { "LJ12_6" , F_LJ, { "sigma", "epsilon" },
        "4*epsilon*((sigma/r)^12 - (sigma/r)^6)", "" } },
    { Potential::LJ12_6_4,
      { "LJ12_6_4" , -1, { "sigma", "epsilon", "gamma" },
        "4*epsilon*((sigma/r)^12 - (sigma/r)^6) - gamma*(1/r)^4", "" } },
    { Potential::LJ8_6,
      {  "LJ8_6", -1, { "sigma", "epsilon" },
         "epsilon*(3*(sigma/r)^8 - 4*(sigma/r)^6)", ""} },
    { Potential::LJ14_7,
      { "LJ14_7", -1, { "sigma", "epsilon", "gamma", "delta" },
        "select(epsilon*sigma,( epsilon*( ( (1+ delta)/((r/sigma)+ delta))^7 ) * ( ( (1+ gamma)/(((r/sigma)^7) +gamma )  ) -2 ) ),0)", "" } },
    { Potential::BORN_MAYER,
      { "BORN_MAYER", -1, { "aexp", "bexp", "vsite" },
        // Implement xor function for Kronecker combination rule
        "-aexp*exp(-bexp*r); bexp=0.5*(bexp1+bexp2); aexp=select(vsite1*vsite2+(1-vsite1)*(1-vsite2),0,sqrt(aexp1*aexp2))",
        "" } },
    { Potential::SLATER_ISA,
      { "SLATER_ISA", -1, { "aexp", "bexp" },
        "aexp*exp(-br)*(br*br/3 + br + 1)", 
        "br=bexp*r; bexp=0.5*(bexp1+bexp2); aexp=sqrt(aexp1*aexp2)" } },
    { Potential::MACDANIEL_SCHMIDT,
      { "MACDANIEL_SCHMIDT", -1, { "a1dexp", "a2dexp", "bdexp" },
        "(Ab-Aa)*exp(-b*r)",
        "b=0.5*(bdexp1+bdexp2); Aa=sqrt(a1dexp1*a1dexp2); Ab=sqrt(a2dexp1*a2dexp2)" } },
    { Potential::WANG_BUCKINGHAM,
      { "WANG_BUCKINGHAM", F_WBHAM, { "sigma", "epsilon", "gamma" }, 
        "select(epsilon*sigma,(((2*epsilon)/(1-(3/(gamma+3)))) * (1.0/(1.0+(r/sigma)^6)) * ((3/(gamma+3))*exp(gamma*(1-(r/sigma)))-1)),0)", "" } },
    { Potential::BUCKINGHAM,
      { "BUCKINGHAM", -1, { "Abh", "bbh", "c6bh" },
        "0", "" } },
    { Potential::TANG_TOENNIES,
      { "TANG_TOENNIES", -1, { "Att", "btt", "c6tt", "c8tt", "c10tt" },
        "0", "" } },
    { Potential::SLATER_ISA_TT,
      { "SLATER_ISA_TT", -1, { "A", "bExch", "bDisp", "c6", "c8", "c10" },
        "A*exp(-bExchR)*((bExchR^2)/3 + bExchR + 1) - ((1-ebDr*sum6)*c6n + (1-ebDr*sum8)*c8n + (1-ebDr*sum10)*c10n)",
        "bExchR=bExch*r; sum10=sum8+br9/362880+br10/3628800; sum8=sum6+br7/5040+br8/40320; sum6=br0+br1+br2/2+br3/6+br4/24+br5/120+br6/720 ; ebDr=exp(-br1); c6n=c6/(r^6); c8n=c8/(r^8); c10n=c10/(r^10); br10=br5*br5; br9=br5*br4; br8=br4*br4; br7=br4*br3; br6=br3*br3; br5=br3*br2; br4=br2*br2; br3=br1*br2; br2=br1*br1; br1=(bDisp*r); br0=1;" } },
    { Potential::TT2b,
      { "TT2b", -1, { "Att2b", "bExchtt2b", "bDisptt2b", "c6tt2b", "c8tt2b", "c10tt2b" },
        "Att2b*exp(-bExchtt2b*r) - ((1-ebDr*sum6)*c6n + (1-ebDr*sum8)*c8n + (1-ebDr*sum10)*c10n)", 
        "sum10=sum8+br9/362880+br10/3628800; sum8=sum6+br7/5040+br8/40320; sum6=br0+br1+br2/2+br3/6+br4/24+br5/120+br6/720 ; ebDr=exp(-br1); c6n=c6tt2b/(r^6); c8n=c8tt2b/(r^8); c10n=c10tt2b/(r^10); br10=br5*br5; br9=br5*br4; br8=br4*br4; br7=br4*br3; br6=br3*br3; br5=br3*br2; br4=br2*br2; br3=br1*br2; br2=br1*br1; br1=(bDisptt2b*r); br0=1;" } },
    { Potential::GENERALIZED_BUCKINGHAM,
      { "GENERALIZED_BUCKINGHAM", F_GBHAM, { "rmin", "epsilon", "gamma", "delta" },
        "select(epsilon*rmin*gamma,( epsilon*((delta + 2*gamma + 6)/(2*gamma)) * (1/(1+((r/rmin)^6))) * (  ((6+delta)/(delta + 2*gamma + 6)) * exp(gamma*(1-(r/rmin))) -1 ) - (epsilon/(1+(r/rmin)^delta)) ),0)", "" } },
    { Potential::COULOMB_GAUSSIAN, { "COULOMB_GAUSSIAN", -1, { "zeta", "zeta2" }, "0", "" } },
    { Potential::COULOMB_POINT, { "COULOMB_POINT", F_COUL_SR, { "zeta", "zeta2" }, "0", "" } },
    { Potential::COULOMB_SLATER, { "COULOMB_SLATER", -1, { "zeta", "zeta2" }, "0", "" } },
    { Potential::HARMONIC_BONDS, { "HARMONIC_BONDS", F_BONDS, { "kb", "bondlength", "bondenergy" }, "0", "" } },
    { Potential::CUBIC_BONDS,
      { "CUBIC_BONDS", F_CUBICBONDS, { "bondlength", "rmax", "kb", "De" },
        "(kb*(r-bondlength)^2 * (rmax-r) - D_e)*step(2*rmax+bondlength-3*r) + step(3*r-2*rmax-bondlength)*(-D_e - (4*kb/27)*(bondlength-rmax)^2)", ""} },
    { Potential::HARMONIC_ANGLES, { "HARMONIC_ANGLES", F_ANGLES, { "kt", "angle" }, "0", "" } },
    { Potential::MORSE_BONDS,
      { "MORSE_BONDS", F_MORSE, { "beta", "De", "D0", "bondlength" },
        "(De*((1 - exp(-beta*(r-bondlength)))^2-1)+D0)", "" } },
    { Potential::HUA_BONDS,
      { "HUA_BONDS", F_HUA, { "bondlength", "De", "b", "c" },
        "(De*(((1-myexp)/(1-c*myexp))^2 -1));myexp=exp(-b*(r-bondlength))", "" } },
    { Potential::LINEAR_ANGLES, { "LINEAR_ANGLES", F_LINEAR_ANGLES, { "a", "klin" }, "0", "" } },
    { Potential::UREY_BRADLEY_ANGLES, { "UREY_BRADLEY_ANGLES", F_UREY_BRADLEY, { "kt", "ub_angle", "r13", "kub" }, "0", "" } },
    { Potential::HARMONIC_DIHEDRALS, { "HARMONIC_DIHEDRALS", F_IDIHS, { "kimp" }, "0", "" } },
    { Potential::FOURIER_DIHEDRALS, { "FOURIER_DIHEDRALS", F_FOURDIHS, { "c0", "c1", "c2", "c3", "c4", "c5" }, "0", "" } },
    { Potential::PROPER_DIHEDRALS, { "PROPER_DIHEDRALS", F_PDIHS, { "phi0", "kp", "mult" }, "0", "" } },
    { Potential::POLARIZATION, { "POLARIZATION", F_POLARIZATION, { "alpha", "rhyper", "fchyper" }, "0", "" } },
    { Potential::VSITE1, { "VSITE1", F_VSITEN, { "vs1a" }, "0", "" } },
    { Potential::VSITE2, { "VSITE2", F_VSITE2, { "vs2a" }, "0", "" } },
    { Potential::VSITE2FD, { "VSITE2FD", F_VSITE2FD, { "vs2fd_a" }, "0", "" } },
    { Potential::VSITE3, { "VSITE3", F_VSITE3, { "vs3a", "vs3b" }, "0", "" } },
    { Potential::VSITE3S, { "VSITE3S", F_VSITE3, { "vs3sa" }, "0", "" } },
    { Potential::VSITE3FD, { "VSITE3FD", F_VSITE3FD, { "vs3fd_a", "vs3fd_b" }, "0", "" } },
    { Potential::VSITE3FAD, { "VSITE3FAD", F_VSITE3FAD, { "vs3fad_a", "vs3fad_b" }, "0", "" } },
    { Potential::VSITE3OUT, { "VSITE3OUT", F_VSITE3OUT, { "vs3out_a", "vs3out_b", "vs3out_c" }, "0", "" } },
    { Potential::VSITE4, { "VSITE4", -1, { "vs4a", "vs4b", "vs4c" }, "0", "" } },
    { Potential::VSITE4S, { "VSITE4S", -1, { "vs4sa", "vs4sb" }, "0", "" } },
    { Potential::VSITE3OUTS, { "VSITE3OUTS", F_VSITE3OUTS, { "vs3outs_a", "vs3outs_c" }, "0", "" } }
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

const std::string &potentialToPreFactor(Potential p)
{
    auto pp = potprops.find(p);
    if (potprops.end() == pp)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No support for potential %d", static_cast<int>(p))));
    }
    return pp->second.prefactor;
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

//! \brief Map from ChargeDistributionType to Potential for Coulomb interactions
static const std::map<ChargeDistributionType, Potential> cp = {
    { ChargeDistributionType::Point, Potential::COULOMB_POINT },
    { ChargeDistributionType::Gaussian, Potential::COULOMB_GAUSSIAN },
    { ChargeDistributionType::Slater, Potential::COULOMB_SLATER }
};

Potential chargeDistributionTypeToPotential(ChargeDistributionType c)
{
    return cp.find(c)->second;
}

ChargeDistributionType potentialToChargeDistributionType(Potential p)
{
    ChargeDistributionType ct = ChargeDistributionType::Point;
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
