/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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
#include "combinationrules.h"

#include <cmath>
#include <cstring>
#include <string>

#include "act/basics/mutability.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/potential.h"
#include "gromacs/math/functions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

const std::map<CombRule, const std::string> combRuleName = 
    {
        { CombRule::Geometric, "Geometric" },
        { CombRule::Arithmetic, "Arithmetic" },
        { CombRule::HogervorstEpsilon, "HogervorstEpsilon" },
        { CombRule::HogervorstSigma, "HogervorstSigma" },
        { CombRule::Yang, "Yang" },
        { CombRule::WaldmanSigma, "WaldmanSigma" },
        { CombRule::WaldmanEpsilon, "WaldmanEpsilon" },
        { CombRule::QiSigma, "QiSigma" },
        { CombRule::QiEpsilon, "QiEpsilon" },
        { CombRule::MasonGamma, "MasonGamma" },
        { CombRule::Volumetric, "Volumetric" },
        { CombRule::InverseSquare, "InverseSquare" },
        { CombRule::HalgrenEpsilon, "HalgrenEpsilon" },
        { CombRule::Kronecker, "Kronecker" },
        { CombRule::GeneralizedMean, "GeneralizedMean" }
    };

const std::string &combinationRuleName(CombRule c)
{
    const auto crfind = combRuleName.find(c);
    if (crfind == combRuleName.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Unsupported combination rule %d",
                                                       static_cast<int>(c)).c_str()));
    }
    return crfind->second;
}

bool combinationRuleRule(const std::string &name, CombRule *cr)
{
    for(const auto &m : combRuleName)
    {
        if (m.second == name)
        {
            *cr = m.first;
            return true;
        }
    }
    return false;
}

ParamCombRule::ParamCombRule(const std::string &rule)
{
    CombRule c;
    if (!combinationRuleRule(rule, &c))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown combination rule '%s' (note these are case-sensitive)", rule.c_str())));
    }
    rule_ = c;
}

// The formulas are in combinationkernels.h, here the input is checked

double combineTwo(CombRule comb, double x1, double x2)
{
    bool   ok    = false;
    double value = combination::two<combination::StdMath>(comb, x1, x2, &ok);
    if (!ok)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Unknown combination rule %s",
                                                       combRuleName.find(comb)->second.c_str()).c_str()));
    }
    return value;
}

double combineHogervorstSigma(double e1, double e2, double g1, double g2, double s1, double s2)
{
    if (g1 <= 6 || g2 <= 6)
    {
        GMX_THROW(gmx::InvalidInputError("Combination rule HogervorstSigma not defined if gamma1 or gamma2 <= 6"));
    }
    return combination::hogervorstSigma<combination::StdMath>(e1, e2, g1, g2, s1, s2);
}

double combineWaldmanEpsilon(double e1, double e2, double s1, double s2)
{
    return combination::waldmanEpsilon<combination::StdMath>(e1, e2, s1, s2);
}

double combineGeneralizedMean(double x1, double x2, double exponent)
{
    return combination::generalizedMean<combination::StdMath>(x1, x2, exponent);
}

//! \return the gamma of Mason1955a, see combination::masonGamma
static double combineMasonGamma(double g1, double g2, double s1, double s2)
{
    return combination::masonGamma<combination::StdMath>(g1, g2, s1, s2);
}

void evalCombinationRule(Potential                     ftype,
                         const CombRuleSet            &combrule,
                         const ForceFieldParameterMap &ivdw,
                         const ForceFieldParameterMap &jvdw,
                         bool                          includePair,
                         ForceFieldParameterMap       *pmap)
{
    // Fudge unit
    const std::string unit("kJ/mol");

    // We use dependent mutability to show these are not independent params
    auto mutd = Mutability::Dependent;

    for(const auto &param : ivdw)
    {
        if (jvdw.end() == jvdw.find(param.first))
        {
            GMX_THROW(gmx::InternalError("Van der Waals parameters of different types cannot be combined"));
        }
        if (combrule.end() == combrule.find(param.first))
        {
            std::string allrules;
            for(const auto& c: combrule)
            {
                allrules += " " + c.first;
            }
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Parameter %s not found. There are combination rules for: %s.", param.first.c_str(), allrules.c_str()).c_str()));
        }
        const auto   crule = combrule.find(param.first)->second;
        double value = 0;
        if (Potential::BORN_MAYER == ftype ||
            Potential::MACDANIEL_SCHMIDT == ftype ||
            Potential::BUCKINGHAM == ftype ||
            Potential::TANG_TOENNIES == ftype ||
            Potential::TT2b == ftype ||
            Potential::SLATER_ISA_TT == ftype ||
            Potential::MORSE_BONDS == ftype)
        {
            if (CombRule::Kronecker == crule.rule())
            {
                if (includePair)
                {
                    value = combineTwo(CombRule::Arithmetic,
                                       ivdw.find(param.first)->second.value(),
                                       jvdw.find(param.first)->second.value());
                }
                else
                {
                    value = 0;
                }
            }
            else if (CombRule::GeneralizedMean == crule.rule())
            {
                value = combineGeneralizedMean(ivdw.find(param.first)->second.value(),
                                               jvdw.find(param.first)->second.value(),
                                               crule.ffplConst().value());
            }
            else if (CombRule::HogervorstSigma == crule.rule() ||
                     CombRule::WaldmanEpsilon == crule.rule() ||
                     CombRule::QiEpsilon == crule.rule() ||
                     CombRule::MasonGamma == crule.rule())
            {
                // These need the sigma and epsilon of the potential
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Combination rule %s for %s needs sigma and epsilon, which potential %s does not have",
                                                                   combinationRuleName(crule.rule()).c_str(),
                                                                   param.first.c_str(),
                                                                   potentialToString(ftype).c_str()).c_str()));
            }
            else
            {
                value = combineTwo(crule.rule(),
                                   ivdw.find(param.first)->second.value(),
                                   jvdw.find(param.first)->second.value());
            }
        }
        else
        {
            // Defining some strings that we may or may not need
            auto vdwname = potentialToParameterName(ftype);
            std::string cdist, cepsilon, cgamma;
            switch (ftype)
            {
            case Potential::GENERALIZED_BUCKINGHAM:
                cdist    = vdwname[gbhRMIN];
                cepsilon = vdwname[gbhEPSILON];
                cgamma   = vdwname[gbhGAMMA];
                break;
            case Potential::WANG_BUCKINGHAM:
                cdist    = vdwname[wbhSIGMA];
                cepsilon = vdwname[wbhEPSILON];
                cgamma   = vdwname[wbhGAMMA];
                break;
            case Potential::LJ14_7:
                cdist    = vdwname[lj14_7SIGMA];
                cepsilon = vdwname[lj14_7EPSILON];
                cgamma   = vdwname[lj14_7GAMMA];
                break;
            case Potential::LJ12_6:
                cdist    = vdwname[lj12_6SIGMA];
                cepsilon = vdwname[lj12_6EPSILON];
                break;
            case Potential::LJ12_6_4:
                cdist    = vdwname[lj12_6_4SIGMA];
                cepsilon = vdwname[lj12_6_4EPSILON];
                cgamma   = vdwname[lj12_6_4GAMMA];
                break;
            case Potential::LJ8_6:
                cdist    = vdwname[lj8_6SIGMA];
                cepsilon = vdwname[lj8_6EPSILON];
                break;
            default:
                gmx_fatal(FARGS, "Please implement support for potential '%s' ftype %d",
                          potentialToString(ftype).c_str(),
                          static_cast<int>(ftype));   
            }
            auto ieps = ivdw.find(cepsilon)->second.value();
            auto jeps = jvdw.find(cepsilon)->second.value();
            double igam = 0;
            double jgam = 0;
            if (ivdw.end() != ivdw.find(cgamma))
            {
                igam = ivdw.find(cgamma)->second.value();
            }
            if (jvdw.end() != jvdw.find(cgamma))
            {
                jgam = jvdw.find(cgamma)->second.value();
            }
            double isig = 0, jsig = 0;
            if (ivdw.find(cdist) == ivdw.end() || jvdw.find(cdist) == jvdw.end())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Parameter %s missing trying to apply combination rule %s for potential %s",
                                                               cdist.c_str(),
                                                               combinationRuleName(crule.rule()).c_str(),
                                                               potentialToString(ftype).c_str()).c_str()));
            }
            isig = ivdw.find(cdist)->second.value();
            jsig = jvdw.find(cdist)->second.value();
            switch (crule.rule())
            {
            case CombRule::GeneralizedMean:
                value = combineGeneralizedMean(ivdw.find(param.first)->second.value(),
                                               jvdw.find(param.first)->second.value(),
                                               crule.ffplConst().value());
                break;
            case CombRule::HogervorstSigma:
                value = combineHogervorstSigma(ieps, jeps, igam, jgam, isig, jsig);
                break;
            case CombRule::WaldmanEpsilon:
            case CombRule::QiEpsilon:
                // Qi et al. use the epsilon of Waldman and Hagler, Eqn. 3 of Qi2016a
                value = combineWaldmanEpsilon(ieps, jeps, isig, jsig);
                break;
            case CombRule::MasonGamma:
                value = combineMasonGamma(igam, jgam, isig, jsig);
                break;
            default: // perform correct default action, replace by all values
                value = combineTwo(crule.rule(),
                                   ivdw.find(param.first)->second.value(),
                                   jvdw.find(param.first)->second.value());
                break;
            }
        }
        pmap->insert_or_assign(param.first,
                               ForceFieldParameter(unit, value, 0, 1, value, value, mutd, true, true));
    }
}

} // namespace alexandria
