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

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "act/basics/mutability.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/potential.h"
#include "act/utility/units.h"
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

combination::PairCombination pairCombination(Potential ftype)
{
    combination::PairCombination pc;
    switch (ftype)
    {
    case Potential::BORN_MAYER:
    case Potential::MACDANIEL_SCHMIDT:
    case Potential::BUCKINGHAM:
    case Potential::TANG_TOENNIES:
    case Potential::TT2b:
    case Potential::SLATER_ISA_TT:
    case Potential::MORSE_BONDS:
        pc.simple = true;
        return pc;
    case Potential::GENERALIZED_BUCKINGHAM:
        pc = { false, gbhRMIN, gbhEPSILON, gbhGAMMA };
        return pc;
    case Potential::WANG_BUCKINGHAM:
        pc = { false, wbhSIGMA, wbhEPSILON, wbhGAMMA };
        return pc;
    case Potential::LJ14_7:
        pc = { false, lj14_7SIGMA, lj14_7EPSILON, lj14_7GAMMA };
        return pc;
    case Potential::LJ12_6:
        pc = { false, lj12_6SIGMA, lj12_6EPSILON, -1 };
        return pc;
    case Potential::LJ12_6_4:
        pc = { false, lj12_6_4SIGMA, lj12_6_4EPSILON, lj12_6_4GAMMA };
        return pc;
    case Potential::LJ8_6:
        pc = { false, lj8_6SIGMA, lj8_6EPSILON, -1 };
        return pc;
    default:
        gmx_fatal(FARGS, "Please implement support for potential '%s' ftype %d",
                  potentialToString(ftype).c_str(),
                  static_cast<int>(ftype));
    }
    return pc;
}

void evalCombinationRule(Potential                     ftype,
                         const CombRuleSet            &combrule,
                         const ForceFieldParameterMap &ivdw,
                         const ForceFieldParameterMap &jvdw,
                         bool                          includePair,
                         ForceFieldParameterMap       *pmap)
{
    // We use dependent mutability to show these are not independent params
    auto mutd = Mutability::Dependent;

    // The parameters of both particle types in the order of ivdw, with
    // their rules, see combination::pairParameters. They are combined in
    // the units of the program, such that particle types may have their
    // parameters in other units, and so are the parameters of the pair.
    size_t                   n = ivdw.size();
    std::vector<std::string> names, units;
    std::vector<double>      vi, vj, exponents, out(n, 0.0);
    std::vector<CombRule>    rules;
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
        const auto &crule = combrule.find(param.first)->second;
        names.push_back(param.first);
        units.push_back(gromacsUnit(param.second.unit()));
        vi.push_back(param.second.internalValue());
        vj.push_back(jvdw.find(param.first)->second.internalValue());
        rules.push_back(crule.rule());
        exponents.push_back(crule.ffplConst().value());
    }
    // The positions of sigma, epsilon and gamma among the parameters
    auto pc = pairCombination(ftype);
    if (!pc.simple)
    {
        auto vdwname = potentialToParameterName(ftype);
        auto where   = [&names, &vdwname](int index)
        {
            if (index < 0)
            {
                return -1;
            }
            auto found = std::find(names.begin(), names.end(), vdwname[index]);
            return found == names.end() ? -1 : static_cast<int>(found - names.begin());
        };
        int dist = where(pc.dist);
        if (dist < 0 && n > 0)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Parameter %s missing trying to apply combination rule %s for potential %s",
                                                           vdwname[pc.dist],
                                                           combinationRuleName(rules[0]).c_str(),
                                                           potentialToString(ftype).c_str()).c_str()));
        }
        int epsilon = where(pc.epsilon);
        if (epsilon < 0 && n > 0)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Parameter %s missing trying to apply combination rule %s for potential %s",
                                                           vdwname[pc.epsilon],
                                                           combinationRuleName(rules[0]).c_str(),
                                                           potentialToString(ftype).c_str()).c_str()));
        }
        pc.dist    = dist;
        pc.epsilon = epsilon;
        pc.gamma   = where(pc.gamma);
    }
    int  bad    = 0;
    auto status = combination::pairParameters<combination::StdMath>(pc, static_cast<int>(n), rules.data(),
                                                                      exponents.data(), vi.data(), vj.data(),
                                                                      includePair, out.data(), &bad);
    switch (status)
    {
    case combination::PairStatus::Ok:
        break;
    case combination::PairStatus::UnknownRule:
        GMX_THROW(gmx::InternalError(gmx::formatString("Unknown combination rule %s",
                                                       combinationRuleName(rules[bad]).c_str()).c_str()));
    case combination::PairStatus::NeedsSigma:
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Combination rule %s for %s needs sigma and epsilon, which potential %s does not have",
                                                           combinationRuleName(rules[bad]).c_str(),
                                                           names[bad].c_str(),
                                                           potentialToString(ftype).c_str()).c_str()));
    case combination::PairStatus::GammaTooSmall:
        GMX_THROW(gmx::InvalidInputError("Combination rule HogervorstSigma not defined if gamma1 or gamma2 <= 6"));
    }
    for (size_t k = 0; k < n; k++)
    {
        pmap->insert_or_assign(names[k],
                               ForceFieldParameter(units[k], out[k], 0, 1, out[k], out[k], mutd, true, true));
    }
}

} // namespace alexandria
