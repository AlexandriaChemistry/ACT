/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

#define sqr(x) (x*x)

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
        { CombRule::Kronecker, "Kronecker" }
    };

const std::string &combinationRuleName(CombRule c)
{
    auto crfind = combRuleName.find(c);
    if (crfind == combRuleName.end())
    {
        GMX_THROW(gmx::InternalError("Unsupported combination rule"));
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

double combineTwo(CombRule comb, double x1, double x2)
{
    switch (comb)
    {
    case CombRule::Geometric:
        return std::sqrt(x1*x2);
    case CombRule::Arithmetic:
        return 0.5 * (x1 + x2);
    case CombRule::HogervorstEpsilon:
        return (2.0 * x1 * x2)/( x1 + x2 );
    case CombRule::Yang:
        if (x1*x2 == 0)
        {
            return 0;
        }
        else
        {
            return ( x1 * x2 ) * ( x1 + x2 ) / (sqr(x1) + sqr(x2) );
        }
    case CombRule::WaldmanSigma:
        return std::pow(0.5*(std::pow(x1, 6) + std::pow(x2, 6)), (1.0/6.0));
    case CombRule::Volumetric:
        // Kriz2024a Eqn. 20
        return std::pow(0.5*(x1*x1*x1 + x2*x2*x2), (1.0/3.0));
    case CombRule::InverseSquare:
        // Kriz2024a Eqn. 21
        if (x1 == 0 || x2 == 0)
        {
            return 0;
        }
        else
        {
            double denom = 1/sqr(x1) + 1/sqr(x2);
            return std::sqrt(2/denom);
        }
    case CombRule::QiSigma:
        if (x1 == 0 && x2 == 0)
        {
            return 0;
        }
        // Qi2016a Eqn. 2 
        return (x1*x1*x1+x2*x2*x2)/(x1*x1+x2*x2);
    case CombRule::HalgrenEpsilon:
        // Halgren1992a, Eqn. 14
        return 4*x1*x2/sqr(std::sqrt(x1) + std::sqrt(x2));
    default: // throws
        GMX_THROW(gmx::InternalError(gmx::formatString("Unknown combination rule %s",
                                                       combRuleName.find(comb)->second.c_str()).c_str()));
    }
    return 0;
}

double combineHogervorstSigma(double e1, double e2, double g1, double g2, double s1, double s2)
{
    if (g1 <= 6 || g2 <= 6)
    {
        GMX_THROW(gmx::InvalidInputError("Combination rule HogervorstSigma not defined if gamma1 or gamma2 <= 6"));
    }
    double tempi = std::abs(e1 * g1 * (std::pow(s1, 6)) /(g1 - 6 ));
    double tempj = std::abs(e2 * g2 * (std::pow(s2, 6)) /(g2 - 6  ));
    double gam12 = combineTwo(CombRule::Arithmetic, g1, g2);
    double eps12 = combineTwo(CombRule::HogervorstEpsilon, e1, e2);
    if (0 == eps12)
    {
        return 0;
    }
    return std::pow((std::sqrt( tempi * tempj ) )* std::abs(gam12 - 6) / (gam12 * eps12), 1.0/6.0);
}

double combineWaldmanEpsilon(double e1, double e2, double s1, double s2)
{
    // Qi2106 Eqn 3.
    double s13 = s1*s1*s1;
    double s23 = s2*s2*s2;
    if (s13 == 0 && s23 == 0)
    {
        return 0;
    }
    return std::sqrt( e1 * e2 ) * ( 2 * s13 * s23 / (sqr(s13) + sqr(s23)));
}

double combineMasonGamma(double g1, double g2, double s1, double s2)
{
    double sigmaIJ = combineTwo(CombRule::Geometric, s1, s2);
    return sigmaIJ * (0.5*((g1/s1)+(g2/s2)));
}

std::map<const std::string, CombRule> getCombinationRule(const ForceFieldParameterList &vdw)
{
    std::map<const std::string, CombRule> myCombRule;
    std::string oldCombRule("combination_rule");
    if (vdw.optionExists(oldCombRule))
    {
        GMX_THROW(gmx::InvalidInputError("Old style combination_rule not supported anymore"));
    }
    else
    {
        for(const auto &opt : vdw.combinationRules())
        {
            CombRule cr;
            if (combinationRuleRule(opt.second, &cr))
            {
                myCombRule.insert({ opt.first, cr });
            }
            else
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Invalid combination rule name %s for parameter %s", opt.second.c_str(), opt.first.c_str()).c_str()));
            }
        }
    }
    return myCombRule;
}

void evalCombinationRule(Potential                                    ftype,
                         const std::map<const std::string, CombRule> &combrule,
                         const ForceFieldParameterMap                &ivdw,
                         const ForceFieldParameterMap                &jvdw,
                         bool                                         includePair,
                         ForceFieldParameterMap                      *pmap)
{
    // Fudge unit
    std::string unit("kJ/mol");

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
            for(auto c: combrule)
            {
                allrules += " " + c.first;
            }
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Parameter %s not found. There are combination rules for: %s.", param.first.c_str(), allrules.c_str()).c_str()));
        }
        auto   crule = combrule.find(param.first)->second;
        double value = 0;
        if (Potential::EXPONENTIAL == ftype ||
            Potential::DOUBLEEXPONENTIAL == ftype ||
            Potential::BUCKINGHAM == ftype ||
            Potential::TANG_TOENNIES == ftype ||
            Potential::TT2b == ftype ||
            Potential::MORSE_BONDS == ftype)
        {
            if (CombRule::Kronecker == crule)
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
            else
            {
                value = combineTwo(crule,
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
            case Potential::LJ8_6:
                cdist    = vdwname[lj8_6SIGMA];
                cepsilon = vdwname[lj8_6EPSILON];
                break;
            default:
                gmx_fatal(FARGS, "Please implement support for potential %s",
                          potentialToString(ftype).c_str());             
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
                                                               combinationRuleName(crule).c_str(),
                                                               potentialToString(ftype).c_str()).c_str()));
            }
            isig = ivdw.find(cdist)->second.value();
            jsig = jvdw.find(cdist)->second.value();
            switch (crule)
            {
            case CombRule::HogervorstSigma:
                value = combineHogervorstSigma(ieps, jeps, igam, jgam, isig, jsig);
                break;
            case CombRule::WaldmanEpsilon:
                value = combineWaldmanEpsilon(ieps, jeps, isig, jsig);
                break;
            case CombRule::MasonGamma:
                value = combineMasonGamma(igam, jgam, isig, jsig);
                break;
            default: // perform correct default action, replace by all values
                value = combineTwo(crule,
                                   ivdw.find(param.first)->second.value(),
                                   jvdw.find(param.first)->second.value());
                break;
            }
        }
        pmap->insert_or_assign(param.first,
                               ForceFieldParameter(unit, value, 0, 1, value, value, mutd, true, true));
    }
}

static void generateParameterPairs(ForceField      *pd,
                                   InteractionType  itype,
                                   bool             force)
{
    // Do not crash if e.g. there is no VDWCORRECTION.
    if (!pd->interactionPresent(itype))
    {
        return;
    }
    auto forcesVdw = pd->findForces(itype);
    auto comb_rule = getCombinationRule(*forcesVdw);

    // We temporarily store the new parameters here
    ForceFieldParameterListMap *parm = forcesVdw->parameters();;
    
    // Now do the double loop
    for (auto &ivdw : *forcesVdw->parameters())
    {
        auto &iid    = ivdw.first;
        // Check whether this is a single atom parameter
        if (iid.atoms().size() > 1)
        {
            continue;
        }
        auto &iparam = ivdw.second;
        bool iupdated = false;
        for(const auto &ip : iparam)
        {
            iupdated = iupdated || ip.second.updated();
        }
        for (auto &jvdw : *forcesVdw->parameters())
        {
            auto &jid    = jvdw.first;
            // Check whether this is a single atom parameter and
            // whether this is is larger or equal to iid.
            if (jid.atoms().size() > 1 || jid.id() < iid.id())
            {
                continue;
            }
            // Test whether or not to include this pair.
            // This will only be used in case the potential is exponential.
            bool includePair = true;
            if (Potential::EXPONENTIAL == forcesVdw->potential())
            {
                auto ai = iid.atoms()[0];
                auto aj = jid.atoms()[0];
                if (pd->hasParticleType(ai) && pd->hasParticleType(aj))
                {
                    auto pti = pd->findParticleType(ai)->apType();
                    auto ptj = pd->findParticleType(aj)->apType();
                    // The pair will be included only if one particle is a vsite
                    // and the other an atom.
                    includePair = ((ActParticle::Vsite == pti && ActParticle::Atom == ptj) ||
                                   (ActParticle::Vsite == ptj && ActParticle::Atom == pti));
                }
            }
            auto &jparam = jvdw.second;
            bool jupdated = false;
            for(const auto &jp : jparam)
            {
                jupdated = jupdated || jp.second.updated();
            }
            if (!(iupdated || jupdated || force))
            {
                continue;
            }
            // Fill the parameters, potential dependent
            ForceFieldParameterMap pmap;
            evalCombinationRule(forcesVdw->potential(),
                                comb_rule, ivdw.second, jvdw.second, includePair, &pmap);

            parm->insert_or_assign(Identifier({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes),
                                   std::move(pmap));
        }
    }
    // Phew, we're done!
}

static void generateCoulombParameterPairs(ForceField *pd, bool force)
{
    auto forcesCoul = pd->findForces(InteractionType::ELECTROSTATICS);
    
    // Fudge unit
    std::string unit("kJ/mol");
    
    // We use dependent mutability to show these are not independent params
    auto mutd = Mutability::Dependent;
    auto cname = potentialToParameterName(Potential::COULOMB_GAUSSIAN);
    auto zeta = cname[coulZETA];
    // Finally add the new parameters to the exisiting list
    auto fold = forcesCoul->parameters();
    // Now do the double loop
    int nid = 0;
    for (auto &icoul : *forcesCoul->parameters())
    {
        auto &iid    = icoul.first;
        // Check whether this is a single atom parameter
        if (iid.atoms().size() > 1)
        {
            continue;
        }
        double izeta = icoul.second[zeta].internalValue();
        for (auto &jcoul : *forcesCoul->parameters())
        {
            auto &jid    = jcoul.first;
            // Check whether this is a single atom parameter and
            // whether this is is larger or equal to iid.
            if (jid.atoms().size() > 1 || jid.id() < iid.id() ||
                (!icoul.second[zeta].updated() && !jcoul.second[zeta].updated() && !force))
            {
                continue;
            }
            double     jzeta  = jcoul.second[zeta].internalValue();
            Identifier pairID({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes);
            nid += 1;
            auto       oldfp  = fold->find(pairID);
            if (oldfp == fold->end())
            {
                ForceFieldParameterMap ffpm = {
                    { cname[coulZETA], 
                      ForceFieldParameter(unit, izeta, 0, 0, izeta, izeta, mutd, false, true) },
                    { cname[coulZETA2],
                      ForceFieldParameter(unit, jzeta, 0, 0, jzeta, jzeta, mutd, false, true) }
                };
                fold->insert({pairID, ffpm});
            }
            else
            {
                auto &pi = oldfp->second.find(cname[coulZETA])->second;
                pi.forceSetValue(izeta);
                auto &pj = oldfp->second.find(cname[coulZETA2])->second;
                pj.forceSetValue(jzeta);
            }
        }
    }
    if (debug)
    {
        int np = forcesCoul->parameters()->size();
        fprintf(debug, "Made %d/%d identifiers in generateCoulombParameterPairs\n", nid, np);
    }
    // Phew, we're done!
}

void generateDependentParameter(ForceField *pd, bool force)
{
    generateParameterPairs(pd, InteractionType::VDW, force);
    generateParameterPairs(pd, InteractionType::VDWCORRECTION, force);
    generateParameterPairs(pd, InteractionType::INDUCTIONCORRECTION, force);
    generateCoulombParameterPairs(pd, force);
}

} // namespace alexandria
