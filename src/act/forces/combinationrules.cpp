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
#include "combinationrules.h"

#include <cmath>
#include <cstring>
#include <string>

#include "act/basics/mutability.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/potential.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/md_enums.h"
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
    default:
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

std::map<const std::string, CombRule> oldCombinationRule(const std::string &vdw_comb,
                                                         Potential          ftype)
{
    std::map<const std::string, CombRule> myCombRule;
    int  i;
    for(i = 0; i < eCOMB_NR; i++)
    {
        if (vdw_comb.compare(ecomb_names[i]) == 0)
        {
            break;
        }
    }
    GMX_RELEASE_ASSERT(i < eCOMB_NR, gmx::formatString("Cannot find combination rule %s in GROMACS",
                                                       vdw_comb.c_str()).c_str());
    std::string cdist;
    if (Potential::GENERALIZED_BUCKINGHAM == ftype)
    {
        cdist = gbh_name[gbhRMIN];
    }
    else
    {
        cdist = lj14_7_name[lj14_7SIGMA];
    }
    const std::string cepsilon(lj14_7_name[lj14_7EPSILON]);
    const std::string cgamma(lj14_7_name[lj14_7GAMMA]);
    const std::string cdelta(lj14_7_name[lj14_7DELTA]);
    bool haveDelta = Potential::GENERALIZED_BUCKINGHAM == ftype || Potential::LJ14_7 ==ftype;
    bool haveGamma = haveDelta || Potential::WANG_BUCKINGHAM == ftype;
    switch (i)
    {
    case eCOMB_GEOMETRIC:
        myCombRule.insert({ cepsilon, CombRule::Geometric });
        myCombRule.insert({ cdist,    CombRule::Geometric });
        if (haveGamma)
        {
            myCombRule.insert({ cgamma,   CombRule::Geometric });
        }
        if (haveDelta)
        {
            myCombRule.insert({ cdelta,   CombRule::Geometric });
        }
        break;
    case eCOMB_ARITHMETIC:
        myCombRule.insert({ cepsilon, CombRule::Arithmetic });
        myCombRule.insert({ cdist,    CombRule::Arithmetic });
        if (haveGamma)
        {
            myCombRule.insert({ cgamma,   CombRule::Arithmetic });
        }
        if (haveDelta)
        {
            myCombRule.insert({ cdelta,   CombRule::Arithmetic });
        }
        break;
    case eCOMB_LORENTZ_BERTHELOT:
        myCombRule.insert({ cdist,    CombRule::Arithmetic });
        myCombRule.insert({ cepsilon, CombRule::Geometric });
        break;
    case eCOMB_KONG_MASON:
        // Kong, C. L. Combining Rules for Intermolecular Potential Parameters. II. Rules for the Lennard-Jones (12−6) Potential
        // and the Morse Potential. J. Chem. Phys. 1973, 59.
        myCombRule.insert({ cdist,    CombRule::Geometric });
        myCombRule.insert({ cepsilon, CombRule::HogervorstEpsilon });
        myCombRule.insert({ cgamma,   CombRule::MasonGamma });
        if (haveDelta)
        {
            myCombRule.insert({ cdelta, CombRule::Yang });
        }
        break;
    case eCOMB_HOGERVORST:
        // Hogervorst, Physica, Volume: 51, Page: 77, Year: 1971. Combination rules for Buckingham.
        myCombRule.insert({ cdist,    CombRule::HogervorstSigma });
        myCombRule.insert({ cepsilon, CombRule::HogervorstEpsilon });
        myCombRule.insert({ cgamma,   CombRule::Arithmetic });
        if (haveDelta)
        {
            myCombRule.insert({ cdelta, CombRule::Yang });
        }
        break;   
    case eCOMB_YANG:
        // Yang, JPhysChemA, Volume: 122, Page: 1672, Year: 2018. Combination rules for Morse.
        myCombRule.insert({ cdist,    CombRule::Yang });
        myCombRule.insert({ cepsilon, CombRule::HogervorstEpsilon });
        myCombRule.insert({ cgamma,   CombRule::Yang });
        if (haveDelta)
        {
            myCombRule.insert({ cdelta, CombRule::Yang });
        }
        break;
    case eCOMB_QI:
        // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. Combination rules for Buf-14-7.
        // Cubic-mean for sigma, and Waldman-Hagler for epsilon.
        myCombRule.insert({ cdist,    CombRule::QiSigma });
        myCombRule.insert({ cepsilon, CombRule::WaldmanEpsilon });
        myCombRule.insert({ cgamma,   CombRule::Arithmetic });
        if (haveDelta)
        {
            if (Potential::GENERALIZED_BUCKINGHAM == ftype)
            {
               myCombRule.insert({ cdelta, CombRule::Yang });
            }
            else if (Potential::LJ14_7 == ftype)
            {
                myCombRule.insert({ cdelta, CombRule::Geometric });
            }
        }
        break;
    case eCOMB_QI_2:
        // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. Combination rules for Buf-14-7.
        // Cubic-mean for sigma, and Waldman-Hagler for epsilon.
        // 2023 testing, Kriz. is almost the same asi Qi.
        myCombRule.insert({ cdist,    CombRule::QiSigma });
        myCombRule.insert({ cepsilon, CombRule::WaldmanEpsilon });
        myCombRule.insert({ cgamma,   CombRule::WaldmanSigma });
        if (haveDelta)
        {
            if (Potential::GENERALIZED_BUCKINGHAM == ftype)
            {
                myCombRule.insert({ cdelta, CombRule::Yang });
            }
            else if (Potential::LJ14_7 == ftype)
            {
                myCombRule.insert({ cdelta, CombRule::Geometric });
            }
            else
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Don't know how to handle delta for function type %s",
                                                               potentialToString(ftype).c_str()).c_str()));
            }
        }
        break;    
    case eCOMB_WALDMAN_HAGLER:
        // Waldman & Hagler, J. Comp. Chem., Year: 1993.
        myCombRule.insert({ cdist,    CombRule::WaldmanSigma });
        myCombRule.insert({ cepsilon, CombRule::WaldmanEpsilon });
        if (haveGamma)
        {
            myCombRule.insert({ cgamma, CombRule::Yang });
        }
        if (haveDelta)
        {
            myCombRule.insert({ cdelta, CombRule::Yang });
        }
        break;
    case eCOMB_QYQY:
        // Waldman & Hagler, J. Comp. Chem., Year: 1993.
        // Kriz: changing to the best rule for GBHAM Qi, Yang, Qi (with the "yang" for delta)
        myCombRule.insert({ cdist,    CombRule::QiSigma });
        myCombRule.insert({ cepsilon, CombRule::WaldmanEpsilon });
        myCombRule.insert({ cgamma,   CombRule::Yang });
        if (haveDelta)
        {
            myCombRule.insert({ cdelta, CombRule::Yang });
        }
        break;
    case eCOMB_QKmQG:
        // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. The best combination rules for Buf-14-7. Cubic-mean for sigma, and Waldman-Hagler for epsilon. Qi /WH for epsilon, KM for gamma (but with geometric sigmaIJ), qi for sigma and geometric for delta
        myCombRule.insert({ cdist,    CombRule::QiSigma });
        myCombRule.insert({ cepsilon, CombRule::WaldmanEpsilon });
        myCombRule.insert({ cgamma,   CombRule::MasonGamma });
        myCombRule.insert({ cdelta,   CombRule::Geometric });
        break;	    
    default:
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Combination rule %s not supported anymore. Sorry.", ECOMBNAME(i)).c_str()));
    }
    return myCombRule;
}

std::map<const std::string, CombRule> getCombinationRule(const ForceFieldParameterList &vdw)
{
    std::map<const std::string, CombRule> myCombRule;
    std::string oldCombRule("combination_rule");
    if (vdw.optionExists(oldCombRule))
    {
        return oldCombinationRule(vdw.optionValue(oldCombRule), vdw.potential());
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

ForceFieldParameterMap evalCombinationRule(Potential                                    ftype,
                                           const std::map<const std::string, CombRule> &combrule,
                                           const ForceFieldParameterMap                &ivdw,
                                           const ForceFieldParameterMap                &jvdw,
                                           bool                                         same)
{
    // Fudge unit
    std::string unit("kJ/mol");
    
    // We use dependent mutability to show these are not independent params
    auto mutd = Mutability::Dependent;

    ForceFieldParameterMap pmap;

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
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Parameter %s not found. There are combination rules for:%s.", param.first.c_str(), allrules.c_str()).c_str()));
        }
        double value = 0;
        if (Potential::EXPONENTIAL == ftype)
        {
            auto crule = combrule.find(param.first)->second;
            if (CombRule::Kronecker == crule)
            {
                if (same)
                {
                    value = 0;
                }
                else
                {
                    value = combineTwo(CombRule::Geometric,
                                       ivdw.find(param.first)->second.value(),
                                       jvdw.find(param.first)->second.value());
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
            const std::string cepsilon(lj14_7_name[lj14_7EPSILON]);
            const std::string cgamma(lj14_7_name[lj14_7GAMMA]);
            std::string       cdist;
            if (ftype == Potential::GENERALIZED_BUCKINGHAM)
            {
                cdist = gbh_name[gbhRMIN];
            }
            else
            {
                cdist = lj14_7_name[lj14_7SIGMA];
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
                GMX_THROW(gmx::InternalError(gmx::formatString("Parameter %s missing", cdist.c_str()).c_str()));
            }
            isig = ivdw.find(cdist)->second.value();
            jsig = jvdw.find(cdist)->second.value();
            auto   crule    = combrule.find(param.first)->second;
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
            default:
                value = combineTwo(crule,
                                   ivdw.find(param.first)->second.value(),
                                   jvdw.find(param.first)->second.value());
                break;
            }
        }
        std::string pij = gmx::formatString("%s_ij", param.first.c_str());
        pmap.insert({ pij, ForceFieldParameter(unit, value, 0, 1, value, value, mutd, true, true)});
    }
    return pmap;
}

static void generateParameterPairs(ForceField      *pd,
                                   InteractionType  itype,
                                   bool             force)
{
    // Do not crash if e.g. there is no CHARGETRANSFER.
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
        auto iparam = ivdw.second;
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
            bool same = jid.id() == iid.id();
            auto jparam = jvdw.second;
            bool jupdated = false;
            for(const auto &jp : jparam)
            {
                jupdated = jupdated || jp.second.updated();
            }
            if ((!iupdated || !iupdated) && !force)
            {
                continue;
            }
            // Fill the parameters, potential dependent
            auto pmap = evalCombinationRule(forcesVdw->potential(),
                                            comb_rule, ivdw.second, jvdw.second, same);

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
    auto zeta = coul_name[coulZETA];
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
                    { coul_name[coulZETAI], 
                      ForceFieldParameter(unit, izeta, 0, 0, izeta, izeta, mutd, false, true) },
                    { coul_name[coulZETAJ],
                      ForceFieldParameter(unit, jzeta, 0, 0, jzeta, jzeta, mutd, false, true) }
                };
                fold->insert({pairID, ffpm});
            }
            else
            {
                auto &pi = oldfp->second.find(coul_name[coulZETAI])->second;
                pi.forceSetValue(izeta);
                auto &pj = oldfp->second.find(coul_name[coulZETAJ])->second;
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
    generateParameterPairs(pd, InteractionType::CHARGETRANSFER, force);
    generateCoulombParameterPairs(pd, force);
}

} // namespace alexandria
