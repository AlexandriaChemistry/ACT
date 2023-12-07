/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include "act/basics/mutability.h"
#include "act/forcefield/forcefield_parametername.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/fatalerror.h"

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
        { CombRule::HalgrenEpsilon, "HalgrenEpsilon" }
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
    return std::pow((std::sqrt( tempi * tempj ) )* abs(gam12 - 6) / (gam12 * eps12), 1.0/6.0);
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
                                                         int                ftype)
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
    const std::string csigma(lj14_7_name[lj14_7SIGMA]);
    const std::string crmin(gbh_name[gbhRMIN]);
    const std::string cepsilon(lj14_7_name[lj14_7EPSILON]);
    const std::string cgamma(lj14_7_name[lj14_7GAMMA]);
    const std::string cdelta(lj14_7_name[lj14_7DELTA]);
    bool haveDelta = F_GBHAM == ftype || F_LJ14_7 ==ftype;
    bool haveGamma = haveDelta || F_WBHAM == ftype;
    switch (i)
    {
    case eCOMB_GEOMETRIC:
        myCombRule.insert({ cepsilon, CombRule::Geometric });
        if (F_GBHAM == ftype)
        {
            myCombRule.insert({ crmin,    CombRule::Geometric });
        }
        else
        {
            myCombRule.insert({ csigma,   CombRule::Geometric });
        }
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
        if (F_GBHAM == ftype)
        {
            myCombRule.insert({ crmin,    CombRule::Arithmetic });
        }
        else
        {
            myCombRule.insert({ csigma,   CombRule::Arithmetic });
        }
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
        myCombRule.insert({ csigma,   CombRule::Arithmetic });
        myCombRule.insert({ cepsilon, CombRule::Geometric });
        break;
    case eCOMB_KONG_MASON:
        // Kong, C. L. Combining Rules for Intermolecular Potential Parameters. II. Rules for the Lennard-Jones (12−6) Potential
        // and the Morse Potential. J. Chem. Phys. 1973, 59.
        myCombRule.insert({ csigma,   CombRule::Geometric });
        myCombRule.insert({ cepsilon, CombRule::HogervorstEpsilon });
        myCombRule.insert({ cgamma,   CombRule::MasonGamma });
        if (haveDelta)
        {
            myCombRule.insert({ cdelta, CombRule::Yang });
        }
        break;
    case eCOMB_HOGERVORST:
        // Hogervorst, Physica, Volume: 51, Page: 77, Year: 1971. Combination rules for Buckingham.
        myCombRule.insert({ csigma,   CombRule::HogervorstSigma });
        myCombRule.insert({ cepsilon, CombRule::HogervorstEpsilon });
        myCombRule.insert({ cgamma,   CombRule::Arithmetic });
        if (haveDelta)
        {
            myCombRule.insert({ cdelta, CombRule::Yang });
        }
        break;   
    case eCOMB_YANG:
        // Yang, JPhysChemA, Volume: 122, Page: 1672, Year: 2018. Combination rules for Morse.
        myCombRule.insert({ csigma,   CombRule::Yang });
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
        myCombRule.insert({ csigma,   CombRule::QiSigma });
        myCombRule.insert({ cepsilon, CombRule::WaldmanEpsilon });
        myCombRule.insert({ cgamma,   CombRule::Arithmetic });
        if (haveDelta)
        {
            if (F_GBHAM == ftype)
            {
               myCombRule.insert({ cdelta, CombRule::Yang });
            }
            else if (F_LJ14_7 == ftype)
            {
                myCombRule.insert({ cdelta, CombRule::Geometric });
            }
        }
        break;
    case eCOMB_QI_2:
        // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. Combination rules for Buf-14-7.
        // Cubic-mean for sigma, and Waldman-Hagler for epsilon.
        // 2023 testing, Kriz. is almost the same asi Qi.
        if (F_GBHAM == ftype)
        {
            myCombRule.insert({ crmin,   CombRule::QiSigma });
        }
        else
        {
            myCombRule.insert({ csigma,   CombRule::QiSigma });
        }
        myCombRule.insert({ cepsilon, CombRule::WaldmanEpsilon });
        myCombRule.insert({ cgamma,   CombRule::WaldmanSigma });
        if (haveDelta)
        {
            if (F_GBHAM == ftype)
            {
                myCombRule.insert({ cdelta, CombRule::Yang });
            }
            else if (F_LJ14_7 == ftype)
            {
                myCombRule.insert({ cdelta, CombRule::Geometric });
            }
        }
        break;    
    case eCOMB_WALDMAN_HAGLER:
        // Waldman & Hagler, J. Comp. Chem., Year: 1993.
        myCombRule.insert({ csigma,   CombRule::WaldmanSigma });
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
        myCombRule.insert({ csigma,   CombRule::QiSigma });
        myCombRule.insert({ cepsilon, CombRule::WaldmanEpsilon });
        myCombRule.insert({ cgamma,   CombRule::Yang });
        if (haveDelta)
        {
            myCombRule.insert({ cdelta, CombRule::Yang });
        }
        break;
    case eCOMB_QKmQG:
        // Qi, Bioorg. & Med. Chem., Volume: 24, Page: 4911, Year: 2016. The best combination rules for Buf-14-7. Cubic-mean for sigma, and Waldman-Hagler for epsilon. Qi /WH for epsilon, KM for gamma (but with geometric sigmaIJ), qi for sigma and geometric for delta
        myCombRule.insert({ csigma,   CombRule::QiSigma });
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
    const std::string oldCombRule = "combination_rule";
    if (vdw.optionExists(oldCombRule))
    {
        return oldCombinationRule(vdw.optionValue(oldCombRule), vdw.gromacsType());
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

ForceFieldParameterMap evalCombinationRule(int                                          ftype,
                                           const std::map<const std::string, CombRule> &combrule,
                                           const ForceFieldParameterMap                &ivdw,
                                           const ForceFieldParameterMap                &jvdw)
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
            GMX_THROW(gmx::InternalError(gmx::formatString("Parameter %s not found in combination rule", param.first.c_str()).c_str()));
        }
        const std::string csigma(lj14_7_name[lj14_7SIGMA]);
        const std::string crmin(gbh_name[gbhRMIN]);
        const std::string cepsilon(lj14_7_name[lj14_7EPSILON]);
        const std::string cgamma(lj14_7_name[lj14_7GAMMA]);
        auto ieps = ivdw.find(cepsilon)->second.value();
        auto jeps = jvdw.find(cepsilon)->second.value();
        auto igam = ivdw.find(cgamma)->second.value();
        auto jgam = jvdw.find(cgamma)->second.value();
        auto isig = ivdw.find(csigma)->second.value();
        auto jsig = jvdw.find(csigma)->second.value();
        if (F_GBHAM == ftype)
        {
            isig = ivdw.find(crmin)->second.value();
            jsig = jvdw.find(crmin)->second.value();
        }
        double value    = 0;
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
        std::string pij = gmx::formatString("%s_ij", param.first.c_str());
        pmap.insert({ pij, ForceFieldParameter(unit, value, 0, 1, value, value, mutd, true, true)});
    }
    return pmap;
}

static void generateVdwParameterPairs(ForceField *pd)
{
    auto forcesVdw = pd->findForces(InteractionType::VDW);
    auto comb_rule = getCombinationRule(*forcesVdw);

    // We temporarily store the new parameters here
    ForceFieldParameterListMap *parm = forcesVdw->parameters();;
    
    // Now do the double loop
    for (auto &ivdw : *forcesVdw->parameters())
    {
        auto iid    = ivdw.first;
        // Check whether this is a single atom parameter
        if (iid.atoms().size() > 1)
        {
            continue;
        }
        auto iparam = ivdw.second;
        for (auto &jvdw : *forcesVdw->parameters())
        {
            auto jid    = jvdw.first;
            // Check whether this is a single atom parameter and
            // whether this is is larger or equal to iid.
            if (jid.atoms().size() > 1 || jid.id() < iid.id())
            {
                continue;
            }
            auto jparam = jvdw.second;
            // Fill the parameters, potential dependent
            auto pmap = evalCombinationRule(forcesVdw->gromacsType(),
                                            comb_rule, ivdw.second, jvdw.second);

            parm->insert_or_assign(Identifier({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes),
                                   std::move(pmap));
            
        }
    }
    // Phew, we're done!
}

static void generateCoulombParameterPairs(ForceField *pd)
{
    auto forcesCoul = pd->findForces(InteractionType::COULOMB);
    
    // We temporarily store the new parameters here
    ForceFieldParameterList newParams;
    
    // Fudge unit
    std::string unit("kJ/mol");
    
    // We use dependent mutability to show these are not independent params
    auto mutd = Mutability::Dependent;

    // Now do the double loop
    for (auto &icoul : *forcesCoul->parameters())
    {
        auto iid    = icoul.first;
        // Check whether this is a single atom parameter
        if (iid.atoms().size() > 1)
        {
            continue;
        }
        auto iparam = icoul.second;
        double izeta = icoul.second["zeta"].internalValue();
        for (auto &jcoul : *forcesCoul->parameters())
        {
            auto jid    = jcoul.first;
            // Check whether this is a single atom parameter and
            // whether this is is larger or equal to iid.
            if (jid.atoms().size() > 1 || jid.id() < iid.id())
            {
                continue;
            }
            auto jparam = jcoul.second;
            double jzeta  = jcoul.second["zeta"].internalValue();
            Identifier pairID({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes);
            ForceFieldParameter pi(unit, izeta, 0, 0, izeta, izeta, mutd, true, true);
            ForceFieldParameter pj(unit, jzeta, 0, 0, jzeta, jzeta, mutd, true, true);
            newParams.addParameter(pairID, coul_name[coulZETAI], pi);
            newParams.addParameter(pairID, coul_name[coulZETAJ], pj);
        }
    }
    // Finally add the new parameters to the exisiting list
    auto fold = forcesCoul->parameters();
    for(const auto &np : newParams.parametersConst())
    {
        // Remove old copy if it exists
        auto oldfp = fold->find(np.first);
        if (oldfp != fold->end())
        {
            fold->erase(oldfp);
        }
        // Now add the new one
        fold->insert({ np.first, np.second });
    }
    // Phew, we're done!
}

static void generateShellForceConstants(ForceField *pd)
{
    if (!pd->polarizable())
    {
        return;
    }
    auto itype = InteractionType::POLARIZATION;
    auto ffpl  = pd->findForces(itype)->parameters();
    // Loop over particles
    for(const auto &part : pd->particleTypesConst())
    {
        if (part.second.hasOption("poltype"))
        {
            auto shellType = part.second.optionValue("poltype");
            if (ffpl->find(shellType) == ffpl->end())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Missing polarization term for %s", shellType.c_str()).c_str()));
            }
            auto  &parms   = ffpl->find(shellType)->second;
            auto   alpha   = parms.find("alpha")->second.internalValue();
            auto   qshell  = pd->findParticleType(shellType)->charge();
            double kshell  = 0;
            if (alpha > 0 && qshell != 0)
            {
                kshell = gmx::square(qshell)*ONE_4PI_EPS0/alpha;
            }
            std::string fc_name("kshell");
            auto fc_parm = parms.find(fc_name);
            if (parms.end() == fc_parm)
            {
                ForceFieldParameter fc_new("kJ/mol nm2", kshell, 0, 1,
                                           kshell, kshell,
                                           Mutability::Dependent, true, true);
                parms.insert({ fc_name, fc_new });
            }
            else
            {
                fc_parm->second.setMutability(Mutability::Free);
                fc_parm->second.setMinimum(kshell);
                fc_parm->second.setMaximum(kshell);
                fc_parm->second.setValue(kshell);
                fc_parm->second.setMutability(Mutability::Dependent);
            }
        }
    }
}

void generateDependentParameter(ForceField *pd)
{
    generateVdwParameterPairs(pd);
    generateCoulombParameterPairs(pd);
    generateShellForceConstants(pd);
}

} // namespace alexandria
