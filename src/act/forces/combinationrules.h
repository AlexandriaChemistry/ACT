/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2025
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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

#ifndef ACT_COMBINATIONRULES_H
#define ACT_COMBINATIONRULES_H

#include <map>

#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/potential.h"

namespace alexandria
{

//! \brief Class that defines the combination rules supported by ACT
enum class CombRule {
    Geometric, Arithmetic, Volumetric, InverseSquare, 
    HogervorstEpsilon, HogervorstSigma, Yang,
    WaldmanSigma, WaldmanEpsilon, HalgrenEpsilon,
    QiSigma, QiEpsilon, MasonGamma, Kronecker,
    GeneralizedMean
};

//! Map combination rules to strings
extern const std::map<CombRule, const std::string> combRuleName;

/*! Class to combine all combination rule information for one parameter
 * The name of the parameter is stored externally in a std::map
 */
class ParamCombRule
{
private:
    //! Combination rule type
    CombRule            rule_;
    //! Optional float for GeneralizedMean, since this can be trained it uses the ffp.
    ForceFieldParameter ffpl_;
public:
    /*! Constructor
     * \param[in] rule The new combination rule type
     */
    ParamCombRule(CombRule rule) : rule_(rule) {}
    /*! Constructor
     * \param[in] rule String corresponding to combination rule
     * \throws if string is incorrect (does not match a known combination rule)
     */
    ParamCombRule(const std::string &rule);
    /*! \brief Add a force field parameter
     * \param[in] ffpl The new parameter
     */
    void addForceFieldParameter(ForceFieldParameter ffpl)
    {
        ffpl_ = ffpl;
    }
    //! \return the rule
    CombRule rule() const { return rule_; }
    //! \return the ForceFieldParameter
    const ForceFieldParameter &ffplConst() const { return ffpl_; }
    //! \return the ForceFieldParameter for editing
    ForceFieldParameter *ffpl() { return &ffpl_; }
};

//! Packaging all the combination rules
typedef std::map<std::string, ParamCombRule> CombRuleSet;
    
//! \return string corresponding to CombRule c
const std::string &combinationRuleName(CombRule c);

/*! \brief Determine CombRule from name
 * \param[in]  name The combination rule name
 * \param[out] cr   The combination rule
 * \return true if successful, false otherwise
 */
bool combinationRuleRule(const std::string &name, CombRule *cr);

/*! \brief Combination rule using GeneralizedMean
 * Mathematical background at https://en.wikipedia.org/wiki/Generalized_mean
 * \param[in] x1       First value
 * \param[in] x2       Second value
 * \param[in] exponent Power to use
 * \return the combined value
 */
double combineGeneralizedMean(double x1, double x2, double exponent);

/*! \brief Execute a simple combination rule
 * \param[in] comb The combination rule to apply
 * \param[in] x1   First value
 * \param[in] x2   Second value
 * \return The combined value
 */
double combineTwo(CombRule comb, double x1, double x2);

/*! \brief Execute a combination rule according to Waldman & Hagler https://doi.org/10.1002/jcc.540140909
 * \param[in] e1 First epsilon
 * \param[in] e2 Second epsilon
 * \param[in] s1 First sigma
 * \param[in] s2 Second sigma
 * \return The combined value
 */
double combineWaldmanEpsilon(double e1, double e2, double s1, double s2);

/*! \brief Execute a combination rule according to Hogervorst1971a https://doi.org/10.1016/0031-8914(71)90138-8
 * \param[in] e1 First epsilon
 * \param[in] e2 Second epsilon
 * \param[in] g1 First gamma
 * \param[in] g2 Second gamma
 * \param[in] s1 First sigma
 * \param[in] s2 Second sigma
 * \return The combined value
 */
double combineHogervorstSigma(double e1, double e2, double g1, double g2, double s1, double s2);

/*! \brief Generate combined force field parameter map
 * \param[in] ftype    The force function ACT style
 * \param[in] combrule Map of combination rules per parameter
 * \param[in] ivdw     Parameters for particle i
 * \param[in] jvdw     Parameters for particle j
 * \param[in] same     Should be true if i and j are the same particle
 * \param[out] pmap    Force Field Parameter Map with pair entries
 */
void evalCombinationRule(Potential                     ftype,
                         const CombRuleSet            &combrules,
                         const ForceFieldParameterMap &ivdw,
                         const ForceFieldParameterMap &jvdw,
                         bool                          same,
                         ForceFieldParameterMap       *pmap);

} // namespace alexandria

#endif
