
    /*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022,2023
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

#include "act/forcefield/forcefield.h"

namespace alexandria
{

    //! \brief Class that defines the combination rules supported by ACT
    enum class CombRule 
    {
        Geometric, Arithmetic, Volumetric, InverseSquare, 
            HogervorstEpsilon, HogervorstSigma, Yang,
            WaldmanSigma, WaldmanEpsilon, HalgrenEpsilon,
            QiSigma, QiEpsilon, MasonGamma
            };

    //! Map combination rules to strings            
    extern const std::map<CombRule, const std::string> combRuleName;

    //! \return string corresponding to CombRule c
    const std::string &combinationRuleName(CombRule c);

    /*! \brief Determine CombRule from name
     * \param[in]  name The combination rule name
     * \param[out] cr   The combination rule
     * \return true if successful, false otherwise
     */
    bool combinationRuleRule(const std::string &name, CombRule *cr);

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

    /*! \brief Execute a combination rule according to Mason1955a https://doi.org/10.1063/1.1740561
     * \param[in] g1 First gamma
     * \param[in] g2 Second gamma
     * \param[in] s1 First sigma
     * \param[in] s2 Second sigma
     * \return The combined value
     */
    double combineMasonGamma(double g1, double g2, double s1, double s2);

    /*! \brief Extract a map of combination rules for each parameter
     * \param[in] vdw_comb Old style string designating the combination rule
     * \param[in] ftype    Gromacs function type
     * \return the map
     */
    std::map<const std::string, CombRule> oldCombinationRule(const std::string &vdw_comb,
                                                             int                ftype);
    /*! \brief Extract a map of combination rules for each parameter
     * \param[in] vdw Van der Waals list of ff params
     * \return the map
     */
    std::map<const std::string, CombRule> getCombinationRule(const ForceFieldParameterList &vdw);
    
    /*! \brief Generate combined force field parameter map
     * \param[in] combrule Map of combination rules per parameter
     * \param[in] ivdw     Parameters for particle i
     * \param[in] jvdw     Parameters for particle j
     * \return a Force Field Parameter Map with pair entries
     */
    ForceFieldParameterMap evalCombinationRule(const std::map<const std::string, CombRule> &combrule,
                                               const ForceFieldParameterMap                &ivdw,
                                               const ForceFieldParameterMap                &jvdw);

    /*! \brief Generate nonbonded parameters for pairs of atoms
     * as well as force constants force shells.
     * \param[inout] pd The force field structure
     */
    void generateDependentParameter(ForceField *pd);

} // namespace alexandria

#endif
