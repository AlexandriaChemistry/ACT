
    /*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

#include "act/forcefield/forcefield.h"

namespace alexandria
{
/*! \brief Combine epsilon and sigma into c6 and c12 for LJ
 * \param[in]  CombinationRule The combination rule used
 * \param[in]  sigmaI          The first sigma
 * \param[in]  sigmaJ          The second sigma
 * \param[in]  epsilonI        The first epsilon
 * \param[in]  epsilonJ        The second epsilon
 * \param[out] c6              The LJ c6
 * \param[out] c12             The LJ c12
 */
void CombineLJ(int     CombinationRule,
               double  sigmaI,
               double  sigmaJ,
               double  epsilonI,
               double  epsilonJ,
               double *c6,
               double *c12);

void CombineLJ_86(int     CombinationRule,
               double  sigmaI,
               double  sigmaJ,
               double  epsilonI,
               double  epsilonJ,
               double *c6,
               double *c8);


/*! \brief Extract the GROMACS style combination rule
 * \param[in] vdw Van der Waals list of ff params
 * \return index in GROMACS list of comb rules
 */
int getCombinationRule(const ForceFieldParameterList &vdw);

/*! \brief Combine epsilon, sigma and gamm.
 * \param[in]  CombinationRule The combination rule used
 * \param[in]  sigmaI          The first sigma
 * \param[in]  sigmaJ          The second sigma
 * \param[in]  epsilonI        The first epsilon
 * \param[in]  epsilonJ        The second epsilon
 * \param[in]  gammaI          The first gamma
 * \param[in]  gammaJ          The second gamma
 * \param[out] sigmaIJ         The combined sigma
 * \param[out] epsilonIJ       The combined epsilon
 * \param[out] gammaIJ         The combined gamma
 */
void CombineBham(int     CombinationRule,
                 double  sigmaI,
                 double  sigmaJ,
                 double  epsilonI,
                 double  epsilonJ,
                 double  gammaI,
                 double  gammaJ,
                 double *sigmaIJ,
                 double *epsilonIJ,
                 double *gammaIJ);
void CombineLJ_147(int     CombinationRule,
                 double  sigmaI,
                 double  sigmaJ,
                 double  epsilonI,
                 double  epsilonJ,
                 double  gammaI,
                 double  gammaJ,
		 double  deltaI,
		 double  deltaJ,
                 double *sigmaIJ,
                 double *epsilonIJ,
                 double *gammaIJ,
		 double *deltaIJ);

void CombineGBham(int     CombinationRule,
                  double  rminI,
                  double  rminJ,
                  double  epsilonI,
                  double  epsilonJ,
                  double  gammaI,
                  double  gammaJ,
                  double  deltaI,
                  double  deltaJ,
                  double *rminIJ,
                  double *epsilonIJ,
                  double *gammaIJ,
                  double *deltaIJ);

/*! \brief Generate nonbonded parameters for pairs of atoms
 * as well as force constants force shells.
 * \param[inout] pd The force field structure
 */
void generateDependentParameter(ForceField *pd);

} // namespace alexandria

#endif
