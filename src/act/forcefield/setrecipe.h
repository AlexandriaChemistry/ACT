/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
/*! \internal \file
 * \brief
 * A recipe for the parameters of an interaction, for the processor and the GPU.
 *
 * See ParameterRecipes in parameterrecipes.h, which makes them. This header
 * only depends on combinationkernels.h, such that GPU code can include it.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_forcefield
 */
#ifndef ACT_FORCEFIELD_SETRECIPE_H
#define ACT_FORCEFIELD_SETRECIPE_H

#include <cstdint>

#include "act/forcefield/combinationkernels.h"

namespace alexandria
{

//! The largest number of parameters of a set, see SetRecipe
constexpr int maxSetParameters = 8;

/*! \brief How to compute one set of parameters from the values
 * The parameters are in the order of the potential, see
 * potentialToParameterName, as the topology has them. The type is plain
 * data, such that the recipes can be copied to the GPU.
 */
struct SetRecipe
{
    //! Whether the parameters are those of a pair of particle types
    uint8_t                      pair    = 0;
    //! The number of parameters
    uint8_t                      n       = 0;
    //! For a pair, the number of values of each particle type
    uint8_t                      nin     = 0;
    //! For a pair, for Kronecker, whether the pair interacts
    uint8_t                      include = 1;
    //! For a pair, the positions of sigma, epsilon and gamma among the values
    combination::PairCombination pc;
    //! For each parameter, or for a pair each value of the first type, its position in the values, -1 if missing
    int32_t                      slotI[maxSetParameters];
    //! For a pair, the positions of the values of the second type
    int32_t                      slotJ[maxSetParameters];
    //! For a pair, the position of the exponent of the rule of each value, -1 if there is none
    int32_t                      exponent[maxSetParameters];
    //! For a pair, the rule of each value
    CombRule                     rules[maxSetParameters];
    //! For a pair, the value that each parameter is, -1 if missing
    int8_t                       order[maxSetParameters];
};

/*! \brief Compute the parameters of a set from the values
 * \tparam M The mathematical functions, see combinationkernels.h
 * \param[in]  recipe  The recipe
 * \param[in]  values  The values of the force field, in their own units
 * \param[in]  factors The factor of the unit of each value
 * \param[out] out     The parameters, recipe.n of them
 * \return The status of combining a pair, Ok for others
 */
template <class M>
inline combination::PairStatus evaluateRecipe(const SetRecipe &recipe,
                                              const double    *values,
                                              const double    *factors,
                                              double          *out)
{
    const auto internal = [values, factors](int32_t slot)
    {
        return slot < 0 ? 0.0 : values[slot]*factors[slot];
    };
    if (!recipe.pair)
    {
        for (int k = 0; k < recipe.n; k++)
        {
            out[k] = internal(recipe.slotI[k]);
        }
        return combination::PairStatus::Ok;
    }
    // Only the first recipe.nin are used, the others are set for compilers
    // that cannot see that, such as GCC with -Wmaybe-uninitialized
    double vi[maxSetParameters]       = { 0 };
    double vj[maxSetParameters]       = { 0 };
    double exps[maxSetParameters]     = { 0 };
    double combined[maxSetParameters] = { 0 };
    for (int q = 0; q < recipe.nin; q++)
    {
        vi[q]   = internal(recipe.slotI[q]);
        vj[q]   = internal(recipe.slotJ[q]);
        exps[q] = recipe.exponent[q] < 0 ? 0.0 : values[recipe.exponent[q]];
    }
    int  bad    = 0;
    auto status = combination::pairParameters<M>(recipe.pc, recipe.nin, recipe.rules, exps,
                                                 vi, vj, recipe.include != 0, combined, &bad);
    for (int k = 0; k < recipe.n; k++)
    {
        out[k] = recipe.order[k] < 0 ? 0.0 : combined[recipe.order[k]];
    }
    return status;
}

} // namespace alexandria

#endif
