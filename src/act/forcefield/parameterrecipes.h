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
 * Recipes for the parameters of interactions in the units of the program,
 * from the values of the force field.
 *
 * The parameters of an interaction, as the topology has them, follow from
 * the values of the force field in their own units: either directly, a
 * value times the factor of its unit, see ForceFieldParameter::internalValue,
 * or for pairs of particle types by combining the values of both types with
 * the combination rules, see evalCombinationRule. A recipe says how, for
 * one set of parameters, such as those of a bond type or of a pair of van
 * der Waals types, in terms of the positions of the values in one array.
 * ParameterRecipes makes the recipes and the array, and the parameters can
 * then be computed from the array alone, also on the GPU, with
 * evaluateRecipe, which only uses the functions of combinationkernels.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_forcefield
 */
#ifndef ACT_FORCEFIELD_PARAMETERRECIPES_H
#define ACT_FORCEFIELD_PARAMETERRECIPES_H

#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/combinationkernels.h"
#include "act/forcefield/setrecipe.h"

namespace alexandria
{

class ForceField;
class ForceFieldParameter;

/*! \brief The recipes of the parameters of interactions, with the values they use
 * The recipes point to the values of the force field they were made from,
 * which update() copies again, such that the force field must not change
 * its structure afterwards, only its values.
 */
class ParameterRecipes
{
public:
    /*! \brief The set of parameters of an entry of the force field, such as a bond type
     * \param[in] pd    The force field
     * \param[in] itype The interaction type
     * \param[in] id    The identifier of the entry
     * \return The index of the set, -1 if the force field has no such entry
     */
    int directSet(const ForceField *pd, InteractionType itype, const Identifier &id);

    /*! \brief The set of parameters of a pair of particle types, from the combination rules
     * \param[in] pd    The force field
     * \param[in] itype The interaction type, such as InteractionType::VDW
     * \param[in] pair  The identifier of the pair, with two particle types
     * \return The index of the set, -1 if the force field has not both types
     */
    int pairSet(const ForceField *pd, InteractionType itype, const Identifier &pair);

    //! Copy the values of the force field again, after they changed
    void update();

    //! \return the recipes
    const std::vector<SetRecipe> &recipes() const { return recipes_; }

    //! \return the values, in their own units
    const std::vector<double> &values() const { return values_; }

    //! \return the factors of the units of the values
    const std::vector<double> &factors() const { return factors_; }

    /*! \brief Compute the parameters of a set on the processor
     * \param[in] set The index of the set
     * \return The parameters, in the units of the program
     * \throws if the combination rules cannot be applied
     */
    std::vector<double> evaluate(int set) const;

    /*! \brief The position of a value
     * The value is copied by update as those of the recipes, for instance
     * for the charges of the particles, see gpu::ChargeRecipe.
     * \param[in] p The parameter of the force field
     * \return The position in values and factors, added if new
     */
    int32_t slot(const ForceFieldParameter *p);

    /*! \brief The position of a value, if it has one
     * \param[in] p The parameter of the force field
     * \return The position in values, -1 if none
     */
    int32_t findSlot(const ForceFieldParameter *p) const;

private:

    //! The recipes
    std::vector<SetRecipe>                                     recipes_;
    //! The values of the force field, in their own units
    std::vector<double>                                        values_;
    //! The factors of their units
    std::vector<double>                                        factors_;
    //! The parameters of the force field of the values
    std::vector<const ForceFieldParameter *>                   params_;
    //! The position of each parameter
    std::map<const ForceFieldParameter *, int32_t>             slots_;
    //! The set of each entry of the force field
    std::map<std::pair<InteractionType, Identifier>, int>     direct_;
    //! The set of each pair of particle types
    std::map<std::pair<InteractionType, Identifier>, int>     pairs_;
};

} // namespace alexandria

#endif
