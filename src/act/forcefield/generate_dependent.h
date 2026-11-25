/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2025,2026
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
 
/*! \internal \brief
 * Implements generation of parameter matrices using combination rules.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#ifndef FFUTIL_H
#define FFUTIL_H

#include <set>
#include <string>
#include <vector>

#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"

namespace alexandria
{

class ForceField;
class MsgHandler;

//! The name of the force constant of a shell in the polarization parameters
extern const char *kshellName;

/*! \brief A shell whose force constant is fixed in the force field
 *
 * The force constant of a shell is \f$k = q_s^2/(4\pi\varepsilon_0\alpha)\f$
 * for shell charge \f$q_s\f$ and polarizability \f$\alpha\f$. Normally it
 * follows from these two. If the polarization parameters of the shell type
 * have a parameter kshell that is not Dependent, the force constant is
 * given instead, and the one of alpha and the shell charge that is
 * Dependent follows from it and the other one.
 */
struct ShellForceConstant
{
    //! The identifier of the polarization parameters, the particle type of the shell
    std::string shellType;
    //! The force constant, in the unit of the force field
    double      kshell         = 0;
    //! Its unit
    std::string kshellUnit;
    //! The polarizability, in the unit of the force field
    double      alpha          = 0;
    //! Its unit
    std::string alphaUnit;
    //! The charge of the shell (e)
    double      charge         = 0;
    //! Whether alpha follows from the force constant, otherwise the charge does
    bool        alphaDependent = true;
};

/*! \brief The shells whose force constant is fixed, with their parameters
 * \param[in] pd The force field
 * \return The shells, sorted by type
 */
std::vector<ShellForceConstant> fixedShellForceConstants(const ForceField *pd);

//! What generateShellParameters changed
struct ShellParametersChanged
{
    //! A polarizability
    bool alpha  = false;
    //! A charge of a shell
    bool charge = false;
};

/*! \brief Compute the parameters of shells that follow from a fixed force constant
 *
 * For each shell type with a fixed force constant (see ShellForceConstant),
 * either \f$\alpha = q_s^2/(4\pi\varepsilon_0 k)\f$, if alpha is Dependent,
 * or \f$q_s = \pm\sqrt{4\pi\varepsilon_0 k\alpha}\f$, if the charge of the
 * shell particle type is, keeping the sign of the charge (negative if it
 * is zero).
 * \param[inout] pd The force field
 * \return What changed
 * \throws gmx::InvalidInputError if the force constant is not positive, or
 *         if not exactly one of alpha and the shell charge is Dependent
 */
ShellParametersChanged generateShellParameters(ForceField *pd);

//! Which parameter of a shell follows from its fixed force constant
enum class ShellFollows
{
    //! The polarizability
    Alpha,
    //! The charge of the shell
    Charge
};

/*! \brief The shell types of a force field
 * \param[in] pd The force field
 * \return The identifiers of the polarization parameters with alpha that
 *         are particle types of shells, sorted
 */
std::vector<Identifier> shellTypes(const ForceField *pd);

/*! \brief Fix the force constant of a shell type
 *
 * The parameter kshell becomes Fixed, the one that follows Dependent,
 * and the other one, if it was Dependent, Fixed at its value. The
 * dependent parameter is computed, see generateShellParameters.
 * \param[inout] pd        The force field
 * \param[in]    shellType A shell type, see shellTypes
 * \param[in]    follows   Which parameter follows from the force constant
 * \param[in]    kshell    The force constant in kJ/mol nm2, or, if it is
 *                         not above zero, that of the current
 *                         polarizability and charge
 * \throws gmx::InvalidInputError if shellType is not a shell type, or if
 *         the current force constant is not defined (alpha or the charge
 *         is zero)
 */
void fixShellForceConstant(ForceField       *pd,
                           const Identifier &shellType,
                           ShellFollows      follows,
                           double            kshell);

/*! \brief Let the force constant of a shell type follow from alpha and the charge again
 *
 * The parameter kshell, if there is one, becomes Dependent with the force
 * constant of the current values, and a dependent polarizability or charge
 * becomes Fixed at its value.
 * \param[inout] pd        The force field
 * \param[in]    shellType A shell type, see shellTypes
 * \return Whether the force constant was fixed
 * \throws gmx::InvalidInputError if shellType is not a shell type
 */
bool releaseShellForceConstant(ForceField       *pd,
                               const Identifier &shellType);

/*! \brief Generate nonbonded parameters for pairs of atoms
 * as well as the shell parameters that follow from fixed force
 * constants, see generateShellParameters. Only pair parameters
 * for which one of the atomic parameters were updated will be
 * generated, unless force is true.
 * \param[in]    msghandler For debugging info
 * \param[inout] pd         The force field structure
 * \param[in]    force      Update all parameters in the matrices
 * \param[in]    changed    If not nullptr, only the parameters of these
 *                          interaction types may have been updated, such
 *                          that the other ones can be skipped
 */
void generateDependentParameter(MsgHandler                      *msghandler,
                                ForceField                      *pd,
                                bool                             force   = false,
                                const std::set<InteractionType> *changed = nullptr);

}
#endif
