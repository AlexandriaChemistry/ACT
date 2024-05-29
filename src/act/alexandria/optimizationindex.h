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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef ALEXANDRIA_OPTIMIZATIONINDEX_H
#define ALEXANDRIA_OPTIMIZATIONINDEX_H

#include <map>
#include <set>
#include <string>
#include <vector>

#include "act/ga/genome.h"
#include "act/utility/communicationrecord.h"
#include "act/basics/mutability.h"
#include "molgen.h"

namespace alexandria
{
class ForceFieldParameter;
/*! \brief Convenience storage of parameters to optimize
 */
class OptimizationIndex
{
private:
    //! The particle type. If empty, an interaction type is used.
    std::string          particleType_;
    //! The interaction type
    InteractionType      iType_ = InteractionType::CHARGE;
    //! The name of the parameter matching forcefield
    Identifier           parameterId_;
    //! The type of parameter, eg. sigma, epsilon
    std::string          parameterType_;
    //! Pointer to relevant force field parameter
    ForceFieldParameter *ffp_ = nullptr;
public:
    //! Default constructor
    OptimizationIndex() {}

    /*! \brief Constructor
     * \param[in] iType         The interaction type
     * \param[in] parameterId   The identifier
     * \param[in] parameterType The type
     */
    OptimizationIndex(InteractionType    iType,
                      Identifier         parameterId,
                      const std::string &parameterType) :
        iType_(iType), parameterId_(parameterId), parameterType_(parameterType) {}

    /*! \brief Constructor
     * \param[in] pType         The particle type
     * \param[in] parameterType The type
     */
    OptimizationIndex(const std::string  pType,
                      const std::string &parameterType) :
        particleType_(pType), parameterType_(parameterType) {}

    /*! Lookup the force field parameter
     * \param[in] pd The force field
     */
    void findForceFieldParameter(ForceField *pd);

    //! \return The force field parameter
    ForceFieldParameter *forceFieldParameter() { return ffp_; }

    //! \return the interaction type
    InteractionType iType() const { return iType_; }

    //! \return the id
    Identifier id() const { return parameterId_; }

    //! \return particle type
    const std::string &particleType() const { return particleType_; }

    //! \return the type
    const std::string &parameterType() const { return parameterType_; }

    //! \return a compound string representing the index
    std::string name() const;

    //! \return a compound string representing the atomindex
    std::string parameterName() const
    {
        return parameterId_.id().c_str();
    }

    /*! \brief Send an OptimizationIndex
     * \param[in] cr   The communication information
     * \param[in] dest The destination processor
     */
    CommunicationStatus send(const CommunicationRecord *cr,
                             int                        dest);

    /*! \brief Send an OptimizationIndex
     * \param[in] cr   The communication information
     * \param[in] dest The destination processor
     */
    CommunicationStatus receive(const CommunicationRecord *cr,
                                int                        src);
};

} // namespace

#endif
