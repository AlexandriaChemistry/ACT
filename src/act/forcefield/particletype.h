/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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
 */
#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H

#include <cstdlib>

#include <map>
#include <string>

#include "act/basics/act_particle.h"
#include "act/basics/identifier.h"
#include "act/forcefield/forcefield_parameter.h"

namespace alexandria
{

/*! \brief Produce a string. Will crash if there is no corresponding string.
 * \param[in] itype The interaction type
 * \return a topology option corresponding to the input. 
 */
const std::string &interactionTypeToParticleSubtype(InteractionType itype);

/*! \brief
 * Contains all the information realted to
 * alexandria force field partice types.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class ParticleType
{
 public:

    //! Empty constructor
    ParticleType() {}
        
    /*! \brief
     * ParticleType constructor
     *
     * \param[in] id              Identifier
     * \param[in] desc            Description
     * \param[in] gmxParticleType GROMACS particle type
     */
    ParticleType(const Identifier  &id,
                 const std::string &desc,
                 ActParticle        apType) :
    id_(id), desc_(desc), apType_(apType) {}

    /*! \brief Return the identifier
     */
    const Identifier &id() const { return id_; }

    /*! \brief Return the decription of particle
     */
    const std::string &description() const { return desc_; }

    /*! \brief Return the type of particle
     */
    ActParticle apType() const { return apType_; }

    /*! \brief Set optional values
     * 
     * Optional key/value pairs that are interpreted by the routine
     * \param[in] key   Key
     * \param[in] value Value
     * \throws if the option is not supported
     */
    void setOption(const std::string &key,
                   const std::string &value);
    
    /*! \brief Return true if option is present
     * \param[in] type The key value of the option
     * \return true or false
     */
    bool hasOption(const std::string &type) const
    {
        return option_.find(type) != option_.end();
    }
    
    /*! \brief Get option value
     * 
     * \param[in] type   Key
     * \throws if option does not exist
     * \return value
     */
    const std::string &optionValue(const std::string &type) const;
    
    /*! \brief Return the option map
     * \return All the options in a map data structure
     */
    const std::map<std::string,std::string> optionsConst() const { return option_; }
    
    /*! \brief Return whether the interaction type is supported
     */
    bool hasInteractionType(InteractionType itype) const;
    /*! \brief Get identifier corresponding to interactiontype
     * 
     * \param[in] itype The interaction type to look for
     * \throws if option does not exist for this particle.
     */
    Identifier interactionTypeToIdentifier(InteractionType itype) const;
    
    /*! \brief Add a force field parameter to the atype
     * \param[in] type  The type of parameter
     * \param[in] param The value etc. in a param structure
     */
    void addForceFieldParameter(const std::string         &type,
                                const ForceFieldParameter &param)
    {
        parameterMap_.insert({type, param});
        if (type.compare("charge") != 0)
        {
            parameterMap_.find(type)->second.setNonNegative();
        }
    }
    
    /*! \brief Return all parameters
     * \return Editable map of parameters
     */
    ForceFieldParameterMap *parameters() { return &parameterMap_; }

    /*! \brief Return all parameters
     * \return Constant map of parameters
     */
    const ForceFieldParameterMap parametersConst() const { return parameterMap_; }

    /*! \brief Return true when this parameter exists
     * \param[in] type The parameter name
     * \return whether the type is present
     */
    bool hasParameter(const std::string &type) const
    {
        return parameterMap_.find(type) != parameterMap_.end();
    }
    /*! \brief Return parameter value
     * \param[in] type The parameter name
     * \return the value
     * \throws if not present
     */
    double paramValue(const std::string &type) const;
    
    //! Return force field parameter corresponding to type
    const ForceFieldParameter &parameterConst(const std::string &type) const;

    //! Return mutable force field parameter corresponding to type
    ForceFieldParameter *parameter(const std::string &type);

    /*! \brief Return mass
     * \return mass of the particle if present, zero otherwise
     */
    double mass() const;
        
    /*! \brief Return charge
     * \return charge of the particle if present, zero otherwise
     */
    double charge() const;
        
    /*! \brief Return atomnumber
     * \return atomnumber of the particle or 0 if not available
     */
    int atomnumber() const;
    
    /*! \brief Return row in periodic table
     * \return row of the particle or 0 if not available
     */
    int row() const;
    
    /*! \brief Return element if present
     * \return element or empty string
     */
    std::string element() const;

    /*! \brief Distribute my data
     * \param[in] cr   communication record
     * \param[in] dest Destination processor
     */
    CommunicationStatus Send(const CommunicationRecord *cr, int dest);
    
    /*! \brief Broadcast my data
     * \param[in] cr   communication record
     * \param[in] root The MPI root
     * \param[in] comm MPI communicator
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);

    /*! \brief Receive my data
     * \param[in] cr  communication record
     * \param[in] src Source processor
     */
    CommunicationStatus Receive(const CommunicationRecord *cr, int src);

 private:
    //! My identifier
    Identifier                         id_;
    //! String with particle type description
    std::string                        desc_;
    //! Particle type (Atom,  etc.)
    ActParticle                        apType_;
    //! Map of options including sub types
    std::map<std::string, std::string> option_;
    //! The force field parameters associated with this particle type
    ForceFieldParameterMap             parameterMap_;
};

} // namespace

#endif
