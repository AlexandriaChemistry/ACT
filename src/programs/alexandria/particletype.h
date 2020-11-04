/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H

#include <cstdlib>

#include <map>
#include <string>

#include "forcefieldparameter.h"
#include "identifier.h"

namespace alexandria
{

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
                 int               gmxParticleType) :
    id_(id), desc_(desc), gmxParticleType_(gmxParticleType) {}
    
    /*! \brief Return the identifier
     */
    const Identifier &id() const { return id_; }

    /*! \brief Return the decription of particle
     */
    const std::string &description() const { return desc_; }

    /*! \brief Return the type of particle
     */
    int gmxParticleType() const { return gmxParticleType_; }

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
     * \param[in] key   Key
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
     * \param[in] key itype
     * \return value will be empty if option does not exist an
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
    }
    
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
    
    //! Return pointer to force field parameter corresponding to type
    const ForceFieldParameter &parameter(const std::string &type) const;
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
    
    /*! \brief Return reference enthalpy for this particle
     * \return ref_enthalpy of the particle or 0 if not available
     */
    double refEnthalpy() const;
    
    /*! \brief Return element if present
     * \return element or empty string
     */
    std::string element() const;

    /*! \brief Distribute my data
     * \param[in] cr   GROMACS communication record
     * \param[in] dest Destination processor
     */
    CommunicationStatus Send(const t_commrec *cr, int dest);
    
    /*! \brief Receive my data
     * \param[in] cr  GROMACS communication record
     * \param[in] src Source processor
     */
    CommunicationStatus Receive(const t_commrec *cr, int src);

 private:
    //! My identifier
    Identifier                         id_;
    //! String with particle type description
    std::string                        desc_;
    //! GROMACS particle type (eptAtom, eptShell etc.)
    int                                gmxParticleType_;
    //! Map of options including sub types
    std::map<std::string, std::string> option_;
    //! The force field parameters associated with this particle type
    ForceFieldParameterMap             parameterMap_;
};

//! Iterator over vector of these
using ParticleTypeIterator      = std::vector<ParticleType>::iterator;

//! Const iterator over vector of these
using ParticleTypeConstIterator = std::vector<ParticleType>::const_iterator;

} // namespace

#endif
