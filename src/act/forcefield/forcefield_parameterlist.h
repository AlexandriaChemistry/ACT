/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2025
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef FORCEFIELDPARAMETERLIST_H
#define FORCEFIELDPARAMETERLIST_H

#include <algorithm>
#include <map>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/potential.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"

namespace alexandria
{

//! \brief Shortcut for recurring declaration
typedef std::map<Identifier, ForceFieldParameterMap> ForceFieldParameterListMap;

/*! \brief Class to hold the parameters for an interaction type
 *
 * This class hold a list of parameters of a certain type for a
 * force field. It is important that the class does not hold any
 * implicit information.
 */
class ForceFieldParameterList
{
 public:
    //! Empty constructor for helper nodes
    ForceFieldParameterList() {};
    
    /*! \brief Constructor
     *
     * \param[in] function The function for which parameters are stored. This may be an empty variable.
     * \param[in] canSwap  Can atom/bond identifiers be swapped in order
     * \throw if function is not recognized.
     */
    ForceFieldParameterList(const std::string &function,
                            CanSwap            canSwap);
   
    //! \brief Return the function name
    Potential potential() const { return pot_; }

    //! \brief Set the function name
    void setFunction(const std::string &function);

    //! \brief Return whether or not identifiers can be swapped
    CanSwap canSwap() const { return canSwap_; }

    //! \brief Set whether or not identifiers can be swapped
    void setCanSwap(CanSwap cs) { canSwap_ = cs; }

    /*! \brief Add function specific options
     *
     * \param[in] option  Name of the option
     * \param[in] value   Value of the opton
     */
    void addOption(const std::string &option, const std::string &value);
    
    /*! \brief check whether an option exists
     *
     * \param[in] option  Name of the option
     * \return bool value 
     */
    bool optionExists(const std::string &option) const
    {
        return options_.find(option) != options_.end();
    }

    /*! \brief Remove an option if it exists
     *
     * \param[in] option  Name of the option
     */
    void removeOption(const std::string &option)
    {
        auto myopt = options_.find(option);
        if (myopt != options_.end())
        {
            options_.erase(myopt);
        }
    }

    /*! \brief Extract the value corresponding to an option
     *
     * \param[in] option Name of the option
     * \return Value
     * \throws gmx::InvalidInputError if the option is non-existent
     */
    const std::string &optionValue(const std::string &option) const
    {
        auto ov = options_.find(option);
        if (ov == options_.end())
        {
            auto buf = gmx::formatString("Unknown option %s", option.c_str());
            GMX_THROW(gmx::InvalidInputError(buf.c_str()));
        }
        return ov->second;
    }

    //! Return the options map
    const std::map<std::string, std::string> &option() const { return options_; }

    /*! \brief Add function specific combination rule
     *
     * \param[in] param  Name of the parameter
     * \param[in] rule   String for the rule
     */
    void addCombinationRule(const std::string &param, const std::string &rule)
    {
        combrules_[param] = rule;
    }
    
    /*! \brief check whether a combination rule exists
     *
     * \param[in] param  Name of the parameter
     * \return bool value 
     */
    bool combinationRuleExists(const std::string &param) const
    {
        return combrules_.find(param) != combrules_.end();
    }

    /*! \brief Extract the rule corresponding to an parameter
     *
     * \param[in] param Name of the parameter
     * \return Value
     * \throws gmx::InvalidInputError if the param is non-existent
     */
    const std::string &combinationRule(const std::string &param) const
    {
        auto ov = combrules_.find(param);
        if (ov == combrules_.end())
        {
            auto buf = gmx::formatString("Unknown parameter %s", param.c_str());
            GMX_THROW(gmx::InvalidInputError(buf.c_str()));
        }
        return ov->second;
    }

    //! Return the combination rule map
    const std::map<std::string, std::string> &combinationRules() const { return combrules_; }

    //! Return the parameters map as a const variable
    const ForceFieldParameterListMap &parametersConst() const { return parameters_; };
    
    //! Return the parameters map as an editable variable
    ForceFieldParameterListMap *parameters() { return &parameters_; };
    
    //! Clear the parameter map
    void clearParameters() { parameters_.clear(); }

    /*! \brief Add one parameter to the map
     *
     * Add one parameter for the given identifier (e.g. atomtype or bond).
     * If the parameter type has been set already an exception will be raised 
     * (this check has not been implemented yet).
     * \param[in] identifier Name of the atomtype or bond
     * \param[in] type       The parameter type
     * \param[in] param      The force field parameter structure
     */
    void addParameter(const Identifier          &identifier,
                      const std::string         &type,
                      const ForceFieldParameter &param);
       
    
    /*! \brief check whether any parameter exists for identifier
     *
     * \param[in] identifier  Name of the bond or atomtype
     * \return bool value 
     */
    bool parameterExists(const Identifier &identifier) const;

    /*! \brief Get number of parameters in this list
     * \return number of parameters
     */
    size_t numberOfParameters() const { return parameters_.size(); }

    //! \brief return whether or not there are any parameters
    bool empty() const { return parameters_.size() == 0; }
    
    /*! \brief check whether any parameter exists for identifier
     *
     * \param[in] identifier  Name of the bond or atomtype
     * \return Parameter index
     */
    size_t parameterId(const Identifier &identifier) const;

    /*! \brief Look up map of parameters
     * Will throw an exception when identifier is not found
     * \param[in] identifier Name of the atomtype or bond
     * \return map of parameters
     */
    const ForceFieldParameterMap &findParametersConst(const Identifier &identifier) const;
    
    /*! \brief Look up map of parameters
     * \param[in] identifier Name of the atomtype or bond
     * \return map of parameters that maybe empty if not found
     */
    const ForceFieldParameterMap &findParameterMapConst(const Identifier &identifier) const;
    
    /*! \brief Find map of parameters for editing
     * Will throw an exception when identifier is not found
     * \param[in] identifier Name of the atomtype or bond
     * \return vector of parameters
     */
    const ForceFieldParameterMap *findParametersPtrConst(const Identifier &identifier) const;

    /*! \brief Find map of parameters for editing
     * Will throw an exception when identifier is not found
     * \param[in] identifier Name of the atomtype or bond
     * \return vector of parameters
     */
    ForceFieldParameterMap *findParameters(const Identifier &identifier);

    /*! \brief Convenience function to look up specific parameter
     *
     * \param[in] identifier String with atoms or bonds
     * \param[in] type       The parameter type, e.g. sigma or bondlength
     * \throws if not found
     */
    const ForceFieldParameter &findParameterTypeConst(const Identifier  &identifier,
                                                      const std::string &type) const;

    /*! \brief Convenience function to look up specific parameter for editing
     *
     * \param[in] identifier String with atoms or bonds
     * \param[in] type       The parameter type, e.g. sigma or bondlength
     * \throws if not found
     */
    ForceFieldParameter *findParameterType(const Identifier  &identifier,
                                           const std::string &type);

    /*! \brief
     * Clear the parameter map. Used in bastat to rebuild the parameters.
     */
    void eraseParameter()
    {
        parameters_.clear();
        counter_ = 0;
    }

    /*! \brief Dump contents to a file
     * \param[in] fp File pointer
     */
    void dump(FILE *fp) const;
    
    /*! \brief Send the contents to another processor
     * \param[in] cr   Communication data structure
     * \param[in] dest Processor id to send the data to
     */
    CommunicationStatus Send(const CommunicationRecord *cr, int dest) const;

    /*! \brief Broadcast contents from another processor
     * \param[in] cr   Communication data structure
     * \param[in] root The MPI root
     * \param[in] comm MPI communication structure
     * \return The status of the whole thing
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);

    /*! \brief Receive contents from another processor
     * \param[in] cr  Communication data structure
     * \param[in] src Processor id to receive the data from
     */
    CommunicationStatus Receive(const CommunicationRecord *cr, int src);

    //! \return The counter \p counter_ for index
    int counter() const { return counter_; };

 private:
    //! The potential
    Potential    pot_ = Potential::NONE;

    //! Whether or not swapping is allowed
    CanSwap      canSwap_;

    //! Map structure for the options associated with the parameter list
    std::map<std::string, std::string> options_;
        
    //! Map structure for combination rules associated with the parameter list
    std::map<std::string, std::string> combrules_;
        
    //! Map of parameters 
    ForceFieldParameterListMap parameters_;

    //! Counter for index
    int counter_ = 0;
};

} // namespace

#endif
