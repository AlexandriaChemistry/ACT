/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020 
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

#ifndef FORCEFIELDPARAMETERLIST_H
#define FORCEFIELDPARAMETERLIST_H

#include <map>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "forcefieldparameter.h"

namespace alexandria
{

/*! \brief Class to hold the parameters for an interaction type
 *
 * This class hold a list of parameters of a certain type for a
 * force field. It is important that the class does not hold any
 * implicit information.
 */
class ForceFieldParameterList
{
 public:
    /*! \brief Constructor
     *
     * \param[in] function The function for which parameters are stored
     */
    ForceFieldParameterList(const std::string &function) : function_(function) {}
   
    //! \brief Return the function name
    const std::string &function() const { return function_; }

    /*! \brief Add function specific options
     *
     * \param[in] option  Name of the option
     * \param[in] value   Value of the opton
     */
    void addOption(const std::string &option, const std::string &value)
    {
        options_.insert({option, value});
    }
    
    /*! \brief check whether an option exists
     *
     * \param[in] option  Name of the option
     * \return bool value 
     */
    bool optionExists(const std::string &option) const
    {
        return options_.find(option) != options_.end();
    }
    /*! \brief Extract the value corresponding to an option
     *
     * \param[in] option Name of the option
     * \return Value
     * \throws gmx::InvalidInputError if the option is non-existent
     */
    const std::string &optionValue(const std::string option) const
    {
        auto ov = options_.find(option);
        if (ov == options_.end())
        {
            auto buf = gmx::formatString("Unknown option %s", option.c_str());
            GMX_THROW(gmx::InvalidInputError(buf.c_str()));
        }
        return ov->second;
    }

    /*! \brief Add one parameter to the map
     *
     * Add one parameter for the given identifier (e.g. atomtype or bond).
     * If the parameter type has been set already an exception will be raised 
     * (this check has not been implemented yet).
     * \param[in] identifier Name of the atomtype or bond
     * \param[in] param      The force field parameter structure
     */
    void addParameter(const std::string         &identifier,
                      const ForceFieldParameter &param);
       
    
    /*! \brief check whether any parameter exists for identifier
     *
     * \param[in] identifier  Name of the bond or atomtype
     * \return bool value 
     */
    bool parameterExists(const std::string &identifier) const
    {
        return parameters_.find(identifier) != parameters_.end();
    }

    /*! \brief Look up vector of parameters
     * Will throw an exception when identifier is not found
     * \param[in] identifier Name of the atomtype or bond
     * \return vector of parameters
     */
    const std::vector<ForceFieldParameter> &searchParameter(const std::string &identifier) const;
    
    /*! \brief Find vector of parameters for editing
     * Will throw an exception when identifier is not found
     * \param[in] identifier Name of the atomtype or bond
     * \return vector of parameters
     */
    std::vector<ForceFieldParameter> &findParameter(const std::string &identifier);
    
 private:
    //! The function type
    std::string function_;

    //! Map structure for the options associated with the parameter list
    std::map<const std::string, const std::string>   options_;
        
    //! Map of parameters 
    std::map<const std::string, std::vector<ForceFieldParameter> > parameters_;
};

} // namespace

#endif
