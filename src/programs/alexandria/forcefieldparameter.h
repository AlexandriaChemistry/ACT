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

#ifndef FORCEFIELDPARAMETER_H
#define FORCEFIELDPARAMETER_H

#include <algorithm>
#include <string>

namespace alexandria
{

//! \brief Enum determining whether a parameter can be changed
enum class Mutability
{
    Fixed,
    Bounded,
    Free
};

/*! \brief Return a string corresponding to the mutability
 */
const std::string &mutabilityName(Mutability mutability);

/*! \brief Lookup a string and return mutability value
 *
 * \param[in]  name       String
 * \param[out] mutability Point to Mutability encoded in the string
 * \return true if the name could be converted to a Mutability, false otherwise
 */
bool nameToMutability(const std::string &name, Mutability *mutability);

/*! \brief Class to hold one force field parameter
 */
class ForceFieldParameter
{
 public:
    //! \brief default constructor
    ForceFieldParameter() {}
    
    /*! \brief Constructor initiating all parameters.
     *
     * \param[in] identifier  Id linking the parameter to e.g. an atomtype
     * \param[in] type        Type of parameter
     * \param[in] value       Actual value of the parameter
     * \param[in] uncertainty Uncertainty in the value
     * \param[in] minimum     Minimum allowed value
     * \param[in] maximum     Maximum allowed value
     * \param[in] mutability  In what way this parameter may be changed
     * \param[in] strict      Throw an exception in case of value errors
     */
    ForceFieldParameter(const std::string &identifier,
                        const std::string &type,
                        double             value,
                        double             uncertainty,
                        double             minimum,
                        double             maximum,
                        Mutability         mutability,
                        bool               strict) : 
    identifier_(identifier), type_(type), value_(value), originalValue_ (value),
        uncertainty_(uncertainty), originalUncertainty_(uncertainty),
        minimum_(minimum), maximum_(maximum), mutability_(mutability), strict_(strict) {}
        
    //! \brief Return parameter identifier
    const std::string &identifier() const { return identifier_; }
    
    //! \brief Return type of parameter
    const std::string &type() const { return type_; }
    
    //! \brief Return current parameter value
    double value() const { return value_; }
    
    //! \brief Return original parameter value
    double originalValue() const { return originalValue_; }
    
    /*! \brief Set the parameter value if not fixed
     *
     * \param[in] value The new value
     * Will throw an exception if the strict flag is true and the value
     * is out of range, or the variable is alltogher fixed,
     */
    void setValue(double value);
    
    //! \brief Return the current uncertainty in this value
    double uncertainty() const { return uncertainty_; }
    
    //! \brief Return the original uncertainty in this value
    double originalUncertainty() const { return originalUncertainty_; }
    
    /*! \brief Set the uncertainty in this value if not fixed
     * \param[in] uncertainty  The new uncertainty
     * Will throw an exception if the strict flag is true and the
     * variable is alltogher fixed,
     */
    void setUncertainty(double uncertainty);
    
    //! \brief Return minimum allowed value
    double minimum() const { return minimum_; }
    
    //! \brief Return maximum allowed value
    double maximum() const { return maximum_; }
    
    //! \brief Return how this parameter may be changed
    Mutability mutability() const { return mutability_; }
    
    //! \brief Return whether or not to throw on value errors
    bool strict() const { return strict_; }
 private:
    //! The identifier of the parameter
    std::string identifier_;
    //! The type of the parameter
    std::string type_;
    //! The current value of the parameter
    double      value_               = 0;
    //! The original value of the parameter
    double      originalValue_       = 0;
    //! The current uncertainty in the parameter
    double      uncertainty_         = 0;
    //! The original value of the uncertainty
    double      originalUncertainty_ = 0;
    //! Minimum allowed value for the parameter
    double      minimum_             = 0;
    //! Maximum allowed value for the parameter
    double      maximum_             = 0;
    //! In what way this parameter is mutable
    Mutability  mutability_          = Mutability::Free;
    /*! Whether or not to throw an exception in case value or
     * uncertainty is set incorrectly
     */
    bool        strict_              = false;
};

} // namespace alexandria

#endif
