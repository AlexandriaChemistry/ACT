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
 
/*! \brief Class to hold one force field parameter
 */
class ForceFieldParameter
{
 public:
    //! \brief default constructor
    ForceFieldParameter() {}
    
    /*! \brief Constructor initiating all parameters.
     *
     * \param[in] name        Description of the parameter
     * \param[in] index       Index in a parameter list
     * \param[in] value       Actual value of the parameter
     * \param[in] uncertainty Uncertainty in the value
     * \param[in] minimum     Minimum allowed value
     * \param[in] maximum     Maximum allowed value
     * \param[in] mutability  In what way this parameter may be changed
     * \param[in] strict      Throw an exception in case of value errors
     */
    ForceFieldParameter(const std::string &name,
                        int        index,
                        double     value,
                        double     uncertainty,
                        double     minimum,
                        double     maximum,
                        Mutability mutability,
                        bool       strict) : 
    name_(name), index_(index), value_(value), uncertainty_(uncertainty),
        minimum_(minimum), maximum_(maximum), mutability_(mutability), strict_(strict) {}
        
    //! \brief Return parameter name
    const std::string &name() const { return name_; }
    
    //! \brief Return parameter index
    int index() const { return index_; }
    
    //! \brief Return parameter value
    double value() const { return value_; }
    
    /*! \brief Set the parameter value if not fixed
     *
     * \param[in] value The new value
     * Will throw an exception if the strict flag is true and the value
     * is out of range, or the variable is alltogher fixed,
     */
    void setValue(double value);
    
    //! \brief Return the uncertainty in this value
    double uncertainty() const { return uncertainty_; }
    
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
    //! The name of the parameter
    std::string name_;
    int         index_       = 0;
    //! The index of the parameter
    double      value_       = 0;
    double      uncertainty_ = 0;
    double      minimum_     = 0;
    double      maximum_     = 0;
    Mutability  mutability_  = Mutability::Free;
    bool        strict_      = false;
};

} // namespace alexandria

#endif
