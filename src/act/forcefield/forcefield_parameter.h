/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2022,2026
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

#ifndef FORCEFIELDPARAMETER_H
#define FORCEFIELDPARAMETER_H

#include <algorithm>
#include <map>
#include <string>

#include "gromacs/math/units.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "act/utility/communicationrecord.h"
#include "act/utility/units.h"
#include "act/basics/mutability.h"

/*! \brief Convert string to boolean
 * \param[in] str The string to convert
 * \return boolean
 */
bool stringToBoolean(const std::string &str);

namespace alexandria
{

/*! \brief Class to hold one force field parameter
 */
class ForceFieldParameter
{
 public:
    //! \brief default constructor
    ForceFieldParameter() = default;
    
    /*! \brief Constructor initiating all parameters.
     *
     * If min and max are equal, mutability will be set to Fixed
     * automatically.
     * \todo Check unit
     * \param[in] unit        Physical unit of parameter, e.g. nm or fs
     * \param[in] value       Actual value of the parameter
     * \param[in] uncertainty Uncertainty in the value
     * \param[in] ntrain      Number of data points used for training
     * \param[in] minimum     Minimum allowed value
     * \param[in] maximum     Maximum allowed value
     * \param[in] mutability  In what way this parameter may be changed
     * \param[in] strict      Throw an exception in case of value errors
     * \param[in] nonNegative Prevent parameter from obtaining a negative value
     * \param[in] resetValue  Will reset a value within bounds if true and if
     *                        the value is outside the bounds
     */
    ForceFieldParameter(const std::string &unit,
                        double             value,
                        double             uncertainty,
                        unsigned int       ntrain,
                        double             minimum,
                        double             maximum,
                        Mutability         mutability,
                        bool               strict,
                        bool               nonNegative,
                        bool               resetValue=false)
                         : 
    unit_(unit), value_(value), originalValue_ (value),
    uncertainty_(uncertainty), originalUncertainty_(uncertainty),
    ntrain_(ntrain), originalNtrain_(ntrain),
    minimum_(minimum), maximum_(maximum), mutability_(mutability),
    strict_(strict), nonNegative_(nonNegative)
    {
        if (resetValue)
        {
            // The value is set to the nearest bound, which is recorded such
            // that the force field can tell, see wasReset and
            // validateForceField, rather than written for each parameter
            if (mutability_ == Mutability::Bounded &&
                (value_ < minimum_ || value_ > maximum_))
            {
                reset_       = true;
                valueBefore_ = value_;
                value_       = (value_ < minimum_) ? minimum_ : maximum_;
            }
        }
        else if (value_ < minimum_ || value_ > maximum_)
        {
            GMX_THROW(gmx::InvalidInputError( 
                      gmx::formatString("Value for force field parameter %g (%s) not within bounds [%g, %g]",
                      value_, unit_.c_str(), minimum_, maximum_).c_str()));
        }
        if (minimum_ == maximum_ && mutability_ == Mutability::Bounded)
        {
            mutability_ = Mutability::Fixed;
        }
        calculateInternalValue();
    }
        
    /*! \brief Whether the value was outside the bounds and was set to the nearest bound
     * This happens when a force field is read, see the constructor.
     */
    bool wasReset() const { return reset_; }

    //! \return the value before it was set to the nearest bound, see wasReset
    double valueBeforeReset() const { return valueBefore_; }

    /*! \brief Give a parameter that depends on others a new value
     * The same as assigning ForceFieldParameter(unit(), value, 0, 1, value,
     * value, Mutability::Dependent, true, true), without copying the unit.
     * \param[in] value         The value
     * \param[in] internalValue The value in the units of the program
     */
    void setDependent(double value, double internalValue)
    {
        value_               = value;
        internalValue_       = internalValue;
        originalValue_       = value;
        reset_               = false;
        valueBefore_         = 0;
        uncertainty_         = 0;
        originalUncertainty_ = 0;
        ntrain_              = 1;
        originalNtrain_      = 1;
        minimum_             = value;
        maximum_             = value;
        mutability_          = Mutability::Dependent;
        strict_              = true;
        nonNegative_         = true;
        index_               = 0;
        updated_             = false;
    }

    //! \brief Return unit of parameter
    const std::string &unit() const { return unit_; }
    
    //! \brief Set unit of parameter
    void setUnit(const std::string &unit) { unit_ = unit; }
    
    //! \brief Return index (an externally determined identifier)
    size_t index() const { return index_; }

    /*! \brief Set the index
     * \param[in] index The index
     */
    void setIndex(size_t index) { index_ = index; }
    
    //! \brief Return current parameter value
    double value() const { return value_; }
    
    //! \brief Return current parameter value
    double internalValue() const { return internalValue_; }
    
    //! \brief Return original parameter value
    double originalValue() const { return originalValue_; }
    
    /*! \brief Set the parameter value if not fixed
     *
     * \param[in] value The new value
     * Will throw an exception if the strict flag is true and the value
     * is out of range, or the variable is alltogher fixed,
     */
    void setValue(double value);
    /*! \brief Set the value ignoring mutability and bounds
     * \param[in] value The new value
     */    
    void forceSetValue(double value);
    //! \brief Return the current uncertainty in this value
    double uncertainty() const { return uncertainty_; }
    
    //! \brief Return the original uncertainty in this value
    double originalUncertainty() const { return originalUncertainty_; }
    
    //! \brief Return the number of training points used
    unsigned int ntrain() const { return ntrain_; }

    /*!\brief Set the number of training point used for this value if not fixed
     * \param[in] ntrain  The number of data points
     * \throws an exception if the strict flag is true and the
     * variable is not mutable.
     */ 
    void setNtrain(unsigned int ntrain);
    
    /*! \brief Add one to the number of training points
     * \throws when the parameter is not mutable and strict flag is true
     */
    void incrementNtrain() { setNtrain(ntrain_+1); }

    /*! \brief Set the uncertainty in this value if not fixed
     * \param[in] uncertainty  The new uncertainty
     * \throws an exception if the strict flag is true and the
     * variable is not mutable.
     */
    void setUncertainty(double uncertainty);

    //! \brief Return minimum allowed value
    double minimum() const { return minimum_; }
    
    /*! \brief Set minimum allowed value irrespective of mutability
     * \param[in] minimum the new value
     * \return true if successful, false otherwise
     */
    bool setMinimum(double minimum)
    { 
        if (!nonNegative_ || minimum >= 0)
        {
            minimum_ = minimum;
            return true;
        }
        return false;
    }
    
    /*! \brief Force non-negative values
     * The value is also set to the nearest bound if it is outside the
     * bounds, which is recorded, see wasReset.
     */
    void setNonNegative()
    {
        double before = value_;
        minimum_ = std::max(0.0, minimum_);
        maximum_ = std::max(minimum_, maximum_);
        value_   = std::max(minimum_, value_);
        value_   = std::min(maximum_, value_);
        if (value_ != before && !reset_)
        {
            reset_       = true;
            valueBefore_ = before;
        }
        nonNegative_ = true;
        calculateInternalValue();
    }

    //! Allow negative values again, the bounds and the value are not changed
    void clearNonNegative() { nonNegative_ = false; }

    //! \return true if parameter should not be negative
    bool nonNegative() const { return nonNegative_; }
    
    //! \brief Return maximum allowed value
    double maximum() const { return maximum_; }
    
    //! \brief Set maximum allowed value irrespective of mutability
    void setMaximum(double maximum) { maximum_ = maximum; }

    //! \return whether the parameter was recently updated
    bool updated() const { return updated_; }

    /*! \brief Change the updated state
     * \param[in] updated The new value
     */
    void setUpdated(bool updated) { updated_ = updated; }

    //! \brief Return how this parameter may be changed
    Mutability mutability() const { return mutability_; }

    //! \brief Change the mutability of this parameter
    void setMutability(Mutability m) { mutability_ = m; }

    //! \brief Return whether this parameter is mutable at all
    bool isMutable() const
    { 
        return (mutability_ != Mutability::Fixed &&
                mutability_ != Mutability::ACM &&
                mutability_ != Mutability::Dependent); 
    }

    //! \brief Return whether or not to throw on value errors
    bool strict() const { return strict_; }

    /*! \brief Copy all values except the index
     * \param[in] src Source to take the values from
     */
    void copy(const ForceFieldParameter &src);

    /*! \brief Send the contents to another processor
     * \param[in] cr   Communication data structure
     * \param[in] dest Processor id to send the data to
     */
    CommunicationStatus Send(const CommunicationRecord *cr, int dest)  const;

    /*! \brief Bcast contents to and from other processors
     * \param[in] cr   Communication data structure
     * \param[in] root The MPI root
     * \param[in] comm Communication structure
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);

    /*! \brief Receive contents from another processor
     * \param[in] cr  Communication data structure
     * \param[in] src Processor id to receive the data from
     */
    CommunicationStatus Receive(const CommunicationRecord *cr, int src);

 private:
    //! Update the internal value by converting units
    void calculateInternalValue() { internalValue_ = convertToGromacs(value_, unit_); }
    //! The unit of the parameter
    std::string  unit_;
    //! The current value of the parameter
    double       value_               = 0;
    //! The current internal value of the parameter
    double       internalValue_       = 0;
    //! The original value of the parameter
    double       originalValue_       = 0;
    //! Whether the value was set to the nearest bound when it was read, see wasReset
    bool         reset_               = false;
    //! The value before that, see wasReset
    double       valueBefore_         = 0;
    //! The current uncertainty in the parameter
    double       uncertainty_         = 0;
    //! The original value of the uncertainty
    double       originalUncertainty_ = 0;
    //! The number of training points used
    unsigned int ntrain_              = 0;
    //! The original number of training points used
    unsigned int originalNtrain_      = 0;
    //! Minimum allowed value for the parameter
    double       minimum_             = 0;
    //! Maximum allowed value for the parameter
    double       maximum_             = 0;
    //! In what way this parameter is mutable
    Mutability   mutability_          = Mutability::Fixed;
    /*! Whether or not to throw an exception in case value or
     * uncertainty is set incorrectly
     */
    bool         strict_              = false;
    //! Whether this parameter should be strictly positive
    bool         nonNegative_         = false;
    //! Externally determined index
    size_t       index_               = 0;
    //! Whether or not the parameter has been recently updated
    bool         updated_             = false;
};

//! Mapping from string to force field parameters
typedef std::map<std::string, ForceFieldParameter> ForceFieldParameterMap;

} // namespace alexandria

#endif
