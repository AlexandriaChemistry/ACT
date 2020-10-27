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
#include <map>
#include <string>

#include "gromacs/math/units.h"
#include "gromacs/mdtypes/commrec.h"

#include "communication.h"
#include "gmx_simple_comm.h"
#include "mutability.h"

namespace alexandria
{

/*! \brief Class to hold one force field parameter
 */
class ForceFieldParameter
{
 public:
    //! \brief default constructor
    ForceFieldParameter() {}
    
    /*! \brief Constructor initiating all parameters.
     *
     * TODO: Check unit
     * \param[in] unit        Physical unit of parameter, e.g. nm or fs
     * \param[in] value       Actual value of the parameter
     * \param[in] uncertainty Uncertainty in the value
     * \param[in] ntrain      Number of data points used for training
     * \param[in] minimum     Minimum allowed value
     * \param[in] maximum     Maximum allowed value
     * \param[in] mutability  In what way this parameter may be changed
     * \param[in] strict      Throw an exception in case of value errors
     */
    ForceFieldParameter(const std::string &unit,
                        double             value,
                        double             uncertainty,
                        uint64_t           ntrain,
                        double             minimum,
                        double             maximum,
                        Mutability         mutability,
                        bool               strict)
                         : 
    unit_(unit), value_(value), originalValue_ (value),
    uncertainty_(uncertainty), originalUncertainty_(uncertainty),
    ntrain_(ntrain), originalNtrain_(ntrain),
    minimum_(minimum), maximum_(maximum), mutability_(mutability),
    strict_(strict) {}
        
    //! \brief Return unit of parameter
    const std::string &unit() const { return unit_; }
    
    //! \brief Return GROMACS unit of parameter
    int unitGromacs() const { return string2unit(unit_.c_str()); }

    //! \brief Return index (an externally determined identifier)
    size_t index() const { return index_; }

    /*! \brief Set the index
     * \param[in] index The index
     */
    void setIndex(size_t index) { index_ = index; }
    
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
    
    //! \brief Return the number of training points used
    uint64_t ntrain() const { return ntrain_; }
    
    /*!\brief Set the number of training point used for this value if not fixed
     * \param[in] ntrain  The number of data points
     * \throws an exception if the strict flag is true and the
     * variable is not mutable.
     */ 
    void setNtrain(uint64_t ntrain);
    
    /*! \brief Add one to the number of training points
     * \throws when the parameter is not mutable and strict flag is true
     */
    void incrementNtrain() { setNtrain(ntrain_+1); }

    /*! \brief Subtract one from the number of training points
     * \throws when the parameter is not mutable and strict flag is true
     */
    void decrementNtrain() { setNtrain(ntrain_-1); }

    /*! \brief Set the uncertainty in this value if not fixed
     * \param[in] uncertainty  The new uncertainty
     * \throws an exception if the strict flag is true and the
     * variable is not mutable.
     */
    void setUncertainty(double uncertainty);
    
    //! \brief Return minimum allowed value
    double minimum() const { return minimum_; }
    
    //! \brief Return maximum allowed value
    double maximum() const { return maximum_; }
    
    //! \brief Return how this parameter may be changed
    Mutability mutability() const { return mutability_; }
    
    //! \brief Return whether this parameter is mutable at all
    bool isMutable() const { return mutability_ == Mutability::Free || mutability_ == Mutability::Bounded; }
    
    //! \brief Return whether or not to throw on value errors
    bool strict() const { return strict_; }

    /*! \brief Send the contents to another processor
     * \param[in] cr   Communication data structure
     * \param[in] dest Processor id to send the data to
     */
    CommunicationStatus Send(const t_commrec *cr, int dest)  const;

    /*! \brief Receive contents from another processor
     * \param[in] cr  Communication data structure
     * \param[in] src Processor id to receive the data from
     */
    CommunicationStatus Receive(const t_commrec *cr, int src);

 private:
    //! The unit of the parameter
    std::string unit_;
    //! The current value of the parameter
    double      value_               = 0;
    //! The original value of the parameter
    double      originalValue_       = 0;
    //! The current uncertainty in the parameter
    double      uncertainty_         = 0;
    //! The original value of the uncertainty
    double      originalUncertainty_ = 0;
    //! The number of training points used
    uint64_t    ntrain_              = 0;
    //! The original number of training points used
    uint64_t    originalNtrain_      = 0;
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
    //! Externally determined index
    size_t      index_               = 0;
};

typedef std::map<std::string, ForceFieldParameter> ForceFieldParameterMap;

} // namespace alexandria

#endif
