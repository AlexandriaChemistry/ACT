/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#ifndef MOLPROPOBSERVABLE_H
#define MOLPROPOBSERVABLE_H

#include <array>
#include <map>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"

#include "act/utility/communicationrecord.h"
#include "phase.h"

namespace alexandria
{

/*! \brief
 * Enumerated type holding the types of observables stored in MolProp
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum class MolPropObservable {
    //! Electrostatic potential
    POTENTIAL,
    //! Dipole
    DIPOLE,
    //! Quadrupole
    QUADRUPOLE,
    //! Octupole
    OCTUPOLE,
    //! Hexadecapole
    HEXADECAPOLE,
    //! Polarizability
    POLARIZABILITY,
    //! Harmonic frequencies 
    FREQUENCY,
    //! IR intensity
    INTENSITY,
    //! Hartree fock energy
    HF,
    //! QM energy of molecule minus atomic contribution
    DELTAE0,
    //! Interaction energy between fragments
    INTERACTIONENERGY,
    //! Delta H formation
    DHFORM,
    //! Delta G formation
    DGFORM,
    //! Delta S formation
    DSFORM,
    //! Standard entropy
    ENTROPY,
    //! Translational entropy
    STRANS,
    //! Rotational entropy
    SROT,
    //! Vibrational entropy
    SVIB,
    //! Heat capacity at constant pressure
    CP,
    //! Heat capacity at constant volume
    CV,
    //! Zero point energy
    ZPE,
    //! Charge
    CHARGE,
    //! Coordinates
    COORDINATES
};

/*! \brief
 * Enum to select either QM or Experimental data or either
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum iqmType {
    //! Experimental data only
    Exp,
    //! Both experimental and quantum-chemical
    Both,
    //! Quantum-chemistry only
    QM
};

/*! \brief
 * Strings describing the MolPropObservable enum elements
 * \param[in] MPO The observable
 * \return the string describing the unit
 */
const char *mpo_name(MolPropObservable MPO);

/*! \brief
 * Strings describing the MolPropObservable enum units
 * \param[in] MPO The observable
 * \return the string describing the unit
 */
const char *mpo_unit2(MolPropObservable MPO);

/*! \brief Deduct MolPropObservable from string
 * \param[in]   str The string to use
 * \param[[out] mpo The corresponding MolPropObservable
 * \return false if there is no corresponding MolPropObservable, true if there is
 */
bool stringToMolPropObservable(const std::string &str, MolPropObservable *mpo);

#define crash() GMX_THROW(gmx::InternalError("This should not be called"))

/*! \brief
 * Generic molecular property base clase
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class GenericProperty
{
private:
    /*! \brief
     * The type of property. Upon entering data, the values should be converted
     * to internal units. When the data should be exported it can be converted
     * back. Conversion should use the functions in act/utility/units.h.
     */
    MolPropObservable mpo_ = MolPropObservable::POTENTIAL;
    //! Subtype of this property
    std::string       type_;
    //! Internal unit of this property
    std::string       unit_;
    //! Original unit of this property
    std::string       inputUnit_;
    /*! \brief
     * Temperature at which the property is measured or computed.
     */
    double            T_ = 0;
    /*! \brief
     * Phase in which the property is measured or computed: e.x. Gas, Liquid, and Solid
     */
    ePhase            eP_ = ePhase::PLASMA;
public:
    //! Default constructor
    GenericProperty() {};
    
    /*! \brief
     * Creates a new GenericProperty object.
     *
     * \param[in] mpo       Type of the property
     * \param[in] type      The subtype of this property
     * \param[in] inputUnit The unit for this property
     * \param[in] T         Temperature
     * \param[in] ep        The phase
     */
    GenericProperty(MolPropObservable  mpo,
                    const std::string &type,
                    const std::string &inputUnit,
                    double             T,
                    ePhase             ep) : 
        mpo_(mpo), type_(type), inputUnit_(inputUnit), T_(T), eP_(ep)
    {
    }

    //! Destructor
    virtual ~GenericProperty() {};

    GenericProperty(const GenericProperty& copyFrom) = default;
    GenericProperty& operator=(const GenericProperty& copyFrom) = default;
    GenericProperty(GenericProperty &&) = default;
    GenericProperty& operator=(GenericProperty &&) = default;

    /*! \brief
     * Return the property type
     */
    const char *getType() const { return type_.c_str(); }
    
    /*! \brief
     * Return the internal unit of the property
     */
    const char *getUnit() const { return unit_.c_str(); }
    
    /*! \brief
     * Return the input unit of the property
     */
    const char *getInputUnit() const { return inputUnit_.c_str(); }
    
    /*! \brief
     * Return the temperature
     */
    double getTemperature() const { return T_; }
    
    /*! \brief
     * Return the phase
     */
    ePhase getPhase() const { return eP_; }
    
    /*! \brief
     * Set the internal unit of the property
     *
     **\param[in] unit The unit
     */
    void setUnit(const std::string &unit) { unit_ = unit; }
    
    /*! \brief
     * Set the temperature of the property
     *
     **\param[in] T Temperature
     */
    void setTemperature(double T) { T_ = T; }
    
    /*! \brief
     * Set the phase of the property
     *
     **\param[in] ep Phase of the property
     */
    void setPhase(ePhase ep) { eP_ = ep; }
    
    //! Return my type of property
    MolPropObservable mpo() const { return mpo_; }

    /*! \brief Dump (excerpt of) content tot a file
     * \param[in] fp File pointer
     */
    virtual void Dump(FILE *fp) const;
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   Data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    virtual CommunicationStatus Send(const CommunicationRecord *cr,
                                     int                        dest) const;
    
    /*! \brief
     * Broadcast this object over an MPI connection
     *
     * \param[in] cr   Data structure for MPI communication
     * \param[in] root The MPI root
     * \param[in] comm MPI Communicator
     * \return the CommunicationStatus of the operation
     */
    virtual CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                          int                        root,
                                          MPI_Comm                   comm);

    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  Data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    virtual CommunicationStatus Receive(const CommunicationRecord *cr,
                                        int                        src);

    virtual double getValue() const = 0;
    
    virtual double getError() const = 0;
    
    virtual const std::vector<double> &getVector() const = 0;
    
};

typedef std::vector<GenericProperty>::const_iterator GPConstIterator;

/*! \brief
 * Elements of the electrostatic multipole
 */
class MolecularMultipole : public GenericProperty
{
private:
    //! Average value
    double              average_ = 0;
    //! Error
    double              error_   = 0;
    //! The values of the electric multipole
    std::vector<double> values_;
public:
    //! Default constructor
    MolecularMultipole() {}
    
    /*! Constructor initiating all elements of the multipole
     * \param[in] type      The calculation type
     * \param[in] inputUnit The unit of the values
     * \param[in] T         The temperature
     * \param[in] mpo       The multipole type
     * \throws if mpo is not a multipole
     */
    MolecularMultipole(const std::string &type,
                       const std::string &inputUnit,
                       double             T,
                       MolPropObservable  mpo);
                       
    //! Destructor
    ~MolecularMultipole() {}
    
    MolecularMultipole(const MolecularMultipole& copyFrom) = default;
    MolecularMultipole& operator=(const MolecularMultipole& copyFrom) = default;
    MolecularMultipole(MolecularMultipole &&) = default;
    MolecularMultipole& operator=(MolecularMultipole &&) = default;
    
    /*! \brief Check whether a certain id is present
     * \param[in] id    The name of the parameter, e.g. "xy"
     * \return true if present
     */
    bool hasId(const std::string &id);
    
    /*! \brief Set a value based on the name
     * \param[in] id    The name of the parameter, e.g. "xy", "average", "error"
     * \param[in] value The value
     * \throws if unknown id
     */
    void setValue(const std::string &id, double value);

    void Dump(FILE *fp) const;
    
    double getValue() const
    {
        return average_;
    }
    
    double getError() const
    {
        return error_;
    }
    
    //! \return the values
    const std::vector<double> &getVector() const { return values_; }
    
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   Data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr,
                             int                        dest) const;
    
    /*! \brief
     * Broadcast this object over an MPI connection
     *
     * \param[in] cr   Data structure for MPI communication
     * \param[in] root The MPI root
     * \param[in] comm MPI Communicator
     * \return the CommunicationStatus of the operation
     */
    virtual CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                          int                        root,
                                          MPI_Comm                   comm);

    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  Data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(const CommunicationRecord *cr,
                                int                        src);
};

/*! \brief
 * Harmonic frequencies or intensities. Distinction between the
 * two is through the MolPropObservable
 */
class Harmonics : public GenericProperty
{
private:
    //! The frequencies
    std::vector<double> values_;
public:
    //! Default constructor
    Harmonics() {}
    
    /*! Constructor with units
     * \param[in] unit The unit of the frequency/intensity
     * \param[in] T    The temperature
     * \param[in] mpo  The property 
     * \throw if mpo is not correct
     */
    Harmonics(const std::string &unit,
              double             T,
              MolPropObservable  mpo);

    //! Destructor    
    ~Harmonics() {}

    Harmonics(const Harmonics& copyFrom) = default;
    Harmonics& operator=(const Harmonics& copyFrom) = default;
    Harmonics(Harmonics &&) = default;
    Harmonics& operator=(Harmonics &&) = default;

    /*! \brief Add a value and convert it to internal units.
     * \param[in] frequency The frequency
     * \param[in] intensity The intensity value
     */
    void addValue(double value);
    
    //! \return the frequencies
    const std::vector<double> &getVector() const { return values_; }
    
    double getValue() const
    {
        crash();
    }
    
    double getError() const
    {
        crash();
    }
    
    void Dump(FILE *fp) const;
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   Data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr,
                             int                        dest) const;
    
    /*! \brief
     * Broadcast this object over an MPI connection
     *
     * \param[in] cr  Data structure for MPI communication
     * \param[in] root The MPI root
     * \param[in] comm MPI Communicator
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);
                                
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  Data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(const CommunicationRecord *cr,
                                int                        src);
};

/*! \brief
 * Contains the elements of the molecular polarizability tensor
 *
 * The six elements of the upper diagonal of a polarizability tensor are stored
 * along with the average molecular polarizability and the error if known.
 * The values are dependent on the orientation of the molecule.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularPolarizability : public GenericProperty
{
private:
    double average_, error_;
    //! The polarizability tensor
    tensor alpha_ = {{ 0 }};
public:
    //! Default constructor
    MolecularPolarizability() {}
    
    /*! Constructor initiating all elements of the quadrupole tensor
     * \param[in] type      The calculation type
     * \param[in] inputUnit The unit of the values
     * \param[in] T         The temperature
     * \param[in] xx        The corresponding value
     * \param[in] yy        The corresponding value
     * \param[in] zz        The corresponding value
     * \param[in] xy        The corresponding value
     * \param[in] xz        The corresponding value
     * \param[in] yz        The corresponding value
     * \param[in] average   The average value
     * \param[in] error     The (experimental) error
     */
    MolecularPolarizability(const std::string &type,
                            const std::string &inputUnit,
                            double T,
                            double xx, double yy, double zz,
                            double xy, double xz, double yz,
                            double average, double error);
    
    ~MolecularPolarizability() {}
    
    MolecularPolarizability(const MolecularPolarizability& copyFrom) = default;
    MolecularPolarizability& operator=(const MolecularPolarizability& copyFrom) = default;
    MolecularPolarizability(MolecularPolarizability &&) = default;
    MolecularPolarizability& operator=(MolecularPolarizability &&) = default;
    
    //! Set all the elements of the polarizablity tensor and converts them to internal units
    void Set(double xx, double yy, double zz, double xy, double xz, double yz);
    
    double getValue() const;
    
    double getError() const { return error_; }
    
    const std::vector<double> &getVector() const { crash(); }

    const tensor &getTensor() const
    {
        return alpha_;
    }

    void Dump(FILE *fp) const;
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   Data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr,
                             int                        dest) const;
    
    /*! \brief
     * Broadcast this object over an MPI connection
     *
     * \param[in] cr  Data structure for MPI communication
     * \param[in] root The MPI root
     * \param[in] comm MPI Communicator
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);

    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  Data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(const CommunicationRecord *cr,
                                int                        src);
};

/*! \brief
 * Contains a molecular energy
 *
 * Different energy terms associated with a molecule can be stored based
 * on quantum chemistry calculations or experimental data.
 * For example the Enthalpy of formation at different temperatures.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularEnergy : public GenericProperty
{
private: 
    double average_ = 0;
    double error_  = 0;
public:
    //! Default constructor needed for receiving over a network
    MolecularEnergy() {};

    /*! Constructor storing all properties related to this energy term
     * \param[in] mpo       The energy type
     * \param[in] type      The calculation type
     * \param[in] inputUnit The unit of the values
     * \param[in] T         The temperature
     * \param[in] ep        The physical phase of the measurement
     * \param[in] average   The average value
     * \param[in] error     The (experimental) error
     */
    MolecularEnergy(MolPropObservable mpo,
                    const std::string &type,
                    const std::string &inputUnit,
                    double T,
                    ePhase ep,
                    double average,
                    double error);
                    
    ~MolecularEnergy() {}
    
    MolecularEnergy(const MolecularEnergy& copyFrom) = default;
    MolecularEnergy& operator=(const MolecularEnergy& copyFrom) = default;
    MolecularEnergy(MolecularEnergy &&) = default;
    MolecularEnergy& operator=(MolecularEnergy &&) = default;

    //! Set the average and error for the energy
    void Set(double average, double error) 
    { 
        average_ = average; 
        error_   = error; 
    };

    double getValue() const { return average_; }
    
    double getError() const { return error_; }

    const std::vector<double> &getVector() const { crash(); }

    void Dump(FILE *fp) const;
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   Data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr,
                             int                        dest) const;
    
    /*! \brief
     * Broadcast this object over an MPI connection
     *
     * \param[in] cr   Data structure for MPI communication
     * \param[in] root The MPI root
     * \param[in] comm MPI Communicator
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);
                                
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  Data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(const CommunicationRecord *cr,
                                int                        src);
};

/*! \brief
 * Contains the electrostatic potential in a coordinate close to a molecule.
 *
 * The electrostatic potential (ESP) can be computed using quantum chemistry and
 * stored in this class.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class ElectrostaticPotential : public GenericProperty
{
private:
    //! Internal coordinate unit
    std::string xyzUnit_;
    //! Input coordinate unit
    std::string xyzInputUnit_;
    //! Internal potential unit
    std::string vUnit_;
    //! Input potential unit
    std::string vInputUnit_;
    //! Electrostatic potential ids
    std::vector<int>       espID_;
    //! Grid coordinates
    std::vector<gmx::RVec> xyz_;
    //! Potential values
    std::vector<double>    V_;
public:
    //! Default constructor
    ElectrostaticPotential() {}
    
    /*! Constructor that set the units of coordinates and potential
     * \param[in] xyzInputUnit Unit of the coordinate
     * \param[in] vInputUnit   Unit of the potential
     */
    ElectrostaticPotential(const std::string &xyzInputUnit,
                           const std::string &vInputUnit);
    
    ~ElectrostaticPotential() {}
    
    ElectrostaticPotential(const ElectrostaticPotential& copyFrom) = default;
    ElectrostaticPotential& operator=(const ElectrostaticPotential& copyFrom) = default;
    ElectrostaticPotential(ElectrostaticPotential &&) = default;
    ElectrostaticPotential& operator=(ElectrostaticPotential &&) = default;
    
    double getValue() const { crash(); }
    
    double getError() const { crash(); }

    const std::vector<double> &getVector() const { crash(); }

    /*! Fill the contents of the ESP
     * Set the ESP id, the coordinates and the potential itself.
     * \param[in] espid        Unique id
     * \param[in] x            X coordinate
     * \param[in] y            Y coordinate
     * \param[in] z            Z coordinate
     * \param[in] V            Potential
     */
    void addPoint(int                espid, 
                  double             x,
                  double             y,
                  double             z,
                  double             V);

    //! Return the unit of the coordinates
    const std::string &getXYZunit() const { return xyzUnit_; }
    
    //! Return the unit of the potential
    const std::string &getVunit() const { return vUnit_; }
    
    //! Return the ESP ida from the original calculation
    const std::vector<int> &espid() const { return espID_; }
    
    //! Return the coordinate of the ESP points
    const std::vector<gmx::RVec> &xyz() const { return xyz_; }
    
    //! Return the electrostatic potential at this point in space
    const std::vector<double> V() const { return V_; }

    void Dump(FILE *fp) const;
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   Data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr, int dest) const;
    
    /*! \brief
     * Broadcast this object over an MPI connection
     *
     * \param[in] cr  Data structure for MPI communication
     * \param[in] root The MPI root
     * \param[in] comm MPI Communicator
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);

    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  Data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(const CommunicationRecord *cr, int src);
};

} // namespace alexandria

#endif
