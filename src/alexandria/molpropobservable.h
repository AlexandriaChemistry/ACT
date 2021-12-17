/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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

#ifndef MOLPROPOBSERVABLE_H
#define MOLPROPOBSERVABLE_H

#include <map>
#include <string>
#include <vector>

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"

#include "communication.h"
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
    //! Polarizability
    POLARIZABILITY,
    //! Hartree fock energy
    HF,
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
    //! Zero point energy
    ZPE,
    //! Molecular energy
    EMOL,
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
const char *mpo_unit(MolPropObservable MPO);

/*! \brief Deduct MolPropObservable from string
 * \param[in] str The string to use
 * \return the corresponding MolPropObservable
 * \throws if there is no corresponding MolPropObservable
 */
MolPropObservable stringToMolPropObservable(const std::string &str);

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
     * The type of property
     */
    MolPropObservable mpo_;
    /*! \brief
     * Temperature at which the property is measured or computed.
     */
    double      T_;
    /*! \brief
     * Phase in which the property is measured or computed: e.x. Gas, Liquid, and Solid
     */
    ePhase      eP_;
public:
    //! Default constructor
    GenericProperty() { T_ = 0; eP_ = epNR; };
    
    /*! \brief
     * Creates a new GenericProperty object.
     *
     * \param[in] mpo Type of the property
     * \param[in] T   Temperature
     * \param[in] ep  The phase
     */
    GenericProperty(MolPropObservable mpo,
                    double             T,
                    ePhase             ep) : 
            mpo_(mpo), T_(T), eP_(ep)
    {}

    /*! \brief
     * Return the property type
     */
    const char *getType() const { return mpo_name(mpo_); }
    
    /*! \brief
     * Return the unit of the property
     */
    const char *getUnit() const { return mpo_unit(mpo_); }
    
    /*! \brief
     * Return the temperature
     */
    double getTemperature() const { return T_; }
    
    /*! \brief
     * Return the phase
     */
    ePhase getPhase() const { return eP_; }
    
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
    
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   GROMACS data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(t_commrec *cr, int dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  GROMACS data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(t_commrec *cr,
                                int        src);

    virtual double getValue() const = 0;
    
    virtual double getError() const = 0;
    
    virtual const std::vector<double> &getVector() const = 0;
    
    virtual const tensor &getTensor() const = 0;
};

typedef std::vector<GenericProperty>::const_iterator GPConstIterator;

/*! \brief
 * Contains the elements of the molecular quadrupole
 *
 * The six elements of the upper diagonal of a quadrupole tensor are stored.
 * The values are dependent on the orientation of the molecule and are relative
 * to the center of charge in case that the total charge or total dipole of
 * the molecules are non-zero.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularQuadrupole : public GenericProperty
{
private:
    //! The quadrupole moment
    tensor quad_ = {{ 0 }};
public:
    //! Default constructor
    MolecularQuadrupole() {}
    
    //! Constructor initiating all elements of the quadrupole tensor
    MolecularQuadrupole(double T,
                        double xx, double yy, double zz,
                        double xy, double xz, double yz) :
        GenericProperty(MolPropObservable::QUADRUPOLE, T, epGAS)
    { Set(xx, yy, zz, xy, xz, yz); };
    
    //! Set all the elements of the qudrupole tensor
    void Set(double xx, double yy, double zz, double xy, double xz, double yz)
    { 
        quad_[XX][XX] = xx;
        quad_[YY][YY] = yy;
        quad_[ZZ][ZZ] = zz;
        quad_[XX][YY] = xy;
        quad_[XX][ZZ] = xz; 
        quad_[YY][ZZ] = yz;
    };
    
    double getValue() const { crash(); }
    
    double getError() const { crash(); }
    
    const std::vector<double> &getVector() const { crash(); }
    
    const tensor &getTensor() const { return quad_; }
    
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   GROMACS data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(t_commrec *cr,
                             int        dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  GROMACS data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(t_commrec *cr,
                                int        src);
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
class MolecularPolarizability : public MolecularQuadrupole
{
private:
    double average_, error_;
public:
    //! Default constructor
    MolecularPolarizability() {}
    
    //! Constructor initiating all elements of the quadrupole tensor
    MolecularPolarizability(double T,
                            double xx, double yy, double zz,
                            double xy, double xz, double yz,
                            double average, double error)
        : MolecularQuadrupole(T, xx, yy, zz, xy, xz, yz)
    { 
        average_ = average;
        error_   = error;
    }
    
    double getValue() const;
    
    double getError() const { return error_; }
    
    const std::vector<double> &getVector() const { crash(); }

    const tensor &getTensor() const
    {
        return MolecularQuadrupole::getTensor();
    }
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   GROMACS data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(t_commrec *cr,
                             int        dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  GROMACS data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(t_commrec *cr,
                                int        src);
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

    //! Constructor storing all properties related to this energy term
    MolecularEnergy(MolPropObservable mpo,
                    double T,
                    ePhase ep,
                    double average,
                    double error)
        : GenericProperty(mpo, T, ep)
    { 
        Set(average, error);
    };
    
    //! Set the average and error for the energy
    void Set(double average, double error) 
    { 
        average_ = average; 
        error_   = error; 
    };

    double getValue() const { return average_; }
    
    double getError() const { return error_; }

    const std::vector<double> &getVector() const { crash(); }

    const tensor &getTensor() const { crash(); }
    
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   GROMACS data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(t_commrec *cr,
                             int        dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  GROMACS data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(t_commrec *cr,
                                int        src);
};

/*! \brief
 * Contains the dipole vector
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularDipole : public GenericProperty
{
private:
    std::vector<double> mu_;
    double              average_;
    double              error_;
public:
    //! Default constructor
    MolecularDipole() {}
    
    //! Constructor storing all properties related to this dipole
    MolecularDipole(double T,
                    double x, double y, double z, double aver, double error)
        : GenericProperty(MolPropObservable::DIPOLE, T, epGAS)
    { Set(x, y, z, aver, error); }

    //! Set all properties related to this dipole
    void Set(double x, double y, double z, double aver, double error)
    {
        mu_.clear();
        mu_.push_back(x);
        mu_.push_back(y);
        mu_.push_back(z);
        average_ = aver;
        error_   = error;
    };

    //! Return all properties of this dipole
    const std::vector<double> &getVector() const { return mu_; }

    //! Return the average dipole value
    double getValue() const { return average_; }

    //! Return the error in the average dipole
    double getError() const { return error_; }
    
    const tensor &getTensor() const { crash(); }

    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   GROMACS data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(t_commrec *cr,
                             int        dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  GROMACS data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(t_commrec *cr, int src);
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
class ElectrostaticPotential
{
private:
    std::string xyzUnit_, vUnit_;
    int         espID_;
    double      x_, y_, z_, V_;
public:
    //! Default constructor
    ElectrostaticPotential() {}
    
    //! Constructor that set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
    ElectrostaticPotential(const std::string &xyz_unit,
                           const std::string &V_unit,
                           int espid, double x, double y, double z, double V)
    { Set(xyz_unit, V_unit, espid, x, y, z, V); };
    
    //! Constructor that set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
    ElectrostaticPotential(const char *xyz_unit, const char *V_unit,
                           int espid, double x, double y, double z, double V)
    { Set(xyz_unit, V_unit, espid, x, y, z, V); };
    
    //! Set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
    void Set(const std::string &xyz_unit,
             const std::string &V_unit,
             int espid, double x, double y, double z, double V)
    { xyzUnit_ = xyz_unit; vUnit_ = V_unit; espID_ = espid; x_ = x; y_ = y; z_ = z; V_ = V; };

    //! Set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
    void Set(const char *xyz_unit, const char *V_unit, int espid, double x, double y, double z, double V)
    {
        std::string x_(xyz_unit), V_(V_unit); Set(x_, V_, espid, x, y, z, V);
    }
    
    //! Return the units of coordinates and potential, the ESP id, the coordinates and the potential itself
    void get(char **xyz_unit, char **V_unit, int *espid, double *x, double *y, double *z, double *V) const
    {
        *xyz_unit = strdup(xyzUnit_.c_str()); *V_unit = strdup(vUnit_.c_str()); *espid = espID_; *x = x_; *y = y_; *z = z_; *V = V_;
    };
    
    //! Return the unit of the coordinates
    const std::string &getXYZunit() const { return xyzUnit_; }
    
    //! Return the unit of the potential
    const std::string &getVunit() const { return vUnit_; }
    
    //! Return the ESP id from the original calculation
    int getEspid() const { return espID_; }
    
    //! Return the X coordinate of the ESP point
    double getX() const { return x_; }
    
    //! Return the Y coordinate of the ESP point
    double getY() const { return y_; }
    
    //! Return the Z coordinate of the ESP point
    double getZ() const { return z_; }
    
    //! Return the electrostatic potential at this point in space
    double getV() const { return V_; }

    //! Set the potential for this instance of the class
    void SetV(double V) { V_ = V; }
    
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   GROMACS data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(t_commrec *cr, int dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  GROMACS data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(t_commrec *cr, int src);
};

} // namespace alexandria

#endif
