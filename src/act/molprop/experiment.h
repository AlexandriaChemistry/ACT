/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <map>
#include <string>
#include <vector>

#include "act/molprop/molpropobservable.h"
#include "act/qgen/qtype.h"

namespace alexandria
{

/*! \brief
 * Enum describing the type of the QM job computed by the Gaussian software
 */
enum class JobType {
    OPT       = 0,
    POP       = 1,
    POLAR     = 2,
    G2        = 3,
    G3        = 4,
    G4        = 5,
    CBSQB3    = 6,
    W1U       = 7,
    W1BD      = 8,
    SP        = 9,
    UNKNOWN   = 10
};

/*! \brief
 * Return string corresponding to the job type
 * \param[in] jType the Job type
 */
const char *jobType2string(JobType jType);

/*! \brief
 * Strings describing the job type
 * \param[in] str the String
 * \return the JobType
 * \throws if invalid string is passed
 */
JobType string2jobType(const std::string &str);

/*! \brief
 * Enum describing the source of the data
 */
enum DataSource {
    dsExperiment, dsTheory
};

/*! \brief
 * Return string corresponding to the data source
 */
const char *dataSourceName(DataSource ds);

/*! \brief
 * Return DataSource corresponding to the data source string
 */
DataSource dataSourceFromName(const std::string &name);

/*! \brief
 * Contains data on an atom based on a calculation.
 * Coordinates and Forces are store in the specified units.
 * This class coordinates, name, atom type and an array of
 * Charge values based on different methods.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class CalcAtom
{
private:
    std::string             name_, obType_, coord_unit_, force_unit_, residueName_;
    double                  x_ = 0, y_ = 0, z_ = 0;
    double                  fx_ = 0, fy_ = 0, fz_ = 0;
    int                     atomID_ = 0, residueNumber_ = 0, chainId_ = 0;
    char                    chain_ = ' ';
    std::map<qType, double> q_;
public:
    //! Default constructor
    CalcAtom() {}
    
    //! Constructor initiating the name, type and atomid
    CalcAtom(const char *name, const char *obtype, int atomid)
    {
        name_.assign(name); obType_.assign(obtype); atomID_ = atomid; residueName_.assign(nullptr);
    };
    
    //! Constructor initiating the name, type and atomid
    CalcAtom(const std::string &name,
             const std::string &obtype,
             int                atomid)
    {
        name_ = name; obType_ = obtype; atomID_ = atomid;
    };
    
    
    //! Function returning true if the two atoms are equal
    bool Equal(CalcAtom ca);
    
    /*! \brief Add an AtomicCharge element to the atom
     * \param[in] type Charge type
     * \param[in] q    The charge
     */
    void AddCharge(qType type, double q)
    {
        q_.insert(std::pair<qType, double>(type, q));
    }
    
    /*! \brief Return whether the charge type is present
     * \param[in] type Charge type
     * \return true if found
     */
    bool hasCharge(qType type) const
    {
        return q_.find(type) != q_.end();
    }
    
    /*! \brief Return charge of a certain type
     * \param[in] type Charge type
     * \return charge if found, or crash otherwise
     */
    double charge(qType type) const
    {
        return q_.find(type)->second;
    }
    
    //! \brief Return the whole charge map
    const std::map<qType, double> &chargesConst() const { return q_; }
    
    //! Return the atom id of the atom
    int getAtomid() const { return atomID_; }
    
    //! Return the name of the atom
    const std::string &getName() const { return name_; }
    
    //! Return the OpenBabel type of the atom
    const std::string &getObtype() const { return obType_; }
    
    //! Set the OpenBabel type of the atom
    void setObtype(const std::string &obtype) { obType_ = obtype; }
    
    //! Return the unit of the coordinates of the atom
    const std::string &coordUnit() const { return coord_unit_; }
    
    //! Return the unit of the coordinates of the atom
    const std::string &forceUnit() const { return force_unit_; }
    
    //! Return the name of residue
    const std::string &ResidueName() const { return residueName_; }
    
    //! Return the residue number
    int ResidueNumber() const { return residueNumber_; }

    //! Set the unit of the coordinates of the atom
    void setCoordUnit(const std::string &unit);
    
    //! Set the unit of the forces on the atom
    void setForceUnit(const std::string &unit);
    
    //! Set the residue name for the atom
    void SetResidue(const std::string &residueName, int residueNumber)
    {
        residueName_   = residueName;
        residueNumber_ = residueNumber;
    }
    
    //! Set the chain information (from pdb files)
    void SetChain(int chainId, char chain)
    {
        chainId_ = chainId;
        chain_   = chain;
    }
    
    //! Return chain ID
    int chainId() const { return chainId_; }
    
    //! Return chain
    char chain() const { return chain_; }
    
    //! Set the coordinates of the atom
    void setCoords(double x, double y, double z) { x_ = x; y_ = y; z_ = z; }
    
    //! Return all the coordinates of the atom
    void coords(double *x, double *y, double *z) const
    { *x = x_; *y = y_; *z = z_; }

    //! Set the forces on the atom
    void setForces(double fx, double fy, double fz) { fx_ = fx; fy_ = fy; fz_ = fz; }
    
    //! Return all the coordinates of the atom
    void forces(double *fx, double *fy, double *fz) const
    { *fx = fx_; *fy = fy_; *fz = fz_; }

    //! Return the X coordinate of the atom
    double getX() const { return x_; }
    
    //! Return the Y coordinate of the atom
    double getY() const { return y_; }
    
    //! Return the Z coordinate of the atom
    double getZ() const { return z_; }
    
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   GROMACS data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr,
                             int                        dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  GROMACS data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(const CommunicationRecord *cr,
                                int                        src);
};

//! Iterator over CalcAtom items
using  CalcAtomIterator      = typename std::vector<CalcAtom>::iterator;

//! Const iterator over CalcAtom items
using  CalcAtomConstIterator = typename std::vector<CalcAtom>::const_iterator;

/*! \brief
 * Contains molecular data based on experiments
 *
 * This is a composite class holding the results from experiment, either
 * the dipole, polarizability, energy or the quadrupole. A reference to the
 * publication (or handbook) containing the data is stored, as well as the
 * conformation of the molecule (if known).
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Experiment
{
public:
    //! Default constructor
    Experiment() { }
    
    //! Constructor initiating an Experiment with reference and conformation
    Experiment(const std::string &reference,
               const std::string &conformation) :
        dataSource_(dsExperiment), reference_(reference),
        conformation_(conformation), jobtype_(JobType::UNKNOWN)
    {}
    
    //! Constructor initiating a Calculation
    Experiment(const std::string &program,
               const std::string &method,
               const std::string &basisset,
               const std::string &reference,
               const std::string &conformation,
               const std::string &datafile,
               JobType            jtype);
    
    //! Return the type of data
    DataSource dataSource() const { return dataSource_; }
    
    //! Dump the contents of this object to a file
    void Dump(FILE *fp) const;

    /*! \brief Whether a property is present
     * \param[in] mpo Theproperty
     * \return true if found
     */
    bool hasProperty(MolPropObservable mpo) const
    {
        return property_.find(mpo) != property_.end();
    }
    /*! Add property data
     * \param[in] mpo
     * \param[in] gp
     */  
    void addProperty(MolPropObservable mpo, GenericProperty *gp);
    
    const std::map<MolPropObservable, std::vector<GenericProperty *> > propertyConst() const { return property_; }

    //! Return the molecular conformation
    const std::string &getConformation() const { return conformation_; }
    
    //! Return the literature reference
    const std::string &getReference() const { return reference_; }
    
    //! Return the type of calculation
    const JobType &getJobtype() const { return jobtype_; }
    
    //! Add a CalcAtom object to the list of atoms
    void AddAtom(CalcAtom ca);
    
    //! Return the number of atoms
    int NAtom() const { return catom_.size(); }
    
    //! Return const vector of calcatom
    const std::vector<CalcAtom> &calcAtomConst() const { return catom_; }
    
    //! Return mutable vector of calcatom
    std::vector<CalcAtom> *calcAtom() { return &catom_; }
    
    //! Return iterator pointint to this particular atom or EndAtom() if not found
    CalcAtomIterator searchAtom(CalcAtom ca);
    
    //! Return a complete coordinate array
    const std::vector<gmx::RVec> &getCoordinates() const { return coordinates_; }
    
    //! Return a complete forces array
    const std::vector<gmx::RVec> &getForces() const { return forces_; }
    
    //! Add ElectrostaticPotential element to the array
    void AddPotential(ElectrostaticPotential ep) { potential_.push_back(ep); }
    
    //! Return the number of potential points
    int NPotential() const { return potential_.size(); };
    
    //! Return const vector of electrostaticpotentials
    const std::vector<ElectrostaticPotential> &electrostaticPotentialConst() const { return potential_; }
    
    //! Return the program used to perform the calculation
    const std::string &getProgram() const { return program_; }
    
    //! Return the basis set used to perform the calculation
    const std::string &getBasisset() const { return basisset_; }
    
    //! Return the method used to perform the calculation
    const std::string &getMethod() const { return method_; }
    
    //! Return the datafile from which the calculation output was extracted
    const std::string &getDatafile() const { return datafile_; }
    
    /*! \brief
     * Function that fetches charges from this QM
     * \param[out] q         The charges
     * \param[in]  qtype     The charge type
     * \param[out] reference The literature reference (if nullptr this is ignored)
     * \param[out] lot       The level of theory (if nullptr this is ignored)
     * \return true if successful
     */
    bool getCharges(std::vector<double> *q,
                    qType                qtype,
                    std::string         *reference,
                    std::string         *lot) const;

    //! Merge in another object. Return number of warnings.
    int Merge(const Experiment *src);
    
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   GROMACS data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr,
                             int                        dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  GROMACS data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(const CommunicationRecord *cr,
                                int                        src);
    
    bool hasMolPropObservable(MolPropObservable mpo) const
    {
        return property_.end() != property_.find(mpo);
    }
    
    const std::vector<GenericProperty *> &propertyConst(MolPropObservable mpo) const
    {
        return property_.find(mpo)->second;
    }
private:
    /*! \brief Utility to set reference and lot
     * \param[out] reference Reference to literature, if not nullptr
     * \param[out] lot       Level of theory, if not nullptr
     */
    void getReferenceLot(std::string *reference,
                         std::string *lot) const;
    DataSource                           dataSource_;
    std::string                          reference_;
    std::string                          conformation_;
    std::string                          program_;
    std::string                          method_;
    std::string                          basisset_;
    std::string                          datafile_;
    JobType                              jobtype_;
    std::vector<CalcAtom>                catom_;
    
    std::map<MolPropObservable, std::vector<GenericProperty *> > property_;
    
    std::vector<ElectrostaticPotential>  potential_;
    std::vector<gmx::RVec>               coordinates_;
    std::vector<gmx::RVec>               forces_;
};

//! Iterator over Experiment items
using  ExperimentIterator      = typename std::vector<Experiment>::iterator;

//! Const iterateor over Experiment items
using  ExperimentConstIterator = typename std::vector<Experiment>::const_iterator;

} // namespace alexandria

#endif
