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

#ifndef ACT_MOLPROP_EXPERIMENT_H
#define ACT_MOLPROP_EXPERIMENT_H

#include <map>
#include <string>
#include <vector>

#include "act/molprop/composition.h"
#include "act/molprop/molpropobservable.h"

namespace alexandria
{

/*! \brief
 * Enum describing the type of the QM job computed by the Gaussian software
 */
enum class JobType { OPT, ESP, SP, TOPOLOGY, UNKNOWN };

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
        reference_(reference),
        conformation_(conformation)
    {}
    
    //! Constructor initiating a Calculation
    Experiment(const std::string &program,
               const std::string &method,
               const std::string &basisset,
               const std::string &reference,
               const std::string &conformation,
               const std::string &datafile,
               JobType            jtype);
    
    //#define EXPERIMENT_DESTRUCTOR 1
#ifdef EXPERIMENT_DESTRUCTOR
    //! Destructor
    ~Experiment();
    
    Experiment(const Experiment& other);

    Experiment& operator=(const Experiment& copyFrom) = default;
    
    Experiment(Experiment &&other);

    Experiment& operator=(Experiment &&) = default;
#endif

    //! Return the type of data
    DataSource dataSource() const { return dataSource_; }
    
    //! Dump the contents of this object to a file
    void Dump(FILE *fp) const;

    //! Set the experiment id
    void setId(int id) { id_ = id; }
    
    //! Get the id
    int id() const { return id_; }

    /*! \brief Whether a property is present
     * \param[in] mpo Theproperty
     * \return true if found
     */
    bool hasProperty(MolPropObservable mpo) const
    {
        return property_.find(mpo) != property_.end();
    }
    
    /*! Return a property vector for editing
     * \param[in] mpo The property to look for
     * \return The property stuctures
     */
    std::vector<GenericProperty *> *property(MolPropObservable mpo)
    {
        return &property_[mpo];
    }
    /*! Return a const property vector
     * \param[in] mpo The property to look for
     * \return The property stuctures
     */
    const std::vector<GenericProperty *> &propertyConst(MolPropObservable mpo)
    {
        return property_[mpo];
    }
    /*! Add property data
     * \param[in] mpo
     * \param[in] gp
     */  
    void addProperty(MolPropObservable mpo, GenericProperty *gp);
    
    //! \return map of properties for inspection only
    const std::map<MolPropObservable, std::vector<GenericProperty *> > &propertiesConst() const { return property_; }

    //! \return map of properties for editing
    std::map<MolPropObservable, std::vector<GenericProperty *> > *properties() { return &property_; }

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
    
    //! Return iterator pointing to this particular atom or EndAtom() if not found
    CalcAtomIterator searchAtom(CalcAtom ca);
    
    //! Return a complete coordinate array
    const std::vector<gmx::RVec> &getCoordinates() const { return coordinates_; }
    
    //! Copy the CalcAtom data to a separate array that can be extracted as above
    void setCoordinates();
    
    //! Return a complete forces array
    const std::vector<gmx::RVec> &getForces() const { return forces_; }
    
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
     * \param[in] cr   data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr,
                             int                        dest) const;
    
    /*! \brief
     * Broadcast this object over an MPI connection
     *
     * \param[in] cr   data structure for MPI communication
     * \param[in] root The MPI root
     * \param[in] comm MPI communicator
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  data structure for MPI communication
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
    /*! \brief Utility to get reference and lot
     * \param[out] reference Reference to literature, if not nullptr
     * \param[out] lot       Level of theory, if not nullptr
     */
    void getReferenceLot(std::string *reference,
                         std::string *lot) const;
    DataSource                           dataSource_ = dsExperiment;
    std::string                          reference_;
    std::string                          conformation_;
    std::string                          program_;
    std::string                          method_;
    std::string                          basisset_;
    std::string                          datafile_;
    JobType                              jobtype_ = JobType::UNKNOWN;
    std::vector<CalcAtom>                catom_;
    
    //! Experiment ID
    int                                  id_ = 0;
    std::map<MolPropObservable, std::vector<GenericProperty *> > property_;
    
    // std::vector<ElectrostaticPotential>  potential_;
    std::vector<gmx::RVec>               coordinates_;
    std::vector<gmx::RVec>               forces_;
};

//! Iterator over Experiment items
using  ExperimentIterator      = typename std::vector<Experiment>::iterator;

//! Const iterateor over Experiment items
using  ExperimentConstIterator = typename std::vector<Experiment>::const_iterator;

} // namespace alexandria

#endif
