/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#ifndef MOLPROP_H
#define MOLPROP_H

#include <string.h>

#include <map>
#include <string>
#include <vector>

#include "act/alexandria/topology.h"
#include "act/molprop/composition.h"
#include "act/molprop/experiment.h"
#include "act/molprop/fragment.h"
#include "act/molprop/molpropobservable.h"
#include "act/molprop/phase.h"
#include "act/forcefield/forcefield.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/real.h"

/*! \brief
 * Contains all classes related to alexandria force field tools
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
namespace alexandria
{

class CommunicationRecord;

/*! \brief
 * Contains molecular properties from a range of sources.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolProp
{
private:
    int                               index_        = -1;
    //! Total number of atoms in this compound
    int                               natom_        = 0;
    //! Whether or not all atoms have proper atom types
    bool                              hasAllAtomTypes_ = false;
    std::string                       molname_, iupac_, cas_, cid_, inchi_;
    std::vector<std::string>          category_;
    //! Elemental composition of the compound
    std::map<const char *,int>        composition_;
    std::vector<Fragment>             fragment_;
    std::vector<Experiment>           exper_;
    std::vector<Bond>                 bond_;

public:
    //! Constructor for a MolProp object
    MolProp() {}
    
    /*! \brief
     * Check the internal consistency of this object
     *
     * \todo Implement this
     */
    void checkConsistency() {};
    
    //! Set the index number for sorting
    void setIndex(int index) { index_ = index; }
    
    //! Return the index number for sorting
    int getIndex() const { return index_; }
    
    /*! \brief
     * Merge the content of another MolProp into this one
     *
     * \param[in] mpi The object to be merged into the present one
     * \return Number of warnings
     * \todo Check and double check
     */
    int Merge(const MolProp *mpi);
    
    //! Dump the contents of this object to a file
    void Dump(FILE *fp) const;
    
    //! \return the LaTeX formula
    std::string texFormula() const;
    
    //! \brief Generate the elemental composition
    void generateComposition();
    
    //! \return The elemental composition
    const std::map<const char *,int> &composition() const { return composition_; }
    
    //! \return Formula for compound or complex
    std::string formula() const;

    //! Set the molname
    void SetMolname(const std::string &molname) { molname_ = molname; }
    
    //! Return the molname
    const std::string &getMolname() const { return molname_; }
    
    //! Set the IUPAC name
    void SetIupac(const std::string &iupac) { iupac_ = iupac; }
    
    //! Return IUPAC name or, if not found, the molname
    const std::string &getIupac() const
    {
        if (iupac_.size() > 0)
        {
            return iupac_;
        }
        else
        {
            return molname_;
        }
    }
    
    //! Set the CAS (Chemical Abstract Service) identifier, see http://www.cas.org/
    void SetCas(const std::string &cas) { cas_.assign(cas); }
    
    //! Return the CAS (Chemical Abstract Service) identifier, see http:://www.cas.org
    const std::string &getCas() const { return cas_; }
    
    //! Set the CID (Chemspider identifier) see http:://www.chemspider.com
    void SetCid(const std::string &cid) { cid_ = cid; }
    
    //! Return the CID (Chemspider identifier) see http:://www.chemspider.com
    const std::string &getCid() const { return cid_; }
    
    //! Set the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
    void SetInchi(const std::string &inchi) { inchi_ = inchi; }
    
    //! Return the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
    const std::string &getInchi() const { return inchi_; }
    
    /*! \brief Find a QM property according to JobType below
     * \param[in] mpo The property of interest
     * \param[in] T   The temperature at which the property
     *                should have been measured. If -1 any T
     *                may be returned.
     * \param[in] jt  The JobType, OPT or SP. In the latter case
     *                the first in the list will be returned.
     * \return the property of interest or nullptr if not found
     */
    const GenericProperty *qmProperty(MolPropObservable  mpo, 
                                      double             T,
                                      JobType            jt) const;

    /*! \brief Find an Experimental property
     * \param[in] mpo The property of interest
     * \param[in] T   The temperature at which the property
     *                should have been measured. If -1 any T
     *                may be returned.
     * \return the property of interest or nullptr if not found
     */
    const GenericProperty *expProperty(MolPropObservable  mpo, 
                                       double             T) const;

    //! Add a classification category for this molecule
    void AddCategory(const std::string &category)
    {
        if (!SearchCategory(category)) { category_.push_back(category); }
    }
    
    //! Return the number of categories
    int NCategory() { return category_.size(); }
    
    //! Return category array as const
    const std::vector<std::string> &categoryConst() const { return category_; }

    //! Return true if catname is an existing category
    bool SearchCategory(const std::string &catname) const;
    
    //! Clear the category vector
    void clearCategory() { category_.clear(); }
    
    //! Return number of atoms in the compound
    int NAtom() const { return natom_; }
    
    //! Returns boolean stating whether all atoms have valid atom types
    bool hasAllAtomTypes() const { return hasAllAtomTypes_; }
    
    //! Add a Bond element
    void AddBond(const Bond &b) { bond_.push_back(b); }
    
    //! Check whether a Bond element is present already
    bool BondExists(const Bond &b);
    
    //! Return the number of Bond elements
    int NBond() const { return bond_.size(); }
    
    //! Return const vector of bonds
    const std::vector<Bond> &bondsConst() const { return bond_; }
    
    //! Return pointer to the whole bond array for editing
    std::vector<Bond> *bonds() { return &bond_; }

    //! Return the fragments of this molecule
    const std::vector<Fragment> &fragments() const { return fragment_; }
    
    //! Return the fragments of this molecule
    const std::vector<Fragment> *fragmentPtrConst() const { return &fragment_; }
    
    //! Return editable fragments of this molecule
    std::vector<Fragment> *fragmentPtr() { return &fragment_; }
    
    /*! Add a new fragment
     * \param[in] f The new fragment
     */    
    void addFragment(const Fragment &f) { fragment_.push_back(f); }

    /*! \brief Renumber the residues to one residue per fragment.
     */
    void renumberResidues();

    /*! Generate fragments based on bonds
     * \param[in] pd     The force field needed for looking up atom props
     * \param[in] qtotal The total charge of the molecule
     */    
    void generateFragments(const ForceField *pd,
                           double            qtotal);
    
    //! Clear the fragment information
    void clearFragments() { fragment_.clear(); }
    
    //! Return the total charge by summing over all fragments
    int totalCharge();
    
    //! Return the total multiplicity from all fragments
    int totalMultiplicity();
    
    //! Return the total mass by summing over all fragments
    double totalMass() const;
    
    //! Return the total symmetry number
    int symmetryNumber() const;

    /*! \brief Return the first calculation that matches the specification
     * \param[in] job The JobType, typically JobType::Opt
     * \return the Experiment or nullptr
     */
    const Experiment *findExperimentConst(JobType job) const;
    
    /*! \brief Return the first calculation that matches the specification
     * \param[in] job The JobType, typically JobType::Opt
     * \return the Experiment or nullptr
     */
    Experiment *findExperiment(JobType job);
    
    //! Add an experiment
    void AddExperiment(const Experiment &myexp) { exper_.push_back(myexp); }
    
    //! Return const vector of experiments
    const std::vector<Experiment> &experimentConst() const { return exper_; }

    //! Return mutable vector of experiments
    std::vector<Experiment> *experiment() { return &exper_; }
    
    //! Return end of experiment vector
    ExperimentConstIterator EndExperiment()  const { return exper_.end(); }
    
    //! Return pointer to the last inserted experiment or nullptr if the number of experiments is zero
    Experiment *LastExperiment()
    {
        if (!exper_.empty())
        {
            return &(exper_.back());
        }
        else
        {
            return nullptr;
        }
    }
    
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
     * Receives this object over an MPI connection
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

//! Iterator over MolProp items
using  MolPropIterator      = typename std::vector<MolProp>::iterator;

//! Const iterator over MolProp items
using  MolPropConstIterator = typename std::vector<MolProp>::const_iterator;

/*! \brief Utility to compare temperatures
 *
 * Compares two temperatures
 * \param[in] Tref The reference temperature
 * \param[in] T    The temperature to be tested
 * \return true if the reference T < 0, or the difference between the two
 *              T is negligable.
 */
bool bCheckTemperature(double Tref, double T);

}

#endif
