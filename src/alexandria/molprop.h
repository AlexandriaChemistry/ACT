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

#ifndef MOLPROP_H
#define MOLPROP_H

#include <string.h>

#include <map>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/real.h"

#include "communicationrecord.h"
#include "composition.h"
#include "experiment.h"
#include "molpropobservable.h"
#include "phase.h"
#include "poldata/poldata.h"
#include "topology.h"

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
    double                            mass_         = 0.0;
    int                               charge_       = 0;
    int                               multiplicity_ = 1;
    //! Total number of atoms in this compound
    int                               natom_        = 0;
    //! Whether or not all atoms have proper atom types
    bool                              hasAllAtomTypes_ = false;
        std::string                       formula_, texform_, molname_, iupac_, cas_, cid_, inchi_;
        std::vector<std::string>          category_;
        //! Elemental composition of the compound
        std::map<const char *,int>        composition_;
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
        void CheckConsistency();

        //! Set the index number for sorting
        void SetIndex(int index) { index_ = index; }

        //! Return the index number for sorting
        int getIndex() const { return index_; }

        //! Set the molecular mass
        void SetMass(double mass) { mass_ = mass; }

        //! Return the molecular mass
        double getMass() const { return mass_; }

        //! Set the total charge of the molecule
        void SetTotalCharge(double charge) { charge_ = charge; }

        //! Return the total charge of the molecule
        int totalCharge() const { return charge_; }

        //! Set the multiplicity of the molecule
        void SetMultiplicity(int multiplicity) { multiplicity_ = multiplicity; }

        //! Return the multiplicity  of the molecule
        int getMultiplicity() const { return multiplicity_; }

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

        //! Set the LaTeX formula
        void SetTexFormula(const std::string &formula) { texform_.assign(formula); }

        //! Set the formula
        void SetFormula(const std::string &formula) { formula_.assign(formula); }

        //! \return the formula
        const std::string &formula() const { return formula_; }

        //! \return the LaTeX formula
        const std::string &getTexFormula() const;

        //! \brief Generate the elemental composition
        void generateComposition();
        
        //! \return The elemental composition
        const std::map<const char *,int> &composition() const { return composition_; }
        /*! \brief
         * Generate the chemical formula for this molecule based on atoms
         * present in a calculation
         *
         * \param[in] ap Data structure containing information about atoms
         * \todo Check and double check. If there is no calculation data no
         * formula can be generated
         */
        bool GenerateFormula(gmx_atomprop_t ap);

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

        /*! \brief Find a property according to specifications below
         * \param[in] mpo    The property of interest
         * \param[in] iQM    Whether Exp, QM or Both are acceptable
         * \param[in] T      The temperature at which the property
         *                   should have been measured. If -1 any T
         *                   may be returned.
         * \param[in] method The QM method, may be empty
         * \param[in] basis  The QM basisset, may be empty
         * \param[in] conf   The molecular conformation, may be empty
         * \return the property of interest or nullptr if not found
         */
        const GenericProperty *findProperty(MolPropObservable  mpo, 
                                            iqmType            iQM,
                                            double             T,
                                            const std::string &method,
                                            const std::string &basis,
                                            const std::string &conf) const;
#ifdef OLD
        //! Convenience function
        //! TODO document and clean up
        bool getPropRef(MolPropObservable mpo, iqmType iQM,
                        const std::string &method,
                        const std::string &basis,
                        const std::string &conf,
                        double *value, double *error, double *T,
                        std::string *ref, std::string *mylot,
                        std::vector<double> *vec, tensor quadrupole);

        //! And another one
        //! TODO document and clean up
        bool getProp(MolPropObservable mpo, iqmType iQM,
                     const std::string &method,
                     const std::string &basis,
                     const std::string &conf,
                     double *value, double *error, double *T);
#endif
        //! Returns true if the HF energy of the optimized geometry exists and returns the HF
        bool getOptHF(double *value);

        //! Returns the number of Opt and SP experiments for a molecule in allmols.dat
        int NOptSP();

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
        
        /*! \brief Fetch a calculation according to sepcifications
         * \param[in] method The QM method
         * \param[in] basis  The basisset
         * \param[in] conformation May be empty
         * \return the Experiment or nullptr
         */
        const Experiment *findExperimentConst(const std::string &method,
                                              const std::string &basis,
                                              const std::string &conformation) const;
        //! Add an experiment
        void AddExperiment(const Experiment &myexp) { exper_.push_back(myexp); }

        void Stats()
        {
            printf("%s - %s - %d experiments\n",
                   molname_.c_str(), formula_.c_str(),
                   (int)exper_.size());
        }
        //! Return the number of experiments
        int NExperiment() const { return exper_.size(); }

        //! Return const vector of experiments
        const std::vector<Experiment> &experimentConst() const { return exper_; }

        //! Return mutable vector of experiments
        std::vector<Experiment> *experiment() { return &exper_; }

        //! Return end of experiment vector
        ExperimentConstIterator EndExperiment()  const { return exper_.end(); }
        
        //! Return pointer to the last inserted experiment or nullptr if the number of experiments is zero
        Experiment *LastExperiment()
        {
            if (NExperiment() > 0)
            {
                return &(exper_.back());
            }
            else
            {
                return nullptr;
            }
        }

        /*! \brief Return a calculation iterator
         *
         * Return iterator corresponding to the level of theory
         * as determined by method and basis or EndExperiment in
         * case it is not found. If either method or basis are empty
         * any calculation may be taken. The calculation should
         * hold the requested observable of the type (can be nullptr)
         * \param[in]  method  The QM method
         * \param[in]  basis   The QM basis set
         * \param[out] mylot   The level of theory found
         * \param[in]  mpo     The observable type
         * \param[in]  type    The observable subtype
         */
        ExperimentIterator getCalcPropType(const std::string &method,
                                           const std::string &basis,
                                           std::string       *mylot,
                                           MolPropObservable  mpo,
                                           const char        *type);
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
