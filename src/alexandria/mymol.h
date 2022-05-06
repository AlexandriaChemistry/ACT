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

#ifndef ACT_MYMOL_H
#define ACT_MYMOL_H

#include <map>

#include "act/molprop/molprop.h"
#include "act/poldata/poldata.h"
#include "act/qgen/qgen_acm.h"
#include "act/qgen/qgen_resp.h"
#include "act/qgen/qtype.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/regression.h"
#include "alexandria/fragmenthandler.h"
#include "alexandria/gentop_vsite.h"
#include "alexandria/molselect.h"
#include "alexandria/mymol_low.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

struct gmx_enerdata_t;
struct gmx_shellfc_t;
struct t_forcerec;
struct t_inputrec;
struct gmx_hw_info_t;
struct gmx_vsite_t;
struct t_nrnb;

namespace alexandria
{

//! Determine on which core/processor this compound is supported
enum class eSupport {
    //! Not supported at all
    No,
    //! Supported on this core
    Local,
    //! Supported on another core
    Remote
};

//! Distinguish different sets of coordinates
enum class coordSet {
    //! The original coordinates
    Original,
    //! The minimized coordinates
    Minimized
};

    //! Forward declaration, full declaration is in mymol.cpp
    class MyForceProvider;
    /*! \brief
     * Contains molecular properties from a range of sources.
     * Overloads the regular molprop and adds a lot of functionality.
     * For one thing, it can generate molprop contents from a coordinate
     * file if needed.
     *
     * \inpublicapi
     * \ingroup module_alexandria
     */
    class MyMol : public MolProp
    {
    private:
        // Now class MolHandler can access private members of MyMol
        friend class MolHandler;
        int                             *cgnr_           = nullptr;
        bool                             bHaveShells_    = false;
        bool                             bHaveVSites_    = false;
        bool                             bNeedVsites_    = false;
        //! Tolerance for converged SQE or EEM calculations
        double                           qTolerance_     = 1e-6;
        //! Max iterations for SQE or EEM calculations
        int                              maxQiter_       = 5;
        t_excls                         *excls_          = nullptr;
        std::unique_ptr<gmx_vsite_t>    *vsite_          = nullptr;
        std::unique_ptr<gmx::MDAtoms>   *MDatoms_        = nullptr;
        std::unique_ptr<gmx::MDModules> *mdModules_      = nullptr;
        MyForceProvider                 *myforce_        = nullptr;
        //! How to renumber input atom numbers to numbers with shells
        std::vector<int>                 shellRenumber_;
        //! This points to the atom indices before shells were added
        std::map<int, int>               originalAtomIndex_;
        bool                             gromacsGenerated_ = false;
        gpp_atomtype_t                   gromppAtomtype_;
        double                           atomizationEnergy_ = 0;
        //! The job type corresponding to the coordinates
        JobType                          myJobType_        = JobType::UNKNOWN;
        //! Map of charge type dependent properties
        std::map<qType, QtypeProps>           qProps_;
        //! Center of nuclear charge
        rvec                      CenterOfCharge_ = { 0 };
        //! GROMACS state variable
        t_state                  *state_      = nullptr;
        //! GROMACS force record
        t_forcerec               *fr_         = nullptr;
        //! Function that returns true if a molecule is symmetric
        bool             IsSymmetric(real toler);
        //! Map of vectors to back up coordinates
        std::map<coordSet, std::vector<gmx::RVec> > backupCoordinates_;
        /*! Make a back up of coordinates
         * \param[in] cs Which data set to back up to
         */
        void backupCoordinates(coordSet cs);
        /*! Restored backed up coordinates
         * \param[in] cs Which data set to restore from
         */
        void restoreCoordinates(coordSet cs);
        /*! \brief
         * Generate Atoms based on quantum calculation with specified level of theory.
         * If the requested level of theory is not present, another
         * level can be tried if the strict flag is false.
         * \param[in]  pd     Force field data
         * \param[out] atoms  The structure to update
         * \return The status
         */
        immStatus GenerateAtoms(const Poldata     *pd,
                                t_atoms           *atoms);

        /*! \brief
         * Generate angles, dihedrals, exclusions etc.
         *
         * \param[in] atoms  The atoms structure
         * \param[in] bDihs  Whether to generate dihedrals
         * \param[in] nrexcl Number of exclusions for use in GROMACS code
         */
        void MakeAngles(t_atoms *atoms,
                        bool     bDihs,
                        int      nrexcl);

        /*! \brief
         * Generate virtual sites or linear angles
         *
         * \param[in] pd         Poldata
         * \param[in] atoms      Atoms structure
         */
        //void MakeSpecialInteractions(const Poldata *pd,
        //                             t_atoms       *atoms);
        /*! \brief
         * Add vsites on bonds to hos bond shell particles
         *
         * \param[in]  fp    File to write (debug) information to
         * \param[in]  pd    Data structure containing atomic properties
         * \param[out] atoms Structure to modify with new particles.
         */
        void addBondVsites(FILE          *fp,
                           const Poldata *pd,
                           t_atoms       *atoms);

        /*! \brief
         * Add shell particles
         *
         * \param[in]  fp    File to write (debug) information to
         * \param[in]  pd    Data structure containing atomic properties
         * \param[out] atoms Structure to modify with new particles.
         */
        void addShells(FILE          *fp,
                       const Poldata *pd,
                       t_atoms       *atoms);

        /*! \brief
         * Check whether atom types exist in the force field
         * also check whether the multiplicity is correct.
         *
         * \param[in] pd    The force field structure
         * \param[in] atoms The structure to check
         * \return status code.
         */
        immStatus checkAtoms(const Poldata *pd,
                             const t_atoms *atoms);

        /*! \brief
         * Return true if atom type needs to have virtual site.
         *
         * \param[in] atype  Atom type
         * \param[in] pd     Data structure containing atomic properties
         * \return true if vsite is needed
         */
        bool IsVsiteNeeded(std::string        atype,
                           const Poldata     *pd);

        /*! \brief
         * Find the atoms inside the molcule needed to construct the inplane virtual sites.
         *
         * \param[in]  ca     The index of the central atom
         * \param[out] atoms  Data structure containing atomic properties
         */
        void findInPlaneAtoms(int ca, std::vector<int> &atoms);

        /*! \brief
         * Find the atoms inside the molcule needed to construct the out of plane virtual sites.
         *
         * \param[in]  ca     The index of the central atom
         * \param[out] atoms  Data structure containing atomic properties
         */
        void findOutPlaneAtoms(int ca, std::vector<int> &atoms);

        /*! \brief Compute the atomization energy
         * Will add the reference enthalpy values per atom
         * in the molecule and store it locally.
         * \param[in] pd The force field structure.
         */
        void computeAtomizationEnergy(const Poldata *pd);

        /*! \brief extract frequencies and intensities if present.
         * This will use the optimized structure only.
         */
        void getHarmonics();
        
        //! Energy terms for this compounds
        std::map<MolPropObservable, double> energy_;
        //! Molecular Topology
        Topology                      *topology_;
        gmx_mtop_t                    *mtop_           = nullptr;
        gmx_localtop_t                *ltop_           = nullptr;
        t_symtab                      *symtab_         = nullptr;
        t_inputrec                    *inputrec_       = nullptr;
        gmx_enerdata_t                *enerd_          = nullptr;
        t_fcdata                      *fcd_            = nullptr;
        int                            nRealAtoms_     = 0;
        PaddedVector<gmx::RVec>        f_;
        PaddedVector<gmx::RVec>        optf_;
        // Reference data for devcomputer
        std::vector<double>            ref_frequencies_;
        std::vector<double>            ref_intensities_;
        t_nrnb                         nrnb_;
        gmx_wallcycle_t                wcycle_;
        gmx_shellfc_t                 *shellfc_        = nullptr;
        std::vector<int>               symmetric_charges_;
        std::vector<std::string>       error_messages_;
        eSupport                       eSupp_         = eSupport::Local;
        //! Structure to manage charge generation
        FragmentHandler               *fraghandler_   = nullptr;
   public:
  
        /*! \brief
         * Constructor
         */
        MyMol();

        iMolSelect                     dataset_type_   = iMolSelect::Ignore;

        iMolSelect datasetType() const { return dataset_type_; }

        void set_datasetType(iMolSelect dataset_type) 
        {
            dataset_type_ = dataset_type;
        }

        //! \return the job type corresponding to coordinates
        const JobType &jobType() const { return myJobType_; }

        //! \return true if shells are present
        bool haveShells() const { return nullptr != shellfc_; }

        //! \return whether this is a linear molecule
        bool linearMolecule() const;
        
        //! \return how this compound is supported on this processor
        eSupport support() const { return eSupp_; }

        /*! \brief Set the support type
         * \param[in] esup The support type
         */
        void setSupport(eSupport esup) { eSupp_ = esup; }
        
        //! \return the number of real atoms, i.e. not shells or vsites
        int nRealAtoms() const { return nRealAtoms_; }

        //! \return the atomization energy
        double atomizationEnergy() const { return atomizationEnergy_; }
        
        /*! Return an energy component
         * \param[in]  mpo  The particular term that is requested
         * \param[out] ener The value found
         * \return whether or not the energy is found
         */
        bool energy(MolPropObservable mpo, double *ener) const 
        {
            if (energy_.find(mpo) == energy_.end())
            {
                return false;
            }
            *ener = energy_.find(mpo)->second;
            return true;
        }

        /*! \brief Compute energies and forces for all structures
         * The optimized structure as well as all the excited
         * structures will be used to compute the energy and
         * forces. Store the results in vectors containing two
         * doubles, first the reference, then the calculated one.
         * Forces are stored as a vector of structures and then a
         * vector of atoms. Energies are stored as a 1D vector.
         * \param[out] forceMap The forces
         * \param[out] enerMap  The energies
         */
        void forceEnergyMaps(std::vector<std::vector<std::pair<double, double> > > *forceMap,
                             std::vector<std::pair<double, double> > *enerMap);
        
        //! Return the reference frequencies collected earlier
        const std::vector<double> &referenceFrequencies() const { return ref_frequencies_; }
                
        //! Return the reference intensities collected earlier
        const std::vector<double> &referenceIntensities() const { return ref_intensities_; }
                
        //! \return the center of charge
        const rvec &centerOfCharge() const { return CenterOfCharge_; }

        /*! \brief Return QtypeProps for a charge type
         * \param[in] qt The charge type, e.g. qType::CM5
         * \return the corresponding structure or nullptr 
         */
        QtypeProps *qTypeProps(qType qt);
        
        /*! \brief Return QtypeProps for a charge type
         * \param[in] qt The charge type, e.g. qType::CM5
         * \return the corresponding structure or nullptr 
         */
        const QtypeProps *qTypeProps(qType qt) const;
        
        const std::vector<std::string> &errors() const {return error_messages_;}

        /*! \brief
         * Add the screening factors of the distributed charge to atom structure
         *
         * \param[in] pd     Data structure containing atomic properties
         * \param[out] atoms Structure to fill with force field data
         */
        immStatus zetaToAtoms(const Poldata *pd,
                              t_atoms       *atoms);

        /*! \brief
         * Return the coordinate vector of the molecule
         */
        const gmx::HostVector<gmx::RVec> &x() const { return state_->x; }

        /*! \brief
         * Return the force vector of the molecule
         */
        const PaddedVector<gmx::RVec> &f() const { return f_; }

        //! \return the potential energy of this molecule
        real potentialEnergy() const;

        //! \return a GROMACS style array with energy terms
        const real *energyTerms() const;
        
        /*! \brief
         * \return mdatoms structure
         */
        t_mdatoms *getMdatoms() { return MDatoms_->get()->mdatoms(); }

        /*! \brief
         * \return atoms structure for editing
         */
        t_atoms *atoms();
        
        /*! \brief Return the fragment handler
         */
        const FragmentHandler *fragmentHandler() const { return fraghandler_; }
        
        /*! \brief
         * \return atoms const structure
         */
        const t_atoms &atomsConst() const;
        
        /*! \brief Return the bond order
         * \param[in] ai Atom I
         * \param[in] aj Atom J
         * \return bond order or 0 if not present
         */
        double bondOrder(int ai, int aj) const;

        /*! \brief
         * It generates the atoms structure which will be used to print the topology file.
         *
         * \param[in]  fp      File to write (debug) information to
         * \param[in]  pd      Data structure containing atomic properties
         * \param[in]  missing How to treat missing parameters
         * \return status
         */
        immStatus GenerateTopology(FILE              *fp,
                                   const Poldata     *pd,
                                   missingParameters  missing);

        //! Return the ACT topology structure
        const Topology *topology() const { return topology_; }

        /*! \brief
         *  Computes polarizability tensor in the presence of external
         *  electric field. The result is stored in 
         *  mymol->qTypeProps(qType::Calc).
         *
         * \param[in]  efield   Strenght of the external electric field
         * \returns the result of the calculation, if fine it is immOK
         */
        immStatus CalcPolarizability(double efield);
        
        /*! \brief
         * Generate atomic partial charges
         *
         * \param[in] pd        Data structure containing atomic properties
         * \param[in] fplog     Logger
         * \param[in] cr        Communication parameters
         * \param[in] algorithm The algorithm for determining charges,
         *                      if NONE it is read from the Poldata structure.
         * \param[in] qcustom   Custom (user-provided) charges
         */
        immStatus GenerateCharges(const Poldata             *pd,
                                  const gmx::MDLogger       &fplog,
                                  const CommunicationRecord *cr,
                                  ChargeGenerationAlgorithm  algorithm,
                                  const std::vector<double> &qcustom);
        /*! \brief
         * Generate atomic partial charges using EEM or SQE.
         * If shells are present they will be minimized.
         *
         * \param[in] pd      Data structure containing atomic properties
         */
        immStatus GenerateAcmCharges(const Poldata *pd);
                                     
        /*! \brief Implement charge symmetrization
         *
         * Initiates internal structure for atom charge symmetry
         * (e.g. CH3 with identical charges on H).
         * Must be called before initQgenResp.
         * \param[in] pd                 Data structure containing atomic properties
         * \param[in] bSymmetricCharges  Consider molecular symmetry to calculate partial charge
         * \param[in] symm_string        The type of molecular symmetry
         */
        void symmetrizeCharges(const Poldata  *pd,
                               bool            bSymmetricCharges,
                               const char     *symm_string);

        /*! \brief Calculate the RMSD from ESP for QM charges
         * \param[in]  pd     Poldata structure
         */
        void calcEspRms(const Poldata *pd);

        /*! \brief Make a ESP correlation plot
         *
         * Will compute the electrostatic potential around the compound
         * and make a correlation plot between the QM potential and the
         * Alexandria potential.
         * \param[in] espcorr File name to plot to
         * \param[in] oenv    Gromacs output structure
         */
        void plotEspCorrelation(const char                *espcorr,
                                const gmx_output_env_t    *oenv);

        /*! \brief
         * Collect the properties from QM (Optimized structure) or
         * from experiment.
         *
         * \param[in] iqm      Determine whether to allow exp or QM results or both for each property
         * \param[in] T        The temperature to use. -1 allows any T.
         */
        immStatus getExpProps(const std::map<MolPropObservable, iqmType> &iqm,
                              double             T);

        /*! \brief
         * Print the topology that was generated previously in GROMACS format.
         *
         * \param[in] fn        File name
         * \param[in] bVerbose  Verbose
         * \param[in] pd        Data structure containing atomic properties
         * \param[in] cr        Gromacs communication record
         * \param[in] method    QC method
         * \param[in] basis     QC basis set
         * \param[in] bITP      Whether or not to write an itp file iso top file
         */
        void PrintTopology(const char                *fn,
                           bool                       bVerbose,
                           const Poldata             *pd,
                           const CommunicationRecord *cr,
                           const std::string         &method,
                           const std::string         &basis,
                           bool                       bITP = false);

        //! \brief Update GROMACS data structures
        void updateMDAtoms();
        
        /*! \brief Calculate the forces and energies
         * For a polatizable model the shell positions are minimized.
         * \param[in]  crtmp         Temporary communication record with one core only.
         * \param[out] shellForceRMS Root mean square force on the shells
         * \return immStatus::OK if everything worked fine, error code otherwise.
         */
        immStatus calculateEnergy(const t_commrec *crtmp,
                                  real            *shellForceRMS);

        /*! \brief
         * Relax the shells (if any) or compute the forces in the molecule.
         *
         * \param[out] rmsf                Root mean square force on the shells
         * \return immOK if everything went fine, an error otherwise.
         */
        immStatus computeForces(double *rmsf);

        /*! \brief
         * Return the optimized geometry of the molecule from the data file.
         * The structure returned here corresponds to the one from the input
         * file.
         */
        bool getOptimizedGeometry(rvec *x);

        /*! \brief
         * Update internal structures for bondtype due to changes in pd
         *
         * \param[in] pd           Data structure containing atomic properties
         * \param[in] iType        Interaction type
         */
        void UpdateIdef(const Poldata   *pd,
                        InteractionType  iType);

        /*! \brief
         * Generate GROMACS structures.
         */
        immStatus GenerateGromacs(const gmx::MDLogger       &mdlog,
                                  const CommunicationRecord *cr,
                                  const char                *tabfn,
                                  ChargeType                 iType);

        /*! \brief
         * Generate cube
         *
         * \param[in] pd          Data structure containing atomic properties
         * \param[in] spacing     The grid space
         * \param[in] border      The amount of space around the molecule
         * \param[in] reffn
         * \param[in] pcfn
         * \param[in] pdbdifffn
         * \param[in] potfn
         * \param[in] rhofn
         * \param[in] hisfn
         * \param[in] difffn
         * \param[in] diffhistfn
         * \param[in] oenv
         */
        void GenerateCube(const Poldata          *pd,
                          real                    spacing,
                          real                    border,
                          const char             *reffn,
                          const char             *pcfn,
                          const char             *pdbdifffn,
                          const char             *potfn,
                          const char             *rhofn,
                          const char             *hisfn,
                          const char             *difffn,
                          const char             *diffhistfn,
                          const gmx_output_env_t *oenv);

        /*! \brief
         * Print the coordinates corresponding to topology after adding shell particles and/or vsites.
         *
         * \param[in] fn The filename.
         */
        void PrintConformation(const char *fn);

        /*! \brief
         * set the inputrec of the MyMol object
         *
         * \param[in] ir   GROMACS t_inputrec structure
         */
        void setInputrec(t_inputrec  *ir)
        {
            inputrec_ = ir;
        }
        
        /*! \brief
         * Equal operator for MyMol object
         *
         */
        bool operator==(const MyMol &mol) const
        {
            return (this->getMolname() == mol.getMolname());
        }
        
        /*! \brief Init the class for the RESP algorithm.
         * This must be called before generateCharges even if no RESP
         * is used for charge generation, since ESP points can be used
         * in analysis.
         * \param[in]  pd      Data structure containing atomic properties
         * \param[in]  watoms  Weight for the potential on the atoms in 
         *                     doing the RESP fit. Should be 0 in most cases.
         * \param[in]  maxESP  Percentage of the ESP points to consider (<= 100)
         */
        void initQgenResp(const Poldata     *pd,
                          real               watoms,
                          int                maxESP);

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] cr   GROMACS data structure for MPI communication
         * \param[in] dest Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(const CommunicationRecord *cr, int dest) const;

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] cr  GROMACS data structure for MPI communication
         * \param[in] src Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(const CommunicationRecord *cr, int src);

    };

}

#endif
