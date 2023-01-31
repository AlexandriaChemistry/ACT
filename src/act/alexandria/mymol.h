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

#ifndef ACT_MYMOL_H
#define ACT_MYMOL_H

#include <map>

#include "act/alexandria/fragmenthandler.h"
#include "act/alexandria/gentop_vsite.h"
#include "act/alexandria/molselect.h"
#include "act/alexandria/mymol_low.h"
#include "act/forces/forcecomputer.h"
#include "act/molprop/molprop.h"
#include "act/poldata/poldata.h"
#include "act/qgen/qgen_acm.h"
#include "act/qgen/qgen_resp.h"
#include "act/qgen/qtype.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/regression.h"
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
    // GromacsStuff
    t_excls                         *excls_          = nullptr;
    std::unique_ptr<gmx_vsite_t>    *vsite_          = nullptr;
    std::unique_ptr<gmx::MDAtoms>   *MDatoms_        = nullptr;
    std::unique_ptr<gmx::MDModules> *mdModules_      = nullptr;
    bool                             gromacsGenerated_ = false;
    gpp_atomtype_t                   gromppAtomtype_;
    //! GROMACS state variable
    t_state                         *state_      = nullptr;
    //! GROMACS force record
    t_forcerec                      *fr_         = nullptr;
   //! GROMACS style atoms structure
    t_atoms                         *atoms_;
    //! GROMACS symbol table
    t_symtab                        *symtab_         = nullptr;
    //! GROMACS Topology etc.
    gmx_mtop_t                      *mtop_           = nullptr;
    gmx_localtop_t                  *ltop_           = nullptr;
    t_inputrec                      *inputrec_       = nullptr;
    gmx_enerdata_t                  *enerd_          = nullptr;
    t_fcdata                        *fcd_            = nullptr;

    // The energy storage
    std::map<InteractionType, double> energies_;
    // Optimized coordinates from the input.
    std::vector<gmx::RVec>           optimizedCoordinates_;
    // Older stuff
    bool                             bHaveShells_    = false;
    bool                             bHaveVSites_    = false;
    //! Tolerance for converged SQE or EEM calculations
    double                           qTolerance_     = 1e-6;
    //! Max iterations for SQE or EEM calculations
    int                              maxQiter_       = 5;
    //! How to renumber input atom numbers to numbers with shells
    std::vector<int>                 shellRenumber_;
    //! This points to the atom indices before shells were added
    std::map<int, int>               originalAtomIndex_;
    //! The job type corresponding to the coordinates
    JobType                          myJobType_        = JobType::UNKNOWN;
    //! Map of charge type dependent properties
    std::map<qType, QtypeProps>      qProps_;
    //! Center of nuclear charge
    gmx::RVec                        CenterOfCharge_ = { 0, 0, 0 };
    //! Function that returns true if a molecule is symmetric
    bool IsSymmetric(real toler) const;
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
                       const Poldata     *pd) const;
    
    /*! \brief
     * Find the atoms inside the molcule needed to construct the inplane virtual sites.
     *
     * \param[in]  ca     The index of the central atom
     * \param[out] atoms  Data structure containing atomic properties
     */
    void findInPlaneAtoms(int ca, std::vector<int> *atoms) const;
    
    /*! \brief
     * Find the atoms inside the molcule needed to construct the out of plane virtual sites.
     *
     * \param[in]  ca     The index of the central atom
     * \param[out] atoms  Data structure containing atomic properties
     */
    void findOutPlaneAtoms(int ca, std::vector<int> *atoms) const;
    
    /*! \brief extract frequencies and intensities if present.
     * This will use the optimized structure only.
     */
    void getHarmonics();
    
    //! Energy terms for this compounds
    std::map<MolPropObservable, double> energy_;
    //! Molecular Topology
    Topology                      *topology_;
    //! All the atoms, but not shells or vsites
    std::vector<int>               realAtoms_; 
    // Reference data for devcomputer
    std::vector<double>            ref_frequencies_;
    std::vector<double>            ref_intensities_;
    t_nrnb                         nrnb_;
    gmx_wallcycle_t                wcycle_;
    gmx_shellfc_t                 *shellfc_       = nullptr;
    std::vector<int>               symmetric_charges_;
    std::vector<std::string>       error_messages_;
    eSupport                       eSupp_         = eSupport::Local;
    //! Structure to manage charge generation
    FragmentHandler               *fraghandler_   = nullptr;
    iMolSelect                     dataset_type_  = iMolSelect::Ignore;

public:

    //! \return a GROMACS style array with energy terms
    const real *energyTerms() const;
    
    /*! \brief
     * Return the original coordinate vector of the molecule. If there are shell particles
     * their position may have been optimized, but the atoms are as read from the input.
     */
    const std::vector<gmx::RVec> &xOriginal() const { return optimizedCoordinates_; }

    /*! \brief
     * Constructor
     */
    MyMol();
    
    iMolSelect datasetType() const { return dataset_type_; }
    
    void set_datasetType(iMolSelect dataset_type) 
    {
        dataset_type_ = dataset_type;
    }
    
    //! \return the job type corresponding to coordinates
    const JobType &jobType() const { return myJobType_; }
    
    //! \return true if shells are present
    bool haveShells() const { return bHaveShells_; }
    
    //! \return whether this is a linear molecule
    bool linearMolecule() const;
    
    //! \return how this compound is supported on this processor
    eSupport support() const { return eSupp_; }
    
    /*! \brief Set the support type
     * \param[in] esup The support type
     */
    void setSupport(eSupport esup) { eSupp_ = esup; }
    
    //! \return the number of real atoms, i.e. not shells or vsites
    int nRealAtoms() const { return realAtoms_.size(); }
    
    //! \return the indexex of real atoms
    const std::vector<int> &realAtoms() const { return realAtoms_; }

    //! \return the GROMACS energy data
    const gmx_enerdata_t *enerdata() const { return enerd_; }
    
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
     * vector of atoms. Energies are stored as a 1D vector pair.
     * All energy components are stored as a vector of
     * reference energies paired with a map of ACT energy components. The 
     * InteractionType index points to the energy type in the
     * energyComponentsMap.
     * \param[in]  pd                   The force field structure
     * \param[in]  forceComp            The force computer utility
     * \param[out] forceMap             The forces
     * \param[out] energyMap            The energy components for a single structure
     * \param[out] interactionEnergyMap The interaction energies
     * \param[out] energyComponentsMap  The energy components
     */
    void forceEnergyMaps(const Poldata                                                       *pd,
                         const ForceComputer                                                 *forceComp,
                         std::vector<std::vector<std::pair<double, double> > >               *forceMap,
                         std::vector<std::pair<double, double> >                             *energyMap,
                         std::vector<std::pair<double, double> >                             *interactionEnergyMap,
                         std::vector<std::pair<double, std::map<InteractionType, double> > > *energyComponentMap) const;
    
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
     * \return mdatoms structure
     */
    t_mdatoms *getMdatoms() { return MDatoms_->get()->mdatoms(); }
    
    /*! \brief
     * \return GROMACS atoms structure
     */
    t_atoms *gmxAtoms() const { return atoms_; }
    
    /*! \brief
     * \return atoms data for editing
     */
    std::vector<ActAtom> *atoms() { return topology_->atomsPtr(); }
    
    /*! \brief Return the fragment handler
     */
    const FragmentHandler *fragmentHandler() const { return fraghandler_; }
    
    /*! \brief
     * \return atoms const structure
     */
    const t_atoms *gmxAtomsConst() const { return atoms_; }
    
    /*! \brief
     * \return atoms data
     */
    const std::vector<ActAtom> &atomsConst() const { return topology_->atoms(); }
    
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
     * \param[in]  gromacsSupport Whether or not to include legacy gromacs
     * \return status
     */
    immStatus GenerateTopology(FILE              *fp,
                               const Poldata     *pd,
                               missingParameters  missing,
                               bool               gromacsSupport);
    
    //! Return the ACT topology structure
    const Topology *topology() const { return topology_; }
    
    /*! \brief
     *  Computes polarizability tensor in the presence of external
     *  electric field. The result is stored in 
     *  mymol->qTypeProps(qType::Calc).
     *
     * \param[in] pd        The force field
     * \param[in] forceComp The force computer
     */
    void CalcPolarizability(const Poldata       *pd,
                            const ForceComputer *forceComp);
    
    /*! \brief
     * Generate atomic partial charges
     *
     * \param[in]  pd        Data structure containing atomic properties
     * \param[in]  forceComp Force computer utility
     * \param[in]  fplog     Logger
     * \param[in]  cr        Communication parameters
     * \param[in]  algorithm The algorithm for determining charges,
     *                       if NONE it is read from the Poldata structure.
     * \param[in]  qtype     If algorithm is Read this type of charges will
     *                       be extracted from the molprop structure.
     * \param[in]  qcustom   Custom (user-provided) charges
     * \param[out] coords    The coordinates, will be updated for shells
     * \param[out] forces    This routine will compute energies and forces.
     */
    immStatus GenerateCharges(const Poldata             *pd,
                              const ForceComputer       *forceComp,
                              const gmx::MDLogger       &fplog,
                              const CommunicationRecord *cr,
                              ChargeGenerationAlgorithm  algorithm,
                              qType                      qtype,
                              const std::vector<double> &qcustom,
                              std::vector<gmx::RVec>    *coords,
                              std::vector<gmx::RVec>    *forces);
    /*! \brief
     * Generate atomic partial charges using EEM or SQE.
     * If shells are present they will be minimized.
     *
     * \param[in]  pd        Data structure containing atomic properties
     * \param[in]  forceComp The force computer
     * \param[out] coords    The coordinates, will be updated for shells
     * \param[out] forces    The forces
     */
    immStatus GenerateAcmCharges(const Poldata          *pd,
                                 const ForceComputer    *forceComp,
                                 std::vector<gmx::RVec> *coords,
                                 std::vector<gmx::RVec> *forces);
    
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
     * \param[in] pd     Poldata structure
     * \param[in] coords Coordinates
     */
    void calcEspRms(const Poldata                *pd,
                    const std::vector<gmx::RVec> *coords);
    
    /*! \brief Make a ESP correlation plot
     *
     * Will compute the electrostatic potential around the compound
     * and make a correlation plot between the QM potential and the
     * Alexandria potential.
     * \param[in] pd        Force field
     * \param[in] coords    The coordinates
     * \param[in] espcorr   File name to plot to
     * \param[in] oenv      Gromacs output structure
     * \param[in] forceComp Utility to compute forces
     */
    void plotEspCorrelation(const Poldata                *pd,
                            const std::vector<gmx::RVec> &coords,
                            const char                   *espcorr,
                            const gmx_output_env_t       *oenv,
                            const ForceComputer          *forceComp);
    
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
     * \param[in] forceComp The force computer utility
     * \param[in] cr        Communication record
     * \param[in] coords    The coordinates
     * \param[in] method    QC method
     * \param[in] basis     QC basis set
     * \param[in] bITP      Whether or not to write an itp file iso top file
     */
    void PrintTopology(const char                   *fn,
                       bool                          bVerbose,
                       const Poldata                *pd,
                       const ForceComputer          *forceComp,
                       const CommunicationRecord    *cr,
                       const std::vector<gmx::RVec> &coords,
                       const std::string            &method,
                       const std::string            &basis,
                       bool                          bITP = false);
    
    //! \brief Update GROMACS data structures
    void updateMDAtoms();
    
    /*! \brief Calculate the forces and energies
     * For a polarizable model the shell positions are minimized.
     * This code is maintained only for comparing ACT native energies and forces
     * to the gromacs code. Do not use inproduction code.
     * \param[in]  crtmp         Temporary communication record with one core only.
     * \param[in]  coordinates   The atomic coordinates
     * \param[out] forces        Force array
     * \param[out] energies      The energy components
     * \param[out] shellForceRMS Root mean square force on the shells
     * \return immStatus::OK if everything worked fine, error code otherwise.
     */
    immStatus calculateEnergyOld(const t_commrec                   *crtmp,
                                 std::vector<gmx::RVec>            *coordinates,
                                 PaddedVector<gmx::RVec>           *forces,
                                 std::map<InteractionType, double> *energies,
                                 real                              *shellForceRMS);
    
    /*! \brief Calculate the interaction energies.
     * For a system with multiple fragments this will compute
     * Epot(system) - Sum_f Epot(f) 
     * where f are the fragments.
     * For a polarizable model the shell positions are minimized.
     * \param[in] pd                 The force field
     * \param[in] forceComputer      The code to run the calculations.
     * \param[out] interactionForces The forces on the atoms due to the interacting components
     * \param[inout] coords          Atomic coordinates (shell positions can be updated)
     * \return The interaction energy.
     */
    double calculateInteractionEnergy(const Poldata          *pd,
                                      const ForceComputer    *forceComputer,
                                      std::vector<gmx::RVec> *interactionForces,
                                      std::vector<gmx::RVec> *coords) const;
    
    /*! \brief
     * Update internal structures for bondtype due to changes in pd
     *
     * \param[in] pd     Data structure containing atomic properties
     * \param[in] iTypes Interaction types
     * \param[in] updateZeta Whether to update the atomic zeta as well
     */
    void UpdateIdef(const Poldata                      *pd,
                    const std::vector<InteractionType> &iTypes,
                    bool                                updateZeta);
    
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
     * \param[in] coords      Atomic coordinates
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
                      const std::vector<gmx::RVec> &coords,
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
     * \param[in] fn          The filename.
     * \param[in] coords      The atomic coordinates
     * \param[in] writeShells Whether or not to write the shell particles
     */
    void PrintConformation(const char                   *fn,
                           const std::vector<gmx::RVec> &coords,
                           bool                          writeShells);
    
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
     * \param[in]  coords  The coordinates
     * \param[in]  watoms  Weight for the potential on the atoms in 
     *                     doing the RESP fit. Should be 0 in most cases.
     * \param[in]  maxESP  Percentage of the ESP points to consider (<= 100)
     */
    void initQgenResp(const Poldata                *pd,
                      const std::vector<gmx::RVec> &coords,
                      real                          watoms,
                      int                           maxESP);
    
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr, int dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr   data structure for MPI communication
     * \param[in] root The MPI root
     * \param[in] comm The MPI communicator
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
    CommunicationStatus Receive(const CommunicationRecord *cr, int src);
    
};

}

#endif
