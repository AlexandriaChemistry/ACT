/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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

#ifndef ACT_ACTMOL_H
#define ACT_ACTMOL_H

#include <map>

#include "act/alexandria/fragmenthandler.h"
#include "act/alexandria/molselect.h"
#include "act/basics/msg_handler.h"
#include "act/forces/forcecomputer.h"
#include "act/molprop/molprop.h"
#include "act/forcefield/forcefield.h"
#include "act/qgen/qgen_acm.h"
#include "act/qgen/qgen_resp.h"
#include "act/qgen/qtype.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/regression.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

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

/*! \brief Class to hold corresponding QM (reference) and ACT energies
 */
class ACTEnergy
{
private:
    //! Id of the experiment
    int    experiment_ = 0;
    //! QM energy
    double eqm_        = 0;
    //! ACT energy
    double eact_       = 0;
    //! QM present
    bool haveQM_       = false;
    //! ACT present
    bool haveACT_      = false;
public:
    /*! \brief Constructor
     * \param[in] experiment The id
     * \param[in] eqm        The QM energy
     * \param[in] eact       The ACT energy
     */
    ACTEnergy(int experiment, double eqm, double eact) : experiment_(experiment)
    {
        setQM(eqm);
        setACT(eact);
    }
    /*! \brief Constructor
     * \param[in] experiment The id
     */
    ACTEnergy(int experiment) : experiment_(experiment) {}

    //! \return The experiment ID
    int id() const { return experiment_; }

    //! Is there a QM energy
    bool haveQM() const { return haveQM_; }

    //! \return The QM energy    
    double eqm() const
    {
        if (!haveQM_)
        {
            GMX_THROW(gmx::InvalidInputError("Trying to extract a QM energy when there is none"));
        }
        return eqm_;
    }

    //! set QM value
    void setQM(double eqm)
    {
        eqm_ = eqm;
        haveQM_ = true;
    }

    //! Is there an ACT energy
    bool haveACT() const { return haveACT_; }

    //! \return The ACT energy
    double eact() const
    {
        GMX_RELEASE_ASSERT(haveACT_, "Trying to extract an ACT energy when there is none");
        return eact_;
    }

    //! set ACT value
    void setACT(double eact)
    {
        eact_ = eact;
        haveACT_ = true;
    }
};

typedef std::map<InteractionType, ACTEnergy> ACTEnergyMap;
typedef std::vector<ACTEnergyMap>            ACTEnergyMapVector;

/*! \brief Class to hold corresponding QM (reference) adn ACT props
 */
class ACTQprop
{
private:
    //! Whether or not a particle is a real atom
    std::vector<bool>  isAtom_;
    //! Reference Qprops
    QtypeProps         qPqm_;
    //! Calculated Qprops
    QtypeProps         qPact_;
    //! Resp calculation structure
    QgenResp           QgenResp_;
public:
    //! Default constructor
    ACTQprop() {}

    /*! \brief Constructor initializing both QtypeProps
     * \param[in] atoms The atoms
     * \param[in] x     The coordinates
     */
    ACTQprop(const std::vector<ActAtom>   &atoms,
             const std::vector<gmx::RVec> &x);

    //! \return constant QM property    
    const QtypeProps &qPqmConst() const { return qPqm_; }

    //! \return constant ACT property    
    const QtypeProps &qPactConst() const { return qPact_; }

    //! \return QM property for editing
    QtypeProps *qPqm() { return &qPqm_; }

    //! \return ACT property for editing
    QtypeProps *qPact() { return &qPact_; }

    /*! \brief Return internal structure
     * \return the QgenResp_ data structure
     */
    QgenResp *qgenResp() { return &QgenResp_; }
    
    /*! \brief Return internal structure
     * \return the QgenResp_ data structure
     */
    const QgenResp &qgenRespConst() const { return QgenResp_; }

    //! Copy the charges from Resp to the QtypeProps
    void copyRespQ();
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
class ACTMol : public MolProp
{
private:
    //! Periodic box.
    matrix                            box_;
    // The energy storage
    std::map<InteractionType, double> energies_;
    // Optimized coordinates from the input.
    //std::vector<gmx::RVec>           optimizedCoordinates_;
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
    //! Vector of charge type dependent properties
    std::vector<ACTQprop>            qProps_;
    //! Is this a linear molecule? Will be checked first when generating a topology.
    bool                             isLinear_ = false;
    /*! \brief
     * Add vsites on bonds to hos bond shell particles
     *
     * \param[in]  fp    File to write (debug) information to
     * \param[in]  pd    Data structure containing atomic properties
     */
    void addBondVsites(FILE             *fp,
                       const ForceField *pd);
    
    /*! \brief
     * Add shell particles
     *
     * \param[in]  fp    File to write (debug) information to
     * \param[in]  pd    Data structure containing atomic properties
     * \param[out] atoms Structure to modify with new particles.
     */
    void addShells(FILE             *fp,
                   const ForceField *pd);
    
    /*! \brief
     * Check whether atom types exist in the force field
     * also check whether the multiplicity is correct.
     *
     * \param[in] msghandler Message handler
     * \param[in] pd    The force field structure
     */
    void checkAtoms(MsgHandler       *msghandler,
                    const ForceField *pd);
    
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
    Topology                      *topology_      = nullptr;
    //! All the atoms, but not shells or vsites
    std::vector<int>               realAtoms_; 
    //! Reference data for devcomputer
    std::vector<double>            ref_frequencies_;
    std::vector<double>            ref_intensities_;
    //! QM method used for data
    std::string                    method_;
    //! QM basis set used for data
    std::string                    basis_;
    //! Charge symmetry 
    std::vector<int>               symmetric_charges_;
    eSupport                       eSupp_         = eSupport::Local;
    //! Structure to manage charge generation
    FragmentHandler               *fraghandler_   = nullptr;
    iMolSelect                     dataset_type_  = iMolSelect::Ignore;
    //! Internal storage for original coords
    std::vector<gmx::RVec>         xOriginal_;
public:

    //! \return a GROMACS style array with energy terms
    const real *energyTerms() const;
    
    /*! \brief
     * Return a coordinate vector of the molecule corresponding to the first experiment
     * with Jobtype Opt or Topology or SP. The array includes shells and/or vsites.
     * \throws if not suitable experiment is present.
     */
    std::vector<gmx::RVec> xOriginal() const;

    /*! \brief
     * Constructor
     */
    ACTMol();

    //! \return the molprop object
    MolProp *molProp() { return MolProp::self(); }

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
    bool linearMolecule() const { return isLinear_; }
    
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
    //const gmx_enerdata_t *enerdata() const { return enerd_; }
    
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
     * \param[in]  msghandler                 A message handler
     * \param[in]  pd                         The force field structure
     * \param[in]  forceComp                  The force computer utility
     * \param[out] forceMap                   The forces
     * \param[out] energyMap                  The energy components for a single structure
     * \param[out] interactionEnergyMap       The interaction energies
     * \param[out] energyComponentsMap        The energy components
     * \param[in] separateInductionCorrection Whether to store InductionCorrection separately or add it to Induction
     */
    void forceEnergyMaps(MsgHandler                                                        *msghandler,
                         const ForceField                                                  *pd,
                         const ForceComputer                                               *forceComp,
                         std::vector<std::vector<std::pair<double, double> > >             *forceMap,
                         std::vector<ACTEnergy>                                            *energyMap,
                         ACTEnergyMapVector                                                *interactionEnergyMap,
                         std::vector<std::pair<double, std::map<InteractionType, double>>> *energyComponentMap,
                         bool                                                               separateInductionCorrection) const;
    
    //! Return the reference frequencies collected earlier
    const std::vector<double> &referenceFrequencies() const { return ref_frequencies_; }
    
    //! Return the reference intensities collected earlier
    const std::vector<double> &referenceIntensities() const { return ref_intensities_; }

    /*! \brief
     * \return atoms data for editing
     */
    std::vector<ActAtom> *atoms() { return topology_->atomsPtr(); }

    /*! \brief
     * \return atoms data for reading only
     */
    const std::vector<ActAtom> &atomsConst() { return topology_->atoms(); }

    //! \return the QtypeProps vector for editing
    std::vector<ACTQprop> *qProps() { return &qProps_; }

    //! \return the QtypeProps vector for reading
    const std::vector<ACTQprop> &qPropsConst() const { return qProps_; }

    //! \return the level of theory used for data
    std::string levelOfTheory();

    /*! \brief Set QM basis set 
     * \param[in] basis Will be new basis set if not empty
     */
    void setBasisset(const std::string &basis) { if (!basis.empty()) basis_ = basis; }

    /*! \brief Set QM method 
     * \param[in] method Will be new method if not empty
     */
    void setMethod(const std::string &method) { if (!method.empty()) method_ = method; }

    /*! Whether data is present
     * \param[in] mpo The observable to look for
     * \return true if any data of the type is present 
     */
    bool hasMolPropObservable(MolPropObservable mpo) const;

    /*! \brief Return the fragment handler
     */
    FragmentHandler *fragmentHandler() const { return fraghandler_; }
    
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
     * \param[in]  msghandler Message handler
     * \param[in]  pd         Data structure containing atomic properties
     * \param[in]  missing    How to treat missing parameters
     */
    void GenerateTopology(MsgHandler        *msghandler,
                          ForceField        *pd,
                          missingParameters  missing);
    
    //! Return the ACT topology structure
    const Topology *topology() const { return topology_; }
    
    //! Return the ACT topology structure for editing
    Topology *topologyPtr() { return topology_; }
    /*! \brief
     * Update internal structures with electric moments etc.
     * \param[in]  pd        Data structure containing atomic properties
     * \param[in]  forceComp Force computer utility
     * \param[out] forces    This routine will compute energies and forces.
     */
    void updateQprops(const ForceField          *pd,
                      const ForceComputer       *forceComp,
                      std::vector<gmx::RVec>    *forces);
    /*! \brief
     * Generate atomic partial charges
     *
     * \param[in]  msghandler Message Handler
     * \param[in]  pd        Data structure containing atomic properties
     * \param[in]  forceComp Force computer utility
     * \param[in]  algorithm The algorithm for determining charges,
     *                       if NONE it is read from the ForceField structure.
     * \param[in]  qtype     If algorithm is Read this type of charges will
     *                       be extracted from the molprop structure.
     * \param[in]  qcustom   Custom (user-provided) charges
     * \param[out] coords    The coordinates, will be updated for shells
     * \param[out] forces    This routine will compute energies and forces.
     * \param[in]  updateQprops Whether or not to update the qprops (dipoles, quadrupoles etc.)
     */
    void GenerateCharges(MsgHandler                *msghandler,
                         const ForceField          *pd,
                         const ForceComputer       *forceComp,
                         ChargeGenerationAlgorithm  algorithm,
                         qType                      qtype,
                         const std::vector<double> &qcustom,
                         std::vector<gmx::RVec>    *coords,
                         std::vector<gmx::RVec>    *forces,
                         bool                       updateQProps = false);
    /*! \brief
     * Generate atomic partial charges using EEM or SQE.
     * If shells are present they will be minimized.
     *
     * \param[in]  msghandler Message and status handler 
     * \param[in]  pd         Data structure containing atomic properties
     * \param[in]  forceComp  The force computer
     * \param[out] coords     The coordinates, will be updated for shells
     * \param[out] forces     The forces
     */
    ACTMessage GenerateAcmCharges(MsgHandler             *msg_handler,
                                  const ForceField       *pd,
                                  const ForceComputer    *forceComp,
                                  std::vector<gmx::RVec> *coords,
                                  std::vector<gmx::RVec> *forces);
    
    /*! \brief
     * Collect the properties from QM (Optimized structure) or
     * from experiment.
     *
     * \param[in] msghandler Message and status handler 
     * \param[in] pd   The force field   
     * \param[in] iqm  Determine whether to allow exp or QM results or both for each property
     * \param[in]  watoms  Weight for the potential on the atoms in 
     *                     doing the RESP fit. Should be 0 in most cases.
     * \param[in]  maxESP  Percentage of the ESP points to consider (<= 100)
     */
    void getExpProps(MsgHandler                                 *msghandler,
                     const ForceField                           *pd,
                     const std::map<MolPropObservable, iqmType> &iqm,
                     real                                        watoms = 0,
                     int                                         maxESP = 100);


    /*!\brief Generate molecule info
     *
     * \param[in] pd          The force field
     * \param[in] forceComp   The force computer utility
     * \param[in] coords      The coordinates
     * \return a text about the compound including properties and citation infor
     */
    std::vector<std::string> generateCommercials(const ForceField             *pd,
                                                 const ForceComputer          *forceComp,
                                                 const std::vector<gmx::RVec> &coords);

    /*! \brief
     * Print the topology that was generated previously in GROMACS format.
     *
     * \param[in] msg_handler For warning and status
     * \param[in] fn          File name
     * \param[in] pd          Data structure containing atomic properties
     * \param[in] forceComp   The force computer utility
     * \param[in] coords      The coordinates
     * \param[in] bITP        Whether or not to write an itp file iso top file
     */
    void PrintTopology(MsgHandler                   *msg_handler,
                       const char                   *fn,
                       const ForceField             *pd,
                       const ForceComputer          *forceComp,
                       const std::vector<gmx::RVec> &coords,
                       bool                          bITP = false);
    
    //! \brief Update GROMACS data structures
    void updateMDAtoms();
    
    /*! \brief Calculate the forces and energies
     * For a polarizable model the shell positions are minimized.
     * This code is maintained only for comparing ACT native energies and forces
     * to the gromacs code. Do not use in production code.
     * \param[in]  crtmp         Temporary communication record with one core only.
     * \param[in]  coordinates   The atomic coordinates
     * \param[out] forces        Force array
     * \param[out] energies      The energy components
     * \param[out] shellForceRMS Root mean square force on the shells
     * \return ACTMessage::OK if everything worked fine, error code otherwise.
     */
    ACTMessage calculateEnergyOld(const t_commrec                   *crtmp,
                                  std::vector<gmx::RVec>            *coordinates,
                                  PaddedVector<gmx::RVec>           *forces,
                                  std::map<InteractionType, double> *energies,
                                  real                              *shellForceRMS);
    
    /*! \brief Calculate the interaction energies.
     * For a system with multiple fragments this will compute
     * Epot(system) - Sum_f Epot(f) 
     * where f are the fragments. Importantly the individual components
     * of the energy are stored.
     * For a polarizable model the shell positions are minimized.
     * \param[in] msghandler                  A message handler
     * \param[in] pd                          The force field
     * \param[in] forceComputer               The code to run the calculations.
     * \param[out] einter                     The interaction energy components
     * \param[out] interactionForces          The forces on the atoms due to the interacting components
     * \param[inout] coords                   Atomic coordinates (shell positions can be updated)
     * \param[in] separateInductionCorrection Whether to store InductionCorrection separately or add it to Induction
     */
    void calculateInteractionEnergy(MsgHandler                        *msghandler,
                                    const ForceField                  *pd,
                                    const ForceComputer               *forceComputer,
                                    std::map<InteractionType, double> *einter,
                                    std::vector<gmx::RVec>            *interactionForces,
                                    std::vector<gmx::RVec>            *coords,
                                    bool                               separateInductionCorrection) const;
    
    /*! \brief
     * Update internal structures for bondtype due to changes in pd
     *
     * \param[in] pd     Data structure containing atomic properties
     * \param[in] iTypes Interaction types
     * \param[in] updateZeta Whether to update the atomic zeta as well
     */
    void UpdateIdef(const ForceField                      *pd,
                    const std::vector<InteractionType> &iTypes,
                    bool                                updateZeta);
    
    /*! \brief
     * Generate cube
     *
     * \param[in] msghandler  For debugging and information
     * \param[in] pd          Data structure containing atomic properties
     * \param[in] coords      Atomic coordinates
     * \param[in] spacing     The grid space
     * \param[in] border      The amount of space around the molecule
     * \param[in] reffn
     * \param[in] pcfn
     * \param[in] potfn
     * \param[in] rhofn
     * \param[in] hisfn
     * \param[in] difffn
     * \param[in] diffhistfn
     * \param[in] oenv
     */
    void GenerateCube(MsgHandler                   *msghandler,
                      const ForceField             *pd,
                      const std::vector<gmx::RVec> &coords,
                      const ForceComputer          *forceComp,
                      real                          spacing,
                      real                          border,
                      const char                   *reffn,
                      const char                   *pcfn,
                      const char                   *pdbdifffn,
                      const char                   *potfn,
                      const char                   *rhofn,
                      const char                   *hisfn,
                      const char                   *difffn,
                      const char                   *diffhistfn,
                      const gmx_output_env_t       *oenv);
    
    /*! \brief
     * Print the coordinates corresponding to topology after adding shell particles and/or vsites.
     *
     * \param[in] fn          The filename.
     * \param[in] coords      The atomic coordinates
     * \param[in] writeShells Whether or not to write the shell particles
     * \param[in] box         The simulation box / unitcell
     */
    void PrintConformation(const char                   *fn,
                           const std::vector<gmx::RVec> &coords,
                           bool                          writeShells,
                           const matrix                  box);
    
    /*! \brief
     * Equal operator for ACTMol object
     *
     */
    bool operator==(const ACTMol &mol) const
    {
        return (this->getMolname() == mol.getMolname());
    }
    
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
