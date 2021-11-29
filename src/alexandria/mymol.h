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

#ifndef MYMOL_H
#define MYMOL_H

#include <map>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
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

#include "gentop_core.h"
#include "gentop_vsite.h"
#include "molprop.h"
#include "molselect.h"
#include "mymol_low.h"
#include "poldata.h"
#include "qgen_acm.h"
#include "qgen_resp.h"
#include "qtype.h"

struct gmx_enerdata_t;
struct gmx_shellfc_t;
struct t_commrec;
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
        /*! \brief
         * Gromacs structures
         */
        int                             *cgnr_           = nullptr;
        bool                             bHaveShells_    = false;
        bool                             bHaveVSites_    = false;
        bool                             bNeedVsites_    = false;
        double                           ref_enthalpy_   = 0;
        double                           polarizability_ = 0;
        double                           sig_pol_        = 0;
        t_excls                         *excls_          = nullptr;
        std::unique_ptr<gmx_vsite_t>    *vsite_          = nullptr;
        std::unique_ptr<gmx::MDAtoms>   *MDatoms_        = nullptr;
        std::unique_ptr<gmx::MDModules> *mdModules_      = nullptr;
        MyForceProvider                 *myforce_        = nullptr;
        GentopVsites                     gvt_;
        std::string                      forcefield_;
        double                           isoPol_elec_      = 0;
        double                           isoPol_calc_      = 0;
        bool                             gromacsGenerated_ = false;
        gpp_atomtype_t                   gromppAtomtype_;
        //! Store the bond order for an atom pair
        std::map<std::pair<int, int>, double> bondOrder_;

        //! Map of charge type dependent properties
        std::map<qType, QtypeProps>           qProps_;
        //! Center of nuclear charge
        rvec                      CenterOfCharge_ = { 0 };
        //! Experimental dipole
        double                    dip_exp_    = 0;
        //! Error in experimental dipole
        double                    dip_err_    = 0;
        //! Weighting factor for dipole????
        double                    dip_weight_ = 0;
        //! GROMACS state variable
        t_state                  *state_      = nullptr;
        //! GROMACS force record
        t_forcerec               *fr_         = nullptr;
        //! Function that returns true if a molecule is symmetric
        bool             IsSymmetric(real toler);
        //! Vector to back up coordinates
        std::vector<gmx::RVec>    backupCoordinates_;
        //! Make a back up of coordinates
        void backupCoordinates();
        //! Restored backed up coordinates
        void restoreCoordinates();
        /*! \brief
         * Generate Atoms based on quantum calculation with specified level of theory
         *
         * \param[in]  pd     Force field data
         * \param[out] atoms  The structure to update
         * \param[in]  method Method used for QM calculations
         * \param[in]  basis  Basis set used for QM calculations
         * \param[out] mylot  Level of theory used
         */
        immStatus GenerateAtoms(const Poldata     *pd,
                                t_atoms           *atoms,
                                const std::string &method,
                                const std::string &basis,
                                std::string       *mylot);

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
        void MakeSpecialInteractions(const Poldata *pd,
                                     t_atoms       *atoms);
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
         *
         * \param[in] pd    The force field structure
         * \param[in] atoms The structure to check
         */
        immStatus checkAtoms(const Poldata *pd,
                             const t_atoms *atoms);

        /*! \brief
         * Return true if atom type neesd to have virtual site.
         *
         * \param[in] atype  Atom type
         * \param[in] pd     Data structure containing atomic properties
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

        /*! \brief Utility function to construct an identifier from GROMACS
         * topology indices.
         *
         * \param[in] pd     Poldata
         * \param[in] iType  The interaction type
         * \param[in] btype  Vector of all bond types for the atoms in this mymol
         * \param[in] natoms Number of atoms in this interaction
         * \param[in] iatoms The atom indices
         * \return An identifier
         */
        Identifier getIdentifier(const Poldata                  *pd,
                                 InteractionType                 iType,
                                 const std::vector<std::string> &btype,
                                 int                             natoms,
                                 const int                      *iatoms);

    public:
        double                         chieq_         = 0;
        // Enthalpy of formation (experimental) for this compound
        double                         Hform_         = 0;
        // Target molecule energy for this compound
        double                         Emol_          = 0;
        double                         anisoPol_elec_ = 0;
        double                         anisoPol_calc_ = 0;
        double                         mpad_          = 0; //molecular polarizability anisotropy difference (mpad)
        tensor                         alpha_elec_    = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        tensor                         alpha_calc_    = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        eSupport                       eSupp_         = eSupport::Local;
        PaddedVector<gmx::RVec>        f_;
        PaddedVector<gmx::RVec>        optf_;
        std::vector<int>               symmetric_charges_;
        QgenAcm                       *QgenAcm_        = nullptr;
        std::vector<PlistWrapper>      plist_;
        gmx_mtop_t                    *mtop_           = nullptr;
        gmx_localtop_t                *ltop_           = nullptr;
        gmx_shellfc_t                 *shellfc_        = nullptr;
        t_symtab                      *symtab_         = nullptr;
        t_inputrec                    *inputrec_       = nullptr;
        gmx_enerdata_t                *enerd_          = nullptr;
        t_fcdata                      *fcd_            = nullptr;
        t_nrnb                         nrnb_;
        gmx_wallcycle_t                wcycle_;
        std::vector<std::string>       error_messages_;

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

        /*! \brief
         * \return the center of nuclear charge
         */
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
         * Rotate the molcular dipole vector onto a reference vector
         *
         * \param[in] mu             Molecular dipole vector
         * \param[in] muReference    The reference vector
         */
        void rotateDipole(rvec mu, rvec muReference);

        /*! \brief
         * Return experimental dipole
         */
        double dipExper() const { return dip_exp_; }

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
         * \return mdatoms structure
         */
        t_mdatoms *getMdatoms() { return MDatoms_->get()->mdatoms(); }

        /*! \brief
         * \return atoms structure for editing
         */
        t_atoms *atoms();
        
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
         * Return mtop structure
         */
        const gmx_mtop_t *mtop() const { return mtop_; }
        /*! \brief
         * Return grompp Atomtype structure
         */
        const gpp_atomtype_t *gromppAtomtype() const { return &gromppAtomtype_; }
        /*! \brief
         * Return mol state
         */
        t_state *molState() const { return state_; }

        /*! \brief
         * It generates the atoms structure which will be used to print the topology file.
         *
         * \param[in] fp          File to write (debug) information to
         * \param[in] pd          Data structure containing atomic properties
         * \param[in] method      Method used for QM calculation
         * \param[in] basis       Basis set used for QM calculation
         * \param[out] mylot      Level of theory
         * \param[in] missing     How to treat missing parameters
         * \param[in] tabfn       Table function file for table potentials
         */
        immStatus GenerateTopology(FILE              *fp,
                                   const Poldata     *pd,
                                   const std::string &method,
                                   const std::string &basis,
                                   std::string       *mylot,
                                   missingParameters  missing,
                                   const char        *tabfn);

        /*! \brief
         *  Computes isotropic polarizability at the presence of external
         *  electric field (under construction!!!)
         *
         * \param[in]  efield   Strenght of the external electric field
         * \param[in]  fplog
         * \param[in]  cr
         * \returns the result of the calculation, if fine it is immOK
         */
        immStatus CalcPolarizability(double efield, t_commrec *cr, FILE *fplog);
      
        /*! \brief set the electronic polarizability
         * \param[in] isoPol The isotropic polarizability
         */
        void SetElectronicPolarizability(double isoPol) { isoPol_elec_ = isoPol; }
        /*! \brief get the electronic polarizability
         */
        double ElectronicPolarizability() { return isoPol_elec_; }
        /*! \brief get the calculated polarizability
         */
        double CalculatedPolarizability() { return isoPol_calc_; }
        /*! \brief get the electronic anisotropic polarizability
         */
        double ElectronicAnisoPolarizability() { return anisoPol_elec_; }
        /*! \brief get the calculated anisotropic polarizability
         */
        double CalculatedAnisoPolarizability() { return anisoPol_calc_; }
        /*! \brief Return difference between calculated and electronic isotropy polarizability 
         */
        double PolarizabilityDeviation() const { return isoPol_calc_ - isoPol_elec_; }
        /*! \brief Return difference between calculated and electronic anisotropy polarizability 
         */
        double AnisoPolarizabilityDeviation() const { return anisoPol_calc_ - anisoPol_elec_; }
        /*! \brief Return difference between calculated and electronic polarizability tensor
         */
        double PolarizabilityTensorDeviation() const;

        /*! \brief
         * Generate atomic partial charges
         *
         * \param[in] pd                             Data structure containing atomic properties
         * \param[in] fplog                          Logger
         * \param[in] cr      Communication parameters
         * \param[in] tabfn   Table function
         * \param[in] qcycle  Number of cycles for computing charges
         * \param[in] qtol    Convergence of charges tolerance
         * \param[in] algorithm The algorithm for determining charges,
         *                    if NONE it is read from the Poldata structure.
         * \param[in] qcustom Custom (user-provided) charges
         * \param[in] lot     Level of theory when using QM charges
         */
        immStatus GenerateCharges(const Poldata             *pd,
                                  const gmx::MDLogger       &fplog,
                                  t_commrec                 *cr,
                                  const char                *tabfn,
                                  int                        qcycle,
                                  real                       qtol,
                                  ChargeGenerationAlgorithm  algorithm,
                                  const std::vector<double> &qcustom,
                                  const std::string         &lot);
        /*! \brief
         * Generate atomic partial charges using EEM or SQE.
         * If shells are present they will be minimized.
         *
         * \param[in] pd      Data structure containing atomic properties
         * \param[in] cr      Communication parameters
         * \param[in] qcycle  Number of cycles for computing charges
         * \param[in] qtol    Convergence of charges tolerance
         */
        immStatus GenerateAcmCharges(const Poldata *pd,
                                     t_commrec     *cr,
                                     int            qcycle,
                                     real           qtol);
                                     
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
         * \param[in] cr      Gromacs communication record
         */
        void plotEspCorrelation(const char             *espcorr,
                                const gmx_output_env_t *oenv,
                                t_commrec              *cr);

        /*! \brief
         * Collect the experimental properties
         *
         * \param[in] iqm      Determine whether to allow exp or QM results or both
         * \param[in] bZero    Allow zero dipoles
         * \param[in] bZPE     Use zero point energies
         * \param[in] bDHform  Whether to use the enthalpy of formation
         * \param[in] method   Method used for QM calculation
         * \param[in] basis    Basis set used for QM calculation
         * \param[in] pd       Force field structure
         */
        immStatus getExpProps(iqmType            iqm,
                              gmx_bool           bZero,
                              gmx_bool           bZPE,
                              gmx_bool           bDHform,
                              const std::string &method,
                              const std::string &basis,
                              const Poldata     *pd);

        /*! \brief
         * Print the topology that was generated previously in GROMACS format.
         *
         * \param[in] fn        File name to open
         * \param[in] bVerbose  Verbose output
         * \param[in] pd        Data structure containing atomic properties
         * \param[in] cr        Gromacs communication record
         * \param[in] method    QC method
         * \param[in] basis     QC basis set
         */
        void PrintTopology(const char        *fn,
                           bool               bVerbose,
                           const Poldata     *pd,
                           t_commrec         *cr,
                           const std::string &method,
                           const std::string &basis);

        /*! \brief
         * Print the topology that was generated previously in GROMACS format.
         *
         * \param[in] fp        File pointer opened previously.
         * \param[in] bVerbose  Verbose
         * \param[in] pd        Data structure containing atomic properties
         * \param[in] bITP      Whether or not to write an itp file iso top file
         * \param[in] cr        Gromacs communication record
         * \param[in] method    QC method
         * \param[in] basis     WC basis set
         */
        void PrintTopology(FILE                    *fp,
                           bool                     bVerbose,
                           const Poldata           *pd,
                           bool                     bITP,
                           t_commrec               *cr,
                           const std::string       &method,
                           const std::string       &basis);

        /*! \brief
         *  Compute or derive global info about the molecule
         *
         * \param[in] pd   Data structure containing atomic properties
         */
        void CalcQPol(const Poldata *pd);

        /*! \brief
         * Compute the anisotropic polarizability from the polarizability tensor
         *
         * \param[in]  polar     Tensor of polarizability
         * \param[out] anisoPol  Anisotropic polarizability
         */
        void CalcAnisoPolarizability(tensor polar, double *anisoPol);


        /*! \brief
         * Relax the shells (if any) or compute the forces in the molecule
         *
         * \param[in]  fplog Log file pointer
         * \param[in]  cr    Communication record
         * \param[out] rmsf  Root mean square force on the shells
         * \return immOK if everything went fine, an error otherwise.
         */
        immStatus computeForces(FILE *fplog, t_commrec *cr, double *rmsf);

        /*! \brief
         * Change the coordinate of the molecule based
         * on the coordinate of the conformation stored
         * in molprop experiment class.
         *
         * \param[in] ei     Experiment
         * \param[in] bpolar Whether or not there are shells
         */
        void changeCoordinate(const Experiment &ei, gmx_bool bpolar);

        /*! \brief
         * Return the optimized geometry of the molecule from the data file.
         */
        bool getOptimizedGeometry(rvec *x);

        /*! \brief
         * Set the force field name.
         */
        void SetForceField(const char *ff) { forcefield_.assign(ff); }

        /*! \brief
         * Update internal structures for bondtype due to changes in pd
         *
         * \param[in] pd           Data structure containing atomic properties
         * \param[in] iType        Interaction type
         */
        void UpdateIdef(const Poldata   *pd,
                        InteractionType  iType);

        /*! \brief
         * Get the force field
         *
         * \param[out] forcefield_     Force field
         */
        std::string getForceField() { return forcefield_; }

        /*! \brief
         * Generate Charge Groups
         *
         * \param[in] ecg
         * \param[in] bUsePDBcharge
         */
        immStatus GenerateChargeGroups(eChargeGroup ecg, bool bUsePDBcharge);

        /*! \brief
         * Generate GROMACS structures.
         */
        immStatus GenerateGromacs(const gmx::MDLogger      &mdlog,
                                  t_commrec                *cr,
                                  const char               *tabfn,
                                  ChargeType               iType);

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
         * \param[in] fn A File pointer opened previously.
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
         * \param[in]  method  Method used for QM
         * \param[in]  basis   Basis set used for QM
         * \param[in]  watoms  Weight for the potential on the atoms in 
         *                     doing the RESP fit. Should be 0 in most cases.
         * \param[in]  maxESP  Percentage of the ESP points to consider (<= 100)
         */
        void initQgenResp(const Poldata     *pd,
                          const std::string &method,
                          const std::string &basis,
                          real               watoms,
                          int                maxESP);

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

}

#endif
