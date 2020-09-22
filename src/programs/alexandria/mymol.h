/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
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


struct gmx_enerdata_t;
struct gmx_shellfc_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct gmx_hw_info_t;

enum eDih {
    edihNo,
    edihOne,
    edihAll,
    edihNR
};

enum eSupport {
    eSupportNo,
    eSupportLocal,
    eSupportRemote,
    eSupportNR
};

struct gmx_vsite_t;
struct t_nrnb;

namespace alexandria
{
/*! \brief Enumerated type to differentiate the charge types */
enum qType {
    qtCalc      = 0,
    qtESP       = 1,
    qtMulliken  = 2,
    qtHirshfeld = 3,
    qtCM5       = 4,
    qtElec      = 5,
    qtNR        = 6
};

/*! \brief return string corresponding to charge type */
const char *qTypeName(qType qt);

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
        int                              nexcl_;
        int                             *cgnr_;
        bool                             bHaveShells_;
        bool                             bHaveVSites_;
        bool                             bNeedVsites_;
        double                           ref_enthalpy_;
        double                           polarizability_;
        double                           sig_pol_;
        double                           EspRms_;
        double                           EemRms_;
        t_excls                         *excls_;
        immStatus                        immAtoms_;
        immStatus                        immCharges_;
        immStatus                        immTopology_;
        std::unique_ptr<gmx_vsite_t>    *vsite_     = nullptr;
        std::unique_ptr<gmx::MDAtoms>   *MDatoms_   = nullptr;
        std::unique_ptr<gmx::MDModules> *mdModules_ = nullptr;
        MyForceProvider                 *myforce_;
        GentopVsites                     gvt_;
        std::string                      forcefield_;
        

        //! Array of dipole vectors
        rvec                      mu_qm_[qtNR];
        //! Array of quadrupole tensors
        tensor                    Q_qm_[qtNR];
        //! Array of vectors of charges
        std::vector<double>       charge_QM_[qtNR];
        //! Experimental dipole
        double                    dip_exp_    = 0;
        //! Error in experimental dipole
        double                    dip_err_    = 0;
        //! Weighting factor for dipole????
        double                    dip_weight_ = 0;
        //! Center of charge
        rvec                      coc_        = {0, 0, 0};
        //! GROMACS state variable
        t_state                  *state_;
        //! GROMACS force record
        t_forcerec               *fr_;
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
         * \param[in] ap      Gromacs atom properties
         * \param[in] method  Method used for QM calculations
         * \param[in] basis   Basis set used for QM calculations
         * \param[out] mylot  Level of theory used
         */
        immStatus GenerateAtoms(gmx_atomprop_t     ap,
                                const std::string &method,
                                const std::string &basis,
                                std::string       *mylot);

        /*! \brief
         * Generate angles, dihedrals, exclusions etc.
         *
         * \param[in] bPairs
         * \param[in] bDihs
         */
        void MakeAngles(bool bPairs,
                        bool bDihs);

        /*! \brief
         * Generate virtual sites or linear angles
         *
         * \param[in] bUseVsites
         */
        void MakeSpecialInteractions(const Poldata *pd,
                                     bool           bUseVsites);

        /*! \brief
         * Add shell particles
         *
         * \param[in] pd       Data structure containing atomic properties
         */
        void addShells(const Poldata *pd);

        /*! \brief
         * Check whether atom types exist in the force field
         *
         * \param[in] pd
         */
        immStatus checkAtoms(const Poldata *pd);

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

        

        /*! \brief Extract charges and electric moments and store them.
         *
         * \param[in] qt     Charge type to store in
         * \param[in] natom  Number of atoms
         * \param[in] q      The charges
         */
        void setQandMoments(qType qt, int natom, real q[]);
        
        

        double                         isoPol_elec_   = 0;
        double                         isoPol_calc_   = 0;
        bool                           gromacsGenerated_ = false;
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
        eSupport                       eSupp_;
        PaddedVector<gmx::RVec>        f_;
        PaddedVector<gmx::RVec>        optf_;
        std::vector<int>               symmetric_charges_;
        QgenResp                      *QgenResp_ = nullptr;
        QgenAcm                       *QgenAcm_ = nullptr;
        std::vector<PlistWrapper>      plist_;
        gmx_mtop_t                    *mtop_;
        gmx_localtop_t                *ltop_;
        gpp_atomtype_t                 atype_;
        gmx_shellfc_t                 *shellfc_;
        t_symtab                      *symtab_;
        t_inputrec                    *inputrec_;
        gmx_enerdata_t                *enerd_;
        t_atoms                       *atoms_;
        t_fcdata                      *fcd_;
        t_nrnb                         nrnb_;
        gmx_wallcycle_t                wcycle_;
        std::vector<std::string>       error_messages_;

        /*! \brief
         * Constructor
         */
        MyMol();

        /*! \brief
         * Return QM dipole corresponding to charge type qt
         */
        const rvec &muQM(qType qt) const { return mu_qm_[qt]; }

        rvec &muQM(qType qt) { return mu_qm_[qt]; }

        /*! \brief
         * Return QM quadrupole corresponding to charge type qt.
         */
        const tensor &QQM(qType qt) const { return Q_qm_[qt]; }

        /*! \brief
         * Return Charge vector corresponding to charge type qt.
         */
        const std::vector<double> &chargeQM(qType qt) const { return charge_QM_[qt]; }
        
        
        const std::vector<std::string> &errors() const {return error_messages_;}

        /*! \brief Store dipole in appropriate vector
         *
         * \param[in] qt The charge type
         * \param[in] mu The dipole to be stored
         */
        void set_muQM(qType qt, rvec mu) { copy_rvec(mu, mu_qm_[qt]); }

        /*! \brief
         * Rotate the molcular dipole vector onto a reference vector
         *
         * \param[in] mu             Molecular dipole vector
         * \param[in] muReference    The reference vector
         */
        void rotateDipole(rvec mu, rvec muReference);

        /*! \brief Store quadrupole in appropriate tensor
         *
         * \param[in] qt The charge type
         * \param[in] Q  The quadrupole to be stored
         */
        void set_QQM(qType qt, tensor Q) { copy_mat(Q, Q_qm_[qt]); }

        /*! \brief
         * Return computed dipole for charge type qt.
         */
        double dipQM(qType qt) const { return norm(mu_qm_[qt]); }

        /*! \brief
         * Return experimental dipole
         */
        double dipExper() const { return dip_exp_; }

        /*! \brief
         * Add the screening factors of the distributed charge to atom structure
         *
         * \param[in] pd      Data structure containing atomic properties
         */
        immStatus zeta2atoms(const Poldata *pd);

        /*! \brief
         * Return the coordinate vector of the molecule
         */
        const gmx::HostVector<gmx::RVec> &x() const { return state_->x; }

        /*! \brief
         * \return mdatoms structure
         */
        t_mdatoms *getMdatoms() { return MDatoms_->get()->mdatoms(); }

        /*! \brief
         * Return mtop structure
         */
        const gmx_mtop_t *mtop() { return mtop_; }
        /*! \brief
         * Return mol state
         */
        t_state *molState() const { return state_; }

        /*! \brief
         * It generates the atoms structure which will be used to print the topology file.
         *
         * \param[in] ap          Gromacs atom properties
         * \param[in] pd          Data structure containing atomic properties
         * \param[in] method      Method used for QM calculation
         * \param[in] basis       Basis set used for QM calculation
         * \param[in] nexcl       Number of Exclusions
         * \param[in] bUseVsites  Add virtual sites to the topology structure
         * \param[in] bPairs      Add pairs to the topology structure
         * \param[in] bDih        Add dihedrals to the topology structure
         */
        immStatus GenerateTopology(gmx_atomprop_t     ap,
                                   const Poldata     *pd,
                                   const std::string &method,
                                   const std::string &basis,
                                   std::string       *mylot,
                                   bool               bUseVsites,
                                   bool               bPairs,
                                   bool               bDih,
                                   bool               bBASTAT,
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
         * \param[in] watoms
         * \param[in] method                         Method used for QM calculation
         * \param[in] basis                          Basis set used for QM calculation
         * \param[in] cr
         * \param[in] tabfn
         */
        immStatus GenerateCharges(const Poldata          *pd,
                                  const gmx::MDLogger    &fplog,
                                  t_commrec              *cr,
                                  const char             *tabfn,
                                  gmx_hw_info_t          *hwinfo,
                                  int                     qcycle,
                                  real                    qtol);
                                  
        /*! \brief Implement charge symmetrization
         *
         * Initiates internal structure for atom charge symmetry
         * (e.g. CH3 with identical charges on H).
         * Must be called before initQgenResp.
         * \param[in] pd                 Data structure containing atomic properties
         * \param[in] ap                 Gromacs atom properties
         * \param[in] bSymmetricCharges  Consider molecular symmetry to calculate partial charge
         * \param[in] symm_string        The type of molecular symmetry
         */
        void symmetrizeCharges(const Poldata  *pd,
                               gmx_atomprop_t  ap,
                               bool            bSymmetricCharges,
                               const char     *symm_string);
        /*! \brief Init the class for the RESP algorithm.
         * This must be called before generateCharges even if no RESP
         * is used for charge generation, since ESP points can be used
         * in analysis.
         * \param[in]  pd                 Data structure containing atomic properties
         * \param[in]  method  Method used for QM
         * \param[in]  basis   Basis set used for QM
         * \param[out] mylot   Combined level of theory, may be nullptr
         * \param[in]  watoms  Weight for the potential on the atoms in 
         *                     doing the RESP fit. Should be 0 in most cases.
         * \param[in]  maxESP  Percentage of the ESP points to consider (<= 100)
         */
        void initQgenResp(const Poldata     *pd,
                          const std::string &method,
                          const std::string &basis,
                          std::string       *mylot,
                          real               watoms,
                          int                maxESP);

        /*! \brief Make a ESP correlation plot
         *
         * Will compute the electrostatic potential around the compound
         * and make a correlation plot between the QM potential and the
         * Alexandria potential.
         * \param[in] espcorr File name to plot to
         * \param[in] oenv    Gromacs output structure
         */
        void plotEspCorrelation(const char             *espcorr,
                                const gmx_output_env_t *oenv);
        /*! \brief
         * Return the root-mean square deviation of
         * the generated ESP from the QM ESP.
         *
         */
        double espRms() const { return EspRms_; }

        /*! \brief
         * Collect the experimental properties
         *
         * \param[in] bQM
         * \param[in] bZero
         * \param[in] method   Method used for QM calculation
         * \param[in] basis    Basis set used for QM calculation
         * \param[in] gap      Gaussian atom property
         */
        immStatus getExpProps(gmx_bool           bQM,
                              gmx_bool           bZero,
                              gmx_bool           bZPE,
                              gmx_bool           bDHform,
                              const std::string &method,
                              const std::string &basis,
                              const Poldata     *pd);

        /*! \brief
         * Print the topology that was generated previously in GROMACS format.
         *
         * \param[in] fn        A File pointer opened previously.
         * \param[in] iModel    The distrbution model of charge (e.x. point charge, gaussian, and slater models)
         * \param[in] bVerbose  Verobse
         * \param[in] pd        Data structure containing atomic properties
         * \param[in] aps       Gromacs atom properties
         */
        void PrintTopology(const char        *fn,
                           bool               bVerbose,
                           const Poldata     *pd,
                           gmx_atomprop_t     aps,
                           t_commrec         *cr,
                           double             efield,
                           const std::string &method,
                           const std::string &basis);

        /*! \brief
         * Print the topology that was generated previously in GROMACS format.
         *
         * \param[in] fn        A File pointer opened previously.
         * \param[in] iModel    The distrbution model of charge (e.x. point charge, gaussian, and slater models)
         * \param[in] bVerbose  Verbose
         * \param[in] pd        Data structure containing atomic properties
         * \param[in] aps       Gromacs atom properties
         */
        void PrintTopology(FILE                    *fp,
                           bool                     bVerbose,
                           const Poldata           *pd,
                           gmx_atomprop_t           aps,
                           bool                     bITP,
                           t_commrec               *cr,
                           double                   efield,
                           const std::string       &method,
                           const std::string       &basis);

        /*! \brief
         *  Compute or derive global info about the molecule
         *
         * \param[in] pd   Data structure containing atomic properties
         * \param[out] mu  The dipole
         */
        void CalcQPol(const Poldata *pd, rvec mu);

        /*! \brief
         * Compute the molecular dipole vector
         */
        void CalcDipole();

        /*! \brief
         * Compute the molecular dipole vector and store them in the vector mu
         */
        void CalcDipole(rvec mu);

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
         * Set the force field
         *
         * \param[in] ff   Force field
         */

        /*! \brief
         * Change the coordinate of the molecule based
         * on the coordinate of the conformation stored
         * in molprop experiment class.
         *
         * \param[in] ei   ExperimentIterator
         */
        void changeCoordinate(ExperimentIterator ei, gmx_bool bpolar);

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
         * \param[in] modifiedOnly Only update parameters that were modified.
         *                         Experimental speed-up for tune_fc.
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
         * Calculate quadrupole tensor
         */
        void CalcQuadrupole();

        /*! \brief Calculates dipole components, and quadrupoles.
         *
         * Compute moments using QM-based charges like
         * Mulliken, Hirshfeld, CM5, etc. Since there is no Shell particle in
         * QM calculations, it loops over eptAtoms, only.
         * \param[in]  q  Array of charges
         * \param[out] mu Dipole vector
         * \param[out] Q  Quadrupole tensor
         */
        void CalcQMbasedMoments(real *q, rvec mu, tensor Q);

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
                                  gmx_hw_info_t            *hwinfo,
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
};

const char *immsg(immStatus imm);

}

#endif
