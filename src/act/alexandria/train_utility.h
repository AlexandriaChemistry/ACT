/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#ifndef TRAIN_UTILITY_H
#define TRAIN_UTILITY_H

#include <cstdio>
#include <vector>

#include "act/alexandria/actmol.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/forces/forcecomputer.h"
#include "act/forcefield/forcefield.h"
#include "act/statistics/statistics.h"
#include "act/utility/jsontree.h"
#include "gromacs/commandline/pargs.h"

/*! \brief Utility function to merge command line arguments
 * \param[inout] pargs The complete list of arguments
 * \param[in]    npa   The nmber of new elements
 * \param[in]    pa    The new elements
 */
void doAddOptions(std::vector<t_pargs> *pargs, size_t npa, t_pargs pa[]);

namespace gmx
{
class TextWriter;
}

namespace alexandria
{

class MsgHandler;

//! \brief Short-cut for property-dependent statistics
using qtStats = std::map<qPropertyType, gmx_stats>;


/*! \brief Class to manage output from force field training
 */
class TrainForceFieldPrinter
{
private:
    //! Can make computations for molecules, such as the Hessian and energy minimization (for now).
    MolHandler molHandler_;
    //! Tolerance (kJ/mol e) for marking ESP as an outlier in the log file
    real esp_toler_           = 30;
    //! Tolerance (Debye) for marking dipole as an outlier in the log file
    real dip_toler_           = 0.5;
    //! Tolerance (Buckingham) for marking quadrupole as an outlier in the log file
    real quad_toler_          = 5;
    //! Tolerance (Debye A^2) for marking octupole as an outlier in the log file
    real oct_toler_          = 5;
    //! Tolerance (Debye A^3) for marking hexadecapole as an outlier in the log file
    real hex_toler_          = 5;
    //! Tolerance (A^3) for marking diagonal elements of the polarizability tensor as an outlier in the log file
    real alpha_toler_         = 3;
    //! Tolerance (A^3) for marking isotropic polarizability as an outlier in the log file
    real isopol_toler_        = 2;
    //! Fit regression analysis of results to y = ax+b instead of y = ax
    bool useOffset_           = false;
    //! Print all information from all SP calculation
    bool printSP_             = true;
    //! Perform energy minimization and compute vibrational frequencies for each molecule (after optimizing the force field if -optimize is enabled)
    bool calcFrequencies_     = false;
    //! Do special analysis of diatomic compounds
    bool diatomic_            = false;
    //! Dump outliers to xyz files if larger or equal to zero
    real dumpOutliers_        = -1;
    //! Data structures for storing energies per molecule
    std::map<std::string, std::vector<ACTEnergy> > molEnergyMap_;
    //! Data structures for storing energy terms per dimer
    std::map<std::string, ACTEnergyMapVector>      molInteractionEnergyMap_;

    //! Data structures for storing results of interaction energies
    std::map<InteractionType, std::map<iMolSelect, qtStats>> lsq_einter_;
    //! Comparison of epot, isoPol, anisoPol and alpha
    std::map<iMolSelect, qtStats>      lsq_epot_, lsq_isoPol_, lsq_anisoPol_, lsq_alpha_;
    //! Statistics for RMSF and frequencies.
    std::map<iMolSelect, gmx_stats>    lsq_rmsf_;
    gmx_stats                          lsq_freq_;
    //! Command line stuff
    std::map<MolPropObservable, std::string> mpoOpt_;
    std::map<MolPropObservable, std::string> mpoFnm_;
    //! Statistics for multipoles
    std::map<MolPropObservable, std::map<iMolSelect, qtStats> > lsq_multi_;
    //! Interaction energy terms
    std::vector<InteractionType> terms_;

    //! \brief Analyse polarizability, add to statistics and print
    void analysePolarisability(gmx::TextWriter     *tw,
                               const ForceField    *pd,
                               alexandria::ACTMol  *mol,
                               iMolSelect           ims,
                               const ForceComputer *forceComp);
    
    //! \brief Analyses dipoles, quadrupoles, etc.
    void analyse_multipoles(MsgHandler                                      *msg_handler,
                            const std::vector<alexandria::ACTMol>::iterator &mol,
                            std::map<MolPropObservable, double>              toler,
                            const ForceField                                *pd,
                            const ForceComputer                             *forceComputer);
    //! \brief And the atoms.
    void printAtoms(gmx::TextWriter              *tw,
                    alexandria::ACTMol           *mol,
                    const std::vector<gmx::RVec> &coords,
                    const std::vector<gmx::RVec> &forces);
    
    /*! \brief Write ACT energies to a molprop file
     * \param[in]    msghandler For warnings etc.
     * \param[in]    mpout      File name to write to
     * \param[in]    pd         A force field
     * \param[in]    forceComp  The force computer
     * \param[inout] mols       Vector of molecules, will be modified
     */
    void writeMolpropsEnergies(MsgHandler          *msghandler,
                               const char          *mpout,
                               const ForceField    *pd,
                               const ForceComputer *forceComp,
                               std::vector<ACTMol> *mols);
    /*! \brief do part of the printing, add to statistics
     */
    void printEnergyForces(MsgHandler                          *msghandler,
                           const ForceField                    *pd,
                           const ForceComputer                 *forceComp,
                           const std::map<eRMS, FittingTarget> &targets,
                           const AtomizationEnergy             &atomenergy,
                           alexandria::ACTMol                  *mol,
                           iMolSelect                          ims,
                           const gmx_output_env_t              *oenv,
                           bool                                 printAll);
    /*! \brief Print data on outliers.
     */
    void printOutliers(gmx::TextWriter                       *tw,
                       iMolSelect                             ims,
                       double                                 sigma,
                       bool                                   bIntermolecular,
                       InteractionType                        itype,
                       const std::vector<alexandria::ACTMol> *actmol);
public:
    TrainForceFieldPrinter();
    
    /*! \brief Add my options to the list of command line arguments
     * \param[out] pargs The vector to add to
     */
    void addOptions(std::vector<t_pargs> *pargs);
    
    /*! \brief Add my files to the list of command line arguments
     * \param[out] filenm The vector to add to
     */
    void addFileOptions(std::vector<t_filenm> *filenm);

    /*! \brief Do the force field info printing to the log file.
     * By default only the interactions/energies/observables that are
     * used in optimization will be printed.
     * \param[in] msghandler For printing and status
     * \param[in] sii        Information on the training
     * \param[in] actmol     The compounds/dimers
     * \param[in] oenv       For printing xvg files
     * \param[in] filenm     Filenames for additional output files
     * \param[in] printAll   Tell the printer to print all observables.
     */
    void print(MsgHandler                  *msghandler,
               StaticIndividualInfo        *sii,
               std::vector<ACTMol>         *actmol,
               const gmx_output_env_t      *oenv,
               const std::vector<t_filenm> &filenm,
               bool                         printAll);
};

/*! \brief Print header and command line arguments
 *
 * \param[in] tw      TextWriter pointer, if nullptr
 *                    will not do anything
 * \param[in] pargs   The command line arguments
 * \param[in] filenms The filenms used
 */
void print_header(gmx::TextWriter             *tw,
                  const std::vector<t_pargs>  &pargs,
                  const std::vector<t_filenm> &filenms);
                  
/*! \brief Do an analysis of frequencies compared to reference if present.
 * \param[in]    pd           The force field
 * \param[in]    mol          Molecule data
 * \param[in]    molhandler   The molecule handler
 * \param[in]    forceComp    The force computer
 * \param[in]    coords       Minimized coordinates
 * \param[in]    atomenergy   The atomization energy data
 * \param[out]   lsq_freq_all Statistics structure, may be nullptr
 * \param[inout] jtree        Output structure
 * \param[in]    spectrumFileName If not nullptr, a simulated IR spectrum will be written to this file
 * \param[in]    lineWidth    The Lorentzian line width for printing a spectrum
 * \param[in]    oenv         Structure to print xvg files
 * \param[in]    debugNMA     Will provide excessive printing statements
 */
void doFrequencyAnalysis(const ForceField         *pd,
                         alexandria::ACTMol       *mol,
                         const MolHandler         &molhandler,
                         const ForceComputer      *forceComp,
                         std::vector<gmx::RVec>   *coords,
                         const AtomizationEnergy  &atomenergy,
                         gmx_stats                *lsq_freq_all,
                         JsonTree                 *jtree,
                         const char               *spectrumFileName,
                         double                    lineWidth,
                         gmx_output_env_t         *oenv,
                         bool                      debugNMA);
                            
} // namespace alexandria

#endif
