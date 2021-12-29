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

#ifndef MYMOL_LOW_H
#define MYMOL_LOW_H

#include <assert.h>

#include <cstdio>
#include <cstring>

#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/pbcutil/pbc.h"

#include "chargemodel.h"

struct gpp_atomtype;
struct t_atoms;
struct t_blocka;
struct t_excls;
struct t_symtab;

namespace alexandria
{

class Poldata;
class QgenAcm;
class Topology;

/*! \brief
 * Class to determine how to treat missing parameters.
 */
enum class missingParameters 
{ 
    //! Will generate an error if a parameter is missing
    Error,
    //! Missing parameters will be ignored
    Ignore, 
    //! Missing parameters will be generated
    Generate 
}; 

/*! \brief
 * Class with error codes from mymol.
 */
enum class immStatus {
    //! Unknown error, this should not happen
    Unknown,
    //! No error, all OK.
    OK,
    //! Zero dipole found in a molecule
    ZeroDip,
    //! No Quadrupole present
    NoQuad,
    //! The compound is charged
    Charged,
    //! Cannot determine the atom types for one or more atoms
    AtomTypes,
    //! Cannot determine the atom number for one or more atoms
    AtomNumber,
    //! Incorrect multiplicity
    Multiplicity,
    //! Problem converting data from the underlying molprop structure
    MolpropConv,
    //! Incorrect or unknown bondorder
    BondOrder,
    //! Initializing the Restrained Electrostatic Potential algorithm
    RespInit,
    //! Problem generating charges
    ChargeGeneration,
    //! Shell minimization did not converge
    ShellMinimization,
    //! Requested level of theory missing
    LOT,
    //! QM Inconsistency (ESP dipole does not match Electronic)
    QMInconsistency,
    //! Compound not in training set
    Test,
    //! No experimental data
    NoData,
    //! Problem generating shells
    GenShells,
    //! Problem generating bonds
    GenBonds,
    //! Problem communicating MolProp between processors
    CommProblem,
    //! Charge distribution width zeta is zero unexpectedly
    ZeroZeta,
    //! The number of data is lower than mindata
    InsufficientDATA,
    //! No dipole moment
    NoDipole,
    //! Not a supported bond
    NotSupportedBond,
    //! A not supported angle was found
    NotSupportedAngle,
    //! A not supported LinearAngle was found
    NotSupportedLinearAngle,
    //! A not supported Dihedral was found
    NotSupportedDihedral
};

const char *immsg(immStatus imm);
 
bool is_planar(const rvec   xi,  const rvec xj,
               const rvec   xk,  const rvec xl,
               const t_pbc *pbc, real phi_toler);

bool is_linear(const rvec xi, const rvec xj,
               const rvec xk, const t_pbc *pbc,
               real th_toler);

void add_excl(t_excls *excls, int e);

void add_excl_pair(t_excls excls[], int e1, int e2);

void remove_excl(t_excls *excls, int remove);

void let_shells_see_shells(t_excls excls[], t_atoms *atoms, struct gpp_atomtype *atype);

void copy_atoms(t_atoms *src, t_atoms *dest);

real calc_r13(const Poldata                  *pd,
              const std::vector<std::string> &atoms,
              const std::vector<double>      &bondOrders,
              const real                      angle);

real calc_relposition(const Poldata                  *pd,
                      const std::vector<std::string> &atoms,
                      const std::vector<double>      &bondOrders);

std::vector<double> getDoubles(const std::string &s);

void nonbondedFromPdToMtop(gmx_mtop_t    *mtop,
                           const Poldata *pd,
                           t_forcerec    *fr);

void excls_to_blocka(int natom, t_excls excls_[], t_blocka *blocka);

void mtop_update_cgs(gmx_mtop_t *mtop);

void put_in_box(int natom, matrix box, rvec x[], real dbox);

void write_zeta_q(FILE       *fp,
                  QgenAcm    *qgen,
                  t_atoms    *atoms,
                  ChargeType  iChargeType);

void write_zeta_q2(QgenAcm             *qgen,
                   struct gpp_atomtype *atype,
                   t_atoms             *atoms,
                   ChargeType           iChargeType);

void write_top(FILE                            *out,
               char                            *molname,
               t_atoms                         *at,
               gmx_bool                         bRTPresname,
               const Topology *topology,
               t_excls                          excls[],
               struct gpp_atomtype             *atype,
               int                             *cgnr,
               const Poldata                   *pd);

void calc_rotmatrix(rvec target_vec, rvec ref_vec, matrix rotmatrix);

gmx_mtop_t *do_init_mtop(const Poldata                   *pd,
                         char                           **molname,
                         t_atoms                         *atoms,
                         t_inputrec                      *ir,
                         t_symtab                        *symtab,
                         const char                      *tabfn);

} // namespace alexandria
#endif
