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
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/pbcutil/pbc.h"

#include "gentop_core.h"
#include "gentop_vsite.h"
#include "poldata.h"
#include "qgen_acm.h"
#include "qgen_resp.h"

namespace alexandria
{

enum class missingParameters { Error, Ignore, Generate }; 

enum class immStatus {
    Unknown,
    OK,
    ZeroDip,
    NoQuad,
    Charged,
    AtomTypes,
    AtomNumber,
    MolpropConv,
    BondOrder,
    RespInit,
    ChargeGeneration,
    ShellMinimization,
    LOT,
    QMInconsistency,
    Test,
    NoData,
    GenShells,
    GenBonds,
    CommProblem,
    ZeroZeta,
    InsufficientDATA,
    NoDipole,
    NotSupportedBond,
    NotSupportedAngle,
    NotSupportedDihedral
};

const char *immsg(immStatus imm);
 
bool is_planar(rvec   xi,  rvec xj,
               rvec   xk,  rvec xl,
               t_pbc *pbc, real phi_toler);

bool is_linear(rvec xi, rvec xj,
               rvec xk, t_pbc *pbc,
               real th_toler);

void add_excl(t_excls *excls, int e);

void add_excl_pair(t_excls excls[], int e1, int e2);

void remove_excl(t_excls *excls, int remove);

void let_shells_see_shells(t_excls excls[], t_atoms *atoms, gpp_atomtype_t atype);

void copy_atoms(t_atoms *src, t_atoms *dest);

void cp_plist(t_params                   plist[],
              int                        ftype,
              InteractionType            itype,
              std::vector<PlistWrapper> &plist_);

real calc_r13(const Poldata     *pd,
              const std::string &aai,
              const std::string &aaj,
              const std::string &aak,
              const real         angle);

real calc_relposition(const Poldata     *pd,
                      const std::string  aai,
                      const std::string  aaj,
                      const std::string  aak);

/*! \brief Update GROMACS force field parameters
 *
 * The force field parameters in gromacs idef-style structures are updated from the
 * Poldata structure. This can be a heavy calculation for a large system, but it will 
 * only handle the  bonds etc. in the system (a molecule typically).
 * \param[in]  pd      Poldata structure
 * \param[out] plist   The parameter vector
 * \param[in]  atoms   GROMACS atom structure
 * \param[in]  bAllowMissingParameters Boolean determining whether parameters can be missing
 * \param[in]  molname Molecule name
 * \param[out] errors  Error messages are appended to this vector.
 * \return Warning status
 */ 
immStatus updatePlist(const Poldata             *pd,
                      std::vector<PlistWrapper> *plist,
                      const t_atoms             *atoms,
                      missingParameters          missing,
                      const std::string         &molname,
                      std::vector<std::string>  *errors);

std::vector<double> getDoubles(const std::string &s);

void nonbondedFromPdToMtop(gmx_mtop_t    *mtop,
                           const Poldata *pd,
                           t_forcerec    *fr);

void plist_to_mtop(const std::vector<PlistWrapper> &plist,
                   gmx_mtop_t                      *mtop_);

gmx_mtop_t *do_init_mtop(const Poldata                   *pd,
                         char                           **molname,
                         t_atoms                         *atoms,
                         const std::vector<PlistWrapper> &plist,
                         t_inputrec                      *ir,
                         t_symtab                        *symtab,
                         const char                      *tabfn);

void excls_to_blocka(int natom, t_excls excls_[], t_blocka *blocka);

void mtop_update_cgs(gmx_mtop_t *mtop);

void put_in_box(int natom, matrix box, rvec x[], real dbox);

void write_zeta_q(FILE       *fp,
                  QgenAcm    *qgen,
                  t_atoms    *atoms,
                  ChargeType  iChargeType);

void write_zeta_q2(QgenAcm        *qgen,
                   gpp_atomtype_t  atype,
                   t_atoms        *atoms,
                   ChargeType      iChargeType);

int get_subtype(directive d, int ftype);

void print_bondeds(FILE                            *out,
                   directive                        d,
                   int                              plist_ftype,
                   int                              print_ftype,
                   const std::vector<PlistWrapper> &plist);

void write_top(FILE                            *out,
               char                            *molname,
               t_atoms                         *at,
               gmx_bool                         bRTPresname,
               const std::vector<PlistWrapper> &plist,
               t_excls                          excls[],
               gpp_atomtype_t                   atype,
               int                             *cgnr,
               int                              nrexcl,
               const Poldata                   *pd);

void print_top_header(FILE                    *fp,
                      const Poldata           *pd,
                      bool                     bPol,
                      std::vector<std::string> commercials,
                      bool                     bItp);

void calc_rotmatrix(rvec target_vec, rvec ref_vec, matrix rotmatrix);

} // namespace alexandria
#endif
