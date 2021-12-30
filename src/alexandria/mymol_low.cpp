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


#include "mymol_low.h"

#include <assert.h>

#include <cstdio>
#include <cstring>

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"

#include "forcefieldparameter.h"
#include "poldata.h"
#include "qgen_acm.h"
#include "topology.h"
#include "units.h"

namespace alexandria
{

std::map<immStatus, const char *> immMessages = {
    { immStatus::Unknown,                  "Unknown status" },
    { immStatus::OK,                       "OK" },
    { immStatus::ZeroDip,                  "Zero Dipole" },
    { immStatus::NoQuad,                   "No Quadrupole" },
    { immStatus::Charged,                  "Charged" },
    { immStatus::AtomTypes,                "Atom type problem" },
    { immStatus::AtomNumber,               "Atom number problem" },
    { immStatus::MolpropConv,              "Converting from molprop" },
    { immStatus::Multiplicity,             "Number of electrons does not match the multiplicity. Is the total charge correct?" },
    { immStatus::BondOrder,                "Determining bond order" },
    { immStatus::RespInit,                 "RESP Initialization" },
    { immStatus::ChargeGeneration,         "Charge generation" },
    { immStatus::ShellMinimization,        "Shell minimization" },
    { immStatus::LOT,                      "Requested level of theory missing" },
    { immStatus::QMInconsistency,          "QM Inconsistency (ESP dipole does not match Electronic)" },
    { immStatus::Test,                     "Compound not in training set" },
    { immStatus::NoData,                   "No experimental data" },
    { immStatus::GenShells,                "Generating shells" },
    { immStatus::GenBonds,                 "Generating bonds" },
    { immStatus::CommProblem,              "Communicating MolProp" },
    { immStatus::ZeroZeta,                 "Charge distribution width zeta is zero unexpectedly" },
    { immStatus::InsufficientDATA,         "The number of data is lower than mindata" },
    { immStatus::NoDipole,                 "No Dipole moment" },
    { immStatus::NotSupportedBond,         "NotSupportedBond" },
    { immStatus::NotSupportedAngle,        "NotSupportedAngle" },
    { immStatus::NotSupportedLinearAngle,  "NotSupportedLinearAngle" },
    { immStatus::NotSupportedDihedral,     "NotSupportedDihedral" }
};

const char *immsg(immStatus imm)
{
    return immMessages[imm];
}

bool is_planar(const rvec xi,  const rvec xj, const rvec xk,
               const rvec xl, const t_pbc *pbc,
               real phi_toler)
{
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real phi;

    phi = RAD2DEG*dih_angle(xi, xj, xk, xl, pbc, r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);

    return (fabs(phi) < phi_toler);
}

bool is_linear(const rvec xi, const rvec xj,
               const rvec xk, const t_pbc *pbc,
               real th_toler)
{
    int  t1, t2;
    rvec r_ij, r_kj;
    real costh, th;

    th = fabs(RAD2DEG*bond_angle(xi, xj, xk, pbc, r_ij, r_kj, &costh, &t1, &t2));
    if ((th > th_toler) || (th < 180-th_toler))
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Angle is %g, th_toler is %g\n", th, th_toler);
        }
        return true;
    }
    return false;
}

void copy_atoms(t_atoms *src, t_atoms *dest)
{
    int i;

    if (dest->nr < src->nr)
    {
        srenew(dest->atom, src->nr);
        srenew(dest->atomname, src->nr);
        if (nullptr != src->atomtype)
        {
            srenew(dest->atomtype, src->nr);
        }
        else if (nullptr != dest->atomtype)
        {
            sfree(dest->atomtype);
            dest->atomtype = nullptr;
        }
        if (nullptr != src->atomtypeB)
        {
            srenew(dest->atomtypeB, src->nr);
        }
        else if (nullptr != dest->atomtypeB)
        {
            sfree(dest->atomtypeB);
            dest->atomtypeB = nullptr;
        }
    }
    dest->nr = src->nr;
    for (i = 0; (i < src->nr); i++)
    {
        dest->atom[i]      = src->atom[i];
        dest->atomname[i]  = src->atomname[i];
        if (nullptr != src->atomtype)
        {
            dest->atomtype[i]  = src->atomtype[i];
        }
        if (nullptr != src->atomtypeB)
        {
            dest->atomtypeB[i] = src->atomtypeB[i];
        }
    }
    if (dest->nres < src->nres)
    {
        srenew(dest->resinfo, src->nres);
    }

    if (nullptr != src->pdbinfo)
    {
        srenew(dest->pdbinfo, src->nres);
    }
    else if (nullptr != dest->pdbinfo)
    {
        sfree(dest->pdbinfo);
        dest->pdbinfo = nullptr;
    }
    dest->nres = src->nres;
    for (i = 0; (i < src->nres); i++)
    {
        dest->resinfo[i] = src->resinfo[i];
        if (nullptr != src->pdbinfo)
        {
            dest->pdbinfo[i] = src->pdbinfo[i];
        }
    }
}

real calc_r13(const Poldata                  *pd,
              const std::vector<std::string> &atoms,
              const std::vector<double>      &bondOrders,
              const real                      angle)
{
    double r12 = 0, r23 = 0, r13 = 0;

    Identifier aij ({atoms[0], atoms[1] }, { bondOrders[0] }, CanSwap::Yes);
    Identifier akj ({atoms[2], atoms[1] }, { bondOrders[1] }, CanSwap::Yes);

    std::string type("bondlength");
    auto bij = pd->findForcesConst(InteractionType::BONDS).findParameterTypeConst(aij, type);
    auto bkj = pd->findForcesConst(InteractionType::BONDS).findParameterTypeConst(akj, type);

    r12 = convertToGromacs(bij.value(), bij.unit());
    r23 = convertToGromacs(bkj.value(), bkj.unit());

    r13 = std::sqrt((r12*r12) + (r23*r23) - (2*r12*r23*std::cos(DEG2RAD*angle)));

    return r13;
}

real calc_relposition(const Poldata                  *pd,
                      const std::vector<std::string> &atoms,
                      const std::vector<double>      &bondOrders)
{
    double b0 = 0, b1 = 0, relative_position = 0;

    Identifier aij({atoms[0], atoms[1]}, { bondOrders[0] }, CanSwap::Yes);
    Identifier ajk({atoms[1], atoms[2]}, { bondOrders[1] }, CanSwap::Yes);

    std::string type("bondlength");
    auto bij = pd->findForcesConst(InteractionType::BONDS).findParameterTypeConst(aij, type);
    auto bjk = pd->findForcesConst(InteractionType::BONDS).findParameterTypeConst(ajk, type);

    b0 = convertToGromacs(bij.value(), bij.unit());
    b1 = convertToGromacs(bjk.value(), bjk.unit());

    relative_position = (b1/(b0+b1));

    return relative_position;
}


std::vector<double> getDoubles(const std::string &s)
{
    std::vector<double> d;

    for (auto &ss : gmx::splitString(s))
    {
        d.push_back(gmx::doubleFromString(ss.c_str()));
    }
    return d;
}

static int getCombinationRule(const ForceFieldParameterList &vdw)
{
    auto combRule = vdw.optionValue("combination_rule");
    int  i;
    for(i = 0; i < eCOMB_NR; i++)
    {
        if (combRule.compare(ecomb_names[i]) == 0)
        {
            break;
        }
    }
    GMX_RELEASE_ASSERT(i < eCOMB_NR, gmx::formatString("Cannot find combination rule %s in GROMACS",
                                                       combRule.c_str()).c_str());
    return i;
}

static void getSigmaEpsilon(const ForceFieldParameterList &fa,
                            const std::string             &ai,
                            double                        *sigma,
                            double                        *epsilon)
{
    Identifier idI(ai);
    *sigma   = 0;
    *epsilon = 0;
    if (fa.parameterExists(idI))
    {
        *sigma   = fa.findParameterTypeConst(idI, "sigma").value();
        *epsilon = fa.findParameterTypeConst(idI, "epsilon").value();
    }
}

static void getLjParams(const ForceFieldParameterList &fa,
                        const std::string             &ai,
                        const std::string             &aj,
                        double                        *c6,
                        double                        *cn)
{
    double sigmaI, sigmaJ, epsilonI, epsilonJ;
    getSigmaEpsilon(fa, ai, &sigmaI, &epsilonI);
    getSigmaEpsilon(fa, aj, &sigmaJ, &epsilonJ);

    auto comb     = getCombinationRule(fa);     
    switch (comb)
    {
    case eCOMB_GEOMETRIC:
        {
            *c6 = std::sqrt(sigmaI * sigmaJ);
            *cn = std::sqrt(epsilonI * epsilonJ);
        }
        break;
    case eCOMB_ARITHMETIC:
        {
            double sig  = 0.5 * (sigmaI + sigmaJ);
            double eps  = std::sqrt(epsilonI + epsilonJ);
            double sig6 = std::pow(sig, 6.0);
            *c6 = 4*eps*sig6;
            *cn = *c6 * sig6;
        }
        break;
    case eCOMB_GEOM_SIG_EPS:
        {
            double sig  = std::sqrt(sigmaI * sigmaJ);
            double eps  = std::sqrt(epsilonI * epsilonJ);
            double sig6 = std::pow(sig, 6.0);
            *c6 = 4*eps*sig6;
            *cn = *c6 * sig6;
        }
        break;
    case eCOMB_NONE:
    case eCOMB_NR:
        gmx_fatal(FARGS, "Unsupported combination rule %d for Lennard Jones", comb);
    }
}

static void getSigmaEpsilonGamma(const ForceFieldParameterList &fa,
                                 const std::string             &ai,
                                 double                        *sigma,
                                 double                        *epsilon,
                                 double                        *gamma)
{
    Identifier idI(ai);
    *sigma   = 0;
    *epsilon = 0;
    *gamma   = 0;
    if (fa.parameterExists(idI))
    {
        *sigma   = fa.findParameterTypeConst(idI, "sigma").value();
        *epsilon = fa.findParameterTypeConst(idI, "epsilon").value();
        *gamma   = fa.findParameterTypeConst(idI, "gamma").value();
    }
}

static void getBhamParams(const ForceFieldParameterList &fa,
                          const std::string             &ai,
                          const std::string             &aj,
                          double                        *a,
                          double                        *b,
                          double                        *c)
{
    double sigmaI, sigmaJ, epsilonI, epsilonJ, gammaI, gammaJ;
    getSigmaEpsilonGamma(fa, ai, &sigmaI, &epsilonI, &gammaI);
    getSigmaEpsilonGamma(fa, aj, &sigmaJ, &epsilonJ, &gammaJ);

    auto comb     = getCombinationRule(fa);
    switch (comb)
    {
        case eCOMB_GEOMETRIC:
            *a = std::sqrt(sigmaI * sigmaJ);
            *b = std::sqrt(epsilonI * epsilonJ);
            *c = std::sqrt(gammaI * gammaJ);
            break;
        case eCOMB_ARITHMETIC:
            *a = 0.5 * (sigmaI + sigmaJ);
            *b = std::sqrt(epsilonI * epsilonJ);
            *c = 0.5 * (gammaI + gammaJ);
            break;
        case eCOMB_KONG_MASON:
            *a = std::sqrt(sigmaI * sigmaJ);
            *b = 2.0 * (epsilonI * epsilonJ)/(epsilonI + epsilonJ);
            *c = *a * (0.5*((gammaI/sigmaI)+(gammaJ/sigmaJ)));
            break;
        case eCOMB_GEOM_SIG_EPS:
        case eCOMB_NONE:
        case eCOMB_NR:
            gmx_fatal(FARGS, "Unsupported combination rule %d for Buckingham", comb);
    }
    if (debug)
    {
        fprintf(debug, "BHAM parameters for %s %s %g %g %g\n",
                ai.c_str(), aj.c_str(), *a, *b, *c);
    }
}

void nonbondedFromPdToMtop(gmx_mtop_t    *mtop,
                           const Poldata *pd,
                           t_forcerec    *fr)
{
    auto ntype  = mtop->ffparams.atnr;
    int  ntype2 = gmx::square(ntype);
    if (mtop->ffparams.functype.empty())
    {
        mtop->ffparams.functype.resize(ntype2, 0);
    }
    if (mtop->ffparams.iparams.empty())
    {
        mtop->ffparams.iparams.resize(ntype2, {});
    }
    auto forcesVdw = pd->findForcesConst(InteractionType::VDW);
    auto ftypeVdW  = forcesVdw.fType();
    typedef struct
    {
        std::string name;
        int         ptype;
    } mytype;
    std::vector<mytype> mytypes;
    mytypes.resize(ntype);
    t_atoms *atoms = &mtop->moltype[0].atoms;
    for (int i = 0; i < atoms->nr; i++)
    {
        auto itype = atoms->atom[i].type;
        GMX_ASSERT(itype >= 0 && itype < ntype, 
                   gmx::formatString("Atom type %d for %s out of range [0, %d]",
                                     itype, *atoms->atomname[i], ntype).c_str());
        auto iname = *atoms->atomtype[i];
        if (mytypes[itype].name.empty())
        {
            mytypes[itype].name.assign(iname);
            mytypes[itype].ptype = atoms->atom[i].ptype;
        }
        else if (mytypes[itype].name.compare(iname) != 0)
        {
            gmx_fatal(FARGS, "Atom type mismatch. Found %s for type %d but expected %s", iname, itype, mytypes[itype].name.c_str());
        }
    }
    if (debug)
    {
        int i = 0;
        for (auto &tn : mytypes)
        {
            fprintf(debug, "type %d name %s ptype %s\n",
                    i++, tn.name.c_str(), ptype_str[tn.ptype]);
        }
    }
    // TODO: Use the symmetry in the matrix
    auto fa = pd->findForcesConst(InteractionType::VDW);
    for (auto i = 0; (i < ntype); i++)
    {
        if (mytypes[i].ptype == eptAtom)
        {
            for (auto j = 0; (j < ntype); j++)
            {
                auto idx = ntype*i+j;
                if (mytypes[j].ptype == eptAtom)
                {
                    mtop->ffparams.functype[idx] = ftypeVdW;
                    switch (ftypeVdW)
                    {
                        case F_LJ:
                        {
                            double c6 = 0, c12 = 0;
                            getLjParams(fa, mytypes[i].name,
                                        mytypes[j].name, &c6, &c12);
                            mtop->ffparams.iparams[idx].lj.c6  = c6;
                            mtop->ffparams.iparams[idx].lj.c12 = c12;
                            if (nullptr != fr)
                            {
                                C6(fr->nbfp, ntype, i, j)  = c6*6.0;
                                C12(fr->nbfp, ntype, i, j) = c12*12.0;
                            }
                        }
                        break;
                        case F_BHAM:
                        {
                            double a = 0, b = 0, c = 0;
                            getBhamParams(fa, mytypes[i].name,
                                          mytypes[j].name, &a, &b, &c);
                            mtop->ffparams.iparams[idx].bham.a = a;
                            mtop->ffparams.iparams[idx].bham.b = b;
                            mtop->ffparams.iparams[idx].bham.c = c;
                            if (nullptr != fr)
                            {
                                BHAMA(fr->nbfp, ntype, i, j) = a;
                                BHAMB(fr->nbfp, ntype, i, j) = b;
                                /* nbfp now includes the 6.0 derivative prefactor */
                                /* the 6.0 derivative prefactor is turned off for the modified BHAM implemented in nb_generic */
                                /* BHAMC(nbfp, atnr, i, j) = idef->iparams[k].bham.c*6.0; */
                                BHAMC(fr->nbfp, ntype, i, j) = c;
                            }
                            if (debug)
                            {
                                fprintf(debug, "idx = %3d a = %10g b = %10g c = %10g\n",
                                        idx,
                                        mtop->ffparams.iparams[idx].bham.a,
                                        mtop->ffparams.iparams[idx].bham.b,
                                        mtop->ffparams.iparams[idx].bham.c);
                            }
                        }
                        break;
                        default:
                            fprintf(stderr, "Invalid van der waals type %s\n",
                                    interaction_function[ftypeVdW].longname);
                    }
                }
            }
        }
    }
}


gmx_mtop_t *do_init_mtop(const Poldata                   *pd,
                         char                           **molname,
                         t_atoms                         *atoms,
                         t_inputrec                      *ir,
                         t_symtab                        *symtab,
                         const char                      *tabfn)
{
    gmx_mtop_t *mtop = new gmx_mtop_t();
    mtop->name     = molname;
    mtop->moltype.resize(1);
    mtop->moltype[0].name = molname;
    mtop->molblock.resize(1);
    mtop->molblock[0].nmol        = 1;
    mtop->molblock[0].type        = 0;
    mtop->groups.grps[egcENER].nr = 1;
    mtop->natoms                  = atoms->nr;
    init_t_atoms(&(mtop->moltype[0].atoms), atoms->nr, false);
    snew(mtop->moltype[0].atoms.atomtype, atoms->nr);
    snew(mtop->moltype[0].atoms.atomname, atoms->nr);
    snew(mtop->moltype[0].atoms.resinfo, atoms->nres+1);
    mtop->moltype[0].atoms.nres = atoms->nres;
    /* Count the number of atom types in the molecule */
    int ntype      = 0;
    for (int i = 0; (i < atoms->nr); i++)
    {
        bool found = false;
        int  itp   = atoms->atom[i].type;
        mtop->moltype[0].atoms.atom[i]     = atoms->atom[i];
        mtop->moltype[0].atoms.atomtype[i] = put_symtab(symtab, *atoms->atomtype[i]);
        mtop->moltype[0].atoms.atomname[i] = put_symtab(symtab, *atoms->atomname[i]);
        GMX_RELEASE_ASSERT(atoms->atom[i].resind <= atoms->nres, gmx::formatString("Residue number inconsistency. Atom %d has res number %d out of %d", i, atoms->atom[i].resind, atoms->nres).c_str());
        GMX_RELEASE_ASSERT(atoms->atom[i].resind >= 0, gmx::formatString("Residue number inconsistency. Atom %d has res number %d out of %d", i, atoms->atom[i].resind, atoms->nres).c_str());
        t_atoms_set_resinfo(&mtop->moltype[0].atoms, i, symtab,
                            *atoms->resinfo[atoms->atom[i].resind].name,
                            atoms->atom[i].resind-1,
                            atoms->resinfo[atoms->atom[i].resind].ic,
                            atoms->resinfo[atoms->atom[i].resind].chainid,
                            atoms->resinfo[atoms->atom[i].resind].chainnum);

        for (int j = 0; !found && (j < i); j++)
        {
            found = (itp == atoms->atom[j].type);
        }
        if (!found)
        {
            //ntype++;
            //What if the types are not correct? Need to check.
            // For instance if the atom types are 0 2 3 there will be 3
            // different types but the arrays can not be indexed with a 3.
            ntype = std::max(ntype, itp+1);
        }
    }

    snew(mtop->groups.grpname, ntype);
    snew(mtop->groups.grps[egcENER].nm_ind, ntype);

    int  ind = 0;
    for (int i = 0; i < atoms->nr; i++)
    {
        char  *atp   = *atoms->atomtype[i];
        bool   found = false;
        for (int j = 0; !found && (j < i); j++)
        {
            found = (strcmp(atp, *atoms->atomtype[j]) == 0);
        }
        if (!found)
        {
            mtop->groups.grpname[ind]              = put_symtab(symtab, atp);
            mtop->groups.grps[egcENER].nm_ind[ind] = ind;
            ind++;
        }
    }

    mtop->ffparams.atnr             = ntype;
    mtop->ffparams.reppow           = 12;

    if (nullptr != tabfn)
    {
        mtop->groups.grps[egcENER].nr   = ntype;
        ir->opts.ngener                 = ntype;
        srenew(ir->opts.egp_flags, ntype*ntype);
        for (int k = 0; k < ntype; k++)
        {
            for (int m = k; m < ntype; m++)
            {
                ir->opts.egp_flags[k*ntype + m] |= EGP_TABLE;
            }
        }
    }

    nonbondedFromPdToMtop(mtop, pd, nullptr);

    gmx_mtop_finalize(mtop);
    /* Create a charge group block */
    stupid_fill_block(&(mtop->moltype[0].cgs), atoms->nr, false);
    
    return mtop;
}

void excls_to_blocka(int natom, t_excls excls_[], t_blocka *blocka)
{
    int i, j, k, nra;

    if (blocka->nr < natom)
    {
        srenew(blocka->index, natom+1);
        for (int i = blocka->nr; i < natom+1; i++)
        {
            blocka->index[i] = 0;
        }
    }
    nra = 0;
    for (i = 0; (i < natom); i++)
    {
        nra += excls_[i].nr;
    }
    snew(blocka->a, nra+1);
    nra = 0;
    for (i = j = 0; (i < natom); i++)
    {
        blocka->index[i] = nra;
        for (k = 0; (k < excls_[i].nr); k++)
        {
            blocka->a[j++] = excls_[i].e[k];
        }
        nra += excls_[i].nr;
    }
    blocka->index[natom] = nra;
    blocka->nr           = natom;
    blocka->nra          = nra;
}

void put_in_box(int natom, matrix box, rvec x[], real dbox)
{
    int  i, m;
    rvec xmin, xmax, xcom;

    clear_rvec(xcom);
    copy_rvec(x[0], xmin);
    copy_rvec(x[0], xmax);
    for (i = 0; (i < natom); i++)
    {
        rvec_inc(xcom, x[i]);
        for (m = 0; (m < DIM); m++)
        {
            if (xmin[m] > x[i][m])
            {
                xmin[m] = x[i][m];
            }
            else if (xmax[m] < x[i][m])
            {
                xmax[m] = x[i][m];
            }
        }
    }
    for (m = 0; (m < DIM); m++)
    {
        xcom[m]  /= natom;
        box[m][m] = (dbox+xmax[m]-xmin[m]);
    }
}

static int get_subtype(directive d, int ftype)
{
    int i;
    for (i = 1; i < 20; i++)
    {
        if (!(d == d_angles && i == 7))
        {
            if (ifunc_index(d, i) == ftype)
            {
                return i;
            }
        }
    }
    return 1;
}

static void print_bondeds(FILE                               *out,
                          directive                           d,
                          int                                 ftype,
                          const std::vector<TopologyEntry *> &entries)
{
    if (entries.empty())
    {
        return;
    }
    fprintf(out, "[ %s ]\n", dir2str(d));
    fprintf(out, ";atom i");
    for (int j = 1; (j < NRAL(ftype)); j++)
    {
        fprintf(out, "  %5c", j+'i');
    }
    fprintf(out, "   type  parameters\n");
    int subtype = get_subtype(d, ftype);
    for (auto &entry : entries)
    {
        for (auto &j : entry->atomIndices())
        {
            fprintf(out, "  %5d", 1+j);
        }
        fprintf(out, "  %5d", subtype);
        for (int j = 0; (j < NRFPA(ftype)); j++)
        {
            fprintf(out, "  %10g", 0.0);
        }
        fprintf(out, "; %s\n", entry->id().id().c_str());
    }
    fprintf(out, "\n");
}


void write_top(FILE            *out,
               char            *molname,
               t_atoms         *at,
               gmx_bool         bRTPresname,
               const Topology  *topology,
               t_excls          excls[],
               struct gpp_atomtype *atype,
               const Poldata   *pd)
{
     std::map<int, directive> toPrint = {
            { F_CONSTR,       d_constraints },
            { F_CONSTRNC,     d_constraints },
            { F_LJ14,         d_pairs },
            { F_CMAP,         d_cmap },
            { F_POLARIZATION, d_polarization },
            { F_THOLE_POL,    d_thole_polarization },
            { F_VSITE2,       d_vsites2 },
            { F_VSITE3,       d_vsites3 },
            { F_VSITE3FD,     d_vsites3 },
            { F_VSITE3FAD,    d_vsites3 },
            { F_VSITE3OUT,    d_vsites3 },
            { F_VSITE4FD,     d_vsites4 },
            { F_VSITE4FDN,    d_vsites4 }
        };
    if (at && atype)
    {
        std::vector<int> cgnr;
        cgnr.resize(at->nr, 0);
        for(int i = 0; i < at->nr; i++)
        {
            cgnr[i] = i+1;
        }
        fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", pd->getNexcl());
        print_atoms(out, atype, at, cgnr.data(), bRTPresname);
        for (auto &fs : pd->forcesConst())
        {
            auto iType = fs.first;
            if (!topology->hasEntry(iType))
            {
                continue;
            }
            auto fType = fs.second.fType();
            if (InteractionType::BONDS == fs.first)
            {
                print_bondeds(out, d_bonds, fType, topology->entry(iType));
            }
            else if (InteractionType::ANGLES == iType || InteractionType::LINEAR_ANGLES == iType)
            {
                print_bondeds(out, d_angles, fType, topology->entry(iType));
            }
            else if (InteractionType::PROPER_DIHEDRALS == iType || InteractionType::IMPROPER_DIHEDRALS == iType)
            {
                print_bondeds(out, d_dihedrals, fType, topology->entry(iType));
            }
            else if (toPrint.end() != toPrint.find(fType))
            {
                print_bondeds(out, toPrint.find(fType)->second, fType, topology->entry(iType));
            }
        }
        print_excl(out, at->nr, excls);
    }
}

void print_top_header(FILE                    *fp,
                      const Poldata           *pd,
                      bool                     bPol,
                      std::vector<std::string> commercials,
                      bool                     bItp)
{
    std::string   gt_old, gt_type;
    auto qt          = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));

    fprintf(fp, ";\n");
    fprintf(fp, "; Topology generated by alexandria gentop.\n");
    fprintf(fp, "; Watch this space for information & commercials.\n");
    for (auto i = commercials.begin(); (i < commercials.end()); ++i)
    {
        fprintf(fp, "; %s\n", i->c_str());
    }
    fprintf(fp, ";\n");
    if (!bItp)
    {
        fprintf(fp, "[ defaults ]\n");
        fprintf(fp, "; nbfunc         comb-rule       gen-pairs       fudgeLJ     fudgeQQ\n");
        auto ftype = pd->findForcesConst(InteractionType::VDW).fType();
        std::string ff;
        if (ftype == F_LJ)
        {
            ff.assign("LJ");
        }
        auto combRule = getCombinationRule(pd->findForcesConst(InteractionType::VDW));
        fprintf(fp, "%-15s  %-15s no           %10g  %10g\n\n",
                ff.c_str(),
                ecomb_names[combRule],
                1.0, 1.0);

        fprintf(fp, "[ atomtypes ]\n");
        fprintf(fp, "%-7s%-6s  %6s  %11s  %10s  %5s %-s  %s\n",
                ";atype ", "btype", "at.num", "mass", "charge", "ptype",
                "Van_der_Waals", "Ref_Enthalpy");

        gt_old = "";

        auto vdw = pd->findForcesConst(InteractionType::VDW);
        for (const auto &aType : pd->particleTypesConst())
        {
            gt_type    = aType.id().id();
            auto btype = aType.interactionTypeToIdentifier(InteractionType::BONDS);
            if ((0 ==  gt_old.size()) || (gt_old.compare(gt_type) != 0))
            {
                auto sgt_type= aType.interactionTypeToIdentifier(InteractionType::POLARIZATION);
                auto vdwtype = aType.interactionTypeToIdentifier(InteractionType::VDW);
                double sigma = 0, epsilon = 0, gamma = 0;
                if (!vdwtype.id().empty())
                {
                    auto myvdw = vdw.findParametersConst(vdwtype);
                    sigma      = myvdw["sigma"].value();
                    epsilon    = myvdw["epsilon"].value();
                    gamma      = myvdw["gamma"].value();
                }
                fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f %s %g %g %g %g\n",
                        gt_type.c_str(), 
                        !btype.id().empty() ? btype.id().c_str() : gt_type.c_str(), 
                        aType.atomnumber(), 
                        aType.mass(), 0.0,
                        ptype_str[aType.gmxParticleType()],
                        sigma, epsilon, gamma,
                        aType.refEnthalpy());
                if (false && bPol)
                {
                    if (strcasecmp(ff.c_str(), "LJ") == 0)
                    {
                        fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  S     0  0\n",
                                sgt_type.id().c_str(),
                                sgt_type.id().c_str(),
                                0, 0.0, 0.0);
                    }
                    else
                    {
                        fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  S     0  0  0\n",
                                sgt_type.id().c_str(),
                                sgt_type.id().c_str(),
                                0, 0.0, 0.0);
                    }
                }
            }
            gt_old = gt_type;
        }
        fprintf(fp, "\n");
        if (iChargeType != ChargeType::Point)
        {
            auto eem = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
            fprintf(fp, "[ distributed_charges ]\n");
            for (const auto &atype : pd->particleTypesConst())
            {
                auto ztype     = atype.interactionTypeToIdentifier(InteractionType::CHARGEDISTRIBUTION);
                auto eep       = eem.findParametersConst(ztype);
                auto shellName = atype.interactionTypeToIdentifier(InteractionType::POLARIZATION);
                if (ChargeType::Slater == iChargeType)
                {
                    fprintf(fp, "%-7s  2  %d  %g\n", atype.id().id().c_str(),
                            atype.row(), eep["zeta"].value());
                }
                else if (ChargeType::Gaussian == iChargeType)
                {
                    fprintf(fp, "%-7s  1  %g\n", atype.id().id().c_str(),
                            eep["zeta"].value());
                }
            }
            fprintf(fp, "\n");
        }
    }
}

void calc_rotmatrix(rvec target_vec, rvec ref_vec, matrix rotmatrix)
{
    rvec au = {0, 0, 0};
    rvec bu = {0, 0, 0};

    svmul((1.0/norm(target_vec)), target_vec, au);
    svmul((1.0/norm(ref_vec)), ref_vec, bu);

    rotmatrix[0][0] = bu[0]*au[0];
    rotmatrix[0][1] = bu[0]*au[1];
    rotmatrix[0][2] = bu[0]*au[2];
    rotmatrix[1][0] = bu[1]*au[0];
    rotmatrix[1][1] = bu[1]*au[1];
    rotmatrix[1][2] = bu[1]*au[2];
    rotmatrix[2][0] = bu[2]*au[0];
    rotmatrix[2][1] = bu[2]*au[1];
    rotmatrix[2][2] = bu[2]*au[2];
}

} // namespace alexandria
