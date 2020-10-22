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

#include "gromacs/gmxpreprocess/convparm.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
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
#include "plistwrapper.h"
#include "poldata.h"
#include "units.h"

namespace alexandria
{

bool is_planar(rvec xi, rvec xj, rvec xk,
               rvec xl, t_pbc *pbc,
               real phi_toler)
{
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real phi;

    phi = RAD2DEG*dih_angle(xi, xj, xk, xl, pbc, r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);

    return (fabs(phi) < phi_toler);
}

bool is_linear(rvec xi, rvec xj,
               rvec xk, t_pbc *pbc,
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

void add_excl(t_excls *excls, int e)
{
    int i;

    for (i = 0; (i < excls->nr); i++)
    {
        if (excls->e[i] == e)
        {
            return;
        }
    }
    srenew(excls->e, excls->nr+1);
    excls->e[excls->nr++] = e;
}

void add_excl_pair(t_excls excls[], int e1, int e2)
{
    if (e1 != e2)
    {
        add_excl(&excls[e1], e2);
        add_excl(&excls[e2], e1);
    }
}

void remove_excl(t_excls *excls, int remove)
{
    int i;

    for (i = remove+1; i < excls->nr; i++)
    {
        excls->e[i-1] = excls->e[i];
    }

    excls->nr--;
}

void let_shells_see_shells(t_excls excls[], t_atoms *atoms, gpp_atomtype_t atype)
{
    int i, k, ak;
    for (i = 0; i < atoms->nr; i++)
    {
        if (get_atomtype_ptype(atoms->atom[i].type, atype) == eptShell)
        {
            for (k = 0; k < excls[i].nr; )
            {
                ak = excls[i].e[k];
                if (get_atomtype_ptype(atoms->atom[ak].type, atype) == eptShell)
                {
                    remove_excl(&(excls[i]), k);
                }
                else
                {
                    k++;
                }
            }
        }
    }
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

void cp_plist(t_params                   plist[],
              int                        ftype,
              InteractionType            itype,
              std::vector<PlistWrapper> &plist_)
{
    if (plist[ftype].nr > 0)
    {
        PlistWrapper pw(itype, ftype);
        for (int i = 0; (i < plist[ftype].nr); i++)
        {
            for (int j = interaction_function[ftype].nratoms; j < MAXATOMLIST; j++)
            {
                plist[ftype].param[i].a[j] = 0;
            }
            for (int j = interaction_function[ftype].nrfpA; j < MAXFORCEPARAM; j++)
            {
                plist[ftype].param[i].c[j] = 0;
            }
            pw.addParam(plist[ftype].param[i]);
        }
        plist_.push_back(pw);
    }
}

real calc_r13(const Poldata     *pd,
              const std::string &aai,
              const std::string &aaj,
              const std::string &aak,
              const real         angle)
{
    double r12 = 0, r23 = 0, r13 = 0;

    Identifier aij ({aai, aaj}, CanSwap::Yes);
    Identifier ajk ({aaj, aak}, CanSwap::Yes);

    std::string type("bondlength");
    auto bij = pd->findForcesConst(eitBONDS).findParameterTypeConst(aij, type);
    auto bjk = pd->findForcesConst(eitBONDS).findParameterTypeConst(ajk, type);

    r12 = convertToGromacs(bij.value(), bij.unit());
    r23 = convertToGromacs(bjk.value(), bjk.unit());

    r13 = std::sqrt((r12*r12) + (r23*r23) - (2*r12*r23*std::cos(DEG2RAD*angle)));

    return r13;
}

real calc_relposition(const Poldata     *pd,
                      const std::string  aai,
                      const std::string  aaj,
                      const std::string  aak)
{
    double b0 = 0, b1 = 0, relative_position = 0;

    Identifier aij({aai, aaj}, CanSwap::Yes);
    Identifier ajk({aaj, aak}, CanSwap::Yes);

    std::string type("bondlength");
    auto bij = pd->findForcesConst(eitBONDS).findParameterTypeConst(aij, type);
    auto bjk = pd->findForcesConst(eitBONDS).findParameterTypeConst(ajk, type);

    b0 = convertToGromacs(bij.value(), bij.unit());
    b1 = convertToGromacs(bjk.value(), bjk.unit());

    relative_position = (b1/(b0+b1));

    return relative_position;
}

immStatus updatePlist(const Poldata             *pd,
                      std::vector<PlistWrapper> *plist,
                      const t_atoms              *atoms,
                      bool                       bBASTAT,
                      const std::string         &molname,
                      std::vector<std::string>  *errors)
{
    for (auto pw = plist->begin(); pw < plist->end(); ++pw)
    {
        auto iType   = pw->getItype();
        auto fs      = pd->findForcesConst(iType);
        auto canSwap = fs.canSwap();
        pw->setFtype(fs.fType());
        auto nratoms = interaction_function[fs.fType()].nratoms;
        int bondOrder_index = 0;
        for (auto pwi = pw->beginParam(); pwi < pw->endParam(); ++pwi)
        {
            std::vector<std::string> bondAtomType;
            bool bondNamesOK = true;
            // Store last atom type in case something goes wrong
            std::string lastAtype;
            for(int i = 0; i < nratoms && bondNamesOK; i++)
            {
                std::string tmp;
                GMX_RELEASE_ASSERT(i < atoms->nr, "BAH");
                lastAtype.assign(*atoms->atomtype[pwi->a[i]]);
                bondNamesOK = pd->atypeToBtype(lastAtype, &tmp);
                bondAtomType.push_back(tmp);
            }
            if (bondNamesOK)
            {
                Identifier bondId;
                if (eitBONDS == iType)
                {
                    bondId = Identifier(bondAtomType,
                                        pw->bondOrder(bondOrder_index++),
                                        canSwap);
                }
                else
                {
                    bondId = Identifier(bondAtomType, canSwap);
                }
                int n = 0;
                if (fs.parameterExists(bondId))
                {
                    for(const auto &fp : fs.findParametersConst(bondId))
                    {
                        pwi->c[n++] = convertToGromacs(fp.second.value(), fp.second.unit());
                    }
                }
                else if (!bBASTAT)
                {
                    errors->push_back(gmx::formatString("Could not find information for %s in %s - ftype %s\n",
                                                        bondId.id().c_str(), molname.c_str(), interaction_function[fs.fType()].longname).c_str());
                    if (iType == eitBONDS)
                    {
                        return immNotSupportedBond;
                    }
                    else if (iType ==  eitANGLES)
                    {
                        return immNotSupportedAngle;
                    }
                    else
                    {
                        return immNotSupportedDihedral;
                    }
                }
            }
            else
            {
                errors->push_back(gmx::formatString("Unsupported atom type: %s!\n", lastAtype.c_str()));
                return immAtomTypes;
            }
        }
    }
    return immOK;
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

static void getLjParams(const Poldata     *pd,
                        const std::string &ai,
                        const std::string &aj,
                        double            *c6,
                        double            *cn)
{
    auto fa       = pd->findForcesConst(eitVDW);
    Identifier idI({ai}, CanSwap::Yes); 
    Identifier idJ({aj}, CanSwap::Yes); 
    auto sigmaI   = fa.findParameterTypeConst(idI, "sigma").value();
    auto sigmaJ   = fa.findParameterTypeConst(idJ, "sigma").value();
    auto epsilonI = fa.findParameterTypeConst(idI, "epsilon").value();
    auto epsilonJ = fa.findParameterTypeConst(idJ, "epsilon").value();

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

static void getBhamParams(const Poldata     *pd,
                          const std::string &ai,
                          const std::string &aj,
                          double            *a,
                          double            *b,
                          double            *c)
{
    auto fa       = pd->findForcesConst(eitVDW);
    Identifier idI({ai}, CanSwap::Yes);
    Identifier idJ({aj}, CanSwap::Yes);
    auto sigmaI   = fa.findParameterTypeConst(idI, "sigma").value();
    auto sigmaJ   = fa.findParameterTypeConst(idJ, "sigma").value();
    auto epsilonI = fa.findParameterTypeConst(idI, "epsilon").value();
    auto epsilonJ = fa.findParameterTypeConst(idJ, "epsilon").value();
    auto gammaI   = fa.findParameterTypeConst(idI, "gamma").value();
    auto gammaJ   = fa.findParameterTypeConst(idJ, "gamma").value();

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
                           t_atoms       *atoms,
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
    auto forcesVdw = pd->findForcesConst(eitVDW);
    auto ftypeVdW  = forcesVdw.fType();
    typedef struct
    {
        std::string name;
        int         ptype;
    } mytype;
    std::vector<mytype> mytypes;
    mytypes.resize(ntype);
    for (int i = 0; i < atoms->nr; i++)
    {
        auto itype = atoms->atom[i].type;
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
                            getLjParams(pd, mytypes[i].name,
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
                            getBhamParams(pd, mytypes[i].name,
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

void polarizabilityFromPdToMtop(gmx_mtop_t     *mtop,
                                gmx_localtop_t *ltop,
                                t_atoms        *atoms,
                                const Poldata  *pd,
                                unsigned int    ft)
{
    for (auto i = 0; i < ltop->idef.il[ft].nr; i += interaction_function[ft].nratoms+1)
    {
        auto   tp      = ltop->idef.il[ft].iatoms[i];
        auto   ai      = ltop->idef.il[ft].iatoms[i+1];
        std::string ptype;
        if (pd->atypeToPtype(*atoms->atomtype[ai], &ptype))
        {
            Identifier idS({ptype}, CanSwap::No);
            auto param = pd->findForcesConst(eitPOLARIZATION).findParameterTypeConst(idS, "alpha");
            auto alpha = param.value();
            auto unit  = param.unit(); 
            mtop->ffparams.iparams[tp].polarize.alpha =
                ltop->idef.iparams[tp].polarize.alpha = convertToGromacs(alpha, unit);
        }
    }
}

void plist_to_mtop(const std::vector<PlistWrapper> &plist,
                   gmx_mtop_t                      *mtop_)
{
    int ffparamsSize = mtop_->ffparams.numTypes();
    for (auto &pw : plist)
    {
        int ftype  = pw.getFtype();
        int nra    = NRAL(ftype);
        int nrfp   = NRFPA(ftype);
        int nratot = pw.nParam()*(1+nra);
        if (nratot > 0 && debug)
        {
            fprintf(debug, "There are %d interactions of type %s\n", nratot/(nra+1),
                    interaction_function[ftype].name);
        }
        //mtop_->moltype[0].ilist[ftype].iatoms.resize(nratot, {0});
        /* For generating pairs */
        for (auto j = pw.beginParam(); (j < pw.endParam()); ++j)
        {
            std::vector<real> c;
            c.resize(MAXFORCEPARAM, 0);
            if (ftype == F_LJ14)
            {
                int ati = mtop_->moltype[0].atoms.atom[j->a[0]].type;
                int atj = mtop_->moltype[0].atoms.atom[j->a[1]].type;
                int tp  = ati*mtop_->ffparams.atnr+atj;
                c[0] = c[nrfp]   = mtop_->ffparams.iparams[tp].lj.c6;
                c[1] = c[nrfp+1] = mtop_->ffparams.iparams[tp].lj.c12;
            }
            else
            {
                for (int l = 0; (l < nrfp); l++)
                {
                    c[l] = j->c[l];
                    if (NOTSET == c[l])
                    {
                        c[l] = 0;
                    }
                    c[l+nrfp] = c[l];
                }
            }
            double reppow = 12.0;
            int    type   = enter_params(&mtop_->ffparams, ftype, c.data(), 0, reppow, ffparamsSize, false);
            mtop_->moltype[0].ilist[ftype].iatoms.push_back(type);
            for (int m = 0; m < interaction_function[ftype].nratoms; m++)
            {
                mtop_->moltype[0].ilist[ftype].iatoms.push_back(j->a[m]);
            }
        }
    }
}

gmx_mtop_t *do_init_mtop(const Poldata                   *pd,
                         char                           **molname,
                         t_atoms                         *atoms,
                         const std::vector<PlistWrapper> &plist,
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


    /*Count the number of atom types in the molecule*/
    int ntype      = 0;
    for (int i = 0; (i < atoms->nr); i++)
    {
        bool found = false;
        int  itp   = atoms->atom[i].type;
        mtop->moltype[0].atoms.atom[i] = atoms->atom[i];
        for (int j = 0; !found && (j < i); j++)
        {
            found = (itp == atoms->atom[j].type);
        }
        if (!found)
        {
            ntype++;
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

    nonbondedFromPdToMtop(mtop, atoms, pd, nullptr);

    gmx_mtop_finalize(mtop);
    /* Create a charge group block */
    stupid_fill_block(&(mtop->moltype[0].cgs), atoms->nr, false);
    plist_to_mtop(plist, mtop);

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

void mtop_update_cgs(gmx_mtop_t *mtop)
{
    for (auto &i : mtop->moltype)
    {
        if (i.atoms.nr > i.cgs.nr)
        {
            i.cgs.nr           = i.atoms.nr;
            i.cgs.nalloc_index = i.atoms.nr+1;
            srenew(i.cgs.index, i.cgs.nr+1);
            for (int j = 0; (j <= i.cgs.nr); j++)
            {
                i.cgs.index[j] = j;
            }
        }
    }
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

void write_zeta_q(FILE       *fp,
                  QgenAcm    *qgen,
                  t_atoms    *atoms,
                  ChargeType  iChargeType)
{
    int    i, ii, k, row;
    double zeta, q;
    bool   bAtom, bTypeSet;

    if (nullptr == qgen)
    {
        return;
    }

    fprintf(fp, "[ charge_spreading ]\n");
    fprintf(fp, "; This section describes additional atom type properties.\n");
    fprintf(fp, "; Spreading type (stype) can be either Gaussian (AXg) or Slater (AXs).\n");
    fprintf(fp, "; The zeta are the same for atoms of the same type, and all but the last\n");
    fprintf(fp, "; charge as well. The final charge is different between atoms however,\n");
    fprintf(fp, "; and it is listed below in the [ atoms ] section.\n");
    fprintf(fp, "; atype stype  nq%s      zeta          q  ...\n",
            (eqtSlater == iChargeType) ? "  row" : "");

    k = -1;
    for (i = 0; (i < atoms->nr); i++)
    {
        bAtom = (atoms->atom[i].ptype == eptAtom);
        if (bAtom)
        {
            k++;
        }
        if (k == -1)
        {
            gmx_fatal(FARGS, "The first atom must be a real atom, not a shell");
        }
        bTypeSet = false;
        for (ii = 0; !bTypeSet && (ii < i); ii++)
        {
            bTypeSet = (atoms->atom[ii].type == atoms->atom[i].type);
        }
        if (!bTypeSet)
        {
            fprintf(fp, "%5s %6s",
                    *atoms->atomtype[i],
                    chargeTypeName(iChargeType).c_str());
            row   = qgen->getRow(k);
            q     = qgen->getQ(k);
            zeta  = qgen->getZeta(k);
            if ((row != NOTSET) && (q != NOTSET) && (zeta != NOTSET))
            {
                atoms->atom[i].q      =
                    atoms->atom[i].qB = q;
                if (!bTypeSet)
                {
                    if (eqtSlater == iChargeType)
                    {
                        fprintf(fp, "  %4d", row);
                    }
                    fprintf(fp, " %10f", zeta);
                    fprintf(fp, " %10f", q);
                }
            }
            if (!bTypeSet)
            {
                fprintf(fp, "\n");
            }
        }
    }
    fprintf(fp, "\n");
}

void write_zeta_q2(QgenAcm        *qgen,
                   gpp_atomtype_t  atype,
                   t_atoms        *atoms,
                   ChargeType      iChargeType)
{
    FILE      *fp;
    int        i, k, row;
    double     zeta, q, qtot;
    gmx_bool   bAtom;

    if (nullptr == qgen)
    {
        return;
    }

    fp = fopen("zeta_q.txt", "w");
    k  = -1;
    for (i = 0; (i < atoms->nr); i++)
    {
        bAtom = (atoms->atom[i].ptype == eptAtom);
        if (bAtom)
        {
            k++;
        }
        if (k == -1)
        {
            gmx_fatal(FARGS, "The first atom must be a real atom, not a shell");
        }
        fprintf(fp, "%6s  %5s", chargeTypeName(iChargeType).c_str(),
                get_atomtype_name(atoms->atom[i].type, atype));
        qtot = 0;
        row   = qgen->getRow(k);
        q     = qgen->getQ(k);
        zeta  = qgen->getZeta(k);
        if ((row != NOTSET) && (q != NOTSET) && (zeta != NOTSET))
        {
            qtot += q;
            fprintf(fp, "%5d %10g %10g", row, zeta, q);
        }
        atoms->atom[i].q = qtot;
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fclose(fp);
}

int get_subtype(directive d, int ftype)
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

void print_bondeds(FILE                            *out,
                   directive                        d,
                   int                              plist_ftype,
                   int                              print_ftype,
                   const std::vector<PlistWrapper> &plist)
{
    auto p = SearchPlist(plist, plist_ftype);

    if (plist.end() == p || p->nParam() == 0)
    {
        return;
    }
    fprintf(out, "[ %s ]\n", dir2str(d));
    fprintf(out, ";atom i");
    for (int j = 1; (j < NRAL(print_ftype)); j++)
    {
        fprintf(out, "  %5c", j+'i');
    }
    fprintf(out, "   type  parameters\n");
    int subtype = get_subtype(d, print_ftype);
    for (auto i = p->beginParam(); (i < p->endParam()); ++i)
    {
        for (int j = 0; (j < NRAL(print_ftype)); j++)
        {
            fprintf(out, "  %5d", 1+i->a[j]);
        }
        fprintf(out, "  %5d", subtype);
        for (int j = 0; (j < NRFPA(print_ftype)); j++)
        {
            fprintf(out, "  %10g", i->c[j]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

void write_top(FILE                            *out,
               char                            *molname,
               t_atoms                         *at,
               gmx_bool                         bRTPresname,
               const std::vector<PlistWrapper> &plist,
               t_excls                          excls[],
               gpp_atomtype_t                   atype,
               int                             *cgnr,
               int                              nrexcl,
               const Poldata                   *pd)
{
    if (at && atype && cgnr)
    {
        fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", nrexcl);
        print_atoms(out, atype, at, cgnr, bRTPresname);
        for (auto &fs : pd->forcesConst())
        {
            auto iType = fs.first;
            auto fType = fs.second.fType();
            if (eitBONDS == fs.first)
            {
                print_bondeds(out, d_bonds, fType, fType, plist);
            }
            else if (eitANGLES == iType || eitLINEAR_ANGLES == iType)
            {
                print_bondeds(out, d_angles, fType, fType, plist);
            }
            else if (eitPROPER_DIHEDRALS == iType || eitIMPROPER_DIHEDRALS == iType)
            {
                print_bondeds(out, d_dihedrals, fType, fType, plist);
            }
        }
        print_bondeds(out, d_constraints, F_CONSTR, F_CONSTR, plist);
        print_bondeds(out, d_constraints, F_CONSTRNC, F_CONSTRNC, plist);
        print_bondeds(out, d_pairs, F_LJ14, F_LJ14, plist);
        print_excl(out, at->nr, excls);
        print_bondeds(out, d_cmap, F_CMAP, F_CMAP, plist);
        print_bondeds(out, d_polarization, F_POLARIZATION, F_POLARIZATION, plist);
        print_bondeds(out, d_thole_polarization, F_THOLE_POL, F_THOLE_POL, plist);
        print_bondeds(out, d_vsites2, F_VSITE2, F_VSITE2, plist);
        print_bondeds(out, d_vsites3, F_VSITE3, F_VSITE3, plist);
        print_bondeds(out, d_vsites3, F_VSITE3FD, F_VSITE3FD, plist);
        print_bondeds(out, d_vsites3, F_VSITE3FAD, F_VSITE3FAD, plist);
        print_bondeds(out, d_vsites3, F_VSITE3OUT, F_VSITE3OUT, plist);
        print_bondeds(out, d_vsites4, F_VSITE4FD, F_VSITE4FD, plist);
        print_bondeds(out, d_vsites4, F_VSITE4FDN, F_VSITE4FDN, plist);
    }
}

void print_top_header(FILE                    *fp,
                      const Poldata           *pd,
                      gmx_atomprop_t           aps,
                      bool                     bPol,
                      std::vector<std::string> commercials,
                      bool                     bItp)
{
    std::string   gt_old, gt_type;
    int           atomnumber;
    real          mass;
    ChargeType    iChargeType = pd->chargeType();

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
        auto ftype = pd->findForcesConst(eitVDW).fType();
        std::string ff;
        if (ftype == F_LJ)
        {
            ff.assign("LJ");
        }
        auto combRule = getCombinationRule(pd->findForcesConst(eitVDW));
        fprintf(fp, "%-15s  %-15s no           %10g  %10g\n\n",
                ff.c_str(),
                ecomb_names[combRule],
                1.0, 1.0);

        fprintf(fp, "[ atomtypes ]\n");
        fprintf(fp, "%-7s%-6s  %6s  %11s  %10s  %5s %-s  %s\n",
                ";atype ", "btype", "at.num", "mass", "charge", "ptype",
                "Van_der_Waals", "Ref_Enthalpy");

        gt_old = "";

        auto vdw = pd->findForcesConst(eitVDW);
        for (auto aType = pd->getAtypeBegin(); aType != pd->getAtypeEnd(); aType++)
        {
            gt_type    = aType->getType();
            auto btype = aType->id(eitBONDS);
            if (gmx_atomprop_query(aps, epropMass, "", aType->getElem().c_str(), &mass))
            {
                atomnumber = gmx_atomprop_atomnumber(aps, aType->getElem().c_str());
                if ((0 ==  gt_old.size()) || (gt_old.compare(gt_type) != 0))
                {
                    auto sgt_type= aType->id(eitPOLARIZATION);
                    auto vdwtype = aType->id(eitVDW);
                    auto myvdw   = vdw.findParametersConst(vdwtype);
                    fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  A   %g %g %g %s\n",
                            gt_type.c_str(), btype.id().c_str(),
                            atomnumber, mass, 0.0,
                            myvdw["sigma"].value(),
                            myvdw["epsilon"].value(),
                            myvdw["gamma"].value(),
                            aType->getRefEnthalpy().c_str());
                    if (bPol)
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
            }
            gt_old = gt_type;
        }
        fprintf(fp, "\n");
        if (iChargeType != eqtPoint)
        {
            auto eem = pd->findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
            fprintf(fp, "[ distributed_charges ]\n");
            for (auto atype = pd->getAtypeBegin(); atype != pd->getAtypeEnd(); atype++)
            {
                auto ztype     = atype->id(eitELECTRONEGATIVITYEQUALIZATION);
                auto eep       = eem.findParametersConst(ztype);
                auto shellName = atype->id(eitPOLARIZATION);
                if (eqtSlater == iChargeType)
                {
                    fprintf(fp, "%-7s  2  %g  %g\n", atype->getType().c_str(),
                            eep["row"].value(), eep["zeta"].value());
                    // TODO: Check whether this is necessary
                    //if (polar)
                    //{
                    //  fprintf(fp, "%-7s  2  %d  %g\n", shellName.id().c_str(),
                    //          eem->getRow(1), eem->getZeta(1));
                    //}
                }
                else if (eqtGaussian == iChargeType)
                {
                    fprintf(fp, "%-7s  1  %g\n", atype->getType().c_str(),
                            eep["zeta"].value());
                    // TODO see above
                    //if (polar)
                    //{
                    //   fprintf(fp, "%-7s  1  %g\n", shellName.c_str(),
                    //          eem->getZeta(1));
                    //}
                }
                else
                {
                    GMX_RELEASE_ASSERT(false, "Incorrect charge distribution.");
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
