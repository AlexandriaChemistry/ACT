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

#include "mymol.h"

#include <cstdio>

#include <map>
#include <random>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/gmxpreprocess/convparm.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringcompare.h"

#include "forcefieldparameter.h"
#include "molprop_util.h"
#include "mymol_low.h"
#include "units.h"

namespace alexandria
{

static void vsiteType_to_atomType(const std::string &vsiteType, std::string *atomType)
{
    std::size_t pos = vsiteType.find("L");
    *atomType       = vsiteType.substr (0, pos);
}

class MyForceProvider : public gmx::IForceProvider
{
    private:
        std::vector<double> efield_;

    public:

        MyForceProvider()
        {
            efield_.resize(DIM, 0);
        }

        // From IForceProvider
        //! \copydoc IForceProvider::calculateForces()
        void calculateForces(const gmx::ForceProviderInput &forceProviderInput,
                             gmx::ForceProviderOutput      *forceProviderOutput) override;

        /*! \brief Set the electric field
         * \param[in] efield Vector of field values
         */
        void setField(const std::vector<double> &efield);
};

void MyForceProvider::setField(const std::vector<double> &efield)
{
    efield_ = efield;
}

void MyForceProvider::calculateForces(const gmx::ForceProviderInput &forceProviderInput,
                                      gmx::ForceProviderOutput      *forceProviderOutput)
{
    const t_commrec cr      = forceProviderInput.cr_;
    const t_mdatoms mdatoms = forceProviderInput.mdatoms_;
    double          t       = forceProviderInput.t_;

    rvec           *f = as_rvec_array(forceProviderOutput->forceWithVirial_.force_.data());

    for (int dim = 0; dim < DIM; dim++)
    {
        double efield = FIELDFAC*efield_[dim];
        if (efield != 0)
        {
            for (int i = 0; i < mdatoms.nr; ++i)
            {
                f[i][dim] += mdatoms.chargeA[i]*efield;
            }
        }
    }
    if (MASTER(&cr) && nullptr != debug)
    {
        fprintf(debug, "Electric Field. t: %4g  Ex: %4g  Ey: %4g  Ez: %4g\n",
                t, efield_[XX], efield_[YY], efield_[ZZ]);
    }
}


MyMol::MyMol() : gvt_(VsiteType::ALL)

{
    myforce_           = new MyForceProvider;
    snew(symtab_, 1);
    open_symtab(symtab_);
    gromppAtomtype_    = init_atomtype();
    state_         = new t_state;
    state_->flags |= (1<<estX);
    state_->flags |= (1<<estV);
    state_->flags |= (1<<estCGP);
    snew(enerd_, 1);
    snew(fcd_, 1);
    clear_mat(state_->box);
    for (int m = 0; m < DIM; m++)
    {
        state_->box[m][m] = 10.0;
    }
    init_enerdata(1, 0, enerd_);
    init_nrnb(&nrnb_);
}

t_atoms *MyMol::atoms()
{
    if (!mtop_ || mtop_->moltype.size() != 1) 
    {
        GMX_THROW(gmx::InternalError("mtop_ structure not yet initialized"));
    }
    return &mtop_->moltype[0].atoms; 
}
        
const t_atoms &MyMol::atomsConst() const
{
    if (!mtop_ || mtop_->moltype.size() != 1) 
    {
        GMX_THROW(gmx::InternalError("mtop_ structure not yet initialized"));
    }
    return mtop_->moltype[0].atoms; 
}        

bool MyMol::IsSymmetric(real toler)
{
    int       i, j, m;
    real      mm, tm;
    rvec      com, test;
    gmx_bool *bSymm, bSymmAll;

    clear_rvec(com);
    tm = 0;
    for (i = 0; i < atomsConst().nr; i++)
    {
        mm  = atomsConst().atom[i].m;
        tm += mm;
        for (m = 0; (m < DIM); m++)
        {
            com[m] += mm*state_->x[i][m];
        }
    }
    if (tm > 0)
    {
        for (m = 0; m < DIM; m++)
        {
            com[m] /= tm;
        }
    }
    for (i = 0; i < atomsConst().nr; i++)
    {
        rvec_dec(state_->x[i], com);
    }

    snew(bSymm, atomsConst().nr);
    for (i = 0; i < atomsConst().nr; i++)
    {
        bSymm[i] = (norm(state_->x[i]) < toler);
        for (j = i+1; (j < atomsConst().nr) && !bSymm[i]; j++)
        {
            rvec_add(state_->x[i], state_->x[j], test);
            if (norm(test) < toler)
            {
                bSymm[i] = true;
                bSymm[j] = true;
            }
        }
    }
    bSymmAll = true;
    for (i = 0; i < atomsConst().nr; i++)
    {
        bSymmAll = bSymmAll && bSymm[i];
    }
    sfree(bSymm);
    for (i = 0; i < atomsConst().nr; i++)
    {
        rvec_inc(state_->x[i], com);
    }

    return bSymmAll;
}

void MyMol::findInPlaneAtoms(int ca, std::vector<int> &atoms)
{
    int bca = 0;
    /*First try to find the atom bound to the central atom (ca).*/
    for (auto &bi : bondsConst())
    {
        if ((ca == (bi.getAj() - 1) ||
             ca == (bi.getAi() - 1)))
        {
            if (ca == (bi.getAi() - 1))
            {
                bca = (bi.getAj() - 1);
                atoms.push_back(bca);
            }
            else
            {
                bca = (bi.getAi() - 1);
                atoms.push_back(bca);
            }
        }
    }
    /*Now try to find atoms bound to bca, except ca.*/
    for (auto bi : bondsConst())
    {
        if ((ca != (bi.getAj() - 1)   &&
             ca != (bi.getAi() - 1))  &&
            (bca == (bi.getAj() - 1)  ||
             bca == (bi.getAi() - 1)))
        {
            if (bca == (bi.getAi() - 1))
            {
                atoms.push_back(bi.getAj() - 1);
            }
            else
            {
                atoms.push_back(bi.getAi() - 1);
            }
        }
    }
}

void MyMol::findOutPlaneAtoms(int ca, std::vector<int> &atoms)
{
    for (auto &bi : bondsConst())
    {
        if (bi.getBondOrder() == 1  &&
            (ca == (bi.getAj() - 1) ||
             ca == (bi.getAi() - 1)))
        {
            if (ca == (bi.getAi() - 1))
            {
                atoms.push_back(bi.getAj() - 1);
            }
            else
            {
                atoms.push_back(bi.getAi() - 1);
            }
        }
    }
}

bool MyMol::IsVsiteNeeded(std::string    atype,
                          const Poldata *pd)
{
    auto vsite = pd->findVsite(atype);
    if (vsite != pd->getVsiteEnd())
    {
        return true;
    }
    else
    {
        return false;
    }
}

/*
 * Make Linear Angles, Improper Dihedrals, and Virtual Sites
 */
void MyMol::MakeSpecialInteractions(const Poldata *pd,
                                    t_atoms       *atoms)
{
    std::vector < std::vector < unsigned int>> bonds;
    std::vector<int> nbonds;
    t_pbc            pbc;
    real             th_toler = 175;
    real             ph_toler = 5;

    rvec            *x = as_rvec_array(state_->x.data());

    set_pbc(&pbc, epbcNONE, state_->box);

    bonds.resize(atoms->nr);
    for (auto &bi : bondsConst())
    {
        // Store bonds bidirectionally to get the number correct
        bonds[bi.getAi() - 1].push_back(bi.getAj() - 1);
        bonds[bi.getAj() - 1].push_back(bi.getAi() - 1);
    }
    nbonds.resize(atoms->nr);
    for (auto i = 0; i < atoms->nr; i++)
    {
        nbonds[i] = bonds[i].size();
    }
    for (auto i = 0; i < atoms->nr; i++)
    {
        /* Now test initial geometry */
        if ((bonds[i].size() == 2) &&
            is_linear(x[bonds[i][0]], x[i], x[bonds[i][1]],
                      &pbc, th_toler))
        {
            if (nullptr != debug)
            {
                fprintf(debug, "found linear angle %s-%s-%s in %s\n",
                        *atoms->atomtype[bonds[i][0]],
                        *atoms->atomtype[i],
                        *atoms->atomtype[bonds[i][1]],
                        getMolname().c_str());
            }
            gvt_.addLinear(bonds[i][0], i, bonds[i][1]);
        }
        else if ((bonds[i].size() == 3) &&
                 is_planar(x[i], x[bonds[i][0]],
                           x[bonds[i][1]], x[bonds[i][2]],
                           &pbc, ph_toler))
        {
            if (nullptr != debug)
            {
                fprintf(debug, "found planar group %s-%s-%s-%s in %s\n",
                        *atoms->atomtype[i],
                        *atoms->atomtype[bonds[i][0]],
                        *atoms->atomtype[bonds[i][1]],
                        *atoms->atomtype[bonds[i][2]],
                        getMolname().c_str());
            }
            gvt_.addPlanar(i, bonds[i][0], bonds[i][1], bonds[i][2],
                           &nbonds[0]);
        }
        const auto atype(*atoms->atomtype[i]);
        if (IsVsiteNeeded(atype, pd))
        {
            std::vector<int> vAtoms;
            auto             vsite = pd->findVsite(atype);
            if (vsite->type() == VsiteType::IN_PLANE)
            {
                vAtoms.push_back(i);
                findInPlaneAtoms(i, vAtoms);
                if (vsite->ncontrolatoms() == static_cast<int>(vAtoms.size()))
                {
                    gvt_.addInPlane(vsite->ncontrolatoms(),
                                    vsite->nvsite(),
                                    vAtoms[0], vAtoms[1],
                                    vAtoms[2], vAtoms[3]);
                }
            }
            else if (vsite->type() == VsiteType::OUT_OF_PLANE)
            {
                vAtoms.push_back(i);
                findOutPlaneAtoms(i, vAtoms);
                if (vsite->ncontrolatoms() == static_cast<int>(vAtoms.size()))
                {
                    gvt_.addOutPlane(vsite->ncontrolatoms(),
                                     vsite->nvsite(),
                                     vAtoms[0], vAtoms[1], vAtoms[2]);
                }
            }
            if (!bNeedVsites_)
            {
                bNeedVsites_ = true;
            }
        }
    }
    auto anr = atoms->nr;
    gvt_.generateSpecial(pd, true, atoms, &x, &plist_,
                         symtab_, gromppAtomtype_, &excls_, state_);
    bHaveVSites_ = (atoms->nr > anr);
}

/*
 * Make Harmonic Angles, Proper Dihedrals
 * This needs the bonds to be F_BONDS.
 */
void MyMol::MakeAngles(t_atoms *atoms,
                       bool     bDihs,
                       int      nrexcl)
{
    t_nextnb nnb;
    t_restp  rtp;
    t_params plist[F_NRE];

    init_plist(plist);
    for (auto &pw : plist_)
    {
        if (F_BONDS == pw.getFtype())
        {
            pr_alloc(pw.nParam(), &(plist[F_BONDS]));
            auto i     = 0;
            for (auto &pi : pw.paramsConst())
            {
                for (auto j = 0; j < MAXATOMLIST; j++)
                {
                    plist[F_BONDS].param[i].a[j] = pi.a[j];
                }
                for (auto j = 0; j < MAXFORCEPARAM; j++)
                {
                    plist[F_BONDS].param[i].c[j] = pi.c[j];
                }
                i++;
            }
            plist[F_BONDS].nr = i;
            break;
        }
    }
    /* Make Harmonic Angles and Proper Dihedrals */
    snew(excls_, atoms->nr);
    init_nnb(&nnb, atoms->nr, nrexcl + 2);
    gen_nnb(&nnb, plist);

    print_nnb(&nnb, "NNB");
    rtp.bKeepAllGeneratedDihedrals    = bDihs;
    rtp.bRemoveDihedralIfWithImproper = bDihs;
    rtp.bGenerateHH14Interactions     = true;
    rtp.nrexcl                        = nrexcl;

    gen_pad(&nnb, atoms, &rtp, plist, excls_, nullptr, false);

    t_blocka *EXCL;
    snew(EXCL, 1);
    if (debug)
    {
        fprintf(debug, "Will generate %d exclusions for %d atoms\n",
                nrexcl, atoms->nr);
    }
    generate_excl(nrexcl, atoms->nr, plist, &nnb, EXCL);
    for (int i = 0; i < EXCL->nr; i++)
    {
        int ne = EXCL->index[i+1]-EXCL->index[i];
        srenew(excls_[i].e, ne);
        excls_[i].nr = 0;
        for (auto j = EXCL->index[i]; j < EXCL->index[i+1]; j++)
        {
            if (EXCL->a[j] != i)
            {
                excls_[i].e[excls_[i].nr++] = EXCL->a[j];
            }
        }
        // Set the rest of the memory to zero
        for (auto j = excls_[i].nr; j < ne; j++)
        {
            excls_[i].e[j] = 0;
        }
    }
    done_blocka(EXCL);
    sfree(EXCL);
    if (nullptr != debug)
    {
        for (auto i = 0; i < atoms->nr; i++)
        {
            fprintf(debug, "excl %d", i);
            for (auto j = 0; j < excls_[i].nr; j++)
            {
                fprintf(debug, "  %2d", excls_[i].e[j]);
            }
            fprintf(debug, "\n");
        }
    }
    done_nnb(&nnb);

    cp_plist(plist, F_ANGLES, InteractionType::ANGLES, plist_);

    if (bDihs)
    {
        cp_plist(plist, F_PDIHS, InteractionType::PROPER_DIHEDRALS, plist_);
    }
    for (auto i = 0; i < F_NRE; i++)
    {
        if (plist[i].nr > 0)
        {
            sfree(plist[i].param);
        }
    }
}

immStatus MyMol::GenerateAtoms(const Poldata     *pd,
                               t_atoms           *atoms,
                               const std::string &method,
                               const std::string &basis,
                               std::string       *mylot)
{
    double              xx, yy, zz;
    int                 natom = 0;
    immStatus           imm   = immStatus::OK;

    auto ci     = getCalc(method, basis, mylot);
    if (ci < EndExperiment())
    {
        t_param nb;
        memset(&nb, 0, sizeof(nb));
        natom = 0;
        init_t_atoms(atoms, ci->NAtom(), false);
        snew(atoms->atomtype, ci->NAtom());
        snew(atoms->atomtypeB, ci->NAtom());
        int res0 = -1;
        atoms->nres = 0;
        for (auto &cai : ci->calcAtomConst())
        {
            auto myunit = cai.getUnit();
            cai.getCoords(&xx, &yy, &zz);
            int resnr = cai.ResidueNumber();
            if (resnr != res0)
            {
                res0  = resnr;
                atoms->nres += 1;
            }
            state_->x[natom][XX] = convertToGromacs(xx, myunit);
            state_->x[natom][YY] = convertToGromacs(yy, myunit);
            state_->x[natom][ZZ] = convertToGromacs(zz, myunit);

            atoms->atom[natom].q      =
                atoms->atom[natom].qB = 0;
            atoms->atom[natom].resind = resnr;
            t_atoms_set_resinfo(atoms, natom, symtab_, cai.ResidueName().c_str(),
                                atoms->atom[natom].resind, ' ', 
                                cai.chainId(), cai.chain());
            atoms->atomname[natom]    = put_symtab(symtab_, cai.getName().c_str());

            // First set the atomtype
            atoms->atomtype[natom]      =
                atoms->atomtypeB[natom] = put_symtab(symtab_, cai.getObtype().c_str());
            if (pd->hasParticleType(cai.getObtype()))
            {
                auto atype = pd->findParticleType(cai.getObtype());
                atoms->atom[natom].m      =
                    atoms->atom[natom].mB = atype->mass();
                atoms->atom[natom].atomnumber = atype->atomnumber();
                strncpy(atoms->atom[natom].elem, atype->element().c_str(), sizeof(atoms->atom[natom].elem)-1);
                
                natom++;
            }
            else
            {
                error_messages_.push_back(gmx::formatString("Cannot find atomtype %s (atom %d) in poldata\n", 
                                                            cai.getObtype().c_str(), natom));

                return immStatus::AtomTypes;
            }
        }
        for (auto i = 0; i < natom; i++)
        {
            atoms->atom[i].type      =
                atoms->atom[i].typeB = add_atomtype(gromppAtomtype_, symtab_,
                                                    &(atoms->atom[i]),
                                                    *atoms->atomtype[i],
                                                    &nb, 0,
                                                    atoms->atom[i].atomnumber);
        }
        GMX_RELEASE_ASSERT(atoms->nr == natom, "Inconsitency numbering atoms");
    }
    else
    {
        imm = immStatus::LOT;
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Tried to convert %s to gromacs. LOT is %s/%s. Natoms is %d\n",
                getMolname().c_str(),
                method.c_str(), basis.c_str(), natom);
    }

    return imm;
}

immStatus MyMol::checkAtoms(const Poldata *pd,
                            const t_atoms *atoms)
{
    auto nmissing = 0;
    for (auto i = 0; i < atoms->nr; i++)
    {
        const auto atype(*atoms->atomtype[i]);
        if (!pd->hasParticleType(atype))
        {
            printf("Could not find a force field entry for atomtype %s atom %d in compound '%s'\n",
                   *atoms->atomtype[i], i+1,
                   getMolname().c_str());
            nmissing++;
        }
    }
    if (nmissing > 0)
    {
        return immStatus::AtomTypes;
    }
    return immStatus::OK;
}

immStatus MyMol::zetaToAtoms(const Poldata *pd,
                             t_atoms       *atoms)
{
    /* The first time around we add zeta for the core and addShells will
     * take care of the zeta for the shells.
     * For later calls during optimization of zeta also the
     * zeta on the shells will be set. 
     */
    auto qt        = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    const char *ct = "chargetype";
    if (!qt.optionExists(ct))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No option %s in in the force field file", ct).c_str()));
    }
    auto eqtModel = name2ChargeType(qt.optionValue(ct));
    if (eqtModel == ChargeType::Point)
    {
        return immStatus::OK;
    }
    
    auto iType = InteractionType::CHARGEDISTRIBUTION;
    for (auto i = 0; i < atoms->nr; i++)
    {
        auto atype = pd->findParticleType(*atoms->atomtype[i]);
        auto ztype = atype->interactionTypeToIdentifier(iType);
        GMX_RELEASE_ASSERT(qt.parameterExists(ztype), gmx::formatString("No such parameter '%s' for atom type %s in %s",
                                                                        ztype.id().c_str(),
                                                                        atype->id().id().c_str(),
                                                                        interactionTypeToString(iType).c_str()).c_str());

        auto eep   = qt.findParametersConst(ztype);
        const char *zzz =  "zeta";
        if (eep.find(zzz) ==  eep.end())
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Cannot find parameter type %s for atomtype %s",
                                                               zzz, atype->id().id().c_str()).c_str()));
        }
        auto zeta  = eep["zeta"].value();
        
        if (zeta == 0 && eqtModel != ChargeType::Point)
        {
            return immStatus::ZeroZeta;
        }
        
        atoms->atom[i].zetaA     =
            atoms->atom[i].zetaB = zeta;
        atoms->atom[i].row       = eep["row"].value();
    }
    return immStatus::OK;
}

static void fill_atom(t_atom *atom,
                      real m, real q, real mB, real qB, int atomnumber,
                      int atomtype, int ept,
                      real zetaA, real zetaB, int row, int resind)
{
    atom->m             = m;
    atom->mB            = mB;
    atom->q             = q;
    atom->qB            = qB;
    atom->atomnumber    = atomnumber;
    atom->type          = atomtype;
    atom->typeB         = atomtype;
    atom->ptype         = ept;
    atom->zetaA         = zetaA;
    atom->zetaB         = zetaB;
    atom->row           = row;
    atom->resind        = resind;
}
                       
immStatus MyMol::GenerateTopology(FILE              *fp,
                                  const Poldata     *pd,
                                  const std::string &method,
                                  const std::string &basis,
                                  std::string       *mylot,
                                  missingParameters  missing,
                                  const char        *tabfn)
{
    immStatus   imm = immStatus::OK;
    std::string btype1, btype2;

    if (nullptr != debug)
    {
        fprintf(debug, "Generating topology for %s\n", getMolname().c_str());
    }
    generateComposition();
    if (NAtom() <= 0)
    {
        imm = immStatus::AtomTypes;
    }
    t_atoms *atoms = nullptr;
    if (immStatus::OK == imm)
    {
        snew(atoms, 1);
        state_change_natoms(state_, NAtom());
        imm = GenerateAtoms(pd, atoms, method, basis, mylot);
    }
    if (immStatus::OK == imm)
    {
        imm = checkAtoms(pd, atoms);
    }
    if (immStatus::OK == imm)
    {
        imm = zetaToAtoms(pd, atoms);
    }
    /* Store bonds in harmonic potential list first, update type later */
    if (immStatus::OK == imm)
    {
        if (NBond() > 0)
        {
            int  ftb = F_BONDS;
            auto fs  = pd->findForcesConst(InteractionType::BONDS);
            for (auto &bi : bondsConst())
            {
                t_param b;
                memset(&b, 0, sizeof(b));
                b.a[0] = bi.getAi() - 1;
                b.a[1] = bi.getAj() - 1;
                pd->atypeToBtype(*atoms->atomtype[b.a[0]], &btype1);
                pd->atypeToBtype(*atoms->atomtype[b.a[1]], &btype2);
                Identifier bondId({btype1, btype2}, bi.getBondOrder(), CanSwap::Yes);
                // Store the bond order for later usage.
                bondOrder_.insert({std::make_pair(b.a[0],b.a[1]), bi.getBondOrder()});
                // We add the parameter with zero parameters, they will be
                // set further down. However it is important to set the
                // bondorder.
                add_param_to_plist(&plist_, ftb, InteractionType::BONDS,
                                   b, bi.getBondOrder());
            }
            auto pw = SearchPlist(plist_, ftb);
            if (plist_.end() == pw || pw->nParam() == 0)
            {
                imm = immStatus::GenBonds;
            }
        }
    }
    // Check whether we have dihedrals in the force field.
    bool bDih = false;
    if (pd->interactionPresent(InteractionType::PROPER_DIHEDRALS))
    {
        bDih = !pd->findForcesConst(InteractionType::PROPER_DIHEDRALS).empty();
    }
    
    if (immStatus::OK == imm)
    {
        MakeAngles(atoms, bDih, pd->getNexcl());

        MakeSpecialInteractions(pd, atoms);
        imm = updatePlist(pd, &plist_, atoms, missing,
                          getMolname(), &error_messages_);
    }
    if (immStatus::OK == imm)
    {
        /* Center of charge */
        auto atntot = 0;
        rvec coc    = { 0 };
        for (auto i = 0; i < atoms->nr; i++)
        {
            auto atn = atoms->atom[i].atomnumber;
            atntot  += atn;
            for (auto m = 0; m < DIM; m++)
            {
                coc[m] += state_->x[i][m]*atn;
            }
        }
        /* Center of charge */
        svmul((1.0/atntot), coc, CenterOfCharge_);
        for (auto &qp : qProps_)
        {
            qp.second.setCenterOfCharge(CenterOfCharge_);
        }
        addBondVsites(fp, pd, atoms);
        if (pd->polarizable())
        {
            addShells(debug, pd, atoms);
        }
        char **molnameptr = put_symtab(symtab_, getMolname().c_str());
        // Generate mtop
        mtop_ = do_init_mtop(pd, molnameptr, atoms, plist_,
                             inputrec_, symtab_, tabfn);
        excls_to_blocka(atoms->nr, excls_, &(mtop_->moltype[0].excls));
        if (pd->polarizable())
        {
            // Update mtop internals to account for shell type
            srenew(mtop_->atomtypes.atomnumber,
                   get_atomtype_ntypes(gromppAtomtype_));
            for (auto i = 0; i < get_atomtype_ntypes(gromppAtomtype_); i++)
            {
                mtop_->atomtypes.atomnumber[i] = get_atomtype_atomnumber(i, gromppAtomtype_);
            }
            mtop_->ffparams.atnr = get_atomtype_ntypes(gromppAtomtype_);
            // Generate shell data structure
            shellfc_             = init_shell_flexcon(debug, mtop_, 0, 1, false);
        }
        if (nullptr == ltop_)
        {
            // Generate ltop from mtop
            ltop_ = gmx_mtop_generate_local_top(mtop_, false);
        }
    }
    if (immStatus::OK == imm && missing != missingParameters::Generate)
    {
        UpdateIdef(pd, InteractionType::BONDS);
        UpdateIdef(pd, InteractionType::ANGLES);
        if (pd->polarizable())
        {
            UpdateIdef(pd, InteractionType::POLARIZATION);
        }
        UpdateIdef(pd, InteractionType::LINEAR_ANGLES);
        UpdateIdef(pd, InteractionType::IMPROPER_DIHEDRALS);
        if (bDih)
        {
            UpdateIdef(pd, InteractionType::PROPER_DIHEDRALS);
        }
        UpdateIdef(pd, InteractionType::VSITE2);
    }
    if (immStatus::OK != imm && debug)
    {
        for(const auto &emsg : error_messages_)
        {
            fprintf(debug, "%s\n", emsg.c_str());
        }
    }
    qProps_.insert(std::pair<qType, QtypeProps>(qType::Calc, QtypeProps(qType::Calc)));
    qProps_.insert(std::pair<qType, QtypeProps>(qType::Elec, QtypeProps(qType::Elec)));
    return imm;
}

void MyMol::addBondVsites(FILE          *fp,
                          const Poldata *pd,
                          t_atoms       *atoms)
{
    int     atomNrOld = atoms->nr;
    t_param p         = { { 0 } };
    // First add virtual sites for bond shells if needed.
    auto    vs2       = pd->findForcesConst(InteractionType::VSITE2);
    auto    qt        = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    // TODO: add a flag to turn that on or off?
    t_atom vsite_atom = { 0 };
    for (auto b: bondsConst())
    {
        int ai = b.getAi()-1;
        int aj = b.getAj()-1;
        std::string aTypei(*atoms->atomtype[ai]);
        std::string aTypej(*atoms->atomtype[aj]);
        Identifier  bsId({ aTypei, aTypej }, b.getBondOrder(), CanSwap::No);
        // Check whether a vsite is defined for this bond
        std::string v2("v2");
        if (vs2.parameterExists(bsId) && pd->hasParticleType(v2))
        {
            // Yes! We need to add a virtual site
            p.a[0] = atoms->nr;
            p.a[1] = ai;
            p.a[2] = aj;
            auto param = vs2.findParameterTypeConst(bsId, "v2_a");
            p.c[0] = convertToGromacs(param.value(), param.unit());
            add_param_to_plist(&plist_, F_VSITE2, InteractionType::VSITE2, p);
            // Now add the particle
            add_t_atoms(atoms, 1, 0);
            // Add exclusion for the vsite and its constituting atoms
            srenew(excls_, atoms->nr);
            excls_[atoms->nr-1].nr = 0;
            //snew(excls_[atoms->nr-1].e, 2);
            //excls_[atoms->nr-1].e[0] = ai;
            //excls_[atoms->nr-1].e[1] = aj;
        
            auto ptype     = pd->findParticleType(v2);
            auto m         = ptype->paramValue("mass");
            auto q         = ptype->paramValue("charge");
            auto vsid      = ptype->interactionTypeToIdentifier(InteractionType::POLARIZATION);
            auto vstype    = pd->findParticleType(vsid.id());
            auto vszetaid  = vstype->interactionTypeToIdentifier(InteractionType::CHARGEDISTRIBUTION);
            auto vseep     = qt.findParametersConst(vszetaid);
            // Now fill the newatom
            auto zeta      = vseep["zeta"].value();
            auto vsgpp     = add_atomtype(gromppAtomtype_, symtab_, &vsite_atom,
                                          ptype->id().id().c_str(), &p, 0, 0);
            fill_atom(&atoms->atom[atoms->nr-1],
                      m, q, m, q,
                      ptype->atomnumber(),
                      vsgpp, eptVSite,
                      zeta, zeta, ptype->row(), 
                      atoms->atom[ai].resind);
            atoms->atomname[atoms->nr-1] = put_symtab(symtab_, v2.c_str());
            atoms->atomtype[atoms->nr-1] = put_symtab(symtab_, v2.c_str());
            atoms->atomtypeB[atoms->nr-1] = put_symtab(symtab_, v2.c_str());
        }
    }
    if (atoms->nr > atomNrOld && fp)
    {
        fprintf(fp, "Added %d bond vsite(s) for %s\n",
                atoms->nr - atomNrOld, getMolname().c_str());
    }
}

void MyMol::addShells(FILE          *fp,
                      const Poldata *pd,
                      t_atoms       *atoms)
{
    int                    shell  = 0;
    int                    nshell = 0;
    std::map<int, int>     inv_renum;
    std::vector<std::string> newname;
    t_atoms               *newatoms;
    t_excls               *newexcls;
    std::vector<gmx::RVec> newx;
    t_param                p = { { 0 } };
    std::vector<int>       shellRenumber;
    
    /* Calculate the total number of Atom and Vsite particles and
     * generate the renumbering array.
     */
    shellRenumber.resize(atoms->nr, 0);
    for (int i = 0; i < atoms->nr; i++)
    {
        auto atype                  = pd->findParticleType(*atoms->atomtype[i]);
        shellRenumber[i]            = i + nshell;
        inv_renum[shellRenumber[i]] = i;
        if (atype->hasInteractionType(InteractionType::POLARIZATION))
        {
            // TODO: Update if particles can have more than one shell
            nshell++;
        }
    }
    if (fp)
    {
        fprintf(fp, "Found %d shells to be added\n", nshell);
    }
    int nParticles = atoms->nr+nshell;
    state_change_natoms(state_, nParticles);
    
    /* Add Polarization to the plist. */
    auto qt  = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    auto fs  = pd->findForcesConst(InteractionType::POLARIZATION);
    
    
    // Atoms first
    for (int i = 0; i < atoms->nr; i++)
    {
        std::string atomtype(*atoms->atomtype[i]);
        if (atoms->atom[i].ptype == eptAtom ||
            atoms->atom[i].ptype == eptVSite)
        {
            if (pd->hasParticleType(atomtype))
            {
                // TODO: Allow adding multiple shells
                auto fa = pd->findParticleType(atomtype);
                if (fa->hasInteractionType(InteractionType::POLARIZATION))
                {
                    auto ptype = fa->interactionTypeToIdentifier(InteractionType::POLARIZATION);
                    auto param = fs.findParameterTypeConst(ptype, "alpha");
                    auto pol   = convertToGromacs(param.value(), param.unit());
                    if (pol <= 0)
                    {
                        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Polarizability should be positive for %s", fa->id().id().c_str()).c_str()));
                    }
                    // TODO Multiple shell support
                    p.a[0] = shellRenumber[i];
                    p.a[1] = shellRenumber[i]+1;
                    if (bHaveVSites_)
                    {
                        auto vsite = pd->findVsite(atomtype);
                        if (vsite != pd->getVsiteEnd())
                        {
                            pol /= vsite->nvsite();
                        }
                    }
                    p.c[0] = pol;
                    add_param_to_plist(&plist_, F_POLARIZATION,
                                       InteractionType::POLARIZATION, p);
                }
            }
            else
            {
                error_messages_.push_back(gmx::formatString("Cannot find atomtype %s in poldata\n", 
                                                            atomtype.c_str()));
            }
        }
    }
    /* Make new atoms and x arrays. */
    snew(newatoms, 1);
    init_t_atoms(newatoms, nParticles, true);
    snew(newatoms->atomtype, nParticles);
    snew(newatoms->atomtypeB, nParticles);
    snew(newatoms->atomname, nParticles);
    newatoms->nres = atoms->nres;
    newx.resize(newatoms->nr);
    newname.resize(newatoms->nr);

    /* Make a new exclusion array and put the shells in it. */
    snew(newexcls, newatoms->nr);
    /* Add exclusion for F_POLARIZATION. */
    auto pw = SearchPlist(plist_, F_POLARIZATION);
    if (plist_.end() != pw)
    {
        /* Exclude the vsites and the atoms from their own shell. */
        if (pd->getNexcl() >= 0)
        {
            for (auto &j : pw->paramsConst())
            {
                add_excl_pair(newexcls, j.a[0], j.a[1]);
            }
        }

        // Make a copy of the exclusions of the Atom or Vsite for the shell.
        for (auto &j : pw->paramsConst())
        {
            // We know that the Atom or Vsite is 0 as we added it to plist as such.
            int  i0 = inv_renum[j.a[0]];
            for (auto j0 = 0; j0 < excls_[i0].nr; j0++)
            {
                add_excl_pair(newexcls, j.a[0], shellRenumber[excls_[i0].e[j0]]);
                add_excl_pair(newexcls, j.a[1], shellRenumber[excls_[i0].e[j0]]);
            }
        }
        for (auto &j : pw->paramsConst())
        {
            for (auto j0 = 0; j0 < newexcls[j.a[0]].nr; j0++)
            {
                add_excl_pair(newexcls, j.a[1], newexcls[j.a[0]].e[j0]);
            }
        }
    }

    /* Copy the old atoms to the new structures. */
    for (int i = 0; i < atoms->nr; i++)
    {
        newatoms->atom[shellRenumber[i]]      = atoms->atom[i];
        newatoms->atomname[shellRenumber[i]]  = put_symtab(symtab_, *atoms->atomname[i]);
        newatoms->atomtype[shellRenumber[i]]  = put_symtab(symtab_, *atoms->atomtype[i]);
        newatoms->atomtypeB[shellRenumber[i]] = put_symtab(symtab_, *atoms->atomtypeB[i]);
        copy_rvec(state_->x[i], newx[shellRenumber[i]]);
        newname[shellRenumber[i]].assign(*atoms->atomtype[i]);
        int resind = atoms->atom[i].resind;
        t_atoms_set_resinfo(newatoms, shellRenumber[i], symtab_,
                            *atoms->resinfo[resind].name,
                            atoms->resinfo[resind].nr,
                            atoms->resinfo[resind].ic, 
                            atoms->resinfo[resind].chainnum, 
                            atoms->resinfo[resind].chainid);
    }
    t_atom shell_atom = { 0 };
    for (int i = 0; i < atoms->nr; i++)
    {
        if (atoms->atom[i].ptype == eptAtom ||
            atoms->atom[i].ptype == eptVSite)
        {
            std::string atomtype;
            // Shell sits next to the Atom or Vsite
            // TODO make this more precise.
            auto j            = 1+shellRenumber[i];
            auto atomtypeName = get_atomtype_name(atoms->atom[i].type, gromppAtomtype_);
            auto fa           = pd->findParticleType(atomtypeName);
            auto shellid      = fa->interactionTypeToIdentifier(InteractionType::POLARIZATION);
            auto shelltype    = pd->findParticleType(shellid.id());
            auto shellzetaid  = shelltype->interactionTypeToIdentifier(InteractionType::CHARGEDISTRIBUTION);
            auto shelleep     = qt.findParametersConst(shellzetaid);
            // Now fill the newatom
            newatoms->atom[j]               = atoms->atom[i];
            newatoms->atom[j].m             =
                newatoms->atom[j].mB            = shelltype->mass();
            // Shell has no core
            newatoms->atom[j].atomnumber    = 0;
            shell                           = add_atomtype(gromppAtomtype_, symtab_, &shell_atom, shellid.id().c_str(), &p, 0, 1);
            newatoms->atom[j].type          = shell;
            newatoms->atom[j].typeB         = shell;
            newatoms->atomtype[j]           = put_symtab(symtab_, shellid.id().c_str());
            newatoms->atomtypeB[j]          = put_symtab(symtab_, shellid.id().c_str());
            newatoms->atomname[j]           = put_symtab(symtab_, atomtypeName);
            newatoms->atom[j].ptype         = eptShell;
            newatoms->atom[j].zetaA         = shelleep["zeta"].value();
            newatoms->atom[j].zetaB         = newatoms->atom[j].zetaA;
            newatoms->atom[j].row           = shelltype->row();
            newatoms->atom[j].resind        = atoms->atom[i].resind;
            copy_rvec(state_->x[i], newx[j]);

            newatoms->atom[j].q      =
                newatoms->atom[j].qB = shelltype->charge();
            if (bHaveVSites_)
            {
                if (atoms->atom[i].ptype == eptVSite)
                {
                    vsiteType_to_atomType(atomtypeName, &atomtype);
                }
                auto vsite = pd->findVsite(atomtype);
                if (vsite != pd->getVsiteEnd())
                {
                    newatoms->atom[j].q /= vsite->nvsite();
                    newatoms->atom[j].qB = newatoms->atom[j].q;
                }
            }
        }
    }
    /* Copy newatoms to atoms */
    copy_atoms(newatoms, atoms);

    for (int i = 0; i < newatoms->nr; i++)
    {
        copy_rvec(newx[i], state_->x[i]);
        atoms->atomtype[i] = 
            atoms->atomtypeB[i] = put_symtab(symtab_, *newatoms->atomtype[i]);
    }

    /* Update the bond orders */
    auto  boCopy = bondOrder_;
    bondOrder_.clear();
    for(const auto &boc : boCopy)
    {
        bondOrder_.insert({std::make_pair(shellRenumber[boc.first.first],
                                          shellRenumber[boc.first.second]), boc.second});
    }
    /* Copy exclusions, empty the original first */
    sfree(excls_);
    excls_ = newexcls;

    /* Now renumber atoms in all other plist interaction types */
    for (auto i = plist_.begin(); i < plist_.end(); ++i)
    {
        if (i->getFtype() != F_POLARIZATION)
        {
            auto mypar = i->params();
            for (auto j = mypar->begin(); j < mypar->end(); ++j)
            {
                for (int k = 0; k < NRAL(i->getFtype()); k++)
                {
                    j->a[k] = shellRenumber[j->a[k]];
                }
            }
        }
    }
    bHaveShells_ = true;
}

double MyMol::bondOrder(int ai, int aj) const
{
    if (bondOrder_.find({ai, aj}) != bondOrder_.end())
    {
        return bondOrder_.find({ai, aj})->second;
    }
    else if (bondOrder_.find({aj, ai}) != bondOrder_.end())
    {
        return bondOrder_.find({aj, ai})->second;
    }
    return 0;
}

immStatus MyMol::GenerateGromacs(const gmx::MDLogger       &mdlog,
                                 t_commrec                 *cr,
                                 const char                *tabfn,
                                 gmx_hw_info_t             *hwinfo,
                                 ChargeType                 ieqt)
{
    if (gromacsGenerated_)
    {
        return immStatus::OK;
    }
    GMX_RELEASE_ASSERT(nullptr != mtop_, "mtop_ == nullptr. You forgot to call GenerateTopology");

    if (!fr_)
    {
        fr_ = mk_forcerec();
    }
    if (tabfn)
    {
        inputrec_->coulombtype = eelUSER;
    }
    mdModules_ = new std::unique_ptr<gmx::MDModules>(new gmx::MDModules());
    mdModules_->get()->assignOptionsToModules(*(inputrec_->params), nullptr);
    fr_->forceProviders = mdModules_->get()->initForceProviders();
    fr_->forceProviders->addForceProvider(myforce_);
    // Tell gromacs to use the generic kernel only.
    gmx_nonbonded_setup(fr_, true);

    gmx::ArrayRef<const std::string>  tabbfnm;
    init_forcerec(nullptr, mdlog, fr_, nullptr, inputrec_, mtop_, cr,
                  state_->box, tabfn, tabfn, tabbfnm, *hwinfo, nullptr, false, true, -1);
    gmx_omp_nthreads_set(emntBonded, 1);
    init_bonded_threading(nullptr, 1, &fr_->bondedThreading);
    setup_bonded_threading(fr_->bondedThreading, atomsConst().nr, false, ltop_->idef);
    wcycle_    = wallcycle_init(debug, 0, cr);

    MDatoms_  = new std::unique_ptr<gmx::MDAtoms>(new gmx::MDAtoms());
    *MDatoms_ = gmx::makeMDAtoms(nullptr, *mtop_, *inputrec_, false);
    atoms2md(mtop_, inputrec_, -1, nullptr, atomsConst().nr, MDatoms_->get());
    auto mdatoms = MDatoms_->get()->mdatoms();
    f_.resizeWithPadding(state_->natoms);

    if (nullptr != shellfc_)
    {
        make_local_shells(cr, mdatoms, shellfc_);
    }
    if (ChargeType::Slater != ieqt)
    {
        for (auto i = 0; i < mtop_->natoms; i++)
        {
            mdatoms->row[i] = 0;
        }
    }
    gromacsGenerated_ = true;
    return immStatus::OK;
}

immStatus MyMol::computeForces(FILE *fplog, t_commrec *cr, double *rmsf)
{
    auto mdatoms = MDatoms_->get()->mdatoms();
    auto atoms   = atomsConst();
    for (auto i = 0; i < mtop_->natoms; i++)
    {
        mdatoms->chargeA[i] = atoms.atom[i].q;
        mdatoms->typeA[i]   = atoms.atom[i].type;
        mdatoms->zetaA[i]   = atoms.atom[i].zetaA;
        if (mdatoms->zetaB)
        {
            mdatoms->zetaB[i]   = atoms.atom[i].zetaB;
        }
        if (nullptr != debug)
        {
            fprintf(debug, "QQQ Setting q[%d] to %g\n", i, mdatoms->chargeA[i]);
        }
    }
    if (!vsite_)
    {
        vsite_  = new std::unique_ptr<gmx_vsite_t>(new gmx_vsite_t());
        *vsite_ = initVsite(*mtop_, cr);
    }
    unsigned long force_flags = ~0;
    double        t           = 0;
    rvec          mu_tot      = { 0, 0, 0 };
    tensor        force_vir;
    clear_mat (force_vir);
    for (int i = 0; i < mtop_->natoms; i++)
    {
        clear_rvec(f_[i]);
    }
    for (int i = 0; i < F_NRE; i++)
    {
        enerd_->term[i] = 0;
    }
    for (int i = 0; i < enerd_->grpp.nener; i++)
    {
        for (int j = 0; j < egNR; j++)
        {
            enerd_->grpp.ener[j][i] = 0;
        }
    }
    restoreCoordinates();
    constructVsitesGlobal(*mtop_, state_->x);
    immStatus imm =  immStatus::OK;
    if (nullptr != shellfc_)
    {
        if (debug)
        {
            fprintf(debug, "cr->nodide %d mol %s alpha %g\n", cr->nodeid,
                    getMolname().c_str(),
                    mtop_->ffparams.iparams[mtop_->moltype[0].ilist[F_POLARIZATION].iatoms[0]].polarize.alpha);
        }
        auto nnodes = cr->nnodes;
        cr->nnodes  = 1;
        real force2 = 0;
        try
        {
            force2 = relax_shell_flexcon(nullptr, cr, nullptr, false,
                                         nullptr, 0, inputrec_,
                                         true, force_flags, ltop_, nullptr,
                                         enerd_, fcd_, state_,
                                         f_.arrayRefWithPadding(), force_vir, mdatoms,
                                         &nrnb_, wcycle_, nullptr,
                                         &(mtop_->groups), shellfc_,
                                         fr_, t, mu_tot, vsite_->get(),
                                         DdOpenBalanceRegionBeforeForceComputation::no,
                                         DdCloseBalanceRegionAfterForceComputation::no);
        }
        catch (gmx::SimulationInstabilityError &ex)
        {
            fprintf(stderr, "Something wrong minimizing shells for %s. Error code %d\n",
                    getMolname().c_str(), ex.errorCode());
            imm = immStatus::ShellMinimization;
        }
        *rmsf = std::sqrt(force2);
        if (force2 > inputrec_->em_tol && fplog)
        {
            for (int i = 0;  i<F_NRE; i++)
            {
                auto ei = enerd_->term[i];
                if (ei != 0)
                {
                    fprintf(fplog, "E[%s] = %g\n", interaction_function[i].name,
                            ei);
                }
            }
            fprintf(fplog, "Shell minimization did not converge in %d steps for %s. RMS Force = %g.\n",
                    inputrec_->niter, getMolname().c_str(),
                    *rmsf);
            pr_rvecs(fplog, 0, "f", f_.rvec_array(), mtop_->natoms);
            imm = immStatus::ShellMinimization;
        }
        cr->nnodes = nnodes;
    }
    else
    {
        do_force(fplog, cr, nullptr, inputrec_, nullptr, 0,
                 &nrnb_, wcycle_, ltop_,
                 &(mtop_->groups),
                 state_->box, state_->x.arrayRefWithPadding(), nullptr,
                 f_.arrayRefWithPadding(), force_vir, mdatoms,
                 enerd_, fcd_,
                 state_->lambda, nullptr,
                 fr_, vsite_->get(), mu_tot, t,
                 nullptr,
                 force_flags,
                 DdOpenBalanceRegionBeforeForceComputation::no,
                 DdCloseBalanceRegionAfterForceComputation::no);
        *rmsf = 0;
    }
    return imm;
}

void MyMol::symmetrizeCharges(const Poldata  *pd,
                              bool            bSymmetricCharges,
                              const char     *symm_string)
{
    if (bSymmetricCharges)
    {
        symmetric_charges_.clear();
        ConstPlistWrapperIterator bonds = SearchPlist(plist_, InteractionType::BONDS);
        if (plist_.end() != bonds)
        {
            symmetrize_charges(bSymmetricCharges, atoms(), bonds,
                               pd, symm_string, &symmetric_charges_);
        }
    }
    else
    {
        for (auto i = 0; i < mtop_->natoms; i++)
        {
            symmetric_charges_.push_back(i);
        }
    }
}

immStatus MyMol::GenerateAcmCharges(const Poldata       *pd,
                                    t_commrec           *cr,
                                    int                  maxiter,
                                    real                 tolerance)
{
    if (QgenAcm_ == nullptr)
    {
        QgenAcm_ = new QgenAcm(pd, atoms(), totalCharge());
    }
    std::vector<double> qq;
    for (auto i = 0; i < mtop_->natoms; i++)
    {
        qq.push_back(QgenAcm_->getQ(i));
    }
    immStatus imm       = immStatus::OK;
    int       iter      = 0;
    bool      converged = false;
    double    EemRms    = 0;
    do
    {
        if (eQgen::OK == QgenAcm_->generateCharges(debug,
                                                   getMolname().c_str(),
                                                   pd,
                                                   atoms(),
                                                   state_->x,
                                                   bondsConst()))
        {
            if (nullptr != shellfc_)
            {
                double rmsf;
                auto imm = computeForces(nullptr, cr, &rmsf);
                if (imm != immStatus::OK)
                {
                    return imm;
                }
            }
            EemRms = 0;
            for (int i = 0; i < mtop_->natoms; i++)
            {
                auto q_i = QgenAcm_->getQ(i);
                EemRms  += gmx::square(qq[i] - q_i);
                qq[i]    = q_i;
            }
            EemRms   /= mtop_->natoms;
            converged = (EemRms < tolerance) || (nullptr == shellfc_);
            iter++;
        }
        else
        {
            imm = immStatus::ChargeGeneration;
        }
    }
    while (imm == immStatus::OK && (!converged) && (iter < maxiter));
    if (!converged)
    {
        printf("Alexandria Charge Model did not converge to %g. rms: %g\n",
               tolerance, sqrt(EemRms));
    }
    auto qcalc = qTypeProps(qType::Calc);
    qcalc->setQ(atoms());
    return imm;
}

immStatus MyMol::GenerateCharges(const Poldata             *pd,
                                 const gmx::MDLogger       &mdlog,
                                 t_commrec                 *cr,
                                 const char                *tabfn,
                                 gmx_hw_info_t             *hwinfo,
                                 int                        maxiter,
                                 real                       tolerance,
                                 ChargeGenerationAlgorithm  algorithm,
                                 const std::vector<double> &qcustom,
                                 const std::string         &lot)
{
    immStatus imm         = immStatus::OK;
    bool      converged   = false;
    int       iter        = 0;
    auto      qt          = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    auto      iChargeType = name2ChargeType(qt.optionValue("chargetype"));

    GenerateGromacs(mdlog, cr, tabfn, hwinfo, iChargeType);
    if (backupCoordinates_.size() == 0)
    {
        backupCoordinates();
    }
    if (algorithm == ChargeGenerationAlgorithm::Custom)
    {
        GMX_RELEASE_ASSERT(0 == atoms()->nr - qcustom.size(),
                           gmx::formatString("Number of custom charges %d does not match the number of atoms %d", static_cast<int>(qcustom.size()), atoms()->nr).c_str());
    }
    else if (algorithm == ChargeGenerationAlgorithm::NONE)
    {
        algorithm = pd->chargeGenerationAlgorithm();
        // Check whether there are free charges
        bool allFixed = true;
        auto myatoms = atoms();
        for (auto i = 0; i < myatoms->nr; i++)
        {
            auto atype = *myatoms->atomtype[i];
            auto ptype = pd->findParticleType(atype);
            auto qff = ptype->parameterConst("charge");
            if (qff.mutability() == Mutability::ACM)
            {
                allFixed = false;
            }
        }    
        if (allFixed)
        {
            algorithm = ChargeGenerationAlgorithm::NONE;
        }
    }
    switch (algorithm)
    {
    case ChargeGenerationAlgorithm::NONE:
        {
            if (debug)
            {
                fprintf(debug, "WARNING! Using fixed charges for %s!\n",
                        getMolname().c_str());
            }
            auto myatoms = atoms();
            for (auto i = 0; i < myatoms->nr; i++)
            {
                auto atype = *myatoms->atomtype[i];
                auto ptype = pd->findParticleType(atype);
                auto qval  = ptype->parameterConst("charge").value();
                myatoms->atom[i].q  = myatoms->atom[i].qB = qval;
            }
            // Now if we have shells, we still have to minimize them!
            if (nullptr != shellfc_)
            {
                double rmsf;
                auto imm = computeForces(nullptr, cr, &rmsf);
                if (imm != immStatus::OK)
                {
                    return imm;
                }
                auto qcalc = qTypeProps(qType::Calc);
                qcalc->setQ(myatoms);
                qcalc->setX(state_->x);
            }
            
            return immStatus::OK;
        }
    case ChargeGenerationAlgorithm::CM5:
    case ChargeGenerationAlgorithm::Hirshfeld:
    case ChargeGenerationAlgorithm::Mulliken:
        {
            auto myatoms = atoms();
            for (auto exper : experimentConst())
            {
                std::string mylot(exper.getMethod());
                mylot += "/" + exper.getBasisset();
                if (lot == mylot)
                {
                    std::string qalg = chargeGenerationAlgorithmName(algorithm);
                    int i = 0;
                    for (auto &ca : exper.calcAtomConst())
                    {
                        if (ca.hasCharge(qalg))
                        {
                            myatoms->atom[i].q  = myatoms->atom[i].qB = ca.charge(qalg);
                            i++;
                        }
                        else
                        {
                            gmx_fatal(FARGS, "No charge type %s for %s", qalg.c_str(), getMolname().c_str());
                        }
                    }
                }
            }
            return immStatus::OK;
        }
    case ChargeGenerationAlgorithm::Custom:
        {
            auto myatoms = atoms();
            for (auto i = 0; i < mtop_->natoms; i++)
            {
                myatoms->atom[i].q  = myatoms->atom[i].qB = qcustom[i];
            }
            return immStatus::OK;
        }
    case ChargeGenerationAlgorithm::ESP:
        {
            double chi2[2]   = {1e8, 1e8};
            int    cur       = 0;
            iter             = 0;
            
            // Init Qgresp should be called before this!
            auto qcalc   = qTypeProps(qType::Calc);
            GMX_RELEASE_ASSERT(qcalc != nullptr, "qType::Calc is not initialized");
            qcalc->qgenResp()->optimizeCharges(pd->getEpsilonR());
            qcalc->qgenResp()->calcPot(pd->getEpsilonR());
            qcalc->copyRespQ();
            
            if (debug)
            {
                fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
            }
            auto myatoms = atoms();
            do
            {
                auto qq = qcalc->charge();
                GMX_RELEASE_ASSERT(myatoms->nr-qq.size() == 0,
                                   gmx::formatString("Number of particles (%d) differs from number of charges (%d)", 
                                                     myatoms->nr, static_cast<int>(qq.size())).c_str());
                for (auto i = 0; i < myatoms->nr; i++)
                {
                    myatoms->atom[i].q = myatoms->atom[i].qB = qq[i];
                }
                if (nullptr != shellfc_)
                {
                    double rmsf;
                    auto imm = computeForces(nullptr, cr, &rmsf);
                    if (imm != immStatus::OK)
                    {
                        return imm;
                    }
                    qcalc->setX(state_->x);
                }
                qcalc->qgenResp()->optimizeCharges(pd->getEpsilonR());
                qcalc->qgenResp()->calcPot(pd->getEpsilonR());
                qcalc->copyRespQ();
                if (debug)
                {
                    fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
                }
                converged = (fabs(chi2[cur] - chi2[1-cur]) < tolerance) || (nullptr == shellfc_);
                cur       = 1-cur;
                iter++;
            }
            while ((!converged) && (iter < maxiter));
            auto qq = qcalc->charge();
            for (auto i = 0; i < myatoms->nr; i++)
            {
                myatoms->atom[i].q      =
                    myatoms->atom[i].qB = qq[i];
            }
        }
        break;
    case ChargeGenerationAlgorithm::EEM:
    case ChargeGenerationAlgorithm::SQE:
        {
            return GenerateAcmCharges(pd, cr, maxiter, tolerance);
        }
        break;
    }
    return imm;
}

void MyMol::changeCoordinate(const Experiment &ei, gmx_bool bpolar)
{
    const std::vector<gmx::RVec> &x = ei.getCoordinates();

    if (bpolar)
    {
        for (size_t i = 0; i < x.size(); i++)
        {
            copy_rvec(x[i], state_->x[2*i]);
            copy_rvec(x[i], state_->x[2*i+1]);
        }
    }
    else
    {
        for (size_t i = 0; i < x.size(); i++)
        {
            copy_rvec(x[i], state_->x[i]);
        }
    }
}

bool MyMol::getOptimizedGeometry(rvec *x)
{
    bool    bopt = false;

    for (auto &ei : experimentConst())
    {
        if (JobType::OPT == ei.getJobtype())
        {
            const std::vector<gmx::RVec> &xxx = ei.getCoordinates();
            for (size_t i = 0; i < xxx.size(); i++)
            {
                copy_rvec(xxx[i], x[i]);
            }
            bopt = true;
            break;
        }
    }
    return bopt;
}

void MyMol::CalcQPol(const Poldata *pd)
{
    int     np;
    double  poltot, sptot, ereftot;

    poltot  = 0;
    sptot   = 0;
    ereftot = 0;
    np      = 0;
    auto eep = pd->findForcesConst(InteractionType::POLARIZATION);
    auto myatoms = atomsConst();
    for (int i = 0; i < myatoms.nr; i++)
    {
        std::string ptype;
        auto atype = pd->findParticleType(*myatoms.atomtype[i]);
        auto idP   = atype->interactionTypeToIdentifier(InteractionType::POLARIZATION);
        if (!idP.id().empty() && eep.parameterExists(idP))
        {
            auto param  = eep.findParameterTypeConst(idP, "alpha");
            poltot += param.value();
            sptot  += gmx::square(param.uncertainty());
            np++;
        }
        ereftot += atype->refEnthalpy();
    }
    qProps_.find(qType::Calc)->second.calcMoments();
    ref_enthalpy_   = ereftot;
    if (np > 0)
    {
        polarizability_ = poltot;
        sig_pol_        = sqrt(sptot/np);
    }
}

void MyMol::CalcAnisoPolarizability(tensor polar, double *anisoPol)
{
    auto a = gmx::square(polar[XX][XX] - polar[YY][YY]);
    auto b = gmx::square(polar[XX][XX] - polar[ZZ][ZZ]);
    auto c = gmx::square(polar[ZZ][ZZ] - polar[YY][YY]);
    auto d = 6 * (gmx::square(polar[XX][YY]) + gmx::square(polar[XX][ZZ]) + gmx::square(polar[ZZ][YY]));

    *anisoPol = sqrt(1/2.0) * sqrt(a + b + c + d);
}

double MyMol::PolarizabilityTensorDeviation() const
{
    double delta2 = 0;
    for (auto i = 0; i < DIM; i++)
    {
        for (auto j = 0; j < DIM; j++)
        {
            delta2 += gmx::square(alpha_calc_[i][j] - alpha_elec_[i][j]);
        }
    }
    return delta2;
}

void MyMol::backupCoordinates()
{
    backupCoordinates_.resize(mtop_->natoms);
    for(int i = 0; i < mtop_->natoms; i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            backupCoordinates_[i][m] = state_->x[i][m];
        }
    }
}

void MyMol::restoreCoordinates()
{
    if (static_cast<int>(backupCoordinates_.size()) == mtop_->natoms)
    {
        for(int i = 0; i < mtop_->natoms; i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                state_->x[i][m] = backupCoordinates_[i][m];
            }
        }
    }
}

immStatus MyMol::CalcPolarizability(double     efield,
                                    t_commrec *cr,
                                    FILE      *fplog)
{
    const double        POLFAC = 29.957004; /* C.m**2.V*-1 to Å**3 */
    std::vector<double> field;
    rvec                mu_ref;
    immStatus           imm = immStatus::OK;
    double              rmsf;
    QtypeProps          qtp(qType::Calc);

    backupCoordinates();
    field.resize(DIM, 0);
    myforce_->setField(field);
    imm          = computeForces(fplog, cr, &rmsf);
    qtp.setQ(atoms());
    qtp.setX(state_->x);
    isoPol_calc_ = 0;
    qtp.calcMoments();
    copy_rvec(qtp.mu(), mu_ref);
    for (auto m = 0; imm == immStatus::OK && m < DIM; m++)
    {
        field[m] = efield;
        myforce_->setField(field);
        imm = computeForces(fplog, cr, &rmsf);
        qtp.setX(state_->x);
        field[m] = 0;
        myforce_->setField(field);
        if (imm == immStatus::OK)
        {
            rvec mu_tot;
            qtp.calcMoments();
            copy_rvec(qtp.mu(), mu_tot);
            for (auto n = 0; n < DIM; n++)
            {
                alpha_calc_[n][m] = ((mu_tot[n]-mu_ref[n])/efield)*(POLFAC);
            }
            isoPol_calc_ += alpha_calc_[m][m]/DIM;
        }
    }
    if (immStatus::OK == imm)
    {
        CalcAnisoPolarizability(alpha_calc_, &anisoPol_calc_);
    }
    restoreCoordinates();
    return imm;
}

void MyMol::PrintConformation(const char *fn)
{
    char title[STRLEN];

    put_in_box(mtop_->natoms, state_->box,
               as_rvec_array(state_->x.data()), 0.3);
    sprintf(title, "%s processed by alexandria", getMolname().c_str());
    write_sto_conf(fn, title, atoms(), as_rvec_array(state_->x.data()), nullptr, epbcNONE, state_->box);
}

void MyMol::PrintTopology(const char        *fn,
                          bool               bVerbose,
                          const Poldata     *pd,
                          t_commrec         *cr,
                          const std::string &method,
                          const std::string &basis)
{
    FILE  *fp   = gmx_ffopen(fn, "w");
    bool   bITP = (fn2ftp(fn) == efITP);

    PrintTopology(fp, bVerbose, pd, bITP, cr, method, basis);

    fclose(fp);
}

static void add_tensor(std::vector<std::string> *commercials,
                       const char *title, const tensor &Q)
{
    char buf[256];
    snprintf(buf, sizeof(buf), "%s:\n"
             "; ( %6.2f %6.2f %6.2f )\n"
             "; ( %6.2f %6.2f %6.2f )\n"
             "; ( %6.2f %6.2f %6.2f )\n",
             title,
             Q[XX][XX], Q[XX][YY], Q[XX][ZZ],
             Q[YY][XX], Q[YY][YY], Q[YY][ZZ],
             Q[ZZ][XX], Q[ZZ][YY], Q[ZZ][ZZ]);
    commercials->push_back(buf);
}

void MyMol::PrintTopology(FILE                   *fp,
                          bool                    bVerbose,
                          const Poldata          *pd,
                          bool                    bITP,
                          t_commrec              *cr,
                          const std::string      &method,
                          const std::string      &basis)
{
    char                     buf[256];
    t_mols                   printmol;
    std::vector<std::string> commercials;
    rvec                     vec;
    tensor                   myQ;
    double                   value = 0, error = 0, T = -1;
    std::string              myref;
    auto qt          = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));
    std::string              mylot       = makeLot(method, basis);

    if (fp == nullptr)
    {
        return;
    }

    CalcQPol(pd);
    if (getMolname().size() > 0)
    {
        printmol.name = strdup(getMolname().c_str());
    }
    else if (formula().size() > 0)
    {
        printmol.name = strdup(formula().c_str());
    }
    else
    {
        printmol.name = strdup("Unknown");
    }

    printmol.nr = 1;

    snprintf(buf, sizeof(buf), "Total Mass = %.3f (Da)", getMass());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Reference_Enthalpy = %.3f (kJ/mol)", ref_enthalpy_);
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Total Charge = %d (e)", totalCharge());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Charge Type  = %s\n",
             chargeTypeName(iChargeType).c_str());
    commercials.push_back(buf);
    
    auto qcalc = qTypeProps(qType::Calc);
    qcalc->setQ(atoms());
    qcalc->setX(state_->x);
    qcalc->calcMoments();
    snprintf(buf, sizeof(buf), "Alexandria Dipole Moment (Debye):\n"
             "; ( %.2f %6.2f %6.2f ) Total= %.2f\n",
             qcalc->mu()[XX], qcalc->mu()[YY], qcalc->mu()[ZZ],
             qcalc->dipole());
    commercials.push_back(buf);

    T = -1;
    const char *qm_type = "electronic";
    const char *qm_conf = "minimum";
    if (getPropRef(MolPropObservable::DIPOLE, iqmType::QM, method, basis, qm_conf,
                   qm_type, &value, &error,
                   &T, &myref, &mylot, vec, myQ))
    {
        qProps_.find(qType::Elec)->second.setMu(vec);
        snprintf(buf, sizeof(buf), "%s Dipole Moment (Debye):\n"
                 "; ( %.2f %6.2f %6.2f ) Total= %.2f\n",
                 mylot.c_str(),
                 vec[XX], vec[YY], vec[ZZ],
                 qProps_.find(qType::Elec)->second.dipole());
        commercials.push_back(buf);
    }
    else if (bVerbose)
    {
        printf("WARNING: QM dipole of type %s not found for lot %s\n",
               qm_type, mylot.c_str());
    }
    tensor aquad;
    copy_mat(qProps_.find(qType::Calc)->second.quad(), aquad);
    add_tensor(&commercials,
               "Alexandria Traceless Quadrupole Moments (Buckingham)",
               aquad);

    T = -1;
    if (getPropRef(MolPropObservable::QUADRUPOLE, iqmType::QM, method, basis, qm_conf,
                   qm_type, &value, &error,
                   &T, &myref, &mylot, vec, myQ))
    {
        qProps_.find(qType::Elec)->second.setQuadrupole(myQ);
        // TODO rotate_tensor(Q_qm_[qType::Elec], Q_qm_[qType::Calc]);
        snprintf(buf, sizeof(buf), "%s Traceless Quadrupole Moments (Buckingham)", mylot.c_str());
        add_tensor(&commercials, buf, myQ);
    }
    else if (bVerbose)
    {
        printf("WARNING: QM quadrupole of type %s not found for lot %s\n",
               qm_type, mylot.c_str());
    }

    // snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability (Additive Law): %.2f +/- %.2f (A^3)\n", polarizability_, sig_pol_);
    // commercials.push_back(buf);

    double efield = 0.1;
    if (nullptr != cr)
    {
        auto imm = CalcPolarizability(efield, cr, debug);
        if (imm == immStatus::OK)
        {
            add_tensor(&commercials, "Alexandria Polarizability components (A^3)", alpha_calc_);
            
            snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability: %.2f (A^3)\n", isoPol_calc_);
            commercials.push_back(buf);
            
            snprintf(buf, sizeof(buf), "Alexandria Anisotropic Polarizability: %.2f (A^3)\n", anisoPol_calc_);
            commercials.push_back(buf);

            T = -1;
            if (getPropRef(MolPropObservable::POLARIZABILITY, iqmType::QM, method, basis, "",
                                      (char *)"electronic", &isoPol_elec_, &error,
                                      &T, &myref, &mylot, vec, alpha_elec_))
            {
                CalcAnisoPolarizability(alpha_elec_, &anisoPol_elec_);
                snprintf(buf, sizeof(buf), "%s + Polarizability components (A^3)", mylot.c_str());
                add_tensor(&commercials, buf, alpha_elec_);
                
                snprintf(buf, sizeof(buf), "%s Isotropic Polarizability: %.2f (A^3)\n", mylot.c_str(), isoPol_elec_);
                commercials.push_back(buf);
                snprintf(buf, sizeof(buf), "%s Anisotropic Polarizability: %.2f (A^3)\n", mylot.c_str(), anisoPol_elec_);
                commercials.push_back(buf);
            }
        }
        else
        {
            commercials.push_back("Could not minimize shells. Cannot compute interactive polarizability.");
        }
    }

    print_top_header(fp, pd, bHaveShells_, commercials, bITP);
    write_top(fp, printmol.name, atoms(), false,
              plist_, excls_, gromppAtomtype_, cgnr_, pd);
    if (!bITP)
    {
        print_top_mols(fp, printmol.name, getForceField().c_str(), nullptr, 0, nullptr, 1, &printmol);
    }
    if (bVerbose)
    {
        for (auto &p : plist_)
        {
            if (p.nParam() > 0)
            {
                printf("There are %4d %s interactions\n", p.nParam(),
                       interaction_function[p.getFtype()].name);
            }
        }
        for (auto i = commercials.begin(); (i < commercials.end()); ++i)
        {
            printf("%s\n", i->c_str());
        }
    }

    sfree(printmol.name);
}

immStatus MyMol::GenerateChargeGroups(eChargeGroup ecg, bool bUsePDBcharge)
{
    real qtot, mtot;
    if ((cgnr_ = generate_charge_groups(ecg, atoms(), plist_, bUsePDBcharge,
                                        &qtot, &mtot)) == nullptr)
    {
        return immStatus::ChargeGeneration;
    }
    return immStatus::OK;
}

void MyMol::GenerateCube(const Poldata          *pd,
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
                         const gmx_output_env_t *oenv)
{
    auto qt          = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));

    if (potfn || hisfn || rhofn || difffn || pdbdifffn)
    {
        char     *gentop_version = (char *)"gentop v0.99b";
        auto qc = qProps_.find(qType::Calc);
        GMX_RELEASE_ASSERT(qc != qProps_.end(), "Cannot find alexandria charge information");
        qc->second.setQ(atoms());
        qc->second.setX(state_->x);
        qc->second.qgenResp()->calcPot(pd->getEpsilonR());
        qc->second.qgenResp()->potcomp(pcfn, atoms(),
                                       as_rvec_array(state_->x.data()),
                                       pdbdifffn, oenv);

        /* This has to be done before the grid is f*cked up by
           writing a cube file */
        QgenResp qCalc(qc->second.qgenResp());
        QgenResp grref(qc->second.qgenResp());

        if (reffn)
        {
            grref.setAtomInfo(atoms(), pd, state_->x, totalCharge());
            grref.setAtomSymmetry(symmetric_charges_);
            grref.readCube(reffn, FALSE);
        }
        else
        {
            qCalc.makeGrid(spacing, border, state_->x);
        }
        if (rhofn)
        {
            std::string buf = gmx::formatString("Electron density generated by %s based on %s charges",
                                                gentop_version,
                                                chargeTypeName(iChargeType).c_str());
            qCalc.calcRho();
            qCalc.writeRho(rhofn, buf, oenv);
        }
        if (potfn)
        {
            std::string buf = gmx::formatString("Potential generated by %s based on %s charges",
                                                gentop_version,
                                                chargeTypeName(iChargeType).c_str());
            qCalc.calcPot(pd->getEpsilonR());
            qCalc.writeCube(potfn, buf, oenv);
        }
        if (hisfn)
        {
            std::string buf = gmx::formatString("Potential generated by %s based on %s charges",
                                                gentop_version,
                                                chargeTypeName(iChargeType).c_str());
            qCalc.writeHisto(hisfn, buf, oenv);
        }
        if (difffn || diffhistfn)
        {
            std::string buf = gmx::formatString("Potential difference generated by %s based on %s charges",
                                                gentop_version,
                                                chargeTypeName(iChargeType).c_str());

            qCalc.writeDiffCube(&grref, difffn, diffhistfn, buf, oenv, 0);
        }
    }
}

void MyMol::rotateDipole(rvec mu, rvec muReference)
{
    if (norm(mu) < 0.04 or norm(muReference) < 0.04)
    {
        return;
    }
    return;
    matrix rotmatrix;
    rvec   tmpvec;
    calc_rotmatrix(mu, muReference, rotmatrix);
    mvmul(rotmatrix, mu, tmpvec);
    copy_rvec(tmpvec, mu);
}

void MyMol::calcEspRms(const Poldata *pd)
{
    int natoms = 0;
    for (int i = 0; i < atoms()->nr; i++)
    {
        if (atoms()->atom[i].ptype == eptAtom)
        {
            natoms++;
        }
    }
    t_atoms myatoms;
    gmx::HostVector<gmx::RVec> myx(natoms);
    init_t_atoms(&myatoms, natoms, FALSE);
    snew(myatoms.atomtype, natoms);
    natoms = 0;
    for (int i = 0; i < atoms()->nr; i++)
    {
        if (atoms()->atom[i].ptype == eptAtom)
        {
            myatoms.atom[natoms]     = atoms()->atom[i];
            myatoms.atomtype[natoms] = atoms()->atomtype[i];
            copy_rvec(x()[i], myx[natoms]);
            natoms++;
        }
    }
    
    auto qcalc   = qTypeProps(qType::Calc);
    auto qgrcalc = qcalc->qgenResp();
    for(auto &i : qProps_)
    {
        auto qi = i.first;
        if (qType::Calc == qi)
        {
            qgrcalc->setAtomInfo(atoms(), pd, x(), totalCharge());
            qgrcalc->updateAtomCharges(atoms());
            qgrcalc->calcPot(pd->getEpsilonR());
        }
        else if (qType::Elec != qi)
        {
            QgenResp *qgr = i.second.qgenResp();
            qgr->setChargeType(ChargeType::Point);
            qgr->setAtomInfo(&myatoms, pd, myx, totalCharge());
            qgr->updateAtomCharges(i.second.charge());
            for (const auto &ep : qgrcalc->espPoint())
            {
                auto r = ep.esp();
                qgr->addEspPoint(r[XX], r[YY], r[ZZ], ep.v());
            }
            qgr->calcPot(pd->getEpsilonR());
        }
    }
    done_atom(&myatoms);
}

immStatus MyMol::getExpProps(gmx_bool           bQM,
                             gmx_bool           bZero,
                             gmx_bool           bZPE,
                             gmx_bool           bDHform,
                             const std::string &method,
                             const std::string &basis,
                             const Poldata     *pd)
{
    int          ia    = 0;
    int          natom = 0;
    unsigned int nwarn = 0;
    double       value = 0;
    double       ZPE   = 0;
    double       error = 0;
    double       T     = -1;
    rvec         vec   = { 0 };
    tensor       quadrupole = {{ 0 }};
    tensor       polar      = {{ 0 }};
    std::string  myref;
    std::string  mylot;

    auto myatoms = atomsConst();
    gmx::HostVector<gmx::RVec> xatom(myatoms.nr);
    for (auto i = 0; i < myatoms.nr; i++)
    {
        if (myatoms.atom[i].ptype == eptAtom ||
            myatoms.atom[i].ptype == eptNucleus)
        {
            copy_rvec(state_->x[i], xatom[natom]);
            natom++;
        }
    }
    xatom.resizeWithPadding(natom);
    for(auto &i : qTypes())
    {
        std::vector<real> q;
        q.resize(natom, 0.0);
        qType qi = i.first;
        if (getPropRef(MolPropObservable::CHARGE, iqmType::QM,
                       method, basis, "",
                       qTypeName(qi),
                       &value, &error, &T,
                       &myref, &mylot, q.data(), quadrupole))
        {
            auto qp = qProps_.find(qi);
            if (qp == qProps_.end())
            {
                qProps_.insert(std::pair<qType, QtypeProps>(qi, QtypeProps(qi)));
                qp = qProps_.find(qi);
                GMX_RELEASE_ASSERT(qp != qProps_.end(), "Could not insert a new QtypeProps in qProps_");
            }
            qp->second.setQ(q);
            qp->second.setX(xatom);
            qp->second.calcMoments();
        }
    }
    T = 298.15;
    immStatus imm = immStatus::OK;
    if (bDHform &&
        getProp(MolPropObservable::ENERGY, (bQM ? iqmType::QM : iqmType::Exp),
                           method, basis, "",
                           (char *)"DeltaHform", &value, &error, &T))
    {
        Hform_ = value;
        Emol_  = value;
        for (ia = 0; ia < myatoms.nr; ia++)
        {
            if (myatoms.atom[ia].ptype == eptAtom ||
                myatoms.atom[ia].ptype == eptNucleus)
            {
                auto atype = pd->findParticleType(*myatoms.atomtype[ia]);
                Emol_ -= atype->refEnthalpy();
            }
        }
        if (bZPE)
        {

            if (getProp(MolPropObservable::ENERGY, iqmType::Both,
                                   method, basis, "",
                                   (char *)"ZPE", &ZPE, &error, &T))
            {
                Emol_ -= ZPE;
            }
            else
            {
                fprintf(stderr, "No zero-point energy for molecule %s.\n",
                        getMolname().c_str());
                imm = immStatus::NoData;
            }
        }
        if (ia < myatoms.nr)
        {
            imm = immStatus::NoData;
        }
    }
    if (imm == immStatus::OK)
    {
        T = -1;
        if (getPropRef(MolPropObservable::DIPOLE, iqmType::QM,
                       method, basis, "",
                       (char *)"electronic",
                       &value, &error, &T, &myref, &mylot,
                       vec, quadrupole))
        {
            dip_exp_  = value;
            dip_err_  = error;
            qProps_.find(qType::Elec)->second.setMu(vec);
            
            if (error <= 0)
            {
                if (debug)
                {
                    fprintf(debug, "WARNING: Error for %s is %g, assuming it is 10%%.\n",
                            getMolname().c_str(), error);
                }
                nwarn++;
                error = 0.1*value;
            }
            dip_weight_ = gmx::square(1.0/error);

            if (!bZero && fabs(qProps_.find(qType::Elec)->second.dipole()) < 1e-2)
            {
                imm = immStatus::ZeroDip;
            }
        }
        else
        {
            imm = immStatus::NoDipole;
        }
    }
    if (immStatus::OK == imm)
    {
        T = -1;
        if (getPropRef(MolPropObservable::QUADRUPOLE, iqmType::QM,
                       method, basis, "",
                       (char *)"electronic",
                       &value, &error, &T, &myref, &mylot,
                       vec, quadrupole))
        {
            qProps_.find(qType::Elec)->second.setQuadrupole(quadrupole);
        }
        T = -1;
        if (getPropRef(MolPropObservable::POLARIZABILITY, iqmType::QM,
                       method, basis, "",
                       (char *)"electronic",
                       &isoPol_elec_, &error, &T,
                       &myref, &mylot, vec, polar))
        {
            copy_mat(polar, alpha_elec_);
            CalcAnisoPolarizability(alpha_elec_, &anisoPol_elec_);
        }
    }
    return imm;
}

Identifier MyMol::getIdentifier(const Poldata                  *pd,
                                InteractionType                 iType,
                                const std::vector<std::string> &btype,
                                int                             natoms,
                                const int                      *iatoms)
{
    // A rather dirty hack.
    // TODO: make it less dirty.
    if (iType == InteractionType::VSITE2)
    {
        natoms -= 1;
        iatoms += 1;
    }
    std::vector<std::string> batoms;
    for (int j = 0; j < natoms; j++)
    {
        // Some atom types can be empty, e.g. polarizability of shells
        // since this a property of the atom.
        if (!btype[iatoms[j]].empty())
        {
            batoms.push_back(btype[iatoms[j]]);
        }
    }
    auto fs = pd->findForcesConst(iType);
    if (natoms == 2 && 
        (iType == InteractionType::BONDS ||
         iType == InteractionType::VSITE2))
    {
        auto bb = std::make_pair(iatoms[0], iatoms[1]);
        auto bo = bondOrder_.find(bb);
        if (bo == bondOrder_.end())
        {
            bb = std::make_pair(iatoms[1], iatoms[0]);
            bo = bondOrder_.find(bb);
        }
        if (bo != bondOrder_.end())
        {
            return Identifier(batoms, bo->second, fs.canSwap());
        }
        else
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Cannot find bond order for %s-%s",
                                                               btype[iatoms[0]].c_str(),
                                                               btype[iatoms[1]].c_str()).c_str()));
        }
    }
    return Identifier(batoms, fs.canSwap());
}

void MyMol::UpdateIdef(const Poldata   *pd,
                       InteractionType  iType)
{
    if (debug)
    {
        fprintf(debug, "UpdateIdef for %s\n", interactionTypeToString(iType).c_str());
    }
    if (iType == InteractionType::VDW)
    {
        nonbondedFromPdToMtop(mtop_, pd, fr_);
        if (debug)
        {
            pr_ffparams(debug, 0, "UpdateIdef Before", &mtop_->ffparams, false);
        }
    }
    else
    {
        // Update other iTypes
        auto fs    = pd->findForcesConst(iType);
        auto ftype = fs.fType();
        // Small optimization. This assumes angles etc. use the same
        // types as bonds. For this reason there is some special
        // treatment for polarization and vsites.
        std::vector<std::string> btype;
        auto myatoms = atomsConst();
        for(int j = 0; j < myatoms.nr; j++)
        {
            if (iType == InteractionType::VSITE2)
            {
                btype.push_back(*myatoms.atomtype[j]);
            }
            else
            {
                auto atype = pd->findParticleType(*myatoms.atomtype[j]);
                if (iType == InteractionType::POLARIZATION)
                {
                    btype.push_back(atype->interactionTypeToIdentifier(iType).id());
                }
                else
                {
                    btype.push_back(atype->interactionTypeToIdentifier(InteractionType::BONDS).id());
                }
            }
        }
        for (auto i = 0; i < ltop_->idef.il[ftype].nr; i += interaction_function[ftype].nratoms+1)
        {
            auto tp     = ltop_->idef.il[ftype].iatoms[i];
            auto bondId = getIdentifier(pd, iType, btype,
                                        interaction_function[ftype].nratoms,
                                        &ltop_->idef.il[ftype].iatoms[i+1]);
            
            switch (ftype)
            {
            case F_MORSE:
                {
                    auto fp = fs.findParameterTypeConst(bondId, "bondlength");
                    mtop_->ffparams.iparams[tp].morse.b0A         =
                        mtop_->ffparams.iparams[tp].morse.b0B     =
                        ltop_->idef.iparams[tp].morse.b0A     =
                        ltop_->idef.iparams[tp].morse.b0B =
                        convertToGromacs(fp.value(), fp.unit());
                    
                    fp = fs.findParameterTypeConst(bondId, "Dm");
                    mtop_->ffparams.iparams[tp].morse.cbA         =
                        mtop_->ffparams.iparams[tp].morse.cbB     =
                        ltop_->idef.iparams[tp].morse.cbA     =
                        ltop_->idef.iparams[tp].morse.cbB =
                        convertToGromacs(fp.value(), fp.unit());
                        
                    fp = fs.findParameterTypeConst(bondId, "beta");
                    mtop_->ffparams.iparams[tp].morse.betaA         =
                        mtop_->ffparams.iparams[tp].morse.betaB     =
                        ltop_->idef.iparams[tp].morse.betaA     =
                        ltop_->idef.iparams[tp].morse.betaB =
                        convertToGromacs(fp.value(), fp.unit());
                }
                break;
            case F_ANGLES:
                {
                    auto fp = fs.findParameterTypeConst(bondId, "angle");
                    mtop_->ffparams.iparams[tp].harmonic.rA         =
                        mtop_->ffparams.iparams[tp].harmonic.rB     =
                        ltop_->idef.iparams[tp].harmonic.rA     =
                        ltop_->idef.iparams[tp].harmonic.rB =
                        fp.value();
                        
                    fp = fs.findParameterTypeConst(bondId, "kt");
                    mtop_->ffparams.iparams[tp].harmonic.krA         =
                        mtop_->ffparams.iparams[tp].harmonic.krB     =
                        ltop_->idef.iparams[tp].harmonic.krA     =
                        ltop_->idef.iparams[tp].harmonic.krB =
                        convertToGromacs(fp.value(), fp.unit());
                }
                break;
            case F_POLARIZATION:
                {
                    auto fp = fs.findParameterTypeConst(bondId, "alpha");
                    mtop_->ffparams.iparams[tp].polarize.alpha =
                        ltop_->idef.iparams[tp].polarize.alpha =
                        convertToGromacs(fp.value(), fp.unit());
                }
                break;
            case F_UREY_BRADLEY:
                {
                    auto fp = fs.findParameterTypeConst(bondId, "angle");
                    mtop_->ffparams.iparams[tp].u_b.thetaA         =
                        mtop_->ffparams.iparams[tp].u_b.thetaB     =
                        ltop_->idef.iparams[tp].u_b.thetaA     =
                        ltop_->idef.iparams[tp].u_b.thetaB =
                        fp.value();
                        
                    fp = fs.findParameterTypeConst(bondId, "r13");
                    mtop_->ffparams.iparams[tp].u_b.r13A         =
                        mtop_->ffparams.iparams[tp].u_b.r13B     =
                        ltop_->idef.iparams[tp].u_b.r13A     =
                        ltop_->idef.iparams[tp].u_b.r13B =
                        convertToGromacs(fp.value(), fp.unit());
                        
                    fp = fs.findParameterTypeConst(bondId, "kt");
                    mtop_->ffparams.iparams[tp].u_b.kthetaA         =
                        mtop_->ffparams.iparams[tp].u_b.kthetaB     =
                        ltop_->idef.iparams[tp].u_b.kthetaA     =
                        ltop_->idef.iparams[tp].u_b.kthetaB =
                        convertToGromacs(fp.value(), fp.unit());
                        
                    fp = fs.findParameterTypeConst(bondId, "kub");
                    mtop_->ffparams.iparams[tp].u_b.kUBA         =
                        mtop_->ffparams.iparams[tp].u_b.kUBB     =
                        ltop_->idef.iparams[tp].u_b.kUBA     =
                        ltop_->idef.iparams[tp].u_b.kUBB =
                        convertToGromacs(fp.value(), fp.unit());
                }
                break;
            case F_LINEAR_ANGLES:
                {
                    // TODO: Check whether this is still needed!
                    //double relative_position = calc_relposition(pd, atoms[0], atoms[1], atoms[2]);
                    auto fp = fs.findParameterTypeConst(bondId, "a");
                    mtop_->ffparams.iparams[tp].linangle.aA         =
                        mtop_->ffparams.iparams[tp].linangle.aB     =
                        ltop_->idef.iparams[tp].linangle.aA     =
                        ltop_->idef.iparams[tp].linangle.aB =
                        convertToGromacs(fp.value(), fp.unit());
                    fp = fs.findParameterTypeConst(bondId, "klin");
                    mtop_->ffparams.iparams[tp].linangle.klinA         =
                        mtop_->ffparams.iparams[tp].linangle.klinB     =
                        ltop_->idef.iparams[tp].linangle.klinA     =
                        ltop_->idef.iparams[tp].linangle.klinB =
                        convertToGromacs(fp.value(), fp.unit());
                    // TODO: fix r13 in LINEAR_ANGLES                        
                    mtop_->ffparams.iparams[tp].linangle.r13A         =
                        mtop_->ffparams.iparams[tp].linangle.r13B     =
                        ltop_->idef.iparams[tp].linangle.r13A     =
                        ltop_->idef.iparams[tp].linangle.r13B = 0;
                        
                    mtop_->ffparams.iparams[tp].linangle.kUBA         =
                        mtop_->ffparams.iparams[tp].linangle.kUBB     =
                        ltop_->idef.iparams[tp].linangle.kUBA     =
                        ltop_->idef.iparams[tp].linangle.kUBB = 0;
                }
                break;
            case F_FOURDIHS:
                {
                    auto newparam = &mtop_->ffparams.iparams[tp];
                    std::vector<double> parameters = {
                        fs.findParameterTypeConst(bondId, "c0").value(),
                        fs.findParameterTypeConst(bondId, "c1").value(),
                        fs.findParameterTypeConst(bondId, "c2").value(),
                        fs.findParameterTypeConst(bondId, "c3").value()
                    };
                    newparam->rbdihs.rbcA[0] = parameters[1]+0.5*(parameters[0]+parameters[2]);
                    newparam->rbdihs.rbcA[1] = 0.5*(3.0*parameters[2]-parameters[0]);
                    newparam->rbdihs.rbcA[2] = 4.0*parameters[3]-parameters[1];
                    newparam->rbdihs.rbcA[3] = -2.0*parameters[2];
                    newparam->rbdihs.rbcA[4] = -4.0*parameters[3];
                    newparam->rbdihs.rbcA[5] = 0.0;
                    for (int k = 0; k < NR_RBDIHS; k++)
                    {
                        ltop_->idef.iparams[tp].rbdihs.rbcA[k]     =
                            ltop_->idef.iparams[tp].rbdihs.rbcB[k] =
                            newparam->rbdihs.rbcB[k]           =
                            newparam->rbdihs.rbcA[k];
                    }
                    break;
                }
            case F_PDIHS:
                {
                    auto fp = fs.findParameterTypeConst(bondId, "phi");
                    mtop_->ffparams.iparams[tp].pdihs.phiA         =
                        mtop_->ffparams.iparams[tp].pdihs.phiB     =
                        ltop_->idef.iparams[tp].pdihs.phiA     =
                        ltop_->idef.iparams[tp].pdihs.phiB =
                        convertToGromacs(fp.value(), fp.unit());
                        
                    fp = fs.findParameterTypeConst(bondId, "cp");
                    mtop_->ffparams.iparams[tp].pdihs.cpA         =
                        mtop_->ffparams.iparams[tp].pdihs.cpB     =
                        ltop_->idef.iparams[tp].pdihs.cpA     =
                        ltop_->idef.iparams[tp].pdihs.cpB =
                        convertToGromacs(fp.value(), fp.unit());
                    
                    int mult = fs.findParameterTypeConst(bondId, "mult").value();
                    mtop_->ffparams.iparams[tp].pdihs.mult =
                        ltop_->idef.iparams[tp].pdihs.mult = mult;
                }
                break;
            case F_IDIHS:
                {
                    auto fp = fs.findParameterTypeConst(bondId, "phi");
                    mtop_->ffparams.iparams[tp].harmonic.rA         =
                        mtop_->ffparams.iparams[tp].harmonic.rB     =
                        ltop_->idef.iparams[tp].harmonic.rA     =
                        ltop_->idef.iparams[tp].harmonic.rB =
                        convertToGromacs(fp.value(), fp.unit());
                        
                    fp = fs.findParameterTypeConst(bondId, "kimp");
                    mtop_->ffparams.iparams[tp].harmonic.krA         =
                        mtop_->ffparams.iparams[tp].harmonic.krB     =
                        ltop_->idef.iparams[tp].harmonic.krA     =
                        ltop_->idef.iparams[tp].harmonic.krB =
                        convertToGromacs(fp.value(), fp.unit());
                }
                break;
            case F_VSITE2:
                {
                    auto fp = fs.findParameterTypeConst(bondId, "v2_a");
                    mtop_->ffparams.iparams[tp].vsite.a         =
                        convertToGromacs(fp.value(), fp.unit());
                }
                break;
            default:
                GMX_THROW(gmx::InternalError(gmx::formatString("Do not know what to do for %s",
                                                               interaction_function[ftype].longname).c_str()));
                break;
            }
        }
    }
}

void MyMol::initQgenResp(const Poldata     *pd,
                         const std::string &method,
                         const std::string &basis,
                         real               watoms,
                         int                maxESP)
{
    std::string        mylot;
    auto qt          = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));
    auto qp          = qTypeProps(qType::Calc);
    QgenResp *qgr = qp->qgenResp();
    qgr->setChargeType(iChargeType);
    qgr->setAtomInfo(atoms(), pd, state_->x, totalCharge());
    qp->setQ(atoms());
    qp->setX(state_->x);
    qgr->setAtomSymmetry(symmetric_charges_);
    qgr->summary(debug);

    auto myatoms = atoms();
    int natoms = 0;
    for (int aa = 0; aa < myatoms->nr; aa++)
    {
        if (myatoms->atom[aa].ptype == eptAtom ||
            myatoms->atom[aa].ptype == eptNucleus)
        {
            natoms++;
        }
    }
    std::random_device               rd;
    std::mt19937                     gen(rd());  
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    double                           cutoff = 0.01*maxESP;
 
    auto ci = getCalcPropType(method, basis, &mylot, MolPropObservable::POTENTIAL, nullptr);
    if (ci != EndExperiment())
    {
        int iesp = 0;
        for (auto &epi : ci->electrostaticPotentialConst())
        {
            auto val = uniform(gen);
            if ((iesp >= natoms || watoms > 0) && val <= cutoff)
            {
                auto xu = epi.getXYZunit();
                auto vu = epi.getVunit();
                qgr->addEspPoint(convertToGromacs(epi.getX(), xu),
                                 convertToGromacs(epi.getY(), xu),
                                 convertToGromacs(epi.getZ(), xu),
                                 convertToGromacs(epi.getV(), vu));
            }
            iesp++;
        }
        if (debug)
        {
            fprintf(debug, "%s added %zu out of %zu ESP points to the RESP structure.\n",
                    getMolname().c_str(), qgr->nEsp(),
                    ci->electrostaticPotentialConst().size());
        }
    }
}

QtypeProps *MyMol::qTypeProps(qType qt)
{
    auto qp = qProps_.find(qt);
    if (qp != qProps_.end())
    {
        return &qp->second;
    }
    return nullptr;
}

const QtypeProps *MyMol::qTypeProps(qType qt) const
{
    const auto qp = qProps_.find(qt);
    if (qp != qProps_.end())
    {
        return &qp->second;
    }
    return nullptr;
}

void MyMol::plotEspCorrelation(const char             *espcorr,
                               const gmx_output_env_t *oenv,
                               t_commrec              *cr)
{
    if (espcorr && oenv)
    {
        auto qgr   = qTypeProps(qType::Calc)->qgenResp();
        qgr->updateAtomCharges(atoms());
        qgr->updateAtomCoords(state_->x);
        double rmsf = 0;
        if (immStatus::OK == computeForces(nullptr, cr, &rmsf))
        {
            qgr->calcPot(1.0);
            qgr->plotLsq(oenv, espcorr);
        }
    }
}

CommunicationStatus MyMol::Send(t_commrec *cr, int dest) const
{
    auto cs = MolProp::Send(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, static_cast<int>(dataset_type_));
    }
    return cs;
}

CommunicationStatus MyMol::Receive(t_commrec *cr, int src)
{
    auto cs = MolProp::Receive(cr, src);
    if (CS_OK == cs)
    {
        set_datasetType(static_cast<iMolSelect>(gmx_recv_int(cr, src)));
    }
    return cs;
}

} // namespace alexandria
