/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
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

#include "act/poldata/forcefieldparameter.h"
#include "gromacs_top.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/multipole_names.h"
#include "mymol_low.h"
#include "symmetrize_charges.h"
#include "act/utility/units.h"

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


MyMol::MyMol() //: gvt_(VsiteType::ALL)

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
        if (ca == bi.aJ() ||
            ca == bi.aI())
        {
            if (ca == bi.aI())
            {
                bca = bi.aJ();
                atoms.push_back(bca);
            }
            else
            {
                bca = bi.aI();
                atoms.push_back(bca);
            }
        }
    }
    /*Now try to find atoms bound to bca, except ca.*/
    for (auto bi : bondsConst())
    {
        if ((ca != bi.aJ()   &&
             ca != bi.aI())  &&
            (bca == bi.aJ()  ||
             bca == bi.aI()))
        {
            if (bca == bi.aI())
            {
                atoms.push_back(bi.aJ());
            }
            else
            {
                atoms.push_back(bi.aI());
            }
        }
    }
}

void MyMol::findOutPlaneAtoms(int ca, std::vector<int> &atoms)
{
    for (auto &bi : bondsConst())
    {
        if (bi.bondOrder() == 1  &&
            (ca == bi.aJ() ||
             ca == bi.aI()))
        {
            if (ca == bi.aI())
            {
                atoms.push_back(bi.aJ());
            }
            else
            {
                atoms.push_back(bi.aI());
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

immStatus MyMol::GenerateAtoms(const Poldata     *pd,
                               t_atoms           *atoms,
                               const std::string &method,
                               const std::string &basis,
                               bool               strict=true)
{
    double                    xx, yy, zz;
    int                       natom = 0;
    immStatus                 imm   = immStatus::OK;
    std::vector<std::string>  confs = { "minimum" };
    const Experiment         *ci    = nullptr;
    
    for(size_t iconf = 0; iconf < confs.size(); iconf++)
    { 
        ci = findExperimentConst(method, basis, confs[iconf]);
        if (ci && ci->NAtom() > 0)
        {
            break;
        }
    }
    if (!ci && !strict)
    {
        if (debug)
        {
            fprintf(debug, "Trying to find calculation without known method/basisset for %s\n", getMolname().c_str());
        }
        ci = findExperimentConst("", "", "");
    }
    if (ci)
    {
        if (ci->NAtom() == 0)
        {
            return immStatus::NoAtoms;
        }
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
                error_messages_.push_back(gmx::formatString("Cannot find atomtype %s (atom %d) in poldata, there are %d atomtypes.\n", 
                                                            cai.getObtype().c_str(), natom, pd->nParticleTypes()));

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
    int nmissing        = 0;
    int atomnumberTotal = 0;
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
        else
        {
            atomnumberTotal += pd->findParticleType(atype)->atomnumber();
        }
    }
    if (nmissing > 0)
    {
        return immStatus::AtomTypes;
    }
    // Check multiplicity
    int multOk = atomnumberTotal + getMultiplicity() + totalCharge();
    if (multOk % 2 == 0)
    {
        return immStatus::Multiplicity;
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

static void setTopologyIdentifiers(Topology      *top,
                                   const Poldata *pd,
                                   const t_atoms *myatoms)
{
    auto entries = top->entries(); 
    for(auto &entry : *entries)
    {
        auto fs = pd->findForcesConst(entry.first);
        for(auto &topentry : entry.second)
        {
            std::vector<std::string> btype;
            for(auto &jj : topentry->atomIndices())
            {
                auto atype = pd->findParticleType(*myatoms->atomtype[jj]);
                switch (entry.first)
                {
                case InteractionType::VSITE2:
                {
                    btype.push_back(*myatoms->atomtype[jj]);
                    break;
                }
                case InteractionType::POLARIZATION:
                {
                    btype.push_back(atype->interactionTypeToIdentifier(entry.first).id());
                    break;
                }
                default:
                {
                    btype.push_back(atype->interactionTypeToIdentifier(InteractionType::BONDS).id());
                    break;
                }
                }
            }
            topentry->setId(Identifier(btype, topentry->bondOrders(), fs.canSwap()));   
        }
    }
}

static void UpdateIdefEntry(const ForceFieldParameterList &fs,
                            const Identifier              &bondId,
                            int                            gromacsType,
                            gmx_mtop_t                    *mtop,
                            gmx_localtop_t                *ltop)
{
    double myval = 0;
    switch (fs.fType())
    {
    case F_MORSE:
        {
            auto fp = fs.findParameterTypeConst(bondId, "bondlength");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].morse.b0A         =
                mtop->ffparams.iparams[gromacsType].morse.b0B     = myval;
            if (ltop)
            { 
                ltop->idef.iparams[gromacsType].morse.b0A     =
                ltop->idef.iparams[gromacsType].morse.b0B = myval;
            }
                    
            fp = fs.findParameterTypeConst(bondId, "Dm");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].morse.cbA         =
                mtop->ffparams.iparams[gromacsType].morse.cbB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].morse.cbA     =
                ltop->idef.iparams[gromacsType].morse.cbB = myval;
            }   
                        
            fp = fs.findParameterTypeConst(bondId, "beta");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].morse.betaA         =
                mtop->ffparams.iparams[gromacsType].morse.betaB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].morse.betaA     =
                ltop->idef.iparams[gromacsType].morse.betaB = myval;
            }
        }
        break;
    case F_ANGLES:
        {
            auto fp = fs.findParameterTypeConst(bondId, "angle");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].harmonic.rA         =
                mtop->ffparams.iparams[gromacsType].harmonic.rB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].harmonic.rA     =
                ltop->idef.iparams[gromacsType].harmonic.rB = myval;
            }
                        
            fp = fs.findParameterTypeConst(bondId, "kt");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].harmonic.krA         =
                mtop->ffparams.iparams[gromacsType].harmonic.krB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].harmonic.krA     =
                ltop->idef.iparams[gromacsType].harmonic.krB = myval;
            }
        }
        break;
    case F_POLARIZATION:
        {
            auto fp = fs.findParameterTypeConst(bondId, "alpha");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].polarize.alpha = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].polarize.alpha = myval;
            }
        }
        break;
    case F_UREY_BRADLEY:
        {
            auto fp = fs.findParameterTypeConst(bondId, "angle");
            double angle = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].u_b.thetaA         =
                mtop->ffparams.iparams[gromacsType].u_b.thetaB     = angle;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].u_b.thetaA     =
                ltop->idef.iparams[gromacsType].u_b.thetaB = angle;
            }

            fp = fs.findParameterTypeConst(bondId, "kt");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].u_b.kthetaA         =
                mtop->ffparams.iparams[gromacsType].u_b.kthetaB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].u_b.kthetaA     =
                ltop->idef.iparams[gromacsType].u_b.kthetaB = myval;
            }

            fp = fs.findParameterTypeConst(bondId, "r13");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].u_b.r13A         =
                mtop->ffparams.iparams[gromacsType].u_b.r13B     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].u_b.r13A     =
                ltop->idef.iparams[gromacsType].u_b.r13B = myval;
            }                     
                         
            fp = fs.findParameterTypeConst(bondId, "kub");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].u_b.kUBA         =
                mtop->ffparams.iparams[gromacsType].u_b.kUBB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].u_b.kUBA     =
                ltop->idef.iparams[gromacsType].u_b.kUBB = myval;
            }
        }
        break;
    case F_LINEAR_ANGLES:
        {
            // TODO: Check whether this is still needed!
            //double relative_position = calc_relposition(pd, atoms[0], atoms[1], atoms[2]);
            auto fp = fs.findParameterTypeConst(bondId, "a");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].linangle.aA         =
                mtop->ffparams.iparams[gromacsType].linangle.aB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].linangle.aA     =
                ltop->idef.iparams[gromacsType].linangle.aB = myval;
            }
            
            fp = fs.findParameterTypeConst(bondId, "klin");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].linangle.klinA         =
                mtop->ffparams.iparams[gromacsType].linangle.klinB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].linangle.klinA     =
                ltop->idef.iparams[gromacsType].linangle.klinB = myval;
            }

            fp = fs.findParameterTypeConst(bondId, "r13lin");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].linangle.r13A         =
                mtop->ffparams.iparams[gromacsType].linangle.r13B     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].linangle.r13A     =
                ltop->idef.iparams[gromacsType].linangle.r13B = myval;
            }
                        
            fp = fs.findParameterTypeConst(bondId, "kublin");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].linangle.kUBA         =
                mtop->ffparams.iparams[gromacsType].linangle.kUBB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].linangle.kUBA     =
                ltop->idef.iparams[gromacsType].linangle.kUBB = myval;
            }
        }
        break;
    case F_FOURDIHS:
        {
            auto newparam = &mtop->ffparams.iparams[gromacsType];
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
                if (ltop)
                {
                    ltop->idef.iparams[gromacsType].rbdihs.rbcA[k]     =
                    ltop->idef.iparams[gromacsType].rbdihs.rbcB[k] = newparam->rbdihs.rbcA[k];
                }
                newparam->rbdihs.rbcB[k] = newparam->rbdihs.rbcA[k];
                    
            }
            break;
        }
    case F_PDIHS:
        {
            auto fp = fs.findParameterTypeConst(bondId, "phi");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].pdihs.phiA         =
                mtop->ffparams.iparams[gromacsType].pdihs.phiB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].pdihs.phiA     =
                ltop->idef.iparams[gromacsType].pdihs.phiB = myval;
            }
                        
            fp = fs.findParameterTypeConst(bondId, "cp");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].pdihs.cpA         =
                mtop->ffparams.iparams[gromacsType].pdihs.cpB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].pdihs.cpA     =
                ltop->idef.iparams[gromacsType].pdihs.cpB = myval;
            }
                    
            int mult = fs.findParameterTypeConst(bondId, "mult").value();
            mtop->ffparams.iparams[gromacsType].pdihs.mult = mult;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].pdihs.mult = mult;
            }
        }
        break;
    case F_IDIHS:
        {
            auto fp = fs.findParameterTypeConst(bondId, "phi");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].harmonic.rA         =
                mtop->ffparams.iparams[gromacsType].harmonic.rB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].harmonic.rA     =
                ltop->idef.iparams[gromacsType].harmonic.rB = myval;
            }
                        
            fp = fs.findParameterTypeConst(bondId, "kimp");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].harmonic.krA         =
                mtop->ffparams.iparams[gromacsType].harmonic.krB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].harmonic.krA     =
                ltop->idef.iparams[gromacsType].harmonic.krB = myval;
            }
        }
        break;
    case F_VSITE2:
        {
            auto fp = fs.findParameterTypeConst(bondId, "v2_a");
            myval = convertToGromacs(fp.value(), fp.unit());
            mtop->ffparams.iparams[gromacsType].vsite.a         = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].vsite.a = myval;
            }
        }
        break;
    default:
        GMX_THROW(gmx::InternalError(gmx::formatString("Do not know what to do for %s",
                                                        interaction_function[fs.fType()].longname).c_str()));
        break;
    }
}

static void TopologyToMtop(Topology       *top,
                           const Poldata  *pd,
                           gmx_mtop_t     *mtop)
{
    int ffparamsSize = mtop->ffparams.numTypes();
    auto entries = top->entries();
    for(auto &entry : *entries)
    {
        std::map<Identifier, int> idToGromacsType;
        auto fs = pd->findForcesConst(entry.first);
        for(auto &topentry: entry.second)
        {
            int gromacsType = -1;
            // First check whether we have this gromacsType already
            auto &bondId = topentry->id();
            GMX_RELEASE_ASSERT(!bondId.id().empty(), "Empty bondId");
            if (idToGromacsType.end() != idToGromacsType.find(bondId))
            {
                gromacsType = idToGromacsType.find(bondId)->second;
            }
            else
            {
                mtop->ffparams.functype.push_back(fs.fType());
                t_iparams ip = { { 0 } };
                mtop->ffparams.iparams.push_back(ip);
                gromacsType = mtop->ffparams.numTypes()-1;
                UpdateIdefEntry(fs, topentry->id(), gromacsType, mtop, nullptr);
                
            }
            if (gromacsType >= ffparamsSize)
            {
                topentry->setGromacsType(gromacsType);
                idToGromacsType.insert(std::pair<Identifier, int>(bondId, gromacsType));
            }
            else
            {
                GMX_THROW(gmx::InternalError("Could not add a force field parameter to the gromacs structure"));
            }
            // One more consistency check
            if (interaction_function[fs.fType()].nratoms !=
                static_cast<int>(topentry->atomIndices().size()))
            {
                GMX_THROW(gmx::InternalError(
                           gmx::formatString("Inconsistency in number of atoms. Expected %d, got %d for %s",
                           interaction_function[fs.fType()].nratoms, static_cast<int>(topentry->atomIndices().size()),
                           interaction_function[fs.fType()].name).c_str()));
            }
            // Now fill the gromacs structure
            mtop->moltype[0].ilist[fs.fType()].iatoms.push_back(gromacsType);
            for (auto &i : topentry->atomIndices())
            {
                mtop->moltype[0].ilist[fs.fType()].iatoms.push_back(i);
            }
        }
    }
}
                 
immStatus MyMol::GenerateTopology(FILE              *fp,
                                  const Poldata     *pd,
                                  const std::string &method,
                                  const std::string &basis,
                                  missingParameters  missing,
                                  const char        *tabfn,
                                  bool               strict)
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
        imm = GenerateAtoms(pd, atoms, method, basis, strict);
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
        topology_ = new Topology(bondsConst());
    }
    // Check whether we have dihedrals in the force field.
    bool bDih = false;
    if (pd->interactionPresent(InteractionType::PROPER_DIHEDRALS))
    {
        bDih = !pd->findForcesConst(InteractionType::PROPER_DIHEDRALS).empty();
    }
    
    if (immStatus::OK == imm)
    {
        topology_->makeAngles(state_->x, 175.0);
        topology_->makeImpropers(state_->x, 5.0);
        topology_->makePropers();
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
        else
        {
            snew(excls_, atoms->nr);
            excls_->nr = 0;
        }
        char **molnameptr = put_symtab(symtab_, getMolname().c_str());
        // Generate mtop
        mtop_ = do_init_mtop(pd, molnameptr, atoms,
                             inputrec_, symtab_, tabfn);
        // First create the identifiers for topology entries
        setTopologyIdentifiers(topology_, pd, atoms);
        // Now generate the mtop fields
        if (missing != missingParameters::Generate)
        {
            TopologyToMtop(topology_, pd, mtop_);
        }
        if (excls_)
        {
            excls_to_blocka(atoms->nr, excls_, &(mtop_->moltype[0].excls));
        }
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
        }
        if (nullptr == ltop_)
        {
            // Generate ltop from mtop
            ltop_ = gmx_mtop_generate_local_top(mtop_, false);
        }
    }
    if (immStatus::OK == imm && missing != missingParameters::Generate)
    {
        for(auto &entry : *(topology_->entries()))
        {
            UpdateIdef(pd, entry.first);
        }
    }
    if (immStatus::OK == imm && pd->polarizable())
    {
        shellfc_ = init_shell_flexcon(debug, mtop_, 0, 1, false);
    }
    if (immStatus::OK == imm)
    {
        imm = checkAtoms(pd, atoms);
    }
    if (immStatus::OK != imm && debug)
    {
        for(const auto &emsg : error_messages_)
        {
            fprintf(debug, "%s\n", emsg.c_str());
        }
    }
    if (immStatus::OK == imm)
    {
        qProps_.insert(std::pair<qType, QtypeProps>(qType::Calc, QtypeProps(qType::Calc)));
        qProps_.insert(std::pair<qType, QtypeProps>(qType::Elec, QtypeProps(qType::Elec)));
    }
    return imm;
}

void MyMol::addBondVsites(FILE          *fp,
                          const Poldata *pd,
                          t_atoms       *atoms)
{
    int     atomNrOld = atoms->nr;
    // First add virtual sites for bond shells if needed.
    auto    vs2       = pd->findForcesConst(InteractionType::VSITE2);
    auto    qt        = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    // TODO: add a flag to turn that on or off?
    t_atom  vsite_atom = { 0 };
    t_param nb         = { { 0 } };
    auto bonds = topology_->entry(InteractionType::BONDS);
    std::vector<TopologyEntry *> vsites;
    for (auto &b: bonds)
    {
        int ai = b->atomIndex(0);
        int aj = b->atomIndex(1);
        std::string aTypei(*atoms->atomtype[ai]);
        std::string aTypej(*atoms->atomtype[aj]);
        Identifier  bsId({ aTypei, aTypej }, b->bondOrders(), CanSwap::No);
        // Check whether a vsite is defined for this bond
        std::string v2("v2");
        if (vs2.parameterExists(bsId) && pd->hasParticleType(v2))
        {
            // Yes! We need to add a virtual site
            auto vsite = new TopologyEntry();
            vsite->addAtom(atoms->nr);
            vsite->addAtom(ai);
            vsite->addAtom(aj);
            vsite->addBondOrder(0.0);
            vsite->addBondOrder(0.0);
            vsites.push_back(vsite);
             // Now add the particle
            add_t_atoms(atoms, 1, 0);
            // Add exclusion for the vsite and its constituting atoms
            srenew(excls_, atoms->nr);
            excls_[atoms->nr-1].nr = 2;
            snew(excls_[atoms->nr-1].e, 2);
            excls_[atoms->nr-1].e[0] = ai;
            excls_[atoms->nr-1].e[1] = aj;
        
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
                                          ptype->id().id().c_str(), &nb, 0, 0);
            fill_atom(&atoms->atom[atoms->nr-1],
                      m, q, m, q,
                      ptype->atomnumber(),
                      vsgpp, eptVSite,
                      zeta, zeta, ptype->row(), 
                      atoms->atom[ai].resind);
            atoms->atomname[atoms->nr-1]  = put_symtab(symtab_, v2.c_str());
            atoms->atomtype[atoms->nr-1]  = put_symtab(symtab_, v2.c_str());
            atoms->atomtypeB[atoms->nr-1] = put_symtab(symtab_, v2.c_str());
        }
    }
    topology_->addEntry(InteractionType::VSITE2, vsites);
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
    std::vector<TopologyEntry *> pols;
    
    /* Calculate the total number of Atom and Vsite particles and
     * generate the renumbering array.
     */
    shellRenumber.resize(atoms->nr, 0);
    for (int i = 0; i < atoms->nr; i++)
    {
        auto atype                  = pd->findParticleType(*atoms->atomtype[i]);
        shellRenumber[i]            = i + nshell;
        originalAtomIndex_.insert(std::pair<int, int>(shellRenumber[i], i));
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
                    if (bHaveVSites_)
                    {
                        auto vsite = pd->findVsite(atomtype);
                        if (vsite != pd->getVsiteEnd())
                        {
                            pol /= vsite->nvsite();
                        }
                    }
                    auto pp = new TopologyEntry();
                    pp->addAtom(shellRenumber[i]);
                    pp->addAtom(shellRenumber[i]+1);
                    pp->addBondOrder(1.0);
                    pols.push_back(pp);
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
            // Add an exclusion for the shell
            snew(newexcls[shellRenumber[i]].e, 1);
            newexcls[shellRenumber[i]].e[0] = j;
            newexcls[shellRenumber[i]].nr = 1;
            auto atomtypeName = get_atomtype_name(atoms->atom[i].type, gromppAtomtype_);
            auto fa           = pd->findParticleType(atomtypeName);
            auto shellid      = fa->interactionTypeToIdentifier(InteractionType::POLARIZATION);
            if (shellid.id().empty())
            {
                // This particle has no shell.
                continue;
            }
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

    /* Copy exclusions, empty the original first */
    sfree(excls_);
    excls_ = newexcls;
    topology_->renumberAtoms(shellRenumber);
    topology_->addEntry(InteractionType::POLARIZATION, pols);
    bHaveShells_ = true;
}

double MyMol::bondOrder(int ai, int aj) const
{
    return topology_->findBond(ai, aj)->bondOrder();
}

immStatus MyMol::GenerateGromacs(const gmx::MDLogger       &mdlog,
                                 const CommunicationRecord *cr,
                                 const char                *tabfn,
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
    init_forcerec(nullptr, mdlog, fr_, nullptr, inputrec_, mtop_, cr->commrec(),
                  state_->box, tabfn, tabfn, tabbfnm, false, true, -1);
    gmx_omp_nthreads_set(emntBonded, 1);
    init_bonded_threading(nullptr, 1, &fr_->bondedThreading);
    setup_bonded_threading(fr_->bondedThreading, atomsConst().nr, false, ltop_->idef);
    wcycle_    = wallcycle_init(debug, 0, cr->commrec());

    MDatoms_  = new std::unique_ptr<gmx::MDAtoms>(new gmx::MDAtoms());
    *MDatoms_ = gmx::makeMDAtoms(nullptr, *mtop_, *inputrec_);
    atoms2md(mtop_, inputrec_, -1, nullptr, atomsConst().nr, MDatoms_->get());
    auto mdatoms = MDatoms_->get()->mdatoms();
    f_.resizeWithPadding(state_->natoms);

    if (nullptr != shellfc_)
    {
        make_local_shells(mdatoms, shellfc_);
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

real MyMol::potentialEnergy() const 
{
    return enerd_->term[F_EPOT];
}

immStatus MyMol::computeForces(FILE *fplog, double *rmsf)
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
    t_commrec *crtmp = init_commrec();
    crtmp->nnodes = 1;
    crtmp->nodeid = 0;
    if (!vsite_)
    {
        vsite_  = new std::unique_ptr<gmx_vsite_t>(new gmx_vsite_t());
        *vsite_ = initVsite(*mtop_, crtmp);
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
            fprintf(debug, "mol %s alpha %g\n", 
                    getMolname().c_str(),
                    mtop_->ffparams.iparams[mtop_->moltype[0].ilist[F_POLARIZATION].iatoms[0]].polarize.alpha);
        }
        real force2 = 0;
        try
        {
            force2 = relax_shell_flexcon(nullptr, crtmp, nullptr, false,
                                         0, inputrec_,
                                         true, force_flags, ltop_,
                                         enerd_, fcd_, state_,
                                         f_.arrayRefWithPadding(), force_vir, mdatoms,
                                         &nrnb_, wcycle_, nullptr,
                                         &(mtop_->groups), shellfc_,
                                         fr_, t, mu_tot, vsite_->get());
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
    }
    else
    {
        do_force(fplog, crtmp, nullptr, inputrec_, 0,
                 &nrnb_, wcycle_, ltop_,
                 &(mtop_->groups),
                 state_->box, state_->x.arrayRefWithPadding(), nullptr,
                 f_.arrayRefWithPadding(), force_vir, mdatoms,
                 enerd_, fcd_,
                 state_->lambda, nullptr,
                 fr_, vsite_->get(), mu_tot, t,
                 force_flags);
        *rmsf = 0;
    }
    done_commrec(crtmp);
    return imm;
}

void MyMol::symmetrizeCharges(const Poldata  *pd,
                              bool            bSymmetricCharges,
                              const char     *symm_string)
{
    if (bSymmetricCharges)
    {
        symmetric_charges_.clear();
        symmetrize_charges(bSymmetricCharges, atoms(), topology_->entry(InteractionType::BONDS),
                           pd, symm_string, &symmetric_charges_);
    }
    else
    {
        for (auto i = 0; i < mtop_->natoms; i++)
        {
            symmetric_charges_.push_back(i);
        }
    }
}

immStatus MyMol::GenerateAcmCharges(const Poldata             *pd,
                                    const CommunicationRecord *cr,
                                    int                        maxiter,
                                    real                       tolerance)
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
            if (haveShells())
            {
                double rmsf;
                auto imm = computeForces(nullptr, &rmsf);
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
            converged = (EemRms < tolerance) || !haveShells();
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
                                 const CommunicationRecord *cr,
                                 const char                *tabfn,
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

    GenerateGromacs(mdlog, cr, tabfn, iChargeType);
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
                auto imm = computeForces(nullptr, &rmsf);
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
            std::map<ChargeGenerationAlgorithm, qType> qtmap = {
                { ChargeGenerationAlgorithm::CM5, qType::CM5 },
                { ChargeGenerationAlgorithm::Hirshfeld, qType::Hirshfeld },
                { ChargeGenerationAlgorithm::Mulliken, qType::Mulliken }
            };
            auto myatoms = atoms();
            for (auto exper : experimentConst())
            {
                std::string mylot(exper.getMethod());
                mylot += "/" + exper.getBasisset();
                if (lot == mylot)
                {
                    int i = 0;
                    for (auto &ca : exper.calcAtomConst())
                    {
                        if (ca.hasCharge(qtmap[algorithm]))
                        {
                            myatoms->atom[i].q  = myatoms->atom[i].qB = ca.charge(qtmap[algorithm]);
                            i++;
                        }
                        else
                        {
                            gmx_fatal(FARGS, "No charge type %s for %s",
                                      qTypeName(qtmap[algorithm]).c_str(), getMolname().c_str());
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
                    auto imm = computeForces(nullptr, &rmsf);
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

immStatus MyMol::CalcPolarizability(double                     efield,
                                    const CommunicationRecord *cr,
                                    FILE                      *fplog)
{
    const double        POLFAC = 29.957004; /* C.m**2.V*-1 to Å**3 */
    std::vector<double> field;
    immStatus           imm = immStatus::OK;
    double              rmsf;
    QtypeProps          qtp(qType::Calc);

    backupCoordinates();
    field.resize(DIM, 0);
    myforce_->setField(field);
    imm          = computeForces(fplog, &rmsf);
    qtp.setQ(atoms());
    qtp.setX(state_->x);
    isoPol_calc_ = 0;
    qtp.calcMoments();
    auto mpo = MolPropObservable::DIPOLE;
    if (!qtp.hasMultipole(mpo))
    {
        GMX_THROW(gmx::InternalError("No dipole to compute."));
    }
    auto mu_ref = qtp.getMultipole(mpo);
    for (auto m = 0; imm == immStatus::OK && m < DIM; m++)
    {
        field[m] = efield;
        myforce_->setField(field);
        imm = computeForces(fplog, &rmsf);
        qtp.setX(state_->x);
        field[m] = 0;
        myforce_->setField(field);
        if (imm == immStatus::OK)
        {
            qtp.calcMoments();
            auto qmu = qtp.getMultipole(mpo);
            for (auto n = 0; n < DIM; n++)
            {
                alpha_calc_[n][m] = ((qmu[n]-mu_ref[n])/efield)*(POLFAC);
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

void MyMol::PrintTopology(const char                *fn,
                          bool                       bVerbose,
                          const Poldata             *pd,
                          const CommunicationRecord *cr,
                          const std::string         &method,
                          const std::string         &basis)
{
    FILE  *fp   = gmx_ffopen(fn, "w");
    bool   bITP = (fn2ftp(fn) == efITP);

    PrintTopology(fp, bVerbose, pd, bITP, cr, method, basis);

    fclose(fp);
}

static void add_tensor(std::vector<std::string> *commercials,
                       const char *title, const std::vector<double> &Q)
{
    char buf[256];
    snprintf(buf, sizeof(buf), "%s:\n"
             "; ( %6.2f %6.2f %6.2f )\n"
             "; (       %6.2f %6.2f )\n"
             "; (             %6.2f )\n",
             title,
             Q[0], Q[1], Q[2], Q[3], Q[4], Q[5]);
    commercials->push_back(buf);
}

void MyMol::PrintTopology(FILE                      *fp,
                          bool                       bVerbose,
                          const Poldata             *pd,
                          bool                       bITP,
                          const CommunicationRecord *cr,
                          const std::string         &method,
                          const std::string         &basis)
{
    char                     buf[256];
    t_mols                   printmol;
    std::vector<std::string> commercials;
    std::vector<double>      vec;
    double                   T = -1;
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
    auto qelec = qTypeProps(qType::Elec);
    qcalc->setQ(atoms());
    qcalc->setX(state_->x);
    qcalc->calcMoments();
    
    T = -1;
    const char *qm_type = "electronic";
    const char *qm_conf = "minimum";
    for(auto &mpo : mpoMultiPoles)
    {
        auto gp = findProperty(mpo, iqmType::QM, T, method, basis, qm_conf);
        if (gp)
        {
            auto vec = gp->getVector();
            qelec->setMultipole(mpo, vec);
            auto mymu = qelec->getMultipole(mpo);
            commercials.push_back(gmx::formatString("%s %s (%s)\n",
                                                    mylot.c_str(), mpo_name(mpo), mpo_unit(mpo)));
            for(auto &fmp : formatMultipole(mpo, mymu))
            {
                commercials.push_back(fmp);
            }
        }
        if (qcalc->hasMultipole(mpo))
        {
            auto mymu = qcalc->getMultipole(mpo);
            commercials.push_back(gmx::formatString("Alexandria %s (%s)\n", mpo_name(mpo), mpo_unit(mpo)));
            for(auto &fmp : formatMultipole(mpo, mymu))
            {
                commercials.push_back(fmp);
            }
        }
    }

    double efield = 0.1;
    if (nullptr != cr)
    {
        auto imm = CalcPolarizability(efield, cr, debug);
        if (imm == immStatus::OK)
        {
            std::vector<double> ac = { alpha_calc_[XX][XX], alpha_calc_[XX][YY], alpha_calc_[XX][ZZ],
                alpha_calc_[YY][YY], alpha_calc_[YY][ZZ], alpha_calc_[ZZ][ZZ] };
            add_tensor(&commercials, "Alexandria Polarizability components (A^3)", ac);
            
            snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability: %.2f (A^3)\n", isoPol_calc_);
            commercials.push_back(buf);
            
            snprintf(buf, sizeof(buf), "Alexandria Anisotropic Polarizability: %.2f (A^3)\n", anisoPol_calc_);
            commercials.push_back(buf);

            T = -1;
            auto gp = findProperty(MolPropObservable::POLARIZABILITY, iqmType::QM, T, method, basis, "");
            if (gp)
            {
                //&isoPol_elec_, &error, &T, &myref, &mylot, &vec, alpha_elec_))
                auto alpha_elec = gp->getVector();
                CalcAnisoPolarizability(alpha_elec_, &anisoPol_elec_);
                snprintf(buf, sizeof(buf), "%s + Polarizability components (A^3)", mylot.c_str());
                add_tensor(&commercials, buf, alpha_elec);
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

    // TODO write a replacement for this function
    print_top_header(fp, pd, bHaveShells_, commercials, bITP);
    write_top(fp, printmol.name, atoms(), false,
              topology_, excls_, gromppAtomtype_, pd);
    if (!bITP)
    {
        print_top_mols(fp, printmol.name, getForceField().c_str(), nullptr, 0, nullptr, 1, &printmol);
    }
    if (bVerbose)
    {
        for (auto &entry : *(topology_->entries()))
        {
            int ftype = pd->findForcesConst(entry.first).fType();
            if (entry.second.size() > 0)
            {
                printf("There are %4d %s interactions\n", static_cast<int>(entry.second.size()),
                       interaction_function[ftype].name);
            }
        }
        for (auto i = commercials.begin(); (i < commercials.end()); ++i)
        {
            printf("%s\n", i->c_str());
        }
    }

    sfree(printmol.name);
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
        QgenResp qCalc(*qc->second.qgenResp());
        QgenResp grref(*qc->second.qgenResp());

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

const real *MyMol::energyTerms() const
{
    return enerd_->term;
}
        
immStatus MyMol::getExpProps(const std::map<MolPropObservable, iqmType> &iqm,
                             gmx_bool                                    bZero,
                             const std::string                          &method,
                             const std::string                          &basis,
                             const Poldata                              *pd)
{
    int                 ia    = 0;
    int                 natom = 0;
    unsigned int        nwarn = 0;
    std::vector<double> vec;
    std::string         myref;
    std::string         mylot;
    immStatus           imm        = immStatus::OK;
    
    auto myatoms = atomsConst();
    GMX_RELEASE_ASSERT(myatoms.nr > 0, "No atoms!");
    
    // Make a copy of the coordinates without shells
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
    
    for (const auto &miq : iqm)
    {
        auto mpo = miq.first;
        switch (mpo)
        {
        case MolPropObservable::CHARGE:
            {
                std::string conf;
                auto ei = findExperimentConst(method, basis, conf);
                if (ei)
                {
                    for(auto &i : qTypes())
                    {
                        qType               qi = i.first;
                        std::string         reference;
                        std::string         lot;
                        std::vector<double> q;
                        if (ei->getCharges(&q, qi, &reference, &lot))
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
                }
                else
                {
                    imm = immStatus::NoData;
                }
            }
            break;
        case MolPropObservable::DHFORM:
        case MolPropObservable::DGFORM:
        case MolPropObservable::ZPE:
            {
                double    T   = 298.15;
                auto gp = static_cast<const MolecularEnergy *>(findProperty(mpo, miq.second, T, method, basis, ""));
                if (gp)
                {
                    energy_.insert(std::pair<MolPropObservable, double>(mpo, gp->getValue()));
                }
                else
                {
                    imm = immStatus::NoData;
                }
            }
            break;
        case MolPropObservable::DIPOLE:
        case MolPropObservable::QUADRUPOLE:
        case MolPropObservable::OCTUPOLE:
        case MolPropObservable::HEXADECAPOLE:
            {
                double T = -1;
                auto gp = static_cast<const MolecularMultipole *>(findProperty(mpo, miq.second, T, method, basis, ""));
                if (gp)
                {
                    qProps_.find(qType::Elec)->second.setMultipole(mpo, gp->getVector());
                }
                else
                {
                    imm = immStatus::NoData;
                }
            }
            break;
        case MolPropObservable::POLARIZABILITY:
            {
                double T = -1;
                auto gp = static_cast<const MolecularPolarizability *>(findProperty(mpo, miq.second, T, method, basis, ""));
                if (gp)
                {
                    copy_mat(gp->getTensor(), alpha_elec_);
                    CalcAnisoPolarizability(alpha_elec_, &anisoPol_elec_);
                }
                else
                {
                    imm = immStatus::NoData;
                }
            }
            break;
        default:
            break;
        }
    }

    if (immStatus::OK == imm &&
        energy_.find(MolPropObservable::DHFORM) != energy_.end())
    {
        double Emol = energy_[MolPropObservable::DHFORM];
        for (ia = 0; ia < myatoms.nr; ia++)
        {
            if (myatoms.atom[ia].ptype == eptAtom ||
                myatoms.atom[ia].ptype == eptNucleus)
            {
                auto atype = pd->findParticleType(*myatoms.atomtype[ia]);
                Emol -= atype->refEnthalpy();
            }
        }
        energy_.insert(std::pair<MolPropObservable, double>(MolPropObservable::EMOL, Emol));
        if (energy_.find(MolPropObservable::ZPE) != energy_.end())
        {
            energy_[MolPropObservable::EMOL] -= energy_[MolPropObservable::ZPE];
        }
    }
    return imm;
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
        if (topology_->hasEntry(iType))
        {
            auto fs    = pd->findForcesConst(iType);
            auto entry = topology_->entry(iType);
            for (size_t i = 0; i < entry.size(); i++)
            {
                auto &bondId = entry[i]->id();
                auto  gromacsType     = entry[i]->gromacsType();
                if (gromacsType < 0 || gromacsType >= mtop_->ffparams.numTypes())
                {
                    GMX_THROW(gmx::InternalError(gmx::formatString("gromacsType = %d should be >= 0 and < %d for %s %s", gromacsType, mtop_->ffparams.numTypes(), interactionTypeToString(iType).c_str(), bondId.id().c_str()).c_str()));
                }
                UpdateIdefEntry(fs, bondId, gromacsType, mtop_, ltop_);
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
 
    auto ci = findExperimentConst(method, basis, "");
    if (ci)
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

void MyMol::plotEspCorrelation(const char                *espcorr,
                               const gmx_output_env_t    *oenv,
                               const CommunicationRecord *cr)
{
    if (espcorr && oenv)
    {
        auto qgr   = qTypeProps(qType::Calc)->qgenResp();
        qgr->updateAtomCharges(atoms());
        qgr->updateAtomCoords(state_->x);
        double rmsf = 0;
        if (immStatus::OK == computeForces(nullptr, &rmsf))
        {
            qgr->calcPot(1.0);
            qgr->plotLsq(oenv, espcorr);
        }
    }
}

CommunicationStatus MyMol::Send(const CommunicationRecord *cr, int dest) const
{
    auto cs = MolProp::Send(cr, dest);
    if (CommunicationStatus::OK == cs)
    {
        cr->send_int(dest, static_cast<int>(dataset_type_));
    }
    return cs;
}

CommunicationStatus MyMol::Receive(const CommunicationRecord *cr, int src)
{
    auto cs = MolProp::Receive(cr, src);
    if (CommunicationStatus::OK == cs)
    {
        set_datasetType(static_cast<iMolSelect>(cr->recv_int(src)));
    }
    return cs;
}

} // namespace alexandria
