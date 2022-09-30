/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#include "mymol.h"

#include <cstdio>

#include <map>
#include <random>
#include <string>

#include "act/molprop/molprop_util.h"
#include "act/molprop/multipole_names.h"
#include "act/poldata/forcefieldparameter.h"
#include "act/poldata/forcefieldparametername.h"
#include "act/utility/regression.h"
#include "act/utility/units.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/paddedvector.h"
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
#include "gromacs_top.h"
#include "mymol_low.h"
#include "symmetrize_charges.h"

namespace alexandria
{

static void vsiteType_to_atomType(const std::string &vsiteType, std::string *atomType)
{
    std::size_t pos = vsiteType.find("L");
    *atomType       = vsiteType.substr (0, pos);
}

MyMol::MyMol() //: gvt_(VsiteType::ALL)

{
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

bool MyMol::IsSymmetric(real toler) const
{
    real  tm;
    rvec  com, test;
    bool  bSymmAll;
    std::vector<bool> bSymm;
    
    clear_rvec(com);
    tm = 0;
    const auto &myatoms = topology_->atoms();
    std::vector<gmx::RVec> myx = optimizedCoordinates_;
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        real mm  = myatoms[i].mass();
        tm += mm;
        for (int m = 0; (m < DIM); m++)
        {
            com[m] += mm*myx[i][m];
        }
    }
    if (tm > 0)
    {
        for (int m = 0; m < DIM; m++)
        {
            com[m] /= tm;
        }
    }
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        rvec_dec(myx[i], com);
    }

    bSymm.resize(myatoms.size());
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        bSymm[i] = (norm(myx[i]) < toler);
        for (size_t j = i+1; (j < myatoms.size()) && !bSymm[i]; j++)
        {
            rvec_add(myx[i], myx[j], test);
            if (norm(test) < toler)
            {
                bSymm[i] = true;
                bSymm[j] = true;
            }
        }
    }
    bSymmAll = true;
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        bSymmAll = bSymmAll && bSymm[i];
    }

    return bSymmAll;
}

void MyMol::findInPlaneAtoms(int ca, std::vector<int> *atoms) const
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
                atoms->push_back(bca);
            }
            else
            {
                bca = bi.aI();
                atoms->push_back(bca);
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
                atoms->push_back(bi.aJ());
            }
            else
            {
                atoms->push_back(bi.aI());
            }
        }
    }
}

void MyMol::findOutPlaneAtoms(int ca, std::vector<int> *atoms) const
{
    for (auto &bi : bondsConst())
    {
        if (bi.bondOrder() == 1  &&
            (ca == bi.aJ() ||
             ca == bi.aI()))
        {
            if (ca == bi.aI())
            {
                atoms->push_back(bi.aJ());
            }
            else
            {
                atoms->push_back(bi.aI());
            }
        }
    }
}

bool MyMol::IsVsiteNeeded(std::string    atype,
                          const Poldata *pd) const
{
    auto vsite = pd->findVsite(atype);
    return vsite != pd->getVsiteEnd();
}

immStatus MyMol::GenerateAtoms(const Poldata     *pd,
                               t_atoms           *atoms)
{
    double                    xx, yy, zz;
    int                       natom = 0;
    immStatus                 imm   = immStatus::OK;
    const Experiment         *ci    = nullptr;
    
    ci = findExperimentConst(JobType::OPT);
    if (!ci)
    {
        ci = findExperimentConst(JobType::TOPOLOGY);
    }
    if (ci)
    {
        if (ci->NAtom() == 0)
        {
            return immStatus::NoAtoms;
        }
        myJobType_ = ci->getJobtype();
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
            auto myunit = cai.coordUnit();
            cai.coords(&xx, &yy, &zz);
            int resnr = cai.ResidueNumber();
            if (resnr != res0)
            {
                res0  = resnr;
                atoms->nres += 1;
            }
            gmx::RVec xxx = { convertToGromacs(xx, myunit),
                              convertToGromacs(yy, myunit),
                              convertToGromacs(zz, myunit) };
            copy_rvec(xxx, state_->x[natom]);
            optimizedCoordinates_.push_back(xxx);
            
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
                atoms->atom[natom].q      =
                    atoms->atom[natom].qB = atype->charge();
                
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
        GMX_RELEASE_ASSERT(atoms->nr == natom, "Inconsistency numbering atoms");
    }
    else
    {
        imm = immStatus::Topology;
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Tried to convert %s to gromacs. LOT is %s/%s. Natoms is %d\n",
                getMolname().c_str(),
                ci->getMethod().c_str(),
                ci->getBasisset().c_str(), natom);
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
    int multOk = atomnumberTotal + totalMultiplicity() + totalCharge();
    if (multOk % 2 == 0)
    {
        fprintf(stderr, "atomnumberTotal %d, totalMultiplicity %d, totalCharge %d\n",
                atomnumberTotal, totalMultiplicity(), totalCharge());
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
    auto &qt       = pd->findForcesConst(InteractionType::COULOMB);
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
    
    auto iType = InteractionType::COULOMB;
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

void MyMol::forceEnergyMaps(const Poldata                                           *pd,
                            const ForceComputer                                     *forceComp,
                            std::vector<std::vector<std::pair<double, double> > >   *forceMap,
                            std::vector<std::pair<double, double> >                 *enerMap,
                            std::vector<std::pair<double, std::map<InteractionType, double> > > *enerAllMap)
{
    auto       myatoms = topology_->atoms();
    t_commrec *crtmp   = init_commrec();
    crtmp->nnodes      = 1;
    forceMap->clear();
    enerMap->clear();
    enerAllMap->clear();
    for (auto &ei : experimentConst())
    {
        // TODO: no need to recompute the energy if we just have
        // done that. Check for OPT being the first calculation.
        const std::vector<gmx::RVec> &xxx = ei.getCoordinates();
        std::vector<gmx::RVec> coords(myatoms.size());
        int j = 0;
        for(size_t i = 0; i < myatoms.size(); i++)
        {
            if (myatoms[i].pType() == eptAtom)
            {
                copy_rvec(xxx[j], coords[i]);
                copy_rvec(xxx[j], state_->x[i]);
                j += 1;
            }
            // TODO initiate shells?
        }
        std::map<InteractionType, double> energies;
        std::vector<gmx::RVec> forces(myatoms.size());
        if (forceComp)
        {
            (void) forceComp->compute(pd, topology_, &coords, &forces, &energies);
        }
        else
        {
            real shellForceRMS;
            PaddedVector<gmx::RVec> gmxforces;
            gmxforces.resizeWithPadding(myatoms.size());
            calculateEnergyOld(crtmp, &coords, &gmxforces, &energies, &shellForceRMS);
            for(size_t i = 0; i < myatoms.size(); i++)
            {
                copy_rvec(gmxforces[i], forces[i]);
            }
        }
        if (ei.hasProperty(MolPropObservable::INTERACTIONENERGY))
        {
            auto eprops = ei.propertyConst(MolPropObservable::INTERACTIONENERGY);
            if (eprops.size() > 1)
            {
                gmx_fatal(FARGS, "Multiple interaction energies for this experiment");
            }
            else if (eprops.size() == 1)
            {
                if (forceComp)
                {
                    double Einter = calculateInteractionEnergy(pd, forceComp);
                    
                    enerMap->push_back({ eprops[0]->getValue(), Einter });
                }
            }
        }
        else if (ei.hasProperty(MolPropObservable::DELTAE0))
        {
            auto eprops = ei.propertyConst(MolPropObservable::DELTAE0);
            if (eprops.size() > 1)
            {
                gmx_fatal(FARGS, "Multiple energies for this experiment");
            }
            else if (eprops.size() == 1)
            {
                enerMap->push_back({ eprops[0]->getValue(), energies[InteractionType::EPOT] });
                enerAllMap->push_back({ eprops[0]->getValue(), std::move(energies) });
            }
        }
        const std::vector<gmx::RVec> &fff = ei.getForces();
        if (!fff.empty())
        {
            size_t ifff = 0;
            std::vector<std::pair<double, double> > thisForce;
            for (size_t i = 0; i < myatoms.size(); i++)
            {
                if (myatoms[i].pType() == eptAtom)
                {
                    if (ifff >= fff.size())
                    {
                        GMX_THROW(gmx::InternalError(gmx::formatString("Inconsistency: there are %lu atoms and shells, but only %zu forces", myatoms.size(), fff.size())));
                    }
                    for(int m = 0; m < DIM; m++)
                    {
                        thisForce.push_back({ fff[ifff][m], forces[i][m] });
                    }
                    ifff += 1;
                }
            }
            forceMap->push_back(std::move(thisForce));
        }
    }
    done_commrec(crtmp);
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
            auto fp = fs.findParameterTypeConst(bondId, morse_name[morseLENGTH]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].morse.b0A         =
                mtop->ffparams.iparams[gromacsType].morse.b0B     = myval;
            if (ltop)
            { 
                ltop->idef.iparams[gromacsType].morse.b0A     =
                ltop->idef.iparams[gromacsType].morse.b0B = myval;
            }
                    
            fp = fs.findParameterTypeConst(bondId, morse_name[morseDE]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].morse.cbA         =
                mtop->ffparams.iparams[gromacsType].morse.cbB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].morse.cbA     =
                ltop->idef.iparams[gromacsType].morse.cbB = myval;
            }   

            fp = fs.findParameterTypeConst(bondId, morse_name[morseBETA]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].morse.betaA         =
                mtop->ffparams.iparams[gromacsType].morse.betaB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].morse.betaA     =
                ltop->idef.iparams[gromacsType].morse.betaB = myval;
            }
            
            fp = fs.findParameterTypeConst(bondId, morse_name[morseD0]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].morse.D0A         =
                mtop->ffparams.iparams[gromacsType].morse.D0B     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].morse.D0A     =
                        ltop->idef.iparams[gromacsType].morse.D0B = myval;
            }
        }
        break;
    case F_BONDS:
        {
            auto fp = fs.findParameterTypeConst(bondId, bond_name[bondLENGTH]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].harmonic.rA         =
                mtop->ffparams.iparams[gromacsType].harmonic.rB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].harmonic.rA     =
                ltop->idef.iparams[gromacsType].harmonic.rB = myval;
            }
                        
            fp = fs.findParameterTypeConst(bondId, bond_name[bondKB]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].harmonic.krA         =
                mtop->ffparams.iparams[gromacsType].harmonic.krB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].harmonic.krA     =
                ltop->idef.iparams[gromacsType].harmonic.krB = myval;
            }
        }
        break;
    case F_ANGLES:
        {
            auto fp = fs.findParameterTypeConst(bondId, angle_name[angleANGLE]);
            // GROMACS uses degrees internally, therefore no conversion
            // myval = fp.internalValue();
            myval = fp.value();
            mtop->ffparams.iparams[gromacsType].harmonic.rA         =
                mtop->ffparams.iparams[gromacsType].harmonic.rB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].harmonic.rA     =
                ltop->idef.iparams[gromacsType].harmonic.rB = myval;
            }
                        
            fp = fs.findParameterTypeConst(bondId, angle_name[angleKT]);
            myval = fp.internalValue();
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
            auto fp = fs.findParameterTypeConst(bondId, pol_name[polALPHA]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].polarize.alpha = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].polarize.alpha = myval;
            }
        }
        break;
    case F_UREY_BRADLEY:
        {
            auto fp = fs.findParameterTypeConst(bondId, ub_name[ubANGLE]);
            // GROMACS uses degrees internally, therefore no conversion
            // double angle = fp.internalValue();
            double angle = fp.value();
            mtop->ffparams.iparams[gromacsType].u_b.thetaA         =
                mtop->ffparams.iparams[gromacsType].u_b.thetaB     = angle;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].u_b.thetaA     =
                ltop->idef.iparams[gromacsType].u_b.thetaB = angle;
            }

            fp = fs.findParameterTypeConst(bondId, ub_name[ubKT]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].u_b.kthetaA         =
                mtop->ffparams.iparams[gromacsType].u_b.kthetaB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].u_b.kthetaA     =
                ltop->idef.iparams[gromacsType].u_b.kthetaB = myval;
            }

            fp = fs.findParameterTypeConst(bondId, ub_name[ubR13]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].u_b.r13A         =
                mtop->ffparams.iparams[gromacsType].u_b.r13B     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].u_b.r13A     =
                ltop->idef.iparams[gromacsType].u_b.r13B = myval;
            }                     
                         
            fp = fs.findParameterTypeConst(bondId, ub_name[ubKUB]);
            myval = fp.internalValue();
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
            auto fp = fs.findParameterTypeConst(bondId, linang_name[linangA]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].linangle.aA         =
                mtop->ffparams.iparams[gromacsType].linangle.aB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].linangle.aA     =
                ltop->idef.iparams[gromacsType].linangle.aB = myval;
            }
            
            fp = fs.findParameterTypeConst(bondId, linang_name[linangKLIN]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].linangle.klinA         =
                mtop->ffparams.iparams[gromacsType].linangle.klinB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].linangle.klinA     =
                ltop->idef.iparams[gromacsType].linangle.klinB = myval;
            }

            mtop->ffparams.iparams[gromacsType].linangle.r13A         =
                mtop->ffparams.iparams[gromacsType].linangle.r13B     = 0;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].linangle.r13A     =
                ltop->idef.iparams[gromacsType].linangle.r13B = myval;
            }
                        
            mtop->ffparams.iparams[gromacsType].linangle.kUBA         =
                mtop->ffparams.iparams[gromacsType].linangle.kUBB     = 0;
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
            std::vector<double> parameters;
            for(int i = 0; i < fdihNR; i++)
            {
                auto p = fs.findParameterTypeConst(bondId, fdih_name[i]).value();
                parameters.push_back(p);
            }
            newparam->rbdihs.rbcA[0] =  parameters[0];
            newparam->rbdihs.rbcA[1] = -parameters[1];
            newparam->rbdihs.rbcA[2] =  parameters[2];
            newparam->rbdihs.rbcA[3] = -parameters[3];
            newparam->rbdihs.rbcA[4] =  0.0;
            newparam->rbdihs.rbcA[5] =  0.0;
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
            auto fp = fs.findParameterTypeConst(bondId, pdih_name[pdihANGLE]);
            // GROMACS uses degrees internally, therefore no conversion
            // double angle = fp.internalValue();
            myval   = fp.value();
            mtop->ffparams.iparams[gromacsType].pdihs.phiA         =
                mtop->ffparams.iparams[gromacsType].pdihs.phiB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].pdihs.phiA     =
                ltop->idef.iparams[gromacsType].pdihs.phiB = myval;
            }
                        
            fp    = fs.findParameterTypeConst(bondId, pdih_name[pdihKP]);
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].pdihs.cpA         =
                mtop->ffparams.iparams[gromacsType].pdihs.cpB     = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].pdihs.cpA     =
                ltop->idef.iparams[gromacsType].pdihs.cpB = myval;
            }
                    
            int mult = fs.findParameterTypeConst(bondId, pdih_name[pdihMULT]).value();
            mtop->ffparams.iparams[gromacsType].pdihs.mult = mult;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].pdihs.mult = mult;
            }
        }
        break;
    case F_IDIHS:
        {
            mtop->ffparams.iparams[gromacsType].harmonic.rA         =
                mtop->ffparams.iparams[gromacsType].harmonic.rB     = 0;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].harmonic.rA     =
                ltop->idef.iparams[gromacsType].harmonic.rB = myval;
            }
            
            auto fp    = fs.findParameterTypeConst(bondId, idih_name[idihKPHI]);
            auto myval = fp.internalValue();
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
            myval = fp.internalValue();
            mtop->ffparams.iparams[gromacsType].vsite.a         = myval;
            if (ltop)
            {
                ltop->idef.iparams[gromacsType].vsite.a = myval;
            }
        }
        break;
    case F_LJ:
    case F_COUL_SR:
    case F_BHAM:
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
        auto &fs = pd->findForcesConst(entry.first);
        for(auto &topentry: entry.second)
        {
            int gromacsType = -1;
            // First check whether we have this gromacsType already
            // TODO check multiple
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
            if (interaction_function[fs.fType()].nratoms > 0 &&
                interaction_function[fs.fType()].nratoms !=
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
                                  missingParameters  missing,
                                  bool               gromacsSupport)
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
    if (immStatus::OK == imm)
    {
        snew(atoms_, 1);
        state_change_natoms(state_, NAtom());
        imm = GenerateAtoms(pd, atoms_);
    }
    if (immStatus::OK == imm)
    {
        imm = checkAtoms(pd, atoms_);
    }
    /* Store bonds in harmonic potential list first, update type later */
    if (immStatus::OK == imm)
    {
        topology_ = new Topology(bondsConst());
    }
    
    if (immStatus::OK == imm)
    {
        topology_->build(pd, optimizedCoordinates_, 175.0, 5.0, missing);
        excls_ = topology_->gromacsExclusions();
    }
    if (immStatus::OK == imm)
    {
        qProps_.insert({ qType::Calc, QtypeProps(qType::Calc) });
        qProps_.insert({ qType::Elec, QtypeProps(qType::Elec) });
    }
    if (immStatus::OK == imm)
    {
        /* Center of charge */
        auto atntot = 0;
        rvec coc    = { 0 };
        for (auto i = 0; i < atoms_->nr; i++)
        {
            auto atn = atoms_->atom[i].atomnumber;
            atntot  += atn;
            for (auto m = 0; m < DIM; m++)
            {
                coc[m] += optimizedCoordinates_[i][m]*atn;
            }
        }
        /* Center of charge */
        svmul((1.0/atntot), coc, CenterOfCharge_);
        for (auto &qp : qProps_)
        {
            qp.second.setCenterOfCharge(CenterOfCharge_);
        }
        addBondVsites(fp, pd, atoms_);
        if (pd->polarizable())
        {
            addShells(debug, pd, atoms_);
            topology_->addShellPairs();
        }
        // Now we can add the atom structures, whether or not the FF
        // is polarizable.
        topology_->setAtoms(atoms_);
        nRealAtoms_ = 0;
        for(int i = 0; i < atoms_->nr; i++)
        {
            if (atoms_->atom[i].ptype == eptAtom)
            {
                nRealAtoms_ += 1;
            }
        }
        char **molnameptr = put_symtab(symtab_, getMolname().c_str());
        if (gromacsSupport)
        {
            // Generate mtop
            mtop_ = do_init_mtop(pd, molnameptr, atoms_,
                                 inputrec_, symtab_, nullptr);
        }
        // First create the identifiers for topology entries
        topology_->setIdentifiers(pd);
        if (missing != missingParameters::Generate)
        {
            // Fill the parameters
            topology_->fillParameters(pd);
            if (gromacsSupport)
            {
                // Now generate the mtop fields
                TopologyToMtop(topology_, pd, mtop_);
            }
        }
        if (excls_ && gromacsSupport)
        {
            excls_to_blocka(atoms_->nr, excls_, &(mtop_->moltype[0].excls));
        }
        
        if (pd->polarizable() && gromacsSupport)
        {
            // Update mtop internals to account for shell type
            srenew(mtop_->atomtypes.atomnumber,
                   get_atomtype_ntypes(gromppAtomtype_));
            for (auto i = 0; i < get_atomtype_ntypes(gromppAtomtype_); i++)
            {
                mtop_->atomtypes.atomnumber[i] = get_atomtype_atomnumber(i, gromppAtomtype_);
            }
            mtop_->ffparams.atnr = get_atomtype_ntypes(gromppAtomtype_);
        }
        if (nullptr == ltop_ && gromacsSupport)
        {
            // Generate ltop from mtop
            ltop_ = gmx_mtop_generate_local_top(mtop_, false);
        }
    }
    if (immStatus::OK == imm && missing != missingParameters::Generate)
    {
        std::vector<InteractionType> itUpdate;
        for(auto &entry : *(topology_->entries()))
        {
            itUpdate.push_back(entry.first);
        }
        UpdateIdef(pd, itUpdate, true);
    }
    if (immStatus::OK == imm && pd->polarizable() && gromacsSupport)
    {
        // Generate shell data structure
        shellfc_ = init_shell_flexcon(debug, mtop_, 0, 1, false);
    }
    if (immStatus::OK == imm)
    {
        imm = checkAtoms(pd, atoms_);
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
        fraghandler_ = new FragmentHandler(pd, optimizedCoordinates_, topology_->atoms(), 
                                           bondsConst(), 
                                           fragmentPtr(), shellRenumber_, missing);
        // Finally, extract frequencies etc.
        getHarmonics();
    }
    return imm;
}

void MyMol::addBondVsites(FILE          *fp,
                          const Poldata *pd,
                          t_atoms       *atoms)
{
    int     atomNrOld = atoms->nr;
    // First add virtual sites for bond shells if needed.
    auto   &vs2       = pd->findForcesConst(InteractionType::VSITE2);
    auto   &qt        = pd->findForcesConst(InteractionType::COULOMB);
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
            auto vszetaid  = vstype->interactionTypeToIdentifier(InteractionType::COULOMB);
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

bool MyMol::linearMolecule() const
{
    const auto myatoms = topology_->atoms();
    std::vector<int> core;
    for(size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == eptAtom)
        {
            core.push_back(i);
        }
    }
    if (core.size() <= 2)
    {
        return true;
    }
    t_pbc *pbc    = nullptr;
    real th_toler = 175; // Degrees
    bool linear   = true;
    for(size_t c = 2; c < core.size(); c++)
    {
        linear = linear && is_linear(optimizedCoordinates_[core[c-2]],
                                     optimizedCoordinates_[core[c-1]], 
                                     optimizedCoordinates_[core[c]], pbc, th_toler);
        if (!linear)
        {
            break;
        }
    }
    return linear;
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
    std::vector<TopologyEntry *> pols;
    
    /* Calculate the total number of Atom and Vsite particles and
     * generate the renumbering array.
     */
    shellRenumber_.resize(atoms->nr, 0);
    for (int i = 0; i < atoms->nr; i++)
    {
        auto atype                  = pd->findParticleType(*atoms->atomtype[i]);
        shellRenumber_[i]           = i + nshell;
        originalAtomIndex_.insert(std::pair<int, int>(shellRenumber_[i], i));
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
    auto &qt  = pd->findForcesConst(InteractionType::COULOMB);
    auto &fs  = pd->findForcesConst(InteractionType::POLARIZATION);
    
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
                    pp->addAtom(shellRenumber_[i]);
                    pp->addAtom(shellRenumber_[i]+1);
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
    init_t_atoms(newatoms, nParticles, false);
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
        newatoms->atom[shellRenumber_[i]]      = atoms->atom[i];
        newatoms->atomname[shellRenumber_[i]]  = put_symtab(symtab_, *atoms->atomname[i]);
        newatoms->atomtype[shellRenumber_[i]]  = put_symtab(symtab_, *atoms->atomtype[i]);
        newatoms->atomtypeB[shellRenumber_[i]] = put_symtab(symtab_, *atoms->atomtypeB[i]);
        copy_rvec(optimizedCoordinates_[i], newx[shellRenumber_[i]]);
        newname[shellRenumber_[i]].assign(*atoms->atomtype[i]);
        int resind = atoms->atom[i].resind;
        t_atoms_set_resinfo(newatoms, shellRenumber_[i], symtab_,
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
            auto j            = 1+shellRenumber_[i];
            // Add an exclusion for the shell
            snew(newexcls[shellRenumber_[i]].e, 1);
            newexcls[shellRenumber_[i]].e[0] = j;
            newexcls[shellRenumber_[i]].nr = 1;
            auto atomtypeName = get_atomtype_name(atoms->atom[i].type, gromppAtomtype_);
            auto fa           = pd->findParticleType(atomtypeName);
            auto shellid      = fa->interactionTypeToIdentifier(InteractionType::POLARIZATION);
            if (shellid.id().empty())
            {
                // This particle has no shell.
                continue;
            }
            auto shelltype    = pd->findParticleType(shellid.id());
            auto shellzetaid  = shelltype->interactionTypeToIdentifier(InteractionType::COULOMB);
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
            copy_rvec(optimizedCoordinates_[i], newx[j]);

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
    optimizedCoordinates_ = newx;

    /* Copy exclusions, empty the original first */
    sfree(excls_);
    excls_ = newexcls;
    topology_->renumberAtoms(shellRenumber_);
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
    // Tell gromacs to use the generic kernel only.
    gmx_nonbonded_setup(fr_, true);

    gmx::ArrayRef<const std::string>  tabbfnm;
    init_forcerec(nullptr, mdlog, fr_, nullptr, inputrec_, mtop_, cr->commrec(),
                  state_->box, tabfn, tabfn, tabbfnm, false, true, -1);
    gmx_omp_nthreads_set(emntBonded, 1);
    init_bonded_threading(nullptr, 1, &fr_->bondedThreading);
    setup_bonded_threading(fr_->bondedThreading, mtop_->natoms, false, ltop_->idef);
    wcycle_    = wallcycle_init(debug, 0, cr->commrec());

    MDatoms_  = new std::unique_ptr<gmx::MDAtoms>(new gmx::MDAtoms());
    *MDatoms_ = gmx::makeMDAtoms(nullptr, *mtop_, *inputrec_);
    atoms2md(mtop_, inputrec_, -1, nullptr, mtop_->natoms, MDatoms_->get());
    auto mdatoms = MDatoms_->get()->mdatoms();

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

void MyMol::updateMDAtoms()
{
    auto mdatoms = MDatoms_->get()->mdatoms();
    auto myatoms = topology_->atoms();
    for (auto i = 0; i < atoms_->nr; i++)
    {
        mdatoms->chargeA[i] = myatoms[i].charge();
        mdatoms->typeA[i]   = atoms_->atom[i].type;
        mdatoms->zetaA[i]   = atoms_->atom[i].zetaA;
        if (mdatoms->zetaB)
        {
            mdatoms->zetaB[i]   = atoms_->atom[i].zetaB;
        }
        if (nullptr != debug)
        {
            fprintf(debug, "QQQ Setting q[%d] to %g\n", i, mdatoms->chargeA[i]);
        }
    }
}

static void reset_f_e(int                      natoms, 
                      PaddedVector<gmx::RVec> *f_,
                      gmx_enerdata_t          *enerd_)
{
    for (int i = 0; i < natoms; i++)
    {
        clear_rvec((*f_)[i]);
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
}

double MyMol::calculateInteractionEnergy(const Poldata       *pd,
                                         const ForceComputer *forceComputer)
{
    auto &tops = fraghandler_->topologies();
    if (tops.size() <= 1)
    {
        return 0;
    }
    // First, compute the total energy
    std::vector<gmx::RVec> coords = xOriginal();
    std::vector<gmx::RVec> ftotal(coords.size());
    std::map<InteractionType, double> etot;
    (void) forceComputer->compute(pd, topology_, &coords, &ftotal, &etot);
    // Now compute interaction energies if there are fragments
    double Einter  = etot[InteractionType::EPOT];
    auto   &astart = fraghandler_->atomStart();
    for(size_t ff = 0; ff < tops.size(); ff++)
    {
        int natom = tops[ff].atoms().size();
        std::vector<gmx::RVec> forces(natom);
        std::vector<gmx::RVec> myx(natom);
        int j = 0;
        for (size_t i = astart[ff]; i < astart[ff+1]; i++)
        {
            copy_rvec(coords[i], myx[j]);
            clear_rvec(forces[j]);
            j++;
        }
        std::map<InteractionType, double> energies;
        (void) forceComputer->compute(pd, &tops[ff], &myx, &forces, &energies);
        Einter -= energies[InteractionType::EPOT];
        if (debug)
        {
            fprintf(debug, "%s Fragment %zu Epot %g\n", getMolname().c_str(), ff, 
                    energies[InteractionType::EPOT]);
        }
    }
    return Einter;
}

immStatus MyMol::calculateEnergyOld(const t_commrec                   *crtmp,
                                    std::vector<gmx::RVec>            *coordinates,
                                    PaddedVector<gmx::RVec>           *forces,
                                    std::map<InteractionType, double> *energies,
                                    real                              *shellForceRMS)
{
    auto          imm         = immStatus::OK;
    unsigned long force_flags = ~0;
    double        t           = 0;
    rvec          mu_tot      = { 0, 0, 0 };
    tensor        force_vir   = { { 0 } };
    auto          mdatoms     = MDatoms_->get()->mdatoms();
    updateMDAtoms();

    for(size_t i = 0; i < coordinates->size(); i++)
    {
        copy_rvec((*coordinates)[i], state_->x[i]);
    }
    // (Re)create v-sites if needed
    if (!vsite_)
    {
        vsite_  = new std::unique_ptr<gmx_vsite_t>(new gmx_vsite_t());
        *vsite_ = initVsite(*mtop_, crtmp);
    }
    constructVsitesGlobal(*mtop_, state_->x);
    
    // Set force and energy to zero
    reset_f_e(mtop_->natoms, forces, enerd_);
    
    if (nullptr != shellfc_)
    {
        if (debug)
        {
            fprintf(debug, "mol %s alpha %g\n",
                    getMolname().c_str(),
                    mtop_->ffparams.iparams[mtop_->moltype[0].ilist[F_POLARIZATION].iatoms[0]].polarize.alpha);
        }
        try
        {
            *shellForceRMS = relax_shell_flexcon(nullptr, crtmp, nullptr, false,
                                                 0, inputrec_,
                                                 true, force_flags, ltop_,
                                                 enerd_, fcd_, state_,
                                                 forces->arrayRefWithPadding(), force_vir, mdatoms,
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
        if (*shellForceRMS > inputrec_->em_tol && debug)
        {
            for (int i = 0;  i<F_NRE; i++)
            {
                auto ei = enerd_->term[i];
                if (ei != 0)
                {
                    fprintf(debug, "E[%s] = %g\n", interaction_function[i].name,
                                ei);
                }
            }
            fprintf(debug, "Shell minimization did not converge in %d steps for %s. RMS Force = %g.\n",
                    inputrec_->niter, getMolname().c_str(),
                    *shellForceRMS);
            pr_rvecs(debug, 0, "f", forces->rvec_array(), mtop_->natoms);
            imm = immStatus::ShellMinimization;
        }
    }
    else
    {
        do_force(debug, crtmp, nullptr, inputrec_, 0,
                 &nrnb_, wcycle_, ltop_,
                 &(mtop_->groups),
                 state_->box, state_->x.arrayRefWithPadding(), nullptr,
                 forces->arrayRefWithPadding(), force_vir, mdatoms,
                 enerd_, fcd_,
                 state_->lambda, nullptr,
                 fr_, vsite_->get(), mu_tot, t,
                 force_flags);
        *shellForceRMS = 0;
    }
    std::map<int, InteractionType> ifmap = {
    { F_BONDS,         InteractionType::BONDS              },
    { F_MORSE,         InteractionType::BONDS              },
    { F_ANGLES,        InteractionType::ANGLES             },
    { F_LINEAR_ANGLES, InteractionType::LINEAR_ANGLES      },
    { F_LJ,            InteractionType::VDW                },
    { F_BHAM,          InteractionType::VDW                },
    { F_COUL_SR,       InteractionType::COULOMB            },
    { F_POLARIZATION,  InteractionType::POLARIZATION       },
    { F_IDIHS,         InteractionType::IMPROPER_DIHEDRALS },
    { F_PDIHS,         InteractionType::PROPER_DIHEDRALS   },
    { F_FOURDIHS,      InteractionType::PROPER_DIHEDRALS   },
    { F_UREY_BRADLEY,  InteractionType::ANGLES             },
    { F_EPOT,          InteractionType::EPOT               }
    };

    for(const auto &ifm : ifmap)
    {
        if (enerd_->term[ifm.first] != 0)
        {
            energies->insert({ifm.second, enerd_->term[ifm.first]});
        }
    }
    auto xOri = xOriginal();
    for(size_t i = 0; i < coordinates->size(); i++)
    {
        copy_rvec(xOri[i], state_->x[i]);
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
        symmetrize_charges(bSymmetricCharges, topology_,
                           pd, symm_string, &symmetric_charges_);
    }
    else
    {
        auto natoms = atomsConst().size();
        for (size_t i = 0; i < natoms; i++)
        {
            symmetric_charges_.push_back(i);
        }
    }
}

immStatus MyMol::GenerateAcmCharges(const Poldata          *pd,
                                    const ForceComputer    *forceComp,
                                    std::vector<gmx::RVec> *forces)
{
    std::vector<double> qold;
    fraghandler_->fetchCharges(&qold);
    if (qold.size() != atomsConst().size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot fetch old charges for %s. #atom %lu #qold %zu",
                                                       getMolname().c_str(), atomsConst().size(), qold.size()).c_str()));
    }
    immStatus imm       = immStatus::OK;
    int       iter      = 0;
    bool      converged = false;
    double    EemRms    = 0;
    auto      natom     = atomsConst().size();
    // Make copy of the coordinates
    std::vector<gmx::RVec> coords = optimizedCoordinates_;
    std::map<InteractionType, double> energies;
    do
    {
        if (eQgen::OK == fraghandler_->generateCharges(debug, getMolname(),
                                                       coords, pd, atoms()))
        {
            (void) forceComp->compute(pd, topology_, &coords, forces, &energies);
            EemRms = 0;
            std::vector<double> qnew;
            fraghandler_->fetchCharges(&qnew);
            GMX_RELEASE_ASSERT(qold.size()==qnew.size(), "Cannot fetch new charges");
            for (size_t i = 0; i < qnew.size(); i++)
            {
                EemRms  += gmx::square(qnew[i] - qold[i]);
                qold[i]  = qnew[i];
            }
            EemRms   /= natom;
            converged = (EemRms < qTolerance_) || !haveShells();
            iter++;
        }
        else
        {
            imm = immStatus::ChargeGeneration;
        }
    }
    while (imm == immStatus::OK && (!converged) && (iter < maxQiter_));
    if (!converged)
    {
        printf("Alexandria Charge Model did not converge to %g. rms: %g\n",
               qTolerance_, sqrt(EemRms));
    }
    auto myatoms = atoms();
    for(size_t i = 0; i < natom; i++)
    {
        (*myatoms)[i].setCharge(qold[i]);
    }
    auto qcalc = qTypeProps(qType::Calc);
    qcalc->setQ(atomsConst());
    qcalc->setX(x());
    return imm;
}

immStatus MyMol::GenerateCharges(const Poldata             *pd,
                                 const ForceComputer       *forceComp,
                                 const gmx::MDLogger       &mdlog,
                                 const CommunicationRecord *cr,
                                 ChargeGenerationAlgorithm  algorithm,
                                 const std::vector<double> &qcustom,
                                 std::vector<gmx::RVec>    *forces)
{
    immStatus imm         = immStatus::OK;
    bool      converged   = false;
    auto      &qt         = pd->findForcesConst(InteractionType::COULOMB);
    auto      iChargeType = name2ChargeType(qt.optionValue("chargetype"));

    if (nullptr != mtop_)
    {
        GenerateGromacs(mdlog, cr, nullptr, iChargeType);
    }
    // TODO check whether this needed
    std::map<InteractionType, double> energies;
    auto myatoms = atoms();
    if (algorithm == ChargeGenerationAlgorithm::Custom)
    {
        GMX_RELEASE_ASSERT(atomsConst().size() == qcustom.size(),
                           gmx::formatString("Number of custom charges %lu does not match the number of atoms %lu", qcustom.size(), atomsConst().size()).c_str());
    }
    else if (algorithm == ChargeGenerationAlgorithm::NONE)
    {
        algorithm = pd->chargeGenerationAlgorithm();
        // Check whether there are free charges
        bool allFixed = true;
        for (size_t i = 0; i < myatoms->size(); i++)
        {
            auto atype = (*myatoms)[i].ffType();
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
            for (size_t i = 0; i < myatoms->size(); i++)
            {
                auto atype = (*myatoms)[i].ffType();
                auto ptype = pd->findParticleType(atype);
                auto qval  = ptype->parameterConst("charge").value();
                (*myatoms)[i].setCharge(qval);
            }
            // If we have shells, we still have to minimize them,
            // but we may want to know the energies anyway.
            (void) forceComp->compute(pd, topology_, &optimizedCoordinates_,
                                      forces, &energies);
            if (haveShells())
            {
                auto qcalc = qTypeProps(qType::Calc);
                qcalc->setQ(*myatoms);
                qcalc->setX(optimizedCoordinates_);
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
            for (auto exper : experimentConst())
            {
                int i = 0;
                for (auto &ca : exper.calcAtomConst())
                {
                    if (ca.hasCharge(qtmap[algorithm]))
                    {
                        (*myatoms)[i].setCharge(ca.charge(qtmap[algorithm]));
                        i++;
                    }
                    else
                    {
                        gmx_fatal(FARGS, "No charge type %s for %s",
                                  qTypeName(qtmap[algorithm]).c_str(), getMolname().c_str());
                    }
                }
            }
            // TODO check this. Copy charges to topology
            // topology_->setAtoms(myatoms);
            return immStatus::OK;
        }
    case ChargeGenerationAlgorithm::Custom:
        {
            for (size_t i = 0; i < myatoms->size(); i++)
            {
                (*myatoms)[i].setCharge(qcustom[i]);
            }

            return immStatus::OK;
        }
    case ChargeGenerationAlgorithm::ESP:
        {
            double chi2[2] = {1e8, 1e8};
            int    cur     = 0;
            int    maxiter = 5;
            int    iter    = 0;
            
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
            do
            {
                auto qq = qcalc->charge();
                GMX_RELEASE_ASSERT(myatoms->size() == qq.size(),
                                   gmx::formatString("Number of particles (%lu) differs from number of charges (%lu)", 
                                                     myatoms->size(),
                                                     qq.size()).c_str());
                for (size_t i = 0; i < myatoms->size(); i++)
                {
                    (*myatoms)[i].setCharge(qq[i]);
                }
                // Copy charges to topology
                chi2[cur] = forceComp->compute(pd, topology_, &optimizedCoordinates_,
                                               forces, &energies);
                qcalc->setX(x());
                qcalc->qgenResp()->optimizeCharges(pd->getEpsilonR());
                qcalc->qgenResp()->calcPot(pd->getEpsilonR());
                qcalc->copyRespQ();
                if (debug)
                {
                    fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
                }
                converged = (fabs(chi2[cur] - chi2[1-cur]) < qTolerance_) || !pd->polarizable();
                cur       = 1-cur;
                iter++;
            }
            while ((!converged) && (iter < maxiter));
            auto qq = qcalc->charge();
            for (size_t i = 0; i < myatoms->size(); i++)
            {
                (*myatoms)[i].setCharge(qq[i]);
            }
            // Copy charges to topology
            // topology_->setAtoms(myatoms);
        }
        break;
    case ChargeGenerationAlgorithm::EEM:
    case ChargeGenerationAlgorithm::SQE:
        {
            imm = GenerateAcmCharges(pd, forceComp, forces);
        }
        break;
    }
    return imm;
}

void MyMol::CalcPolarizability(const Poldata       *pd,
                               const ForceComputer *forceComp)
{
    auto natoms = atomsConst().size();
    std::vector<gmx::RVec> coordinates(natoms);
    for(size_t i = 0; i < natoms; i++)
    {
        copy_rvec(state_->x[i], coordinates[i]);
    }
    forceComp->calcPolarizability(pd, topology_, &coordinates,
                                  qTypeProps(qType::Calc));
}

void MyMol::PrintConformation(const char *fn)
{
    char title[STRLEN];
    
    sprintf(title, "%s processed by ACT - The Alexandria Chemistry Tookit",
            getMolname().c_str());
    write_sto_conf(fn, title, gmxAtoms(), as_rvec_array(x().data()),
                   nullptr, epbcNONE, state_->box);
}

static void add_tensor(std::vector<std::string> *commercials,
                       const char               *title,
                       const char               *unit,
                       const std::vector<double> &Q)
{
    double fac = convertFromGromacs(1.0, unit);
    char buf[256];
    snprintf(buf, sizeof(buf), "%s:\n"
             "; ( %6.2f %6.2f %6.2f )\n"
             "; (       %6.2f %6.2f )\n"
             "; (             %6.2f )\n",
             title,
             fac*Q[0], fac*Q[1], fac*Q[2],
             fac*Q[3], fac*Q[4], fac*Q[5]);
    commercials->push_back(buf);
}

void MyMol::PrintTopology(const char                *fn,
                          bool                       bVerbose,
                          const Poldata             *pd,
                          const ForceComputer       *forceComp,
                          const CommunicationRecord *cr,
                          const std::string         &method,
                          const std::string         &basis,
                          bool                       bITP)
{
    char                     buf[256];
    t_mols                   printmol;
    std::vector<std::string> commercials;
    std::vector<double>      vec;
    double                   T = -1;
    std::string              myref;
    auto &qt         = pd->findForcesConst(InteractionType::COULOMB);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));
    std::string              mylot       = makeLot(method, basis);

    FILE *fp = gmx_ffopen(fn, "w");
    
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

    snprintf(buf, sizeof(buf), "Total Mass = %.3f (Da)", totalMass());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Total Charge = %d (e)", totalCharge());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Charge Type  = %s\n",
             chargeTypeName(iChargeType).c_str());
    commercials.push_back(buf);
    
    auto qcalc = qTypeProps(qType::Calc);
    auto qelec = qTypeProps(qType::Elec);
    qcalc->setQ(atomsConst());
    qcalc->setX(x());
    qcalc->calcMoments();
    
    T = -1;
    for(auto &mpo : mpoMultiPoles)
    {
        auto gp = qmProperty(mpo, T, JobType::OPT);
        if (gp)
        {
            auto vec = gp->getVector();
            qelec->setMultipole(mpo, vec);
            auto mymu = qelec->getMultipole(mpo);
            commercials.push_back(gmx::formatString("%s %s (%s)\n",
                                                    mylot.c_str(), mpo_name(mpo), gp->getUnit()));
            for(auto &fmp : formatMultipole(mpo, mymu))
            {
                commercials.push_back(fmp);
            }
            if (qcalc->hasMultipole(mpo))
            {
                auto mymu = qcalc->getMultipole(mpo);
                commercials.push_back(gmx::formatString("Alexandria %s (%s)\n", mpo_name(mpo), gp->getUnit()));
                for(auto &fmp : formatMultipole(mpo, mymu))
                {
                    commercials.push_back(fmp);
                }
            }
        }
    }

    if (nullptr != cr)
    {
        CalcPolarizability(pd, forceComp);
        auto qcalc = qTypeProps(qType::Calc);
        auto acalc = qcalc->polarizabilityTensor();
        std::vector<double> ac = { acalc[XX][XX], acalc[XX][YY], acalc[XX][ZZ],
                                   acalc[YY][YY], acalc[YY][ZZ], acalc[ZZ][ZZ] };
        auto unit = mpo_unit2(MolPropObservable::POLARIZABILITY);
        add_tensor(&commercials, "Alexandria Polarizability components (A^3)", unit, ac);
        
        snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability: %.2f (A^3)\n",
                 qcalc->isotropicPolarizability());
        commercials.push_back(buf);
        
        snprintf(buf, sizeof(buf), "Alexandria Anisotropic Polarizability: %.2f (A^3)\n",
                 qcalc->anisotropicPolarizability());
        commercials.push_back(buf);
        
        T = -1;
        auto gp = qmProperty(MolPropObservable::POLARIZABILITY,
                             T, JobType::OPT);
        if (gp)
        {
            auto qelec = qTypeProps(qType::Elec);
            auto aelec = qelec->polarizabilityTensor();
            std::vector<double> ae = { aelec[XX][XX], aelec[XX][YY], aelec[XX][ZZ],
                                       aelec[YY][YY], aelec[YY][ZZ], aelec[ZZ][ZZ] };
            snprintf(buf, sizeof(buf), "%s + Polarizability components (A^3)", mylot.c_str());
            add_tensor(&commercials, buf, unit, ae);
            snprintf(buf, sizeof(buf), "%s Isotropic Polarizability: %.2f (A^3)\n",
                     mylot.c_str(), qelec->isotropicPolarizability());
            commercials.push_back(buf);
            snprintf(buf, sizeof(buf), "%s Anisotropic Polarizability: %.2f (A^3)\n",
                     mylot.c_str(), qelec->anisotropicPolarizability());
            commercials.push_back(buf);
        }
    }

    // TODO write a replacement for this function
    print_top_header(fp, pd, bHaveShells_, commercials, bITP);
    write_top(fp, printmol.name, gmxAtoms(),
              topology_, excls_, gromppAtomtype_, pd);
    if (!bITP)
    {
        print_top_mols(fp, printmol.name, pd->filename().c_str(),
                       nullptr, 0, nullptr, 1, &printmol);
    }
    if (bVerbose)
    {
        for (auto &entry : *(topology_->entries()))
        {
            auto &fs = pd->findForcesConst(entry.first);
            int ftype = fs.fType();
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
    
    gmx_ffclose(fp);
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
    auto &qt         = pd->findForcesConst(InteractionType::COULOMB);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));

    if (potfn || hisfn || rhofn || difffn || pdbdifffn)
    {
        char     *gentop_version = (char *)"gentop v0.99b";
        auto qc = qProps_.find(qType::Calc);
        GMX_RELEASE_ASSERT(qc != qProps_.end(), "Cannot find alexandria charge information");
        qc->second.setQ(atomsConst());
        qc->second.setX(x());
        qc->second.qgenResp()->calcPot(pd->getEpsilonR());
        qc->second.qgenResp()->potcomp(pcfn, atomsConst(),
                                       as_rvec_array(state_->x.data()),
                                       pdbdifffn, oenv);

        /* This has to be done before the grid is f*cked up by
           writing a cube file */
        QgenResp qCalc(*qc->second.qgenResp());
        QgenResp grref(*qc->second.qgenResp());

        if (reffn)
        {
            grref.setAtomInfo(atomsConst(), pd, x(), totalCharge());
            grref.setAtomSymmetry(symmetric_charges_);
            grref.readCube(reffn, FALSE);
        }
        else
        {
            qCalc.makeGrid(spacing, border, x());
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

void MyMol::calcEspRms(const Poldata *pd)
{
    int   natoms  = 0;
    auto &myatoms = atomsConst();
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == eptAtom)
        {
            natoms++;
        }
    }
    gmx::HostVector<gmx::RVec> myx(nRealAtoms());
    std::vector<ActAtom> realAtoms;
    size_t ii = 0;
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == eptAtom)
        {
            realAtoms.push_back(myatoms[i]);
            copy_rvec(optimizedCoordinates_[i], myx[ii]);
            ii++;
        }
    }
    
    auto qcalc   = qTypeProps(qType::Calc);
    auto qgrcalc = qcalc->qgenResp();
    for(auto &i : qProps_)
    {
        auto qi = i.first;
        if (qType::Calc == qi)
        {
            //qgrcalc->setAtomInfo(atomsConst(), pd, x(), totalCharge());
            qgrcalc->updateAtomCharges(atomsConst());
            qgrcalc->calcPot(pd->getEpsilonR());
        }
        else if (qType::Elec != qi)
        {
            QgenResp *qgr = i.second.qgenResp();
            qgr->setChargeType(ChargeType::Point);
            qgr->setAtomInfo(realAtoms, pd, x(), totalCharge());
            qgr->updateAtomCharges(i.second.charge());
            for (size_t j = realAtoms.size(); j < qgrcalc->nEsp(); j++)
            {
                auto &ep = qgrcalc->espPoint(j);
                auto &r  = ep.esp();
                qgr->addEspPoint(r[XX], r[YY], r[ZZ], ep.v());
            }
            qgr->calcPot(pd->getEpsilonR());
        }
    }
}

void MyMol::getHarmonics()
{
    for(auto &mpo : { MolPropObservable::FREQUENCY, 
                     MolPropObservable::INTENSITY })
    {
        std::vector<GenericProperty *> harm;
        for (auto &ee : experimentConst())
        {
            if (ee.hasMolPropObservable(mpo))
            {
                harm = ee.propertyConst(mpo);
                break;
            }
        }
        if (!harm.empty())
        {
            for(auto &ff : harm[0]->getVector())
            {
                if (mpo == MolPropObservable::FREQUENCY)
                {
                    ref_frequencies_.push_back(ff);
                }
                else
                {
                    ref_intensities_.push_back(ff);
                }
            }
        }
    }
}

const real *MyMol::energyTerms() const
{
    return enerd_->term;
}
        
immStatus MyMol::getExpProps(const std::map<MolPropObservable, iqmType> &iqm,
                             double                                      T)
{
    int                 natom = 0;
    std::vector<double> vec;
    std::string         myref;
    std::string         mylot;
    immStatus           imm        = immStatus::OK;
    
    auto &myatoms = atomsConst();
    GMX_RELEASE_ASSERT(myatoms.size() > 0, "No atoms!");
    
    // Make a copy of the coordinates without shells
    std::vector<gmx::RVec> xatom(myatoms.size());
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == eptAtom ||
            myatoms[i].pType() == eptNucleus)
        {
            copy_rvec(x()[i], xatom[natom]);
            natom++;
        }
    }
    xatom.resize(natom);
    bool foundNothing = true;
    for (const auto &miq : iqm)
    {
        auto mpo = miq.first;
        switch (mpo)
        {
        case MolPropObservable::CHARGE:
        case MolPropObservable::POTENTIAL:
            {
                std::string conf;
                auto ei = findExperimentConst(JobType::OPT);
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
                    foundNothing = false;
                }
            }
            break;
        case MolPropObservable::DELTAE0:
        case MolPropObservable::DHFORM:
        case MolPropObservable::DGFORM:
        case MolPropObservable::ZPE:
            {
                auto gp = static_cast<const MolecularEnergy *>(qmProperty(mpo, T, JobType::OPT));
                if (gp)
                {
                    energy_.insert(std::pair<MolPropObservable, double>(mpo, gp->getValue()));
                    foundNothing = false;
                }
            }
            break;
        case MolPropObservable::DIPOLE:
        case MolPropObservable::QUADRUPOLE:
        case MolPropObservable::OCTUPOLE:
        case MolPropObservable::HEXADECAPOLE:
            {
                auto gp = static_cast<const MolecularMultipole *>(qmProperty(mpo, T, JobType::OPT));
                if (gp)
                {
                    qProps_.find(qType::Elec)->second.setMultipole(mpo, gp->getVector());
                    foundNothing = false; 
                }
            }
            break;
        case MolPropObservable::POLARIZABILITY:
            {
                auto gp = static_cast<const MolecularPolarizability *>(qmProperty(mpo, T, JobType::OPT));
                if (gp)
                {
                    auto qelec = qTypeProps(qType::Elec);
                    qelec->setPolarizabilityTensor(gp->getTensor());
                    foundNothing = false;
                }
            }
            break;
        default:
            break;
        }
    }
    if (foundNothing)
    {
        imm = immStatus::NoData;
    }
    return imm;
}

void MyMol::UpdateIdef(const Poldata                      *pd,
                       const std::vector<InteractionType> &iTypes,
                       bool                                updateZeta)
{
    topology_->fillParameters(pd);
    // TODO Check whether this is sufficient for updating the particleTypes
    if (updateZeta)
    {
        // Update the electronegativity parameters
        zetaToAtoms(pd, gmxAtoms());
    }
    if (mtop_ == nullptr)
    {
        return;
    }
    for(const auto &iType : iTypes)
    {
        if (debug)
        {
            fprintf(debug, "UpdateIdef for %s\n", interactionTypeToString(iType).c_str());
        }
        // The rest may not be needed anymore!
        
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
                auto &fs   = pd->findForcesConst(iType);
                auto entry = topology_->entry(iType);
                for (size_t i = 0; i < entry.size(); i++)
                {
                    // TODO check multiple
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
}

void MyMol::initQgenResp(const Poldata     *pd,
                         real               watoms,
                         int                maxESP)
{
    std::string        mylot;
    auto &qt         = pd->findForcesConst(InteractionType::COULOMB);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));
    auto qp          = qTypeProps(qType::Calc);
    QgenResp *qgr = qp->qgenResp();
    qgr->setChargeType(iChargeType);
    qgr->setAtomInfo(atomsConst(), pd, x(), totalCharge());
    qp->setQ(atomsConst());
    qp->setX(x());
    qgr->setAtomSymmetry(symmetric_charges_);
    qgr->summary(debug);

    int natoms = nRealAtoms();

    std::random_device               rd;
    std::mt19937                     gen(rd());  
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    double                           cutoff = 0.01*maxESP;
 
    auto ci = findExperimentConst(JobType::OPT);
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

void MyMol::plotEspCorrelation(const Poldata          *pd,
                               const char             *espcorr,
                               const gmx_output_env_t *oenv,
                               const ForceComputer    *forceComp)
{
    if (espcorr && oenv)
    {
        auto qgr   = qTypeProps(qType::Calc)->qgenResp();
        qgr->updateAtomCharges(atomsConst());
        qgr->updateAtomCoords(x());
        std::vector<gmx::RVec> forces(atomsConst().size());
        std::map<InteractionType, double> energies;
        (void) forceComp->compute(pd, topology_, &optimizedCoordinates_,
                                  &forces, &energies);
        qgr->calcPot(1.0);
        qgr->plotLsq(oenv, espcorr);
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
