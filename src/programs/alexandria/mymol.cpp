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
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
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

static const double A2CM = E_CHARGE*1.0e-10;        /* e Angstrom to Coulomb meter */

static const double CM2D = SPEED_OF_LIGHT*1.0e+24;  /* Coulomb meter to Debye */

static inline int delta(int a, int b) { return ( a == b ) ? 1 : 0; }

static inline double e2d(double a) {return a*ENM2DEBYE; }

static void vsiteType_to_atomType(const std::string &vsiteType, std::string *atomType)
{
    std::size_t pos = vsiteType.find("L");
    *atomType        = vsiteType.substr (0, pos);
}

const char *qTypeName(qType qt)
{
    switch (qt)
    {
        case qtESP:       return "ESP";
        case qtMulliken:  return "Mulliken";
        case qtHirshfeld: return "Hirshfeld";
        case qtCM5:       return "CM5";
        case qtCalc:      return "Calculated";
        case qtElec:      return "Electronic";
        default:
            return "Unknown charge type";
    }
}

const char *immsg(immStatus imm)
{
    static const char *msg[immNR] = {
        "Unknown status",
        "OK",
        "Zero Dipole",
        "No Quadrupole",
        "Charged",
        "Atom type problem",
        "Atom number problem",
        "Converting from molprop",
        "Determining bond order",
        "RESP Initialization",
        "Charge generation",
        "Shell minimization",
        "Requested level of theory missing",
        "QM Inconsistency (ESP dipole does not match Elec)",
        "Not in training set",
        "No experimental data",
        "Generating shells",
        "Generating bonds",
        "Communicating MolProp",
        "Zeta is zero",
        "The number of data is lower than mindata",
        "No Dipole moment",
        "NotSupportedBond",
        "NotSupportedAngle",
        "NotSupportedDihedral"
    };
    return msg[imm];
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


MyMol::MyMol() : gvt_(evtALL)

{
    bHaveShells_       = false;
    bHaveVSites_       = false;
    bNeedVsites_       = false;
    cgnr_              = nullptr;
    immAtoms_          = immOK;
    immTopology_       = immOK;
    immCharges_        = immOK;
    shellfc_           = nullptr;
    vsite_             = nullptr;
    myforce_           = new MyForceProvider;
    snew(symtab_, 1);
    open_symtab(symtab_);
    gromppAtomtype_    = init_atomtype();
    mtop_          = nullptr;
    fr_            = nullptr;
    ltop_          = nullptr;
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
    for (int j = 0; j < qtNR; j++)
    {
        clear_rvec(mu_qm_[j]);
        clear_mat(Q_qm_[j]);
    }
}

bool MyMol::IsSymmetric(real toler)
{
    int       i, j, m;
    real      mm, tm;
    rvec      com, test;
    gmx_bool *bSymm, bSymmAll;

    clear_rvec(com);
    tm = 0;
    for (i = 0; i < atoms_->nr; i++)
    {
        mm  = atoms_->atom[i].m;
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
    for (i = 0; i < atoms_->nr; i++)
    {
        rvec_dec(state_->x[i], com);
    }

    snew(bSymm, atoms_->nr);
    for (i = 0; i < atoms_->nr; i++)
    {
        bSymm[i] = (norm(state_->x[i]) < toler);
        for (j = i+1; (j < atoms_->nr) && !bSymm[i]; j++)
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
    for (i = 0; i < atoms_->nr; i++)
    {
        bSymmAll = bSymmAll && bSymm[i];
    }
    sfree(bSymm);
    for (i = 0; i < atoms_->nr; i++)
    {
        rvec_inc(state_->x[i], com);
    }

    return bSymmAll;
}

void MyMol::findInPlaneAtoms(int ca, std::vector<int> &atoms)
{
    int bca = 0;
    /*First try to find the atom bound to the central atom (ca).*/
    for (auto &bi : bondConst())
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
    for (auto bi : bondConst())
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
    for (auto &bi : bondConst())
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
                                    bool           bUseVsites)
{
    std::vector < std::vector < unsigned int>> bonds;
    std::vector<int> nbonds;
    t_pbc            pbc;
    real             th_toler = 175;
    real             ph_toler = 5;

    rvec            *x = as_rvec_array(state_->x.data());

    set_pbc(&pbc, epbcNONE, state_->box);

    bonds.resize(atoms_->nr);
    for (auto &bi : bondConst())
    {
        // Store bonds bidirectionally to get the number correct
        bonds[bi.getAi() - 1].push_back(bi.getAj() - 1);
        bonds[bi.getAj() - 1].push_back(bi.getAi() - 1);
    }
    nbonds.resize(atoms_->nr);
    for (auto i = 0; i < atoms_->nr; i++)
    {
        nbonds[i] = bonds[i].size();
    }
    for (auto i = 0; i < atoms_->nr; i++)
    {
        /* Now test initial geometry */
        if ((bonds[i].size() == 2) &&
            is_linear(x[i], x[bonds[i][0]], x[bonds[i][1]],
                      &pbc, th_toler))
        {
            if (nullptr != debug)
            {
                fprintf(debug, "found linear angle %s-%s-%s in %s\n",
                        *atoms_->atomtype[bonds[i][0]],
                        *atoms_->atomtype[i],
                        *atoms_->atomtype[bonds[i][1]],
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
                        *atoms_->atomtype[i],
                        *atoms_->atomtype[bonds[i][0]],
                        *atoms_->atomtype[bonds[i][1]],
                        *atoms_->atomtype[bonds[i][2]],
                        getMolname().c_str());
            }
            gvt_.addPlanar(i, bonds[i][0], bonds[i][1], bonds[i][2],
                           &nbonds[0]);
        }
        if (bUseVsites)
        {
            const auto atype(*atoms_->atomtype[i]);
            if (IsVsiteNeeded(atype, pd))
            {
                std::vector<int> atoms;
                auto             vsite = pd->findVsite(atype);
                if (vsite->type() == evtIN_PLANE)
                {
                    atoms.push_back(i);
                    findInPlaneAtoms(i, atoms);
                    if (vsite->ncontrolatoms() == static_cast<int>(atoms.size()))
                    {
                        gvt_.addInPlane(vsite->ncontrolatoms(),
                                        vsite->nvsite(),
                                        atoms[0], atoms[1],
                                        atoms[2], atoms[3]);
                    }
                }
                else if (vsite->type() == evtOUT_OF_PLANE)
                {
                    atoms.push_back(i);
                    findOutPlaneAtoms(i, atoms);
                    if (vsite->ncontrolatoms() == static_cast<int>(atoms.size()))
                    {
                        gvt_.addOutPlane(vsite->ncontrolatoms(),
                                         vsite->nvsite(),
                                         atoms[0], atoms[1],
                                         atoms[2]);
                    }
                }
                if (!bNeedVsites_)
                {
                    bNeedVsites_ = true;
                }
            }
        }
    }
    auto anr = atoms_->nr;
    gvt_.generateSpecial(pd, bUseVsites, atoms_,
                         &x, plist_, symtab_, gromppAtomtype_, &excls_, state_);
    bHaveVSites_ = (atoms_->nr > anr);
}

/*
 * Make Harmonic Angles, Proper Dihedrals, and 14 Pairs.
 * This needs the bonds to be F_BONDS.
 */
void MyMol::MakeAngles(bool bPairs,
                       bool bDihs)
{
    t_nextnb                            nnb;
    t_restp                             rtp;
    t_params                            plist[F_NRE];

    init_plist(plist);
    for (auto &pw : plist_)
    {
        if (F_BONDS == pw.getFtype())
        {
            pr_alloc(pw.nParam(), &(plist[F_BONDS]));
            auto i = 0;
            for (auto pi = pw.beginParam(); (pi < pw.endParam()); ++pi)
            {
                for (auto j = 0; j < MAXATOMLIST; j++)
                {
                    plist[F_BONDS].param[i].a[j] = pi->a[j];
                }
                for (auto j = 0; j < MAXFORCEPARAM; j++)
                {
                    plist[F_BONDS].param[i].c[j] = pi->c[j];
                }
                i++;
            }
            plist[F_BONDS].nr = i;
            break;
        }
    }
    /* Make Harmonic Angles and Proper Dihedrals */
    snew(excls_, atoms_->nr);
    init_nnb(&nnb, atoms_->nr, nexcl_ + 2);
    gen_nnb(&nnb, plist);

    print_nnb(&nnb, "NNB");
    rtp.bKeepAllGeneratedDihedrals    = bDihs;
    rtp.bRemoveDihedralIfWithImproper = bDihs;
    rtp.bGenerateHH14Interactions     = bPairs;
    rtp.nrexcl                        = nexcl_;

    gen_pad(&nnb, atoms_, &rtp, plist, excls_, nullptr, false);

    t_blocka *EXCL;
    snew(EXCL, 1);
    if (debug)
    {
        fprintf(debug, "Will generate %d exclusions for %d atoms\n",
                nexcl_, atoms_->nr);
    }
    generate_excl(nexcl_, atoms_->nr, plist, &nnb, EXCL);
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
        for (auto i = 0; i < atoms_->nr; i++)
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

    cp_plist(plist, F_ANGLES, eitANGLES, plist_);

    if (bDihs)
    {
        cp_plist(plist, F_PDIHS, eitPROPER_DIHEDRALS, plist_);
    }
    if (bPairs)
    {
        /* Make 1-4 table */
        cp_plist(plist, F_LJ14, eitLJ14, plist_);
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
                               const std::string &method,
                               const std::string &basis,
                               std::string       *mylot)
{
    double              xx, yy, zz;
    int                 natom = 0;
    immStatus           imm   = immOK;

    ExperimentIterator  ci = getCalc(method, basis, mylot);
    if (ci < EndExperiment())
    {
        t_param nb;
        memset(&nb, 0, sizeof(nb));
        natom = 0;
        init_t_atoms(atoms_, ci->NAtom(), false);
        snew(atoms_->atomtype, ci->NAtom());
        snew(atoms_->atomtypeB, ci->NAtom());
        int res0 = -1;
        int nres =  0;
        for (auto &cai : ci->calcAtomConst())
        {
            auto myunit = cai.getUnit();
            cai.getCoords(&xx, &yy, &zz);
            int resnr = cai.ResidueNumber();
            if (resnr != res0)
            {
                res0  = resnr;
                nres += 1;
            }
            state_->x[natom][XX] = convertToGromacs(xx, myunit);
            state_->x[natom][YY] = convertToGromacs(yy, myunit);
            state_->x[natom][ZZ] = convertToGromacs(zz, myunit);

            double q = 0;
            for (auto &qi : cai.atomicChargeConst())
            {
                // TODO Clean up this mess.
                if (qi.getType().compare("ESP") == 0)
                {
                    myunit = string2unit((char *)qi.getUnit().c_str());
                    q      = convertToGromacs(qi.getQ(), myunit);
                    break;
                }
            }
            atoms_->atom[natom].q      =
                atoms_->atom[natom].qB = q;
            atoms_->atom[natom].resind = resnr;
            t_atoms_set_resinfo(atoms_, natom, symtab_, cai.ResidueName().c_str(),
                                atoms_->atom[natom].resind, ' ', 
                                cai.chainId(), cai.chain());
            atoms_->atomname[natom]        = put_symtab(symtab_, cai.getName().c_str());

            // First set the atomtype
            atoms_->atomtype[natom]      =
                atoms_->atomtypeB[natom] = put_symtab(symtab_, cai.getObtype().c_str());
            auto atype = pd->findAtype(cai.getObtype());
            if (atype != pd->getAtypeEnd())
            {
                atoms_->atom[natom].m      =
                    atoms_->atom[natom].mB = atype->mass();
                atoms_->atom[natom].atomnumber = atype->atomnumber();
                strncpy(atoms_->atom[natom].elem, atype->getElem().c_str(), sizeof(atoms_->atom[natom].elem)-1);
                
                natom++;
            }
            else
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find atomtype %s in poldata", cai.getObtype().c_str()).c_str()));
            }
        }
        for (auto i = 0; i < natom; i++)
        {
            atoms_->atom[i].type      =
                atoms_->atom[i].typeB = add_atomtype(gromppAtomtype_, symtab_,
                                                     &(atoms_->atom[i]),
                                                     *atoms_->atomtype[i],
                                                     &nb, 0,
                                                     atoms_->atom[i].atomnumber);
        }
        atoms_->nr   = natom;
        atoms_->nres = nres;
        //printf("natom %d nres %d\n", natom, nres);
    }
    else
    {
        imm = immLOT;
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Tried to convert %s to gromacs. LOT is %s/%s. Natoms is %d\n",
                getMolname().c_str(),
                method.c_str(), basis.c_str(), natom);
    }

    return imm;
}

immStatus MyMol::checkAtoms(const Poldata *pd)
{
    auto nmissing = 0;
    for (auto i = 0; i < atoms_->nr; i++)
    {
        const auto atype(*atoms_->atomtype[i]);
        auto       fa = pd->findAtype(atype);
        if (fa == pd->getAtypeEnd())
        {
            printf("Could not find a force field entry for atomtype %s atom %d in compound '%s'\n",
                   *atoms_->atomtype[i], i+1,
                   getMolname().c_str());
            nmissing++;
        }
    }
    if (nmissing > 0)
    {
        return immAtomTypes;
    }
    return immOK;
}

immStatus MyMol::zeta2atoms(const Poldata *pd)
{
    /* The first time around we add zeta for the core and addShells will
     * take care of the zeta for the shells.
     * For later calls during optimization of zeta also the
     * zeta on the shells will be set. 
     */
    auto eqtModel = pd->chargeType();
    auto eem      = pd->findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
    for (auto i = 0; i < atoms_->nr; i++)
    {
        auto atype = pd->findAtype(*atoms_->atomtype[i]);
        auto ztype = atype->id(eitELECTRONEGATIVITYEQUALIZATION);
        auto eep   = eem.findParametersConst(ztype);
        auto zeta  = eep["zeta"].value();
        
        if (zeta == 0 && eqtModel != eqtPoint)
        {
            return immZeroZeta;
        }
        
        atoms_->atom[i].zetaA     =
            atoms_->atom[i].zetaB = zeta;
        atoms_->atom[i].row       = eep["row"].value();
    }
    return immOK;
}

immStatus MyMol::GenerateTopology(const Poldata     *pd,
                                  const std::string &method,
                                  const std::string &basis,
                                  std::string       *mylot,
                                  bool               bUseVsites,
                                  bool               bPairs,
                                  bool               bDih,
                                  bool               bBASTAT,
                                  const char        *tabfn)
{
    immStatus   imm = immOK;
    std::string btype1, btype2;

    if (nullptr != debug)
    {
        fprintf(debug, "Generating topology for %s\n", getMolname().c_str());
    }
    nexcl_ = pd->getNexcl();
    GenerateComposition();
    if (NAtom() <= 0)
    {
        imm = immAtomTypes;
    }
    if (immOK == imm)
    {
        snew(atoms_, 1);
        state_change_natoms(state_, NAtom());
        imm = GenerateAtoms(pd, method, basis, mylot);
    }
    if (immOK == imm)
    {
        imm = checkAtoms(pd);
    }
    if (immOK == imm)
    {
        imm = zeta2atoms(pd);
    }
    /* Store bonds in harmonic potential list first, update type later */
    if (immOK == imm)
    {
        int  ftb = F_BONDS;
        auto fs  = pd->findForcesConst(eitBONDS);
        for (auto &bi : bondConst())
        {
            t_param b;
            memset(&b, 0, sizeof(b));
            b.a[0] = bi.getAi() - 1;
            b.a[1] = bi.getAj() - 1;
            pd->atypeToBtype(*atoms_->atomtype[b.a[0]], &btype1);
            pd->atypeToBtype(*atoms_->atomtype[b.a[1]], &btype2);
            Identifier bondId({btype1, btype2}, bi.getBondOrder(), CanSwap::Yes);
            // Store the bond order for later usage.
            bondOrder_.insert({std::make_pair(b.a[0],b.a[1]), bi.getBondOrder()});
            // We add the parameter with zero parameters, they will be
            // set further down. However it is important to set the
            // bondorder.
            add_param_to_plist(plist_, ftb, eitBONDS, b, bi.getBondOrder());
        }
        auto pw = SearchPlist(plist_, ftb);
        if (plist_.end() == pw || pw->nParam() == 0)
        {
            imm = immGenBonds;
        }
    }
    if (immOK == imm)
    {
        MakeAngles(bPairs, bDih);

        MakeSpecialInteractions(pd, bUseVsites);
        imm = updatePlist(pd, &plist_, atoms_, bBASTAT, getMolname(), &error_messages_);
    }
    if (immOK == imm)
    {
        /* Center of charge */
        auto atntot = 0;
        for (auto i = 0; i < atoms_->nr; i++)
        {
            auto atn = atoms_->atom[i].atomnumber;
            atntot  += atn;
            for (auto m = 0; m < DIM; m++)
            {
                coc_[m] += state_->x[i][m]*atn;
            }
        }
        svmul((1.0/atntot), coc_, coc_);
        /* Center of charge */

        bool bAddShells = pd->polarizable();
        if (bAddShells)
        {
            addShells(pd);
        }
        char **molnameptr = put_symtab(symtab_, getMolname().c_str());
        // Generate mtop
        mtop_ = do_init_mtop(pd, molnameptr, atoms_, plist_,
                             inputrec_, symtab_, tabfn);
        excls_to_blocka(atoms_->nr, excls_, &(mtop_->moltype[0].excls));
        if (bAddShells)
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
    if (immOK == imm && !bBASTAT)
    {
        UpdateIdef(pd, eitBONDS);
        UpdateIdef(pd, eitANGLES);
        UpdateIdef(pd, eitIMPROPER_DIHEDRALS);
        if (bDih)
        {
            UpdateIdef(pd, eitPROPER_DIHEDRALS);
        }
    }
    return imm;
}

void MyMol::addShells(const Poldata *pd)
{
    int                    shell  = 0;
    int                    nshell = 0;
    std::vector<int>       renum;
    std::map<int, int>     inv_renum;
    std::vector<std::string> newname;
    t_atoms               *newatoms;
    t_excls               *newexcls;
    std::vector<gmx::RVec> newx;
    t_param                p;

    /* Calculate the total number of Atom and Vsite particles and
     * generate the renumbering array.
     */
    renum.resize(atoms_->nr, 0);
    for (int i = 0; i < atoms_->nr; i++)
    {
        auto atype          = pd->findAtype(*atoms_->atomtype[i]);
        renum[i]            = i + nshell;
        inv_renum[renum[i]] = i;
        if (atype->hasId(eitPOLARIZATION))
        {
            nshell++;
        }
    }
    int nParticles = atoms_->nr+nshell;
    state_change_natoms(state_, nParticles);

    /* Add Polarization to the plist. */
    memset(&p, 0, sizeof(p));
    auto eem = pd->findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
    for (int i = 0; i < atoms_->nr; i++)
    {
        std::string atomtype(*atoms_->atomtype[i]);
        if (atoms_->atom[i].ptype == eptAtom ||
            atoms_->atom[i].ptype == eptVSite)
        {
            auto fa = pd->findAtype(atomtype);
            if (pd->getAtypeEnd() != fa && fa->hasId(eitPOLARIZATION))
            {
                auto ptype = fa->id(eitPOLARIZATION);
                auto param = pd->findForcesConst(eitPOLARIZATION).findParameterTypeConst(ptype, "alpha");
                auto pol   = convertToGromacs(param.value(), param.unit());
                if (pol > 0)
                {
                    p.a[0] = renum[i];
                    p.a[1] = renum[i]+1;
                    if (bHaveVSites_)
                    {
                        auto vsite = pd->findVsite(atomtype);
                        if (vsite != pd->getVsiteEnd())
                        {
                            pol /= vsite->nvsite();
                        }
                    }
                    p.c[0] = pol;
                    add_param_to_plist(plist_, F_POLARIZATION, eitPOLARIZATION, p);
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
    newatoms->nres = atoms_->nres;
    newx.resize(newatoms->nr);
    newname.resize(newatoms->nr);

    /* Make a new exclusion array and put the shells in it. */
    snew(newexcls, newatoms->nr);
    /* Add exclusion for F_POLARIZATION. */
    auto pw = SearchPlist(plist_, F_POLARIZATION);
    if (plist_.end() != pw)
    {
        /* Exclude the vsites and the atoms from their own shell. */
        if (nexcl_ >= 0)
        {
            for (auto j = pw->beginParam(); (j < pw->endParam()); ++j)
            {
                add_excl_pair(newexcls, j->a[0], j->a[1]);
            }
        }

        // Make a copy of the exclusions of the Atom or Vsite for the shell.
        for (auto j = pw->beginParam(); (j < pw->endParam()); ++j)
        {
            // We know that the Atom or Vsite is 0 as we added it to plist as such.
            int  i0 = inv_renum[j->a[0]];
            for (auto j0 = 0; j0 < excls_[i0].nr; j0++)
            {
                add_excl_pair(newexcls, j->a[0], renum[excls_[i0].e[j0]]);
                add_excl_pair(newexcls, j->a[1], renum[excls_[i0].e[j0]]);
            }
        }
        for (auto j = pw->beginParam(); j < pw->endParam(); ++j)
        {
            for (auto j0 = 0; j0 < newexcls[j->a[0]].nr; j0++)
            {
                add_excl_pair(newexcls, j->a[1], newexcls[j->a[0]].e[j0]);
            }
        }
    }

    /* Copy the old atoms to the new structures. */
    for (int i = 0; i < atoms_->nr; i++)
    {
        newatoms->atom[renum[i]]      = atoms_->atom[i];
        newatoms->atomname[renum[i]]  = put_symtab(symtab_, *atoms_->atomname[i]);
        newatoms->atomtype[renum[i]]  = put_symtab(symtab_, *atoms_->atomtype[i]);
        newatoms->atomtypeB[renum[i]] = put_symtab(symtab_, *atoms_->atomtypeB[i]);
        copy_rvec(state_->x[i], newx[renum[i]]);
        newname[renum[i]].assign(*atoms_->atomtype[i]);
        int resind = atoms_->atom[i].resind;
        t_atoms_set_resinfo(newatoms, renum[i], symtab_,
                            *atoms_->resinfo[resind].name,
                            atoms_->resinfo[resind].nr,
                            atoms_->resinfo[resind].ic, 
                            atoms_->resinfo[resind].chainnum, 
                            atoms_->resinfo[resind].chainid);
    }
    t_atom shell_atom = { 0 };
    for (int i = 0; i < atoms_->nr; i++)
    {
        if (atoms_->atom[i].ptype == eptAtom ||
            atoms_->atom[i].ptype == eptVSite)
        {
            std::string atomtype;
            // Shell sits next to the Atom or Vsite
            auto        j            = 1+renum[i];
            auto        atomtypeName = get_atomtype_name(atoms_->atom[i].type, gromppAtomtype_);
            auto fa                  = pd->findAtype(atomtypeName);
            auto ztype               = fa->id(eitELECTRONEGATIVITYEQUALIZATION);
            auto eep                 = eem.findParametersConst(ztype);
            auto shellid             = fa->id(eitPOLARIZATION);
            auto shelltype           = pd->findAtype(shellid.id());
            auto shellzetaid         = shelltype->id(eitELECTRONEGATIVITYEQUALIZATION);
            auto shelleep            = eem.findParametersConst(shellzetaid);
            // Now fill the newatom
            newatoms->atom[j]               = atoms_->atom[i];
            newatoms->atom[j].m             =
                newatoms->atom[j].mB            = shelltype->mass();
            // Shell has no core
            newatoms->atom[j].atomnumber    = 0;
            shell                           = add_atomtype(gromppAtomtype_, symtab_, &shell_atom, shellid.id().c_str(), &p, 0, 0);
            newatoms->atom[j].type          = shell;
            newatoms->atom[j].typeB         = shell;
            newatoms->atomtype[j]           = put_symtab(symtab_, shellid.id().c_str());
            newatoms->atomtypeB[j]          = put_symtab(symtab_, shellid.id().c_str());
            newatoms->atomname[j]           = put_symtab(symtab_, atomtypeName);
            newatoms->atom[j].ptype         = eptShell;
            newatoms->atom[j].zetaA         = shelleep["zeta"].value();
            newatoms->atom[j].zetaB         = newatoms->atom[j].zetaA;
            newatoms->atom[j].row           = shelleep["row"].value();
            newatoms->atom[j].resind        = atoms_->atom[i].resind;
            copy_rvec(state_->x[i], newx[j]);

            newatoms->atom[j].q      =
                newatoms->atom[j].qB = shelleep["charge"].value();
            if (bHaveVSites_)
            {
                if (atoms_->atom[i].ptype == eptVSite)
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
    copy_atoms(newatoms, atoms_);

    for (int i = 0; i < newatoms->nr; i++)
    {
        copy_rvec(newx[i], state_->x[i]);
        atoms_->atomtype[i] = 
            atoms_->atomtypeB[i] = put_symtab(symtab_, *newatoms->atomtype[i]);
    }

    /* Update the bond orders */
    auto  boCopy = bondOrder_;
    bondOrder_.clear();
    for(const auto &boc : boCopy)
    {
        bondOrder_.insert({std::make_pair(renum[boc.first.first],
                                          renum[boc.first.second]), boc.second});
    }
    /* Copy exclusions, empty the original first */
    sfree(excls_);
    excls_ = newexcls;

    /*Now renumber atoms in all other plist interaction types */
    for (auto i = plist_.begin(); i < plist_.end(); ++i)
    {
        if (i->getFtype() != F_POLARIZATION)
        {
            for (auto j = i->beginParam(); j < i->endParam(); ++j)
            {
                for (int k = 0; k < NRAL(i->getFtype()); k++)
                {
                    j->a[k] = renum[j->a[k]];
                }
            }
        }
    }
    bHaveShells_ = true;
}

immStatus MyMol::GenerateGromacs(const gmx::MDLogger       &mdlog,
                                 t_commrec                 *cr,
                                 const char                *tabfn,
                                 gmx_hw_info_t             *hwinfo,
                                 ChargeType                 ieqt)
{
    if (gromacsGenerated_)
    {
        return immOK;
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
    setup_bonded_threading(fr_->bondedThreading, atoms_->nr, false, ltop_->idef);
    wcycle_    = wallcycle_init(debug, 0, cr);

    MDatoms_  = new std::unique_ptr<gmx::MDAtoms>(new gmx::MDAtoms());
    *MDatoms_ = gmx::makeMDAtoms(nullptr, *mtop_, *inputrec_, false);
    atoms2md(mtop_, inputrec_, -1, nullptr, atoms_->nr, MDatoms_->get());
    auto mdatoms = MDatoms_->get()->mdatoms();
    f_.resizeWithPadding(state_->natoms);

    if (nullptr != shellfc_)
    {
        make_local_shells(cr, mdatoms, shellfc_);
    }
    if (eqtSlater != ieqt)
    {
        for (auto i = 0; i < mtop_->natoms; i++)
        {
            mdatoms->row[i] = 0;
        }
    }
    gromacsGenerated_ = true;
    return immOK;
}

immStatus MyMol::computeForces(FILE *fplog, t_commrec *cr, double *rmsf)
{
    auto mdatoms = MDatoms_->get()->mdatoms();
    for (auto i = 0; i < mtop_->natoms; i++)
    {
        mdatoms->chargeA[i] = mtop_->moltype[0].atoms.atom[i].q;
        mdatoms->typeA[i]   = mtop_->moltype[0].atoms.atom[i].type;
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
    immStatus imm =  immOK;
    if (nullptr != shellfc_)
    {
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
            imm = immShellMinimization;
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
            imm = immShellMinimization;
        }
        cr->nnodes = nnodes;
    }
    else
    {
        do_force(fplog, cr, nullptr, inputrec_, nullptr, nullptr, 0,
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
        ConstPlistWrapperIterator bonds = SearchPlist(plist_, eitBONDS);
        if (plist_.end() != bonds)
        {
            symmetrize_charges(bSymmetricCharges, atoms_, bonds,
                               pd, symm_string, &symmetric_charges_);
        }
    }
    else
    {
        for (auto i = 0; i < atoms_->nr; i++)
        {
            symmetric_charges_.push_back(i);
        }
    }
}

void MyMol::initQgenResp(const Poldata     *pd,
                         const std::string &method,
                         const std::string &basis,
                         std::string       *mylot,
                         real              watoms,
                         int               maxESP)
{
    auto iChargeType = pd->chargeType();

    GMX_RELEASE_ASSERT(QgenResp_ == nullptr,
                       "QgenResp_ already initialized");
    QgenResp_ = new QgenResp();
    QgenResp_->setChargeType(iChargeType);
    QgenResp_->setAtomWeight(watoms);
    QgenResp_->setAtomInfo(atoms_, pd, x(), getCharge());
    QgenResp_->setAtomSymmetry(symmetric_charges_);
    QgenResp_->setMolecularCharge(getCharge());
    QgenResp_->summary(debug);

    std::random_device               rd;
    std::mt19937                     gen(rd());  
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    double                           cutoff = 0.01*maxESP;
 
    auto ci = getCalcPropType(method, basis, mylot, MPO_POTENTIAL, nullptr);
    if (ci != EndExperiment())
    {
        int iesp = 0;
        for (auto &epi : ci->electrostaticPotentialConst())
        {
            auto val = uniform(gen);
            if (QgenResp_->myWeight(iesp) > 0 && val < cutoff)
            {
                auto xu = epi.getXYZunit();
                auto vu = epi.getVunit();
                QgenResp_->addEspPoint(convertToGromacs(epi.getX(), xu),
                                       convertToGromacs(epi.getY(), xu),
                                       convertToGromacs(epi.getZ(), xu),
                                       convertToGromacs(epi.getV(), vu));
            }
            iesp++;
        }
        if (debug)
        {
            fprintf(debug, "Added %zu ESP points to the RESP structure.\n", QgenResp_->nEsp());
        }
    }
}

immStatus MyMol::GenerateCharges(const Poldata       *pd,
                                 const gmx::MDLogger &mdlog,
                                 t_commrec           *cr,
                                 const char          *tabfn,
                                 gmx_hw_info_t       *hwinfo,
                                 int                  maxiter,
                                 real                 tolerance)
{
    immStatus           imm          = immOK;
    bool                converged    = false;
    int                 iter         = 0;
    auto                iChargeType  = pd->chargeType();

    GenerateGromacs(mdlog, cr, tabfn, hwinfo, iChargeType);
    if (backupCoordinates_.size() == 0)
    {
        backupCoordinates();
    }
    switch (pd->chargeGenerationAlgorithm())
    {
    case eqgNONE:
        if (debug)
        {
            fprintf(debug, "WARNING! Using zero charges for %s!\n",
                    getMolname().c_str());
        }
        for (auto i = 0; i < atoms_->nr; i++)
        {
            atoms_->atom[i].q  = atoms_->atom[i].qB = 0;
        }
        return immOK;
    case eqgESP:
        {
            double chi2[2]   = {1e8, 1e8};
            real   rrms      = 0;
            int    cur       = 0;
            real   cosangle  = 0;
            EspRms_          = 0;
            iter             = 0;
            
            // Init Qgresp should be called before this!
            QgenResp_->optimizeCharges(pd->getEpsilonR());
            QgenResp_->calcPot(pd->getEpsilonR());
            EspRms_ = chi2[cur] = QgenResp_->getRms(&rrms, &cosangle);
            if (debug)
            {
                fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
            }
            do
            {
                for (auto i = 0; i < atoms_->nr; i++)
                {
                    mtop_->moltype[0].atoms.atom[i].q      =
                        mtop_->moltype[0].atoms.atom[i].qB = QgenResp_->getCharge(i);
                }
                if (nullptr != shellfc_)
                {
                    double rmsf;
                    auto imm = computeForces(nullptr, cr, &rmsf);
                    if (imm != immOK)
                    {
                        return imm;
                    }
                    QgenResp_->updateAtomCoords(state_->x);
                }
                QgenResp_->optimizeCharges(pd->getEpsilonR());
                QgenResp_->calcPot(pd->getEpsilonR());
                real cosangle = 0;
                EspRms_ = chi2[cur] = QgenResp_->getRms(&rrms, &cosangle);
                if (debug)
                {
                    fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
                }
                converged = (fabs(chi2[cur] - chi2[1-cur]) < tolerance) || (nullptr == shellfc_);
                cur       = 1-cur;
                iter++;
            }
            while ((!converged) && (iter < maxiter));
            for (auto i = 0; i < atoms_->nr; i++)
            {
                atoms_->atom[i].q      =
                    atoms_->atom[i].qB = QgenResp_->getCharge(i);
            }
        }
        break;
    case eqgEEM:
    case eqgSQE:
        {
            if (QgenAcm_ == nullptr)
            {
                QgenAcm_ = new QgenAcm(pd, atoms_, getCharge());
            }

            auto q                 = QgenAcm_->q();
            std::vector<double> qq = q;
            iter                   = 0;
            do
            {
                if (eQGEN_OK == QgenAcm_->generateCharges(debug,
                                                          getMolname().c_str(),
                                                          pd,
                                                          atoms_,
                                                          state_->x,
                                                          bonds()))
                {
                    for (auto i = 0; i < mtop_->natoms; i++)
                    {
                        mtop_->moltype[0].atoms.atom[i].q      =
                            mtop_->moltype[0].atoms.atom[i].qB = atoms_->atom[i].q;
                    }
                    if (nullptr != shellfc_)
                    {
                        double rmsf;
                        auto imm = computeForces(nullptr, cr, &rmsf);
                        if (imm != immOK)
                        {
                            return imm;
                        }
                    }
                    EemRms_ = 0;
                    for (size_t i = 0; i < q.size(); i++)
                    {
                        EemRms_  += gmx::square(qq[i] - q[i]);
                        qq[i]     = q[i];
                    }
                    EemRms_  /= q.size();
                    converged = (EemRms_ < tolerance) || (nullptr == shellfc_);
                    iter++;
                }
                else
                {
                    imm = immChargeGeneration;
                }
            }
            while (imm == immOK && (!converged) && (iter < maxiter));
            for (auto i = 0; i < mtop_->natoms; i++)
            {
                mtop_->moltype[0].atoms.atom[i].q      =
                    mtop_->moltype[0].atoms.atom[i].qB = atoms_->atom[i].q;
            }
            if (!converged)
            {
                printf("Alexandria Charge Model did not converge to %g. rms: %g\n", tolerance, sqrt(EemRms_));
            }
        }
        break;
    default:
        gmx_fatal(FARGS, "ChargeGenerationAlgorithm %s is not implemented",
                  chargeGenerationAlgorithmName(pd->chargeGenerationAlgorithm()).c_str());
        break;
    }
    return imm;
}

void MyMol::plotEspCorrelation(const char             *espcorr,
                               const gmx_output_env_t *oenv)
{
    if (espcorr && oenv)
    {
        QgenResp_->updateAtomCharges(atoms_);
        QgenResp_->calcPot(1.0);
        QgenResp_->plotLsq(oenv, espcorr);
    }
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
        if (JOB_OPT == ei.getJobtype())
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

void MyMol::CalcDipole()
{
    rvec mu;
    CalcDipole(mu);
    set_muQM(qtCalc, mu);
}

void MyMol::CalcDipole(rvec mu)
{
    rvec   r; /* distance of atoms to center of charge */
    clear_rvec(mu);
    for (auto i = 0; i < atoms_->nr; i++)
    {
        rvec_sub(state_->x[i], coc_, r);
        auto q = e2d(atoms_->atom[i].q);
        for (auto m = 0; m < DIM; m++)
        {
            mu[m] += r[m]*q;
        }
    }
}

void MyMol::CalcQuadrupole()
{
    real   r2, q;
    rvec   r; /* distance of atoms to center of charge */
    tensor Q;
    clear_mat(Q);
    for (auto i = 0; i < atoms_->nr; i++)
    {
        rvec_sub(state_->x[i], coc_, r);
        r2   = iprod(r, r);
        q    = atoms_->atom[i].q;
        for (auto m = 0; m < DIM; m++)
        {
            for (auto n = 0; n < DIM; n++)
            {
                Q[m][n] += 0.5*q*(3.0*r[m]*r[n] - r2*delta(m, n))*NM2A*A2CM*CM2D*10;
            }
        }
    }
    set_QQM(qtCalc, Q);
}

void MyMol::CalcQMbasedMoments(real *q, rvec mu, tensor Q)
{
    int   i, j;
    real  r2;
    rvec  r;  /* distance of atoms to center of charge */

    clear_rvec(mu);
    clear_mat(Q);
    for (i = j = 0; i < atoms_->nr; i++)
    {
        if (atoms_->atom[i].ptype == eptAtom ||
            atoms_->atom[i].ptype == eptNucleus)
        {
            rvec_sub(state_->x[i], coc_, r);
            r2       = iprod(r, r);
            for (auto m = 0; m < DIM; m++)
            {
                mu[m] += (r[m]*e2d(q[j]));
                for (auto n = 0; n < DIM; n++)
                {
                    Q[m][n] += 0.5*q[j]*(3.0*r[m]*r[n] - r2*delta(m, n))*NM2A*A2CM*CM2D*10;
                }
            }
            j++;
        }
    }
}

void MyMol::CalcQPol(const Poldata *pd, rvec mu)

{
    int     np;
    double  poltot, sptot, ereftot, eref;

    poltot  = 0;
    sptot   = 0;
    ereftot = 0;
    np      = 0;
    auto eep = pd->findForcesConst(eitPOLARIZATION);
    for (int i = 0; i < atoms_->nr; i++)
    {
        std::string ptype;
        auto atype = pd->findAtype(*atoms_->atomtype[i]);
        auto idP   = atype->id(eitPOLARIZATION);
        if (eep.parameterExists(idP))
        {
            auto param  = eep.findParameterTypeConst(idP, "alpha");
            poltot += param.value();
            sptot  += gmx::square(param.uncertainty());
            np++;
        }
        if (pd->getAtypeRefEnthalpy(*atoms_->atomtype[i], &eref))
        {
            ereftot += eref;
        }
    }
    CalcDipole(mu);
    CalcQuadrupole();
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
    GMX_RELEASE_ASSERT(backupCoordinates_.size() == 0, "Can only backup coordinates once");
    backupCoordinates_.resize(atoms_->nr);
    for(int i = 0; i < atoms_->nr; i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            backupCoordinates_[i][m] = state_->x[i][m];
        }
    }
}

void MyMol::restoreCoordinates()
{
    if (static_cast<int>(backupCoordinates_.size()) == atoms_->nr)
    {
        for(int i = 0; i < atoms_->nr; i++)
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
    immStatus           imm = immOK;
    double              rmsf;

    field.resize(DIM, 0);
    myforce_->setField(field);
    imm          = computeForces(fplog, cr, &rmsf);
    isoPol_calc_ = 0;
    CalcDipole(mu_ref);
    for (auto m = 0; imm == immOK && m < DIM; m++)
    {
        field[m] = efield;
        myforce_->setField(field);
        imm = computeForces(fplog, cr, &rmsf);
        field[m] = 0;
        myforce_->setField(field);
        if (imm == immOK)
        {
            rvec mu_tot;
            CalcDipole(mu_tot);
            for (auto n = 0; n < DIM; n++)
            {
                alpha_calc_[n][m] = ((mu_tot[n]-mu_ref[n])/efield)*(POLFAC);
            }
            isoPol_calc_ += alpha_calc_[m][m]/DIM;
        }
    }
    if (immOK == imm)
    {
        CalcAnisoPolarizability(alpha_calc_, &anisoPol_calc_);
    }
    restoreCoordinates();
    return imm;
}

void MyMol::PrintConformation(const char *fn)
{
    char title[STRLEN];

    put_in_box(atoms_->nr, state_->box, as_rvec_array(state_->x.data()), 0.3);
    sprintf(title, "%s processed by alexandria", getMolname().c_str());
    write_sto_conf(fn, title, atoms_, as_rvec_array(state_->x.data()), nullptr, epbcNONE, state_->box);
}

void MyMol::PrintTopology(const char        *fn,
                          bool               bVerbose,
                          const Poldata     *pd,
                          gmx_atomprop_t     aps,
                          t_commrec         *cr,
                          double             efield,
                          const std::string &method,
                          const std::string &basis)
{
    FILE  *fp   = gmx_ffopen(fn, "w");
    bool   bITP = (fn2ftp(fn) == efITP);

    PrintTopology(fp, bVerbose, pd, aps, bITP, cr, efield, method, basis);

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

static void rotate_tensor(tensor Q, tensor Qreference)
{
    matrix rotmatrix;
    rvec   tmpvec;
    // TODO: this code is not correct!
    // the whole tensor should be taken into account, not
    // just the components. All vectors should be transformed
    // by the same matrix.
    for (int m = 0; m < DIM; m++)
    {
        if (norm(Q[m]) > 0 && norm(Qreference[m]) > 0)
        {
            calc_rotmatrix(Q[m], Qreference[m], rotmatrix);
            mvmul(rotmatrix, Q[m], tmpvec);
            copy_rvec(tmpvec, Q[m]);
        }
    }
}

void MyMol::PrintTopology(FILE                   *fp,
                          bool                    bVerbose,
                          const Poldata          *pd,
                          gmx_atomprop_t          aps,
                          bool                    bITP,
                          t_commrec              *cr,
                          double                  efield,
                          const std::string      &method,
                          const std::string      &basis)
{
    char                     buf[256];
    t_mols                   printmol;
    std::vector<std::string> commercials;
    rvec                     vec, mu;
    tensor                   myQ;
    double                   value = 0, error = 0, T = -1;
    std::string              myref;
    auto                     iChargeType = pd->chargeType();
    std::string              mylot       = makeLot(method, basis);

    if (fp == nullptr)
    {
        return;
    }

    CalcQPol(pd, mu);
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
    snprintf(buf, sizeof(buf), "Total Charge = %d (e)", getCharge());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Charge Type  = %s\n",
             chargeTypeName(iChargeType).c_str());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Alexandria Dipole Moment (Debye):\n"
             "; ( %.2f %6.2f %6.2f ) Total= %.2f\n",
             mu[XX], mu[YY], mu[ZZ],
             norm(mu));
    commercials.push_back(buf);

    T = -1;
    const char *qm_type = "electronic";
    const char *qm_conf = "minimum";
    if (getPropRef(MPO_DIPOLE, iqmQM, method, basis, qm_conf,
                   qm_type, &value, &error,
                   &T, &myref, &mylot, vec, myQ))
    {
        set_muQM(qtElec, vec);
        if (value > 0)
        {
            rotateDipole(mu_qm_[qtElec], mu);
        }
        snprintf(buf, sizeof(buf), "%s Dipole Moment (Debye):\n"
                 "; ( %.2f %6.2f %6.2f ) Total= %.2f\n",
                 mylot.c_str(),
                 mu_qm_[qtElec][XX], mu_qm_[qtElec][YY], mu_qm_[qtElec][ZZ],
                 norm(mu_qm_[qtElec]));
        commercials.push_back(buf);
    }
    else if (bVerbose)
    {
        printf("WARNING: QM dipole of type %s not found for lot %s\n",
               qm_type, mylot.c_str());
    }

    add_tensor(&commercials,
               "Alexandria Traceless Quadrupole Moments (Buckingham)",
               QQM(qtCalc));

    T = -1;
    if (getPropRef(MPO_QUADRUPOLE, iqmQM, method, basis, qm_conf,
                   qm_type, &value, &error,
                   &T, &myref, &mylot, vec, myQ))
    {
        set_QQM(qtElec, myQ);
        rotate_tensor(Q_qm_[qtElec], Q_qm_[qtCalc]);
        snprintf(buf, sizeof(buf), "%s Traceless Quadrupole Moments (Buckingham)", mylot.c_str());
        add_tensor(&commercials, buf, Q_qm_[qtElec]);
    }
    else if (bVerbose)
    {
        printf("WARNING: QM quadrupole of type %s not found for lot %s\n",
               qm_type, mylot.c_str());
    }

    snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability (Additive Law): %.2f +/- %.2f (A^3)\n", polarizability_, sig_pol_);
    commercials.push_back(buf);

    if (efield > 0 && nullptr != cr)
    {
        auto imm = CalcPolarizability(efield, cr, debug);
        if (imm == immOK)
        {
            add_tensor(&commercials, "Alexandria Polarizability components (A^3)", alpha_calc_);
            
            snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability (Interactive): %.2f (A^3)\n", isoPol_calc_);
            commercials.push_back(buf);
            
            snprintf(buf, sizeof(buf), "Alexandria Anisotropic Polarizability: %.2f (A^3)\n", anisoPol_calc_);
            commercials.push_back(buf);

            T = -1;
            if (getPropRef(MPO_POLARIZABILITY, iqmQM, method, basis, "",
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

    print_top_header(fp, pd, aps, bHaveShells_, commercials, bITP);
    write_top(fp, printmol.name, atoms_, false,
              plist_, excls_, gromppAtomtype_, cgnr_, nexcl_, pd);
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
    if ((cgnr_ = generate_charge_groups(ecg, atoms_, plist_, bUsePDBcharge,
                                        &qtot, &mtot)) == nullptr)
    {
        return immChargeGeneration;
    }
    if (ecg != ecgAtom)
    {
        //sort_on_charge_groups(cgnr_, atoms_,
        //                    plist_, x_, excls_, ndxfn, nmol);
    }
    return immOK;
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
    ChargeType iChargeType = pd->chargeType();

    if (potfn || hisfn || rhofn || difffn || pdbdifffn)
    {
        char     *gentop_version = (char *)"gentop v0.99b";
        QgenResp *grref;

        QgenResp_->updateAtomCharges(atoms_);
        QgenResp_->calcPot(pd->getEpsilonR());
        QgenResp_->potcomp(pcfn, atoms_, 
                        as_rvec_array(state_->x.data()), pdbdifffn, oenv);

        /* This has to be done before the grid is f*cked up by
           writing a cube file */
        grref = new QgenResp(QgenResp_);

        if (reffn)
        {
            grref->setAtomInfo(atoms_, pd, state_->x, getCharge());
            grref->setAtomSymmetry(symmetric_charges_);
            grref->readCube(reffn, FALSE);
            delete QgenResp_; 
            QgenResp_ = new QgenResp(grref);
        }
        else
        {
            QgenResp_->makeGrid(spacing, border, as_rvec_array(state_->x.data()));
        }
        if (rhofn)
        {
            std::string buf = gmx::formatString("Electron density generated by %s based on %s charges",
                                                gentop_version,
                                                chargeTypeName(iChargeType).c_str());
            QgenResp_->calcRho();
            QgenResp_->writeRho(rhofn, buf, oenv);
        }
        if (potfn)
        {
            std::string buf = gmx::formatString("Potential generated by %s based on %s charges",
                                                gentop_version,
                                                chargeTypeName(iChargeType).c_str());
            QgenResp_->calcPot(pd->getEpsilonR());
            QgenResp_->writeCube(potfn, buf, oenv);
        }
        if (hisfn)
        {
            std::string buf = gmx::formatString("Potential generated by %s based on %s charges",
                                                gentop_version,
                                                chargeTypeName(iChargeType).c_str());
            QgenResp_->writeHisto(hisfn, buf, oenv);
        }
        if (difffn || diffhistfn)
        {
            std::string buf = gmx::formatString("Potential difference generated by %s based on %s charges",
                                                gentop_version,
                                                chargeTypeName(iChargeType).c_str());

            QgenResp_->writeDiffCube(grref, difffn, diffhistfn, buf, oenv, 0);
        }
    }
}

void MyMol::rotateDipole(rvec mu, rvec muReference)
{
    if (norm(mu) < 0.04 or norm(muReference) < 0.04)
    {
        return;
    }
    matrix rotmatrix;
    rvec   tmpvec;
    calc_rotmatrix(mu, muReference, rotmatrix);
    mvmul(rotmatrix, mu, tmpvec);
    copy_rvec(tmpvec, mu);
}

void MyMol::setQandMoments(qType qt, int natom, real q[])
{
    int i, j;

    if (natom > 0)
    {
        charge_QM_[qt].resize(natom);
        for (i = j = 0; i < atoms_->nr; i++)
        {
            if (atoms_->atom[i].ptype == eptAtom ||
                atoms_->atom[i].ptype == eptNucleus)
            {
                charge_QM_[qt][j] = q[j];
                j++;
            }
        }
        CalcQMbasedMoments(q, mu_qm_[qt], Q_qm_[qt]);
    }
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
    double       Hatom = 0;
    double       ZPE   = 0;
    double       error = 0;
    double       T     = -1;
    rvec         vec;
    tensor       quadrupole = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    tensor       polar      = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::string  myref;
    std::string  mylot;
    bool         esp_dipole_found  = false;

    for (auto i = 0; i < atoms_->nr; i++)
    {
        if (atoms_->atom[i].ptype == eptAtom ||
            atoms_->atom[i].ptype == eptNucleus)
        {
            natom++;
        }
    }
    real q[natom];
    if (getPropRef(MPO_CHARGE, iqmQM,
                              method, basis, "",
                              (char *)"ESP charges",
                              &value, &error, &T,
                              &myref, &mylot, q, quadrupole))
    {
        setQandMoments(qtESP, natom, q);
        esp_dipole_found = true;
    }
    T = -1;
    if (getPropRef(MPO_CHARGE, iqmQM,
                              method, basis, "",
                              (char *)"Mulliken charges",
                              &value, &error, &T,
                              &myref, &mylot, q, quadrupole))
    {
        setQandMoments(qtMulliken, natom, q);
        if (esp_dipole_found && dipQM(qtMulliken) > 0)
        {
            rotate_tensor(Q_qm_[qtMulliken], Q_qm_[qtESP]);
        }
    }
    T = -1;
    if (getPropRef(MPO_CHARGE, iqmQM,
                              method, basis, "",
                              (char *)"Hirshfeld charges",
                              &value, &error, &T,
                              &myref, &mylot, q, quadrupole))
    {
        setQandMoments(qtHirshfeld, natom, q);
        if (esp_dipole_found && dipQM(qtHirshfeld) > 0)
        {
            rotate_tensor(Q_qm_[qtHirshfeld], Q_qm_[qtESP]);
        }

    }
    T = -1;
    if (getPropRef(MPO_CHARGE, iqmQM,
                              method, basis, "",
                              (char *)"CM5 charges",
                              &value, &error, &T,
                              &myref, &mylot, q, quadrupole))
    {
        setQandMoments(qtCM5, natom, q);
        if (esp_dipole_found && dipQM(qtCM5) > 0)
        {
            rotate_tensor(Q_qm_[qtCM5], Q_qm_[qtESP]);
        }

    }
    T = 298.15;
    immStatus imm = immOK;
    if (bDHform &&
        getProp(MPO_ENERGY, (bQM ? iqmQM : iqmExp),
                           method, basis, "",
                           (char *)"DeltaHform", &value, &error, &T))
    {
        Hform_ = value;
        Emol_  = value;
        for (ia = 0; ia < atoms_->nr; ia++)
        {
            if (atoms_->atom[ia].ptype == eptAtom ||
                atoms_->atom[ia].ptype == eptNucleus)
            {
                if (pd->getAtypeRefEnthalpy(*atoms_->atomtype[ia], &Hatom))
                {
                    Emol_ -= Hatom;
                }
                else
                {
                    fprintf(stderr, "WARNING: NO reference enthalpy for molecule %s.\n",
                            getMolname().c_str());
                    Emol_ = 0;
                    imm   = immNoData;
                    break;
                }
            }
        }
        if (bZPE)
        {

            if (getProp(MPO_ENERGY, iqmBoth,
                                   method, basis, "",
                                   (char *)"ZPE", &ZPE, &error, &T))
            {
                Emol_ -= ZPE;
            }
            else
            {
                fprintf(stderr, "No zero-point energy for molecule %s.\n",
                        getMolname().c_str());
                imm = immNoData;
            }
        }
        if (ia < atoms_->nr)
        {
            imm = immNoData;
        }
    }
    if (imm == immOK)
    {
        T = -1;
        if (getPropRef(MPO_DIPOLE, iqmQM,
                                  method, basis, "",
                                  (char *)"electronic",
                                  &value, &error, &T, &myref, &mylot,
                                  vec, quadrupole))
        {
            dip_exp_  = value;
            dip_err_  = error;
            set_muQM(qtElec, vec);

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

            if (!bZero && (dipQM(qtElec) - 0.0) < 1e-2)
            {
                imm = immZeroDip;
            }
            if (immOK == imm && esp_dipole_found)
            {
                rotateDipole(mu_qm_[qtElec], mu_qm_[qtESP]);
            }
        }
        else
        {
            imm = immNoDipole;
        }
    }
    if (immOK == imm)
    {
        T = -1;
        if (getPropRef(MPO_QUADRUPOLE, iqmQM,
                       method, basis, "",
                       (char *)"electronic",
                       &value, &error, &T, &myref, &mylot,
                       vec, quadrupole))
        {
            set_QQM(qtElec, quadrupole);
            if (immOK == imm && esp_dipole_found && norm(mu_qm_[qtElec]) > 0)
            {
                rotate_tensor(Q_qm_[qtElec], Q_qm_[qtESP]);
            }
        }
        T = -1;
        if (getPropRef(MPO_POLARIZABILITY, iqmQM,
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
    std::vector<std::string> batoms;
    for (int j = 0; j < natoms; j++)
    {
        batoms.push_back(btype[iatoms[j]]);
    }
    auto fs = pd->findForcesConst(iType);
    if (iType == eitBONDS && natoms == 2)
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
    std::string params;

    if (debug)
    {
        fprintf(debug, "UpdateIdef for %s\n", interactionTypeToString(iType).c_str());
    }
    if (iType == eitVDW)
    {
        nonbondedFromPdToMtop(mtop_, atoms_, pd, fr_);
        if (debug)
        {
            pr_ffparams(debug, 0, "UpdateIdef Before", &mtop_->ffparams, false);
        }
    }
    else if (iType == eitPOLARIZATION)
    {
        auto pw = SearchPlist(plist_, iType);
        if (plist_.end() != pw)
        {
            auto ft = pw->getFtype();
            polarizabilityFromPdToMtop(mtop_, ltop_, atoms_, pd, ft);           
        }
    }
    else
    {
        // Update other iTypes
        auto fs    = pd->findForcesConst(iType);
        auto ftype = fs.fType();
        // Small optimization. This assumes angles etc. use the same
        // types as bonds.
        std::vector<std::string> btype;
        for(int j = 0; j < atoms_->nr; j++)
        {
            auto atype = pd->findAtype(*atoms_->atomtype[j]);
            btype.push_back(atype->id(eitBONDS).id());
        }
        for (auto i = 0; i < ltop_->idef.il[ftype].nr; i += interaction_function[ftype].nratoms+1)
        {
            auto  tp     = ltop_->idef.il[ftype].iatoms[i];
            auto  bondId = getIdentifier(pd, iType, btype, 
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
                        convertToGromacs(fp.value(), fp.unit());
                        
                    fp = fs.findParameterTypeConst(bondId, "kt");
                    mtop_->ffparams.iparams[tp].harmonic.krA         =
                        mtop_->ffparams.iparams[tp].harmonic.krB     =
                        ltop_->idef.iparams[tp].harmonic.krA     =
                        ltop_->idef.iparams[tp].harmonic.krB =
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
                        convertToGromacs(fp.value(), fp.unit());
                        
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
            default:
                GMX_THROW(gmx::InternalError(gmx::formatString("Do not know what to do for %s",
                                                               interaction_function[ftype].longname).c_str()));
                break;
            }
        }
    }
}

} // namespace alexandria
