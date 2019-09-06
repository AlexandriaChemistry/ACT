/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019
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

#include <assert.h>

#include <cstdio>
#include <cstring>

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
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
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

#include "mymol_low.h"

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
    atype_         = init_atomtype();
    mtop_          = nullptr;
    fr_            = nullptr;
    ltop_          = nullptr;
    mp_            = new MolProp;
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
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        mm  = topology_->atoms.atom[i].m;
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
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        rvec_dec(state_->x[i], com);
    }

    snew(bSymm, topology_->atoms.nr);
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        bSymm[i] = (norm(state_->x[i]) < toler);
        for (j = i+1; (j < topology_->atoms.nr) && !bSymm[i]; j++)
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
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        bSymmAll = bSymmAll && bSymm[i];
    }
    sfree(bSymm);
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        rvec_inc(state_->x[i], com);
    }

    return bSymmAll;
}

void MyMol::findInPlaneAtoms(int ca, std::vector<int> &atoms)
{
    int bca = 0;
    /*First try to find the atom bound to the central atom (ca).*/
    for (auto bi = molProp()->BeginBond();
         bi < molProp()->EndBond(); bi++)
    {
        if ((ca == (bi->getAj() - 1) ||
             ca == (bi->getAi() - 1)))
        {
            if (ca == (bi->getAi() - 1))
            {
                bca = (bi->getAj() - 1);
                atoms.push_back(bca);
            }
            else
            {
                bca = (bi->getAi() - 1);
                atoms.push_back(bca);
            }
        }
    }
    /*Now try to find atoms bound to bca, except ca.*/
    for (auto bi = molProp()->BeginBond();
         bi < molProp()->EndBond(); bi++)
    {
        if ((ca != (bi->getAj() - 1)   &&
             ca != (bi->getAi() - 1))  &&
            (bca == (bi->getAj() - 1)  ||
             bca == (bi->getAi() - 1)))
        {
            if (bca == (bi->getAi() - 1))
            {
                atoms.push_back(bi->getAj() - 1);
            }
            else
            {
                atoms.push_back(bi->getAi() - 1);
            }
        }
    }
}

void MyMol::findOutPlaneAtoms(int ca, std::vector<int> &atoms)
{
    for (auto bi = molProp()->BeginBond();
         bi < molProp()->EndBond(); bi++)
    {
        if (bi->getBondOrder() == 1  &&
            (ca == (bi->getAj() - 1) ||
             ca == (bi->getAi() - 1)))
        {
            if (ca == (bi->getAi() - 1))
            {
                atoms.push_back(bi->getAj() - 1);
            }
            else
            {
                atoms.push_back(bi->getAi() - 1);
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

    bonds.resize(topology_->atoms.nr);
    for (auto bi = molProp()->BeginBond(); (bi < molProp()->EndBond()); bi++)
    {
        // Store bonds bidirectionally to get the number correct
        bonds[bi->getAi() - 1].push_back(bi->getAj() - 1);
        bonds[bi->getAj() - 1].push_back(bi->getAi() - 1);
    }
    nbonds.resize(topology_->atoms.nr);
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        nbonds[i] = bonds[i].size();
    }
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        /* Now test initial geometry */
        if ((bonds[i].size() == 2) &&
            is_linear(x[i], x[bonds[i][0]], x[bonds[i][1]],
                      &pbc, th_toler))
        {
            if (nullptr != debug)
            {
                fprintf(debug, "found linear angle %s-%s-%s in %s\n",
                        *topology_->atoms.atomtype[bonds[i][0]],
                        *topology_->atoms.atomtype[i],
                        *topology_->atoms.atomtype[bonds[i][1]],
                        molProp()->getMolname().c_str());
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
                        *topology_->atoms.atomtype[i],
                        *topology_->atoms.atomtype[bonds[i][0]],
                        *topology_->atoms.atomtype[bonds[i][1]],
                        *topology_->atoms.atomtype[bonds[i][2]],
                        molProp()->getMolname().c_str());
            }
            gvt_.addPlanar(i, bonds[i][0], bonds[i][1], bonds[i][2],
                           &nbonds[0]);
        }
        if (bUseVsites)
        {
            const auto atype(*topology_->atoms.atomtype[i]);
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
    auto anr = topology_->atoms.nr;
    gvt_.generateSpecial(pd, bUseVsites, &topology_->atoms,
                         &x, plist_, symtab_, atype_, &excls_, state_);
    bHaveVSites_ = (topology_->atoms.nr > anr);
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
    snew(excls_, topology_->atoms.nr);
    init_nnb(&nnb, topology_->atoms.nr, nexcl_ + 2);
    gen_nnb(&nnb, plist);

    print_nnb(&nnb, "NNB");
    rtp.bKeepAllGeneratedDihedrals    = bDihs;
    rtp.bRemoveDihedralIfWithImproper = bDihs;
    rtp.bGenerateHH14Interactions     = bPairs;
    rtp.nrexcl                        = nexcl_;

    gen_pad(&nnb, &(topology_->atoms), &rtp, plist, excls_, nullptr, false);

    t_blocka *EXCL;
    snew(EXCL, 1);
    if (debug)
    {
        fprintf(debug, "Will generate %d exclusions for %d atoms\n",
                nexcl_, topology_->atoms.nr);
    }
    generate_excl(nexcl_, topology_->atoms.nr, plist, &nnb, EXCL);
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
        for (auto i = 0; i < topology_->atoms.nr; i++)
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

immStatus MyMol::GenerateAtoms(gmx_atomprop_t            ap,
                               const char               *lot)
{
    int                 myunit;
    double              xx, yy, zz;
    int                 natom = 0;
    immStatus           imm   = immOK;

    ExperimentIterator  ci = molProp()->getLot(lot);
    if (ci < molProp()->EndExperiment())
    {
        t_param nb;
        memset(&nb, 0, sizeof(nb));
        natom = 0;
        init_t_atoms(&(topology_->atoms), ci->NAtom(), false);
        snew(topology_->atoms.atomtype, ci->NAtom());
        snew(topology_->atoms.atomtypeB, ci->NAtom());

        for (auto cai = ci->BeginAtom(); (cai < ci->EndAtom()); cai++)
        {
            myunit = string2unit((char *)cai->getUnit().c_str());
            if (myunit == -1)
            {
                gmx_fatal(FARGS, "Unknown length unit '%s' for atom coordinates",
                          cai->getUnit().c_str());
            }
            cai->getCoords(&xx, &yy, &zz);

            state_->x[natom][XX] = convert2gmx(xx, myunit);
            state_->x[natom][YY] = convert2gmx(yy, myunit);
            state_->x[natom][ZZ] = convert2gmx(zz, myunit);

            double q = 0;
            for (auto qi = cai->BeginQ(); (qi < cai->EndQ()); qi++)
            {
                // TODO Clean up this mess.
                if (qi->getType().compare("ESP") == 0) 
                {
                    myunit = string2unit((char *)qi->getUnit().c_str());
                    q      = convert2gmx(qi->getQ(), myunit);
                    break;
                }
            }
            topology_->atoms.atom[natom].q      =
                topology_->atoms.atom[natom].qB = q;

            t_atoms_set_resinfo(&(topology_->atoms), natom, symtab_, ((cai->ResidueName().c_str() != nullptr) ? cai->ResidueName().c_str() : molProp()->getMolname().c_str()), 1, ' ', 1, ' ');
            topology_->atoms.atomname[natom]        = put_symtab(symtab_, cai->getName().c_str());
            topology_->atoms.atom[natom].atomnumber = gmx_atomprop_atomnumber(ap, cai->getName().c_str());

            real mass = 0;
            if (!gmx_atomprop_query(ap, epropMass, "???", cai->getName().c_str(), &mass))
            {
                fprintf(stderr, "Could not find mass for %s\n", cai->getName().c_str());
            }
            topology_->atoms.atom[natom].m      =
                topology_->atoms.atom[natom].mB = mass;

            strcpy(topology_->atoms.atom[natom].elem, gmx_atomprop_element(ap, topology_->atoms.atom[natom].atomnumber));

            topology_->atoms.atom[natom].resind = 0;
            // First set the atomtype
            topology_->atoms.atomtype[natom]      =
                topology_->atoms.atomtypeB[natom] = put_symtab(symtab_, cai->getObtype().c_str());

            natom++;
        }
        for (auto i = 0; i < natom; i++)
        {
            topology_->atoms.atom[i].type      =
                topology_->atoms.atom[i].typeB = add_atomtype(atype_, symtab_,
                                                              &(topology_->atoms.atom[i]),
                                                              *topology_->atoms.atomtype[i],
                                                              &nb, 0,
                                                              topology_->atoms.atom[i].atomnumber);
        }
        topology_->atoms.nr   = natom;
        topology_->atoms.nres = 1;
    }
    else
    {
        imm = immLOT;
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Tried to convert %s to gromacs. LOT is %s. Natoms is %d\n",
                molProp()->getMolname().c_str(), lot, natom);
    }

    return imm;
}

immStatus MyMol::checkAtoms(const Poldata *pd)
{
    auto nmissing = 0;
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        const auto atype(*topology_->atoms.atomtype[i]);
        auto       fa = pd->findAtype(atype);
        if (fa == pd->getAtypeEnd())
        {
            printf("Could not find a force field entry for atomtype %s atom %d in compound '%s'\n",
                   *topology_->atoms.atomtype[i], i+1,
                   MolProp().getMolname().c_str());
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
    /* Here, we add zeta for the core. addShells will
       take care of the zeta for the shells later. */
    auto zeta = 0.0;
    auto row  = 0;
    auto eqdModel = pd->getChargeModel();
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        zeta = pd->getZeta(*topology_->atoms.atomtype[i], 0);
        row  = pd->getRow(*topology_->atoms.atomtype[i], 0);
        if (zeta == 0 && getEemtypeDistributed(eqdModel))
        {
            return immZeroZeta;
        }
        topology_->atoms.atom[i].zetaA     =
            topology_->atoms.atom[i].zetaB = zeta;
        topology_->atoms.atom[i].row       =  row;
    }
    return immOK;
}

immStatus MyMol::GenerateTopology(gmx_atomprop_t  ap,
                                  const Poldata  *pd,
                                  const char     *lot,
                                  bool            bUseVsites,
                                  bool            bPairs,
                                  bool            bDih,
                                  bool            bBASTAT,
                                  const char     *tabfn)
{
    int         ftb;
    t_param     b;
    immStatus   imm = immOK;
    std::string btype1, btype2;
    ChargeModel iChargeModel = pd->getChargeModel();
                                  
    if (nullptr != debug)
    {
        fprintf(debug, "Generating topology for %s\n", molProp()->getMolname().c_str());
    }
    nexcl_ = pd->getNexcl();
    molProp()->GenerateComposition(pd);
    if (molProp()->NAtom() <= 0)
    {
        imm = immAtomTypes;
    }
    if (immOK == imm)
    {
        snew(topology_, 1);
        init_top(topology_);
        state_change_natoms(state_, molProp()->NAtom());
        imm = GenerateAtoms(ap, lot);
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
    ftb = F_BONDS;
    if (immOK == imm)
    {
        for (auto fs = pd->forcesBegin(); fs != pd->forcesEnd(); fs++)
        {
            if (eitBONDS == fs->iType())
            {
                ListedForceConstIterator force;
                auto                     lengthUnit = string2unit(fs->unit().c_str());
                if (-1 == lengthUnit)
                {
                    gmx_fatal(FARGS, "No such length unit '%s' for bonds", fs->unit().c_str());
                }
                memset(&b, 0, sizeof(b));
                for (auto bi = molProp()->BeginBond(); bi < molProp()->EndBond(); bi++)
                {
                    b.a[0] = bi->getAi() - 1;
                    b.a[1] = bi->getAj() - 1;
                    pd->atypeToBtype(*topology_->atoms.atomtype[b.a[0]], btype1);
                    pd->atypeToBtype(*topology_->atoms.atomtype[b.a[1]], btype2);
                    std::vector<std::string> atoms = {btype1, btype2};
                    if (pd->findForce(atoms, &force))
                    {
                        std::string         pp = force->params();
                        std::vector<double> dd = getDoubles(pp);
                        int                 ii = 0;
                        b.c[ii++] = convert2gmx(force->refValue(), lengthUnit);
                        for (auto &d : dd)
                        {
                            b.c[ii++] = d;
                        }
                        add_param_to_plist(plist_, ftb, eitBONDS, b, bi->getBondOrder());
                    }
                    else
                    {
                        // Insert a dummy bond with a bond order of 1 to be replaced later
                        for (auto i = 0; i < MAXFORCEPARAM; i++)
                        {
                            b.c[i] = 0;
                        }
                        add_param_to_plist(plist_, ftb, eitBONDS, b, 1);
                    }
                }
            }
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

        imm = updatePlist(pd, plist_, topology_, bBASTAT, molProp()->getMolname());
    }
    if (immOK == imm)
    {
        /*Center of charge*/
        auto atntot = 0;
        for (auto i = 0; i < topology_->atoms.nr; i++)
        {
            auto atn = topology_->atoms.atom[i].atomnumber;
            atntot  += atn;
            for (auto m = 0; m < DIM; m++)
            {
                coc_[m] += state_->x[i][m]*atn;
            }
        }
        svmul((1.0/atntot), coc_, coc_);
        /*Center of charge*/

        bool bAddShells = getEemtypePolarizable(iChargeModel);
        if (bAddShells)
        {
            addShells(pd);
        }
        char **molnameptr = put_symtab(symtab_, molProp()->getMolname().c_str());
        mtop_ = do_init_mtop(pd, molnameptr, &topology_->atoms, plist_, inputrec_, symtab_, tabfn); // Generate mtop
        excls_to_blocka(topology_->atoms.nr, excls_, &(mtop_->moltype[0].excls));
        if (bAddShells)
        {
            srenew(mtop_->atomtypes.atomnumber, get_atomtype_ntypes(atype_)); // Update mtop internals to account for shell type
            for (auto i = 0; i < get_atomtype_ntypes(atype_); i++)
            {
                mtop_->atomtypes.atomnumber[i] = get_atomtype_atomnumber(i, atype_);
            }
            mtop_->ffparams.atnr = get_atomtype_ntypes(atype_);
            shellfc_             = init_shell_flexcon(debug, mtop_, 0, 1, false); // Generate shell data structure
        }
        if (nullptr == ltop_)
        {
            ltop_ = gmx_mtop_generate_local_top(mtop_, false); // Generate ltop from mtop
        }
    }
    return imm;
}

void MyMol::addShells(const Poldata *pd)
{
    int                    shell  = 0;
    int                    nshell = 0;
    double                 pol    = 0;
    double                 sigpol = 0;

    std::vector<int>       renum;
    std::vector<int>       inv_renum;

    char                   buf[32];
    char                 **newname;
    t_atoms               *newatoms;
    t_excls               *newexcls;
    std::vector<gmx::RVec> newx;
    t_param                p;

    auto                   polarUnit = string2unit(pd->getPolarUnit().c_str());
    if (-1 == polarUnit)
    {
        gmx_fatal(FARGS, "No such polarizability unit '%s'", pd->getPolarUnit().c_str());
    }

    /*Calculate the total number of Atom and Vsite particles.*/
    auto nParticles = topology_->atoms.nr;
    for (int i = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptVSite)
        {
            nParticles++; // We add 1 shell particle per Atom and Vsite particles
        }
    }
    state_change_natoms(state_, nParticles);
    inv_renum.resize(nParticles, -1);
    renum.resize(topology_->atoms.nr, 0);

    /*Renumber the atoms.*/
    for (int i = 0; i < topology_->atoms.nr; i++)
    {
        renum[i]              = i + nshell;
        inv_renum[i + nshell] = i;
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptVSite)
        {
            nshell++;
        }
    }

    /*Add Polarization to the plist.*/
    memset(&p, 0, sizeof(p));
    for (int i = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptVSite)
        {
            std::string atomtype;
            vsiteType_to_atomType(*topology_->atoms.atomtype[i], &atomtype);
            auto        fa    = pd->findAtype(atomtype);
            if (pd->getAtypeEnd() != fa)
            {
                auto ztype = fa->getZtype();
                if (pd->getAtypePol(atomtype, &pol, &sigpol) && (pol > 0) &&
                    (pd->getNzeta(ztype) == 2))
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

                    if (nexcl_ == 0)
                    {
                        /*
                           This is an ugly hack to turn off the Polarize
                           routine. If we can make nexcl_ = 0 work, it should
                           be implemenetd in a nice way.
                         */
                        pol = 1e+10;
                    }
                    p.c[0] = convert2gmx(pol, polarUnit);
                    add_param_to_plist(plist_, F_POLARIZATION, eitPOLARIZATION, p);
                }
                else
                {
                    gmx_fatal(FARGS, "Polarizability is %f for %s ztype %s btype %s ptype %s.\n",
                              pol, *topology_->atoms.atomtype[i], ztype.c_str(),
                              fa->getBtype().c_str(), fa->getPtype().c_str());
                }
            }
            else
            {
                printf("Cannot find atomtype %s in poldata\n", atomtype.c_str());
            }
        }
    }

    t_atom *shell_atom;
    snew(shell_atom, 1);
    shell_atom->ptype = eptShell;

    /* Make new atoms and x arrays. */
    snew(newatoms, 1);
    init_t_atoms(newatoms, nParticles, true);
    snew(newatoms->atomtype, nParticles);
    snew(newatoms->atomtypeB, nParticles);
    newatoms->nres = topology_->atoms.nres;
    newx.resize(newatoms->nr);
    snew(newname, newatoms->nr);

    /* Make a new exclusion array and put the shells in it. */
    snew(newexcls, newatoms->nr);

    /* Add exclusion for F_POLARIZATION. */
    auto pw = SearchPlist(plist_, F_POLARIZATION);
    if (plist_.end() != pw)
    {
        /*
           Exclude the vsites and the atoms from their own shell.
           This step will be done if the number of exclusions is
           bigger than zero, otherwsie, the vsite or the core will
           interact with its own shell.
         */
        if (nexcl_ > 0)
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
            char buf[256];
            snprintf(buf, sizeof(buf), "Uninitialized inv_renum entry for atom %d (%d) shell %d (%d)",
                     j->a[0], inv_renum[j->a[0]],
                     j->a[1], inv_renum[j->a[1]]);
            GMX_RELEASE_ASSERT(i0 >= 0, buf);
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
    for (int i = 0; i < topology_->atoms.nr; i++)
    {
        newatoms->atom[renum[i]]      = topology_->atoms.atom[i];
        newatoms->atomname[renum[i]]  = put_symtab(symtab_, *topology_->atoms.atomname[i]);
        newatoms->atomtype[renum[i]]  = put_symtab(symtab_, *topology_->atoms.atomtype[i]);
        newatoms->atomtypeB[renum[i]] = put_symtab(symtab_, *topology_->atoms.atomtypeB[i]);
        copy_rvec(state_->x[i], newx[renum[i]]);
        newname[renum[i]] = *topology_->atoms.atomtype[i];
        t_atoms_set_resinfo(newatoms, renum[i], symtab_,
                            *topology_->atoms.resinfo[topology_->atoms.atom[i].resind].name,
                            topology_->atoms.atom[i].resind, ' ', 1, ' ');
    }
    for (int i = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptVSite)
        {
            std::string atomtype;
            auto        iat = renum[i]; // Atom or Vsite
            auto        j   = iat + 1;  // Shell sits next to the Atom or Vsite

            auto        atomtypeName = get_atomtype_name(topology_->atoms.atom[i].type, atype_);
            vsiteType_to_atomType(atomtypeName, &atomtype);

            newatoms->atom[j]               = topology_->atoms.atom[i];
            newatoms->atom[j].m             = 0;
            newatoms->atom[j].mB            = 0;
            newatoms->atom[j].atomnumber    = topology_->atoms.atom[i].atomnumber;
            sprintf(buf, "%s_s", topology_->atoms.atom[i].elem);
            newatoms->atomname[j]           = put_symtab(symtab_, buf);
            sprintf(buf, "%s_s", atomtype.c_str());
            newname[j]                      = strdup(buf);
            shell                           = add_atomtype(atype_, symtab_, shell_atom, buf, &p, 0, 0);
            newatoms->atom[j].type          = shell;
            newatoms->atom[j].typeB         = shell;
            newatoms->atomtype[j]           = put_symtab(symtab_, buf);
            newatoms->atomtypeB[j]          = put_symtab(symtab_, buf);
            newatoms->atom[j].ptype         = eptShell;
            newatoms->atom[j].zetaA         = pd->getZeta(atomtype, 1);
            newatoms->atom[j].zetaB         = newatoms->atom[j].zetaA;
            newatoms->atom[j].row           = pd->getRow(atomtype, 1);
            newatoms->atom[j].resind        = topology_->atoms.atom[i].resind;
            copy_rvec(state_->x[i], newx[j]);

            newatoms->atom[j].q      =
                newatoms->atom[j].qB = pd->getQ(atomtype, 1);
            if (bHaveVSites_)
            {
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
    copy_atoms(newatoms, &topology_->atoms);

    for (int i = 0; i < newatoms->nr; i++)
    {
        copy_rvec(newx[i], state_->x[i]);
        topology_->atoms.atomtype[i] = put_symtab(symtab_, newname[i]);
    }
    sfree(newname);

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
    sfree(shell_atom);
}

immStatus MyMol::GenerateGromacs(const gmx::MDLogger       &mdlog,
                                 t_commrec                 *cr,
                                 const char                *tabfn,
                                 gmx_hw_info_t             *hwinfo,
                                 ChargeModel    ieqd)
{
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
    setup_bonded_threading(fr_->bondedThreading, topology_->atoms.nr, false, ltop_->idef);
    wcycle_    = wallcycle_init(debug, 0, cr);

    MDatoms_  = new std::unique_ptr<gmx::MDAtoms>(new gmx::MDAtoms());
    *MDatoms_ = gmx::makeMDAtoms(nullptr, *mtop_, *inputrec_, false);
    atoms2md(mtop_, inputrec_, -1, nullptr, topology_->atoms.nr, MDatoms_->get());
    auto mdatoms = MDatoms_->get()->mdatoms();
    f_.resizeWithPadding(state_->natoms);

    if (nullptr != shellfc_)
    {
        make_local_shells(cr, mdatoms, shellfc_);
    }
    if (!getEemtypeSlater(ieqd))
    {
        for (auto i = 0; i < mtop_->natoms; i++)
        {
            mdatoms->row[i] = 0;
        }
    }
    return immOK;
}

void MyMol::computeForces(FILE *fplog, t_commrec *cr)
{
    auto mdatoms = MDatoms_->get()->mdatoms();
    if (mdatoms->typeA[0] == 0)
    {
        for (auto i = 0; i < mtop_->natoms; i++)
        {
            mdatoms->chargeA[i] = mtop_->moltype[0].atoms.atom[i].q;
            mdatoms->typeA[i]   = mtop_->moltype[0].atoms.atom[i].type;
            if (nullptr != debug)
            {
                fprintf(debug, "QQQ Setting q[%d] to %g\n", i, mdatoms->chargeA[i]);
            }
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
    if (nullptr != shellfc_)
    {
        auto nnodes = cr->nnodes;
        cr->nnodes  = 1;
        relax_shell_flexcon(fplog, cr, nullptr, false,
                            nullptr, 0, inputrec_,
                            true, force_flags, ltop_, nullptr,
                            enerd_, fcd_, state_,
                            f_.arrayRefWithPadding(), force_vir, mdatoms,
                            &nrnb_, wcycle_, nullptr,
                            &(mtop_->groups), shellfc_,
                            fr_, t, mu_tot, vsite_->get(),
                            DdOpenBalanceRegionBeforeForceComputation::no,
                            DdCloseBalanceRegionAfterForceComputation::no);
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
    }
}

void MyMol::initQgresp(const Poldata             *pd,
                       const char                *lot,
                       real                       watoms,
                       int                        maxESP)
{
    auto iChargeModel = pd->getChargeModel();
                       
    Qgresp_.setChargeModel(iChargeModel);
    Qgresp_.setAtomWeight(watoms);
    Qgresp_.setAtomInfo(&topology_->atoms, pd, state_->x, molProp()->getCharge());
    Qgresp_.setAtomSymmetry(symmetric_charges_);
    Qgresp_.setMolecularCharge(molProp()->getCharge());
    Qgresp_.summary(debug);

    auto ci = molProp()->getLotPropType(lot, MPO_POTENTIAL, nullptr);
    if (ci != molProp()->EndExperiment())
    {
        int mod  = 100/maxESP;
        int iesp = 0;
        for (auto epi = ci->BeginPotential(); epi < ci->EndPotential(); ++epi, ++iesp)
        {
            if (Qgresp_.myWeight(iesp) == 0 || ((iesp-ci->NAtom()) % mod) != 0)
            {
                continue;
            }
            auto xu = string2unit(epi->getXYZunit().c_str());
            auto vu = string2unit(epi->getVunit().c_str());
            if (-1 == xu)
            {
                gmx_fatal(FARGS, "No such length unit '%s' for potential",
                          epi->getXYZunit().c_str());
            }
            if (-1 == vu)
            {
                gmx_fatal(FARGS, "No such potential unit '%s' for potential",
                          epi->getVunit().c_str());
            }
            Qgresp_.addEspPoint(convert2gmx(epi->getX(), xu),
                                convert2gmx(epi->getY(), xu),
                                convert2gmx(epi->getZ(), xu),
                                convert2gmx(epi->getV(), vu));
        }
        if (debug)
        {
            fprintf(debug, "Added %zu ESP points to the RESP structure.\n", Qgresp_.nEsp());
        }
    }
}

immStatus MyMol::GenerateCharges(const Poldata             *pd,
                                 const gmx::MDLogger       &mdlog,
                                 gmx_atomprop_t             ap,
                                 real                       watoms,
                                 real                       hfac,
                                 const char                *lot,
                                 bool                       bSymmetricCharges,
                                 const char                *symm_string,
                                 t_commrec                 *cr,
                                 const char                *tabfn,
                                 gmx_hw_info_t             *hwinfo,
                                 int                        maxiter,
                                 int                        maxESP,
                                 real                       tolerance,
                                 const gmx_output_env_t    *oenv,
                                 const char                *ESPcorr)
{
    std::vector<double> qq;
    immStatus           imm         = immOK;
    bool                converged   = false;
    int                 iter        = 0;
    auto                iChargeModel = pd->getChargeModel();
                                 
    GenerateGromacs(mdlog, cr, tabfn, hwinfo, iChargeModel);
    if (bSymmetricCharges)
    {
        symmetric_charges_.clear();
        ConstPlistWrapperIterator bonds = SearchPlist(plist_, eitBONDS);
        if (plist_.end() != bonds)
        {
            symmetrize_charges(bSymmetricCharges, &topology_->atoms, bonds,
                               pd, ap, symm_string, symmetric_charges_);
        }
    }
    else
    {
        for (auto i = 0; i < topology_->atoms.nr; i++)
        {
            symmetric_charges_.push_back(i);
        }
    }
    switch (chargeGenerationAlgorithm(pd->getChargeModel()))
    {
        case eqgNONE:
            if (debug)
            {
                fprintf(debug, "WARNING! Using zero charges for %s!\n",
                        molProp()->getMolname().c_str());
            }
            for (auto i = 0; i < topology_->atoms.nr; i++)
            {
                topology_->atoms.atom[i].q  = topology_->atoms.atom[i].qB = 0;
            }
            return immOK;
        case eqgESP:
        {
            double chi2[2]   = {1e8, 1e8};
            real   rrms      = 0;
            real   wtot      = 0;
            int    cur       = 0;
            real   cosangle  = 0;
            EspRms_          = 0;
            iter             = 0;

            initQgresp(pd, lot, watoms, maxESP);
            Qgresp_.optimizeCharges();
            Qgresp_.calcPot();
            EspRms_ = chi2[cur] = Qgresp_.getRms(&wtot, &rrms, &cosangle);
            if (debug)
            {
                fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
            }
            do
            {
                for (auto i = 0; i < topology_->atoms.nr; i++)
                {
                    mtop_->moltype[0].atoms.atom[i].q      =
                        mtop_->moltype[0].atoms.atom[i].qB = Qgresp_.getAtomCharge(i);
                }
                if (nullptr != shellfc_)
                {
                    computeForces(nullptr, cr);
                    Qgresp_.updateAtomCoords(state_->x);
                }
                Qgresp_.optimizeCharges();
                Qgresp_.calcPot();
                real cosangle = 0;
                EspRms_ = chi2[cur] = Qgresp_.getRms(&wtot, &rrms, &cosangle);
                if (debug)
                {
                    fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
                }
                converged = (fabs(chi2[cur] - chi2[1-cur]) < tolerance) || (nullptr == shellfc_);
                cur       = 1-cur;
                iter++;
            }
            while ((!converged) && (iter < maxiter));
            for (auto i = 0; i < topology_->atoms.nr; i++)
            {
                topology_->atoms.atom[i].q      =
                    topology_->atoms.atom[i].qB = Qgresp_.getAtomCharge(i);
            }
            if (ESPcorr && oenv)
            {
                Qgresp_.plotLsq(oenv, ESPcorr);
            }
        }
        break;
        case eqgACM:
        {
            Qgacm_.setInfo(pd,
                           &topology_->atoms,
                           hfac,
                           molProp()->getCharge(),
                           bHaveShells_);

            auto q     = Qgacm_.q();
            auto natom = Qgacm_.natom();

            qq.resize(natom + 1);
            for (auto i = 0; i < natom + 1; i++)
            {
                qq[i] = q[i][0];
            }
            iter = 0;
            do
            {
                if (eQGEN_OK == Qgacm_.generateCharges(debug,
                                                       molProp()->getMolname().c_str(),
                                                       pd,
                                                       &topology_->atoms,
                                                       state_->x))
                {
                    for (auto i = 0; i < mtop_->natoms; i++)
                    {
                        mtop_->moltype[0].atoms.atom[i].q      =
                            mtop_->moltype[0].atoms.atom[i].qB = topology_->atoms.atom[i].q;
                    }
                    if (nullptr != shellfc_)
                    {
                        computeForces(nullptr, cr);
                    }
                    q       = Qgacm_.q();
                    EemRms_ = 0;
                    for (auto i = 0; i < natom + 1; i++)
                    {
                        EemRms_  += gmx::square(qq[i] - q[i][0]);
                        qq[i]     = q[i][0];
                    }
                    EemRms_  /= natom;
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
                    mtop_->moltype[0].atoms.atom[i].qB = topology_->atoms.atom[i].q;
            }
            if (!converged)
            {
                printf("Alexandria Charge Model did not converge to %g. rms: %g\n", tolerance, sqrt(EemRms_));
            }
        }
        break;
        default:
            gmx_fatal(FARGS, "Not implemented");
            break;
    }
    return imm;
}

void MyMol::changeCoordinate(ExperimentIterator ei, gmx_bool bpolar)
{
    const std::vector<gmx::RVec> &x = ei->getCoordinates();

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

    for (auto ei = molProp()->BeginExperiment();
         (!bopt) && (ei < molProp()->EndExperiment()); ++ei)
    {
        if (JOB_OPT == ei->getJobtype())
        {
            const std::vector<gmx::RVec> &xxx = ei->getCoordinates();
            for (size_t i = 0; i < xxx.size(); i++)
            {
                copy_rvec(xxx[i], x[i]);
            }
            bopt = true;
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
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        rvec_sub(state_->x[i], coc_, r);
        auto q = e2d(topology_->atoms.atom[i].q);
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
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        rvec_sub(state_->x[i], coc_, r);
        r2   = iprod(r, r);
        q    = topology_->atoms.atom[i].q;
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
    for (i = j = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptNucleus)
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
    int     i, np;
    double  poltot, pol, sigpol, sptot, ereftot, eref;

    poltot  = 0;
    sptot   = 0;
    ereftot = 0;
    np      = 0;
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        if (pd->getAtypePol(*topology_->atoms.atomtype[i], &pol, &sigpol))
        {
            np++;
            poltot += pol;
            sptot  += gmx::square(sigpol);
        }
        if (pd->getAtypeRefEnthalpy(*topology_->atoms.atomtype[i], &eref))
        {
            ereftot += eref;
        }
    }
    CalcDipole(mu);
    CalcQuadrupole();
    ref_enthalpy_   = ereftot;
    polarizability_ = poltot;
    sig_pol_        = sqrt(sptot/topology_->atoms.nr);
}

void MyMol::CalcAnisoPolarizability(tensor polar, double *anisoPol)
{
    auto a = gmx::square(polar[XX][XX] - polar[YY][YY]);
    auto b = gmx::square(polar[XX][XX] - polar[ZZ][ZZ]);
    auto c = gmx::square(polar[ZZ][ZZ] - polar[YY][YY]);
    auto d = 6 * (gmx::square(polar[XX][YY]) + gmx::square(polar[XX][ZZ]) + gmx::square(polar[ZZ][YY]));

    *anisoPol = sqrt(1/2.0) * sqrt(a + b + c + d);
}

void MyMol::CalcPolarizability(double     efield,
                               t_commrec *cr,
                               FILE      *fplog)
{
    const double        POLFAC = 29.957004; /* C.m**2.V*-1 to Å**3 */
    std::vector<double> field;
    rvec                mu_ref, mu_tot;

    field.resize(DIM, 0);
    myforce_->setField(field);
    CalcDipole(mu_ref);
    computeForces(fplog, cr);
    for (auto m = 0; m < DIM; m++)
    {
        field[m] = efield;
        myforce_->setField(field);
        computeForces(fplog, cr);
        CalcDipole(mu_tot);
        for (auto n = 0; n < DIM; n++)
        {
            alpha_calc_[n][m] = ((mu_tot[n]-mu_ref[n])/efield)*(POLFAC);
        }
        isoPol_calc_     += alpha_calc_[m][m];
        field[m]          = 0.0;
    }
    isoPol_calc_ /= DIM;
    CalcAnisoPolarizability(alpha_calc_, &anisoPol_calc_);
}

void MyMol::PrintConformation(const char *fn)
{
    char title[STRLEN];

    put_in_box(topology_->atoms.nr, state_->box, as_rvec_array(state_->x.data()), 0.3);
    sprintf(title, "%s processed by alexandria", molProp()->getMolname().c_str());
    write_sto_conf(fn, title, &topology_->atoms, as_rvec_array(state_->x.data()), nullptr, epbcNONE, state_->box);
}

void MyMol::PrintTopology(const char     *fn,
                          bool            bVerbose,
                          const Poldata  *pd,
                          gmx_atomprop_t  aps,
                          t_commrec      *cr,
                          double          efield,
                          const char     *lot)
{
    FILE  *fp   = gmx_ffopen(fn, "w");
    bool   bITP = (fn2ftp(fn) == efITP);

    PrintTopology(fp, bVerbose, pd, aps, bITP, cr, efield, lot);

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
                          const char             *lot)
{
    char                     buf[256];
    t_mols                   printmol;
    std::vector<std::string> commercials;
    rvec                     vec, mu;
    tensor                   myQ;
    double                   value = 0, error = 0, T = -1;
    std::string              myref, mylot;
    auto                     iChargeModel = pd->getChargeModel();
                          
    if (fp == nullptr)
    {
        return;
    }

    CalcQPol(pd, mu);
    if (molProp()->getMolname().size() > 0)
    {
        printmol.name = strdup(molProp()->getMolname().c_str());
    }
    else if (molProp()->formula().size() > 0)
    {
        printmol.name = strdup(molProp()->formula().c_str());
    }
    else
    {
        printmol.name = strdup("Unknown");
    }

    printmol.nr = 1;

    snprintf(buf, sizeof(buf), "Total Mass = %.3f (Da)", molProp()->getMass());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Reference_Enthalpy = %.3f (kJ/mol)", ref_enthalpy_);
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Total Charge = %d (e)", molProp()->getCharge());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Charge Type  = %s\n", getEemtypeName(iChargeModel));
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Alexandria Dipole Moment (Debye):\n"
             "; ( %.2f %6.2f %6.2f ) Total= %.2f\n",
             mu[XX], mu[YY], mu[ZZ],
             norm(mu));
    commercials.push_back(buf);

    T = -1;
    const char *qm_type = "electronic";
    const char *qm_conf = "minimum";
    if (molProp()->getPropRef(MPO_DIPOLE, iqmQM, lot, qm_conf,
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
                 lot,
                 mu_qm_[qtElec][XX], mu_qm_[qtElec][YY], mu_qm_[qtElec][ZZ],
                 norm(mu_qm_[qtElec]));
        commercials.push_back(buf);
    }
    else
    {
        printf("WARNING: QM dipole of type %s not found for lot %s\n", qm_type, lot);
    }

    add_tensor(&commercials,
               "Alexandria Traceless Quadrupole Moments (Buckingham)",
               QQM(qtCalc));

    T = -1;
    if (molProp()->getPropRef(MPO_QUADRUPOLE, iqmQM, lot, qm_conf,
                              qm_type, &value, &error,
                              &T, &myref, &mylot, vec, myQ))
    {
        set_QQM(qtElec, myQ);
        rotate_tensor(Q_qm_[qtElec], Q_qm_[qtCalc]);
        snprintf(buf, sizeof(buf), "%s Traceless Quadrupole Moments (Buckingham)", lot);
        add_tensor(&commercials, buf, Q_qm_[qtElec]);
    }
    else
    {
        printf("WARNING: QM quadrupole of type %s not found for lot %s\n", qm_type, lot);
    }


    snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability (Additive Law): %.2f +/- %.2f (A^3)\n", polarizability_, sig_pol_);
    commercials.push_back(buf);

    if (efield > 0 && nullptr != cr)
    {
        CalcPolarizability(efield, cr, debug);
        add_tensor(&commercials, "Alexandria Polarizability components (A^3)", alpha_calc_);

        snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability (Interactive): %.2f (A^3)\n", isoPol_calc_);
        commercials.push_back(buf);

        snprintf(buf, sizeof(buf), "Alexandria Anisotropic Polarizability: %.2f (A^3)\n", anisoPol_calc_);
        commercials.push_back(buf);

        T = -1;
        if (molProp()->getPropRef(MPO_POLARIZABILITY, iqmQM, lot, "",
                                  (char *)"electronic", &isoPol_elec_, &error,
                                  &T, &myref, &mylot, vec, alpha_elec_))
        {
            CalcAnisoPolarizability(alpha_elec_, &anisoPol_elec_);
            snprintf(buf, sizeof(buf), "%s + Polarizability components (A^3)", lot);
            add_tensor(&commercials, buf, alpha_elec_);

            snprintf(buf, sizeof(buf), "%s Isotropic Polarizability: %.2f (A^3)\n", lot, isoPol_elec_);
            commercials.push_back(buf);
            snprintf(buf, sizeof(buf), "%s Anisotropic Polarizability: %.2f (A^3)\n", lot, anisoPol_elec_);
            commercials.push_back(buf);
        }
    }

    print_top_header(fp, pd, aps, bHaveShells_, commercials, bITP);
    write_top(fp, printmol.name, &topology_->atoms, false,
              plist_, excls_, atype_, cgnr_, nexcl_, pd);
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
    if ((cgnr_ = generate_charge_groups(ecg, &topology_->atoms,
                                        plist_,
                                        bUsePDBcharge,
                                        &qtot, &mtot)) == nullptr)
    {
        return immChargeGeneration;
    }
    if (ecg != ecgAtom)
    {
        //sort_on_charge_groups(cgnr_, &topology_->atoms,
        //                    plist_, x_, excls_, ndxfn, nmol);
    }
    return immOK;
}

void MyMol::GenerateCube(const Poldata          *pd,
                         real                    spacing,
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
    ChargeModel iChargeModel = pd->getChargeModel();
                         
    if ((nullptr  != potfn) || (nullptr != hisfn) || (nullptr != rhofn) ||
        ((nullptr != difffn) && (nullptr != reffn)))
    {
        char      buf[256];
        char     *gentop_version = (char *)"v0.99b";
        QgenResp  grref;

        Qgresp_.potcomp(pcfn, pdbdifffn, oenv);

        /* This has to be done before the grid is f*cked up by
           writing a cube file */
        grref = Qgresp_;

        sprintf(buf, "Potential generated by %s based on %s charges",
                gentop_version,
                getEemtypeName(iChargeModel));

        if (nullptr != difffn)
        {
            grref.setAtomInfo(&topology_->atoms, pd, state_->x, molProp()->getCharge());
            grref.setAtomSymmetry(symmetric_charges_);
            grref.readCube(reffn, FALSE);
            Qgresp_ = grref;
        }
        else
        {
            Qgresp_.makeGrid(spacing, state_->box, as_rvec_array(state_->x.data()));
        }
        if (nullptr != rhofn)
        {
            sprintf(buf, "Electron density generated by %s based on %s charges",
                    gentop_version, getEemtypeName(iChargeModel));
            Qgresp_.calcRho();
            Qgresp_.writeRho(rhofn, buf);
        }
        sprintf(buf, "Potential generated by %s based on %s charges",
                gentop_version, getEemtypeName(iChargeModel));
        if (nullptr != potfn)
        {
            Qgresp_.calcPot();
            Qgresp_.writeCube(potfn, buf);
        }
        if (nullptr != hisfn)
        {
            Qgresp_.writeHisto(hisfn, buf, oenv);
        }
        if ((nullptr != difffn) || (nullptr != diffhistfn))
        {
            sprintf(buf, "Potential difference generated by %s based on %s charges",
                    gentop_version,
                    getEemtypeName(iChargeModel));

            Qgresp_.writeDiffCube(grref, difffn, diffhistfn, buf, oenv, 0);
        }
    }
}

void MyMol::rotateDipole(rvec mu, rvec muReference)
{
    if (norm(mu) == 0 or norm(muReference) == 0)
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
        for (i = j = 0; i < topology_->atoms.nr; i++)
        {
            if (topology_->atoms.atom[i].ptype == eptAtom ||
                topology_->atoms.atom[i].ptype == eptNucleus)
            {
                charge_QM_[qt][j] = q[j];
                j++;
            }
        }
        CalcQMbasedMoments(q, mu_qm_[qt], Q_qm_[qt]);
    }
}

immStatus MyMol::getExpProps(gmx_bool bQM, gmx_bool bZero,
                             gmx_bool bZPE, const char *lot,
                             const Poldata *pd)
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


    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptNucleus)
        {
            natom++;
        }
    }
    real q[natom];
    if (molProp()->getPropRef(MPO_CHARGE, iqmQM,
                              (char *)mylot.c_str(), "",
                              (char *)"ESP charges",
                              &value, &error, &T,
                              &myref, &mylot, q, quadrupole))
    {
        setQandMoments(qtESP, natom, q);
        esp_dipole_found = true;
    }
    T = -1;
    if (molProp()->getPropRef(MPO_CHARGE, iqmQM,
                              (char *)mylot.c_str(), "",
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
    if (molProp()->getPropRef(MPO_CHARGE, iqmQM,
                              (char *)mylot.c_str(), "",
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
    if (molProp()->getPropRef(MPO_CHARGE, iqmQM,
                              (char *)mylot.c_str(), "",
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
    if (molProp()->getProp(MPO_ENERGY, (bQM ? iqmQM : iqmExp),
                           lot, "", (char *)"DeltaHform", &value, &error, &T))
    {
        Hform_ = value;
        Emol_  = value;
        for (ia = 0; ia < topology_->atoms.nr; ia++)
        {
            if (topology_->atoms.atom[ia].ptype == eptAtom ||
                topology_->atoms.atom[ia].ptype == eptNucleus)
            {
                if (pd->getAtypeRefEnthalpy(*topology_->atoms.atomtype[ia], &Hatom))
                {
                    Emol_ -= Hatom;
                }
                else
                {
                    fprintf(stderr, "WARNING: NO reference enthalpy for molecule %s.\n",
                            molProp()->getMolname().c_str());
                    Emol_ = 0;
                    imm   = immNoData;
                    break;
                }
            }
        }
        if (bZPE)
        {

            if (molProp()->getProp(MPO_ENERGY, iqmBoth, lot, "",
                                   (char *)"ZPE", &ZPE, &error, &T))
            {
                Emol_ -= ZPE;
            }
            else
            {
                fprintf(stderr, "No zero-point energy for molecule %s.\n",
                        molProp()->getMolname().c_str());
                imm = immNoData;
            }
        }
        if (ia < topology_->atoms.nr)
        {
            imm = immNoData;
        }
    }
    if (imm == immOK)
    {
        T = -1;
        if (molProp()->getPropRef(MPO_DIPOLE, iqmQM, lot, "",
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
                            molProp()->getMolname().c_str(), error);
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
        if (molProp()->getPropRef(MPO_QUADRUPOLE, iqmQM,
                                  lot, "", (char *)"electronic",
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
        if (molProp()->getPropRef(MPO_POLARIZABILITY, iqmQM,
                                  lot, "", (char *)"electronic",
                                  &isoPol_elec_, &error, &T,
                                  &myref, &mylot, vec, polar))
        {
            copy_mat(polar, alpha_elec_);
            CalcAnisoPolarizability(alpha_elec_, &anisoPol_elec_);
        }
    }
    return imm;
}

void MyMol::UpdateIdef(const Poldata   *pd,
                       InteractionType  iType)
{
    std::string params;

    if (debug)
    {
        fprintf(debug, "UpdateIdef for %s\n", iType2string(iType));
    }
    if (iType == eitVDW)
    {
        nonbondedFromPdToMtop(mtop_, &topology_->atoms, pd, fr_);
        if (debug)
        {
            pr_ffparams(debug, 0, "UpdateIdef Before", &mtop_->ffparams, FALSE);
        }
        return;
    }
    auto fs = pd->findForces(iType);
    if (fs == pd->forcesEnd())
    {
        gmx_fatal(FARGS, "Can not find the force %s to update",
                  iType2string(iType));
    }
    auto ftype = fs->fType();
    // Make a list of bonded types that can be indexed
    // with the atomtype. That should speed up the code
    // below somewhat.
    std::vector<std::string> btype(mtop_->ffparams.atnr);
    for (int i = 0; i < mtop_->natoms; i++)
    {
        std::string bt;
        if (pd->atypeToBtype(*topology_->atoms.atomtype[i], bt))
        {
            btype[topology_->atoms.atom[i].type] = bt;
        }
        else if (!(topology_->atoms.atom[i].ptype == eptShell ||
                   topology_->atoms.atom[i].ptype == eptVSite))
        {
            gmx_fatal(FARGS, "Cannot find bonded type for atomtype %s",
                      *topology_->atoms.atomtype[i]);
        }
    }
    switch (iType)
    {
        case eitBONDS:
        {
            auto lu = string2unit(fs->unit().c_str());
            if (-1 == lu)
            {
                gmx_fatal(FARGS, "Unknown length unit '%s' for bonds",
                          fs->unit().c_str());
            }
            for (auto i = 0; i < ltop_->idef.il[ftype].nr; i += interaction_function[ftype].nratoms+1)
            {
                auto                     tp  = ltop_->idef.il[ftype].iatoms[i];
                std::vector<std::string> atoms;
                for (int j = 1; j < interaction_function[ftype].nratoms+1; j++)
                {
                    atoms.push_back(btype[topology_->atoms.atom[ltop_->idef.il[ftype].iatoms[i+j]].type]);
                }
                double value, sigma;
                size_t ntrain;
                if (pd->searchForce(atoms, params, &value, &sigma, &ntrain, iType))
                {
                    auto bondLength = convert2gmx(value, lu);
                    auto parameters = gmx::splitString(params);
                    switch (ftype)
                    {
                        case F_MORSE:
                        {
                            mtop_->ffparams.iparams[tp].morse.b0A         =
                                mtop_->ffparams.iparams[tp].morse.b0B     =
                                    ltop_->idef.iparams[tp].morse.b0A     =
                                        ltop_->idef.iparams[tp].morse.b0B = bondLength;
                            GMX_RELEASE_ASSERT(parameters.size() == 2, "Need exactly two parameters for Morse bonds");
                            mtop_->ffparams.iparams[tp].morse.cbA         =
                                mtop_->ffparams.iparams[tp].morse.cbB     =
                                    ltop_->idef.iparams[tp].morse.cbA     =
                                        ltop_->idef.iparams[tp].morse.cbB =
                                            gmx::doubleFromString(parameters[0].c_str());
                            mtop_->ffparams.iparams[tp].morse.betaA         =
                                mtop_->ffparams.iparams[tp].morse.betaB     =
                                    ltop_->idef.iparams[tp].morse.betaA     =
                                        ltop_->idef.iparams[tp].morse.betaB =
                                            gmx::doubleFromString(parameters[1].c_str());
                        }
                        break;
                        default:
                            gmx_fatal(FARGS, "Unsupported bondtype %s in UpdateIdef",
                                      fs->function().c_str());
                    }
                }
            }
        }
        break;
        case eitANGLES:
        case eitLINEAR_ANGLES:
        {
            for (auto i = 0; i < ltop_->idef.il[ftype].nr; i += interaction_function[ftype].nratoms+1)
            {
                auto                     tp  = ltop_->idef.il[ftype].iatoms[i];
                std::vector<std::string> atoms;
                for (int j = 1; j < interaction_function[ftype].nratoms+1; j++)
                {
                    atoms.push_back(btype[topology_->atoms.atom[ltop_->idef.il[ftype].iatoms[i+j]].type]);
                }
                double angle, sigma;
                size_t ntrain;
                if (pd->searchForce(atoms, params, &angle, &sigma, &ntrain, iType))
                {
                    auto parameters = gmx::splitString(params);
                    auto r13        = calc_r13(pd, atoms[0], atoms[1], atoms[2], angle);
                    switch (ftype)
                    {
                        case F_ANGLES:
                        {
                            mtop_->ffparams.iparams[tp].harmonic.rA         =
                                mtop_->ffparams.iparams[tp].harmonic.rB     =
                                    ltop_->idef.iparams[tp].harmonic.rA     =
                                        ltop_->idef.iparams[tp].harmonic.rB = angle;
                            GMX_RELEASE_ASSERT(parameters.size() == 1, "Need exactly one parameters for Harmonic angles");
                            mtop_->ffparams.iparams[tp].harmonic.krA         =
                                mtop_->ffparams.iparams[tp].harmonic.krB     =
                                    ltop_->idef.iparams[tp].harmonic.krA     =
                                        ltop_->idef.iparams[tp].harmonic.krB =
                                            gmx::doubleFromString(parameters[0].c_str());
                        }
                        break;
                        case F_UREY_BRADLEY:
                        {
                            mtop_->ffparams.iparams[tp].u_b.thetaA         =
                                mtop_->ffparams.iparams[tp].u_b.thetaB     =
                                    ltop_->idef.iparams[tp].u_b.thetaA     =
                                        ltop_->idef.iparams[tp].u_b.thetaB = angle;

                            mtop_->ffparams.iparams[tp].u_b.r13A         =
                                mtop_->ffparams.iparams[tp].u_b.r13B     =
                                    ltop_->idef.iparams[tp].u_b.r13A     =
                                        ltop_->idef.iparams[tp].u_b.r13B = r13;
                            GMX_RELEASE_ASSERT(parameters.size() == 2, "Need exactly two parameters for Urey-Bradley angles");

                            mtop_->ffparams.iparams[tp].u_b.kthetaA         =
                                mtop_->ffparams.iparams[tp].u_b.kthetaB     =
                                    ltop_->idef.iparams[tp].u_b.kthetaA     =
                                        ltop_->idef.iparams[tp].u_b.kthetaB =
                                            gmx::doubleFromString(parameters[0].c_str());
                            mtop_->ffparams.iparams[tp].u_b.kUBA         =
                                mtop_->ffparams.iparams[tp].u_b.kUBB     =
                                    ltop_->idef.iparams[tp].u_b.kUBA     =
                                        ltop_->idef.iparams[tp].u_b.kUBB =
                                            gmx::doubleFromString(parameters[1].c_str());
                        }
                        break;
                        case F_LINEAR_ANGLES:
                        {
                            double relative_position = calc_relposition(pd, atoms[0], atoms[1], atoms[2]);

                            mtop_->ffparams.iparams[tp].linangle.aA         =
                                mtop_->ffparams.iparams[tp].linangle.aB     =
                                    ltop_->idef.iparams[tp].linangle.aA     =
                                        ltop_->idef.iparams[tp].linangle.aB = relative_position;

                            mtop_->ffparams.iparams[tp].linangle.r13A         =
                                mtop_->ffparams.iparams[tp].linangle.r13B     =
                                    ltop_->idef.iparams[tp].linangle.r13A     =
                                        ltop_->idef.iparams[tp].linangle.r13B = r13;
                            GMX_RELEASE_ASSERT(parameters.size() == 2, "Need exactly two parameters for Linear angles");

                            mtop_->ffparams.iparams[tp].linangle.klinA         =
                                mtop_->ffparams.iparams[tp].linangle.klinB     =
                                    ltop_->idef.iparams[tp].linangle.klinA     =
                                        ltop_->idef.iparams[tp].linangle.klinB =
                                            gmx::doubleFromString(parameters[0].c_str());
                            mtop_->ffparams.iparams[tp].linangle.kUBA         =
                                mtop_->ffparams.iparams[tp].linangle.kUBB     =
                                    ltop_->idef.iparams[tp].linangle.kUBA     =
                                        ltop_->idef.iparams[tp].linangle.kUBB =
                                            gmx::doubleFromString(parameters[1].c_str());
                        }
                        break;
                        default:
                            gmx_fatal(FARGS, "Unsupported angletype %s in UpdateIdef",
                                      fs->function().c_str());
                    }
                }
            }
        }
        break;
        case eitPROPER_DIHEDRALS:
        case eitIMPROPER_DIHEDRALS:
        {
            for (auto i = 0; i < ltop_->idef.il[ftype].nr; i += interaction_function[ftype].nratoms+1)
            {
                auto                     tp  = ltop_->idef.il[ftype].iatoms[i];
                std::vector<std::string> atoms;
                for (int j = 1; j < interaction_function[ftype].nratoms+1; j++)
                {
                    atoms.push_back(btype[topology_->atoms.atom[ltop_->idef.il[ftype].iatoms[i+j]].type]);
                }
                double angle, sigma;
                size_t ntrain;
                if (pd->searchForce(atoms, params, &angle, &sigma, &ntrain, iType))
                {
                    auto parameters = gmx::splitString(params);
                    switch (ftype)
                    {
                        case F_FOURDIHS:
                        {
                            std::vector<double> old;
                            for (auto parm : parameters)
                            {
                                old.push_back(gmx::doubleFromString(parm.c_str())); 
                            }
                            auto newparam = &mtop_->ffparams.iparams[tp];
                            newparam->rbdihs.rbcA[0] = old[1]+0.5*(old[0]+old[2]);
                            newparam->rbdihs.rbcA[1] = 0.5*(3.0*old[2]-old[0]);
                            newparam->rbdihs.rbcA[2] = 4.0*old[3]-old[1];
                            newparam->rbdihs.rbcA[3] = -2.0*old[2];
                            newparam->rbdihs.rbcA[4] = -4.0*old[3];
                            newparam->rbdihs.rbcA[5] = 0.0;
                            for(int k = 0; k < NR_RBDIHS; k++)
                            {
                                ltop_->idef.iparams[tp].rbdihs.rbcA[k] = 
                                    ltop_->idef.iparams[tp].rbdihs.rbcB[k] = 
                                    newparam->rbdihs.rbcB[k] = 
                                    newparam->rbdihs.rbcA[k];
                            }
                            break;
                        }
                        case F_PDIHS:
                        {
                            mtop_->ffparams.iparams[tp].pdihs.phiA         =
                                mtop_->ffparams.iparams[tp].pdihs.phiB     =
                                    ltop_->idef.iparams[tp].pdihs.phiA     =
                                        ltop_->idef.iparams[tp].pdihs.phiB = angle;

                            GMX_RELEASE_ASSERT(parameters.size() == 2, "Need exactly two parameters for proper dihedrals");
                            mtop_->ffparams.iparams[tp].pdihs.cpA         =
                                mtop_->ffparams.iparams[tp].pdihs.cpB     =
                                    ltop_->idef.iparams[tp].pdihs.cpA     =
                                        ltop_->idef.iparams[tp].pdihs.cpB =
                                            gmx::doubleFromString(parameters[0].c_str());
                            /* Multiplicity for Proper Dihedral must be integer
                               This assumes that the second parameter is Multiplicity */
                            mtop_->ffparams.iparams[tp].pdihs.mult =
                                ltop_->idef.iparams[tp].pdihs.mult =
                                    atoi(parameters[1].c_str());
                        }
                        break;
                        case F_IDIHS:
                        {
                            mtop_->ffparams.iparams[tp].harmonic.rA         =
                                mtop_->ffparams.iparams[tp].harmonic.rB     =
                                    ltop_->idef.iparams[tp].harmonic.rA     =
                                        ltop_->idef.iparams[tp].harmonic.rB = angle;

                            GMX_RELEASE_ASSERT(parameters.size() == 1, "Need exactly one parameter for proper dihedrals");
                            mtop_->ffparams.iparams[tp].harmonic.krA         =
                                mtop_->ffparams.iparams[tp].harmonic.krB     =
                                    ltop_->idef.iparams[tp].harmonic.krA     =
                                        ltop_->idef.iparams[tp].harmonic.krB =
                                            gmx::doubleFromString(parameters[0].c_str());
                        }
                        break;
                        default:
                            gmx_fatal(FARGS, "Unsupported dihedral type %s in UpdateIdef",
                                      fs->function().c_str());
                    }
                }
            }
        }
        break;
        case eitLJ14:
        case eitPOLARIZATION:
        {
            auto pw = SearchPlist(plist_, iType);
            if (plist_.end() != pw)
            {
                auto ft = pw->getFtype();
                for (auto i = 0; i < ltop_->idef.il[ft].nr; i += interaction_function[ft].nratoms+1)
                {
                    auto   tp  = ltop_->idef.il[ft].iatoms[i];
                    auto   ai  = ltop_->idef.il[ft].iatoms[i+1];
                    double alpha, sigma;
                    if (pd->getAtypePol(*topology_->atoms.atomtype[ai], &alpha, &sigma))
                    {
                        mtop_->ffparams.iparams[tp].polarize.alpha =
                            ltop_->idef.iparams[tp].polarize.alpha = alpha;
                    }
                }
            }
        }
        break;
        case eitVDW:
        case eitVSITE2:
        case eitVSITE3FAD:
        case eitVSITE3OUT:
        case eitCONSTR:
        case eitNR:
        default:
            break;
    }
}
} // namespace alexandria
