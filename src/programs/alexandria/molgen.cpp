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


#include "molgen.h"

#include <cmath>

#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/smalloc.h"

#include "alex_modules.h"
#include "fill_inputrec.h"
#include "gmx_simple_comm.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata_xml.h"

#define STRLEN 256

const char *rmsName(int e)
{
    static const char *rms[ermsNR] =
    {
        "BOUNDS", "MU", "QUAD", "CHARGE", "ESP",
        "EPOT", "Force2", "Polar", "Penalty", "TOT"
    };
    if (e >= 0 && e < ermsNR)
    {
        return rms[e];
    }
    else
    {
        return "Incorrect index in rmsName";
    }
}

namespace alexandria
{

static void dump_index_count(const IndexCount       *ic,
                             FILE                   *fp,
                             const Poldata          &pd,
                             gmx_bool                bFitZeta)
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "Atom index for this optimization.\n");
    fprintf(fp, "Name  Number  Action   #Zeta\n");
    for (auto i = ic->beginIndex(); i < ic->endIndex(); ++i)
    {
        int nZeta  = pd.getNzeta(i->name());
        int nZopt  = 0;
        for (int j = 0; (j < nZeta); j++)
        {
            if (pd.getZeta(i->name(), j) > 0)
            {
                nZopt++;
            }
        }
        if (i->isConst())
        {
            fprintf(fp, "%-4s  %6d  Constant\n",
                    i->name().c_str(), i->count());
        }
        else
        {
            fprintf(fp, "%-4s  %6d  Optimized %4d%s\n",
                    i->name().c_str(),
                    i->count(), nZopt,
                    bFitZeta ? " optimized" : " constant");
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
}

static void make_index_count(IndexCount                *ic,
                             const Poldata             &pd,
                             char                      *opt_elem,
                             gmx_bool                   bFitZeta)
{
    if (nullptr != opt_elem)
    {
        /* Only members of opt_elem will be optimized */
        std::vector<std::string> opt_elems = gmx::splitString(opt_elem);
        for (auto eep = pd.BeginEemprops(); eep != pd.EndEemprops(); ++eep)
        {
            auto bConst   = true;
            auto name     = eep->getName();
            auto opt_elem = std::find_if(opt_elems.begin(), opt_elems.end(),
                                         [name](const std::string a)
                                         {
                                             return a == name;
                                         });

            if (opt_elems.end() != opt_elem)
            {
                bConst = false;
            }

            ic->addName(eep->getName(), bConst);
        }
    }
    else
    {
        /*All types of Eeemprops will be optimized */
        bool bConst = false;
        for (auto eep = pd.BeginEemprops(); eep != pd.EndEemprops(); ++eep)
        {
            ic->addName(eep->getName(), bConst);
        }
    }
    dump_index_count(ic, debug, pd, bFitZeta);
}

void IndexCount::addName(const std::string &name,
                         bool               bConst)
{
    auto ai = std::find_if(atomIndex_.begin(), atomIndex_.end(),
                           [name](const AtomIndex a)
                           {
                               return a.name().compare(name) == 0;
                           });
    if (atomIndex_.end() == ai)
    {
        AtomIndex aaa(name, bConst);
        atomIndex_.push_back(aaa);
    }
    else
    {
        if (ai->isConst() == bConst)
        {
            gmx_fatal(FARGS, "Trying to add atom %s as both constant and optimized", name.c_str());
        }
        else
        {
            fprintf(stderr, "Trying to add %s twice\n", name.c_str());
            // ai.increment();
        }
    }
}

void IndexCount::sumCount(t_commrec *cr)
{
    totCount_.resize(atomIndex_.size(), 0);
    int i = 0;
    for (const auto &ai : atomIndex_)
    {
        totCount_[i++] = ai.count();
    }
    if (cr->nnodes > 1)
    {
        gmx_sumi(totCount_.size(), totCount_.data(), cr);
    }
}

void IndexCount::incrementName(const std::string &name)
{
    auto ai = findName(name);
    if (ai == atomIndex_.end())
    {
        gmx_fatal(FARGS, "No such atom %s", name.c_str());
    }
    ai->increment();
}

bool IndexCount::isOptimized(const std::string &name)
{
    bool isObtimized = false;
    auto ai          = findName(name);
    if (ai != atomIndex_.end())
    {
        if (!ai->isConst())
        {
            isObtimized = true;
        }
    }
    return isObtimized;
}

void IndexCount::decrementName(const std::string &name)
{
    auto ai = findName(name);
    if (ai == atomIndex_.end())
    {
        gmx_fatal(FARGS, "No such atom %s", name.c_str());
    }
    ai->decrement();
}

int IndexCount::count(const std::string &name)
{
    auto ai = findName(name);
    if (ai == atomIndex_.end())
    {
        return ai->count();
    }
    return 0;
}

int IndexCount::cleanIndex(int   minimum_data,
                           FILE *fp)
{
    int nremove = 0;

    for (auto ai = atomIndex_.begin(); ai < atomIndex_.end(); )
    {
        if (!ai->isConst() && (ai->count() < minimum_data))
        {
            if (fp)
            {
                fprintf(fp, "Not enough support in data set for optimizing %s\n",
                        ai->name().c_str());
            }
            ai = atomIndex_.erase(ai);
            nremove++;
        }
        else
        {
            ++ai;
        }
    }
    return nremove;
}

MolGen::MolGen()
{
    cr_        = nullptr;
    bFinal_    = false;
    bDone_     = false;
    bGenVsite_ = false;
    bOptHfac_  = false;
    qsymm_     = false;
    constrain_ = false;
    bQM_       = false;
    J0_min_    = 5;
    Chi0_min_  = 1;
    zeta_min_  = 2;
    J0_max_    = 30;
    Chi0_max_  = 30;
    zeta_max_  = 30;
    watoms_    = 0;
    qtol_      = 1e-6;
    qcycle_    = 500;
    mindata_   = 3;
    nexcl_     = 2;
    hfac_      = 0;
    maxESP_    = 100;
    fixchi_    = (char *)"";
    lot_       = "B3LYP/aug-cc-pVTZ";
    inputrec_  = new t_inputrec();
    fill_inputrec(inputrec_);
    fc_[ermsTOT] = 1;
}

MolGen::~MolGen()
{
    if (cr_)
    {
        done_commrec(cr_);
    }
}

static void doAddOptions(std::vector<t_pargs> *pargs, size_t npa, t_pargs pa[])
{
    for (size_t i = 0; i < npa; i++)
    {
        pargs->push_back(pa[i]);
    }
}

void MolGen::addOptions(std::vector<t_pargs> *pargs, eTune etune)
{
    t_pargs pa_general[] =
    {
        { "-mindata", FALSE, etINT, {&mindata_},
          "Minimum number of data points to optimize force field parameters" },
        { "-lot",    FALSE, etSTR,  {&lot_},
          "Use this method and level of theory when selecting coordinates and charges. Multiple levels can be specified which will be used in the order given, e.g.  B3LYP/aug-cc-pVTZ:HF/6-311G**" },
        { "-fc_bound",    FALSE, etREAL, {&fc_[ermsBOUNDS]},
          "Force constant in the penalty function for going outside the borders given with the fitting options (see below)." },
        { "-qtol",   FALSE, etREAL, {&qtol_},
          "Tolerance for assigning charge generation algorithm." },
        { "-qcycle", FALSE, etINT, {&qcycle_},
          "Max number of tries for optimizing the charges." },
        { "-qsymm",  FALSE, etBOOL, {&qsymm_},
          "Symmetrize the charges on symmetric groups, e.g. CH3, NH2." },
        { "-constrain",  FALSE, etBOOL, {&constrain_},
          "Perform Box-Constraint optimization" },
        { "-genvsites", FALSE, etBOOL, {&bGenVsite_},
          "Generate virtual sites. Check and double check." },
        { "-nexcl",  FALSE, etINT, {&nexcl_},
          "[HIDDEN]Exclusion number." },
        { "-qm",     FALSE, etBOOL, {&bQM_},
          "[HIDDEN]Use only quantum chemistry results (from the levels of theory below) in order to fit the parameters. If not set, experimental values will be used as reference with optional quantum chemistry results, in case no experimental results are available" }
    };
    t_pargs pa_zeta[] =
    {
        { "-watoms", FALSE, etREAL, {&watoms_},
          "Weight for the atoms when fitting the charges to the electrostatic potential. The potential on atoms is usually two orders of magnitude larger than on other points (and negative). For point charges or single smeared charges use zero. For point+smeared charges 1 is recommended." },
        { "-maxpot", FALSE, etINT, {&maxESP_},
          "Maximum percent of the electrostatic potential points that will be used to fit partial charges." },
        { "-fixchi", FALSE, etSTR,  {&fixchi_},
          "Electronegativity for this atom type is fixed. Set to FALSE if you want this variable as well, but read the help text above (or somewhere)." },
        { "-z0",    FALSE, etREAL, {&zeta_min_},
          "Minimum value that inverse radius (1/nm) can obtain in fitting" },
        { "-z1",    FALSE, etREAL, {&zeta_max_},
          "Maximum value that inverse radius (1/nm) can obtain in fitting" }
    };
    t_pargs pa_eem[] =
    {
        { "-j0",    FALSE, etREAL, {&J0_min_},
          "Minimum value that J0 (eV) can obtain in fitting" },
        { "-j1",    FALSE, etREAL, {&J0_max_},
          "Maximum value that J0 (eV) can obtain in fitting" },
        { "-chi0",    FALSE, etREAL, {&Chi0_min_},
          "Minimum value that Chi0 (eV) can obtain in fitting" },
        { "-chi1",    FALSE, etREAL, {&Chi0_max_},
          "Maximum value that Chi0 (eV) can obtain in fitting" },
        { "-fc_mu",    FALSE, etREAL, {&fc_[ermsMU]},
          "Force constant in the penalty function for the magnitude of the dipole components." },
        { "-fc_quad",  FALSE, etREAL, {&fc_[ermsQUAD]},
          "Force constant in the penalty function for the magnitude of the quadrupole components." },
        { "-fc_esp",   FALSE, etREAL, {&fc_[ermsESP]},
          "Force constant in the penalty function for the magnitude of the electrostatic potential." },
        { "-fc_charge",  FALSE, etREAL, {&fc_[ermsCHARGE]},
          "Force constant in the penalty function for 'unchemical' charges, i.e. negative hydrogens, and positive oxygens." },
        { "-fc_polar",  FALSE, etREAL, {&fc_[ermsPolar]},
          "Force constant in the penalty function for polarizability." },
        { "-hfac",  FALSE, etREAL, {&hfac_},
          "[HIDDEN]Fudge factor to scale the J00 of hydrogen by (1 + hfac * qH). Default hfac is 0, means no fudging." },
        { "-opthfac",  FALSE, etBOOL, {&bOptHfac_},
          "[HIDDEN]Optimize the fudge factor to scale the J00 of hydrogen (see above). If set, then [TT]-hfac[tt] set the absolute value of the largest hfac. Above this, a penalty is incurred." }
    };
    t_pargs pa_fc[] =
    {
        { "-fc_epot",  FALSE, etREAL, {&fc_[ermsEPOT]},
          "Force constant in the penalty function for the magnitude of the potential energy." },
        { "-fc_force",  FALSE, etREAL, {&fc_[ermsForce2]},
          "Force constant in the penalty function for the magnitude of the force." }
    };
    doAddOptions(pargs, asize(pa_general), pa_general);
    switch (etune)
    {
        case etuneEEM:
            doAddOptions(pargs, asize(pa_eem), pa_eem);
            doAddOptions(pargs, asize(pa_zeta), pa_zeta);
            break;
        case etuneZETA:
            doAddOptions(pargs, asize(pa_zeta), pa_zeta);
            break;
        case etuneFC:
            doAddOptions(pargs, asize(pa_fc), pa_fc);
            break;
        default:
            break;
    }
}

void MolGen::optionsFinished()
{
    hfac0_                      = hfac_;
    cr_                         = init_commrec();
    mdlog_                      = gmx::MDLogger {};
    gmx_omp_nthreads_init(mdlog_, cr_, 1, 1, 1, 0, false, false);
    auto pnc                    = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
    hwinfo_                     = gmx_detect_hardware(mdlog_, pnc);
    if (MASTER(cr_))
    {
        printf("There are %d threads/processes.\n", cr_->nnodes);
    }
    // Make sure all the nodes know about the fitting options
    // Not needed anymore
    // gmx_sum(ermsNR, fc_, cr_);
}

immStatus MolGen::check_data_sufficiency(alexandria::MyMol mymol,
                                         IndexCount       *ic)
{
    immStatus imm = immOK;

    for (int i = 0; i < mymol.atoms_->nr; i++)
    {
        if ((mymol.atoms_->atom[i].atomnumber > 0) &&
            (mymol.atoms_->atom[i].ptype == eptAtom))
        {
            auto fa = pd_.findAtype(*(mymol.atoms_->atomtype[i]));
            if (pd_.getAtypeEnd() != fa)
            {
                const std::string &ztype = fa->getZtype();
                auto               ai    = ic->findName(ztype);
                if (ic->endIndex() == ai)
                {
                    if (debug)
                    {
                        fprintf(debug, "Removing %s because of lacking support for atom %s\n",
                                mymol.getMolname().c_str(),
                                ztype.c_str());
                    }
                    imm = immInsufficientDATA;
                }
            }
            else
            {
                imm = immInsufficientDATA;
            }
        }
    }
    if (imm == immOK)
    {
        for (int i = 0; i < mymol.atoms_->nr; i++)
        {
            if ((mymol.atoms_->atom[i].atomnumber > 0) &&
                (mymol.atoms_->atom[i].ptype == eptAtom))
            {
                auto fa = pd_.findAtype(*(mymol.atoms_->atomtype[i]));
                ic->incrementName(fa->getZtype());
            }
        }
    }
    return imm;
}

void MolGen::Read(FILE            *fp,
                  const char      *fn,
                  const char      *pd_fn,
                  gmx_bool         bZero,
                  char            *opt_elem,
                  const MolSelect &gms,
                  gmx_bool         bCheckSupport,
                  bool             bPairs,
                  bool             bDihedral,
                  bool             bZPE,
                  bool             bFitZeta,
                  bool             bDHform,
                  const char      *tabfn)
{
    int                              nwarn    = 0;
    int                              imm_count[immNR];
    immStatus                        imm      = immOK;
    std::vector<alexandria::MolProp> mp;

    atomprop_  = gmx_atomprop_init();
    for (int i = 0; i < immNR; i++)
    {
        imm_count[i] = 0;
    }
    /* Reading Force Field Data from gentop.dat */
    if (MASTER(cr_))
    {
        GMX_RELEASE_ASSERT(nullptr != pd_fn, "Give me a poldata file name");
        try
        {
            alexandria::readPoldata(pd_fn, pd_, atomprop_);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        if (pd_.getNexcl() != nexcl_)
        {
            fprintf(stderr, "Exclusion number changed from %d in gentop.dat to %d read from the command line.\n",
                    pd_.getNexcl(), nexcl_);
            pd_.setNexcl(nexcl_);
        }
    }
    /* Broadcasting Force Field Data from Master to Slave nodes */
    if (PAR(cr_))
    {
        pd_.broadcast(cr_);
    }
    if (nullptr != fp)
    {
        fprintf(fp, "There are %d atom types in the input file %s:\n---\n",
                static_cast<int>(pd_.getNatypes()), pd_fn);
        fprintf(fp, "---\n\n");
    }
    /* Reading Molecules from allmols.dat */
    if (MASTER(cr_))
    {
        MolPropRead(fn, &mp);
        for (auto mpi = mp.begin(); mpi < mp.end(); )
        {
            mpi->CheckConsistency();
            if (false == mpi->GenerateComposition(&pd_) || imsTrain != gms.status(mpi->getIupac()))
            {
                mpi = mp.erase(mpi);
            }
            else
            {
                ++mpi;
            }
        }
        generate_index(&mp);
    }
    /* Sort Molecules based on the number of atoms */
    if (MASTER(cr_))
    {
        std::sort(mp.begin(), mp.end(),
                  [](alexandria::MolProp &mp1,
                     alexandria::MolProp &mp2)
                  {
                      return (mp1.NAtom() < mp2.NAtom());
                  });
    }
    if (bCheckSupport && MASTER(cr_))
    {
        /* Make a index of eemprop types to be either optimized or
         * being kept constant.
         * TODO: This should probably only be done for tune_eem
         */
        make_index_count(&indexCount_,
                         pd_,
                         opt_elem,
                         bFitZeta);
    }
    /* Generate topology for Molecules and distribute them among the nodes */
    std::string      method, basis;
    splitLot(lot_, &method, &basis);
    int              ntopol    = 0;
    std::vector<int> nmolpar;
    int              nlocaltop = 0;
    if (MASTER(cr_))
    {
        for (auto mpi = mp.begin(); mpi < mp.end(); ++mpi)
        {
            if (imsTrain == gms.status(mpi->getIupac()))
            {
                int               dest = (ntopol % cr_->nnodes);
                alexandria::MyMol mymol;
                if (debug)
                {
                    fprintf(debug, "%s\n", mpi->getMolname().c_str());
                }
                mymol.Merge(mpi);
                mymol.setInputrec(inputrec_);
                imm = mymol.GenerateTopology(atomprop_,
                                             &pd_,
                                             method,
                                             basis,
                                             nullptr,
                                             bGenVsite_,
                                             bPairs,
                                             bDihedral,
                                             false,
                                             tabfn);
                if (bCheckSupport && immOK == imm)
                {
                    imm = check_data_sufficiency(mymol, &indexCount_);
                }
                if (immOK == imm)
                {
                    imm = mymol.GenerateCharges(&pd_,
                                                mdlog_,
                                                atomprop_,
                                                watoms_,
                                                hfac_,
                                                method,
                                                basis,
                                                nullptr,
                                                qsymm_,
                                                nullptr,
                                                cr_,
                                                tabfn,
                                                hwinfo_,
                                                qcycle_,
                                                maxESP_,
                                                qtol_,
                                                nullptr,
                                                nullptr);
                    (void) mymol.espRms();
                }
                if (immOK == imm)
                {
                    imm = mymol.GenerateChargeGroups(ecgGroup, false);
                }
                if (immOK == imm)
                {
                    imm = mymol.getExpProps(bQM_, bZero, bZPE, bDHform,
                                            method, basis, &pd_);
                }
                if (immOK == imm)
                {
                    if (dest > 0)
                    {
                        mymol.eSupp_ = eSupportRemote;
                        if (nullptr != debug)
                        {
                            fprintf(debug, "Going to send %s to cpu %d\n", mpi->getMolname().c_str(), dest);
                        }
                        gmx_send_int(cr_, dest, 1);
                        CommunicationStatus cs = mpi->Send(cr_, dest);
                        if (CS_OK != cs)
                        {
                            imm = immCommProblem;
                        }
                        else
                        {
                            imm = (immStatus)gmx_recv_int(cr_, dest);
                        }
                        if (imm != immOK)
                        {
                            fprintf(stderr, "Molecule %s was not accepted on node %d - error %s\n",
                                    mymol.getMolname().c_str(), dest, alexandria::immsg(imm));
                        }
                        else if (nullptr != debug)
                        {
                            fprintf(debug, "Succesfully beamed over %s\n", mpi->getMolname().c_str());
                        }

                    }
                    else
                    {
                        mymol.eSupp_ = eSupportLocal;
                        nlocaltop   += 1;
                    }
                    if (immOK == imm)
                    {
                        mymol_.push_back(std::move(mymol));
                        ntopol += 1;
                        if (nullptr != debug)
                        {
                            fprintf(debug, "Added %s, ntopol = %d\n", mymol.getMolname().c_str(), ntopol);
                        }
                    }
                }
                if ((immOK != imm) && (nullptr != debug))
                {
                    fprintf(debug, "IMM: Dest: %d %s - %s\n", dest, mpi->getMolname().c_str(), immsg(imm));
                }
            }
            else
            {
                imm = immTest;
            }
            imm_count[imm]++;
        }
        /* Send signal done with transferring molecules */
        for (int i = 1; i < cr_->nnodes; i++)
        {
            gmx_send_int(cr_, i, 0);
        }
        nmolpar.push_back(nlocaltop);
        for (int i = 1; i < cr_->nnodes; i++)
        {
            nmolpar.push_back(gmx_recv_int(cr_, i));
        }
        if (fp)
        {
            for (int i = 0; i < cr_->nnodes; i++)
            {
                fprintf(fp, "Node %d has %d compounds\n", i, nmolpar[i]);
            }
        }
    }
    else
    {
        /***********************************************
         *                                             *
         *           S L A V E   N O D E S             *
         *                                             *
         ***********************************************/
        while (gmx_recv_int(cr_, 0) == 1)
        {
            alexandria::MyMol mymol;
            if (nullptr != debug)
            {
                fprintf(debug, "Going to retrieve new compound\n");
            }
            CommunicationStatus cs = mymol.Receive(cr_, 0);
            if (CS_OK != cs)
            {
                imm = immCommProblem;
            }
            else if (nullptr != debug)
            {
                fprintf(debug, "Succesfully retrieved %s\n", mymol.getMolname().c_str());
                fflush(debug);
            }
            mymol.setInputrec(inputrec_);
            imm = mymol.GenerateTopology(atomprop_,
                                         &pd_,
                                         method,
                                         basis,
                                         nullptr,
                                         bGenVsite_,
                                         bPairs,
                                         bDihedral,
                                         false,
                                         tabfn);

            if (immOK == imm)
            {
                imm = mymol.GenerateCharges(&pd_,
                                            mdlog_,
                                            atomprop_,
                                            watoms_,
                                            hfac_,
                                            method,
                                            basis,
                                            nullptr,
                                            qsymm_,
                                            nullptr,
                                            cr_,
                                            tabfn,
                                            hwinfo_,
                                            qcycle_,
                                            maxESP_,
                                            qtol_,
                                            nullptr,
                                            nullptr);
                (void) mymol.espRms();
            }
            if (immOK == imm)
            {
                imm = mymol.GenerateChargeGroups(ecgAtom, false);
            }
            if (immOK == imm)
            {
                imm = mymol.getExpProps(bQM_, bZero, bZPE, bDHform,
                                        method, basis, &pd_);
            }
            mymol.eSupp_ = eSupportLocal;
            imm_count[imm]++;
            if (immOK == imm)
            {
                mymol_.push_back(std::move(mymol));
                nlocaltop += 1;
                if (nullptr != debug)
                {
                    fprintf(debug, "Added molecule %s. Hform = %g Emol = %g\n",
                            mymol.getMolname().c_str(),
                            mymol.Hform_, mymol.Emol_);
                }
            }
            gmx_send_int(cr_, 0, imm);
        }
        gmx_send_int(cr_, 0, nlocaltop);
    }
    if (fp)
    {
        fprintf(fp, "There were %d warnings because of zero error bars.\n", nwarn);
        fprintf(fp, "Made topologies for %d out of %d molecules.\n",
                ntopol, static_cast<int>(mp.size()));

        for (int i = 0; (i < immNR); i++)
        {
            if (imm_count[i] > 0)
            {
                fprintf(fp, "%d molecules - %s.\n", imm_count[i], alexandria::immsg((immStatus)i));
            }
        }
        if (imm_count[immOK] != (int)mp.size())
        {
            fprintf(fp, "Check alexandria.debug for more information.\nYou may have to use the -debug 1 flag.\n\n");
        }
    }
    if (bCheckSupport && MASTER(cr_))
    {
        indexCount_.cleanIndex(mindata_, fp);
    }
    gmx_sumi(1, &nlocaltop, cr_);
    nmol_support_ = nlocaltop;
    if (nmol_support_ == 0)
    {
        gmx_fatal(FARGS, "No support for any molecule!");
    }
}
}
