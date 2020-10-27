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

MolGen::MolGen()
{
    cr_        = nullptr;
    bFinal_    = false;
    bDone_     = false;
    bGenVsite_ = false;
    qsymm_     = false;
    constrain_ = false;
    bQM_       = false;
    watoms_    = 0;
    qtol_      = 1e-6;
    qcycle_    = 500;
    mindata_   = 3;
    nexcl_     = 0;
    nexcl_orig_= nexcl_;
    maxESP_    = 100;
    etune_     = etuneEEM;
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
        { "-fit", FALSE, etSTR, {&fitString_},
          "Quoted list of parameters to fit,  e.g. 'alpha zeta'." },
        { "-nexcl",  FALSE, etINT, {&nexcl_},
          "[HIDDEN]Exclusion number." },
        { "-qm",     FALSE, etBOOL, {&bQM_},
          "[HIDDEN]Use only quantum chemistry results (from the levels of theory below) in order to fit the parameters. If not set, experimental values will be used as reference with optional quantum chemistry results, in case no experimental results are available" }
    };
    doAddOptions(pargs, asize(pa_general), pa_general);
    
    if (etune == etuneEEM)
    {
        t_pargs pa_eem[] =
            {
                { "-watoms", FALSE, etREAL, {&watoms_},
                  "Weight for the atoms when fitting the charges to the electrostatic potential. The potential on atoms is usually two orders of magnitude larger than on other points (and negative). For point charges or single smeared charges use zero. For point+smeared charges 1 is recommended." },
                { "-maxpot", FALSE, etINT, {&maxESP_},
                  "Maximum percent of the electrostatic potential points that will be used to fit partial charges. Note that the input file may have a reduced amount of ESP points compared to the Gaussian output already so do not reduce the amount twice unless you know what you are doing." },
                { "-fc_mu",    FALSE, etREAL, {&fc_[ermsMU]},
                  "Force constant in the penalty function for the magnitude of the dipole components." },
                { "-fc_quad",  FALSE, etREAL, {&fc_[ermsQUAD]},
                  "Force constant in the penalty function for the magnitude of the quadrupole components." },
                { "-fc_esp",   FALSE, etREAL, {&fc_[ermsESP]},
                  "Force constant in the penalty function for the magnitude of the electrostatic potential." },
                { "-fc_charge",  FALSE, etREAL, {&fc_[ermsCHARGE]},
                  "Force constant in the penalty function for 'unchemical' charges, i.e. negative hydrogens, and positive oxygens." },
                { "-fc_polar",  FALSE, etREAL, {&fc_[ermsPolar]},
                  "Force constant in the penalty function for polarizability." }
            };
        doAddOptions(pargs, asize(pa_eem), pa_eem);
    }
    else if (etune == etuneFC)
    {
        t_pargs pa_fc[] =
            {
                { "-fc_epot",  FALSE, etREAL, {&fc_[ermsEPOT]},
                  "Force constant in the penalty function for the magnitude of the potential energy." },
                { "-fc_force",  FALSE, etREAL, {&fc_[ermsForce2]},
                  "Force constant in the penalty function for the magnitude of the force." }
            };
        doAddOptions(pargs, asize(pa_fc), pa_fc);
    }
}

void MolGen::optionsFinished()
{
    cr_                         = init_commrec();
    mdlog_                      = gmx::MDLogger {};
    gmx_omp_nthreads_init(mdlog_, cr_, 1, 1, 1, 0, false, false);
    auto pnc                    = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
    hwinfo_                     = gmx_detect_hardware(mdlog_, pnc);
    for(const auto &toFit : gmx::splitString(fitString_))
    {
        fit_.insert({ toFit, true });
    }
    if (MASTER(cr_))
    {
        printf("There are %d threads/processes and %zu parameter types to optimize.\n", cr_->nnodes, fit_.size());
    }
}

void MolGen::fillIopt()
{
    for(const auto &fit : fit_)
    {
        auto itype = poldata()->typeToInteractionType(fit.first);
        iOpt_.insert({ itype, true });
    }
}

void MolGen::checkDataSufficiency(FILE *fp)
{
    size_t nmol = 0;
    do
    {
        nmol = mymol_.size();
        /* First set the ntrain values for all forces
         * present that should be optimized to zero.
         */
        for(auto &io : iOpt_)
        {
            if (io.second)
            {
                ForceFieldParameterList *fplist = pd_.findForces(io.first);
                for(auto &force : *(fplist->parameters()))
                {
                    for(auto &ff : force.second)
                    {
                        if (ff.second.isMutable())
                        {
                            ff.second.setNtrain(0);
                        }
                    }
                }
            }
        }
        // TODO: Handle bonded interactions
        std::vector<InteractionType> atomicItypes = {
            InteractionType::VDW, InteractionType::POLARIZATION, InteractionType::ELECTRONEGATIVITYEQUALIZATION
        };
        // Now loop over molecules and add interactions
        for(auto &mol : mymol_)
        {
            for(int i = 0; i < mol.atoms_->nr; i++)
            {
                for(auto &itype : atomicItypes)
                {
                    if (optimize(itype))
                    {
                        auto atype  = pd_.findAtype(*mol.atoms_->atomtype[i]);
                        auto fplist = pd_.findForces(itype);
                        auto subId  = atype->id(itype);
                        if (!subId.id().empty())
                        {
                            for(auto &ff : *(fplist->findParameters(subId)))
                            {
                                if (ff.second.isMutable())
                                {
                                    ff.second.incrementNtrain();
                                }
                            }
                        }
                    }
                }
            }
        }
        // Now loop over molecules and remove those without sufficient support
        std::vector<std::string> removeMol;
        for(auto &mol : mymol_)
        {
            bool keep = true;
            for(int i = 0; i < mol.atoms_->nr; i++)
            {
                auto atype = pd_.findAtype(*mol.atoms_->atomtype[i]);
                for(auto &itype : atomicItypes)
                {
                    if (optimize(itype))
                    {
                        auto fplist = pd_.findForces(itype);
                        auto ztype  = atype->id(itype);
                        if (!ztype.id().empty())
                        {
                            for(auto &force : fplist->findParametersConst(ztype))
                            {
                                if (force.second.ntrain() < static_cast<uint64_t>(mindata_))
                                {
                                    keep = false;
                                    break;
                                }
                            }
                        }
                    }
                    if (!keep)
                    {
                        break;
                    }
                }
                if (!keep)
                {
                    break;
                }
            }
            if (!keep)
            {
                removeMol.push_back(mol.getIupac());
            }
        }
        for(auto &rmol : removeMol)
        {
            auto moliter = std::find_if(mymol_.begin(), mymol_.end(),
                                        [rmol](MyMol const &f)
                                        { return (rmol == f.getIupac()); });
            if (moliter == mymol_.end())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find %s in mymol_", rmol.c_str()).c_str()));
            }
            if (debug)
            {
                fprintf(debug, "Removing %s because of lacking support\n",
                        rmol.c_str());
            }
            mymol_.erase(moliter);
        }
    }
    while (mymol_.size() > 0 && mymol_.size() < nmol);
    if (fp)
    {
        fprintf(fp, "There are %zu molecules left to optimize the parameters for.\n", nmol);
    }
}

void MolGen::generateOptimizationIndex(FILE *fp)
{
    for(auto &fs : pd_.forcesConst())
    {
        if (optimize(fs.first))
        {
            for(auto &fpl : fs.second.parametersConst())
            {
                for(auto &param : fpl.second)
                {
                    if (fit(param.first))
                    {
                        if (param.second.isMutable() && param.second.ntrain() >= static_cast<uint64_t>(mindata_))
                        {
                            optIndex_.push_back(OptimizationIndex(fs.first, fpl.first, param.first));
                        }
                    }
                }
            }
        }
    }
    if (fp)
    {
        fprintf(fp, "There are %zu variables to optimize.\n", optIndex_.size());
    }
}

static void incrementImmCount(std::map<immStatus, int> *map, immStatus imm)
{
    if (map->find(imm) == map->end())
    {
        map->insert({imm, 1});
    }
    else
    {
        map->find(imm)->second++;
    }
}

void MolGen::Read(FILE            *fp,
                  const char      *fn,
                  const char      *pd_fn,
                  gmx_bool         bZero,
                  const MolSelect &gms,
                  bool             bZPE,
                  bool             bDHform,
                  const char      *tabfn,
                  iMolSelect       SelectType)
{
    int                              nwarn    = 0;
    std::map<immStatus, int>         imm_count;
    immStatus                        imm      = immStatus::OK;
    std::vector<alexandria::MolProp> mp;

    atomprop_  = gmx_atomprop_init();
    /* Reading Force Field Data from gentop.dat */
    if (MASTER(cr_))
    {
        GMX_RELEASE_ASSERT(nullptr != pd_fn, "Give me a poldata file name");
        try
        {
            alexandria::readPoldata(pd_fn, &pd_);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        if (pd_.getNexcl() != nexcl_ && nexcl_ != nexcl_orig_)
        {
            fprintf(stderr, "WARNING: Changing exclusion number from %d in force field file\n", pd_.getNexcl());
            fprintf(stderr, "         to %d (command line), Please check your output carefully.\n", nexcl_);
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
    //  Now  we have read the poldata and spread it to processors
    fillIopt();
    /* Reading Molecules from allmols.dat */
    if (MASTER(cr_))
    {
        MolPropRead(fn, &mp);
        for (auto mpi = mp.begin(); mpi < mp.end(); )
        {
            mpi->CheckConsistency();
            if (!mpi->GenerateComposition() || SelectType != gms.status(mpi->getIupac()))
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
    /* Generate topology for Molecules and distribute them among the nodes */
    std::string      method, basis;
    splitLot(lot_, &method, &basis);
    std::vector<int> nmolpar;
    int              nlocaltop = 0;
    int              ntopol = 0;
    
    if (MASTER(cr_))
    {
        printf("Generating molecules!\n");
        for (auto mpi = mp.begin(); mpi < mp.end(); ++mpi)
        {
            if (SelectType == gms.status(mpi->getIupac()))
            {
                alexandria::MyMol mymol;
                if (debug)
                {
                    fprintf(debug, "%s\n", mpi->getMolname().c_str());
                }
                mymol.Merge(&(*mpi));
                mymol.setInputrec(inputrec_);
                imm = mymol.GenerateTopology(&pd_,
                                             method,
                                             basis,
                                             nullptr,
                                             bGenVsite_,
                                             false,
                                             optimize(InteractionType::PROPER_DIHEDRALS),
                                             false,
                                             tabfn);
                if (immStatus::OK != imm && debug)
                {
                    fprintf(debug, "Tried to generate topology for %s. Outcome: %s\n",
                            mymol.getMolname().c_str(), immsg(imm));
                }
                if (immStatus::OK == imm)
                {
                    mymol.symmetrizeCharges(&pd_, qsymm_, nullptr);
                    mymol.initQgenResp(&pd_, method, basis, nullptr, 0.0, maxESP_);
                    imm = mymol.GenerateCharges(&pd_,
                                                mdlog_,
                                                cr_,
                                                tabfn,
                                                hwinfo_,
                                                qcycle_,
                                                qtol_);
                    (void) mymol.espRms();
                }
                if (immStatus::OK != imm && debug)
                {
                    fprintf(debug, "Tried to generate charges for %s. Outcome: %s\n",
                            mymol.getMolname().c_str(), immsg(imm));
                }
                if (immStatus::OK == imm)
                {
                    imm = mymol.GenerateChargeGroups(ecgGroup, false);
                }
                if (immStatus::OK == imm)
                {
                    imm = mymol.getExpProps(bQM_, bZero, bZPE, bDHform,
                                            method, basis, &pd_);
                }
                if (immStatus::OK == imm)
                {
                    mymol_.push_back(std::move(mymol));
                }
            }
        }
        checkDataSufficiency(fp);
        generateOptimizationIndex(fp);
        for(auto &mymol : mymol_)
        {
            int dest = (ntopol++ % cr_->nnodes);

            if (dest > 0)
            {
                mymol.eSupp_ = eSupport::Remote;
                if (nullptr != debug)
                {
                    fprintf(debug, "Going to send %s to cpu %d\n", mymol.getMolname().c_str(), dest);
                }
                gmx_send_int(cr_, dest, 1);
                CommunicationStatus cs = mymol.Send(cr_, dest);
                if (CS_OK != cs)
                {
                    imm = immStatus::CommProblem;
                }
                else
                {
                    imm = static_cast<immStatus>(gmx_recv_int(cr_, dest));
                }
                if (imm != immStatus::OK)
                {
                    fprintf(stderr, "Molecule %s was not accepted on node %d - error %s\n",
                            mymol.getMolname().c_str(), dest, alexandria::immsg(imm));
                }
                else if (nullptr != debug)
                {
                    fprintf(debug, "Succesfully beamed over %s\n", mymol.getMolname().c_str());
                }
            }
            else
            {
                mymol.eSupp_ = eSupport::Local;
                nlocaltop   += 1;
            }
            if ((immStatus::OK != imm) && (nullptr != debug))
            {
                fprintf(debug, "IMM: Dest: %d %s - %s\n",
                        dest, mymol.getMolname().c_str(), immsg(imm));
            }
            incrementImmCount(&imm_count, imm);
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
                imm = immStatus::CommProblem;
            }
            else if (nullptr != debug)
            {
                fprintf(debug, "Succesfully retrieved %s\n", mymol.getMolname().c_str());
                fflush(debug);
            }
            mymol.setInputrec(inputrec_);
            imm = mymol.GenerateTopology(&pd_,
                                         method,
                                         basis,
                                         nullptr,
                                         bGenVsite_,
                                         false,
                                         optimize(InteractionType::PROPER_DIHEDRALS),
                                         false,
                                         tabfn);

            if (immStatus::OK == imm)
            {
                mymol.symmetrizeCharges(&pd_, qsymm_, nullptr);
                mymol.initQgenResp(&pd_, method, basis, nullptr, 0.0, maxESP_);
                imm = mymol.GenerateCharges(&pd_,
                                            mdlog_,
                                            cr_,
                                            tabfn,
                                            hwinfo_,
                                            qcycle_,
                                            qtol_);
                (void) mymol.espRms();
            }
            if (immStatus::OK == imm)
            {
                imm = mymol.GenerateChargeGroups(ecgAtom, false);
            }
            if (immStatus::OK == imm)
            {
                imm = mymol.getExpProps(bQM_, bZero, bZPE, bDHform,
                                        method, basis, &pd_);
            }
            mymol.eSupp_ = eSupport::Local;
            incrementImmCount(&imm_count, imm);
            if (immStatus::OK == imm)
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
            gmx_send_int(cr_, 0, static_cast<int>(imm));
        }
        gmx_send_int(cr_, 0, nlocaltop);
    }
    if (fp)
    {
        fprintf(fp, "There were %d warnings because of zero error bars.\n", nwarn);
        fprintf(fp, "Made topologies for %d out of %d molecules.\n",
                ntopol, static_cast<int>(mp.size()));

        for (const auto &imm : imm_count)
        {
            fprintf(fp, "%d molecules - %s.\n", imm.second,
                    alexandria::immsg(imm.first));
        }
        if (imm_count.find(immStatus::OK) != imm_count.end())
        {
            if (imm_count.find(immStatus::OK)->second != static_cast<int>(mp.size()))
            {
                fprintf(fp, "Check alexandria.debug for more information.\nYou may have to use the -debug 1 flag.\n\n");
            }
        }
    }
    gmx_sumi(1, &nlocaltop, cr_);
    nmol_support_ = nlocaltop;
    if (nmol_support_ == 0)
    {
        gmx_fatal(FARGS, "No support for any molecule!");
    }
}
}
