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

#include "molgen.h"

#include <cmath>
#include <map>
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
#include "memory_check.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata_xml.h"

#define STRLEN 256

namespace alexandria
{

  std::map<eRMS, const char *> ermsNames = 
    {
     { eRMS::BOUNDS, "BOUNDS" },
     { eRMS::MU,     "MU" },
     { eRMS::QUAD,   "QUAD" },
     { eRMS::CHARGE, "CHARGE" },
     { eRMS::CM5,    "CM5" },
     { eRMS::ESP,    "ESP" },
     { eRMS::EPOT,   "EPOT" },
     { eRMS::Force2, "Force2" },
     { eRMS::Polar,  "Polar" },
     { eRMS::TOT,    "TOT" }
    };					    

const char *rmsName(eRMS e)
{
  return ermsNames[e];
}

void FittingTarget::print(FILE *fp) const
{
  if (fp != nullptr && chiSquared_ > 0 && numberOfDatapoints_ > 0)
    {
      fprintf(fp, "%-8s  %12.3f  N: %6d  fc: %10g  weighted: %10g  %s\n",
              rmsName(erms_), chiSquared_, numberOfDatapoints_,
              weight_, chiSquaredWeighted(),
              iMolSelectName(ims_));
    }
}
 
MolGen::MolGen()
{
    cr_        = nullptr;
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
    lot_       = "B3LYP/aug-cc-pVTZ";
    inputrec_  = new t_inputrec();
    fill_inputrec(inputrec_);
    for ( const auto &ims : iMolSelectNames()) 
    {
        RmsFittingTarget ft;
        for ( auto &rms : ermsNames )
        {
            FittingTarget    fft(rms.first, ims.first);
            ft.insert(std::pair<eRMS, FittingTarget>(rms.first, std::move(fft)));
        }
        targets_.insert(std::pair<iMolSelect, RmsFittingTarget>(ims.first, std::move(ft)));
        auto etot = target(ims.first, eRMS::TOT);
        GMX_RELEASE_ASSERT(etot != nullptr, "Could not find etot");
        etot->setWeight(1);
    }
}

MolGen::~MolGen()
{
    if (cr_)
    {
        done_commrec(cr_);
    }
}

FittingTarget *MolGen::target(iMolSelect ims, eRMS rms)
{
    auto itarget = targets_.find(ims);
    if (itarget == targets_.end())
    {
        return nullptr;
    }
    auto ift = itarget->second.find(rms);
    if (ift == itarget->second.end())
    {
        return nullptr;
    }
    return &ift->second;
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
        { "-fc_bound",    FALSE, etREAL, {target(iMolSelect::Train, eRMS::BOUNDS)->weightPtr()},
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
                { "-fc_mu",    FALSE, etREAL, {target(iMolSelect::Train, eRMS::MU)->weightPtr()},
                  "Force constant in the penalty function for the magnitude of the dipole components." },
                { "-fc_quad",  FALSE, etREAL, {target(iMolSelect::Train, eRMS::QUAD)->weightPtr()},
                  "Force constant in the penalty function for the magnitude of the quadrupole components." },
                { "-fc_esp",   FALSE, etREAL, {target(iMolSelect::Train, eRMS::ESP)->weightPtr()},
                  "Force constant in the penalty function for the magnitude of the electrostatic potential." },
                { "-fc_charge",  FALSE, etREAL, {target(iMolSelect::Train, eRMS::CHARGE)->weightPtr()},
                  "Force constant in the penalty function for 'unchemical' charges, i.e. negative hydrogens, and positive oxygens." },
                { "-fc_cm5",  FALSE, etREAL, {target(iMolSelect::Train, eRMS::CM5)->weightPtr()},
                  "Force constant in the penalty function for deviation from CM5 charges." },
                { "-fc_polar",  FALSE, etREAL, {target(iMolSelect::Train, eRMS::Polar)->weightPtr()},
                  "Force constant in the penalty function for polarizability." }
            };
        doAddOptions(pargs, asize(pa_eem), pa_eem);
    }
    else if (etune == etuneFC)
    {
        t_pargs pa_fc[] =
            {
                { "-fc_epot",  FALSE, etREAL, {target(iMolSelect::Train, eRMS::EPOT)->weightPtr()},
                  "Force constant in the penalty function for the magnitude of the potential energy." },
                { "-fc_force",  FALSE, etREAL, {target(iMolSelect::Train, eRMS::Force2)->weightPtr()},
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
    if (nullptr != fitString_)
    {
        for(const auto &toFit : gmx::splitString(fitString_))
        {
            fit_.insert({ toFit, true });
        }
    }
    if (MASTER(cr_))
    {
        printf("There are %d threads/processes and %zu parameter types to optimize.\n", cr_->nnodes, fit_.size());
    }
    if (debug)
    {
        fprintf(debug, "optionsFinished: qtol = %g mindata = %d nexcl = %d qcycle = %d\n",
                qtol_, mindata_, nexcl_, qcycle_);
    }
    // Copy the fitting weights from the training set to the other sets
    for(auto &ims : iMolSelectNames())
    {
        if (ims.first != iMolSelect::Train)
        {
            auto ft = fittingTargets(ims.first);
            if (ft != nullptr)
            {
                for(auto &rms : ermsNames)
                {
                    auto w = target(iMolSelect::Train, rms.first)->weight();
                    target(ims.first, rms.first)->setWeight(w);
                }
            }
        }
    }
}

void MolGen::printChiSquared(FILE *fp, iMolSelect ims) const
{
    if (nullptr != fp && MASTER(commrec()))
    {
        fprintf(fp, "Components of fitting function\n");
        for (const auto &ft : fittingTargetsConst(ims))
        {
            ft.second.print(fp);
        }
        fflush(fp);
    }
}

void MolGen::sumChiSquared(bool parallel, iMolSelect ims)
{
    // Now sum over processors
    if (PAR(commrec()) && parallel)
    {
        auto target = targets_.find(ims);
        for( auto &ft : target->second )
        {
            auto chi2 = ft.second.chiSquared();
            gmx_sum(1, &chi2, commrec());
            ft.second.setChiSquared(chi2);
            auto ndp = ft.second.numberOfDatapoints();
            gmx_sumi(1, &ndp, commrec());
            ft.second.setNumberOfDatapoints(ndp);
        }
    }
    auto etot = target(ims, eRMS::TOT);
    GMX_RELEASE_ASSERT(etot != nullptr, "Cannot find etot");
    etot->reset();
    for (const auto &ft : fittingTargetsConst(ims) )
    {
        if (ft.first != eRMS::TOT)
        { 
            etot->increase(1.0, ft.second.chiSquaredWeighted());
        }
    }
    etot->setNumberOfDatapoints(1);
}

void MolGen::fillIopt()
{
    for(const auto &fit : fit_)
    {
        InteractionType itype;
        if (poldata()->typeToInteractionType(fit.first, &itype))
        {
            iOpt_.insert({ itype, true });
        }
    }
}

void MolGen::checkDataSufficiency(FILE *fp)
{
    size_t nmol = 0;
    do
    {
        // We check data sufficiency only for the training set
        nmol = targetSize_.find(iMolSelect::Train)->second;
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
            InteractionType::VDW,
            InteractionType::POLARIZATION,
            InteractionType::CHARGEDISTRIBUTION,
            InteractionType::ELECTRONEGATIVITYEQUALIZATION
        };
        // Now loop over molecules and add interactions
        for(auto &mol : mymol_)
        {
            auto myatoms = mol.atomsConst();
            for(int i = 0; i < myatoms.nr; i++)
            {
                for(auto &itype : atomicItypes)
                {
                    if (optimize(itype))
                    {
                        auto atype  = pd_.findParticleType(*myatoms.atomtype[i]);
                        auto fplist = pd_.findForces(itype);
                        if (atype->hasInteractionType(itype))
                        {
                            auto subId  = atype->interactionTypeToIdentifier(itype);
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
            // Now check bonds and bondcorrections
            auto btype = InteractionType::BONDS;
            auto bonds = pd_.findForces(btype);
            for(int i = 0; i < mol.ltop_->idef.il[bonds->fType()].nr;
                i+= interaction_function[bonds->fType()].nratoms+1)
            {
                // Skip type, which is the first entry in iatoms
                int ai = mol.ltop_->idef.il[bonds->fType()].iatoms[i+1];
                int aj = mol.ltop_->idef.il[bonds->fType()].iatoms[i+2];
                auto bo     = mol.bondOrder(ai, aj);
                if (optimize(btype))
                {
                    auto iPType = pd_.findParticleType(*myatoms.atomtype[ai])->interactionTypeToIdentifier(btype).id();
                    auto jPType = pd_.findParticleType(*myatoms.atomtype[aj])->interactionTypeToIdentifier(btype).id();
                    auto bondId = Identifier({iPType, jPType}, bo, bonds->canSwap());
                    for(auto &ff : *(bonds->findParameters(bondId)))
                    {
                        if (ff.second.isMutable())
                        {
                            ff.second.incrementNtrain();
                        }
                    }
                }
                auto bcctype = InteractionType::BONDCORRECTIONS;
                if (optimize(bcctype) && pd_.interactionPresent(bcctype))
                {
                    auto ztype  = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
                    auto iPType = pd_.findParticleType(*myatoms.atomtype[ai])->interactionTypeToIdentifier(ztype).id();
                    auto jPType = pd_.findParticleType(*myatoms.atomtype[aj])->interactionTypeToIdentifier(ztype).id();
                    auto bcc   = pd_.findForces(bcctype);
                    auto bccId = Identifier({iPType, jPType}, bo, bcc->canSwap());
                    if (!bcc->parameterExists(bccId))
                    {
                        bccId = Identifier({jPType, iPType}, bo, bcc->canSwap());
                        if (!bcc->parameterExists(bccId))
                        {
                            GMX_THROW(gmx::InternalError("Unknown bondcorrection"));
                        }
                    }
                    for(auto &ff : *(bcc->findParameters(bccId)))
                    {
                        if (ff.second.isMutable())
                        {
                            ff.second.incrementNtrain();
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
            auto myatoms = mol.atomsConst();
            for(int i = 0; i < myatoms.nr; i++)
            {
                auto atype = pd_.findParticleType(*myatoms.atomtype[i]);
                for(auto &itype : atomicItypes)
                {
                    if (optimize(itype))
                    {
                        auto fplist = pd_.findForces(itype);
                        if (atype->hasInteractionType(itype))
                        {
                            auto ztype  = atype->interactionTypeToIdentifier(itype);
                            if (!ztype.id().empty())
                            {
                                for(auto &force : fplist->findParametersConst(ztype))
                                {
                                    if (force.second.isMutable() &&
                                        force.second.ntrain() < mindata_)
                                    {
                                        if (fp)
                                        {
                                            fprintf(fp, "No support for %s - %s in %s. Ntrain is %d, should be at least %d.\n",
                                                    ztype.id().c_str(),
                                                    interactionTypeToString(itype).c_str(),
                                                    mol.getMolname().c_str(),
                                                    force.second.ntrain(),
                                                    mindata_);
                                        }
                                        keep = false;
                                        break;
                                    }
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
                GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find %s in mymol", rmol.c_str()).c_str()));
            }
            if (fp)
            {
                fprintf(fp, "Removing %s because of lacking support\n",
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
                        if (param.second.isMutable() && param.second.ntrain() >= mindata_)
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

void MolGen::countTargetSize()
{
    targetSize_.clear();
    for(auto &m : mymol_)
    {
        auto ims = m.datasetType();
        auto ts  = targetSize_.find(ims);
        if (ts == targetSize_.end())
        {
            targetSize_.insert(std::pair<iMolSelect, size_t>(ims, 1));
        }
        else
        {
            ts->second += 1;
        }
    }
}

void MolGen::Read(FILE            *fp,
                  const char      *fn,
                  const char      *pd_fn,
                  gmx_bool         bZero,
                  const MolSelect &gms,
                  bool             bZPE,
                  bool             bDHform,
                  const char      *tabfn)
{
    int                              nwarn    = 0;
    std::map<immStatus, int>         imm_count;
    immStatus                        imm      = immStatus::OK;
    std::vector<alexandria::MolProp> mp;

    print_memory_usage(fp);
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
        print_memory_usage(fp);
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
        print_memory_usage(fp);
        for (auto mpi = mp.begin(); mpi < mp.end(); )
        {
            mpi->CheckConsistency();
            if (!mpi->GenerateComposition())
            {
                /*
                 * Here we remove a compound if we cannot generate it's
                 * composition.
                 */
                mpi = mp.erase(mpi);
            }
            else
            {
                ++mpi;
            }
        }
        generate_index(&mp);

        /* Sort Molecules based on the number of atoms */
        std::sort(mp.begin(), mp.end(),
                  [](alexandria::MolProp &mp1,
                     alexandria::MolProp &mp2)
                  {
                      return (mp1.NAtom() < mp2.NAtom());
                  });
        print_memory_usage(fp);
    }
    /* Generate topology for Molecules and distribute them among the nodes */
    std::string      method, basis;
    splitLot(lot_, &method, &basis);
    std::map<iMolSelect, int> nLocal;
    for(const auto &ims : iMolSelectNames())
    {
        nLocal.insert(std::pair<iMolSelect, int>(ims.first, 0));
    }

    if (MASTER(cr_))
    {
        printf("Generating molecules!\n");
        for (auto mpi = mp.begin(); mpi < mp.end(); ++mpi)
        {
            iMolSelect ims;
            if (gms.status(mpi->getIupac(), &ims) && 
                (iMolSelect::Train == ims || iMolSelect::Test == ims))
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
                                             missingParameters::Error,
                                             tabfn);
                if (immStatus::OK != imm && debug)
                {
                    fprintf(debug, "Tried to generate topology for %s. Outcome: %s\n",
                            mymol.getMolname().c_str(), immsg(imm));
                }
                if (immStatus::OK == imm)
                {
                    mymol.symmetrizeCharges(&pd_, qsymm_, nullptr);
                    mymol.initQgenResp(&pd_, method, basis, 0.0, maxESP_);
                    std::vector<double> dummy;
                    imm = mymol.GenerateCharges(&pd_,
                                                mdlog_,
                                                cr_,
                                                tabfn,
                                                hwinfo_,
                                                qcycle_,
                                                qtol_,
                                                ChargeGenerationAlgorithm::NONE,
                                                dummy,
                                                lot_);
                    //(void) mymol.espRms(qType::Calc);
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
                    mymol.set_datasetType(ims);

                    // mymol_ contains all molecules
                    mymol_.push_back(std::move(mymol));
                }
            }
        }
        print_memory_usage(fp);
        countTargetSize();
        checkDataSufficiency(fp);
        generateOptimizationIndex(fp);
        // Now distribute the molecules over processors.
        // Make sure the master has a bit less work to do
        // than the slaves and that in particular train
        // compounds are distributed equally otherwise. 
        for(auto &ts: targetSize_)
        {
            std::vector<int> ntsNode;
            int              nts  = nTargetSize(ts.first);
            double           ntsD = 0;
            if (cr_->nnodes > 32)
            {
                ntsNode.push_back(nts/(2*cr_->nnodes));
                ntsD  = (1.0*(nts-ntsNode[0]))/(cr_->nnodes-1);
            }
            else
            {
                ntsNode.push_back(nts/cr_->nnodes);
                ntsD = (1.0*nts)/cr_->nnodes;
            }
            double nTot = ntsNode[0]+0.001;
            for(int i = 1; i < cr_->nnodes; i++)
            {
                nTot += ntsD;
                ntsNode.push_back(std::min(nts, int(std::round(nTot))));
            }
            int mymolIndex = 0;
            int dest       = 0;
            for(auto &mymol : mymol_)
            {
                if (mymol.datasetType() != ts.first)
                {
                    continue;
                }
                if (mymolIndex++ >= ntsNode[dest])
                {
                    dest = std::min(cr_->nnodes-1, dest + 1);
                }
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
                    nLocal.find(mymol.datasetType())->second += 1;
                }
                if ((immStatus::OK != imm) && (nullptr != debug))
                {
                    fprintf(debug, "IMM: Dest: %d %s - %s\n",
                            dest, mymol.getMolname().c_str(), immsg(imm));
                }
                incrementImmCount(&imm_count, imm);
            }
        }
        /* Send signal done with transferring molecules */
        for (int i = 1; i < cr_->nnodes; i++)
        {
            gmx_send_int(cr_, i, 0);
        }
        for (int i = 0; i < cr_->nnodes; i++)
        {
            if (fp)
            {
                fprintf(fp, "Node %2d ", i);
            }
            for(const auto &ims : iMolSelectNames())
            {
                int n = nLocal.find(ims.first)->second;
                if (i > 0)
                {
                    n = gmx_recv_int(cr_, i);
                }
                if (fp)
                {
                    fprintf(fp, " %s: %d", ims.second, n);
                }
            }
            if (fp)
            {
                fprintf(fp, " compounds.\n");
            }
        }
        print_memory_usage(fp);
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
                                         missingParameters::Error,
                                         tabfn);

            if (immStatus::OK == imm)
            {
                std::vector<double> dummy;
                mymol.symmetrizeCharges(&pd_, qsymm_, nullptr);
                mymol.initQgenResp(&pd_, method, basis, 0.0, maxESP_);
                imm = mymol.GenerateCharges(&pd_,
                                            mdlog_,
                                            cr_,
                                            tabfn,
                                            hwinfo_,
                                            qcycle_,
                                            qtol_,
                                            ChargeGenerationAlgorithm::NONE,
                                            dummy,
                                            lot_);
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
                
                nLocal.find(mymol.datasetType())->second += 1;
                if (nullptr != debug)
                {
                    fprintf(debug, "Added molecule %s. Hform = %g Emol = %g\n",
                            mymol.getMolname().c_str(),
                            mymol.Hform_, mymol.Emol_);
                }
            }
            gmx_send_int(cr_, 0, static_cast<int>(imm));
        }
        for(const auto &ims : iMolSelectNames())
        {
            gmx_send_int(cr_, 0, nLocal.find(ims.first)->second);
        }
        countTargetSize();
    }
    if (fp)
    {
        fprintf(fp, "There were %d warnings because of zero error bars.\n", nwarn);
        fprintf(fp, "Made topologies for %d out of %d molecules.\n",
                static_cast<int>(nTargetSize(iMolSelect::Train)+nTargetSize(iMolSelect::Test)),
                static_cast<int>(mp.size()));

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
    nmol_support_ = nLocal.find(iMolSelect::Train)->second;
    if (nmol_support_ == 0)
    {
        gmx_fatal(FARGS, "No support for training any molecule!");
    }
}

} // namespace alexandria
