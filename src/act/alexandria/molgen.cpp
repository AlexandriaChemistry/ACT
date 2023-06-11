/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include "molgen.h"

#include <cmath>

#include <map>
#include <set>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/arraysize.h"

#include "alex_modules.h"
#include "act/forces/combinationrules.h"
#include "act/utility/memory_check.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/forcefield/forcefield_xml.h"
#include "tuning_utility.h"

namespace alexandria
{

std::map<eRMS, const char *> ermsNames = {
    { eRMS::BOUNDS,     "BOUNDS"     },
    { eRMS::UNPHYSICAL, "UNPHYSICAL" },
    { eRMS::CHARGE,     "CHARGE"     },
    { eRMS::MU,         "MU"         },
    { eRMS::QUAD,       "QUAD"       },
    { eRMS::OCT,        "OCT"        },
    { eRMS::HEXADEC,    "HEXADEC"    },
    { eRMS::FREQUENCY,  "FREQUENCY"  },
    { eRMS::INTENSITY,  "INTENSITY"  },
    { eRMS::CM5,        "CM5"        },
    { eRMS::ESP,        "ESP"        },
    { eRMS::EPOT,       "EPOT"       },
    { eRMS::Interaction,"Interaction"},
    { eRMS::Force2,     "Force2"     },
    { eRMS::Polar,      "Polar"      },
    { eRMS::TOT,        "TOT"        }
};

const std::map<eRMS, const char *> &geteRMSNames()
{
    return ermsNames;
}

const char *rmsName(eRMS e)
{
  return ermsNames[e];
}

void FittingTarget::print(FILE *fp) const
{
  if (fp != nullptr && chiSquared_ > 0 && numberOfDatapoints_ > 0)
    {
      fprintf(fp, "%-10s  %12.3f  N: %6d  fc: %10g  weighted: %10g  %s\n",
              rmsName(erms_), chiSquared_, numberOfDatapoints_,
              weight_, chiSquaredWeighted(),
              iMolSelectName(ims_));
    }
}

MolGen::MolGen(const CommunicationRecord *cr)
{
    cr_        = cr;
}

void MolGen::addOptions(std::vector<t_pargs>          *pargs,
                        std::map<eRMS, FittingTarget> *targets)
{
    t_pargs pa_general[] =
    {
        { "-mindata", FALSE, etINT, {&mindata_},
          "Minimum number of data points to optimize force field parameters" },
        { "-qsymm",  FALSE, etBOOL, {&qsymm_},
          "Symmetrize the charges on symmetric groups, e.g. CH3, NH2. The list of groups to symmetrize is specified in the force field file." },
        { "-qqm",    FALSE, etSTR,  {&chargeMethod_},
          "Use a method from quantum mechanics that needs to be present in the input file. Either ACM, ESP, Hirshfeld, CM5 or Mulliken may be available., depending on your molprop file." },
        { "-zetadiff", FALSE, etREAL, {&zetaDiff_},
          "Difference between core and shell zeta to count as unphysical" },
        { "-lb",  FALSE, etBOOL, {&loadBalance_},
          "Try to divide the computational load evenly over helpers." },
        { "-fit", FALSE, etSTR, {&fitString_},
          "Quoted list of parameters to fit,  e.g. 'alpha zeta'." },
        { "-fc_bound",    FALSE, etREAL, {targets->find(eRMS::BOUNDS)->second.weightPtr()},
          "Force constant in the penalty function for going outside the borders given with the fitting options (see below)." },
        { "-fc_unphysical",  FALSE, etREAL, {targets->find(eRMS::UNPHYSICAL)->second.weightPtr()},
          "Force constant in the penalty function for unphysical combinations of parameters. In particular this parameter penalized that shell zeta is larger the core zeta." },
        { "-fc_epot",    FALSE, etREAL, {targets->find(eRMS::EPOT)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the potential energy of the compound from the reference." },
        { "-fc_inter",    FALSE, etREAL, {targets->find(eRMS::Interaction)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of interaction energies of multimers from the reference." },
        { "-fc_force",  FALSE, etREAL, {targets->find(eRMS::Force2)->second.weightPtr()},
          "Force constant in the penalty function for the magnitude of the squared force on the atoms. For optimized structures the force on the atoms should be zero." },
        { "-fc_freq",  FALSE, etREAL, {targets->find(eRMS::FREQUENCY)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the vibrational frequencies from the reference." },
        { "-fc_inten",  FALSE, etREAL, {targets->find(eRMS::INTENSITY)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the infrared spectrum intensities from the reference." },
        { "-fc_mu",    FALSE, etREAL, {targets->find(eRMS::MU)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the magnitude of the dipole components from the reference." },
        { "-fc_quad",  FALSE, etREAL, {targets->find(eRMS::QUAD)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the magnitude of the quadrupole components from the reference." },
        { "-fc_oct",  FALSE, etREAL, {targets->find(eRMS::OCT)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the magnitude of the octupole components from the reference." },
        { "-fc_hexadec",  FALSE, etREAL, {targets->find(eRMS::HEXADEC)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the magnitude of the hexedecapole components from the reference." },
        { "-fc_esp",   FALSE, etREAL, {targets->find(eRMS::ESP)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the magnitude of the electrostatic potential from the reference." },
        { "-fc_charge",  FALSE, etREAL, {targets->find(eRMS::CHARGE)->second.weightPtr()},
          "Force constant in the penalty function for 'unchemical' charges as determined by the bounds set in the force field file. This is applied to charges that are determined using the Alexandria Charge Model only." },
        { "-fc_cm5",  FALSE, etREAL, {targets->find(eRMS::CM5)->second.weightPtr()},
          "Force constant in the penalty function for deviation from CM5 charges." },
        { "-fc_polar",  FALSE, etREAL, {targets->find(eRMS::Polar)->second.weightPtr()},
          "Force constant in the penalty function for deviation of the six independent components of the molecular polarizability tensor from the reference." },
    };
    doAddOptions(pargs, asize(pa_general), pa_general);
}

void MolGen::optionsFinished()
{
    if (nullptr != fitString_)
    {
        for(const auto &toFit : gmx::splitString(fitString_))
        {
            fit_.insert({ toFit, true });
        }
    }
    if (cr_->isMaster())
    {
        printf("There are %d threads/processes and %zu parameter types to optimize.\n", cr_->size(), fit_.size());
    }
    if (debug)
    {
        fprintf(debug, "optionsFinished: mindata = %d \n", mindata_);
    }
}

void MolGen::fillIopt(ForceField *pd) // This is called in the read method, the filled structure is used for the optimize() method
{
    for(const auto &fit : fit_)
    {
        InteractionType itype;
        if (pd->typeToInteractionType(fit.first, &itype))
        {
            iOpt_.insert({ itype, true });
            if (debug)
            {
                fprintf(debug, "Adding parameter %s to fitting\n", fit.first.c_str());
            }
        }
        else
        {
            printf("Cannot find parameter '%s' in force field\n", fit.first.c_str()); 
        }
    }
}

void MolGen::checkDataSufficiency(FILE        *fp,
                                  ForceField  *pd) // Called in read method
{
    size_t nmol = 0;
    if (targetSize_.find(iMolSelect::Train) == targetSize_.end())
    {
        return;
    }
    do
    {
        // We check data sufficiency only for the training set
        nmol = targetSize_.find(iMolSelect::Train)->second;
        if (fp)
        {
            fprintf(fp, "Will check data sufficiency for %d training compounds.\n",
                    static_cast<int>(nmol));
        }
        /* First set the ntrain values for all forces
         * present that should be optimized to zero.
         */
        for(auto &io : iOpt_)
        {
            if (io.second)
            {
                ForceFieldParameterList *fplist;
                if (pd->interactionPresent(io.first))
                {
                    // Loop over interactions
                    fplist = pd->findForces(io.first);
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
                else
                {
                    GMX_RELEASE_ASSERT(io.first == InteractionType::CHARGE,
                                       "Death Horror Programming Error.");
                    // Loop over particles to find mutable charges
                    auto pv = pd->particleTypes();
                    std::string ccc("charge");
                    for (auto &pt : *pv)
                    {
                        if (pt.second.hasParameter(ccc))
                        {
                            auto p = pt.second.parameter(ccc);
                            if (p->isMutable())
                            {
                                p->setNtrain(0);
                            }
                        }
                    }
                }
            }
        }
        
        std::vector<InteractionType> atomicItypes = {
            InteractionType::POLARIZATION,
            InteractionType::COULOMB,
            InteractionType::ELECTRONEGATIVITYEQUALIZATION,
            InteractionType::CHARGE
        };
        // If we use a combination rule for Van der Waals, it should be
        // treated as an atomic type. If not, it should be a "bond"
        // potential.
        auto itype_vdw  = InteractionType::VDW;
        auto forces_vdw = pd->findForces(itype_vdw);
        int  comb_rule  = getCombinationRule(*forces_vdw);
        if (comb_rule != eCOMB_NONE)
        {
            atomicItypes.push_back(InteractionType::VDW);
        }
        
        // Now loop over molecules and add interactions
        for(auto &mol : actmol_)
        {
            // We can only consider molecules in the training set here
            // since it is only for those that we will update the
            // parameters. That means that if there is no support in the
            // training set for a compound, we can not compute charges
            // or other parameters for the Test or Ignore set either.
            if (mol.datasetType() != iMolSelect::Train)
            {
                continue;
            }
            auto myatoms = mol.topology()->atoms();
            for(size_t i = 0; i < myatoms.size(); i++)
            {
                for(auto &itype : atomicItypes)
                {
                    if (optimize(itype))
                    {
                        auto atype  = pd->findParticleType(myatoms[i].ffType());
                        if (atype->hasInteractionType(itype))
                        {
                            auto fplist = pd->findForces(itype);
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
                        else if (itype == InteractionType::CHARGE)
                        {
                            std::string ccc("charge");
                            if (atype->hasParameter(ccc))
                            {
                                auto p = atype->parameter(ccc);
                                if (p->isMutable())
                                {
                                    p->incrementNtrain();
                                }
                            }
                        }
                    }
                }
                if (comb_rule == eCOMB_NONE && optimize(itype_vdw))
                {
                    // This is a hack to allow fitting of a whole matrix of Van der Waals parameters.
                    for(size_t j = i+1; j < myatoms.size(); j++)
                    {
                        auto iPType = pd->findParticleType(myatoms[i].ffType())->interactionTypeToIdentifier(itype_vdw).id();
                        auto jPType = pd->findParticleType(myatoms[j].ffType())->interactionTypeToIdentifier(itype_vdw).id();
                        auto vdwId  = Identifier({iPType, jPType}, { 1 }, forces_vdw->canSwap());
                        if (!forces_vdw->parameterExists(vdwId))
                        {
                            GMX_THROW(gmx::InternalError("Unknown Van der Waals pair"));
                        }
                        for(auto &ff : *(forces_vdw->findParameters(vdwId)))
                        {
                            if (ff.second.isMutable())
                            {
                                ff.second.incrementNtrain();
                            }
                        }
                    }
                }
            }
            // Now check bonds and bondcorrections
            auto top   = mol.topology();
            auto btype = InteractionType::BONDS;
            if (pd->interactionPresent(btype))
            {
                auto bonds = pd->findForces(btype);
                if (top->hasEntry(btype))
                {
                    for (const auto &topentry : top->entry(btype))
                    {
                        if (optimize(btype))
                        {
                            // TODO check the loop over multiple ids
                            for(auto &ff : *(bonds->findParameters(topentry->id())))
                            {
                                if (ff.second.isMutable())
                                {
                                    ff.second.incrementNtrain();
                                }
                            }
                        }
                        
                        auto bcctype = InteractionType::BONDCORRECTIONS;
                        if (optimize(bcctype) && pd->interactionPresent(bcctype))
                        {
                            int ai = topentry->atomIndex(0);
                            int aj = topentry->atomIndex(1);
                            auto ztype  = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
                            auto iPType = pd->findParticleType(myatoms[ai].ffType())->interactionTypeToIdentifier(ztype).id();
                            auto jPType = pd->findParticleType(myatoms[aj].ffType())->interactionTypeToIdentifier(ztype).id();
                            auto bcc   = pd->findForces(bcctype);
                            auto bccId = Identifier({iPType, jPType}, topentry->bondOrders(), bcc->canSwap());
                            if (!bcc->parameterExists(bccId))
                            {
                                bccId = Identifier({jPType, iPType}, topentry->bondOrders(), bcc->canSwap());
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
            }
            // Now angles and dihedrals
            std::vector<InteractionType> atypes = {
                InteractionType::ANGLES, InteractionType::LINEAR_ANGLES,
                InteractionType::PROPER_DIHEDRALS,
                InteractionType::IMPROPER_DIHEDRALS,
                InteractionType::VSITE2
            };
            for (const auto &atype : atypes)
            {
                if (fp && InteractionType::VSITE2 == atype)
                {
                    fprintf(fp, "Looking for vsite2 in force field\n");
                }
                if (optimize(atype) && pd->interactionPresent(atype))
                {
                    auto angles = pd->findForces(atype);
                    if (fp && InteractionType::VSITE2 == atype)
                    {
                        fprintf(fp, "Looking for vsite2 in molecular tpology\n");
                        top->dump(fp);
                    }
                    if (top->hasEntry(atype))
                    {
                        for (const auto &topentry : top->entry(atype))
                        {
                            // TODO check multiple ids
                            if (fp && InteractionType::VSITE2 == atype)
                            {
                                fprintf(fp, "Looking for vsite2 %s\n", topentry->id().id().c_str());
                            }
                            for (auto &ff : *(angles->findParameters(topentry->id())))
                            {
                                if (fp && InteractionType::VSITE2 == atype)
                                {
                                    fprintf(fp, "Found vsite2 %s\n", topentry->id().id().c_str());
                                }
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
        for(auto &mol : actmol_)
        {
            // Now we check all molecules, including the Test and Ignore
            // set.
            bool keep = true;
            auto myatoms = mol.topology()->atoms();
            for(size_t i = 0; i < myatoms.size(); i++)
            {
                auto atype = pd->findParticleType(myatoms[i].ffType());
                for(auto &itype : atomicItypes)
                {
                    if (optimize(itype))
                    {
                        if (atype->hasInteractionType(itype))
                        {
                            auto fplist = pd->findForces(itype);
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
                        else if (itype == InteractionType::CHARGE)
                        {
                            std::string ccc("charge");
                            if (atype->hasParameter(ccc))
                            {
                                auto p = atype->parameter(ccc);
                                if (p->isMutable() && p->ntrain() < mindata_)
                                {
                                    if (fp)
                                    {
                                        fprintf(fp, "No support for %s - %s in %s. Ntrain is %d, should be at least %d.\n",
                                                ccc.c_str(),
                                                interactionTypeToString(itype).c_str(),
                                                mol.getMolname().c_str(),
                                                p->ntrain(),
                                                mindata_);
                                    }
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
        if (fp && removeMol.size() > 0)
        {
            fprintf(fp, "Found %d molecules without sufficient support, will remove them.\n",
                    static_cast<int>(removeMol.size()));
        }
        for(auto &rmol : removeMol)
        {
            auto moliter = std::find_if(actmol_.begin(), actmol_.end(),
                                        [rmol](ACTMol const &f)
                                        { return (rmol == f.getIupac()); });
            if (moliter == actmol_.end())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find %s in actmol", rmol.c_str()).c_str()));
            }
            if (fp)
            {
                fprintf(fp, "Removing %s from %s set because of lacking support\n",
                        rmol.c_str(),
                        iMolSelectName(moliter->datasetType()));
            }
            actmol_.erase(moliter);
        }
    }
    while (actmol_.size() > 0 && actmol_.size() < nmol);
    if (fp)
    {
        fprintf(fp, "There are %zu molecules left to optimize the parameters for.\n", nmol);
    }
}

static void incrementImmCount(std::map<immStatus, int> *map, immStatus imm) // Called in read method
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

void MolGen::countTargetSize() // Called in read method
{
    targetSize_.clear();
    for(auto &m : actmol_)
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

static double computeCost(const ACTMol                         *actmol,
                          const std::map<eRMS, FittingTarget> &targets)
{
    // Estimate the cost of computing the chi squared for this compound.
    // Since the number of experiments can vary between molecules,
    // so can the cost. We take the cost to be proportional to the number
    // coulomb calculations.
    double w = 1;
    for(auto &myexp : actmol->experimentConst())
    {
        for(const auto &tg: targets)
        {
            switch(tg.first)
            {
            case eRMS::ESP:
                w += myexp.NPotential();
                break;
            case eRMS::EPOT:
                // All versus all interactions
                w += gmx::square(myexp.NAtom());
                break;
            case eRMS::Interaction:
                // All versus all interactions for dimer and monomers separately
                w += 2*gmx::square(myexp.NAtom());
                break;
            case eRMS::Force2:
                // Since we typically do both energy and force, this one is for free
                break;
            case eRMS::Polar:
                // Does 7 energy calculations per polarizability
                w += 7*gmx::square(myexp.NAtom());
                break;
            case eRMS::BOUNDS:
            case eRMS::UNPHYSICAL:
                // Too small to measure
                break;
            case eRMS::CHARGE:
            case eRMS::MU:
            case eRMS::QUAD:
            case eRMS::OCT:
            case eRMS::HEXADEC:
            case eRMS::FREQUENCY:
            case eRMS::INTENSITY:
            case eRMS::CM5:
            case eRMS::TOT:
                break;
            }
        }
    }
    return w;
}
                     
size_t MolGen::Read(FILE                                *fp,
                    const char                          *fn,
                    ForceField                          *pd,
                    const MolSelect                     &gms,
                    const std::map<eRMS, FittingTarget> &targets,
                    bool                                 verbose)
{
    int                              nwarn    = 0;
    std::map<immStatus, int>         imm_count;
    immStatus                        imm      = immStatus::OK;
    std::vector<alexandria::MolProp> mp;
    
    auto forceComp = new ForceComputer();
    print_memory_usage(debug);

    //  Now  we have read the forcefield and spread it to processors
    fillIopt(pd);
    /* Reading Molecules from allmols.dat */
    if (cr_->isMaster())
    {
        MolPropRead(fn, &mp);
        if (fp)
        {
            fprintf(fp, "Read %d compounds from %s\n", static_cast<int>(mp.size()), fn);
        }
        print_memory_usage(debug);
        for (auto mpi = mp.begin(); mpi < mp.end(); )
        {
            mpi->generateComposition();
            mpi->checkConsistency();
            if (!mpi->hasAllAtomTypes())
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

        print_memory_usage(debug);
    }
    /* Generate topology for Molecules and distribute them among the nodes */
    std::string      method, basis;
    std::map<iMolSelect, int> nLocal;
    for(const auto &ims : iMolSelectNames())
    {
        nLocal.insert(std::pair<iMolSelect, int>(ims.first, 0));
    }
    std::map<MolPropObservable, iqmType> iqmMap = 
        {
            { MolPropObservable::DELTAE0,        iqmType::QM },
            { MolPropObservable::DIPOLE,         iqmType::QM },
            { MolPropObservable::QUADRUPOLE,     iqmType::QM },
            { MolPropObservable::OCTUPOLE,       iqmType::QM },
            { MolPropObservable::HEXADECAPOLE,   iqmType::QM },
            { MolPropObservable::POLARIZABILITY, iqmType::QM },
            { MolPropObservable::CHARGE,         iqmType::QM }
        };
    int  root  = 0;
    auto alg   = pd->chargeGenerationAlgorithm();
    auto qtype = qType::Calc;
    if (nullptr != chargeMethod_)
    {
        qtype = stringToQtype(chargeMethod_);
        alg   = ChargeGenerationAlgorithm::Read;
    }
    if (cr_->isMaster())
    {
        if (fp)
        {
            fprintf(fp, "Trying to generate topologies for %d molecules!\n",
                    static_cast<int>(mp.size()));
        }
        for (auto mpi = mp.begin(); mpi < mp.end(); ++mpi)
        {
            iMolSelect ims;
            if (gms.status(mpi->getIupac(), &ims))
            {
                alexandria::ACTMol actmol;
                if (fp && debug)
                {
                    fprintf(debug, "%s\n", mpi->getMolname().c_str());
                }
                actmol.Merge(&(*mpi));
                imm = actmol.GenerateTopology(fp, pd,
                                             missingParameters::Error);
                if (immStatus::OK != imm)
                {
                    if (verbose && fp)
                    {
                        fprintf(fp, "Tried to generate topology for %s. Outcome: %s\n",
                                actmol.getMolname().c_str(), immsg(imm));
                    }
                    continue;
                }

                std::vector<gmx::RVec> coords = actmol.xOriginal();
                actmol.symmetrizeCharges(pd, qsymm_, nullptr);
                actmol.initQgenResp(pd, coords, 0.0, 100);
                std::vector<double> dummy;
                std::vector<gmx::RVec> forces(actmol.atomsConst().size());
                imm = actmol.GenerateCharges(pd, forceComp, alg,
                                            qtype, dummy, &coords, &forces);

                if (immStatus::OK != imm)
                {
                    if (verbose && fp)
                    {
                        fprintf(fp, "Tried to generate charges for %s. Outcome: %s\n",
                                actmol.getMolname().c_str(), immsg(imm));
                    }
                    continue;
                }
                imm = actmol.getExpProps(iqmMap, 0);
                if (immStatus::OK != imm)
                {
                    if (verbose && fp)
                    {
                        fprintf(fp, "Warning: Tried to extract experimental reference data for %s. Outcome: %s\n",
                                actmol.getMolname().c_str(), immsg(imm));
                    }
                }
                else
                {
                    // Only include the compound if we have all data.
                    actmol.set_datasetType(ims);

                    // actmol_ contains all molecules
                    actmol_.push_back(std::move(actmol));
                }
                incrementImmCount(&imm_count, imm);
            }
            else if (verbose && fp)
            {
                fprintf(fp, "Could not find %s in selection.\n",
                        mpi->getIupac().c_str());
            }
        }
        print_memory_usage(debug);
        countTargetSize();
        // Now distribute the molecules over processors.
        // TODO: Make sure the master has a bit less work to do
        // than the helpers and that in particular train
        // compounds are distributed equally otherwise.
        std::map<int, MPI_Comm> mycomms;
        int start = 1;
        if (cr_->nmiddlemen() > 1)
        {
            start = 0;
        }
        for (int i = start; i < start + 1 + cr_->nhelper_per_middleman(); i++)
        {
            mycomms.insert({ i, cr_->create_column_comm(i) });
        }
        int  actmolIndex = 0;
        auto cs         = CommunicationStatus::OK;
        int  bcint      = 1;
        std::set<int> comm_used;
        std::vector<double> totalCost;
        totalCost.resize(1+cr_->nhelper_per_middleman(), 0);
        for(auto actmol = actmol_.begin(); actmol < actmol_.end() && CommunicationStatus::OK == cs; ++actmol)
        {
            // Local compounds are computed locally
            actmol->setSupport(eSupport::Local);
            // See if we need to send stuff around
            if (mycomms.size() > 0)
            {
                int comm_index;
                if (loadBalance_)
                {
                    comm_index = std::min_element(totalCost.begin(), totalCost.end()) - totalCost.begin();
                }
                else
                {
                    comm_index = actmolIndex % (1+cr_->nhelper_per_middleman());
                }
                totalCost[comm_index] += computeCost(&(*actmol), targets);
                if (cr_->nmiddlemen() > 1)
                {
                    // We always need to send the molecules to the middlemen
                    int myindex = 0;
                    // Send a 1 to signify more is coming
                    cr_->bcast(&bcint, mycomms[myindex]);
                    // Send the index of the target group
                    cr_->bcast(&comm_index, mycomms[myindex]);
                    // Send the molecule
                    cs = actmol->BroadCast(cr_, root, mycomms[myindex]);
                    comm_used.insert(myindex);
                }
                if (comm_index > 0)
                {
                    // We only send relevant molecules to the helpers
                    cr_->bcast(&bcint, mycomms[comm_index]);
                    cr_->bcast(&comm_index, mycomms[comm_index]);
                    cs = actmol->BroadCast(cr_, root, mycomms[comm_index]);
                    comm_used.insert(comm_index);
                    actmol->setSupport(eSupport::Remote);
                }
            }
            if (fp)
            {
                fprintf(fp, "Computing %s on %s nexp = %zu\n", actmol->getMolname().c_str(),
                        eSupport::Local == actmol->support() ? "master" : "helper",
                        actmol->experimentConst().size());
            }

            actmolIndex += 1;
        }
        if (CommunicationStatus::OK != cs)
        {
            GMX_THROW(gmx::InternalError("Error broadcasting molecules"));
        }
        bcint = 0;
        for (auto &cc : comm_used)
        {
            cr_->bcast(&bcint, mycomms[cc]);
        }
        // TODO: Free mycomms
        print_memory_usage(debug);
        // Print cost per helper
        if (fp)
        {
            fprintf(fp, "Computational cost estimate per helper:\n");
            for(size_t i = 0; i < totalCost.size(); i++)
            {
                fprintf(fp,"Helper %2lu  Cost %g\n", i, totalCost[i]);
            }
        }
    }
    else
    {
        /***********************************************
         *                                             *
         *  H E L P E R S & M I D D L E M E N          *
         *                                             *
         ***********************************************/
        // Generate the correct communicator for this helper node.
        auto mycomm     = cr_->create_column_comm(cr_->rank() % (1 + cr_->nhelper_per_middleman()));
        int  bcint      = 0;
        cr_->bcast(&bcint, mycomm);
        while (1 == bcint)
        {
            alexandria::ACTMol actmol;
            if (nullptr != debug)
            {
                fprintf(debug, "Going to retrieve new compound\n");
            }
            int comm_index = 0;
            cr_->bcast(&comm_index, mycomm);
            CommunicationStatus cs = actmol.BroadCast(cr_, root, mycomm);
            if (CommunicationStatus::OK != cs)
            {
                imm = immStatus::CommProblem;
            }
            else if (nullptr != debug)
            {
                fprintf(debug, "Succesfully retrieved %s\n", actmol.getMolname().c_str());
                fflush(debug);
            }
            imm = actmol.GenerateTopology(debug, pd, missingParameters::Error);

            if (immStatus::OK == imm)
            {
                std::vector<double> dummy;
                std::vector<gmx::RVec> coords = actmol.xOriginal();
                actmol.symmetrizeCharges(pd, qsymm_, nullptr);
                actmol.initQgenResp(pd, coords, 0.0, 100);
                std::vector<gmx::RVec> forces(actmol.atomsConst().size());
                imm = actmol.GenerateCharges(pd, forceComp, alg,
                                            qtype, dummy, &coords, &forces);
            }
            if (immStatus::OK == imm)
            {
                imm = actmol.getExpProps(iqmMap, 0);
            }
            if (cr_->isMiddleMan())
            {
                actmol.setSupport(eSupport::Remote);
                if (cr_->nhelper_per_middleman() == 0 ||
                    comm_index % (1+cr_->nhelper_per_middleman()) == 0)
                {
                    actmol.setSupport(eSupport::Local);
                }
            }
            else
            {
                actmol.setSupport(eSupport::Local);
            }
            incrementImmCount(&imm_count, imm);
            if (immStatus::OK == imm)
            {
                if (actmol.support() == eSupport::Local)
                {
                    nLocal.find(actmol.datasetType())->second += 1;
                }
                // TODO Checks for energy should be done only when energy is a target for fitting.
                if (false)
                {
                    double deltaE0;
                    if (!actmol.energy(MolPropObservable::DELTAE0, &deltaE0))
                    {
                        if (nullptr != debug)
                        {
                            fprintf(debug, "No DeltaE0 for %s",
                                    actmol.getMolname().c_str());
                        }
                        imm = immStatus::NoData;
                    }
                    if (immStatus::OK == imm)
                    {
                        double hform;
                        if (!actmol.energy(MolPropObservable::DHFORM, &hform))
                        {
                            if (nullptr != debug)
                            {
                                fprintf(debug, "No DeltaHform for %s",
                                        actmol.getMolname().c_str());
                            }
                            imm = immStatus::NoData;
                        }
                        else if (nullptr != debug)
                        {
                            fprintf(debug, "Added molecule %s. Hform = %g DeltaE0 = %g\n",
                                    actmol.getMolname().c_str(), hform, deltaE0);
                        }
                    }
                }
                actmol_.push_back(std::move(actmol));
                if (cr_->isMiddleMan() && debug)
                {
                    int nexp = actmol_.back().experimentConst().size();
                    fprintf(debug, "Computing %s on %s nexepriment %d\n", actmol_.back().getMolname().c_str(),
                            eSupport::Local == actmol.support() ? "middleman" : "helper", nexp);
                }
            }
            // See whether there is more data coming.
            cr_->bcast(&bcint, mycomm);
        }
        countTargetSize();
    }
    // After all molecules have been sent around let's check whether we have
    // support for them, at least on the middlemen but may be needed everywhere.
    // TODO Check that this is correct.
    checkDataSufficiency(fp, pd);
    
    // Some extra debug printing.
    if (debug)
    {
        std::map<iMolSelect, int> nCount;
        nCount.insert({iMolSelect::Train, 0});
        nCount.insert({iMolSelect::Test, 0});
        for(const auto &m : actmol_)
        {
            if (m.support() == eSupport::Remote)
            {
                auto ims = m.datasetType();
                nCount[ims] += 1;
            }
        }
        fprintf(debug, "Node %d Train: %d Test: %d #mols: %zu\n", cr_->rank(), nCount[iMolSelect::Train],
                nCount[iMolSelect::Test], actmol_.size());
    }
    if (fp)
    {
        fprintf(fp, "There were %d warnings because of zero error bars.\n", nwarn);
        fprintf(fp, "Made topologies for %d out of %d molecules.\n",
                static_cast<int>(nTargetSize(iMolSelect::Train)+nTargetSize(iMolSelect::Test)+
                                 nTargetSize(iMolSelect::Ignore)),
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
    return actmol_.size();
}

} // namespace alexandria
