/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/fetch_charges.h"
#include "act/alexandria/train_utility.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forces/combinationrules.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/utility/memory_check.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"

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
    { eRMS::Electrostatics, "Electrostatics" },
    { eRMS::Exchange,   "Exchange"   },
    { eRMS::Dispersion, "Dispersion" },
    { eRMS::Induction,  "Induction"  },
    { eRMS::AllElec,    "AllElec"    },
    { eRMS::DeltaHF,    "DeltaHF"    },
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
  if (fp != nullptr && chiSquared_ > 0 && totalWeight_ > 0)
    {
      fprintf(fp, "%-15s  %12.3f  TW: %10g  fc: %10g  weighted: %10g  %s\n",
              rmsName(erms_), chiSquared_, totalWeight_,
              weight_, chiSquaredWeighted(),
              iMolSelectName(ims_));
    }
}

MolGen::MolGen(const CommunicationRecord *cr)
{
    cr_ = cr;
}

void MolGen::addFilenames(std::vector<t_filenm> *filenms)
{
    filenms->push_back({ efXML, "-mp",      "allmols", ffREAD   });
    filenms->push_back({ efXML, "-charges", "charges", ffOPTRD  });
}

void MolGen::addOptions(std::vector<t_pargs>          *pargs,
                        std::map<eRMS, FittingTarget> *targets)
{
    std::vector<t_pargs> pa_general =
    {
        { "-mindata", FALSE, etINT, {&mindata_},
          "Minimum number of data points to optimize force field parameters" },
        { "-qsymm",  FALSE, etBOOL, {&qsymm_},
          "Symmetrize the charges on symmetric groups, e.g. CH3, NH2. The list of groups to symmetrize is specified in the force field file." },
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
        { "-fc_elec",    FALSE, etREAL, {targets->find(eRMS::Electrostatics)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the electrostatic component of the interaction energies of multimers from the reference." },
        { "-fc_allelec",    FALSE, etREAL, {targets->find(eRMS::AllElec)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the sum of all energy components involving electrostatics, that is Coulomb, Polarization, Charge Transfer, of the sum of SAPT energy terms Electrostatics and Induction for multimers from the reference." },
        { "-fc_disp",    FALSE, etREAL, {targets->find(eRMS::Dispersion)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the dispersion component of the interaction energies of multimers from the reference." },
        { "-fc_exch",    FALSE, etREAL, {targets->find(eRMS::Exchange)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the exchange component of the interaction energies of multimers from the reference." },
        { "-fc_induc",    FALSE, etREAL, {targets->find(eRMS::Induction)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the induction component of the interaction energies of multimers from the reference." },
        { "-fc_deltahf",  FALSE, etREAL, {targets->find(eRMS::DeltaHF)->second.weightPtr()},
          "Force constant in the penalty function for the deviation of the DeltaHF term that is part of the SAPT induction energy from the reference. If this option is not present, but the DeltaHF term from SAPT is present in the data, this will be added to the second order induction, which means force field parameters related to induction will be trained on the total induction." },
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
        { "-ener_boltz_temp", FALSE, etREAL, { &ener_boltz_temp_},
          "Apply Boltzmann weighting of energies when using either the [TT]-fc_epot[tt] or [TT]-fc_einter[tt] flags. The weight will be determined from the difference in reference energy at the given data point and in the minimum. " },
        { "-maxpot", FALSE, etINT, { &maxpot_},
          "Fraction of ESP points to use when training on it. Number given in percent." },
        { "-watoms", FALSE, etREAL, { &watoms_},
          "Weight for the potential on the atoms in training on ESP. Should be 0 in most cases, except possibly when using Slater distributions." },
        { "-qtype", FALSE, etSTR, {&qTypeString_},
          "Charge type to use" }
    };
    doAddOptions(pargs, pa_general.size(), pa_general.data());
}

bool MolGen::checkOptions(FILE                        *logFile,
                          const std::vector<t_filenm> &filenames,
                          ForceField                  *pd)
{
    std::set iTypeElec  = { InteractionType::ELECTROSTATICS, InteractionType::ELECTRONEGATIVITYEQUALIZATION,
                            InteractionType::BONDCORRECTIONS, InteractionType::POLARIZATION };
    bool Electrostatics = false;

    for(auto toFit = fit_.begin(); toFit != fit_.end(); )
    {
        InteractionType itype;
        if (pd->typeToInteractionType(toFit->first, &itype))
        {
            Electrostatics = Electrostatics || isVsite(itype) || (iTypeElec.end() != iTypeElec.find(itype));
            toFit++;
        }
        else
        {
            if (debug)
            {
                fprintf(debug, "Ignoring unknown training parameter '%s'\n",
                        toFit->first.c_str());
            }
            toFit = fit_.erase(toFit);
        }
    }
    auto charge_fn = opt2fn_null("-charges", filenames.size(), filenames.data());
    if (Electrostatics && charge_fn && strlen(charge_fn) > 0)
    {
        std::string w("ERROR: supplying a chargemap cannot be combined with changing parameters related to electrostatic interactions.");
        fprintf(logFile, "\n%s\n", w.c_str());
        fprintf(stderr, "\n%s\n\n", w.c_str());
        return false;
    }
    return true;
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

void MolGen::fillIopt(ForceField *pd,
                      FILE       *fp)
// This is called in the read method, the filled structure is used for the optimize() method
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
        else if (fp)
        {
            fprintf(fp, "Cannot find parameter '%s' in force field\n",
                    fit.first.c_str()); 
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
            fprintf(fp, "Will check data sufficiency for %zu training compounds.\n",
                    nmol);
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
            InteractionType::ELECTROSTATICS,
            InteractionType::ELECTRONEGATIVITYEQUALIZATION,
            InteractionType::CHARGE
        };
        // If we use a combination rule for Van der Waals or QT, it should be
        // treated as an atomic type. If not, it should be a "bond"
        // potential.
        for(const auto itype : { InteractionType::VDW, InteractionType::VDWCORRECTION,
                                 InteractionType::INDUCTIONCORRECTION })
        {
            if (!pd->interactionPresent(itype))
            {
                continue;
            }
            auto forces    = pd->findForces(itype);
            auto comb_rule = getCombinationRule(*forces);
            if (!comb_rule.empty())
            {
                atomicItypes.push_back(itype);
            }
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
                for(const auto itype : { InteractionType::VDW, InteractionType::VDWCORRECTION,
                                         InteractionType::INDUCTIONCORRECTION })
                {
                    if (!pd->interactionPresent(itype))
                    {
                        continue;
                    }
                    auto forces    = pd->findForces(itype);
                    auto comb_rule = getCombinationRule(*forces);
                    if (!comb_rule.empty() && optimize(itype))
                    {
                        // This is a hack to allow fitting of a whole matrix of Van der Waals parameters.
                        for(size_t j = i+1; j < myatoms.size(); j++)
                        {
                            auto iPType = pd->findParticleType(myatoms[i].ffType())->interactionTypeToIdentifier(itype).id();
                            auto jPType = pd->findParticleType(myatoms[j].ffType())->interactionTypeToIdentifier(itype).id();
                            auto vdwId  = Identifier({iPType, jPType}, { 1 }, forces->canSwap());
                            if (!forces->parameterExists(vdwId))
                            {
                                GMX_THROW(gmx::InternalError(gmx::formatString("Unknown pair %s-%s for %s",
                                                                               myatoms[i].ffType().c_str(),
                                                                               myatoms[j].ffType().c_str(),
                                                                               interactionTypeToString(itype).c_str()).c_str()));
                            }
                            for(auto &ff : *(forces->findParameters(vdwId)))
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
                                    GMX_THROW(gmx::InternalError(gmx::formatString("Unknown bondcorrection %s", bccId.id().c_str()).c_str()));
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
                InteractionType::VSITE1,
                InteractionType::VSITE2,
                InteractionType::VSITE2FD,
                InteractionType::VSITE3,
                InteractionType::VSITE3S,
                InteractionType::VSITE3FD,
                InteractionType::VSITE3FAD,
                InteractionType::VSITE3OUT,
                InteractionType::VSITE3OUTS
            };
            for (const auto &atype : atypes)
            {
                if (optimize(atype) && pd->interactionPresent(atype))
                {
                    auto angles = pd->findForces(atype);
                    if (top->hasEntry(atype))
                    {
                        for (const auto &topentry : top->entry(atype))
                        {
                            // TODO check multiple ids
                            if (!angles->parameterExists(topentry->id()))
                            {
                                if (fp)
                                {
                                    std::string swapped;
                                    if (topentry->id().haveSwapped())
                                    {
                                        swapped = topentry->id().swapped();
                                    }
                                    fprintf(fp, "Cannot find parameter '%s' CanSwap %s swapped '%s' in parameter list %s CanSwap %s\n",
                                            topentry->id().id().c_str(),
                                            canSwapToString(topentry->id().canSwap()).c_str(),
                                            swapped.c_str(),
                                            interactionTypeToString(atype).c_str(),
                                            canSwapToString(angles->canSwap()).c_str());
                                }
                                continue;
                            }
                            for (auto &ff : *(angles->findParameters(topentry->id())))
                            {
                                if (debug && isVsite(atype))
                                {
                                    fprintf(debug, "Found %s %s\n",
                                            interactionTypeToString(atype).c_str(),
                                            topentry->id().id().c_str());
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
            fprintf(fp, "Found %zu molecules without sufficient support, will remove them.\n",
                    removeMol.size());
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
                {
                    auto mpo = MolPropObservable::POTENTIAL;
                    if (myexp.hasProperty(mpo))
                    {
                        auto eee = myexp.propertyConst(mpo);
                        for (size_t k = 0; k < eee.size(); k++)
                        {
                            auto esp = static_cast<const ElectrostaticPotential *>(eee[k]);
                            w += esp->V().size();
                        }
                    }
                }
                break;
            case eRMS::EPOT:
                // All versus all interactions
                w += gmx::square(myexp.NAtom());
                break;
            case eRMS::Interaction:
            case eRMS::Electrostatics:
            case eRMS::Dispersion:
            case eRMS::Induction:
            case eRMS::DeltaHF:
            case eRMS::Exchange:
            case eRMS::AllElec:
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
                    const std::vector<t_filenm>         &filenms,
                    ForceField                          *pd,
                    const MolSelect                     &gms,
                    const std::map<eRMS, FittingTarget> &targets,
                    bool                                 verbose)
{
    int                              nwarn    = 0;
    std::map<immStatus, int>         imm_count;
    immStatus                        imm      = immStatus::OK;
    std::vector<alexandria::MolProp> mp;
    chargeMap                        qmap;
    auto forceComp = new ForceComputer();
    print_memory_usage(debug);

    //  Now  we have read the forcefield and spread it to processors
    fillIopt(pd, cr_->isMaster() ? fp : nullptr);
    /* Reading Molecules from allmols.dat */
    if (cr_->isMaster())
    {
        auto molfn = opt2fn("-mp", filenms.size(),filenms.data());
        MolPropRead(molfn, &mp);
        auto qmapfn = opt2fn_null("-charges", filenms.size(), filenms.data());
        if (qmapfn && strlen(qmapfn) > 0)
        {
            auto qt = qType::ACM;
            if (nullptr != qTypeString_)
            {
                qt = stringToQtype(qTypeString_);
            }
            qmap = fetchChargeMap(pd, forceComp, qmapfn, qt);
        }
        // Even if we did not read a file, we have to tell the other processors
        // about it.
        broadcastChargeMap(cr_, &qmap);
        if (fp)
        {
            fprintf(fp, "Read %zu compounds from %s\n", mp.size(), molfn);
            if (qmap.size() > 0)
            {
                fprintf(fp, "Read chargemap containing %lu entries from %s\n", qmap.size(), qmapfn);
            }
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
    else
    {
        // Other processors receive the qmap from the master.
        broadcastChargeMap(cr_, &qmap);
    }
    /* Generate topology for Molecules and distribute them among the nodes */
    std::string      method, basis;
    std::map<MolPropObservable, iqmType> iqmMap = 
        {
            { MolPropObservable::DELTAE0,           iqmType::QM },
            { MolPropObservable::POTENTIAL,         iqmType::QM },
            { MolPropObservable::INTERACTIONENERGY, iqmType::QM },
            { MolPropObservable::ELECTROSTATICS,    iqmType::QM },
            { MolPropObservable::INDUCTION,         iqmType::QM },
            { MolPropObservable::EXCHANGE,          iqmType::QM },
            { MolPropObservable::DISPERSION,        iqmType::QM },
            { MolPropObservable::DIPOLE,            iqmType::QM },
            { MolPropObservable::QUADRUPOLE,        iqmType::QM },
            { MolPropObservable::OCTUPOLE,          iqmType::QM },
            { MolPropObservable::HEXADECAPOLE,      iqmType::QM },
            { MolPropObservable::POLARIZABILITY,    iqmType::QM },
            { MolPropObservable::CHARGE,            iqmType::QM }
        };
    int  root  = 0;
    auto alg   = pd->chargeGenerationAlgorithm();
    auto qtype = qType::Calc;
    if (cr_->isMaster())
    {
        if (fp)
        {
            fprintf(fp, "Trying to generate topologies for %zu out of %zu molecules!\n",
                    gms.nMol(), mp.size());
        }
        for (auto mpi = mp.begin(); mpi < mp.end(); ++mpi)
        {
            iMolSelect ims;
            if (gms.status(mpi->getIupac(), &ims))
            {
                alexandria::ACTMol actmol;
                if (debug)
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
                imm = actmol.getExpProps(pd, iqmMap, watoms_, maxpot_);
                if (immStatus::OK == imm)
                {
                    auto fragments = actmol.fragmentHandler();
                    if (!qmap.empty())
                    {
                        if (fragments->setCharges(qmap))
                        {
                            // Copy charges to the high-level topology as well
                            fragments->fetchCharges(actmol.atoms());
                        }
                        else
                        {
                            imm = immStatus::NoMolpropCharges;
                        }
                    }
                    else
                    {
                        std::vector<double> dummy;
                        std::vector<gmx::RVec> forces(actmol.atomsConst().size());
                        imm = actmol.GenerateCharges(pd, forceComp, alg,
                                                     qtype, dummy, &coords, &forces, true);
                    }
                }
                if (immStatus::OK != imm)
                {
                    if (verbose && fp)
                    {
                        fprintf(fp, "Tried to generate charges for %s. Outcome: %s\n",
                                actmol.getMolname().c_str(), immsg(imm));
                    }
                    continue;
                }
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
        if (actmol_.size() < static_cast<size_t>(1+cr_->nhelper_per_middleman()))
        {
            if (cr_->size() > 1)
            {
                // Send a 0 to signify more is coming
                bcint = 0;
                cr_->bcast(&bcint, mycomms[0]);
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Fewer molecules (%zu) than number of helpers per middle man (%d). Use a larger dataset or fewer helpers.", actmol_.size(), 1+cr_->nhelper_per_middleman()).c_str()));
            }
            else
            {
                GMX_THROW(gmx::InvalidInputError("No molecules to train on"));
            }
        }
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
                    if (CommunicationStatus::OK != cs)
                    {
                        GMX_THROW(gmx::InternalError(gmx::formatString("Could not broadcast %s to middlemen", actmol->getMolname().c_str()).c_str()));
                    }
                    comm_used.insert(myindex);
                }
                if (comm_index > 0)
                {
                    // We only send relevant molecules to the helpers
                    cr_->bcast(&bcint, mycomms[comm_index]);
                    cr_->bcast(&comm_index, mycomms[comm_index]);
                    cs = actmol->BroadCast(cr_, root, mycomms[comm_index]);
                    if (CommunicationStatus::OK != cs)
                    {
                        GMX_THROW(gmx::InternalError(gmx::formatString("Could not broadcast %s to helpers", actmol->getMolname().c_str()).c_str()));
                    }
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

            std::vector<gmx::RVec> coords = actmol.xOriginal();
            if (immStatus::OK == imm)
            {
                imm = actmol.getExpProps(pd, iqmMap, watoms_, maxpot_);
            }
            if (immStatus::OK == imm)
            {
                auto fragments = actmol.fragmentHandler();
                if (fragments->setCharges(qmap))
                {
                    // Copy charges to the high-level topology as well
                    fragments->fetchCharges(actmol.atoms());
                }
                else
                {
                    std::vector<gmx::RVec> forces(actmol.atomsConst().size());
                    std::vector<double> dummy;
                    imm = actmol.GenerateCharges(pd, forceComp, alg,
                                                 qtype, dummy, &coords, &forces, true);
                }
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
                    fprintf(debug, "Computing %s on %s nexperiment %d\n", actmol_.back().getMolname().c_str(),
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
            if (m.support() == eSupport::Local)
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
        fprintf(fp, "Made topologies for %zu out of %zu molecules.\n",
                nTargetSize(iMolSelect::Train)+nTargetSize(iMolSelect::Test)+
                nTargetSize(iMolSelect::Ignore),
                mp.size());

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
