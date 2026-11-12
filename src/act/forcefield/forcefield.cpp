/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#include "forcefield.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "gromacs/fileio/md5.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/textreader.h"

#include "act/forcefield/act_checksum.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/generate_dependent.h"
#include "act/forcefield/potential.h"
#include "act/basics/msg_handler.h"
#include "act/utility/stringutil.h"

namespace alexandria
{

bool ForceField::verifyCheckSum(MsgHandler        *msgHandler,
                                const std::string &checkSum)
{
    bool match = checkSum == checkSum_;
    if (!match && msgHandler)
    {
        msgHandler->write(gmx::formatString("Checksum mismatch in %s. Expected %s Found %s\n",
                filename_.c_str(),
                checkSum_.c_str(), checkSum.c_str()));
    }
    return match;
}

std::vector<std::string> ForceField::info() const
{
    std::vector<std::string> out;
    out.emplace_back("");
    out.emplace_back("Force field information");
    out.emplace_back("=======================");
    out.push_back(gmx::formatString("Filename:    %s", filename_.c_str()));
    out.push_back(gmx::formatString("CheckSum:    %s", checkSum_.c_str()));
    out.push_back(gmx::formatString("TimeStamp:   %s", timeStamp_.c_str()));
    out.push_back(gmx::formatString("Polarizable: %s", polarizable() ? "True" : "False"));
    out.emplace_back("Interactions:");
    for(const auto &fs : forces_)
    {
        out.push_back(gmx::formatString("  %s function %s #entries %zu", interactionTypeToString(fs.first).c_str(),
                                        potentialToString(fs.second.potential()).c_str(),
                                        fs.second.parametersConst().size()));
        for(const auto &opt : fs.second.option())
        {
            out.push_back(gmx::formatString("    option %s value %s", opt.first.c_str(), opt.second.c_str()));
        }
        for(const auto &cr : fs.second.combinationRules())
        {
            out.push_back(gmx::formatString("    parameter %s combination rule %s parameter %g",
                                            cr.first.c_str(),
                                            combinationRuleName(cr.second.rule()).c_str(),
                                            cr.second.ffplConst().value()));
        }
    }
    out.emplace_back("");
    return out;
}

bool ForceField::verifyCheckSum(MsgHandler *msgHandler)
{
    const std::string checkSum = forcefieldCheckSum(this);
    return verifyCheckSum(msgHandler, checkSum);
}

void ForceField::updateCheckSum()
{
    checkSum_ = forcefieldCheckSum(this);
}

void ForceField::updateTimeStamp()
{
    const auto tnow    = std::chrono::system_clock::now();
    auto ttt     = std::chrono::system_clock::to_time_t(tnow);
    std::tm *now = std::localtime(&ttt);
    timeStamp_   = gmx::formatString("%d-%02d-%02d %02d:%02d:%02d",
                                     now->tm_year+1900, 1+now->tm_mon,
                                     now->tm_mday, now->tm_hour,
                                     now->tm_min, now->tm_sec);
}

void ForceField::setFilename(const std::string &fn2)
{
    GMX_RELEASE_ASSERT((!fn2.empty()),
                       "Trying to set empty ForceField filename");
    if (!filename_.empty())
    {
        fprintf(stderr, "Changing ForceField filename_ from %s to %s\n",
                filename_.c_str(), fn2.c_str());

    }
    filename_ = fn2;
}

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Atom STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

bool ForceField::typeToInteractionType(const std::string &type, 
                                       InteractionType   *itype) const
{
    if (type2Itype_.empty())
    {
        // The parameters of particle types, see ParticleType::chargeVS
        type2Itype_.insert({InteractionType::CHARGE, { "charge", zetaName, chargeVSName, zetaVSName } });
        for(const auto &fs : forcesConst())
        {
            auto iType = fs.first;
            std::set<std::string> params;
            if (iType == InteractionType::VDW)
            {
                params.insert("exponent");
            }
            for(const auto &fp : fs.second.parametersConst())
            {
                for(const auto &myfp : fp.second)
                {
                    params.insert(myfp.first);
                }
            }
            type2Itype_.insert({iType, params});
        }
    }
    size_t colon = type.find(":");
    if (colon != std::string::npos)
    {
        const auto ittt = type.substr(0, colon);
        InteractionType itp;
        if (stringToInteractionType(ittt, &itp))
        {
            const auto iii  = type2Itype_.find(itp);
            if (iii == type2Itype_.end())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such interaction type '%s' in force field", ittt.c_str()).c_str()));
            }
            const auto ptype = type.substr(colon+1, type.size());
            if (iii->second.find(ptype) == iii->second.end())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such parameter '%s' for interaction type '%s' in force field",
                                                                   ptype.c_str(), ittt.c_str()).c_str()));
            }
            *itype = itp;
            return true;
        }
        else
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such interaction type '%s' in force field", ittt.c_str()).c_str()));
        }
    }
    else
    {
        for (const auto &tt : type2Itype_)
        {
            if (tt.second.find(type) != tt.second.end())
            {
                *itype = tt.first;
                return true;
            }
        }
        return false;
    }
}

void ForceField::guessChargeGenerationAlgorithm()
{
    if (interactionPresent(InteractionType::ELECTRONEGATIVITYEQUALIZATION))
    {
        //! \todo Add check for number of interactions?
        if (interactionPresent(InteractionType::BONDCORRECTIONS))
        {
            ChargeGenerationAlgorithm_ = ChargeGenerationAlgorithm::SQE;
        }
        else
        {
            ChargeGenerationAlgorithm_ = ChargeGenerationAlgorithm::EEM;
        }
    }
    else
    {
        ChargeGenerationAlgorithm_ =  ChargeGenerationAlgorithm::NONE;
    } 
}

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Polarizability STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

void ForceField::checkForPolarizability()
{
    const auto f = forces_.find(InteractionType::POLARIZATION);
    
    polarizable_ = (f != forces_.end() && f->second.numberOfParameters() > 0);
}

/*! \brief The parameters that scale the energy of an interaction
 * \param[in] pot The potential
 * \return The indices of the parameters, all of them if empty
 */
static std::vector<int> scaleIndices(Potential pot)
{
    std::vector<int> scale;
    switch (pot)
    {
    case Potential::LJ12_6:
        scale = { lj12_6EPSILON };
        break;
    case Potential::LJ8_6:
        scale = { lj8_6EPSILON };
        break;
    case Potential::LJ14_7:
        scale = { lj14_7EPSILON };
        break;
    case Potential::WANG_BUCKINGHAM:
        scale = { wbhEPSILON };
        break;
    case Potential::GENERALIZED_BUCKINGHAM:
        scale = { gbhEPSILON };
        break;
    case Potential::LJ12_6_4:
        scale = { lj12_6_4EPSILON, lj12_6_4GAMMA };
        break;
    case Potential::BUCKINGHAM:
        scale = { bhA, bhC6 };
        break;
    case Potential::TANG_TOENNIES:
        scale = { ttA, ttC6, ttC8, ttC10 };
        break;
    case Potential::TT2b:
        scale = { tt2bA, tt2bC6, tt2bC8, tt2bC10 };
        break;
    case Potential::SLATER_ISA_TT:
        scale = { SIttA, SIttC6, SIttC8, SIttC10 };
        break;
    case Potential::BORN_MAYER:
    case Potential::SLATER_ISA:
        scale = { expA };
        break;
    case Potential::MACDANIEL_SCHMIDT:
        scale = { dexpA1, dexpA2 };
        break;
    case Potential::MORSE_BONDS:
        // As used for induction corrections
        scale = { morseDE, morseD0 };
        break;
    default:
        break;
    }
    return scale;
}

/*! \brief Whether a parameter scales the energy of an interaction
 * \param[in] pot  The potential
 * \param[in] name The name of the parameter
 * \return true if the energy vanishes when the parameter is zero
 */
static bool scalesEnergy(Potential pot, const std::string &name)
{
    const auto &names = potentialToParameterName(pot);
    const auto  scale = scaleIndices(pot);
    bool        scales = scale.empty();
    for(auto k : scale)
    {
        scales = scales || (static_cast<size_t>(k) < names.size() && names[k] == name);
    }
    return scales;
}

/*! \brief Whether an interaction is zero for any distance, also when training
 * \param[in] pot    The potential
 * \param[in] params Its parameters
 * \return true if the parameters that scale the energy are zero and cannot
 *         be trained
 */
static bool vanishes(Potential pot, const ForceFieldParameterMap &params)
{
    for(const auto &p : params)
    {
        bool scales = scalesEnergy(pot, p.first);
        // Some combination rules give NaN for two zero parameters, which
        // means that neither particle has this interaction.
        if (scales && ((p.second.value() != 0 && !std::isnan(p.second.value())) ||
                       p.second.isMutable()))
        {
            return false;
        }
    }
    return true;
}

void ForceField::setParticleTypesOption(const std::string &key, const std::string &value)
{
    try
    {
        if (key == "nexcl")
        {
            nexcl_ = std::stoi(value);
        }
        else if (key == "chargedistribution")
        {
            chargeDistribution_ = name2ChargeDistributionType(value);
        }
        else if (key == "epsilonr")
        {
            epsilonr_ = std::stod(value);
        }
        else
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown option '%s' of the particle types", key.c_str()).c_str()));
        }
    }
    catch (const std::logic_error &)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Invalid value '%s' of option '%s' of the particle types",
                                                           value.c_str(), key.c_str()).c_str()));
    }
}

std::map<std::string, std::string> ForceField::particleTypesOptions() const
{
    std::map<std::string, std::string> options;
    if (nexcl_)
    {
        options.insert({ "nexcl", std::to_string(*nexcl_) });
    }
    if (chargeDistribution_)
    {
        options.insert({ "chargedistribution", chargeDistributionTypeName(*chargeDistribution_) });
    }
    if (epsilonr_)
    {
        options.insert({ "epsilonr", gmx::formatString("%g", *epsilonr_) });
    }
    return options;
}

void ForceField::convertCoulomb(MsgHandler *msghandler)
{
    const auto conflict = [this](const std::string &what)
          {
              GMX_THROW(gmx::InvalidInputError(gmx::formatString("%s in %s", what.c_str(), filename_.c_str()).c_str()));
          };
    // The number of exclusions of the corrections is not used
    for(auto itype : { InteractionType::VDWCORRECTION, InteractionType::INDUCTIONCORRECTION })
    {
        const auto fs = forces_.find(itype);
        if (forces_.end() != fs)
        {
            fs->second.removeOption("nexcl");
        }
    }
    // The number of exclusions of both COULOMB and VANDERWAALS
    for(auto itype : { InteractionType::ELECTROSTATICS, InteractionType::VDW })
    {
        const auto fs = forces_.find(itype);
        if (forces_.end() == fs || !fs->second.optionExists("nexcl"))
        {
            continue;
        }
        int nexcl = std::stoi(fs->second.optionValue("nexcl"));
        if (nexcl_ && *nexcl_ != nexcl)
        {
            conflict(gmx::formatString("The number of exclusions of %s (%d) differs from that of the particle types or of another interaction (%d)",
                                       interactionTypeToString(itype).c_str(), nexcl, *nexcl_));
        }
        nexcl_ = nexcl;
        fs->second.removeOption("nexcl");
    }
    const auto fq = forces_.find(InteractionType::ELECTROSTATICS);
    if (forces_.end() == fq)
    {
        return;
    }
    const auto &fsCoul = fq->second;
    auto        cd     = potentialToChargeDistributionType(fsCoul.potential());
    if (chargeDistribution_ && *chargeDistribution_ != cd)
    {
        conflict(gmx::formatString("The charge distribution of %s (%s) differs from that of the particle types (%s)",
                                   interactionTypeToString(InteractionType::ELECTROSTATICS).c_str(),
                                   chargeDistributionTypeName(cd).c_str(),
                                   chargeDistributionTypeName(*chargeDistribution_).c_str()));
    }
    chargeDistribution_ = cd;
    if (fsCoul.optionExists("epsilonr"))
    {
        double epsilonr = std::stod(fsCoul.optionValue("epsilonr"));
        if (epsilonr_ && *epsilonr_ != epsilonr)
        {
            conflict(gmx::formatString("The relative dielectric constant of %s (%g) differs from that of the particle types (%g)",
                                       interactionTypeToString(InteractionType::ELECTROSTATICS).c_str(),
                                       epsilonr, *epsilonr_));
        }
        epsilonr_ = epsilonr;
    }
    // The widths of the charge distributions, a copy for each particle type
    const auto &zetatype = interactionTypeToParticleSubtype(InteractionType::ELECTROSTATICS);
    for(auto &pt : alexandria_)
    {
        if (!pt.second.hasOption(zetatype))
        {
            continue;
        }
        const auto ztype = pt.second.interactionTypeToIdentifier(InteractionType::ELECTROSTATICS);
        if (!fsCoul.parameterExists(ztype))
        {
            continue;
        }
        const auto &zeta = fsCoul.findParameterTypeConst(ztype, zetaName);
        if (pt.second.hasParameter(zetaName))
        {
            if (pt.second.parameterConst(zetaName).value() != zeta.value())
            {
                conflict(gmx::formatString("The zeta of particle type %s (%g) differs from that of %s %s (%g)",
                                           pt.first.id().c_str(), pt.second.parameterConst(zetaName).value(),
                                           zetatype.c_str(), ztype.id().c_str(), zeta.value()));
            }
            continue;
        }
        pt.second.addForceFieldParameter(zetaName, zeta);
        if (msghandler)
        {
            msghandler->msg(ACTStatus::Debug,
                            gmx::formatString("Particle type %s gets zeta %g of %s %s",
                                              pt.first.id().c_str(), zeta.value(),
                                              zetatype.c_str(), ztype.id().c_str()));
        }
    }
    // The option zetatype can only be read, and COULOMB is not needed
    for(auto &pt : alexandria_)
    {
        pt.second.removeOption(zetatype);
    }
    forces_.erase(fq);
    type2Itype_.clear();
}

Potential ForceField::potential(InteractionType itype) const
{
    if (InteractionType::ELECTROSTATICS == itype)
    {
        return chargeDistributionTypeToPotential(chargeDistribution());
    }
    const auto fs = forces_.find(itype);
    return forces_.end() == fs ? Potential::NONE : fs->second.potential();
}

void ForceField::convertVsite1(MsgHandler *msghandler)
{
    const auto fv = forces_.find(InteractionType::VSITE1);
    if (forces_.end() == fv)
    {
        return;
    }
    // The pair parameters made with combination rules have to vanish too,
    // which can only be checked if there are combination rules.
    for(auto itype : { InteractionType::VDW, InteractionType::VDWCORRECTION,
                       InteractionType::INDUCTIONCORRECTION })
    {
        const auto fs = forces_.find(itype);
        if (forces_.end() != fs && !fs->second.parametersConst().empty() &&
            fs->second.combinationRules().empty())
        {
            if (msghandler)
            {
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("Keeping virtual sites, since there are no combination rules for %s",
                                                  interactionTypeToString(itype).c_str()));
            }
            return;
        }
    }
    generateDependentParameter(msghandler, this, true);
    // Interaction types in which a virtual site that only carries a charge
    // may occur, and those of which its parameters have to be zero.
    const std::set<InteractionType> allowed    = { InteractionType::VSITE1, InteractionType::VDW,
                                                   InteractionType::VDWCORRECTION,
                                                   InteractionType::INDUCTIONCORRECTION,
                                                   InteractionType::POLARIZATION,
                                                   InteractionType::ELECTROSTATICS };
    const std::set<InteractionType> mustBeZero = { InteractionType::VDW, InteractionType::VDWCORRECTION,
                                                   InteractionType::INDUCTIONCORRECTION };
    // The names by which the parameters of other interaction types refer to a particle type
    const auto names = [](const ParticleType &pt)
          {
              std::set<std::string> nm = { pt.id().id() };
              for(const auto &opt : pt.optionsConst())
              {
                  if (opt.first != "element" && opt.first != "atomnumber" && opt.first != "row" &&
                      opt.first != rowVSName)
                  {
                      nm.insert(opt.second);
                  }
              }
              return nm;
          };
    // Whether an identifier refers to one of the names
    const auto refers = [](const Identifier &id, const std::set<std::string> &nm)
          {
              for(const auto &a : id.atoms())
              {
                  if (nm.count(a) > 0)
                  {
                      return true;
                  }
              }
              return false;
          };
    // The virtual site type of each atom type, if it can be converted
    std::map<std::string, std::string> vsiteOf;
    // The number of VSITE1 entries of each virtual site type
    std::map<std::string, int>         nuse;
    // The virtual sites that are kept, with the reason
    std::map<std::string, std::string> rejected;
    const auto reject = [&rejected](const std::string &vs, const std::string &why)
          {
              rejected.insert({ vs, why });
          };
    for(const auto &entry : fv->second.parametersConst())
    {
        const auto &atoms = entry.first.atoms();
        if (atoms.size() != 2 || !hasParticleType(atoms[0]) || !hasParticleType(atoms[1]))
        {
            continue;
        }
        nuse[atoms[1]] += 1;
        const auto *pa = findParticleType(atoms[0]);
        const auto *pv = findParticleType(atoms[1]);
        if (pa->apType() != ActParticle::Atom || pv->apType() != ActParticle::Vsite ||
            pa->hasChargeVS() || vsiteOf.count(atoms[0]) > 0)
        {
            reject(atoms[1], "not a single virtual site on an atom");
            continue;
        }
        vsiteOf.insert({ atoms[0], atoms[1] });
    }
    // Check that the virtual sites only carry a charge, and that their names
    // are not used by other particle types.
    std::set<std::string> vsites;
    for(const auto &av : vsiteOf)
    {
        vsites.insert(av.second);
    }
    for(const auto &vs : vsites)
    {
        const auto *pv = findParticleType(vs);
        const auto  nm = names(*pv);
        for(const auto &pt : alexandria_)
        {
            if (vsites.count(pt.first.id()) == 0)
            {
                for(const auto &n : names(pt.second))
                {
                    if (nm.count(n) > 0)
                    {
                        reject(vs, "name " + n + " is used by " + pt.first.id());
                    }
                }
            }
        }
        for(const auto &fs : forces_)
        {
            for(const auto &p : fs.second.parametersConst())
            {
                if (!refers(p.first, nm))
                {
                    continue;
                }
                if (allowed.count(fs.first) == 0)
                {
                    reject(vs, "it has " + interactionTypeToString(fs.first) + " parameters");
                }
                else if (InteractionType::POLARIZATION == fs.first)
                {
                    // The polarization parameters only matter for a shell
                    if (pv->hasInteractionType(InteractionType::POLARIZATION))
                    {
                        reject(vs, "it has a shell");
                    }
                }
                else if (mustBeZero.count(fs.first) > 0 && !vanishes(fs.second.potential(), p.second))
                {
                    reject(vs, "it has " + interactionTypeToString(fs.first) + " interactions with " + p.first.id());
                }
            }
        }
    }
    // Now convert
    std::set<std::string> removeNames;
    std::set<std::string> converted;
    if (msghandler)
    {
        for(const auto &r : rejected)
        {
            msghandler->msg(ACTStatus::Verbose,
                            gmx::formatString("Keeping virtual site %s, since %s", r.first.c_str(), r.second.c_str()));
        }
    }
    for(const auto &av : vsiteOf)
    {
        if (rejected.count(av.second) > 0)
        {
            continue;
        }
        auto *pa = findParticleType(av.first);
        const auto *pv = findParticleType(av.second);
        pa->addForceFieldParameter(chargeVSName, pv->parameterConst("charge"));
        ForceFieldParameter zeta("1/nm", 0, 0, 0, 0, 0, Mutability::Fixed, true, true);
        if (pv->hasParameter(zetaName))
        {
            zeta = pv->parameterConst(zetaName);
        }
        pa->addForceFieldParameter(zetaVSName, zeta);
        pa->setOption(rowVSName, std::to_string(pv->row()));
        auto nm = names(*pv);
        removeNames.insert(nm.begin(), nm.end());
        converted.insert(av.second);
        if (msghandler)
        {
            msghandler->msg(ACTStatus::Verbose,
                            gmx::formatString("Replaced virtual site %s on %s by a second charge %g with zeta %g",
                                              av.second.c_str(), av.first.c_str(),
                                              pa->chargeVS(), pa->zetaVS()));
        }
    }
    // Remove the virtual sites that were converted completely
    for(const auto &vs : converted)
    {
        if (nuse[vs] != std::count_if(vsiteOf.begin(), vsiteOf.end(),
                                      [&vs](const auto &av) { return av.second == vs; }))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Virtual site %s was only partly converted", vs.c_str())));
        }
        alexandria_.erase(Identifier(vs));
        ptypeVersion_.update();
    }
    for(auto &fs : forces_)
    {
        auto params = fs.second.parameters();
        for(auto p = params->begin(); p != params->end(); )
        {
            if (refers(p->first, removeNames))
            {
                p = params->erase(p);
            }
            else
            {
                ++p;
            }
        }
    }
    // Without virtual sites on atoms there may be no VSITE1 left
    if (!converted.empty() && forces_.find(InteractionType::VSITE1)->second.parametersConst().empty())
    {
        forces_.erase(InteractionType::VSITE1);
    }
    if (msghandler && !converted.empty())
    {
        size_t natype = std::count_if(vsiteOf.begin(), vsiteOf.end(),
                                      [&converted](const auto &av) { return converted.count(av.second) > 0; });
        msghandler->msg(ACTStatus::Verbose,
                        gmx::formatString("Replaced %zu virtual site types by second charges of %zu atom types.",
                                          converted.size(), natype));
    }
}

void ForceField::expandSecondCharges(MsgHandler *msghandler)
{
    std::vector<std::string> atypes;
    for(const auto &pt : alexandria_)
    {
        if (pt.second.hasChargeVS())
        {
            atypes.push_back(pt.first.id());
        }
    }
    if (atypes.empty())
    {
        return;
    }
    if (forces_.find(InteractionType::VSITE1) == forces_.end())
    {
        addForces(InteractionType::VSITE1,
                  ForceFieldParameterList(potentialToString(Potential::VSITE1), CanSwap::Yes));
    }
    auto *fv = findForces(InteractionType::VSITE1);
    // Whether a name is used by a particle type or by parameters
    const auto used = [this](const std::string &name)
          {
              if (hasParticleType(name))
              {
                  return true;
              }
              for(const auto &fs : forces_)
              {
                  if (fs.second.parameterExists(Identifier(name)) ||
                      fs.second.parameterExists(Identifier(name + "_z")))
                  {
                      return true;
                  }
              }
              return false;
          };
    for(const auto &atype : atypes)
    {
        auto       *pa = findParticleType(atype);
        std::string vs = "v1" + atype;
        for(int k = 2; used(vs); k++)
        {
            vs = gmx::formatString("v1%s_%d", atype.c_str(), k);
        }
        ParticleType pv(Identifier(vs), "", ActParticle::Vsite);
        pv.setOption("atomnumber", "0");
        pv.setOption("element", "VS");
        pv.setOption("row", std::to_string(pa->rowVS()));
        pv.setOption(interactionTypeToParticleSubtype(InteractionType::BONDS), vs);
        pv.addForceFieldParameter("charge", pa->parameterConst(chargeVSName));
        pv.addForceFieldParameter("mass", ForceFieldParameter("Da", 0, 0, 0, 0, 0, Mutability::Fixed, true, true));
        // The distribution of the second charge
        pv.addForceFieldParameter(zetaName, pa->parameterConst(zetaVSName));
        // Interactions that the virtual site does not have, its parameters
        // are those of the atom with the energy scaled to zero.
        for(auto itype : { InteractionType::VDW, InteractionType::VDWCORRECTION,
                           InteractionType::INDUCTIONCORRECTION })
        {
            const auto fs = forces_.find(itype);
            if (forces_.end() == fs || !pa->hasInteractionType(itype))
            {
                continue;
            }
            const auto pid = pa->interactionTypeToIdentifier(itype);
            if (!fs->second.parameterExists(pid))
            {
                continue;
            }
            auto pmap = fs->second.findParameterMapConst(pid);
            for(auto &p : pmap)
            {
                if (scalesEnergy(fs->second.potential(), p.first))
                {
                    p.second.setMinimum(0);
                    p.second.forceSetValue(0);
                    p.second.setMutability(Mutability::Fixed);
                }
                fs->second.addParameter(Identifier(vs), p.first, p.second);
            }
            pv.setOption(interactionTypeToParticleSubtype(itype), vs);
        }
        addParticleType(pv);
        fv->addParameter(Identifier(InteractionType::VSITE1, atype + "!" + vs, CanSwap::Yes), "vs1a",
                         ForceFieldParameter("", 1, 0, 0, 1, 1, Mutability::Fixed, true, false));
        if (msghandler)
        {
            msghandler->msg(ACTStatus::Verbose,
                            gmx::formatString("Replaced second charge %g with zeta %g of %s by virtual site %s",
                                              pa->chargeVS(), pa->zetaVS(), atype.c_str(), vs.c_str()));
        }
        pa->parameters()->erase(chargeVSName);
        pa->parameters()->erase(zetaVSName);
        pa->removeOption(rowVSName);
    }
    generateDependentParameter(msghandler, this, true);
}

void ForceField::addParticleType(const ParticleType &ptp)
{
    if (hasParticleType(ptp.id()))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Trying to add the same particle type %s twice", ptp.id().id().c_str()).c_str()));
    }
    alexandria_.insert({ptp.id(), ptp});
    ptypeVersion_.update();
}

void ForceField::addForces(InteractionType                iType,
                           const ForceFieldParameterList &forces)
{
    const auto f     = forces_.find(iType);
    
    GMX_RELEASE_ASSERT(f == forces_.end(),
                       gmx::formatString("Will not add a second ForceFieldParameterList for %s\n",
                                         interactionTypeToString(iType).c_str()).c_str());
    forces_.insert({iType, forces});
}

bool ForceField::atypeToBtype(const std::string &atype,
                           std::string       *btype) const
{
    auto ai    = findParticleType(atype);
    auto itype = InteractionType::BONDS;
    if (ai && ai->hasInteractionType(itype))
    {
        btype->assign(ai->interactionTypeToIdentifier(itype).id());
        return true;
    }
    return false;
}

/*
 *-+-+-+-+-+-+-+-+
 * COULOMB STUFF
 *-+-+-+-+-+-+-+-+
 */

void ForceField::addSymcharges(const std::string &central,
                            const std::string &attached,
                            int                numattach)
{
    size_t        i;
    Symcharges  * sc;
    for (i = 0; i < symcharges_.size(); i++)
    {
        sc = &(symcharges_[i]);
        if ((strcasecmp(sc->getCentral().c_str(), central.c_str()) == 0)   &&
            (strcasecmp(sc->getAttached().c_str(), attached.c_str()) == 0) &&
            (sc->getNumattach() == numattach))
        {
            break;
        }
    }
    if (i == symcharges_.size())
    {
        const Symcharges symcharges(central, attached, numattach);
        symcharges_.push_back(symcharges);
    }
}

CommunicationStatus ForceField::Send(const CommunicationRecord *cr, int dest)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send(dest, filename_);
        cr->send(dest, checkSum_);
        cr->send(dest, timeStamp_);
        cr->send(dest, static_cast<int>(ChargeGenerationAlgorithm_));
        cr->send(dest, polarizable_ ? 1 : 0);
        const auto ptOptions = particleTypesOptions();
        cr->send(dest, ptOptions.size());
        for(const auto &opt : ptOptions)
        {
            cr->send(dest, opt.first);
            cr->send(dest, opt.second);
        }
        /* Send Ffatype */
        cr->send(dest, alexandria_.size());
        for (auto &alexandria : alexandria_)
        {
            cs = alexandria.second.Send(cr, dest);
            if (CommunicationStatus::OK != cs)
            {
                break;
            }
        }

        /* Force Field Parameter Lists */
        if (CommunicationStatus::OK == cs)
        {
            cr->send(dest, forces_.size());
            for (auto &force : forces_)
            {
                const std::string key(interactionTypeToString(force.first));
                cr->send(dest, key);
                cs = force.second.Send(cr, dest);
                if (CommunicationStatus::OK != cs)
                {
                    break;
                }
            }
        }

        /* Send Symcharges */
        if (CommunicationStatus::OK == cs)
        {
            cr->send(dest, symcharges_.size());
            for (auto &symcharges : symcharges_)
            {
                cs = symcharges.Send(cr, dest);
                if (CommunicationStatus::OK != cs)
                {
                    break;
                }
            }
        }
    }
    return cs;
}

CommunicationStatus ForceField::BroadCast(const CommunicationRecord *cr,
                                       int                        root,
                                       MPI_Comm                   comm,
                                       MsgHandler                *msgHandler)
{
    CommunicationStatus cs = cr->bcast_data(comm);
    if (CommunicationStatus::OK == cs)
    {
        cr->bcast(&filename_, comm);
        cr->bcast(&checkSum_, comm);
        cr->bcast(&timeStamp_, comm);
        int icq = static_cast<int>(ChargeGenerationAlgorithm_);
        cr->bcast(&icq, comm);
        ChargeGenerationAlgorithm_ = static_cast<ChargeGenerationAlgorithm>(icq);
        int pol = polarizable_ ? 1 : 0;
        cr->bcast(&pol, comm);
        polarizable_ = pol == 1;
        {
            auto   ptOptions = particleTypesOptions();
            size_t nopt      = ptOptions.size();
            cr->bcast(&nopt, comm);
            auto   opt       = ptOptions.begin();
            for(size_t i = 0; i < nopt; i++)
            {
                std::string key, value;
                if (cr->rank() == root)
                {
                    key   = opt->first;
                    value = opt->second;
                    ++opt;
                }
                cr->bcast(&key, comm);
                cr->bcast(&value, comm);
                if (cr->rank() != root)
                {
                    setParticleTypesOption(key, value);
                }
            }
        }
        /* Bcast Ffatype */
        size_t asize = alexandria_.size();
        cr->bcast(&asize, comm);
        if (cr->rank() == root)
        {
            for(auto &aa : alexandria_)
            {
                aa.second.BroadCast(cr, root, comm);
            }
        }
        else
        {
            for(size_t as = 0; as < asize; as++)
            {
                ParticleType pt;
                cs = pt.BroadCast(cr, root, comm);
                if (CommunicationStatus::OK == cs)
                {
                    alexandria_.insert({pt.id(), pt});
                    ptypeVersion_.update();
                }
            }
        }

        /* Force Field Parameter Lists */
        if (CommunicationStatus::OK == cs)
        {
            size_t fsize = forces_.size();
            cr->bcast(&fsize, comm);
            if (cr->rank() == root)
            {
                for(auto &ff : forces_)
                {
                    std::string key(interactionTypeToString(ff.first));
                    cr->bcast(&key, comm);
                    ff.second.BroadCast(cr, root, comm);   
                }
            }
            else
            {
                forces_.clear();
            
                for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < fsize); n++)
                {
                    ForceFieldParameterList fs;
                    std::string             key;
                    cr->bcast(&key, comm);
                    InteractionType iType;
                    if (!stringToInteractionType(key, &iType))
                    {
                        GMX_THROW(gmx::InternalError("Incorrect interaction type"));
                    }
                    cs                    = fs.BroadCast(cr, root, comm);
                    if (CommunicationStatus::OK == cs && cr->rank() != root)
                    {
                        forces_.insert({iType, fs});
                    }
                    if (msgHandler)
                    {
                        msgHandler->writeDebug(gmx::formatString("Done Listed force %s\n", key.c_str()));
                    }
                }
            }
        }

        /* Bcast Symcharges */
        if (CommunicationStatus::OK == cs)
        {
            size_t scsize = symcharges_.size();
            cr->bcast(&scsize, comm);
            if (cr->rank() != root)
            {
                symcharges_.resize(scsize);
            }
            for(size_t scs = 0; scs < scsize; scs++)
            {
                cs = symcharges_[scs].BroadCast(cr, root, comm);
                if (CommunicationStatus::OK != cs)
                {
                    break;
                }
            }
        }
    }
    return cs;
}

CommunicationStatus ForceField::Receive(const CommunicationRecord *cr, int src,
                                        MsgHandler *msgHandler)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv(src, &filename_);
        cr->recv(src, &checkSum_);
        cr->recv(src, &timeStamp_);
        int cga;
        cr->recv(src, &cga);
        ChargeGenerationAlgorithm_ = static_cast<ChargeGenerationAlgorithm>(cga);
        cr->recv(src, &polarizable_);
        size_t nopt;
        cr->recv(src, &nopt);
        nexcl_.reset();
        chargeDistribution_.reset();
        epsilonr_.reset();
        for(size_t i = 0; i < nopt; i++)
        {
            std::string key, value;
            cr->recv(src, &key);
            cr->recv(src, &value);
            setParticleTypesOption(key, value);
        }
        /* Rceive Ffatype */
        size_t nalexandria;
        cr->recv(src, &nalexandria);
        alexandria_.clear();
        ptypeVersion_.update();
        for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < nalexandria); n++)
        {
            ParticleType alexandria;
            cs = alexandria.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                alexandria_.insert({alexandria.id(), alexandria});
            }
        }
        if (msgHandler)
        {
            msgHandler->writeDebug("Done receiving atomtypes\n");
        }

        /* Receive Listed Forces */
        size_t nforces;
        cr->recv(src, &nforces);
        forces_.clear();
        for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < nforces); n++)
        {
            ForceFieldParameterList fs;
            std::string             key;
            cr->recv(src, &key);
            InteractionType iType;
            if (!stringToInteractionType(key, &iType))
            {
                GMX_THROW(gmx::InternalError("Incorrect interaction type"));
            }
            cs                    = fs.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                forces_.insert({iType, fs});
            }
            if (msgHandler)
            {
                msgHandler->writeDebug(gmx::formatString("Done Listed force %s\n", key.c_str()));
            }
        }
        if (msgHandler)
        {
            msgHandler->writeDebug("Done Listed forces\n");
        }
        /* Receive Symcharges */
        if (CommunicationStatus::OK == cs)
        {
            size_t nsymcharges;
            cr->recv(src, &nsymcharges);
            symcharges_.clear();
            for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < nsymcharges); n++)
            {
                Symcharges symcharges;
                cs = symcharges.Receive(cr, src);
                if (CommunicationStatus::OK == cs)
                {
                    symcharges_.push_back(symcharges);
                }
            }
        }
    }
    return cs;
}

void ForceField::sendToHelpers(const CommunicationRecord *cr, int root,
                               MsgHandler *msgHandler)
{
    if (msgHandler)
    {
        std::string msg = gmx::formatString("Will send force field from node %d to helper", cr->rank());
        for(const auto &i: cr->helpers())
        {
            msg += gmx::formatString(" %d", i);
        }
        msg += "\n";
        msgHandler->writeDebug(msg);
    }
    // Packed into a single broadcast, the broadcasts in BroadCast use root 0
    cr->startPackedBcast(cr->send_helpers(), 0);
    BroadCast(cr, root, cr->send_helpers(), msgHandler);
    cr->finishPackedBcast();
}

void ForceField::checkUnits() const
{
    // The charge generation, also on the GPU, uses these as they are, in
    // eV, see QgenAcm, the others are converted with internalValue
    const std::set<std::string> eV      = { "eV", "eV/e" };
    const std::set<std::string> perNm   = { "1/nm" };
    std::string                 errors;
    auto check = [&errors](const std::string         &where,
                           const std::string         &name,
                           const ForceFieldParameter &p,
                           const std::set<std::string> &allowed)
    {
        if (allowed.count(p.unit()) == 0)
        {
            std::string units;
            for (const auto &u : allowed)
            {
                units += (units.empty() ? "'" : " or '") + u + "'";
            }
            errors += gmx::formatString("\n  %s of %s has unit '%s', it has to be %s",
                                        name.c_str(), where.c_str(), p.unit().c_str(), units.c_str());
        }
    };
    for (auto itype : { InteractionType::ELECTRONEGATIVITYEQUALIZATION, InteractionType::BONDCORRECTIONS })
    {
        if (!interactionPresent(itype))
        {
            continue;
        }
        for (const auto &entry : findForcesConst(itype).parametersConst())
        {
            for (const auto &p : entry.second)
            {
                if (p.first == "eta" || p.first == "chi" || p.first == "delta_eta" || p.first == "delta_chi")
                {
                    check(interactionTypeToString(itype) + " " + entry.first.id(), p.first, p.second, eV);
                }
            }
        }
    }
    for (const auto &pt : alexandria_)
    {
        for (const char *zname : { zetaName, zetaVSName })
        {
            if (pt.second.hasParameter(zname))
            {
                check("particle type " + pt.first.id(), zname, pt.second.parameterConst(zname), perNm);
            }
        }
    }
    if (!errors.empty())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Force field %s has parameters in units that cannot be used:%s",
                                                           filename().c_str(), errors.c_str())));
    }
}

void ForceField::checkVsitePotentials() const
{
    // The construction follows from the interaction type, see VsiteHandler,
    // the parameters from the potential, which must therefore be the same
    std::string errors;
    for (const auto &fs : forces_)
    {
        if (isVsite(fs.first) &&
            interactionTypeToString(fs.first) != potentialToString(fs.second.potential()))
        {
            errors += gmx::formatString("\n  the virtual sites of type %s have function %s, it has to be %s",
                                        interactionTypeToString(fs.first).c_str(),
                                        potentialToString(fs.second.potential()).c_str(),
                                        interactionTypeToString(fs.first).c_str());
        }
    }
    if (!errors.empty())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Force field %s has virtual sites that cannot be constructed:%s",
                                                           filename().c_str(), errors.c_str())));
    }
}

void ForceField::checkConsistency(MsgHandler *msgHandler) const
{
    int  nerror = 0;
    auto cga    = chargeGenerationAlgorithm();
    if (cga == ChargeGenerationAlgorithm::NONE || cga == ChargeGenerationAlgorithm::ESP)
    {
        return;
    }
    if (interactionPresent(InteractionType::BONDCORRECTIONS) &&
        chargeGenerationAlgorithm() != ChargeGenerationAlgorithm::SQE)
    {
        if (msgHandler)
        {
            msgHandler->msg(ACTStatus::Warning, "Can only have bond corrections when ChargeGenerationAlgorithm = SQE\n");
        }
        nerror += 1;
    }
    auto itype = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
    const auto eem = findForcesConst(itype);
    for (const auto &atp : alexandria_)
    {
        if (!atp.second.hasInteractionType(itype))
        {
            auto &qparam = atp.second.parameterConst("charge");
            if (qparam.mutability() == Mutability::ACM)
            {
                if (msgHandler)
                {
                    msgHandler->msg(ACTStatus::Warning,
                                    gmx::formatString("No %s type for particletype %s\n",
                                                      interactionTypeToParticleSubtype(itype).c_str(),
                                                      atp.second.id().id().c_str()));
                }
                nerror += 1;
            }
            else
            {
                continue;
            }
        }
        else
        {
            const auto acmtype = atp.second.interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            // Check whether the ACM parameters are present
            if (!eem.parameterExists(acmtype))
            {
                if (msgHandler)
                {
                    msgHandler->msg(ACTStatus::Error,
                                    gmx::formatString("No %s parameters for %s of particletype %s in %s",
                                                      interactionTypeToString(itype).c_str(),
                                                      acmtype.id().c_str(),
                                                      atp.second.id().id().c_str(),
                                                      filename().c_str()));
                }
                nerror += 1;
            }
        }
    }
    if (chargeDistribution() == ChargeDistributionType::Point)
    {
        // Check for non-zero zeta
        for(const auto &pt : alexandria_)
        {
            for(const char *zname : { zetaName, zetaVSName })
            {
                if (pt.second.hasParameter(zname) && pt.second.parameterConst(zname).internalValue() != 0)
                {
                    if (msgHandler)
                    {
                        msgHandler->msg(ACTStatus::Warning,
                                        gmx::formatString("Force field specifies Point charges but %s = %g for %s",
                                                          zname, pt.second.parameterConst(zname).internalValue(),
                                                          pt.first.id().c_str()));
                    }
                    nerror += 1;
                }
            }
        }
    }
    if (nerror > 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("ForceField inconsistency.").c_str()));
    }
}

//! \brief Implementation of alexandria::ffOption for std::string
template <> bool ffOption(const ForceField  &pd,
                          InteractionType    itype,
                          const std::string &name,
                          std::string       *value)
{
    if (!pd.interactionPresent(itype))
    {
        return false;
    }
    const auto& fs = pd.findForcesConst(itype);
    if (!fs.optionExists(name))
    {
        return false;
    }
    *value = fs.optionValue(name);
    return true;
}

//! \brief Implementation of alexandria::ffOption for int
template <> bool ffOption(const ForceField  &pd,
                          InteractionType    itype,
                          const std::string &name,
                          int               *value)
{
    std::string vstr;
    if (!ffOption(pd, itype, name, &vstr))
    {
        return false;
    }
    *value = my_atoi(vstr, name);
    return true;
}

//! \brief Implementation of alexandria::ffOption for double
template <> bool ffOption(const ForceField  &pd,
                          InteractionType    itype,
                          const std::string &name,
                          double            *value)
{
    std::string vstr;
    if (!ffOption(pd, itype, name, &vstr))
    {
        return false;
    }
    *value = my_atof(vstr, name);
    return true;
}

} // namespace alexandria
