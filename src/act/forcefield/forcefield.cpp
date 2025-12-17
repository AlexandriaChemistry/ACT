/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
#include <vector>

#include "gromacs/fileio/md5.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/textreader.h"

#include "act/forcefield/act_checksum.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/potential.h"
#include "act/utility/stringutil.h"

namespace alexandria
{

bool ForceField::verifyCheckSum(FILE              *fp,
                                const std::string &checkSum)
{
    bool match = checkSum == checkSum_;
    if (!match && fp)
    {
        fprintf(fp, "Checksum mismatch in %s. Expected %s Found %s\n",
                filename_.c_str(),
                checkSum_.c_str(), checkSum.c_str());
    }
    return match;
}

std::vector<std::string> ForceField::info() const
{
    std::vector<std::string> out;
    out.push_back("Force field information");
    out.push_back("-----------------------------------------------");
    out.push_back(gmx::formatString("Filename:    %s", filename_.c_str()));
    out.push_back(gmx::formatString("CheckSum:    %s", checkSum_.c_str()));
    out.push_back(gmx::formatString("TimeStamp:   %s", timeStamp_.c_str()));
    out.push_back(gmx::formatString("Polarizable: %s", polarizable() ? "True" : "False"));
    out.push_back("Interactions:");
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
    out.push_back("-----------------------------------------------");
    return out;
}

bool ForceField::verifyCheckSum(FILE *fp)
{
    std::string checkSum = forcefieldCheckSum(this);
    return verifyCheckSum(fp, checkSum);
}

void ForceField::updateCheckSum()
{
    checkSum_ = forcefieldCheckSum(this);
}

void ForceField::updateTimeStamp()
{
    auto tnow    = std::chrono::system_clock::now();
    auto ttt     = std::chrono::system_clock::to_time_t(tnow);
    std::tm *now = std::localtime(&ttt);
    timeStamp_   = gmx::formatString("%d-%02d-%02d %02d:%02d:%02d",
                                     now->tm_year+1900, 1+now->tm_mon,
                                     now->tm_mday, now->tm_hour,
                                     now->tm_min, now->tm_sec);
}

void ForceField::setFilename(const std::string &fn2)
{
    GMX_RELEASE_ASSERT((fn2.size() > 0),
                       "Trying to set empty ForceField filename");
    if (0 != filename_.size())
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

gmx_bool ForceField::strcasestrStart(std::string needle, std::string haystack)
{
    std::string ptr;
    ptr = strcasestr(haystack.c_str(), needle.c_str());
    return (ptr == haystack);
}

const std::string ForceField::ztype2elem(const std::string &ztype) const
{
    if (ztype.size() != 0)
    {
        for (auto i : alexandria_)
        {
            if (i.second.interactionTypeToIdentifier(InteractionType::ELECTROSTATICS).id() == ztype)
            {
                return i.second.element();
            }
        }
    }
    gmx_fatal(FARGS, "No such zeta type %s", ztype.c_str());
}

bool ForceField::typeToInteractionType(const std::string &type, 
                                       InteractionType   *itype)
{
    if (type2Itype_.empty())
    {
        type2Itype_.insert({InteractionType::CHARGE, { "charge" } });
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
        auto ittt = type.substr(0, colon);
        InteractionType itp;
        if (stringToInteractionType(ittt, &itp))
        {
            auto iii  = type2Itype_.find(itp);
            if (iii == type2Itype_.end())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such interaction type '%s' in force field", ittt.c_str()).c_str()));
            }
            auto ptype = type.substr(colon+1, type.size());
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
        for (auto tt : type2Itype_)
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
        // TODO: Add check for number of interactions?
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

bool ForceField::atypeToPtype(const std::string &atype,
                           std::string       *ptype) const
{
    auto ai = findParticleType(atype);
    if (ai)
    {
        ptype->assign(ai->interactionTypeToIdentifier(InteractionType::POLARIZATION).id());
        return true;
    }
    return false;
}

void ForceField::checkForPolarizability()
{
    auto f = forces_.find(InteractionType::POLARIZATION);
    
    polarizable_ = (f != forces_.end() && f->second.numberOfParameters() > 0);
}

void ForceField::addParticleType(const ParticleType &ptp)
{
    if (hasParticleType(ptp.id()))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Trying to add the same particle type %s twice", ptp.id().id().c_str()).c_str()));
    }
    alexandria_.insert({ptp.id(), ptp});
}

void ForceField::addForces(InteractionType                iType,
                           const ForceFieldParameterList &forces)
{
    auto f     = forces_.find(iType);
    
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

bool ForceField::atypeToZtype(const std::string &atype,
                           std::string       *ztype) const
{
    auto ai = findParticleType(atype);
    if (ai)
    {
        ztype->assign(ai->interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION).id());
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
        Symcharges symcharges(central, attached, numattach);
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
                std::string key(interactionTypeToString(force.first));
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
                                       MPI_Comm                   comm)
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
                    if (debug)
                    {
                        fprintf(debug, "Done Listed force %s\n", key.c_str());
                        fflush(debug);
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

CommunicationStatus ForceField::Receive(const CommunicationRecord *cr, int src)
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
        /* Rceive Ffatype */
        size_t nalexandria;
        cr->recv(src, &nalexandria);
        alexandria_.clear();
        for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < nalexandria); n++)
        {
            ParticleType alexandria;
            cs = alexandria.Receive(cr, src);
            if (CommunicationStatus::OK == cs)
            {
                alexandria_.insert({alexandria.id(), alexandria});
            }
        }
        if (debug)
        {
            fprintf(debug, "Done receiving atomtypes\n");
            fflush(debug);
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
            if (debug)
            {
                fprintf(debug, "Done Listed force %s\n", key.c_str());
                fflush(debug);
            }
        }
        if (debug)
        {
            fprintf(debug, "Done Listed forces\n");
            fflush(debug);
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

void ForceField::sendParticles(const CommunicationRecord *cr, int dest)
{
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Going to update ForceField::particles on node %d\n", dest);
        }
        for(auto &ax : alexandria_)
        {
            for(auto &p : ax.second.parametersConst())
            {
                auto mut = p.second.mutability();
                if (Mutability::Free    == mut ||
                    Mutability::Bounded == mut)
                {
                    cr->send(dest, 1);
                    cr->send(dest, ax.second.id().id());
                    cr->send(dest, p.first);
                    cr->send(dest, p.second.value());
                }
            }
        }
        cr->send(dest, 0);
    }
    cr->send_done(dest);
}


void ForceField::receiveParticles(const CommunicationRecord *cr, int src)
{
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        /* Receive Particle info */
        int status;
        cr->recv(src, &status);
        while (1 == status)
        {
            std::string axid, paramname;
            double value;
            cr->recv(src, &axid);
            cr->recv(src, &paramname);
            cr->recv(src, &value);
            findParticleType(axid)->parameter(paramname)->setValue(value);
            cr->recv(src, &status);
        }
    }
    else
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Could not update eem properties on node %d\n", cr->rank());
        }
    }
    GMX_RELEASE_ASSERT(CommunicationStatus::DONE == cr->recv_data(src),
                       "Communication did not end correctly");
}

//! Force Field Parameter Lists
static std::vector<InteractionType> eemlist = 
    { InteractionType::BONDCORRECTIONS,
      InteractionType::ELECTROSTATICS,
      InteractionType::POLARIZATION,
      InteractionType::ELECTRONEGATIVITYEQUALIZATION
    };

void ForceField::sendEemprops(const CommunicationRecord *cr, int dest)
{
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Going to update ForceField::eemprop on node %d\n", dest);
        }
        for(auto myeem : eemlist)
        {
            auto fs = forces_.find(myeem);
            if (fs != forces_.end())
            {
                cr->send(dest, 1);
                // TODO do not ignore return value
                (void) fs->second.Send(cr, dest);
            }
            else
            {
                cr->send(dest, 0);
            }
        }
    }
    cr->send_done(dest);
}

void ForceField::receiveEemprops(const CommunicationRecord *cr, int src)
{
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        /* Receive EEMprops and Bond Corrections */
        for(auto myeem : eemlist)
        {
            int nbc;
            cr->recv(src, &nbc);
            if (nbc == 1)
            {
                auto fs = forces_.find(myeem);
                if (fs != forces_.end())
                {
                    forces_.erase(fs);
                }
                ForceFieldParameterList eem;
                eem.Receive(cr, src);
                addForces(myeem, eem);
            }
        }
    }
    else
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Could not update eem properties on node %d\n", cr->rank());
        }
    }
    GMX_RELEASE_ASSERT(CommunicationStatus::DONE == cr->recv_data(src),
                       "Communication did not end correctly");
}

void ForceField::sendToHelpers(const CommunicationRecord *cr, int root, bool bcast)
{
    if (bcast && debug)
    {
        fprintf(debug, "Will send force field from node %d to helper", cr->rank());
        for(const auto &i: cr->helpers())
        {
            fprintf(debug, " %d", i);
        }
        fprintf(debug, "\n");
        fflush(debug);
    }
    if (cr->rank() == root)
    {
        if (bcast)
        {
            BroadCast(cr, root, cr->send_helpers());
        }
        else
        {
            for (auto &dest : cr->helpers())
            {
                if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
                {
                    if (nullptr != debug)
                    {
                        fprintf(debug, "Going to update ForceField on node %d\n", dest);
                    }
                    Send(cr, dest);
                }
                cr->send_done(dest);
            }
        }
    }
    else
    {
        if (bcast)
        {
            BroadCast(cr, root, cr->send_helpers());
        }
        else
        {
            int src = cr->superior();
            if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
            {
                auto cs = Receive(cr, src);
                if (CommunicationStatus::OK == cs)
                {
                    if (nullptr != debug)
                    {
                        fprintf(debug, "ForceField is updated on node %d\n", cr->rank());
                    }
                }
                else
                {
                    if (nullptr != debug)
                    {
                        fprintf(debug, "Could not update ForceField on node %d\n", cr->rank());
                    }
                }
            }
            GMX_RELEASE_ASSERT(CommunicationStatus::DONE == cr->recv_data(src),
                               "Communication did not end correctly");
        }
    }
}

void ForceField::checkConsistency(FILE *fp) const
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
        fprintf(stderr, "Can only have bond corrections when ChargeGenerationAlgorithm = SQE\n");
        nerror += 1;
    }
    auto itype = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
    auto eem = findForcesConst(itype);
    for (const auto &atp : alexandria_)
    {
        if (!atp.second.hasInteractionType(itype))
        {
            auto &qparam = atp.second.parameterConst("charge");
            if (qparam.mutability() == Mutability::ACM)
            {
                fprintf(stderr, "No %s type for particletype %s\n",
                        interactionTypeToParticleSubtype(itype).c_str(),
                        atp.second.id().id().c_str());
                nerror += 1;
            }
            else
            {
                continue;
            }
        }
        else
        {
            auto acmtype = atp.second.interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
            // Check whether zeta types are present
            if (!eem.parameterExists(acmtype))
            {
                fprintf(stderr, "ERROR: No eemprops for %s in ForceField::checkConsistency\n", acmtype.id().c_str());
                nerror += 1;
            }
            else
            {
                auto eep = eem.findParameterMapConst(acmtype);
                double chi0 = eep["chi"].value();
                double J00  = eep["eta"].value();
                if (nullptr != fp)
                {
                    fprintf(fp, "chi0 %g eta %g", chi0, J00);
                }
                double zeta = 0;
                int    row  = eep["row"].value();
                double q    = eep["charge"].value();
                if (nullptr != fp)
                {
                    fprintf(fp, " row %d zeta %g q %g", row, zeta, q);
                }
            }
        }
        if (nullptr != fp)
        {
            fprintf(fp, "\n");
        }
    }
    const auto itq = InteractionType::ELECTROSTATICS;
    if (interactionPresent(itq))
    {
        auto fs  = findForcesConst(itq);
        if (fs.potential() == Potential::COULOMB_POINT)
        {
            // Check for non-zero zeta
            std::string zeta("zeta");
            for(const auto &myparam : fs.parametersConst())
            {
                for(const auto &param : myparam.second)
                {
                    auto myval = param.second.internalValue();
                    if (param.first == zeta && myval != 0)
                    {
                        fprintf(stderr, "Force field specifies Point charges but %s = %g for %s", zeta.c_str(), myval, myparam.first.id().c_str());
                        nerror += 1;
                    }
                }
            }
        }
    }
    if (nerror > 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("ForceField inconsistency.").c_str()));
    }
}

void ForceField::calcDependent()
{
    auto btype = InteractionType::BONDS;
    std::vector<InteractionType> atypes = { InteractionType::LINEAR_ANGLES,
        InteractionType::ANGLES };
    GMX_RELEASE_ASSERT(interactionPresent(btype), "No bond information present");
    auto ffpbonds = findForcesConst(btype);
    for(auto &atype : atypes)
    {
        if (!interactionPresent(atype))
        {
            continue;
        }
        auto ffpl = findForces(atype)->parameters();
        for (auto &fp : *ffpl)
        {
            for (auto &param : fp.second)
            {                               
                if (param.second.mutability() == Mutability::Dependent &&
                    param.first == "r13")
                {
                    const std::vector<std::string> &atoms = fp.first.atoms();
                    const std::vector<double>       bo    = fp.first.bondOrders();
                    Identifier b1({atoms[0], atoms[1]}, { bo[0] }, CanSwap::Yes);
                    Identifier b2({atoms[1], atoms[2]}, { bo[1] }, CanSwap::Yes);
                        
                }
            }
        }
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
    auto fs = pd.findForcesConst(itype);
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
    *value = my_atoi(vstr.c_str(), name.c_str());
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
    *value = my_atof(vstr.c_str(), name.c_str());
    return true;
}

} // namespace alexandria
