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

#include "poldata.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <map>
#include <vector>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"

#include "gmx_simple_comm.h"
#include "stringutil.h"

namespace alexandria
{

void Poldata::setFilename(const std::string &fn2)
{
    GMX_RELEASE_ASSERT((fn2.size() > 0),
                       "Trying to set empty Poldata filename");
    if (0 != filename_.size())
    {
        fprintf(stderr, "Changing Poldata filename_ from %s to %s\n",
                filename_.c_str(), fn2.c_str());

    }
    filename_ = fn2;
}

bool Poldata::yang() const
{
    // Note that this is not a good way of doing things. Filenames may change.
    return (filename_.find("Yang.dat") != std::string::npos);
}

bool Poldata::rappe() const
{
    // Note that this is not a good way of doing things. Filenames may change.
    return (filename_.find("Rappe.dat") != std::string::npos);
}

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Atom STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

gmx_bool Poldata::strcasestrStart(std::string needle, std::string haystack)
{
    std::string ptr;
    ptr = strcasestr(haystack.c_str(), needle.c_str());
    return (ptr == haystack);
}

const std::string Poldata::ztype2elem(const std::string &ztype) const
{
    size_t i;
    if (ztype.size() != 0)
    {
        for (i = 0; i < alexandria_.size(); i++)
        {
            if (alexandria_[i].interactionTypeToIdentifier(InteractionType::CHARGEDISTRIBUTION).id() == ztype)
            {
                return alexandria_[i].element();
            }
        }
    }
    gmx_fatal(FARGS, "No such zeta type %s", ztype.c_str());
}

bool Poldata::typeToInteractionType(const std::string &type, 
                                    InteractionType   *itype)
{
    if (type2Itype_.empty())
    {
        for(const auto &fs : forcesConst())
        {
            auto iType = fs.first;
            for(const auto &fp : fs.second.parametersConst())
            {
                for(const auto &myfp : fp.second)
                {
                    auto mytype = myfp.first;
                    if (type2Itype_.find(mytype) != type2Itype_.end() &&
                        type2Itype_.find(mytype)->second != iType)
                    {
                        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Parameter type %s used for more than one InteractionType (%s and %s).",
                                                                           mytype.c_str(),
                                                                           interactionTypeToString(iType).c_str(),
                                                                           interactionTypeToString(type2Itype_.find(mytype)->second).c_str()).c_str()));
                    }
                    type2Itype_.insert({mytype, iType});
                }
            }
        }
    }
    if (type2Itype_.find(type) == type2Itype_.end())
    {
        printf("Cannot find type %s in force field file %s", type.c_str(), filename_.c_str());
        return false;
    }
    else
    {
        *itype = type2Itype_.find(type)->second;
        return true;
    }
}

ChargeGenerationAlgorithm Poldata::chargeGenerationAlgorithm() const
{
    if (interactionPresent(InteractionType::ELECTRONEGATIVITYEQUALIZATION))
    {
        if (interactionPresent(InteractionType::BONDCORRECTIONS))
        {
            return ChargeGenerationAlgorithm::SQE;
        }
        else
        {
            return ChargeGenerationAlgorithm::EEM;
        }
    }
    else
    {
        return ChargeGenerationAlgorithm::ESP;
    } 
}

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Polarizability STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

bool Poldata::atypeToPtype(const std::string &atype,
                           std::string       *ptype) const
{
    auto ai = findParticleType(atype);
    if (ai != alexandria_.end())
    {
        ptype->assign(ai->interactionTypeToIdentifier(InteractionType::POLARIZATION).id());
        return true;
    }
    return false;
}

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Virtual Site STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

void  Poldata::addVsite(const std::string &atype,
                        const std::string &type,
                        int                number,
                        double             distance,
                        double             angle,
                        int                ncontrolatoms)
{
    size_t i;
    for (i = 0; i < vsite_.size(); i++)
    {
        if (vsite_[i].atype() == atype &&
            vsite_[i].type()  == string2vsiteType(atype.c_str()))
        {
            break;
        }
    }
    if (i == vsite_.size())
    {
        Vsite vs(atype, type, number, distance, angle, ncontrolatoms);
        vsite_.push_back(vs);
    }
    else
    {
        fprintf(stderr, "vsite type %s was already added to Poldata record\n", atype.c_str());
    }
}

void Poldata::checkForPolarizability()
{
    auto f = forces_.find(InteractionType::POLARIZATION);
    
    polarizable_ = (f != forces_.end() && f->second.numberOfParameters() > 0);
}

void Poldata::addParticleType(const ParticleType &ptp)
{
    if (hasParticleType(ptp.id()))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Trying to add the same particle type %s twice", ptp.id().id().c_str()).c_str()));
    }
    alexandria_.push_back(ptp);
}

void Poldata::addForces(const std::string             &interaction,
                        const ForceFieldParameterList &forces)
{
    auto iType = stringToInteractionType(interaction.c_str());
    auto f     = forces_.find(iType);
    
    GMX_RELEASE_ASSERT(f == forces_.end(),
                       gmx::formatString("Will not add a second ForceFieldParameterList for %s\n",
                                         interaction.c_str()).c_str());
    forces_.insert({iType, forces});
}

bool Poldata::atypeToBtype(const std::string &atype,
                           std::string       *btype) const
{
    auto ai = findParticleType(atype);
    if (ai != alexandria_.end())
    {
        btype->assign(ai->interactionTypeToIdentifier(InteractionType::BONDS).id());
        return true;
    }
    return false;
}

bool Poldata::atypeToZtype(const std::string &atype,
                           std::string       *ztype) const
{
    auto ai = findParticleType(atype);
    if (ai != alexandria_.end())
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

void Poldata::addSymcharges(const std::string &central,
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

const Identifier Poldata::atomtypesToZetaIdentifier(const std::vector<std::string> atoms) const
{
    std::vector<std::string> ztypes;
    for (auto &a : atoms)
    {
        std::string ztype;
        if (!atypeToZtype(a, &ztype))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find zeta type for %s", a.c_str()).c_str()));
        }
        ztypes.push_back(ztype);
    }
    return Identifier(ztypes, CanSwap::No);
}

CommunicationStatus Poldata::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &filename_);
        gmx_send_str(cr, dest, &alexandriaVersion_);
        gmx_send_int(cr, dest, nexcl_);
        gmx_send_double(cr, dest, gtEpsilonR_);
        gmx_send_str(cr, dest, &vsite_angle_unit_);
        gmx_send_str(cr, dest, &vsite_length_unit_);
        gmx_send_int(cr, dest, static_cast<int>(ChargeGenerationAlgorithm_));
        gmx_send_int(cr, dest, polarizable_ ? 1 : 0);
        /* Send Ffatype */
        gmx_send_int(cr, dest, alexandria_.size());
        for (auto &alexandria : alexandria_)
        {
            cs = alexandria.Send(cr, dest);
            if (CS_OK != cs)
            {
                break;
            }
        }

        /* Send Vsite */
        if (CS_OK == cs)
        {
            gmx_send_int(cr, dest, vsite_.size());
            for (auto &vsite : vsite_)
            {
                cs = vsite.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }

        /* Force Field Parameter Lists */
        if (CS_OK == cs)
        {
            gmx_send_int(cr, dest, forces_.size());
            for (auto &force : forces_)
            {
                std::string key(interactionTypeToString(force.first));
                gmx_send_str(cr, dest, &key);
                cs = force.second.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }

        /* Send Symcharges */
        if (CS_OK == cs)
        {
            gmx_send_int(cr, dest, symcharges_.size());
            for (auto &symcharges : symcharges_)
            {
                cs = symcharges.Send(cr, dest);
                if (CS_OK != cs)
                {
                    break;
                }
            }
        }
    }
    return cs;
}

CommunicationStatus Poldata::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &filename_);
        gmx_recv_str(cr, src, &alexandriaVersion_);
        nexcl_                = gmx_recv_int(cr, src);
        gtEpsilonR_           = gmx_recv_double(cr, src);
        gmx_recv_str(cr, src, &vsite_angle_unit_);
        gmx_recv_str(cr, src, &vsite_length_unit_);
        ChargeGenerationAlgorithm_ = static_cast<ChargeGenerationAlgorithm>(gmx_recv_int(cr, src));
        polarizable_          = static_cast<bool>(gmx_recv_int(cr, src));
        /* Rceive Ffatype */
        size_t nalexandria = gmx_recv_int(cr, src);
        alexandria_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nalexandria); n++)
        {
            ParticleType alexandria;
            cs = alexandria.Receive(cr, src);
            if (CS_OK == cs)
            {
                alexandria_.push_back(alexandria);
            }
        }
        if (debug)
        {
            fprintf(debug, "Done receiving atomtypes\n");
            fflush(debug);
        }

        /* Receive Vsites */
        size_t nvsite = gmx_recv_int(cr, src);
        vsite_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nvsite); n++)
        {
            Vsite vsite;
            cs = vsite.Receive(cr, src);
            if (CS_OK == cs)
            {
                vsite_.push_back(vsite);
            }
        }
        if (debug)
        {
            fprintf(debug, "Done receiving vsites\n");
            fflush(debug);
        }

        /* Receive Listed Forces */
        size_t nforces           = gmx_recv_int(cr, src);
        forces_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nforces); n++)
        {
            ForceFieldParameterList fs;
            std::string             key;
            gmx_recv_str(cr, src, &key);
            InteractionType iType = stringToInteractionType(key.c_str());
            cs                    = fs.Receive(cr, src);
            if (CS_OK == cs)
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
        if (CS_OK == cs)
        {
            size_t nsymcharges = gmx_recv_int(cr, src);
            symcharges_.clear();
            for (size_t n = 0; (CS_OK == cs) && (n < nsymcharges); n++)
            {
                Symcharges symcharges;
                cs = symcharges.Receive(cr, src);
                if (CS_OK == cs)
                {
                    symcharges_.push_back(symcharges);
                }
            }
        }
    }
    return cs;
}

void Poldata::broadcast_eemprop(const t_commrec *cr)
{
    const int src = 0;
    /* Force Field Parameter Lists */
    std::vector<InteractionType> eemlist = 
        { InteractionType::BONDCORRECTIONS,
          InteractionType::CHARGEDISTRIBUTION,
          InteractionType::POLARIZATION,
          InteractionType::ELECTRONEGATIVITYEQUALIZATION };
    if (MASTER(cr))
    {
        for (auto dest = 1; dest < cr->nnodes; dest++)
        {
            auto cs = gmx_send_data(cr, dest);
            if (CS_OK == cs)
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Going to update Poldata::eemprop on node %d\n", dest);
                }
                for(auto myeem : eemlist)
                {
                    auto fs = forces_.find(myeem);
                    if (fs != forces_.end())
                    {
                        gmx_send_int(cr, dest, 1);
                        cs = fs->second.Send(cr, dest);
                    }
                    else
                    {
                        gmx_send_int(cr, dest, 0);
                    }
                }
            }
            gmx_send_done(cr, dest);
        }
    }
    else
    {
        auto cs = gmx_recv_data(cr, src);
        if (CS_OK == cs)
        {
            /* Receive EEMprops and Bond Corrections */
            for(auto myeem : eemlist)
            {
                int nbc = gmx_recv_int(cr, src);
                if (nbc == 1)
                {
                    auto fs = forces_.find(myeem);
                    if (fs != forces_.end())
                    {
                        forces_.erase(fs);
                    }
                    ForceFieldParameterList eem;
                    eem.Receive(cr, src);
                    addForces(interactionTypeToString(myeem), eem);
                }
            }
        }
        else
        {
            if (nullptr != debug)
            {
                fprintf(debug, "Could not update eem properties on node %d\n", cr->nodeid);
            }
        }
        gmx_recv_data(cr, src);
    }   
}

void Poldata::broadcast(const t_commrec *cr)
{
    const int src = 0;
    if (MASTER(cr))
    {
        for (int dest = 1; dest < cr->nnodes; dest++)
        {
            auto cs = gmx_send_data(cr, dest);
            if (CS_OK == cs)
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Going to update Poldata on node %d\n", dest);
                }
                Send(cr, dest);
            }
            gmx_send_done(cr, dest);
        }
    }
    else
    {
        auto cs = gmx_recv_data(cr, src);
        if (CS_OK == cs)
        {
            auto cs = Receive(cr, src);
            if (CS_OK == cs)
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Poldata is updated on node %d\n", cr->nodeid);
                }
            }
            else
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Could not update Poldata on node %d\n", cr->nodeid);
                }
            }
        }
        gmx_recv_data(cr, src);
    }
}

void Poldata::checkConsistency(FILE *fp) const
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
        fprintf(fp, "Can only have bond corrections when ChargeGenerationAlgorithm = SQE\n");
        nerror += 1;
    }
    auto eem = findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
    for (const auto &atp : alexandria_)
    {
        auto ztype = atp.interactionTypeToIdentifier(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
        if (ztype.id().empty())
        {
            continue;
        }
        // Check whether zeta types are present
        if (!eem.parameterExists(ztype))
        {
            if (fp)
            {
                fprintf(fp, "ERROR: No eemprops for %s in Poldata::checkConsistency\n", ztype.id().c_str());
            }
            nerror += 1;
        }
        else
        {
            auto eep = eem.findParametersConst(ztype);
            double chi0 = eep["chi"].value();
            double J00  = eep["jaa"].value();
            if (nullptr != fp)
            {
                fprintf(fp, "chi0 %g J00 %g", chi0, J00);
            }
            double zeta = 0;//eep["zeta"].value();
            int    row  = eep["row"].value();
            double q    = eep["charge"].value();
            if (nullptr != fp)
            {
                fprintf(fp, " row %d zeta %g q %g", row, zeta, q);
            }
        }
        if (nullptr != fp)
        {
            fprintf(fp, "\n");
        }
    }
    if (nerror > 0)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Poldata inconsistency. Use the -debug 1 flag to find out more").c_str()));
    }
}

} // namespace alexandria
