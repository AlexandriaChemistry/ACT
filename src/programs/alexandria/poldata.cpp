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

void Poldata::addAtype(const std::string &elem,
                       const std::string &desc,
                       const std::string &atype,
                       const std::string &ptype,
                       const std::string &btype,
                       const std::string &ztype,
                       const std::string &refEnthalpy)
{

    size_t i;
    for (i = 0; i < alexandria_.size(); i++)
    {
        if (alexandria_[i].getType().compare(atype) == 0)
        {
            break;
        }
    }
    if (i == alexandria_.size())
    {
        Ffatype sp(desc, atype, ptype, btype, ztype,
                   elem, refEnthalpy);

        alexandria_.push_back(sp);
    }
    else
    {
        fprintf(stderr, "Atom type %s was already added to Poldata record\n", atype.c_str());
    }
}

gmx_bool Poldata::strcasestrStart(std::string needle, std::string haystack)
{
    std::string ptr;
    ptr = strcasestr(haystack.c_str(), needle.c_str());
    return (ptr == haystack);
}

bool Poldata::getAtypeRefEnthalpy(const std::string &atype,
                                  double            *Href) const
{
    auto fa = findAtype(atype);
    if (alexandria_.end() != fa)
    {
        *Href = my_atof(fa->getRefEnthalpy().c_str(), "Href");
        return true;
    }
    return false;
}

const std::string &Poldata::getDesc(const std::string &atype) const
{
    size_t i;
    if (atype.size() != 0)
    {
        for (i = 0; i < alexandria_.size(); i++)
        {
            if (alexandria_[i].getType().compare(atype) == 0)
            {
                return alexandria_[i].getDesc();
            }
        }
    }
    gmx_fatal(FARGS, "No such atomtype %s", atype.c_str());
}

const std::string &Poldata::getElem(const std::string &atype) const
{
    size_t i;
    if (atype.size() != 0)
    {
        for (i = 0; i < alexandria_.size(); i++)
        {
            if (alexandria_[i].getType().compare(atype) == 0)
            {
                return alexandria_[i].getElem();
            }
        }
    }
    gmx_fatal(FARGS, "No such atomtype %s", atype.c_str());
}

const std::string &Poldata::ztype2elem(const std::string &ztype) const
{
    size_t i;
    if (ztype.size() != 0)
    {
        for (i = 0; i < alexandria_.size(); i++)
        {
            if (alexandria_[i].id(eitELECTRONEGATIVITYEQUALIZATION).id() ==
                ztype)
            {
                return alexandria_[i].getElem();
            }
        }
    }
    gmx_fatal(FARGS, "No such zeta type %s", ztype.c_str());
}
#ifdef OLDER
const std::string &Poldata::ztype2atype(const std::string &ztype) const
{
    size_t i;
    if (ztype.size() != 0)
    {
        for (i = 0; i < alexandria_.size(); i++)
        {
            if (alexandria_[i].getZtype().compare(ztype) == 0)
            {
                return alexandria_[i].getType();
            }
        }
    }
    gmx_fatal(FARGS, "No such zeta type %s", ztype.c_str());
}

std::vector<std::string> Poldata::ztype_names() const
{
    std::vector<std::string> ztype_names;   
    for (auto atpi = alexandria_.begin(); atpi != alexandria_.end(); atpi++)
    { 
        if (std::find(ztype_names.begin(), ztype_names.end(), atpi->getZtype()) == ztype_names.end()) 
        {
            ztype_names.push_back(atpi->getZtype());
        }
    }
    return ztype_names;
}
#endif

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Polarizability STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

bool Poldata::atypeToPtype(const std::string &atype,
                           std::string       *ptype) const
{
    auto ai = findAtype(atype);
    if (ai != alexandria_.end())
    {
        ptype->assign(ai->id(eitPOLARIZATION).id());
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
    auto f = forces_.find(eitPOLARIZATION);
    
    polarizable_ = (f != forces_.end() && f->second.numberOfParameters() > 0);
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
    auto ai = findAtype(atype);
    if (ai != alexandria_.end())
    {
        btype->assign(ai->id(eitBONDS).id());
        return true;
    }
    return false;
}

bool Poldata::atypeToZtype(const std::string &atype,
                           std::string       *ztype) const
{
    auto ai = findAtype(atype);
    if (ai != alexandria_.end())
    {
        ztype->assign(ai->id(eitELECTRONEGATIVITYEQUALIZATION).id());
        return true;
    }
    return false;
}

void Poldata::addBtype(const std::string &btype)
{
    size_t i;
    for (i = 0; i < btype_.size(); i++)
    {
        if (btype.compare(btype_[i]) == 0)
        {
            break;
        }
    }
    if (i == btype_.size())
    {
        btype_.push_back(btype);
    }
}


#ifdef OLDER
bool Poldata::searchForce(std::vector<std::string> &atoms,
                          std::string              &params,
                          double                   *refValue,
                          double                   *sigma,
                          size_t                   *ntrain) const
{
    for (auto &f : forces_)
    {
        if (f.searchForce(atoms, params, refValue,
                          sigma, ntrain, false))
        {
            return true;
        }
    }
    return false;
}

bool Poldata::searchForceIType(std::vector<std::string> &atoms,
                               std::string              &params,
                               double                   *refValue,
                               double                   *sigma,
                               size_t                   *ntrain,
                               InteractionType           iType) const
{
    auto f  = findForces(iType);
    if (forcesEnd() != f)
    {
        if (f->searchForce(atoms, params, refValue,
                           sigma, ntrain, iType == eitIMPROPER_DIHEDRALS))
        {
            return true;
        }
    }
    return false;
}

bool Poldata::searchForceBondOrder(std::vector<std::string> &atoms,
                                   std::string              &params,
                                   double                   *refValue,
                                   double                   *sigma,
                                   size_t                   *ntrain,
                                   size_t                    bondOrder) const
{
    for (auto &f : forces_)
    {
        if (f.searchForceBondOrder(atoms, params, refValue,
                                   sigma, ntrain, bondOrder))
        {
            return true;
        }
    }
    return false;
}

bool Poldata::searchForceBondOrderIType(std::vector<std::string> &atoms,
                                        std::string              &params,
                                        double                   *refValue,
                                        double                   *sigma,
                                        size_t                   *ntrain,
                                        size_t                    bondOrder,
                                        InteractionType           iType) const
{
    auto f  = findForces(iType);
    if (forcesEnd() != f)
    {
        if (f->searchForceBondOrder(atoms, params, refValue,
                                    sigma, ntrain, bondOrder))
        {
            return true;
        }
    }
    return false;
}

#endif

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

int Poldata::havePolSupport(const std::string &atype) const
{
    size_t i;
    for (i = 0; i < alexandria_.size(); i++)
    {
        if (atype.compare(alexandria_[i].getType()) == 0)
        {
            return 1;
        }
    }
    return 0;
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
#ifdef OLDER
bool Poldata::haveEemSupport(const std::string &atype,
                             gmx_bool           bAllowZeroParameters) const
{
    auto eep = atypeToEempropsConst(atype);
    return (eep != nullptr &&
            (bAllowZeroParameters || ((eep->getJ0() > 0) && (eep->getChi0() > 0))));
}
#endif
CommunicationStatus Poldata::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &filename_);
        gmx_send_int(cr, dest, nexcl_);
        gmx_send_double(cr, dest, gtEpsilonR_);
        gmx_send_str(cr, dest, &vsite_angle_unit_);
        gmx_send_str(cr, dest, &vsite_length_unit_);
        gmx_send_int(cr, dest, static_cast<int>(ChargeType_));
        gmx_send_int(cr, dest, static_cast<int>(ChargeGenerationAlgorithm_));

        /* Send Ffatype */
        gmx_send_int(cr, dest, alexandria_.size());
        for (auto &alexandria : alexandria_)
        {
            cs = alexandria.Send(cr, dest);
        }

        /* Send Vsite */
        gmx_send_int(cr, dest, vsite_.size());
        for (auto &vsite : vsite_)
        {
            cs = vsite.Send(cr, dest);
        }

        /* Send btype */
        gmx_send_int(cr, dest, btype_.size());
        for (auto &btype : btype_)
        {
            gmx_send_str(cr, dest, &btype);
        }

        /* Force Field Parameter Lists */
        gmx_send_int(cr, dest, forces_.size());
        for (auto &force : forces_)
        {
            std::string key(interactionTypeToString(force.first));
            gmx_send_str(cr, dest, &key);
            cs = force.second.Send(cr, dest);
        }

        /* Send Symcharges */
        gmx_send_int(cr, dest, symcharges_.size());
        for (auto &symcharges : symcharges_)
        {
            cs = symcharges.Send(cr, dest);
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
        nexcl_                = gmx_recv_int(cr, src);
        gtEpsilonR_           = gmx_recv_double(cr, src);
        gmx_recv_str(cr, src, &vsite_angle_unit_);
        gmx_recv_str(cr, src, &vsite_length_unit_);
        ChargeType_                = static_cast<ChargeType>(gmx_recv_int(cr, src));
        ChargeGenerationAlgorithm_ = static_cast<ChargeGenerationAlgorithm>(gmx_recv_int(cr, src));

        /* Rceive Ffatype */
        size_t nalexandria = gmx_recv_int(cr, src);
        alexandria_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nalexandria); n++)
        {
            Ffatype alexandria;
            cs = alexandria.Receive(cr, src);
            if (CS_OK == cs)
            {
                alexandria_.push_back(alexandria);
            }
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

        /* Receive btype */
        size_t nbtype = gmx_recv_int(cr, src);
        btype_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nbtype); n++)
        {
            std::string btype;
            gmx_recv_str(cr, src, &btype);
            if (!btype.empty())
            {
                btype_.push_back(btype);
            }
        }

        /* Receive Listed Forces */
        size_t nforces           = gmx_recv_int(cr, src);
        forces_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nforces); n++)
        {
            ForceFieldParameterList fs;
            std::string             key;
            gmx_recv_str(cr, src, &key);
            InteractionType iType = static_cast<InteractionType>(stringToInteractionType(key.c_str()));
            cs                    = fs.Receive(cr, src);
            if (CS_OK == cs)
            {
                forces_.insert({iType, fs});
            }
        }

        /* Receive Symcharges */
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
    return cs;
}

void Poldata::broadcast_eemprop(const t_commrec *cr)
{
    const int src = 0;
    /* Force Field Parameter Lists */
    std::vector<InteractionType> eemlist = 
        { eitBONDCORRECTIONS,
          eitELECTRONEGATIVITYEQUALIZATION };
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
    if (cga == eqgNONE || cga == eqgESP)
    {
        return;
    }
    auto eem = findForcesConst(eitELECTRONEGATIVITYEQUALIZATION);
    for (auto atp = getAtypeBegin(); atp < getAtypeEnd(); ++atp)
    {
        auto atype = atp->id(eitELECTRONEGATIVITYEQUALIZATION);
        // Check whether zeta types are present
        if (!eem.parameterExists(atype))
        {
            if (fp)
            {
                fprintf(fp, "ERROR: No eemprops for %s in Poldata::checkConsistency\n", atype.id().c_str());
            }
            nerror += 1;
        }
        else
        {
            auto eep = eem.findParametersConst(atype);
            double chi0 = eep["chi"].value();
            double J00  = eep["jaa"].value();
            if (nullptr != fp)
            {
                fprintf(fp, "chi0 %g J00 %g", chi0, J00);
            }
            double zeta = eep["zeta"].value();
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
