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

#include "tune_fc_utils.h"

#include "gromacs/topology/topology.h"
#include "gromacs/utility/strconvert.h"

#include "forcefieldparameter.h"
#include "gmx_simple_comm.h"

namespace alexandria
{

CommunicationStatus ParameterNames::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, ncopies_);
        gmx_send_int(cr, dest, ftype_);
        gmx_send_int(cr, dest, poldataIndex_);
        gmx_send_int(cr, dest, params_.size());
        for (auto &p : params_)
        {
            gmx_send_double(cr, dest, p);
        }
    }
    return cs;
}

CommunicationStatus ParameterNames::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        ncopies_      = gmx_recv_int(cr, src);
        ftype_        = gmx_recv_int(cr, src);
        poldataIndex_ = gmx_recv_int(cr, src);
        int np        = gmx_recv_int(cr, src);
        params_.resize(np, 0.0);
        for (int n = 0; n < np; n++)
        {
            params_[n] = gmx_recv_double(cr, src);
        }
    }
    return cs;
}

CommunicationStatus ForceConstants::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    std::string         itype;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, ftype_);
        itype.assign(interactionTypeToString(itype_));
        gmx_send_str(cr, dest, &itype);
        gmx_send_int(cr, dest, bOpt_);
        gmx_send_int(cr, dest, bn_.size());

        for (auto &bn : bn_)
        {
            bn.first.Send(cr, dest);
            bn.second.Send(cr, dest);
        }
        gmx_send_int(cr, dest, reverseIndex_.size());
        for (auto &ri : reverseIndex_)
        {
            gmx_send_int(cr, dest, ri);
        }
        //gmx_send_int(cr, dest, params_.size());
        //for (auto &param : params_)
        //{
        //   gmx_send_double(cr, dest, param);
        //}
    }
    return cs;
}

CommunicationStatus ForceConstants::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    std::string         itype;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        ftype_        = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &itype);
        itype_        = stringToInteractionType(itype.c_str());
        bOpt_         = gmx_recv_int(cr, src);
        int nbn       = gmx_recv_int(cr, src);

        for (int n = 0; (CS_OK == cs) && (n < nbn); n++)
        {
            Identifier     identifier;
            ParameterNames bn;
            cs = identifier.Receive(cr, src);
            if (CS_OK == cs)
            {
                cs = bn.Receive(cr, src);
            }
            if (CS_OK == cs)
            {
                bn_.insert({identifier, bn});
            }
        }
        int nri       = gmx_recv_int(cr, src);
        for (int n = 0; (CS_OK == cs) && (n < nri); n++)
        {
            auto ri = gmx_recv_int(cr, src);
            reverseIndex_.push_back(ri);
        }
        //int nparam    = gmx_recv_int(cr, src);
        //for (int n = 0; (CS_OK == cs) && (n < nparam); n++)
        //{
        //   auto param = gmx_recv_double(cr, src);
        //   params_.push_back(param);
        //}
    }
    return cs;
}

const ParameterNames &ForceConstants::bondNamesConst(const Identifier &identifier) const
{
    auto p = bn_.find(identifier);
    if (p == bn_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find identifier  %s", identifier.id().c_str()).c_str()));
    }
    return p->second;
}

void ForceConstants::analyzeIdef(const std::vector<MyMol> &mm,
                                 const Poldata            *pd)
{
    if (!bOpt_)
    {
        return;
    }
    for (auto &mymol : mm)
    {
        GMX_THROW(gmx::InternalError("Use UpdateIdef"));
        bool bondsFound = true;
        for (int i = 0; (i < mymol.ltop_->idef.il[ftype_].nr) && bondsFound;
             i += interaction_function[ftype_].nratoms+1)
        {
            std::vector<std::string> atoms;
            auto myatoms = mymol.atomsConst();
            // Loop starts from 1 because the first value is the function type
            for(int j = 1; j <= interaction_function[ftype_].nratoms && bondsFound; j++)
            {
                std::string aa;
                int         ai = mymol.ltop_->idef.il[ftype_].iatoms[i+j];
                if (!pd->atypeToBtype(*myatoms.atomtype[ai], &aa))
                {
                    bondsFound = false;
                }
                else
                {
                    atoms.push_back(aa);
                }
            }
            auto fs = pd->findForcesConst(itype_);
            // TODO: Take bond order into account.
            Identifier bondId(atoms, CanSwap::Yes);
            if (bondsFound && fs.parameterExists(bondId))
            {
                std::vector<double> params;
                for(const auto &f : fs.findParametersConst(bondId))
                {
                    params.push_back(f.second.value());
                }
                auto c = bn_.find(bondId);
                if (c != bn_.end())
                {
                    c->second.inc();
                }
                else
                {
                    ParameterNames bn(1, ftype_, params, 0);
                    addForceConstant(bondId, std::move(bn));
                }
            }
        }
    }
}

void ForceConstants::makeReverseIndex()
{
    int N = 0;
    for (const auto &i : bn_)
    {
        N = std::max(N, i.second.poldataIndex());
    }
    reverseIndex_.resize(N+1, -1);
    int j = 0;
    for (const auto &i : bn_)
    {
        reverseIndex_[i.second.poldataIndex()] = j++;
    }
}

void ForceConstants::dump(FILE *fp) const
{
    if (bOpt_)
    {
        int         ntot = 0;
        fprintf(fp, "Interaction  Bondtypes             Copies Poldata entry\n");
        const char *name = interaction_function[ftype_].name;
        for (const auto &i : bn_)
        {
            fprintf(fp, "%-10s  %-20s  %5d  %5d\n",
                    name, i.first.id().c_str(),
                    i.second.nCopies(), i.second.poldataIndex());
            ntot += i.second.nCopies();
        }
        fprintf(fp, "%d out of %d %s types will be optimized.\n",
                static_cast<int>(bn_.size()), ntot, name);
    }
}

void PoldataUpdate::execute(Poldata *pd)
{
    auto fs = pd->findForces(iType_);
    GMX_RELEASE_ASSERT(fs->parameterExists(identifier_),
                       gmx::formatString("No parameter for %s in %s",
                                         identifier_.id().c_str(),
                                         interactionTypeToString(iType_).c_str()).c_str());
    size_t i = 0;
    auto ptr = fs->findParameters(identifier_);
    for(auto f : *ptr)
    {
        GMX_RELEASE_ASSERT(f.second.mutability() != Mutability::Fixed, 
                           "Fixed listed force parameters should not be here");
        f.second.setValue(parameterValues_[i]);
        i++;
    }
}

void PoldataUpdate::dump(FILE *fp) const
{
    if (nullptr != fp)
    {
        fprintf(fp, "iType: %s identifier: %s parameterValues:",
                interactionTypeToString(iType_).c_str(), identifier_.id().c_str());
        for(const auto &p : parameterValues_)
        {
            fprintf(fp," %g", p);
        }
        fprintf(fp, "\n");
    }
}

CommunicationStatus PoldataUpdate::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        dump(debug);
        gmx_send_int(cr, dest, static_cast<int>(iType_));
        cs = identifier_.Send(cr, dest);
        if (CS_OK == cs)
        {
            gmx_send_int(cr, dest, parameterValues_.size());
            for(const auto &p : parameterValues_)
            {
                gmx_send_double(cr, dest, p);
            }
        }
    }
    return cs;
}

CommunicationStatus PoldataUpdate::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        iType_ = static_cast<InteractionType>(gmx_recv_int(cr, src));
        cs = identifier_.Receive(cr, src);
        if (CS_OK == cs)
        {
            int nparam = gmx_recv_int(cr, src);
            parameterValues_.clear();
            for(int i = 0; i < nparam; i++)
            {
                double p = gmx_recv_double(cr, src);
                parameterValues_.push_back(p);
            }
        }
        dump(debug);
    }
    return cs;
}

} // namespace alexandria
