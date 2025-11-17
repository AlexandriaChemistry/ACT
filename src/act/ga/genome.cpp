/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024
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
#include "genome.h"
    
#include <cstdio>

#include "act/basics/dataset.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace ga
{

std::string Genome::print(const char *name) const
{
    std::string out = gmx::formatString("%s.\nGenome: [ ", name);
    for (auto &ele : genome_)
    {
        out += gmx::formatString("%8g ", ele);
    }
    out += "]; ";
    for (const auto &pair : fitness_)
    {
        out += gmx::formatString("fitness_[%s]: %8g; ", iMolSelectName(pair.first), pair.second);
    }
    out += gmx::formatString("probability_: %8g\n", probability_);
    return out;
}

double Genome::fitness(iMolSelect ims) const
{
    auto ff = fitness_.find(ims);
    if (ff == fitness_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No fitness for data set %s", iMolSelectName(ims)).c_str()));
    }
    return ff->second;
}

void Genome::setFitness(iMolSelect ims, double fitness)
{
    auto ff = fitness_.find(ims);
    if (ff == fitness_.end())
    {
        fitness_.insert({ims, fitness});
    }
    else
    {
        ff->second = fitness;
    }
}

void Genome::unsetFitness(iMolSelect ims)
{
    auto ff = fitness_.find(ims);
    if (ff != fitness_.end())
    {
        fitness_.erase(ff);
    }
}

double Genome::base(size_t index) const
{
    if (index >= genome_.size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Index %zu out of range (should be less than %zu)", index, genome_.size()).c_str()));
    }
    return genome_[index];
}

void Genome::setBase(size_t index, double value)
{
    if (index >= genome_.size())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Index (%zu) out of range (should be < %zu)", index, genome_.size()).c_str()));
    }
    genome_[index] = value;
}

void Genome::Send(const alexandria::CommunicationRecord *cr, int dest) const
{
    cr->send(dest, genome_);
    cr->send(dest, fitness_.size());
    for(auto &f : fitness_)
    {
        std::string imsName(iMolSelectName(f.first));
        cr->send(dest, imsName);
        cr->send(dest, f.second);
    }
    cr->send(dest, probability_);
}
    
void Genome::BroadCast(const alexandria::CommunicationRecord *cr, int root, MPI_Comm comm)
{
    cr->bcast(&genome_, comm);
    size_t nfmap = fitness_.size();
    cr->bcast(&nfmap, comm);
    if (cr->rank() == root)
    {
        for(auto &f : fitness_)
        {
            std::string imsName(iMolSelectName(f.first));
            cr->bcast(&imsName, comm);
            cr->bcast(&f.second, comm);
        }
    }
    else
    {
        fitness_.clear();
        for (size_t i = 0; i < nfmap; i++)
        {
            std::string imsName;
            cr->bcast(&imsName, comm);
            iMolSelect ims;
            if (!name2molselect(imsName, &ims))
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Invalid iMolSelect name %s",
                                                               imsName.c_str()).c_str()));
            }
            double fitness = 0;
            cr->bcast(&fitness, comm);
            fitness_.insert({ims, fitness});
        }
    }
    cr->bcast(&probability_, comm);
}

void Genome::Receive(const alexandria::CommunicationRecord *cr, int src)
{
    cr->recv(src, &genome_);
    size_t nfmap;
    cr->recv(src, &nfmap);
    fitness_.clear();
    for (size_t i = 0; i < nfmap; i++)
    {
        std::string imsName;
        cr->recv(src, &imsName);
        iMolSelect ims;
        if (!name2molselect(imsName, &ims))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Invalid iMolSelect name %s",
                                                           imsName.c_str()).c_str()));
        }
        double fitness;
        cr->recv(src, &fitness);
        fitness_.insert({ims, fitness});
    }
    cr->recv(src, &probability_);
}

} // namespace ga
