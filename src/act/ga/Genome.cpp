/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
#include "Genome.h"
    
#include <cstdio>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "act/utility/communicationrecord.h"

namespace ga
{
    
void Genome::print(const char *name, FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "genome_%s: [ ", name);
    for (auto &ele : genome_)
    {
        fprintf(fp, "%8g ", ele);
    }
    fprintf(fp, "]; ");
    for (const auto &pair : fitness_)
    {
        fprintf(fp, "fitness_[%s]: %8g; ", iMolSelectName(pair.first), pair.second);
    }
    fprintf(fp, "probability_: %8g\n", probability_);
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
    GMX_RELEASE_ASSERT(index < genome_.size(), "Index out of range");
    genome_[index] = value;
}
    
void Genome::Send(const alexandria::CommunicationRecord *cr, int dest) const
{
    cr->send_double_vector(dest, &genome_);
    cr->send_int(dest, fitness_.size());
    for(auto &f : fitness_)
    {
        std::string imsName(iMolSelectName(f.first));
        cr->send_str(dest, &imsName);
        cr->send_double(dest, f.second);
    }
    cr->send_double(dest, probability_);
}
    
void Genome::BroadCast(const alexandria::CommunicationRecord *cr, int root, MPI_Comm comm)
{
    cr->bcast(&genome_, comm);
    int nfmap = fitness_.size();
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
        for (int i = 0; i < nfmap; i++)
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
    cr->recv_double_vector(src, &genome_);
    int nfmap = cr->recv_int(src);
    fitness_.clear();
    for (int i = 0; i < nfmap; i++)
    {
        std::string imsName;
        cr->recv_str(src, &imsName);
        iMolSelect ims;
        if (!name2molselect(imsName, &ims))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Invalid iMolSelect name %s",
                                                           imsName.c_str()).c_str()));
        }
        double fitness = cr->recv_double(src);
        fitness_.insert({ims, fitness});
    }
    probability_ = cr->recv_double(src);
}

} // namespace ga
