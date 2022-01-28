#include "Genome.h"
    
#include <cstdio>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace ga
{
    
void Genome::print(FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "genome_: [ ");
    for (auto &ele : genome_)
    {
        fprintf(fp, "%8g ", ele);
    }
    fprintf(fp, "]; ");
    auto ff = fitness_.find(iMolSelect::Train);
    if (ff != fitness_.end())
    {
        fprintf(fp, "fitness_[Train]: %8g; probability_: %8g\n",
                ff->second, probability_);
    }
    else
    {
        fprintf(fp, "\n");
    }
}

void Genome::print(const char *name, FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "%s", name);
    print(fp);
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
