#include "Genome.h"
    
#include <cstdio>

#include "alexandria/gmx_simple_comm.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/exceptions.h"

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
    GMX_RELEASE_ASSERT(index < genome_.size(), "Index out of range");
    return genome_[index];
}

void Genome::setBase(size_t index, double value)
{
    GMX_RELEASE_ASSERT(index < genome_.size(), "Index out of range");
    genome_[index] = value;
}
    
void Genome::Send(const t_commrec *cr, int dest) const
{
    gmx_send_double_vector(cr, dest, &genome_);
    gmx_send_int(cr, dest, fitness_.size());
    for(auto &f : fitness_)
    {
        std::string imsName(iMolSelectName(f.first));
        gmx_send_str(cr, dest, &imsName);
        gmx_send_double(cr, dest, f.second);
    }
    gmx_send_double(cr, dest, probability_);
}
    
void Genome::Receive(const t_commrec *cr, int src)
{
    gmx_recv_double_vector(cr, src, &genome_);
    int nfmap = gmx_recv_int(cr, src);
    fitness_.clear();
    for (int i = 0; i < nfmap; i++)
    {
        std::string imsName;
        gmx_recv_str(cr, src, &imsName);
        iMolSelect ims;
        if (!name2molselect(imsName, &ims))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Invalid iMolSelect name %s",
                                                           imsName.c_str()).c_str()));
        }
        double fitness = gmx_recv_double(cr, src);
        fitness_.insert({ims, fitness});
    }
    probability_ = gmx_recv_double(cr, src);
}

} // namespace ga
