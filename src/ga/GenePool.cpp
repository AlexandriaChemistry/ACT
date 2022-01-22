#include "GenePool.h"

#include <cstdio>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "Genome.h"

namespace ga
{

void GenePool::print(FILE *fp) const
{
    if (fp)
    {
        fprintf(fp, "Population:\n");
        for (auto &ind : genomes_)
        {
            ind.print(fp);
        }
    }
}

size_t GenePool::findBestIndex() const
{
    return std::min_element(genomes_.begin(), genomes_.end(),
                            [](const Genome &a, const Genome &b) 
                            { return a.fitness(iMolSelect::Train) < b.fitness(iMolSelect::Train); }) - genomes_.begin();
}

void GenePool::addGenome(const std::vector<double> &genome, double fitness)
{
    if (genome.size() != genomeSize_)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("All genes must be the same length. Expected %zu, got %zu.", genomeSize_, genome.size()).c_str()));
    }
    FitnessMap fm = { { iMolSelect::Train, fitness } };
    Genome g(genome, fm);
    genomes_.push_back(g);
}

void GenePool::addGenome(const Genome &genome)
{
    if (genome.nBase() != genomeSize_)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("All genes must be the same length. Expected %zu, got %zu.", genomeSize_, genome.nBase()).c_str()));
    }
    genomes_.push_back(genome);
}

} // namespace ga
