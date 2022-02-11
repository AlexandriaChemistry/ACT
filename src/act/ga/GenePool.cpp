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

void GenePool::sort(iMolSelect ims)
{
    std::sort(genomes_.begin(), genomes_.end(), 
              [ims](const Genome &a, const Genome &b) -> bool
              { 
                  if (a.hasFitness(ims))
                  {
                      if (b.hasFitness(ims))
                      {
                          return a.fitness(ims) < b.fitness(ims); 
                      }
                      else
                      {
                          return false;
                      } 
                  }
                  else
                  {
                      return b.hasFitness(ims);
                  }
              });
}

size_t GenePool::findBestIndex() const
{
    auto ims = iMolSelect::Train;
    return std::min_element(genomes_.begin(), genomes_.end(),
                            [ims](const Genome &a, const Genome &b) 
                            {
                                if (a.hasFitness(ims))
                                {
                                    if (b.hasFitness(ims))
                                    {
                                        return a.fitness(ims) < b.fitness(ims); 
                                    }
                                    else
                                    {
                                        return false;
                                    } 
                                }
                                else
                                {
                                    return b.hasFitness(ims);
                                }
                            }
                            ) - genomes_.begin();
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
