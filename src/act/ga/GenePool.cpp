#include "GenePool.h"

#include <algorithm>
#include <limits>
#include <vector>

#include <cstdio>
#include <cmath>

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
        int i = 0;
        for (auto &ind : genomes_)
        {
            auto istr = gmx::formatString(" %3d ", i++);
            ind.print(istr.c_str(), fp);
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

size_t GenePool::findBestIndex(const iMolSelect ims) const
{
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

const Genome &GenePool::getBest(const iMolSelect ims) const
{
    const size_t index = findBestIndex(ims);
    GMX_RELEASE_ASSERT(
        index < popSize() && genomes_[index].hasFitness(ims),
        "Oh no! The returned index is beyond the vector limits or the genome does not have a fitness entry for the dataset..."
    );
    return genomes_[index];
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

std::vector<double> GenePool::min() const
{
    std::vector<double> vec(genomeSize_, std::numeric_limits<double>::max());
    for (size_t i = 0; i < genomeSize_; i++)
    {
        for (size_t j = 0; j < popSize(); j++)
        {
            vec[i] = genomes_[j].base(i) < vec[i] ? genomes_[j].base(i) : vec[i];
        }
    }
    return vec;
}

std::vector<double> GenePool::max() const
{
    std::vector<double> vec(genomeSize_, std::numeric_limits<double>::min());
    for (size_t i = 0; i < genomeSize_; i++)
    {
        for (size_t j = 0; j < popSize(); j++)
        {
            vec[i] = genomes_[j].base(i) > vec[i] ? genomes_[j].base(i) : vec[i];
        }
    }
    return vec;
}

std::vector<double> GenePool::mean() const
{
    std::vector<double> vec(genomeSize_);
    for (size_t i = 0; i < genomeSize_; i++)
    {
        double sum = 0;
        for (size_t j = 0; j < popSize(); j++)
        {
            sum += genomes_[j].base(i);
        }
        vec[i] = sum / static_cast<double>(popSize());
    }
    return vec;
}

std::vector<double> GenePool::stdev(const std::vector<double> &mean) const
{
    std::vector<double> vec(genomeSize_);
    for (size_t i = 0; i < genomeSize_; i++)
    {
        double sum = 0;
        for (size_t j = 0; j < popSize(); j++)
        {
            sum += (genomes_[j].base(i)-mean[i])*(genomes_[j].base(i)-mean[i]);
        }
        vec[i] = sqrt(sum / static_cast<double>(popSize()));
    }
    return vec;
}

std::vector<double> GenePool::median() const
{
    std::vector<double> vec(genomeSize_);
    size_t popsize = popSize();
    std::vector<double> values(popsize);
    for (size_t i = 0; i < genomeSize_; i++)
    {
        for (size_t j = 0; j < popsize; j++)
        {
            values[j] = genomes_[j].base(i);
        }
        std::sort(values.begin(), values.end());
        // Remember that popsize can be odd in MCMC
        if (popsize % 2 != 0)
        {
            vec[i] = values[popsize/2];
        }
        else
        {
            vec[i] = (values[popsize/2 - 1]+values[popsize/2])/2.0;
        }
    }
    return vec;
}

} // namespace ga
