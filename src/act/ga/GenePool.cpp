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
    std::vector<double> values(popSize());
    for (size_t i = 0; i < genomeSize_; i++)
    {
        for (size_t j = 0; j < popSize(); j++)
        {
            values[j] = genomes_[j].base(i);
        }
        double sum = 0;
        for (auto ele : values)
        {
            sum += ele;
        }
        vec[i] = sum / static_cast<double>(popSize());
    }
    return vec;
}

std::vector<double> GenePool::median() const
{
    std::vector<double> vec(genomeSize_);
    std::vector<double> values(popSize());
    for (size_t i = 0; i < genomeSize_; i++)
    {
        for (size_t j = 0; j < popSize(); j++)
        {
            values[j] = genomes_[j].base(i);
        }
        std::sort(values.begin(), values.end());
        // Remember that popsize can be odd in MCMC
        if (popSize() % 2 != 0)
        {
            vec[i] = values[genomeSize_/2];
        }
        else
        {
            vec[i] = (values[genomeSize_/2 - 1]+values[genomeSize_/2])/2.0;
        }
    }
    return vec;
}

} // namespace ga
