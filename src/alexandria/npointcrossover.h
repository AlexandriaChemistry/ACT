#ifndef ALEXANDRIA_NPOINTCROSSOVER_H
#define ALEXANDRIA_NPOINTCROSSOVER_H


#include "ga/Crossover.h"


namespace alexandria
{


class NPointCrossover : public ga::Crossover
{

private:

    // Random number stuff
    std::random_device  rd;
    std::mt19937        gen;

    //! Amount of crossover cutting points
    size_t order_;
    //! List of all existing indices
    std::vector<size_t> availableIndices_;
    //! List of all crossover points. First element is 0, and last element is \p chromosomeLength_.
    std::vector<size_t> crossoverPoints_;

public:

    /*!
     * Property constructor
     * \param[in] chromosomeLength  the length of the chromosomes
     * \param[in] order             order of the crossover operator (amount of cutting points)
     */
    NPointCrossover(const size_t chromosomeLength,
                    const size_t order)
    : ga::Crossover(chromosomeLength), gen(rd()), order_(order), availableIndices_(chromosomeLength - 1),
      crossoverPoints_(order_ + 2)
    {
        for (size_t i = 1; i < chromosomeLength; i++) availableIndices_[i-1] = i;
        crossoverPoints_[0] = 0;
        crossoverPoints_[crossoverPoints_.size() - 1] = chromosomeLength;
    }

    virtual void offspring(ga::Individual  *parent1,
                           ga::Individual  *parent2,
                           ga::Individual  *child1,
                           ga::Individual  *child2);

};


} //namespace alexandria


#endif //ALEXANDRIA_NPOINTCROSSOVER_H