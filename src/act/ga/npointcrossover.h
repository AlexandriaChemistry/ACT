/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_NPOINTCROSSOVER_H
#define GA_NPOINTCROSSOVER_H


#include "Crossover.h"


namespace ga
{


/*!
 * \brief Performs N-Point Crossover operation
 * See <a href="https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)#One-point_crossover">this</a> for details.
 *
 */
class NPointCrossover : public Crossover
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
     * \param[in] seed              seed for the random number generator
     */
    NPointCrossover(const size_t chromosomeLength,
                    const size_t order,
                    const int seed)
    : ga::Crossover(chromosomeLength, seed), gen(rd()), 
      crossoverPoints_(order + 2)
    {
        gen.seed(seed);
        order_ = order;
        if (chromosomeLength >= 2)
        {
            availableIndices_.resize(chromosomeLength - 1);
            for (size_t i = 1; i < chromosomeLength; i++)
            {
                availableIndices_[i-1] = i;
            }
        }
        crossoverPoints_[0] = 0;
        crossoverPoints_[crossoverPoints_.size() - 1] = chromosomeLength;
    }

    virtual void offspring(ga::Genome *parent1,
                           ga::Genome *parent2,
                           ga::Genome *child1,
                           ga::Genome *child2);

};


} //namespace ga


#endif //ALEXANDRIA_NPOINTCROSSOVER_H
