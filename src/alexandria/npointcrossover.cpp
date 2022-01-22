/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#include "npointcrossover.h"

#include "acmindividual.h"


namespace alexandria
{


void NPointCrossover::offspring(ga::Genome  *parent1,
                                ga::Genome  *parent2,
                                ga::Genome  *child1,
                                ga::Genome  *child2)
{

    // Casting individuals
    //ACMGenome *tmpParent1   = static_cast<ACMGenome*>(parent1);
    //ACMGenome *tmpParent2   = static_cast<ACMGenome*>(parent2);
    //ACMGenome *tmpChild1    = static_cast<ACMGenome*>(child1);
    //ACMGenome *tmpChild2    = static_cast<ACMGenome*>(child2);

    // Iteration variable(s)
    size_t i;
    size_t j;

    // Create a copy of availableIndices_ to sample without replacement
    std::vector<size_t> tmpAvailableIndices(availableIndices_);
    // Fill in the crossover points
    // We start by shuffling the temporal storage for indices
    std::shuffle(tmpAvailableIndices.begin(), tmpAvailableIndices.end(), gen);
    // Now we copy the first order_ elements to crossoverPoints_
    for (i = 0; i < order_; i++)
    {
        crossoverPoints_[i+1] = tmpAvailableIndices[i];
    }
    // Finally, sort the newly added crossover points
    std::sort(crossoverPoints_.begin()+1, crossoverPoints_.end()-1);
    // DONE! Sampled without replacement!

    // We now cross the genes of the individuals
    // Start by the regions that are not swapped
    for (i = 0; i < crossoverPoints_.size() - 1; i += 2)
    {
        for (j = crossoverPoints_[i]; j < crossoverPoints_[i+1]; j++)
        {
            child1->setBase(j, parent1->base(j));
            child2->setBase(j, parent2->base(j));
        }
    }
    // Now the regions that should be swapped
    for (i = 1; i < crossoverPoints_.size() - 1; i += 2)
    {
        for (j = crossoverPoints_[i]; j < crossoverPoints_[i+1]; j++)
        {
            child1->setBase(j, parent2->base(j));
            child2->setBase(j, parent1->base(j));
        }
    }
}


} //namespace alexandria
