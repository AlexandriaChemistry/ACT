#include "npointcrossover.h"

#include "acmindividual.h"


namespace alexandria
{


void NPointCrossover::offspring(ga::Individual  *parent1,
                                ga::Individual  *parent2,
                                ga::Individual  *child1,
                                ga::Individual  *child2)
{

    // Casting individuals
    ACMIndividual *tmpParent1   = static_cast<ACMIndividual*>(parent1);
    ACMIndividual *tmpParent2   = static_cast<ACMIndividual*>(parent2);
    ACMIndividual *tmpChild1    = static_cast<ACMIndividual*>(child1);
    ACMIndividual *tmpChild2    = static_cast<ACMIndividual*>(child2);

    // Iteration variable(s)
    size_t i;
    size_t j;

    // Create a copy of availableIndices_ to sample without replacement
    std::vector<size_t> tmpAvailableIndices(availableIndices_);
    // Fill in the crossover points
    // We start by shuffling the temporal storage for indices
    std::shuffle(tmpAvailableIndices.begin(), tmpAvailableIndices.end(), gen);
    // Now we copy the first order_ elements to crossoverPoints_
    for (i = 0; i < order_; i++) crossoverPoints_[i+1] = tmpAvailableIndices[i];
    // Finally, sort the newly added crossover points
    std::sort(crossoverPoints_.begin()+1, crossoverPoints_.end()-1);
    // DONE! Sampled without replacement!

    // We now cross the genes of the individuals
    // Start by the regions that are not swapped
    for (i = 0; i < crossoverPoints_.size() - 1; i += 2)
    {
        for (j = crossoverPoints_[i]; j < crossoverPoints_[i+1]; j++)
        {
            tmpChild1->setParam(j, tmpParent1->paramAtIndex(j));
            tmpChild2->setParam(j, tmpParent2->paramAtIndex(j));
        }
    }
    // Now the regions that should be swapped
    for (i = 1; i < crossoverPoints_.size() - 1; i += 2)
    {
        for (j = crossoverPoints_[i]; j < crossoverPoints_[i+1]; j++)
        {
            tmpChild1->setParam(j, tmpParent2->paramAtIndex(j));
            tmpChild2->setParam(j, tmpParent1->paramAtIndex(j));
        }
    }

}


} //namespace alexandria
