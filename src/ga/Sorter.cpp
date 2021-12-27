#include "Sorter.h"


namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: MergeSorter                       *
* * * * * * * * * * * * * * * * * * * * * */

void MergeSorter::sort(std::vector<Individual*> *pop)
{

    // Clone the population into tmpPop_
    for (int i = 0; i < pop->size(); i++) tmpPop_[i] = (*pop)[i]->clone();
    topDownSplitMerge(&tmpPop_, 0, pop->size(), pop);

}


void MergeSorter::topDownSplitMerge(      std::vector<Individual*> *popB,
                                    const int                       left,
                                    const int                       right,
                                          std::vector<Individual*> *popA)
{

    if (right - left <= 1) return;

    const int middle = (right + left) / 2;

    topDownSplitMerge(popA, left, middle, popB);
    topDownSplitMerge(popA, middle, right, popB);

    topDownMerge(popB, left, middle, right, popA);

}


void MergeSorter::topDownMerge(      std::vector<Individual*>  *popA,
                               const int                        left,
                               const int                        middle,
                               const int                        right,
                                     std::vector<Individual*>  *popB)
{

    int i = left;
    int j = middle;

    for (int k = left; k < right; k++)
    {
        if ( i < middle && ( j >= right ||
            ( descending_ == ( (*popA)[i]->fitnessTrain() >= (*popA)[j]->fitnessTrain() ) ) ) )
        {
            (*popB)[k] = (*popA)[i]->clone();
            i++;
        }
        else
        {
            (*popB)[k] = (*popA)[j]->clone();
            j++;
        }
    }

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: MergeSorter                         *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: QuickSorter                       *
* * * * * * * * * * * * * * * * * * * * * */

void QuickSorter::sort(std::vector<Individual*> *pop)
{
    const int low = 0;
    const int high = pop->size() - 1;
    quickSort(pop, low, high);
}

void QuickSorter::quickSort(      std::vector<Individual*> *pop,
                            const int                       low,
                            const int                       high)
{
    if (low >= 0 && high >= 0 && low < high)
    {
        if (!descending_) {
            const int p = ascendingPartition(pop, low, high);
            quickSort(pop, low, p - 1);
            quickSort(pop, p + 1, high);
        } else {
            const int p = descendingPartition(pop, low, high);
            quickSort(pop, low, p - 1);
            quickSort(pop, p + 1, high);
        }
    }
}

int QuickSorter::ascendingPartition(      std::vector<Individual*> *pop,
                                    const int                       low,
                                    const int                       high)
{
    const double pivot = (*pop)[high]->fitnessTrain();
    int candidate = low - 1;
    double temp;

    for (int check = low; check <= high; check++)
    {
        if ( (*pop)[check]->fitnessTrain() <= pivot )
        {
            candidate += 1;
            tmpInd_ = (*pop)[candidate];
            (*pop)[candidate] = (*pop)[check];
            (*pop)[check] = tmpInd_;
        }
    }
    return candidate;
}

int QuickSorter::descendingPartition(      std::vector<Individual*> *pop,
                                     const int        low,
                                     const int        high)
{
    const double pivot = (*pop)[high]->fitnessTrain();
    int candidate = low - 1;
    double temp;

    for (int check = low; check <= high; check++)
    {
        if ( (*pop)[check]->fitnessTrain() >= pivot )
        {
            candidate += 1;
            tmpInd_ = (*pop)[candidate];
            (*pop)[candidate] = (*pop)[check];
            (*pop)[check] = tmpInd_;
        }
    }
    return candidate;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: QuickSorter                         *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga
