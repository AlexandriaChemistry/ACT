#include "Sorter.h"

#include "aliases.h"
#include "helpers.h"


MergeSorter::MergeSorter(const int  popSize,
                         const int  chromosomeLength) {

    tmpFitness = vector(popSize);
    tmpPop = allocateMatrix(popSize, chromosomeLength);
    this->chromosomeLength = chromosomeLength;

}

void MergeSorter::sort(matrix&       pop,
                       vector&       fitness,
                       const int     popSize) {

    copyVectorValues(fitness, tmpFitness, 0, popSize);
    copyMatrixValues(pop, tmpPop, 0, popSize, 0, chromosomeLength);
    topDownSplitMerge(tmpPop, tmpFitness, 0, popSize, pop, fitness);

}


void MergeSorter::topDownSplitMerge(matrix&      popB,
                                    vector&      fitB,
                                    const int    left,
                                    const int    right,
                                    matrix&      popA,
                                    vector&      fitA) {

    if (right - left <= 1) return;

    const int middle = (right + left)/2;

    topDownSplitMerge(popA, fitA, left, middle, popB, fitB);
    topDownSplitMerge(popA, fitA, middle, right, popB, fitB);

    topDownMerge(popB, fitB, left, middle, right, popA, fitA);

}


void MergeSorter::topDownMerge(matrix&       popA,
                               vector&       fitA,
                               const int     left,
                               const int     middle,
                               const int     right,
                               matrix&       popB,
                               vector&       fitB) {

    int i = left;
    int j = middle;

    // While there are elements in the left or right runs...
    for (int k = left; k < right; k++) {
        // If left run head exists and is >= existing right run head.
        if (i < middle && (j >= right || fitA[i] >= fitA[j])) {
            fitB[k] = fitA[i];
            copyVectorValues(popA[i], popB[k], 0, chromosomeLength);
            i++;
        } else {
            fitB[k] = fitA[j];
            copyVectorValues(popA[j], popB[k], 0, chromosomeLength);
            j++;
        }
    }

}

QuickSorter::QuickSorter(const int  popSize,
                         const int  chromosomeLength) {

    tmpFitness = vector(popSize);
    tmpPop = allocateMatrix(popSize, chromosomeLength);
    this->chromosomeLength = chromosomeLength;

}

void QuickSorter::sort(matrix&       pop,
                       vector&       fitness,
                       const int     popSize) {
    int low = 0;
    int high = popSize;
    quickSort(fitness, low, high);
}

void QuickSorter::quickSort(vector&       fitness,
               const int     low,
               const int     high) {
    if (low >= 0 && high >= 0 && low < high) {
        int p = partition(fitness, low, high);
        quickSort(fitness, low, p - 1);
        quickSort(fitness, p + 1, high);
    }
}

int QuickSorter::partition(vector&        fitness,
              const int      low,
              const int      high) {
    double pivot = fitness[high];
    int candidate = low - 1;
    double temp;

    for (int check = low; check < high; check++) {
        if (fitness[check] <= pivot) {
            candidate = candidate + 1;
            temp = fitness[candidate];
            fitness[candidate] = fitness[check];
            fitness[check] = temp;
        }
    }

    return candidate;
}





