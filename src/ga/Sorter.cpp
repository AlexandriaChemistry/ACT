#include "Sorter.h"

#include "aliases.h"
#include "helpers.h"


MergeSorter::MergeSorter(const int  popSize,
                         const int  chromosomeLength) {

    tmpFitness = vector(popSize);
    tmpPop = allocateMatrix(popSize, chromosomeLength);
    this->chromosomeLength = chromosomeLength;

}


void MergeSorter::sort(matrix       pop,
                       vector       fitness,
                       const int    popSize) {

    copyVectorValues(fitness, tmpFitness, 0, popSize);
    copyMatrixValues(pop, tmpPop, 0, popSize, 0, chromosomeLength);
    topDownSplitMerge(tmpPop, tmpFitness, 0, popSize, pop, fitness);

}


void MergeSorter::topDownSplitMerge(matrix      popB,
                                    vector      fitB,
                                    const int   left,
                                    const int   right,
                                    matrix      popA,
                                    vector      fitA) {

    if (right - left <= 1) return;

    const int middle = (right + left)/2;

    topDownSplitMerge(popA, fitA, left, middle, popB, fitB);
    topDownSplitMerge(popA, fitA, middle, right, popB, fitB);

    topDownMerge(popB, fitB, left, middle, right, popA, fitA);

}


void MergeSorter::topDownMerge(matrix       popA,
                               vector       fitA,
                               const int    left,
                               const int    middle,
                               const int    right,
                               matrix       popB,
                               vector       fitB) {

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
