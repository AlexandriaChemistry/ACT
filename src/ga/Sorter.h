#ifndef ACT_SORTER_H
#define ACT_SORTER_H

#include "aliases.h"


/*!
 * Abstract class for sorting the population based on fitness
 */
class Sorter {

public:
    /*!
     * Sort individuals (in place) based on fitness
     * @param pop                   the population. Each row is an individual.
     * @param fitness               the fitness of each individual in the population
     * @param popSize               number of individuals in the population
     */
    virtual void sort(matrix pop, vector fitness, const int popSize);

};


/*!
 * Sorter which does nothing.
 */
class EmptySorter : public Sorter {

public:
    void sort(matrix pop, vector fitness, const int popSize) {};

};


/*!
 * Class for Merge-sort
 */
class MergeSorter: public Sorter {

    double chromosomeLength;
    matrix tmpPop;
    vector tmpFitness;

public:
    /*!
     * Create a new MergeSorter object
     * @param popSize               number of individuals in the population
     * @param chromosomeLength      the size of each individual
     */
    MergeSorter(const int popSize, const int chromosomeLength);

    void sort(matrix pop, vector fitness, const int popSize);

    /*!
     * Split <fitA> into 2 runs, sort both runs into <fitB>, merge both runs from <fitB> into <fitA>
     * @param popB      population B
     * @param fitB      fitness B
     * @param left      left index (inclusive)
     * @param right     right index (exclusive)
     * @param popA      population A
     * @param fitA      fitness A
     */
    void topDownSplitMerge(matrix popB, vector fitB, const int left, const int right, matrix popA, vector fitA);

    /*!
     * Left source half is A[left:middle-1].
     * Right source half is A[middle:right-1].
     * Result is B[left:right-1].
     * @param popA          population A
     * @param fitA          fitness A
     * @param left          left index (inclusive)
     * @param middle        middle index
     * @param right         right index (exclusive)
     * @param popB          population B
     * @param fitB          fitness B
     */
    void topDownMerge(matrix popA, vector fitA, const int left, const int middle, const int right, matrix popB,
                      vector fitB);

};

#endif //ACT_SORTER_H
