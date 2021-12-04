#ifndef ACT_SORTER_H
#define ACT_SORTER_H

#include "aliases.h"


namespace ga
{

    /*!
     * Abstract class for sorting the population based on fitness
     */
    class Sorter
    {

    protected:

        //! Whether to sort in descending order of fitness (default false)
        bool descending = false;

        //! Default constructor
        Sorter() {}

        /*!
         * Create a new Sorter
         * @param descending whether to sort in descending order of fitness
         */
        Sorter(const bool descending)
        {
              this->descending = descending;
        }

    public:

        /*!
         * Sort individuals (in place) based on fitness
         * @param pop                   pointer to the population. Each row is an individual.
         * @param fitness               pointer to the fitness of each individual in the population
         * @param popSize               number of individuals in the population
         */
        virtual void sort(      matrix *pop,
                                vector *fitness,
                          const int     popSize) = 0;

    };


    /*!
     * Sorter which does nothing. Will be used when sorting is not required.
     */
    class EmptySorter : public Sorter
    {

    public:
        virtual void sort(      matrix *pop,
                                vector *fitness,
                          const int     popSize) {};

    };


    /*!
     * Class for Merge-sort
     */
    class MergeSorter : public Sorter
    {

        //! Temporal storage for a population
        matrix tmpPop;
        //! Temporal vector for fitness
        vector tmpFitness;

        /*!
         * Split \p fitA into 2 runs, sort both runs into \p fitB, merge both runs from \p fitB into \p fitA
         * @param popB      pointer to population B
         * @param fitB      pointer to fitness B
         * @param left      left index (inclusive)
         * @param right     right index (exclusive)
         * @param popA      pointer to population A
         * @param fitA      pointer to fitness A
         */
        void topDownSplitMerge(      matrix    *popB,
                                     vector    *fitB,
                               const int        left,
                               const int        right,
                                     matrix    *popA,
                                     vector    *fitA);

        /*!
         * Left source half is A[left:middle-1].
         * Right source half is A[middle:right-1].
         * Result is B[left:right-1].
         * @param popA          pointer to population A
         * @param fitA          pointer to fitness A
         * @param left          left index (inclusive)
         * @param middle        middle index
         * @param right         right index (exclusive)
         * @param popB          pointer to population B
         * @param fitB          pointer to fitness B
         */
        void topDownMerge(      matrix     *popA,
                                vector     *fitA,
                          const int         left,
                          const int         middle,
                          const int         right,
                                matrix     *popB,
                                vector     *fitB);

    public:
        /*!
         * Create a new MergeSorter object
         * @param popSize               number of individuals in the population
         * @param chromosomeLength      the size of each individual
         * @param descending            true if sorting in descending order (using fitness),
         *                              false otherwise (using chi2)
         */
        MergeSorter(const int   popSize,
                    const int   chromosomeLength,
                    const bool  descending);

        virtual void sort(      matrix *pop,
                                vector *fitness,
                          const int     popSize);

    };


    /*!
     * Class for Quick-sort
     */
    class QuickSorter : public Sorter
    {

        //! Temporal storage for fitness
        vector tmpFitness;

        /*!
         * Split \p fitness into 2 parts, one left of the pivot element and one to the right of it, and sort both.
         * @param pop       pointer to the population
         * @param fitness   pointer to the fitness vector
         * @param low       the left-most point of the part of the population vector in this recursion
         * @param high      the right-most point of the part of the population vector in this recursion
         */
        void quickSort(      matrix    *pop,
                             vector    *fitness,
                       const int        low,
                       const int        high);

        /*!
         * Find the pivot element and sort everything by comparing with it in an ascending order.
         * @param pop       pointer to the population
         * @param fitness   pointer to the fitness vector
         * @param low       the left-most point of the part of the population vector in this recursion
         * @param high      the right-most point of the part of the population vector in this recursion
         */
        int ascendingPartition(      matrix     *pop,
                                     vector     *fitness,
                               const int         low,
                               const int         high);

         /*!
          * Find the pivot element and sort everything by comparing with it in a descending order.
          * @param pop       pointer to the population
          * @param fitness   pointer to the fitness vector
          * @param low       the left-most point of the part of the population vector in this recursion
          * @param high      the right-most point of the part of the population vector in this recursion
          */
         int descendingPartition(      matrix     *pop,
                                       vector     *fitness,
                                 const int         low,
                                 const int         high);
    public:
        /*!
         * Create a new QuickSorter
         * @param popSize       amount of individuals in the population
         * @param descending    true if sorting in descending order (using fitness),
         *                      false otherwise (using chi2)
         */
        QuickSorter(const int   popSize,
                    const bool  descending);


        virtual void sort(      matrix *pop,
                                vector *fitness,
                          const int     popSize);

    };

}

#endif //ACT_SORTER_H
