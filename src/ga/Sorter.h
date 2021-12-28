/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef GA_SORTER_H
#define GA_SORTER_H


#include <vector>

#include "Individual.h"

#include "gromacs/utility/basedefinitions.h"


namespace ga
{


/*!
 * Abstract class for sorting the population based on fitness
 */
class Sorter
{

protected:

      //! Whether to sort in descending order of fitness (default false)
      bool descending_ = false;

      //! Default constructor
      Sorter() {}

      /*!
       * Property constructor
       * @param descending whether to sort in descending order of fitness
       */
      Sorter(const bool descending)
      : descending_(descending) {}

public:

      /*!
       * Sort individuals (in place) based on fitness
       * \param[in] pop                   pointer to the population
       */
      virtual void sort(std::vector<Individual*> *pop) = 0;

};


/*!
 * Sorter which does nothing. Will be used when sorting is not required.
 */
class EmptySorter : public Sorter
{

public:

      virtual void sort(gmx_unused std::vector<Individual*> *pop) {};

};


/*!
 * Class for Merge-sort
 */
class MergeSorter : public Sorter
{

private:

      //! Temporal storage for a population
      std::vector<Individual*> tmpPop_;

      /*!
       * Split \p fitA into 2 runs, sort both runs into \p fitB, merge both runs from \p fitB into \p fitA
       * @param popB      pointer to population B
       * @param left      left index (inclusive)
       * @param right     right index (exclusive)
       * @param popA      pointer to population A
       */
      void topDownSplitMerge(      std::vector<Individual*>      *popB,
                             const int                            left,
                             const int                            right,
                                   std::vector<Individual*>      *popA);

      /*!
       * Left source half is A[left:middle-1].
       * Right source half is A[middle:right-1].
       * Result is B[left:right-1].
       * @param popA          pointer to population A
       * @param left          left index (inclusive)
       * @param middle        middle index
       * @param right         right index (exclusive)
       * @param popB          pointer to population B
       */
      void topDownMerge(      std::vector<Individual*>     *popA,
                        const int                           left,
                        const int                           middle,
                        const int                           right,
                              std::vector<Individual*>     *popB);

public:

      /*!
       * Create a new MergeSorter object
       * \param[in] popSize               number of individuals in the population
       * \param[in] descending            true if sorting in descending order (using regular fitness),
       *                                  false otherwise (using chi2)
       */
      MergeSorter(const int   popSize,
                  const bool  descending)
      : Sorter(descending), tmpPop_(popSize) {}

      virtual void sort(std::vector<Individual*> *pop);

};


/*!
 * Class for Quick-sort
 */
class QuickSorter : public Sorter
{

private:

      //! Temporal storage for an individual
      Individual *tmpInd_;

      /*!
       * Split \p fitness into 2 parts, one left of the pivot element and one to the right of it, and sort both.
       * @param pop       pointer to the population
       * @param low       the left-most point of the part of the population vector in this recursion
       * @param high      the right-most point of the part of the population vector in this recursion
       */
      void quickSort(      std::vector<Individual*>  *pop,
                     const int                        low,
                     const int                        high);

      /*!
       * Find the pivot element and sort everything by comparing with it in an ascending order.
       * @param pop       pointer to the population
       * @param low       the left-most point of the part of the population vector in this recursion
       * @param high      the right-most point of the part of the population vector in this recursion
       */
      int ascendingPartition(      std::vector<Individual*>      *pop,
                             const int                            low,
                             const int                            high);

      /*!
       * Find the pivot element and sort everything by comparing with it in a descending order.
       * @param pop       pointer to the population
       * @param low       the left-most point of the part of the population vector in this recursion
       * @param high      the right-most point of the part of the population vector in this recursion
       */
      int descendingPartition(      std::vector<Individual*>     *pop,
                              const int                           low,
                              const int                           high);

public:

      /*!
       * Create a new QuickSorter
       * \param[in] descending      true if sorting in descending order (using fitness),
       *                            false otherwise (using chi2)
       */
      QuickSorter(const bool descending)
      : Sorter(descending) {}

      virtual void sort(std::vector<Individual*> *pop);

};


} //namespace ga


#endif //GA_SORTER_H
