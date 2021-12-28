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
 * \brief Abstract class for sorting the population based on the fitness of each Individual
 */
class Sorter
{

protected:

      //! Whether to sort in descending order of fitness (default false)
      bool descending_ = false;

      //! \brief Default constructor
      Sorter() {}

      /*!
       * \brief Property constructor
       * \param[in] descending whether to sort in descending order of fitness
       */
      Sorter(const bool descending)
      : descending_(descending) {}

public:

      /*!
       * \brief Sort the population
       * \param[in] pop pointer to the population
       */
      virtual void sort(std::vector<Individual*> *pop) = 0;

};


/*!
 * \brief Sorter which does nothing. Will be used when sorting is not required.
 */
class EmptySorter : public Sorter
{

public:

      virtual void sort(gmx_unused std::vector<Individual*> *pop) {};

};


/*!
 * \brief Class for Merge-sort
 */
class MergeSorter : public Sorter
{

private:

      //! Temporal storage for a population
      std::vector<Individual*> tmpPop_;

      /*!
       * \brief Split \p fitA into 2 runs, sort both runs into \p fitB, merge both runs from \p fitB into \p fitA
       * \param[in] popB      pointer to population B
       * \param[in] left      left index (inclusive)
       * \param[in] right     right index (exclusive)
       * \param[in] popA      pointer to population A
       */
      void topDownSplitMerge(      std::vector<Individual*>      *popB,
                             const int                            left,
                             const int                            right,
                                   std::vector<Individual*>      *popA);

      /*!
       * \brief Merge both runs from \p fitA into \p fitB <br>
       * Left source half is A[left:middle-1].
       * Right source half is A[middle:right-1].
       * Result is B[left:right-1].
       * \param[in] popA      pointer to population A
       * \param[in] left      left index (inclusive)
       * \param[in] middle    middle index
       * \param[in] right     right index (exclusive)
       * \param[in] popB      pointer to population B
       */
      void topDownMerge(      std::vector<Individual*>     *popA,
                        const int                           left,
                        const int                           middle,
                        const int                           right,
                              std::vector<Individual*>     *popB);

public:

      /*!
       * \brief Constructor
       * \param[in] popSize         number of individuals in the population
       * \param[in] descending      whether we sort in descending order
       */
      MergeSorter(const int   popSize,
                  const bool  descending)
      : Sorter(descending), tmpPop_(popSize) {}

      virtual void sort(std::vector<Individual*> *pop);

};


/*!
 * \brief Class for Quick-sort
 */
class QuickSorter : public Sorter
{

private:

      //! Temporal storage for an individual
      Individual *tmpInd_;

      /*!
       * \brief Split the population into 2 parts, one left of the pivot element and one to the right of it, and sort both.
       * \param[in] pop     pointer to the population
       * \param[in] low     the left-most point of the part of the population vector in this recursive run
       * \param[in] high    the right-most point of the part of the population vector in this recursive run
       */
      void quickSort(      std::vector<Individual*>  *pop,
                     const int                        low,
                     const int                        high);

      /*!
       * \brief Find the pivot element and sort everything by comparing with it in an ascending order.
       * \param[in] pop       pointer to the population
       * \param[in] low       the left-most point of the part of the population vector in this recursive run
       * \param[in] high      the right-most point of the part of the population vector in this recursive run
       */
      int ascendingPartition(      std::vector<Individual*>      *pop,
                             const int                            low,
                             const int                            high);

      /*!
       * \brief Find the pivot element and sort everything by comparing with it in a descending order.
       * \param[in] pop       pointer to the population
       * \param[in] low       the left-most point of the part of the population vector in this recursive run
       * \param[in] high      the right-most point of the part of the population vector in this recursive run
       */
      int descendingPartition(      std::vector<Individual*>     *pop,
                              const int                           low,
                              const int                           high);

public:

      /*!
       * \brief Constructor
       * \param[in] descending whether we sort in descending order
       */
      QuickSorter(const bool descending)
      : Sorter(descending) {}

      virtual void sort(std::vector<Individual*> *pop);

};


} //namespace ga


#endif //GA_SORTER_H
