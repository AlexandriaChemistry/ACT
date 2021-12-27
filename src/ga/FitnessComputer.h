/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julianramon.marradesfurquet.8049@student.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef GA_FITNESSCOMPUTER_H
#define GA_FITNESSCOMPUTER_H

#include "Individual.h"


namespace ga
{


//! Target for the fitness computation
enum class Target
{
    //! Train fitness
    Train,
    //! Test fitness
    Test
};

/*!
 * Abstract class for computing the fitness of an individual
 */
class FitnessComputer
{

public:

    /*!
     * Compute the fitness of an individual
     * \param[in] individual    the individual
     * \param[in] trgtFit       the target for fitness computation. Either Train or Test
     */
    virtual void compute(Individual *individual,
                         Target      trgtFit) = 0;

};


} //namespace ga


#endif //GA_FITNESSCOMPUTER_H
