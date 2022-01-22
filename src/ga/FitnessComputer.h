/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_FITNESSCOMPUTER_H
#define GA_FITNESSCOMPUTER_H

#include "Dataset.h"
#include "Genome.h"

namespace ga
{

/*!
 * \brief Abstract class for computing the fitness of an individual
 */
class FitnessComputer
{

public:

    /*!
     * \brief Compute the fitness of a genome
     * \param[in] genome  The genome
     * \param[in] trgtFit The target for fitness computation. Either Train or Test
     */
    virtual void compute(Genome     *genome,
                         iMolSelect  trgtFit) = 0;

};


} //namespace ga


#endif //GA_FITNESSCOMPUTER_H
