#ifndef GA_TERMINATOR_H
#define GA_TERMINATOR_H


#include <vector>

#include "Individual.h"


namespace ga
{


/*!
 * Abstract class to check for evolution termination conditions
 */
class Terminator
{

public:

    /*!
     * Check whether the evolution should be terminated
     * \param[in] population            the population
     * \param[in] generationNumber      the generation number
     * \return true if we should terminate, false otherwise
     */
    virtual bool terminate(const std::vector<Individual*>  &pop,
                           const int                        generationNumber) = 0;

};

/*!
 * Terminator which stops evolution after a given amount of generations.
 */
class GenerationTerminator : public Terminator
{

private:

    int maxGenerations_;

public:

    /*!
     * Create a new GenerationTerminator object
     * @param maxGenerations    the maximum amount of generations
     */
    GenerationTerminator(const int maxGenerations)
    : maxGenerations_(maxGenerations) {}

    /*!
     * Will return true when \p generationNumber \f$\geq\f$ \p maxGenerations, and false otherwise.
     * \param[in] population            the population
     * \param[in] generationNumber      the generation number
     * \return true if we should terminate, false otherwise
     */
    virtual bool terminate(const std::vector<Individual*>  &pop,
                           const int                        generationNumber)

};


} //namespace ga


#endif //GA_TERMINATOR_H
