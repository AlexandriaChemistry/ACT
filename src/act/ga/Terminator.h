/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_TERMINATOR_H
#define GA_TERMINATOR_H


#include <vector>

#include "GenePool.h"


namespace ga
{


/*!
 * \brief Abstract class to check for evolution termination conditions
 */
class Terminator
{

public:

    /*!
     * \brief Check whether the evolution should be terminated
     * \param[in] pool              The GenePool
     * \param[in] generationNumber  The generation number
     * \return true if we should terminate, false otherwise
     */
    virtual bool terminate(const GenePool *pool,
                           const int       generationNumber) = 0;

};

/*!
 * \brief Terminator which stops evolution after a given amount of generations.
 */
class GenerationTerminator : public Terminator
{

private:

    //! Maximum allowed amount of generations
    int maxGenerations_;

public:

    /*!
     * \brief Constructor
     * @param maxGenerations the maximum allowed amount of generations
     */
    GenerationTerminator(const int maxGenerations)
    : maxGenerations_(maxGenerations) {}

    /*!
     * Will return true when \p generationNumber \f$\geq\f$ \p maxGenerations, and false otherwise.
     * \param[in] pool             The gene pool
     * \param[in] generationNumber The generation number
     * \return true if we should terminate, false otherwise
     */
    virtual bool terminate(const GenePool *pool,
                           const int       generationNumber);

};


} //namespace ga


#endif //GA_TERMINATOR_H
