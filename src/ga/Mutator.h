/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#ifndef GA_MUTATOR_H
#define GA_MUTATOR_H

#include <random>
#include <time.h>

#include "Genome.h"

namespace ga
{
/*!
 * \brief Abstract class for gene mutation of Individual objects
 */
class Mutator
{
private:
    // Ranom number generation
    std::random_device                      rd_base;
    std::mt19937                            gen_base;
    std::uniform_real_distribution<double>  dis_base;

protected:

    /*! \brief Constructor
     * \param[in] seed seed for the random number generator
     */
    Mutator(int seed)
    : gen_base(rd_base()), dis_base(std::uniform_real_distribution<double>(0.0, 1.0))
    {
        gen_base.seed(seed);
    }

    //! \return a random number in \f$ [0, 1] \f$
    double randNum() { return dis_base(gen_base); }

public:

    /*!
     * \brief Mutate genes of an Individual (in place)
     * \param[inout] genome     Pointer to the genome to mutate
     * \param[out]   bestGenome Pointer to the best genome found
     * \param[in]    prMut      Probability of mutating a gene
     */
    virtual void mutate(Genome *genome,
                        Genome *bestGenome,
                        double  prMut) = 0;
    
    virtual void finalize() = 0;
    
    virtual bool foundMinimum() = 0;
};


} //namespace ga


#endif //GA_MUTATOR_H
