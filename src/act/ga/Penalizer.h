#ifndef GA_PENALIZER_H
#define GA_PENALIZER_H

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "GenePool.h"
#include "Initializer.h"

namespace ga
{

/*!
 * \brief Penalizes a population
 */
class Penalizer
{

protected:

    //! File to print stuff to (may be nullptr)
    FILE *logfile_ = nullptr;

    //! \brief Default constructor
    Penalizer(FILE *logfile): logfile_(logfile) {}

public:

    /*!
     * \brief Penalize a population
     * \param[inout] pool       the collection of genomes
     * \param[in]    generation the current generation number
     */
    virtual void penalize(      GenePool *pool,
                          const int       generation) = 0;

}; // class Penalizer

/*!
 * \brief Penalizer that punishes the population if its volume is not large enough
 *
 * We create a parallelogram whose side lengths are defined by the maximum and
 * minimum value of each parameter in the population.
 * If the volume of the parallelogram divided by the total volume of the parameter
 * space is smaller than a given fraction, we reinitialize a given fraction of the
 * worst genomes.
 */
class VolumeFractionPenalizer : public Penalizer
{

private:

    //! Total volume of the parameter space
    double totalVolume_;
    //! Limit of fraction of volume
    double volFracLimit_;
    //! Fraction of the population to reinitialize (worst genomes)
    double popFrac_;
    //! Initializer to randomize genomes
    Initializer *initializer_;

    //! \return the volume of a population
    double getPoolVolume(const GenePool &pool) const;

public:

    /*!
     * \brief Create a new VolumeFractionPenalizer
     * \param[in] totalVolume  total volume of the parameter space
     * \param[in] volFracLimit if the volume of the population divided by
     *                         \p totalVolume is smaller than \p volFracLimit then
     *                         we punish
     * \param[in] popFrac      fraction of the population to penalize (worst genomes)
     * \param[in] initializer  Initializer to randomize genomes
     */
    VolumeFractionPenalizer(      FILE        *logfile,
                            const double       totalVolume,
                            const double       volFracLimit,
                            const double       popFrac,
                                  Initializer *initializer);

    void penalize(      GenePool *pool,
                  const int       generation) override;

}; // class VolumeFactorPenalizer

} // namespace ga


#endif // GA_PENALIZER_H