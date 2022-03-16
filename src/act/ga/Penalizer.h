#ifndef GA_PENALIZER_H
#define GA_PENALIZER_H

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "GenePool.h"
#include "Initializer.h"

#include "gromacs/fileio/xvgr.h"

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

    //! \brief Default base destructor
    virtual ~Penalizer() {}

    /*!
     * \brief Penalize a population
     * \param[inout] pool       the collection of genomes
     * \param[in]    generation the current generation number
     * \return true if the population has been penalized, false otherwise
     */
    virtual bool penalize(      GenePool *pool,
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

    //! Whether we compute the volume in log scale
    bool logVolume_;
    //! Total volume of the parameter space (may be log scale)
    double totalVolume_;
    //! Limit of fraction of volume
    double volFracLimit_;
    //! Fraction of the population to reinitialize (worst genomes)
    double popFrac_;
    //! Initializer to randomize genomes
    Initializer *initializer_;

    //! Output file
    FILE *outfile_ = nullptr;

    //! \return the volume of a population
    double getPoolVolume(const GenePool &pool) const;

public:

    //! \brief Destructor: closes output file
    ~VolumeFractionPenalizer();

    /*!
     * \brief Create a new VolumeFractionPenalizer
     * \param[in] oenv         gromacs output environment
     * \param[in] logVolume    true if volume in log scale, false otherwise
     * \param[in] logfile      file to print log stuff to (may be nullptr)
     * \param[in] totalVolume  total volume of the parameter space
     * \param[in] volFracLimit if the volume of the population divided by
     *                         \p totalVolume is smaller than \p volFracLimit then
     *                         we punish
     * \param[in] popFrac      fraction of the population to penalize (worst genomes)
     * \param[in] initializer  Initializer to randomize genomes
     */
    VolumeFractionPenalizer(      gmx_output_env_t  *oenv,
                            const bool               logVolume,
                                  FILE              *logfile,
                            const double             totalVolume,
                            const double             volFracLimit,
                            const double             popFrac,
                                  Initializer       *initializer);

    bool penalize(      GenePool *pool,
                  const int       generation) override;

}; // class VolumeFactorPenalizer

/*!
 * \brief A catastrophe happens and a fraction of the population (at random) is
 * randomized
 * 
 * Each X generations a fraction Y (at random) of the population is randomized
 * TODO: if needed, make this penalizer print the indices of the individuals it is
 * randomizing to a file
 */
class CatastrophePenalizer : public Penalizer
{

private:

    // Random number stuff (for shuffling)
    std::random_device  rd;
    std::mt19937        gen;

    //! Interval of generations between the penalizer triggers
    int genInterval_;
    //! Fraction of the population to reinitialize
    double popFrac_;
    //! Initializer to randomize genomes
    Initializer *initializer_;
    //! List of all existing indices for genomes in the population (0 to popSize-1)
    std::vector<size_t> availableIndices_;

public:

    /*!
     * \brief Create a new CatastrophePenalizer
     * \param[in] logfile     file to print log stuff to (may be nullptr)
     * \param[in] seed        seed for the random number generation
     * \param[in] genInterval number of generations as interval between penalties
     * \param[in] popFrac     fraction of the population to penalize
     * \param[in] initializer Initializer to randomize genomes
     * \param[in] popSize     Size of the population
     */
    CatastrophePenalizer(      FILE        *logfile,
                         const int          seed,
                         const int          genInterval,
                         const double       popFrac,
                               Initializer *initializer,
                         const size_t       popSize);

    bool penalize(      GenePool *pool,
                  const int       generation) override;

}; // class CatastrophePenalizer

} // namespace ga


#endif // GA_PENALIZER_H