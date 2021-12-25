#ifndef GA_PROBABILITYCOMPUTER_H
#define GA_PROBABILITYCOMPUTER_H


#include <vector>

#include "Individual.h"


namespace ga
{


/*!
 * Abstract class for computing the selection probability of each individual in the population
 */
class ProbabilityComputer
{

public:

    /*!
     * Compute the selection probability of each individual in the population
     * \param[in] pop pointer to the population
     */
    virtual void compute(std::vector<Individual*> *pop) = 0;

};


/*!
 * Class for fitness-based probability computation, that is, the probability is proportional to the fitness.<br>
 * The probability of individual \f$i\f$ can be written as
 * \f[
 *      p_i = \frac { f_i } { \sum \limits_{ j=1 }^{ n } f_j },
 * \f]
 * where \f$n\f$ is the number of individuals in the population.
 */
class FitnessProbabilityComputer : public ProbabilityComputer
{

public:

    virtual void compute(std::vector<Individual*> *pop);

};


/*!
 * Class for Boltzmann temperature probability computation.<br>
 * The probability of individual \f$i\f$ can be written as
 * \f[
 *      p_i = \frac{ e^{ f_i / T } }{ \sum \limits_{ j=1 }^{ n } e^{ f_j / T } },
 * \f]
 * where \f$n\f$ is the number of individuals in the population, and \f$T \in \mathbb{R}\f$ is the Boltzmann
 * temperature parameter.
 */
class BoltzmannProbabilityComputer : public ProbabilityComputer
{

private:

    //! The temperature parameter
    double temperature_;
    //! Stores e^fitness for each individual
    std::vector<double> exponentials_;

public:

    /*!
     * Create a new BoltzmannProbabilityComputer object
     * \param[in] popSize           number of individuals in the population
     * \param[in] temperature       the temperature
     */
    BoltzmannProbabilityComputer(const int      popSize,
                                 const double   temperature)
    : exponentials_(popSize), temperature_(temperature) {}

    virtual void compute(std::vector<Individual*> *pop);

};


/*!
 * Class for Rank probability computation, that is, the probability of an individual is proportional to its rank.
 * Rank 1 means highest fitness, rank 2 second highest, etc.<br>
 * If there are \f$R\f$ ranks, then the probability of the individual at rank \f$i\f$ is
 * \f[
 *      p_i = \frac{R - i + 1}{ R \left( R+1 \right) / 2 }
 * \f]
 */
class RankProbabilityComputer : public ProbabilityComputer
{

private:

    //! Stores the sum of ranks. E.g., if there are 10 individuals, the sum of ranks is 1 + 2  + 3 + ... + 9 + 10
    double sumOfRanks_;

public:

    /*!
     * Create a new RankProbabilityComputer object
     * \param[in] popSize   number of individuals in the population
     */
    RankProbabilityComputer(const int popSize)
    : sumOfRanks_(popSize * (popSize + 1) / 2) {}

    /*!
     * Here we assume that population and fitness are sorted in a descending manner.<br>
     * Compute the probability of each individual in the population
     * @param pop pointer to the population
     */
    virtual void compute(std::vector<Individual*> *pop);

};


} //namespace ga


#endif //GA_PROBABILITYCOMPUTER_H
