/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#ifndef GA_PROBABILITYCOMPUTER_H
#define GA_PROBABILITYCOMPUTER_H

#include <vector>

#include "Genome.h"

namespace ga
{
/*!
 * Abstract class for computing the selection probability of each Genome in the population
 */
class ProbabilityComputer
{

public:

    /*!
     * \brief Compute the selection probability of each genome in the population
     * \param[in] pop pointer to the population
     */
    virtual void compute(std::vector<Genome> *pop) = 0;
};

/*!
 * Class for fitness-based probability computation; that is, the probability is inversely proportional to the fitness.<br>
 * The probability of a genome \f$ i \f$ can be written as
 * \f[
 *      p_i = \frac { 1 / \left( \varepsilon + f_i \right) } { \sum \limits_{ j=1 }^{ n } 1 / \left( \varepsilon + f_j \right) },
 * \f]
 * where \f$ f_i \f$ is the fitness of genome \f$ i \f$ and \f$ n \f$ is the number of genomes in the population.
 */
class FitnessProbabilityComputer : public ProbabilityComputer
{

public:

    virtual void compute(std::vector<Genome> *pop);

};


/*!
 * Class for Boltzmann temperature probability computation.<br>
 * The probability of a genome \f$ i \f$ can be written as
 * \f[
 *      p_i = \frac{ e^{ f_i^* / T } }{ \sum \limits_{ j=1 }^{ n } e^{ f_j^* / T } },
 * \f]
 * where \f$ n \f$ is the number of genomes in the population, \f$ T \in \mathbb{R} \f$ is the Boltzmann
 * temperature parameter, and \f$ f_i^* \f$ is written as
 * \f[
 *      f_i^* = 1 / \left( \varepsilon + f_i \right),
 * \f]
 * where \f$ f_i \f$ is the fitness of genome \f$ i \f$.
 */
class BoltzmannProbabilityComputer : public ProbabilityComputer
{

private:
    //! The temperature parameter
    double temperature_;
    //! Will store \f$ e^{f_i^*} \f$ for each genome \f$ i \f$
    std::vector<double> exponentials_;

public:

    /*!
     * \brief Constructor
     * \param[in] temperature       the temperature
     * \param[in] popSize           amount of genomes in the population
     */
    BoltzmannProbabilityComputer(const double   temperature,
                                 const size_t   popSize)
        : temperature_(temperature), exponentials_(popSize) {}

    virtual void compute(std::vector<Genome> *pop);

};


/*!
 * Class for Rank probability computation, that is, the probability of a genome is proportional to its rank.
 * Rank 1 means highest fitness, rank 2 second highest, etc.<br>
 * If there are \f$ R \f$ ranks, then the probability of the genome at rank \f$ i \f$ is
 * \f[
 *      p_i = \frac{R - i + 1}{ R \left( R+1 \right) / 2 }
 * \f]
 */
class RankProbabilityComputer : public ProbabilityComputer
{

private:

    //! Stores the sum of ranks. E.g., if there are 10 genomes, the sum of ranks is \f$ 1 + 2  + 3 + \dots + 9 + 10 \f$
    double sumOfRanks_;

public:

    /*!
     * \brief
     * \param[in] popSize   number of genomes in the population
     */
    RankProbabilityComputer(const int popSize)
        : sumOfRanks_(popSize * (popSize + 1) / 2) {}

    /*!
     * \brief Here we assume that population and fitness are sorted.<br>
     * Compute the probability of each genome in the population
     * @param pop pointer to the population
     */
    virtual void compute(std::vector<Genome> *pop);

};


} //namespace ga


#endif //GA_PROBABILITYCOMPUTER_H
