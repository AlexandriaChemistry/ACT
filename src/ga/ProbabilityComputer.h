#ifndef ACT_PROBABILITYCOMPUTER_H
#define ACT_PROBABILITYCOMPUTER_H

#include "aliases.h"


/*!
 * Abstract class for computing the selection probability of each individual in the population
 */
class ProbabilityComputer {

public:
    /*!
     * Compute the probability of each individual in the population
     * @param fitness   the fitness of each individual
     * @param prob      structure to store the probabilities
     * @param popSize   size of the population
     */
    virtual void compute(const vector fitness, const vector prob, const int popSize);

};


/*!
 * Class for fitness-based probability computation
 */
class FitnessProbabilityComputer : public ProbabilityComputer {

public:
    void compute(const vector fitness, const vector prob, const int popSize);

};


/*!
 * Class for Boltzmann temperature probability computation
 */
class BoltzmannProbabilityComputer : public ProbabilityComputer {

    double temperature;
    vector exponentials;

public:
    /*!
     * Create a new BoltzmannProbabilityComputer object
     * @param popSize           number of individuals in the population
     * @param temperature       the temperature
     */
    BoltzmannProbabilityComputer(const int popSize, const double temperature);

    void compute(const vector fitness, const vector prob, const int popSize);

};


/*!
 * Class for Rank probability computation
 */
class RankProbabilityComputer : public ProbabilityComputer {

    double sumOfRanks;

public:
    /*!
     * Create a new RankProbabilityComputer object
     * @param popSize       number of individuals in the population
     */
    RankProbabilityComputer(const int popSize);

    /*!
     * Here we assume that population and fitness are sorted in a descending manner
     * @param fitness
     * @param prob
     * @param popSize
     */
    void compute(const vector fitness, const vector prob, const int popSize);

};

#endif //ACT_PROBABILITYCOMPUTER_H
