#ifndef ACT_ALIASES_H
#define ACT_ALIASES_H


#include <tuple>
#include <vector>


//! Regular 1D array of double precision floating point numbers
using vector = std::vector<double>;

//! 2D matrix of double precision floating point numbers
using matrix = std::vector<vector>;


namespace ga
{


/*!
 * Result of the evolution done by the genetic algorithm
 */
typedef struct ga_result
{
    //! Collection of individuals, each being a vector
    matrix  pop;
    //! Fitness score for each individual in the population
    vector  fitness;
    //! Individual with the best fitness
    vector  bestIndividual;
    //! Best fitness
    double  bestFitness;
    //! Amount of generations
    int     generations;
} ga_result_t;


} //namespace ga


#endif //ACT_ALIASES_H
