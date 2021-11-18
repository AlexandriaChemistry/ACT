#ifndef ACT_ALIASES_H
#define ACT_ALIASES_H

#include <tuple>
#include <vector>

/*!
 * Regular 1D array of double precision floating point numbers
 */
using vector = std::vector<double>;

/*!
 * 2D matrix of double precision floating point numbers
 */
using matrix = std::vector<vector>;

namespace ga
{

    /*!
     * Result structure of genetic algorithm
     */
    typedef struct ga_result
    {
        matrix  pop;
        vector  fitness;
        vector  bestIndividual;
        double  bestFitness;
        int     generations;
    } ga_result_t;

}

#endif //ACT_ALIASES_H
