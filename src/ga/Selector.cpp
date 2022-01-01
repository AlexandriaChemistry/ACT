/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#include "Selector.h"

#include <random>


namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: RouletteSelector                  *
* * * * * * * * * * * * * * * * * * * * * */

int RouletteSelector::select(const std::vector<Individual*> &pop)
{
    double num = dis(gen);
    size_t i = 0;
    while (num > 0 and i < pop.size())
    {
        num -= pop[i]->probability();
        i++;
    }
    return i - 1;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: RouletteSelector                    *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga
