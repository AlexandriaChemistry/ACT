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
    int i = 0;
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
