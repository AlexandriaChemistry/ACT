#include "Initializer.h"

#include "aliases.h"


namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: SimpleInitializer                 *
* * * * * * * * * * * * * * * * * * * * * */

void SimpleInitializer::initialize(Individual *individual)
{
    for (int i = 0; i < length; i++) // To do: Get the length from the individual.
    {
        (*individual)[i] = dis(gen);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: SimpleInitializer                   *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga
