#include "Initializer.h"

#include "aliases.h"


namespace ga
{

    void SimpleInitializer::initialize(      vector&    individual,
                                       const int        length)
    {
        for (int i = 0; i < length; i++)
        {
            individual[i] = dis(gen);
        }
    }

}
