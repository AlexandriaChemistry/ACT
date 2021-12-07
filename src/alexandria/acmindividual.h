#ifndef ALEXANDRIA_ACMINDIVIDUAL_H
#define ALEXANDRIA_ACMINDIVIDUAL_H

#include <cstdio>

#include "ga/Individual.h"

namespace alexandria
{


class ACMIndividual : public ga::Individual
{

private:

    //! Whether we are in verbose
    bool verbose_;
    //! Pointer to log file
    FILE* logfile_;

public:

};


}


#endif //ALEXANDRIA_ACMINDIVIDUAL_H