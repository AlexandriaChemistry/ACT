#ifndef ALEXANDRIA_ACMINDIVIDUAL_H
#define ALEXANDRIA_ACMINDIVIDUAL_H


#include "aliases.h"
#include "ga/Initializer.h"


namespace alexandria
{


class ACMInitializer : public ga::Initializer
{

private:

    //! Minimization data
    int mindata_;
    //! Share individual info
    SharedIndividualInfo *sii_;

public: 

    /*!
        * Initialize an individual
        * @param individual pointer to the individual to initialize
        */
    virtual void initialize(ga::Individual *individual);

};


} // namespace alexandria


#endif // ALEXANDRIA_ACMINITIALIZER_H