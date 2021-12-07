#ifndef ALEXANDRIA_ACMINITIALIZER_H
#define ALEXANDRIA_ACMINITIALIZER_H


#include "aliases.h"
#include "ga/Initializer.h"
#include "sharedindividualinfo.h"

#include <time.h>
#include <random>


namespace alexandria
{


class ACMInitializer : public ga::Initializer
{

private:

    //! The minimum amount of data points to consider a parameter
    int mindata_;
    //! Share individual info
    SharedIndividualInfo *sii_;
    //! Random initialization or not.
    bool randInit_;
    // Random number generation
    std::random_device                      rd;
    std::mt19937                            gen;
    std::uniform_real_distribution<double>  dis;
    //! Amount of initialized individuals
    int nCreated_ = 0;
    //! Base name
    std::string outputFile_;

public: 

    ACMInitializer(const int                     mindata,
                         SharedIndividualInfo   *sii,
                   const bool                    randInit,
                   const std::string            &outputFile)
    : gen(rd()), dis(std::uniform_real_distribution<>(0.0, 1.0))
    {
        mindata_  = mindata;
        sii_      = sii;
        randInit_ = randInit;
        outputFile_ = outputFile;
    }

    /*!
        * Initialize an individual
        * @param individual pointer to pointer to the individual to initialize
        */
    virtual void initialize(ga::Individual **individual);



};


} // namespace alexandria


#endif // ALEXANDRIA_ACMINITIALIZER_H