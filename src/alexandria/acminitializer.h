/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 * \author Julian Ramon Marrades Furquet <julianramon.marradesfurquet.8049@student.uu.se>
 */


#ifndef ALEXANDRIA_ACMINITIALIZER_H
#define ALEXANDRIA_ACMINITIALIZER_H


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
    //! Base name plus the inividual part
    std::string outputFile_;

public: 

    ACMInitializer(const int                     mindata,
                         SharedIndividualInfo   *sii,
                   const bool                    randInit,
                   const std::string            &outputFile)
    : gen(rd()), dis(std::uniform_real_distribution<double>(0.0, 1.0))
    {

        gen.seed(::time(NULL));

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