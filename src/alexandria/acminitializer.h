/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef ALEXANDRIA_ACMINITIALIZER_H
#define ALEXANDRIA_ACMINITIALIZER_H


#include "ga/Initializer.h"
#include "sharedindividualinfo.h"

#include <time.h>
#include <random>


namespace alexandria
{

/*!
 * \brief Initializes Individual instances as ACMIndividual objects
 */
class ACMInitializer : public ga::Initializer
{

private:

    //! SharedIndividualInfo pointer
    SharedIndividualInfo *sii_;
    //! Whether we do random initialization or not.
    bool randInit_;
    // Random number generation
    std::random_device                      rd;
    std::mt19937                            gen;
    std::uniform_real_distribution<double>  dis;
    //! Amount of initialized individuals
    int nCreated_ = 0;
    //! Base name for Force Field output file
    std::string outputFile_;

public: 

    /*!
     * \brief Property constructor
     * \param[in] sii           pointer to SharedIndividualInfo instance
     * \param[in] randInit      whether we initialize the force field parameters randomly
     * \param[in] outputFile    base name for Force Field output files
     */
    ACMInitializer(      SharedIndividualInfo   *sii,
                   const bool                    randInit,
                   const std::string            &outputFile)
    : gen(rd()), dis(std::uniform_real_distribution<double>(0.0, 1.0))
    {

        gen.seed(::time(NULL));

        sii_      = sii;
        randInit_ = randInit;
        outputFile_ = outputFile;
    }

    virtual void initialize(ga::Individual **ind);

};


} // namespace alexandria


#endif // ALEXANDRIA_ACMINITIALIZER_H