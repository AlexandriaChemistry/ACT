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
#include "staticindividualinfo.h"

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

    //! StaticIndividualInfo pointer
    StaticIndividualInfo                    *sii_;
    //! Whether we do random initialization or not.
    bool                                     randInit_;
    // Random number generation
    std::random_device                       rd_;
    std::mt19937                             gen_;
    std::uniform_real_distribution<double>   dis_;
    //! Seeds for random number generation
    std::vector<int>                         seeds_;
    //! Base name for Force Field output file
    std::string                              outputFile_;

public: 

    /*!
     * \brief Property constructor
     * \param[in] sii           pointer to StaticIndividualInfo instance
     * \param[in] randInit      whether we initialize the force field parameters randomly
     * \param[in] outputFile    base name for Force Field output files
     * \param[in] seed          Seed for random number initialization,
     *                          if zero it will be generated.
     */
    ACMInitializer(StaticIndividualInfo   *sii,
                   bool                    randInit,
                   const std::string      &outputFile,
                   int                     seed);

    virtual ga::Individual *initialize();

};


} // namespace alexandria


#endif // ALEXANDRIA_ACMINITIALIZER_H
