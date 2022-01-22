/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#ifndef ACT_ACTHELPER_H
#define ACT_ACTHELPER_H

#include "gromacs/mdtypes/commrec.h"

namespace alexandria
{
    // Declare some classes rather than including headers.
    class StaticIndividualInfo;
    class ACMFitnessComputer;
    class MolGen;
   
    /*! \brief Class that runs the deviation calculations only
     */
    class ACTHelper
    {
    private:
        //! Fitness computer
        ACMFitnessComputer  *fitComp_;
    public:
        /*! \brief Constructor
         * \param[in] sii Static information
         * \param[in] mg  Information about this helpers molecules
         */
        ACTHelper(StaticIndividualInfo *sii,
                  MolGen               *mg);
        
        //! \brief Run the helper process
        void run();    
    };
    
} // namespace alexandria

#endif
