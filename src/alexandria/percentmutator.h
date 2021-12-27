#ifndef ALEXANDRIA_PERCENT_MUTATOR
#define ALEXANDRIA_PERCENT_MUTATOR


#include "sharedindividualinfo.h"

#include "ga/Mutator.h"


namespace alexandria
{


/*!
 * Changes values of the genes by a maximum of \p percent * their allowed range
 */
class PercentMutator : public ga::Mutator
{

private:

    //! SharedIndividualInfo pointer
    SharedIndividualInfo *sii_;
    //! The maximum change allowed as percent/100 of the range
    double percent_;

public:

    /*!
     * Constructor
     * \param[in] sii   pointer to SharedIndividualInfo instance
     */
    PercentMutator(      SharedIndividualInfo  *sii,
                   const double                 percent)
    : ga::Mutator(), sii_(sii), percent_(percent) {}                

    virtual void mutate(      ga::Individual   *ind,
                        const double            prMut);

};


} //namespace alexandria


#endif //ALEXANDRIA_PERCENT_MUTATOR