/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#ifndef ALEXANDRIA_PERCENT_MUTATOR
#define ALEXANDRIA_PERCENT_MUTATOR

#include "ga/Mutator.h"

#include "acmindividual.h"

namespace alexandria
{


/*!
 * Changes values of the genes by a maximum of \p percent * their allowed range
 */
class PercentMutator : public ga::Mutator
{

private:
    //! IndividualInfo  pointer
    StaticIndividualInfo *sii_;
    //! The maximum change allowed as percent/100 of the range
    double                percent_;

public:

    /*!
     * Constructor
     * \param[in] sii     Pointer to StaticindividualInfo instance
     * \param[in] percent Maximum allowed change
     */
    PercentMutator(StaticIndividualInfo *sii,
                   double                percent)
    : ga::Mutator(), sii_(sii), percent_(percent) {}                

    /*! \brief Do the actual mutation
     * \param[inout] genome     The genome to mutate
     * \param[out]   bestGenome The best genome found
     * \param[in]    prMut      Probability for mutation
     */
    virtual void mutate(ga::Genome *genome,
                        ga::Genome *bestGenome,
                        double      prMut);

    virtual void finalize() {}
};


} //namespace alexandria


#endif //ALEXANDRIA_PERCENT_MUTATOR
