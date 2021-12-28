/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#include "percentmutator.h"

#include "acmindividual.h"


namespace alexandria
{


void PercentMutator::mutate(      ga::Individual   *ind,
                            const double            prMut)
{

    ACMIndividual *tmpInd = static_cast<ACMIndividual*>(ind);
    const std::vector<double> oldParam = tmpInd->param();
    const std::vector<double> lb = sii_->lowerBound();
    const std::vector<double> ub = sii_->upperBound();
    double newVal;

    for (size_t i = 0; i < tmpInd->nParam(); i++)
    {
        if (randNum() <= prMut)
        {
            newVal = oldParam[i] + percent_*(2*randNum()-1)*(ub[i] - lb[i]);
            newVal = std::max(lb[i], std::min(ub[i], newVal));
            tmpInd->setParam(i, newVal);
        }
    }

}


} //namespace alexandria