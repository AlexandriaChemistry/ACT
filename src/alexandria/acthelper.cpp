#include "acthelper.h"

#include <vector>
    
#include "acmfitnesscomputer.h"
#include "bayes.h"
#include "act/basics/dataset.h"

namespace alexandria
{
    
ACTHelper::ACTHelper(StaticIndividualInfo *sii,
                     MolGen               *mg)
{
    fitComp_ = new ACMFitnessComputer(nullptr, sii, mg, 
                                      false, false, false);
}

void ACTHelper::run()
{
    // H E L P E R   N O D E
    // All variables are set by the master or middlemen, but
    // we have to pass something.
    // If the result is less than zero (-1), we are done.
    std::vector<double> dummy;
    while (fitComp_->calcDeviation(&dummy, CalcDev::Parallel, iMolSelect::Train) >= 0)
    {
        ;
    }
}

} // namespace alexandria
