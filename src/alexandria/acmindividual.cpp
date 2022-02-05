/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#include "acmindividual.h"

#include "gromacs/fileio/xvgr.h"
#include "poldata/poldata_xml.h"

namespace alexandria
{

/* * * * * * * * * * * * * * * * * * * * * *
 * BEGIN: Output stuff                     *
 * * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::printParameters(FILE *fp) const
{
    if (nullptr == fp)
    {
        return;
    }
    for(size_t i = 0; i < genome_.nBase(); i++)
    {
        fprintf(fp, "  %s  %e,", sii_->paramNames()[i].c_str(), genome_.base(i));
    }
    fprintf(fp, "\n");
}

/* * * * * * * * * * * * * * * * * * * * * *
 * END: Output stuff                       *
 * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
 * BEGIN: Poldata stuff                    *
 * * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::copyGenome(const ga::Genome &genome)
{
    genome_ = genome;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Poldata stuff                       *
* * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::addParam(const real val)
{
    initialGenome_.addBase(val);
    genome_.addBase(val);
}

} //namespace alexandria
