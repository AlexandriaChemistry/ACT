/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#include "acmindividual.h"

#include "poldata_xml.h"

#include "gromacs/fileio/xvgr.h"


namespace alexandria
{
/* Constructor */
ACMIndividual::ACMIndividual(const int             id,
                             SharedIndividualInfo *sii,
                             const std::string    &outputFile)
    : ga::Individual()
{
    id_ = id;
    sii_ = sii;
    outputFile_ = "ind" + std::to_string(id_) + "/ind" + std::to_string(id_) + "-" + outputFile;
    
    // Initialize vectors for statistics and bestParam_.
    // initialParam_ and param_ will be initialized later
    size_t nParam = sii_->nParam();
    pmean_.resize(nParam, 0.0);
    psigma_.resize(nParam, 0.0);
    attemptedMoves_.resize(nParam, 0);
    acceptedMoves_.resize(nParam, 0);
    bestParam_.resize(nParam, 0.0);
    
    // Copy targets_ from sii_
    //targets_ = sii_->targets();  // This should make a deep copy if
    // the copy constructors are well made
    
    // Copy poldata from sii_
    //pd_ = sii_->poldataConst();    
}


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: File stuff                        *
* * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::fprintSelf(FILE *fp)
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "id_: %i; ", id_);
    fprintf(fp, "param_: [ ");
    for (const double ele : param_)
    {
        fprintf(fp, "%f ", ele);
    }
    fprintf(fp, "]; ");
    fprintf(fp, "fitnessTrain_: %f; fitnessTest_: %f; probability_: %f\n",
            fitnessTrain_, fitnessTest_, probability_);
}

void ACMIndividual::printHeader(FILE *fp)
{
    if (fp)
    {
        fprintf(fp, "\nIndividual %i\n", id_);
    }
}

void ACMIndividual::openParamConvFiles(const gmx_output_env_t *oenv)
{
    const std::vector<std::string> pClass = sii_->paramClass();
    for (size_t i = 0; i < pClass.size(); i++)
    {
        std::string fileName = "ind" + std::to_string(id_) + "/ind" + std::to_string(id_) + "-" + pClass[i] + "-" + sii_->xvgConv();
        fpc_.push_back(xvgropen(fileName.c_str(),
                                "Parameter convergence",
                                "iteration",
                                "",
                                oenv));

        std::vector<const char*> tmpParamNames;
        for (size_t j = 0; j < sii_->paramNames().size(); j++)
        {
            if ( (sii_->paramClassIndex())[j] == i )
            {
                tmpParamNames.push_back( (sii_->paramNames())[j].c_str() );
            }
        }
        xvgr_legend(fpc_[i], tmpParamNames.size(), tmpParamNames.data(), oenv);
    }
}

void ACMIndividual::openChi2ConvFile(const gmx_output_env_t    *oenv,
                                     const bool                 bEvaluate_testset)
{
    std::string fileName = "ind" + std::to_string(id_) + "/ind" + std::to_string(id_) + "-" + sii_->xvgEpot();
    fpe_ = xvgropen(fileName.c_str(),
                    "Chi squared",
                    "Iteration",
                    "Unknown units",
                    oenv);
    if (bEvaluate_testset)
    {
        std::vector<std::string> legend;
        legend.push_back(iMolSelectName(iMolSelect::Train));
        legend.push_back(iMolSelectName(iMolSelect::Test));
        xvgrLegend(fpe_, legend, oenv);
    }
}

void ACMIndividual::closeConvFiles()
{
    for(FILE *fp: fpc_)  // Close all parameter convergence surveillance files
    {
        xvgrclose(fp);
    }
    if (fpe_ != nullptr)  // Close chi2 surveillance file
    {
        xvgrclose(fpe_);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: File stuff                          *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Chi2 stuff                        *
* * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::printChiSquared(t_commrec *cr, FILE *fp, iMolSelect ims) const
{
    if (nullptr != fp && MASTER(cr))
    {
        fprintf(fp, "\nComponents of fitting function for %s set\n",
                iMolSelectName(ims));
        for (const auto &ft : sii_->targets().find(ims)->second)
        {
            ft.second.print(fp);
        }
        fflush(fp);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Chi2 stuff                          *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Output stuff                      *
* * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::saveState(bool updateCheckSum)
{
    sii_->poldata()->updateTimeStamp();
    if (updateCheckSum)
    {
        sii_->poldata()->updateCheckSum();
    }
    writePoldata(outputFile_, sii_->poldata(), false);
}

void ACMIndividual::printParameters(FILE *fp) const
{
    if (nullptr == fp)
    {
        return;
    }
    for(size_t i = 0; i < param_.size(); i++)
    {
        fprintf(fp, "  %s  %e,", sii_->paramNames()[i].c_str(), param_[i]);
    }
    fprintf(fp, "\n");
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Output stuff                        *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Poldata stuff                     *
* * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::toPoldata(const std::vector<bool> &changed)
{

    size_t n = 0;
    if (psigma_.empty())
    {
        psigma_.resize(param_.size(), 0);
    }
    printParameters(debug);
    for (const auto &optIndex : sii_->optIndex())
    {
        if (changed[n])
        {
            auto                 iType = optIndex.iType();
            ForceFieldParameter *p  = nullptr;
            auto                 pd = sii_->poldata(); 
            if (iType != InteractionType::CHARGE)
            {
                p = pd->findForces(iType)->findParameterType(optIndex.id(), optIndex.parameterType());
            }
            else if (pd->hasParticleType(optIndex.particleType()))
            {
                p = pd->findParticleType(optIndex.particleType())->parameter(optIndex.parameterType());
            }
            GMX_RELEASE_ASSERT(p, gmx::formatString("Could not find parameter %s", optIndex.id().id().c_str()).c_str());
            if (p)
            {
                p->setValue(param_[n]);
                p->setUncertainty(psigma_[n]);
            }
        }
        n++;
    }
    GMX_RELEASE_ASSERT(n == changed.size(),
                       gmx::formatString("n = %zu changed.size() = %zu",
                                         n, changed.size()).c_str());

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Poldata stuff                       *
* * * * * * * * * * * * * * * * * * * * * */

} //namespace alexandria
