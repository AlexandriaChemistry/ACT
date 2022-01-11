/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#include "acmfitnesscomputer.h"


namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMFitnessComputer            *
* * * * * * * * * * * * * * * * * * * */

void ACMFitnessComputer::compute(ga::Individual *ind,
                                 ga::Target      trgtFit)
{
    if (nullptr == ind)
    {
        GMX_THROW(gmx::InternalError("Empty individual"));
    }
    const iMolSelect ims = trgtFit == ga::Target::Train ? iMolSelect::Train : iMolSelect::Test;
    ACMIndividual *tmpInd = static_cast<ACMIndividual*>(ind);
    std::vector<bool> changed;
    // TODO: the middle man should know/decide which parameters have changed.
    changed.resize(sii_->nParam(), true);
    tmpInd->toPoldata(changed);
    const double fitness = calcDeviation(tmpInd->paramPtr(), CalcDev::Parallel, ims);
    if (ims == iMolSelect::Train)
    {
        tmpInd->setFitnessTrain(fitness);
    }
    else
    {
        tmpInd->setFitnessTest(fitness);
    }
    if (debug)
    {
        tmpInd->printParameters(debug);
    }
    if (verbose_ && logfile_)
    {
        tmpInd->printChiSquared(cr_, logfile_, ims);  // Will only be done by MASTER
    }
}

double ACMFitnessComputer::calcDeviation(std::vector<double> *params,
                                         CalcDev              calcDev,
                                         iMolSelect           ims)
{
    // Send / receive parameters
    std::vector<double> *myparams;
    if (MIDDLEMAN(cr_))
    {
        if (PAR(cr_) && calcDev != CalcDev::Master)
        {
            // Send only to my helpers
            for (int dest = cr_->nodeid+1; dest < cr_->nodeid+cr_->nhelper_per_middleman; dest++)
            {
                gmx_send_int(cr_, dest, static_cast<int>(calcDev));
                gmx_send_int(cr_, dest, static_cast<int>(ims));
                gmx_send_double_vector(cr_, dest, params);
            }
        }
        myparams = params;
    }
    else
    {
        // If we have e.g. 3 middlemen with 1 helper each, we have
        // M H M H M H. Now who is my middleman?
        // Do integer division, rounding down, the multiply again.
        int src = (cr_->nodeid / cr_->nhelper_per_middleman) * cr_->nhelper_per_middleman;
        calcDev  = static_cast<CalcDev>(gmx_recv_int(cr_, src));
        ims      = static_cast<iMolSelect>(gmx_recv_int(cr_, src));
        myparams = new std::vector<double>;
        gmx_recv_double_vector(cr_, src, myparams);
    }

    // If final call, return -1
    if (calcDev == CalcDev::Final)
    {
        return -1;
    }

    // Reset the chi2 in FittingTargets for the given dataset in ims
    sii_->resetChiSquared(ims);

    // Gather fitting targets
    std::map<eRMS, FittingTarget> *targets = sii_->fittingTargets(ims);

    // If MIDDLEMAN, penalize out of bounds
    if (MIDDLEMAN(cr_) && bdc_)
    {
        bdc_->calcDeviation(nullptr, targets, sii_->poldata(), *myparams, cr_);
    }

    // If we are running in parallel, spread/receive Poldata properties
    if (PAR(cr_))
    {
        if (calcDev == CalcDev::Parallel)
        {
            sii_->poldata()->broadcast_eemprop(cr_);
            sii_->poldata()->broadcast_particles(cr_);
        }
    }

    // Loop over molecules
    int nmolCalculated = 0;
    for (MyMol &mymol : molgen_->mymols())
    {
        if (ims != mymol.datasetType())
        {
            continue;
        }
        if ((mymol.support() == eSupport::Local) ||
            (calcDev == CalcDev::Master && mymol.support() == eSupport::Remote))
        {
            nmolCalculated += 1;
            for(auto &io : molgen_->iopt())
            {
                if (io.second)
                {
                    // Update the polarizabilities only once before the loop
                    mymol.UpdateIdef(sii_->poldata(), io.first);
                }
            }
            // TODO Check whether this is sufficient for updating the particleTypes
            if (molgen_->fit("zeta"))
            {
                // Update the electronegativity parameters
                mymol.zetaToAtoms(sii_->poldata(), mymol.atoms());
            }
            // Run charge generation including shell minimization
            immStatus imm = mymol.GenerateAcmCharges(sii_->poldata(), cr_,
                                                     molgen_->qcycle(), molgen_->qtol());

            // Check whether we have to disable this compound
            if (immStatus::OK != imm && removeMol_)
            {
                mymol.setSupport(eSupport::No);
                continue;
            }

            computeDiQuad(targets, &mymol);

            for (DevComputer *mydev : devComputers_)
            {
                mydev->calcDeviation(&mymol, targets, sii_->poldata(), *myparams, cr_);
            }
        }
    }
    // Sum the terms of the chi-squared once we have done calculations
    // for all the molecules.
    sii_->sumChiSquared(cr_, calcDev == CalcDev::Parallel, ims);

    numberCalcDevCalled_ += 1;
    
    return (*targets).find(eRMS::TOT)->second.chiSquared();
}

void ACMFitnessComputer::computeDiQuad(std::map<eRMS, FittingTarget> *targets,
                                       MyMol                         *mymol)
{
    QtypeProps *qcalc = mymol->qTypeProps(qType::Calc);
    if ((*targets).find(eRMS::MU)->second.weight() > 0 ||
        (*targets).find(eRMS::QUAD)->second.weight() > 0)
    {
        qcalc->setQ(mymol->atoms());
        qcalc->setX(mymol->x());
        qcalc->calcMoments();
    }
}

void ACMFitnessComputer::fillDevComputers()
{
    if (sii_->target(iMolSelect::Train, eRMS::BOUNDS)->weight() > 0)
        bdc_ = new BoundsDevComputer(logfile_, verbose_, sii_->optIndexPtr());

    if (sii_->target(iMolSelect::Train, eRMS::CHARGE)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::CM5)->weight() > 0)
        devComputers_.push_back(new ChargeCM5DevComputer(logfile_, verbose_));
    if (sii_->target(iMolSelect::Train, eRMS::ESP)->weight() > 0)
        devComputers_.push_back(new EspDevComputer(logfile_, verbose_, molgen_->fit("zeta")));
    if (sii_->target(iMolSelect::Train, eRMS::Polar)->weight() > 0)
        devComputers_.push_back(new PolarDevComputer(logfile_, verbose_, fullQuadrupole_));
    if (sii_->target(iMolSelect::Train, eRMS::QUAD)->weight() > 0)
        devComputers_.push_back(new QuadDevComputer(logfile_, verbose_, fullQuadrupole_));
    if (sii_->target(iMolSelect::Train, eRMS::MU)->weight() > 0)
        devComputers_.push_back(new MuDevComputer(logfile_, verbose_, molgen_->bQM()));
    if (sii_->target(iMolSelect::Train, eRMS::EPOT)->weight() > 0)
        devComputers_.push_back(new EnergyDevComputer(logfile_, verbose_));
}



/* * * * * * * * * * * * * * * * * * * *
* END: ACMFitnessComputer              *
* * * * * * * * * * * * * * * * * * * */


} // namespace alexandria
