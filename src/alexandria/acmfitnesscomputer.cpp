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
    const iMolSelect ims = trgtFit == ga::Target::Train ? iMolSelect::Train : iMolSelect::Test;
    ACMIndividual *tmpInd = static_cast<ACMIndividual*>(ind);
    std::vector<bool> changed;
    changed.resize(tmpInd->nParam(), true);
    tmpInd->toPoldata(changed);
    const double fitness = calcDeviation(tmpInd, false, CalcDev::Parallel, ims);
    if (ims == iMolSelect::Train)
    {
        ind->setFitnessTrain(fitness);
    }
    else
    {
        ind->setFitnessTest(fitness);
    }
}

double ACMFitnessComputer::calcDeviation(ACMIndividual   *ind,
                                         bool             verbose,
                                         CalcDev          calcDev,
                                         iMolSelect       ims)
{

    // Send / receive parameters
    if (MASTER(cr_))
    {
        if (PAR(cr_) && calcDev != CalcDev::Master)
        {
            for (int i = 1; i < cr_->nnodes; i++)
            {
                gmx_send_int(cr_, i, static_cast<int>(calcDev));
                gmx_send_int(cr_, i, static_cast<int>(ims));
            }
        }
    }
    else
    {
        calcDev = static_cast<CalcDev>(gmx_recv_int(cr_, 0));
        ims     = static_cast<iMolSelect>(gmx_recv_int(cr_, 0));
    }

    // If final call, return -1
    if (calcDev == CalcDev::Final) return -1;

    // Reset the chi2 in FittingTargets of individual for the given dataset in ims
    ind->resetChiSquared(ims);

    // Gather fitting targets from the individual
    std::map<eRMS, FittingTarget> *targets = ind->sii()->fittingTargets(ims);

    // If MASTER, penalize out of bounds
    if (MASTER(cr_))
    {
        if (bdc_ != nullptr)
        {
            bdc_->calcDeviation(nullptr, targets, ind->siiConst()->poldata(), ind->param(), cr_);
        }
    }

    // If we are running in parallel, spread/receive Poldata properties
    if (PAR(cr_))
    {
        if (calcDev == CalcDev::Parallel)
        {
            ind->sii()->poldata()->broadcast_eemprop(cr_);
            ind->sii()->poldata()->broadcast_particles(cr_);
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
                    mymol.UpdateIdef(ind->sii()->poldata(), io.first);
                }
            }
            // TODO Check whether this is sufficient for updating the particleTypes
            if (molgen_->fit("zeta"))
            {
                // Update the electronegativity parameters
                mymol.zetaToAtoms(ind->sii()->poldata(), mymol.atoms());
            }
            // Run charge generation including shell minimization
            immStatus imm = mymol.GenerateAcmCharges(ind->sii()->poldata(), cr_,
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
                mydev->calcDeviation(&mymol, targets, ind->sii()->poldata(), ind->param(), cr_);
            }
        }
    }
    // Sum the terms of the chi-squared once we have done calculations
    // for all the molecules.
    ind->sumChiSquared(cr_, calcDev == CalcDev::Parallel, ims);

    if (debug)
    {
        ind->printParameters(debug);
    }

    if (verbose && logfile_)
    {
        ind->printChiSquared(cr_, logfile_, ims);  // Will only be done by MASTER
    }

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
