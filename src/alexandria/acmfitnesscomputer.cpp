/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#include "acmfitnesscomputer.h"

#include "act/basics/dataset.h"
#include "act/ga//Genome.h"

namespace alexandria
{

/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMFitnessComputer            *
* * * * * * * * * * * * * * * * * * * */

void ACMFitnessComputer::compute(ga::Genome *genome,
                                 iMolSelect  trgtFit,
                                 bool        verbose)
{
    if (nullptr == genome)
    {
        GMX_THROW(gmx::InternalError("Empty genome"));
    }
    sii_->updatePoldata(genome);
    double fitness = calcDeviation(genome->basesPtr(),
                                   CalcDev::Parallel, trgtFit);
    genome->setFitness(trgtFit, fitness);
    if (debug)
    {
        // TODO fix printing
        //tmpInd->printParameters(debug);
    }
    if (verbose && logfile_)
    {
        fprintf(logfile_, "Components of fitting function for %s set\n", iMolSelectName(trgtFit));
        for (const auto &ft : sii_->targets().find(trgtFit)->second)
        {
            ft.second.print(logfile_);
        }
    }
}

double ACMFitnessComputer::calcDeviation(std::vector<double> *params,
                                         CalcDev              calcDev,
                                         iMolSelect           ims)
{
    // Send / receive parameters
    std::vector<double> *myparams = nullptr;
    auto cr = sii_->commRec();
    if (cr->isHelper())
    {
        // If we have e.g. 1 overlord and 3 middlemen with 1 helper each, we have
        // O M H M H M H. Now who is my middleman? This is handled by the library
        int src = cr->superior();
        calcDev  = static_cast<CalcDev>(cr->recv_int(src));
        ims      = static_cast<iMolSelect>(cr->recv_int(src));
        myparams = new std::vector<double>;
        cr->recv_double_vector(src, myparams);
        sii_->poldata()->receiveEemprops(cr, src);
        sii_->poldata()->receiveParticles(cr, src);
    }
    else if (calcDev != CalcDev::Master)
    {
        // Send only to my helpers
        for (auto &dest : cr->helpers())
        {
            cr->send_int(dest, static_cast<int>(calcDev));
            cr->send_int(dest, static_cast<int>(ims));
            cr->send_double_vector(dest, params);
            sii_->poldata()->sendEemprops(cr, dest);
            sii_->poldata()->sendParticles(cr, dest);
        }
        myparams = params;
    }
    else
    {
        myparams = params;
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

    // If actMaster or actMiddleMan, penalize out of bounds
    if (cr->isMasterOrMiddleMan() && bdc_)
    {
        bdc_->calcDeviation(nullptr, targets, sii_->poldata(),
                            *myparams, nullptr);
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
            immStatus imm = mymol.GenerateAcmCharges(sii_->poldata());

            // Check whether we have to disable this compound
            if (immStatus::OK != imm && removeMol_)
            {
                mymol.setSupport(eSupport::No);
                continue;
            }

            computeMultipoles(targets, &mymol);

            for (DevComputer *mydev : devComputers_)
            {
                mydev->calcDeviation(&mymol, targets, sii_->poldata(),
                                     *myparams, cr);
            }
        }
    }
    // Sum the terms of the chi-squared once we have done calculations
    // for all the molecules.
    sii_->sumChiSquared(calcDev == CalcDev::Parallel, ims);

    numberCalcDevCalled_ += 1;
    
    return (*targets).find(eRMS::TOT)->second.chiSquared();
}

void ACMFitnessComputer::computeMultipoles(std::map<eRMS, FittingTarget> *targets,
                                           MyMol                         *mymol)
{
    QtypeProps *qcalc = mymol->qTypeProps(qType::Calc);
    if ((*targets).find(eRMS::MU)->second.weight() > 0   ||
        (*targets).find(eRMS::QUAD)->second.weight() > 0 ||
        (*targets).find(eRMS::OCT)->second.weight() > 0  ||
        (*targets).find(eRMS::HEXADEC)->second.weight() > 0)
    {
        qcalc->setQ(mymol->atoms());
        qcalc->setX(mymol->x());
        qcalc->calcMoments();
    }
}

void ACMFitnessComputer::fillDevComputers(const bool verbose)
{
    if (sii_->target(iMolSelect::Train, eRMS::BOUNDS)->weight() > 0)
    {
        bdc_ = new BoundsDevComputer(logfile_, verbose, sii_->optIndexPtr());
    }
    if (sii_->target(iMolSelect::Train, eRMS::CHARGE)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::CM5)->weight() > 0)
    {
        devComputers_.push_back(new ChargeCM5DevComputer(logfile_, verbose));
    }
    if (sii_->target(iMolSelect::Train, eRMS::ESP)->weight() > 0)
    {
        devComputers_.push_back(new EspDevComputer(logfile_, verbose, molgen_->fit("zeta")));
    }
    if (sii_->target(iMolSelect::Train, eRMS::Polar)->weight() > 0)
    {
        devComputers_.push_back(new PolarDevComputer(logfile_, verbose));
    }
    if (sii_->target(iMolSelect::Train, eRMS::MU)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(logfile_, verbose,
                                                         MolPropObservable::DIPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::QUAD)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(logfile_, verbose, 
                                                         MolPropObservable::QUADRUPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::OCT)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(logfile_, verbose, 
                                                         MolPropObservable::OCTUPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::HEXADEC)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(logfile_, verbose, 
                                                         MolPropObservable::HEXADECAPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::EPOT)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Force2)->weight() > 0)
    {
        devComputers_.push_back(new ForceEnergyDevComputer(logfile_, verbose));
    }
    if (sii_->target(iMolSelect::Train, eRMS::FREQUENCY)->weight() > 0)
    {
        devComputers_.push_back(new HarmonicsDevComputer(logfile_, verbose, MolPropObservable::FREQUENCY));
    }
}

/* * * * * * * * * * * * * * * * * * * *
* END: ACMFitnessComputer              *
* * * * * * * * * * * * * * * * * * * */

} // namespace alexandria
