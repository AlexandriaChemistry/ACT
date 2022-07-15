/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#include "acmfitnesscomputer.h"

#include "act/basics/dataset.h"
#include "act/forces/combinationrules.h"
#include "act/ga/Genome.h"

namespace alexandria
{

/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMFitnessComputer            *
* * * * * * * * * * * * * * * * * * * */

void ACMFitnessComputer::compute(ga::Genome    *genome,
                                 iMolSelect     trgtFit,
                                 bool           verbose)
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
        bdc_->calcDeviation(forceComp_, nullptr, targets, sii_->poldata(),
                            *myparams, nullptr);
    }

    generateDependentParameter(sii_->poldata());

    // Loop over molecules
    for (MyMol &mymol : molgen_->mymols())
    {
        if (ims != mymol.datasetType())
        {
            continue;
        }
        if ((mymol.support() == eSupport::Local) ||
            (calcDev == CalcDev::Master && mymol.support() == eSupport::Remote))
        {
            std::vector<InteractionType> itUpdate;
            for(auto &io : molgen_->iopt())
            {
                if (io.second)
                {
                    itUpdate.push_back(io.first);
                }
            }
            // Update the polarizabilities and other params only once before the loop
            mymol.UpdateIdef(sii_->poldata(), itUpdate, 
                             molgen_->fit("zeta"));
            // Run charge generation including shell minimization
            immStatus imm = mymol.GenerateAcmCharges(sii_->poldata(), forceComp_);

            // Check whether we have to disable this compound
            if (immStatus::OK != imm && removeMol_)
            {
                mymol.setSupport(eSupport::No);
                continue;
            }

            computeMultipoles(targets, &mymol);

            for (DevComputer *mydev : devComputers_)
            {
                mydev->calcDeviation(forceComp_, &mymol, targets, sii_->poldata(),
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
        qcalc->setQ(mymol->atomsConst());
        qcalc->setX(mymol->x());
        qcalc->calcMoments();
    }
}

void ACMFitnessComputer::fillDevComputers(const bool verbose)
{
    if (sii_->target(iMolSelect::Train, eRMS::BOUNDS)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::UNPHYSICAL)->weight() > 0)
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
