/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#include "devcomputer.h"

#include "act/utility/communicationrecord.h"
#include "act/utility/units.h"

namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: BoundsDevComputer                 *
* * * * * * * * * * * * * * * * * * * * * */

double BoundsDevComputer::l2_regularizer(      double          x,
                                               double          min,
                                               double          max,
                                         const std::string    &label)
{

    double p = 0;
    if (x < min)
    {
        p = (0.5 * gmx::square(x-min));
    }
    else if (x > max)
    {
        p = (0.5 * gmx::square(x-max));
    }
    if (verbose_ && p != 0.0)
    {
        fprintf(logfile_, "Variable %s is %g, should be within %g and %g\n", label.c_str(), x, min, max);
    }
    return p;

}


void BoundsDevComputer::calcDeviation(gmx_unused MyMol                         *mymol,
                                      gmx_unused std::map<eRMS, FittingTarget> *targets,
                                            Poldata                             *poldata,
                                      const std::vector<double>                 &param,
                                      gmx_unused const CommunicationRecord      *commrec)
{
    double bound = 0;
    size_t n     = 0;
    for (auto &optIndex : *optIndex_)
    {
        InteractionType iType = optIndex.iType();
        ForceFieldParameter p;
        if (iType == InteractionType::CHARGE)
        {
            p = poldata->findParticleType(optIndex.particleType())->parameterConst(optIndex.parameterType());
        }
        else if (poldata->interactionPresent(iType))
        {
            p = poldata->findForcesConst(iType).findParameterTypeConst(optIndex.id(), optIndex.parameterType());
        }
        if (p.mutability() == Mutability::Bounded)
        {
            bound += l2_regularizer(param[n], p.minimum(), p.maximum(), optIndex.name());
        }
        n++;
    }
    (*targets).find(eRMS::BOUNDS)->second.increase(1, bound);
    GMX_RELEASE_ASSERT(n == param.size(),
                        gmx::formatString("Death horror error. n=%zu param.size()=%zu", n, param.size()).c_str());
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: BoundsDevComputer                   *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: ChargeCM5DevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

void ChargeCM5DevComputer::calcDeviation(MyMol                                *mymol,
                                         std::map<eRMS, FittingTarget>        *targets,
                                         Poldata                              *poldata,
                                         gmx_unused const std::vector<double> &param,
                                         gmx_unused const CommunicationRecord *commrec)
{
    double qtot = 0;
    int i = 0;
    const t_atoms myatoms = mymol->atomsConst();
    std::vector<double> qcm5;
    QtypeProps *qp = mymol->qTypeProps(qType::CM5);
    if (qp)
    {
        qcm5 = qp->charge();
        if (debug)
        {
            for (int j = 0; j < myatoms.nr; j++)
            {
                fprintf(debug, "Charge %d. CM5 = %g ACM = %g\n", j, qcm5[j], myatoms.atom[j].q);
            }
        }
    }
    // Iterate over the atoms
    for (int j = 0; j < myatoms.nr; j++)
    {
        if (myatoms.atom[j].ptype == eptShell)
        {
            continue;
        }
        ParticleTypeIterator atype = poldata->findParticleType(*myatoms.atomtype[j]);
        const ForceFieldParameter qparm = atype->parameterConst("charge");
        double qj  = myatoms.atom[j].q;
        double qjj = qj;
        // TODO: only count in real shells
        if (mymol->haveShells() &&
            j < myatoms.nr-1 &&
            myatoms.atom[j+1].ptype == eptShell)
        {
            qjj += myatoms.atom[j+1].q;
        }
        qtot += qjj;
        switch (qparm.mutability())
        {
        case Mutability::Fixed:
            if (qparm.value() != qj)
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Fixed charge for atom %s in %s was changed from %g to %g",
                                                                *myatoms.atomname[j], mymol->getMolname().c_str(), qparm.value(), qj).c_str()));
            }
            break;
        case Mutability::Bounded:
            {
                if ((*targets).find(eRMS::CHARGE)->second.weight() > 0)
                {
                    real dq = 0;
                    if (qj < qparm.minimum())
                    {
                        dq = qparm.minimum() - qj;
                    }
                    else if (qj > qparm.maximum())
                    {
                        dq = qj - qparm.maximum();
                    }
                    (*targets).find(eRMS::CHARGE)->second.increase(1, dq*dq);
                }
            }
            break;
        default:
            break;
        }
        if (qp &&
            qparm.mutability() != Mutability::Fixed &&
            (*targets).find(eRMS::CM5)->second.weight() > 0)
        {
            // TODO: Add charge of shell!
            real dq2 = gmx::square(qjj - qcm5[i]);
            (*targets).find(eRMS::CM5)->second.increase(1, dq2);
        }
        i += 1;
    }
    (*targets).find(eRMS::CHARGE)->second.increase(1, gmx::square(qtot - mymol->totalCharge()));

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: ChargeCM5DevComputer                *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: EspDevComputer                    *
* * * * * * * * * * * * * * * * * * * * * */

void EspDevComputer::calcDeviation(MyMol                                *mymol,
                                   std::map<eRMS, FittingTarget>        *targets,
                                   Poldata                              *poldata,
                                   gmx_unused const std::vector<double> &param,
                                   gmx_unused const CommunicationRecord *commrec)
{

    real rrms     = 0;
    real cosangle = 0;
    QgenResp *qgr = mymol->qTypeProps(qType::Calc)->qgenResp();
    if (mymol->haveShells())
    {
        qgr->updateAtomCoords(mymol->x());
    }
    if (fit_)
    {
        qgr->updateZeta(mymol->atoms(), poldata);
    }
    dumpQX(logfile_, mymol, "ESP");
    qgr->updateAtomCharges(mymol->atoms());
    qgr->calcPot(poldata->getEpsilonR());
    real mae, mse;
    real rms = qgr->getStatistics(&rrms, &cosangle, &mae, &mse);
    double myRms = convertToGromacs(rms, "Hartree/e");
    size_t nEsp = qgr->nEsp();
    (*targets).find(eRMS::ESP)->second.increase(nEsp, gmx::square(myRms)*nEsp);
    if (debug)
    {
        fprintf(debug, "%s ESPrms = %g cosangle = %g\n",
                mymol->getMolname().c_str(),
                myRms, cosangle);
    }

}

void EspDevComputer::dumpQX(FILE *fp, MyMol *mol, const std::string &info)
{
    if (false && fp)
    {
        std::string label = mol->getMolname() + "-" + info;
        fprintf(fp, "%s q:", label.c_str());
        t_mdatoms *md = mol->getMdatoms();
        auto myatoms = mol->atomsConst();
        for (int i = 0; i < myatoms.nr; i++)
        {
            fprintf(fp, " %g (%g)", myatoms.atom[i].q,
                    md->chargeA[i]);
        }
        fprintf(fp, "\n");
        auto top = mol->topology();
        if (top->hasEntry(InteractionType::POLARIZATION))
        {
            fprintf(fp, "%s alpha", label.c_str());
            for(auto &topentry : top->entry(InteractionType::POLARIZATION))
            {
                //fprintf(fp, " %g", mol->ltop_->idef.iparams[topentry->gromacsType()].polarize.alpha);
            }
            fprintf(fp, "\n");
        }
        pr_rvecs(fp, 0, label.c_str(), mol->x().rvec_array(), myatoms.nr);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: EspDevComputer                      *
* * * * * * * * * * * * * * * * * * * * * */


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: PolarDevComputer                  *
* * * * * * * * * * * * * * * * * * * * * */

void PolarDevComputer::calcDeviation(MyMol                                *mymol,
                                     std::map<eRMS, FittingTarget>        *targets,
                                     gmx_unused Poldata                   *poldata,
                                     gmx_unused const std::vector<double> &param,
                                     gmx_unused const CommunicationRecord *commrec)
{
    double diff2 = 0;
    mymol->CalcPolarizability(10, nullptr);
    diff2 = gmx::square(mymol->PolarizabilityDeviation());
    
    if (false && logfile_)
    {
        fprintf(logfile_, "DIFF %s %g\n", mymol->getMolname().c_str(), diff2);
    }
    (*targets).find(eRMS::Polar)->second.increase(1, diff2);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: PolarDevComputer                    *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: MultiPoleDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

void MultiPoleDevComputer::calcDeviation(MyMol                                *mymol,
                                         std::map<eRMS, FittingTarget>        *targets,
                                         gmx_unused Poldata                   *poldata,
                                         gmx_unused const std::vector<double> &param,
                                         gmx_unused const CommunicationRecord *commrec)
{
    if (!(mymol->qTypeProps(qType::Elec)->hasMultipole(mpo_) &&
          mymol->qTypeProps(qType::Calc)->hasMultipole(mpo_)))
    {
        return;
    }
    auto qelec = mymol->qTypeProps(qType::Elec)->getMultipole(mpo_);
    auto qcalc = mymol->qTypeProps(qType::Calc)->getMultipole(mpo_);
    double delta = 0;
    for (size_t mm = 0; mm < qelec.size(); mm++)
    {
        delta += gmx::square(qcalc[mm] - qelec[mm]);
    }
    eRMS rms;
    switch (mpo_)
    {
    case MolPropObservable::DIPOLE:
        rms = eRMS::MU;
        break;
    case MolPropObservable::QUADRUPOLE:
        rms = eRMS::QUAD;
        break;
    case MolPropObservable::OCTUPOLE:
        rms = eRMS::OCT;
        break;
    case MolPropObservable::HEXADECAPOLE:
        rms = eRMS::HEXADEC;
        break;
    default:
        GMX_THROW(gmx::InternalError(gmx::formatString("Not support MolPropObservable %s", mpo_name(mpo_)).c_str()));
    }
    
    (*targets).find(rms)->second.increase(1, delta);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: MultipoleDevComputer                *
* * * * * * * * * * * * * * * * * * * * * */


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: EnergyDevComputer                 *
* * * * * * * * * * * * * * * * * * * * * */

void EnergyDevComputer::calcDeviation(MyMol                                *mymol,
                                      std::map<eRMS, FittingTarget>        *targets,
                                      gmx_unused Poldata                   *poldata,
                                      gmx_unused const std::vector<double> &param,
                                      gmx_unused const CommunicationRecord *commrec)
{
    double emol;
    GMX_RELEASE_ASSERT(mymol->energy(MolPropObservable::EMOL, &emol),
                       gmx::formatString("No molecular energy for %s", mymol->getMolname().c_str()).c_str());
    
    real delta = gmx::square(mymol->potentialEnergy() - emol);
    (*targets).find(eRMS::EPOT)->second.increase(1, delta);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: EnergyDevComputer                       *
* * * * * * * * * * * * * * * * * * * * * */


} // namespace alexandria
