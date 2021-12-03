#include "devcomputer.h"

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


void BoundsDevComputer::calcDeviation(      MyMol                             *mymol,
                                            std::map<eRMS, FittingTarget>     *targets,
                                            Poldata                           *poldata,
                                      const std::vector<double>               &param,
                                            t_commrec                         *commrec)
{
    double bound = 0;
    size_t n     = 0;
    for (auto &optIndex : optIndex_)
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

void ChargeCM5DevComputer::calcDeviation(      MyMol                             *mymol,
                                               std::map<eRMS, FittingTarget>     *targets,
                                               Poldata                           *poldata,
                                         const std::vector<double>               &param,
                                               t_commrec                         *commrec)
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
        if (nullptr != mymol->shellfc_ &&
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


} // namespace alexandria