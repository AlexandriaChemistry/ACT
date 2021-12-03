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

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: ChargeCM5DevComputer                *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: PolarDevComputer                  *
    * * * * * * * * * * * * * * * * * * * * * */

    void PolarDevComputer::calcDeviation(      MyMol                             *mymol,
                                               std::map<eRMS, FittingTarget>     *targets,
                                               Poldata                           *poldata,
                                         const std::vector<double>               &param,
                                               t_commrec                         *commrec)
    {
        double diff2 = 0;
        mymol->CalcPolarizability(10, commrec, nullptr);
        if (bFullQuadrupole_)
        {
            // It is already squared
            diff2 = mymol->PolarizabilityTensorDeviation();
        }
        else
        {
            diff2 = gmx::square(mymol->PolarizabilityDeviation());
        }
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
    * BEGIN: QuadDevComputer                   *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: QuadDevComputer                     *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: MuDevComputer                     *
    * * * * * * * * * * * * * * * * * * * * * */

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: MuDevComputer                       *
    * * * * * * * * * * * * * * * * * * * * * */


}