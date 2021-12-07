#include "acmfitnesscomputer.h"
#include "aliases.h"


namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMFitnessComputer            *
* * * * * * * * * * * * * * * * * * * */

void ACMFitnessComputer::compute(Individual *individual)
{
    // TODO: Implement this function.
}

double ACMFitnessComputer::calcDeviation(Individual   *individual,
                                         CalcDev      calcDev,
                                         iMolSelect   ims)
{
    // TODO: Implement this function.
}

void computeDiQuad(std::map<eRMS, FittingTarget> *targets,
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

void ACFMFitnessComputer::fillDevComputers()
{
    FILE *lf = logFile();
    bool verb = verbose();

    if (target(iMolSelect::Train, eRMS::BOUNDS)->weight() > 0)
        bdc_ = new BoundsDevComputer(lf, verb, &optIndex_);

    if (target(iMolSelect::Train, eRMS::CHARGE)->weight() > 0 ||
        target(iMolSelect::Train, eRMS::CM5)->weight() > 0)
        devComputers_.push_back(new ChargeCM5DevComputer(lf, verb));
    if (target(iMolSelect::Train, eRMS::ESP)->weight() > 0)
        devComputers_.push_back(new EspDevComputer(lf, verb, fit("zeta")));
    if (target(iMolSelect::Train, eRMS::Polar)->weight() > 0)
        devComputers_.push_back(new PolarDevComputer(lf, verb, bFullQuadrupole_));
    if (target(iMolSelect::Train, eRMS::QUAD)->weight() > 0)
        devComputers_.push_back(new QuadDevComputer(lf, verb, bFullQuadrupole_));
    if (target(iMolSelect::Train, eRMS::MU)->weight() > 0)
        devComputers_.push_back(new MuDevComputer(lf, verb, bQM()));
    if (target(iMolSelect::Train, eRMS::EPOT)->weight() > 0)
        devComputers_.push_back(new EnergyDevComputer(lf, verb));
}



/* * * * * * * * * * * * * * * * * * * *
* END: ACMFitnessComputer              *
* * * * * * * * * * * * * * * * * * * */


} // namespace alexandria