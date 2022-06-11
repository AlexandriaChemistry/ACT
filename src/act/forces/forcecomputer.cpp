#include "forcecomputer.h"

#include "alexandria/topology.h"
#include "act/forces/forcecomputerimpl.h"
#include "gromacs/math/vec.h"
#include "act/qgen/qtype.h"

namespace alexandria
{

static double dotProdRvec(const std::vector<bool>      &isShell,
                          const std::vector<gmx::RVec> &rv)
{
    double dpr = 0;
    int    i   = 0;
    for(const auto &rr : rv)
    {
        if (isShell[i++])
        {
            dpr += iprod(rr, rr);
        }
    }
    return dpr;
}

double ForceComputer::compute(const Topology                    *top,
                              std::vector<gmx::RVec>            *coordinates,
                              std::vector<gmx::RVec>            *forces,
                              std::map<InteractionType, double> *energies) const
{
    gmx::RVec field = { 0, 0, 0 };
    // Do first calculation every time.
    computeOnce(top, coordinates, forces, energies, field);
    // Now let's have a look whether we are polarizable
    auto itype = InteractionType::POLARIZATION;
    if (!pd_->polarizable() || !top->hasEntry(itype))
    {
        return 0;
    }
    // Is this particle a shell?
    std::vector<bool>   isShell;
    // One over force constant for this particle
    std::vector<double> fcShell_1;
    auto ffpl = pd_->findForcesConst(itype);
    for(auto &aa : top->atoms())
    {
        bool bIS = aa.pType() == eptShell;
        isShell.push_back(bIS);
        double fc_1 = 0;
        if (bIS)
        {
            Identifier atID(aa.ffType());
            auto fc = ffpl.findParameterTypeConst(atID, "kshell").internalValue();
            if (fc != 0)
            {
                fc_1 = 1.0/fc;
            }
        }
        fcShell_1.push_back(fc_1);
    }
    double msForce  = dotProdRvec(isShell, *forces);
    auto pols       = top->entry(itype);
    // TODO pass the real tolerance
    double toler    = 0.00000001;
    // TODO pass the real maxiter
    int    maxiter  = 25;
    int    iter     = 1;
    // Golden ratio, may be used for overrelaxation
    // double gold     = 0.5*(1+std::sqrt(5.0));
    while (msForce > toler && iter < maxiter)
    {
        // Loop over polarizabilities
        for(const auto &p : pols)
        {
            // Displace the shells according to the force
            // Since the potential is harmonic we use Hooke's law
            // F = k dx -> dx = F / k
            // TODO Optimize this protocol using overrelaxation
            int shell = p->atomIndex(1);
            for(int m = 0; m < DIM; m++)
            {
                (*coordinates)[shell][m] += (*forces)[shell][m] * fcShell_1[shell];
            }
        }
        // Do next calculation
        computeOnce(top, coordinates, forces, energies, field);
        msForce  = dotProdRvec(isShell, *forces);
        iter    += 1;
    }
    return std::sqrt(msForce);
}

void ForceComputer::computeOnce(const Topology                    *top,
                                std::vector<gmx::RVec>            *coordinates,
                                std::vector<gmx::RVec>            *forces,
                                std::map<InteractionType, double> *energies,
                                const gmx::RVec                   &field) const
{
    // Clear energies
    energies->clear();
    // Clear forces
    auto atoms = top->atoms();
    int i = 0;
    for(size_t ff = 0; ff < forces->size(); ++ff)
    {
        for(int m = 0; m < DIM; m++)
        {
            (*forces)[ff][m] = field[m]*atoms[i].charge();
        }
        i += 1;
    }
    double epot = 0;
    for(const auto &entry : top->entries())
    {
        if (entry.second.empty())
        {
            continue;
        }
        // Force field parameter list
        auto ffpl  = pd_->findForcesConst(entry.first);
        // The function we need to do the math
        auto bfc   = getBondForceComputer(ffpl.fType());
        if (nullptr == bfc)
        {
            fprintf(stderr, "Please implement a force function for type %s\n", interaction_function[ffpl.fType()].name);
        }
        else
        {
            // Now do the calculations and store the energy
            double eee = bfc(ffpl, entry.second, top->atoms(), coordinates, forces);
            energies->insert({ entry.first, eee });
            epot += eee;
        }
    }
    energies->insert({ InteractionType::EPOT, epot });
}

void ForceComputer::calcPolarizability(const Topology         *top,
                                       std::vector<gmx::RVec> *coordinates,
                                       QtypeProps             *qtp) const
{
    std::vector<gmx::RVec>            forces(coordinates->size());
    std::map<InteractionType, double> energies;
    gmx::RVec  field = { 0, 0, 0 };
    
    computeOnce(top, coordinates, &forces, &energies, field);
    std::vector<double> q;
    for (auto &at : top->atoms())
    {
        q.push_back(at.charge());
    }
    qtp->setQ(q);
    qtp->setX(*coordinates);
    qtp->calcMoments();
    auto mpo = MolPropObservable::DIPOLE;
    if (!qtp->hasMultipole(mpo))
    {
        GMX_THROW(gmx::InternalError("No dipole to compute."));
    }
    auto mu_ref = qtp->getMultipole(mpo);
    // Convert from e nm2/V to cubic nm
    double enm2_V = E_CHARGE*1e6*1e-18/(4*M_PI*EPSILON0_SI)*1e21;
    tensor alpha  = { { 0 } };
    double efield = 0.1;
    for (auto m = 0; m < DIM; m++)
    {
        field[m] = efield;
        computeOnce(top, coordinates, &forces, &energies, field);
        qtp->setX(*coordinates);
        field[m] = 0;
        qtp->calcMoments();
        auto qmu = qtp->getMultipole(mpo);
        for (auto n = 0; n < DIM; n++)
        {
            alpha[n][m] = enm2_V*((qmu[n]-mu_ref[n])/efield);
        }
    }
    // Store the tensor
    qtp->setPolarizabilityTensor(alpha);
    // Reset energies etc.
    computeOnce(top, coordinates, &forces, &energies, field);
}

int ForceComputer::ftype(InteractionType itype) const
{
    int ftype = F_EPOT;
    if (pd_->interactionPresent(itype))
    {
        ftype = pd_->findForcesConst(itype).fType();
    }
    return ftype;
}

} // namespace alexandria
