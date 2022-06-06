#include "forcecomputer.h"

#include "alexandria/topology.h"
#include "act/forces/forcecomputerimpl.h"
#include "gromacs/math/vec.h"

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

void ForceComputer::compute(const Poldata                     &pd,
                            const Topology                    &top,
                            const std::vector<double>         &charge,
                            std::vector<gmx::RVec>            *coordinates,
                            std::vector<gmx::RVec>            *forces,
                            std::map<InteractionType, double> *energies)
{
    // Do first calculation every time.
    computeOnce(pd, top, coordinates, forces, energies);
    // Now let's have a look whether we are polarizable
    auto itype = InteractionType::POLARIZATION;
    if (!pd.polarizable() || !top.hasEntry(itype))
    {
        return;
    }
    // Is this particle a shell?
    std::vector<bool>   isShell;
    // One over force constant for this particle
    std::vector<double> fcShell_1;
    auto ffpl = pd.findForcesConst(itype);
    int i = 0;
    for(auto &aa : top.atoms())
    {
        bool bIS = aa.pType() == eptShell;
        isShell.push_back(bIS);
        double fc_1 = 0;
        if (bIS)
        {
            Identifier atID({aa.ffType()}, {}, CanSwap::No);
            auto fc = ffpl.findParameterTypeConst(atID, "kshell").value();
            if (fc != 0)
            {
                fc_1 = 1.0/fc;
            }
        }
        fcShell_1.push_back(fc_1);
        i += 1;
    }
    double rmsForce = dotProdRvec(isShell, *forces);
    auto pols       = top.entry(itype);
    // TODO pass the real tolerance
    double toler    = 0.001;
    // TODO pass the real maxiter
    int    maxiter  = 25;
    int    iter     = 1;
    // Golden ratio
    double gold     = 0.5*(1+std::sqrt(5.0));
    while (rmsForce > toler && iter < maxiter)
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
        computeOnce(pd, top, coordinates, forces, energies);
        rmsForce  = dotProdRvec(isShell, *forces);
        iter     += 1;
    }
#undef next
}

void ForceComputer::computeOnce(const Poldata                     &pd,
                                const Topology                    &top,
                                std::vector<gmx::RVec>            *coordinates,
                                std::vector<gmx::RVec>            *forces,
                                std::map<InteractionType, double> *energies)
{
    // Clear energies
    energies->clear();
    // Clear forces
    for(auto ff = forces->begin(); ff < forces->end(); ++ff)
    {
        clear_rvec(*ff);
    }
    for(const auto &entry : top.entries())
    {
        if (entry.second.empty())
        {
            continue;
        }
        // Force field parameter list
        auto ffpl  = pd.findForcesConst(entry.first);
        // The function we need to do the math
        auto bfc   = getBondForceComputer(ffpl.fType());
        if (nullptr == bfc)
        {
            fprintf(stderr, "Please implement a force function for type %s\n", interaction_function[ffpl.fType()].name);
        }
        else
        {
            // Now do the calculations and store the energy
            energies->insert({ entry.first,
                    bfc(ffpl, entry.second, top.atoms(), coordinates, forces) });
        }
    }
}

} // namespace alexandria
