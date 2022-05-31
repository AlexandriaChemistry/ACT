#include "forcecomputer.h"

#include "alexandria/topology.h"
#include "act/forces/forcecomputerimpl.h"
#include "gromacs/math/vec.h"

namespace alexandria
{

void ForceComputer::compute(const Poldata                     &pd,
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
                    bfc(ffpl, entry.second, coordinates, forces) });
        }
    }
}

} // namespace alexandria
