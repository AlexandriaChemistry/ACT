#include "alexandria/topology.h"
    
namespace alexandria
{

void ForceComputer::makePairList(const Topology &top)
{

}

ForceComputer::ForceComputer(const Topology &top,
                             const Poldata  &pd)
{
    makePairList(top);
}
    
double computeBonds(const std::vector<TopologyEntry *> &bonds,
                    const std::vector<gmx::RVec>       *coordinates,
                    std::vector<gmx::RVec>             *forces)
{
}

void ForceComputer::compute(std::vector<gmx::RVec>           *coordinates,
                            std::vector<gmx::RVec>           *forces,
                            std::ma<InteractionType, double> *energies)
{
    // Clear energies
    energies.clear();
    // Clear forces
    for(auto ff = forces->begin(); ff < forces->end(); ++ff)
    {
        clear_rvec(*ff);
    }
    for(const auto &entry : top.entries())
    {
        switch(entry.first)
        {
        case InteractionType::BONDS:
            energies->insert({ entry.first,
                    computeBonds(entry.second, coordinates, forces) });
            break;
        default:
            break;
        }
    }
}

} // namespace alexandria
