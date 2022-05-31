#include <vector>

#include "act/poldata/poldata.h"
#include "alexandria/topology.h"
    
namespace alexandria
{

/*! \brief Class to compute all the forces in a molecule or complex
 */
class ForceComputer
{
 public:
    /*! \brief Constructor
     */
    ForceComputer() {}
    
    /*! Do the actual computations.
     * Will do all the force computations. If shells are present their 
     * positions will be minimized.
     * \param[in]  pd          The force field structure
     * \param[in]  top         The molecular topology
     * \param[in]  coordinates The atomic coordinates. Coordinates of
     *                         shell particles may be changed.
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     */
    void compute(const Poldata                     &pd,
                 const Topology                    &top,
                 std::vector<gmx::RVec>            *coordinates,
                 std::vector<gmx::RVec>            *forces,
                 std::map<InteractionType, double> *energies);
};

} // namespace alexandria
