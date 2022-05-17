#include "alexandria/topology.h"
    
namespace alexandria
{

    /*! \brief Class to compute all the forces in a molecule or complex
     */
class ForceComputer
{
 private:
    std::vector<std::pair<int, int> > pairList_;
    
    void makePairList(const Topology &top);
    
    double computeBonds(const std::vector<TopologyEntry *> &bonds,
                        std::vector<gmx::RVec>             *coordinates,
                        std::vector<gmx::RVec>             *forces);
 public:
    /*! \brief Constructor
     */
    ForceComputer(const Topology &top,
                  const Poldata  &pd);
    
    /*! Do the actual computations.
     * Will do all the force computations. If shells are present their 
     * positions will be minimized.
     * \param[in]  coordinates The atomic coordinates. Coordinates of
     *                         shell particles may be changed.
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     */
    void compute(std::vector<gmx::RVec> *coordinates,
                 std::vector<gmx::RVec> *forces,
                 std::vector<double>    *energies);
};

} // namespace alexandria
