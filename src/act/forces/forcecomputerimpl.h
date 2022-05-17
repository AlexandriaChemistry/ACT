#include <vector>

#include "act/poldata/forcefieldparameterlist.h"
#include "alexandria/topology.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

typedef double (*bondForceComputer)(const ForceFieldParameterList      &ffpl,
                                    const std::vector<TopologyEntry *> &bonds,
                                    const std::vector<gmx::RVec>       *coordinates,
                                    std::vector<gmx::RVec>             *forces);

/*! \brief Return a bonded force computer according to typedef.
 * \param[in] gromacs_index Number corresponding to GROMACS list of energy terms
 * \return pointer to the appropriate function
 * \throws with gmx::InternalError if no such function is implemented.
 */
bondForceComputer getBondForceComputer(int gromacs_index);

} // namespace

