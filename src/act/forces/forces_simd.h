#ifndef ACT_FORCES_SIMD_H
#define ACT_FORCES_SIMD_H

#include <vector>

#include "act/topology/actmol.h"

namespace alexandria
{

class ForceField;
class MsgHandler;

/*! \brief Compute energies and forces for monomers and/or dimers.
 * \param[in] msghandler For debugging and messages
 * \param[in] pd         An ACT force field
 * \param[in] allmols    Vector containing all the molecules
 * \param[in] what       If 0, first monomer, if 1, second monomer, if 2 the dimer
 */
void computeEnergies(MsgHandler          *msghandler,
                     const ForceField    *pd,
                     std::vector<ACTMol> *allmols,
                     size_t               what);

} // namespace

#endif
