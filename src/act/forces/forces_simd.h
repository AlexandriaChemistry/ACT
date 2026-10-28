#ifndef ACT_FORCES_SIMD_H
#define ACT_FORCES_SIMD_H

#include <vector>

#include "act/forcefield/potential.h"
#include "act/molprop/topologyentry.h"
#include "act/topology/actmol.h"

namespace alexandria
{

class ForceField;
class MsgHandler;

/*! \brief Compute energies and forces for monomers and/or dimers.
 * The energies are stored in the conformation data of the molecules, see
 * ACTMol::conformationData, not in the energy maps of the conformations.
 * ConformationEnergies::store copies them there.
 * Unless forces are computed or the force field is polarizable, the energies
 * of an interaction type are taken from the previous call if its parameters
 * (see Topology::parameterVersion), the coordinates (see ConformationData)
 * and the charges did not change.
 * \param[in] msghandler For debugging and messages
 * \param[in] pd         An ACT force field
 * \param[in] allmols    Vector containing all the molecules
 * \param[in] what       If 0, first monomer, if 1, second monomer, if 2 the dimer
 * \param[in] forces     Whether to compute forces. If not, the computed forces
 *                       of the conformations are not meaningful.
 */
void computeEnergies(MsgHandler                  *msghandler,
                     const ForceField            *pd,
                     const std::vector<ACTMol *> &allmols,
                     size_t                       what,
                     bool                         forces = true);

/*! \brief Compute the energies and forces of one potential for all conformations
 * The energy maps and forces of the conformations are cleared first.
 * This is mainly intended for testing.
 * \param[in]    msghandler    For warnings and errors
 * \param[in]    mypot         The potential
 * \param[in]    topentries    The interactions
 * \param[in]    atoms         The atoms
 * \param[inout] coords        The conformations
 * \param[in]    epsilonr      Relative dielectric constant
 * \param[in]    conformations Whether to loop over the conformations in the
 *                             inner loop (vectorized) or the outer loop
 */
void computePotential(MsgHandler                 *msghandler,
                      Potential                   mypot,
                      const TopologyEntryVector  &topentries,
                      const std::vector<ActAtom> &atoms,
                      std::vector<CoordsEnergy>  *coords,
                      double                      epsilonr,
                      bool                        conformations);

} // namespace

#endif
