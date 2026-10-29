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
 * \param[in] skipDone   Whether to leave out the conformations for which the
 *                       shells converged, see CoordsEnergy::done, they keep
 *                       their energies and forces. If false, the done flags
 *                       are reset.
 * \param[in] intermolecular Whether to compute only the intermolecular
 *                       interactions of dimers, which gives the interaction
 *                       energies directly if there is no polarization. The
 *                       monomers of the dimers are then not computed, see
 *                       Topology::separateIntermolecular.
 * \param[in] allForces  Whether the forces on all particles are needed. If not,
 *                       only those on the shells are meaningful, which is
 *                       enough for relaxing the shells.
 */
void computeEnergies(MsgHandler                  *msghandler,
                     const ForceField            *pd,
                     const std::vector<ACTMol *> &allmols,
                     size_t                       what,
                     bool                         forces         = true,
                     bool                         skipDone       = false,
                     bool                         intermolecular = false,
                     bool                         allForces      = true);

/*! \brief Compute the energies and forces of the conformations of one topology
 * The energies are stored in the conformation data, the forces in the
 * conformations. The shells are not relaxed and the forces on virtual
 * sites are not spread, see ForceComputer::computeTotal.
 * \param[in]    msghandler For warnings and errors
 * \param[in]    pd         An ACT force field
 * \param[in]    top        The topology
 * \param[inout] coords     The conformations
 * \param[inout] cdata      The data of the conformations
 * \param[in]    forces     Whether to compute forces
 * \param[in]    skipDone   Whether to leave out the conformations for which the
 *                          shells converged, see computeEnergies
 */
void computeConformations(MsgHandler                *msghandler,
                          const ForceField          *pd,
                          const Topology            *top,
                          std::vector<CoordsEnergy> *coords,
                          ConformationData          *cdata,
                          bool                       forces,
                          bool                       skipDone);

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
