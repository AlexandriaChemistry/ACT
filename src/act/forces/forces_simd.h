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
 * ACTMol::conformationData, the energy maps of the conformations hold the
 * QM energies.
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

/*! \brief Relax the shells of the conformations of all molecules on the GPU
 * As the loops of ForceComputer::relaxMonomer and relaxDimer, after the
 * first calculation, which keeps the energies and forces of the entries
 * without shells, see ConformationData::fixed. Only the entries with
 * shells are computed on the GPU, for all molecules at once, until all
 * shells converged or maxIter is reached. The coordinates, forces,
 * energies and done flags of the conformations are then as those of the
 * loops on the processor, within rounding. Nothing is changed if the GPU
 * cannot do this, such as without GPU, without kept energies, for
 * potentials of entries with shells that the GPU does not have, or with
 * forces on virtual sites to spread.
 * \param[in]    pd               The force field
 * \param[inout] allmols          The molecules
 * \param[in]    what             0 or 1 for a monomer, 2 for the dimers
 * \param[in]    relaxMonomer     For dimers, whether the shells of each monomer are relaxed
 * \param[in]    maxIter          Iterations, including the first calculation
 * \param[in]    msForceToler     Mean square force on the shells for convergence
 * \param[in]    maxShellDistance Largest distance of a shell from its parent, 0 for none
 * \param[in]    minStep          Smallest fraction of the Newton step
 * \param[in]    maxStep          Largest fraction of the Newton step
 * \param[in]    stallIterations  Steps without progress after which shells stop, 0 for never
 * \param[in]    allForces        Whether the forces on all particles are needed
 * \param[out]   nComputed        Number of conformations computed in all iterations
 * \return whether the GPU relaxed the shells
 */
bool relaxShellsGpu(const ForceField            *pd,
                    const std::vector<ACTMol *> &allmols,
                    size_t                       what,
                    const bool                   relaxMonomer[2],
                    int                          maxIter,
                    double                       msForceToler,
                    double                       maxShellDistance,
                    double                       minStep,
                    double                       maxStep,
                    int                          stallIterations,
                    bool                         allForces,
                    size_t                      *nComputed);

/*! \brief Compare the first calculation of the GPU with that of the processor
 * For developers: all topologies of the molecules are copied to the GPU,
 * see gpu::createSystem, all interactions are computed there with the
 * current coordinates and charges, and the energies of each interaction
 * type and conformation are compared with those in the conformation data.
 * The largest differences are written to standard error.
 * \param[in] pd      The force field
 * \param[in] allmols The molecules, after computeEnergies for all their topologies
 * \return false if the GPU cannot compute all interactions
 */
bool checkGpuSystem(const ForceField *pd, const std::vector<ACTMol *> &allmols);

/*! \brief Compute the energies and forces of the conformations of one topology
 * The energies are stored in the conformation data, the forces in its
 * conformations. The shells are not relaxed and the forces on virtual
 * sites are not spread, see ForceComputer::computeTotal.
 * \param[in]    msghandler For warnings and errors
 * \param[in]    pd         An ACT force field
 * \param[in]    top        The topology
 * \param[inout] cdata      The conformations and their data
 * \param[in]    forces     Whether to compute forces
 * \param[in]    skipDone   Whether to leave out the conformations for which the
 *                          shells converged, see computeEnergies
 */
void computeConformations(MsgHandler                *msghandler,
                          const ForceField          *pd,
                          const Topology            *top,
                          ConformationData          *cdata,
                          bool                       forces,
                          bool                       skipDone);

/*! \brief Compute the energies and forces of one potential for all conformations
 * The energies, in ConformationData::energies, and the forces are made
 * anew. This is mainly intended for testing.
 * \param[in]    msghandler    For warnings and errors
 * \param[in]    mypot         The potential
 * \param[in]    topentries    The interactions
 * \param[in]    atoms         The atoms
 * \param[inout] cdata         The conformations, only these are kept, the rest
 *                             of the data for computing them is made anew
 * \param[in]    epsilonr      Relative dielectric constant
 * \param[in]    conformations Whether to loop over the conformations in the
 *                             inner loop (vectorized) or the outer loop
 * \param[in]    forces        Whether to compute forces when looping over the
 *                             conformations in the inner loop, otherwise the
 *                             forces are zero. In the outer loop the forces
 *                             are always computed.
 */
void computePotential(MsgHandler                 *msghandler,
                      Potential                   mypot,
                      const TopologyEntryVector  &topentries,
                      const std::vector<ActAtom> &atoms,
                      ConformationData           *cdata,
                      double                      epsilonr,
                      bool                        conformations,
                      bool                        forces = true);

} // namespace

#endif
