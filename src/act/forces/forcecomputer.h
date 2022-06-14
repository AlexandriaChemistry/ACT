#ifndef ACT_FORCECOMPUTER_H
#define ACT_FORCECOMPUTER_H
    
#include <vector>

#include "act/poldata/poldata.h"
#include "alexandria/topology.h"
    
namespace alexandria
{

class QtypeProps;

/*! \brief Class to compute all the forces in a molecule or complex
 */
class ForceComputer
{
private:
    //! Force field structure
    const Poldata *pd_;
    //! Convergence criterium for minimizing shells: root mean square force
    double         rmsForce_;
    //! Maximum number of iterations to spend on minimizing shells
    int            maxiter_;
    /*! Do one actual computations.
     * Will do one force/energy computation.
     * \param[in]  pd          The force field structure
     * \param[in]  top         The molecular topology
     * \param[in]  coordinates The atomic coordinates. Coordinates of
     *                         shell particles may be changed.
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     * \param[out] field       Electric field
     */
    void computeOnce(const Topology                    *top,
                     std::vector<gmx::RVec>            *coordinates,
                     std::vector<gmx::RVec>            *forces,
                     std::map<InteractionType, double> *energies,
                     const gmx::RVec                   &field) const;
                     
 public:
    /*! \brief Constructor
     * \param[in] pd       Pointer to force field structure
     * \param[in] rmsForce The root mean square force on shells
     * \param[in] maxiter  The maximum number of iterations for shell minimization
     */
    ForceComputer(const Poldata *pd,
                  double         rmsForce = 0.0001,
                  int            maxiter  = 25) : pd_(pd), rmsForce_(rmsForce), maxiter_(maxiter) {}
    
    /*! Do complete energy/force computation.
     * If shells are present their positions will be minimized.
     * \param[in]  top         The molecular topology
     * \param[in]  charges     The charges for all particles
     * \param[in]  coordinates The atomic coordinates. Coordinates of
     *                         shell particles may be changed.
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     * \return The root mean square force on the shells, or zero if not present.
     */
    double compute(const Topology                    *top,
                   std::vector<gmx::RVec>            *coordinates,
                   std::vector<gmx::RVec>            *forces,
                   std::map<InteractionType, double> *energies) const;
                 
    /*! \brief Return the gromacs type used
     * In practice this converts the InteractionType to the ftype
     * used within the force computer.
     * \param[in] itype The interaction type
     * \return the force type
     */
    int ftype(InteractionType itype) const;

    /*! \brief Compute the polarizability tensor
     * \param[in]  top         The topology
     * \param[in]  coordinates The coordinates
     * \param[out] qtp         The charge type properties
     */
    void calcPolarizability(const Topology         *top,
                            std::vector<gmx::RVec> *coordinates,
                            QtypeProps             *qtp) const;
                 
};

} // namespace alexandria

#endif
