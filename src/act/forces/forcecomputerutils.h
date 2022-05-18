#include "gromacs/math/vectypes.h"

/*! \brief Calculate bond-angle. 
 * \param[in] xi The first coordinate i
 * \param[in] xj The second coordinate j
 * \param[in] xk The third coordinate k
 * \param[out] r_ij Distance vector r_j - r_i
 * \param[out] r_kj Distance vector r_j - r_k
 * \param[out] costh The cosine between the bond vectors
 * \return the bond angle
 */
real bond_angle(const rvec  xi, 
                const rvec  xj, 
                const rvec  xk,
                rvec        r_ij, 
                rvec        r_kj, 
                real       *costh);

/*! \brief Calculate dihedral-angle.
 * \param[in] xi The first coordinate i
 * \param[in] xj The second coordinate j
 * \param[in] xk The third coordinate k
 * \param[in] xl The third coordinate l
 * \param[out] r_ij Distance vector r_j - r_i
 * \param[out] r_kj Distance vector r_j - r_k
 * \param[out] r_kl Distance vector r_l - r_k
 * \param[out] m    Normal to plane i-j-k
 * \param[out] n    Normal to plane j-k-l
 * \return the dihedral angle
 */
real dih_angle(const rvec xi,
               const rvec xj,
               const rvec xk,
               const rvec xl,
               rvec       r_ij,
               rvec       r_kj,
               rvec       r_kl,
               rvec       m,
               rvec       n);

