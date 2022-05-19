#include "forcecomputerutils.h"

#include <cmath>

#include "gromacs/math/vec.h"
    
real bond_angle(const rvec xi,
                const rvec xj,
                const rvec xk,
                rvec r_ij,
                rvec r_kj,
                real *costh)
/* Return value is the angle between the bonds i-j and j-k */
{
    real th;

    rvec_sub(xi, xj, r_ij);
    rvec_sub(xk, xj, r_kj);

    *costh = cos_angle(r_ij, r_kj);
    th     = std::acos(*costh);

    return th;
}

real dih_angle(const rvec xi,
               const rvec xj,
               const rvec xk,
               const rvec xl,
               rvec r_ij,
               rvec r_kj,
               rvec r_kl,
               rvec m,
               rvec n)
{
    rvec_sub(xi, xj, r_ij);
    rvec_sub(xk, xj, r_kj);
    rvec_sub(xk, xl, r_kl);

    cprod(r_ij, r_kj, m);
    cprod(r_kj, r_kl, n);
    real phi  = gmx_angle(m, n);
    real ipr  = iprod(r_ij, n);
    real sign = (ipr < 0.0) ? -1.0 : 1.0;

    return sign*phi;
}



