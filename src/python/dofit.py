import math, numpy as np
from jacobi import *

def calc_fit_R(natoms:int, xp, x):
    ndim  = 3
    # Calculate the matrix U
    u     = np.zeros((ndim, ndim))
    for n in range(natoms):
        for c in range(ndim):
           xpc = xp[n][c]
           for r in range(ndim):
               xnr      = x[n][r]
               u[c][r] += xnr*xpc

    # Construct omega
    # Omega is symmetric -> omega==omega'
    omega = np.zeros((2*ndim, 2*ndim))
    for r in range(2*ndim):
        for c in range(0,r+1):
            if r >= ndim and c < ndim:
                omega[r][c] = u[r-ndim][c]
                omega[c][r] = u[r-ndim][c]
            else:
                omega[r][c] = 0
                omega[c][r] = 0

    # Determine h and k
    om    = np.zeros((2*ndim, 2*ndim))
    d     = np.zeros(2*ndim)
    irot = jacobi(omega, 2*ndim, d, om)
    #real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
    # int     natoms = number of rows and columns
    # real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
    # real       **v = v[0..n-1][0..n-1] contains the vectors in columns
    # int      *irot = number of jacobi rotations
    
    if irot == 0:
        print("IROT=0")
        return None

    index = 0
    M_SQRT2 = math.sqrt(2.0)
    vh = np.zeros((ndim,ndim))
    vk = np.zeros((ndim,ndim))
    # Copy only the first ndim-1 eigenvectors
    for j in range(ndim-1):
        max_d = -1000
        for i in range(2*ndim):
            if d[i] > max_d:
                max_d = d[i]
                index = i
        d[index] = -10000
        for i in range(ndim):
            vh[j][i] = M_SQRT2*om[i][index]
            vk[j][i] = M_SQRT2*om[i+ndim][index]

    # Calculate the last eigenvector as the outer-product of the first two.
    # This insures that the conformation is not mirrored and
    # prevents problems with completely flat reference structures.
    vh[2] = np.cross(vh[0], vh[1])
    vk[2] = np.cross(vk[0], vk[1])

    # determine R
    R = np.zeros((3,3))
    for r in range(ndim):
        for c in range(ndim):
            for s in range(ndim):
                R[r][c] += vk[s][r]*vh[s][c]
    return R

