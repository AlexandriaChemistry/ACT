#include <map>

#include "forcecomputerimpl.h"

#include "act/forces/forcecomputerutils.h"

#include "gromacs/topology/ifunc.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"

namespace alexandria
{

static double computeBonds(const ForceFieldParameterList      &ffpl,
                           const std::vector<TopologyEntry *> &bonds,
                           const std::vector<gmx::RVec>       *coordinates,
                           std::vector<gmx::RVec>             *forces)
{
    double ebond = 0;
    auto   x     = *coordinates;
    auto   f     = *forces;
    const  real half = 0.5;
    for (const auto b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto id         = b->id(); 
        auto bondlength = ffpl.findParameterTypeConst(id, "bondlength").gromacsValue();
        auto kb         = ffpl.findParameterTypeConst(id, "kb").gromacsValue();
        // Get the atom indices
        auto indices    = b->atomIndices();
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto dr         = std::sqrt(dr2) - bondlength;
        
        auto fbond      = kb*dr*gmx::invsqrt(dr2);
        ebond          += half*kb*dr*dr;
        
        for (int m = 0; (m < DIM); m++)
        {
            auto fij          = fbond*dx[m];
            (*forces)[indices[0]][m] += fij;
            (*forces)[indices[1]][m] -= fij;
        }
    }
    return ebond;
}

static double computeMorse(const ForceFieldParameterList      &ffpl,
                           const std::vector<TopologyEntry *> &bonds,
                           const std::vector<gmx::RVec>       *coordinates,
                           std::vector<gmx::RVec>             *forces)
{
    double ebond = 0;
    auto   x     = *coordinates;
    auto   f     = *forces;
    for (const auto b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto id         = b->id(); 
        auto bondlength = ffpl.findParameterTypeConst(id, "bondlength").gromacsValue();
        auto beta       = ffpl.findParameterTypeConst(id, "beta").gromacsValue();
        auto De         = ffpl.findParameterTypeConst(id, "De").gromacsValue();
        auto D0         = ffpl.findParameterTypeConst(id, "D0").gromacsValue();
        // Get the atom indices
        auto indices    = b->atomIndices();
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto dr         = std::sqrt(dr2);
        auto temp       = std::exp(-beta*(dr-bondlength));
        auto omtemp     = 1.0-temp;
        auto cbomtemp   = De*omtemp;
        auto vbond      = D0+cbomtemp*omtemp;
        auto fbond      = -2.0*beta*temp*cbomtemp/dr;
        ebond          += vbond-D0;

        for (int m = 0; (m < DIM); m++)
        {
            auto fij          = fbond*dx[m];
            f[indices[0]][m] += fij;
            f[indices[1]][m] -= fij;
        }
    }
    return ebond;
}

// It is not finished yes, needs to be double checked. 
static double computeAngles(const ForceFieldParameterList      &ffpl,
                            const std::vector<TopologyEntry *> &bonds,
                            const std::vector<gmx::RVec>       *coordinates,
                            std::vector<gmx::RVec>             *forces)
{
    double energy = 0, costh = 0;
    auto   x     = *coordinates;
    auto   f     = *forces;
    const  real half = 0.5;
    for (const auto a : angles)
    {
        // Get the parameters. We have to know their names to do this.
        auto id         = a->id(); 
        auto theta0     = ffpl.findParameterTypeConst(id, "theta0").gromacsValue();
        auto ka         = ffpl.findParameterTypeConst(id, "ka").gromacsValue();
        // Get the atom indices
        auto indices    = a->atomIndices();

		  // atom0 -> i
	     // atom1 -> k
		  // atom2 -> j
        rvec r_ij, r_kj;
		  auto theta = bond_angle(x[indices[0]], x[indices[1]], x[indices[2]], 
										  r_ij, r_kj, &costh);

		  // Convert Dergree to Radian
		  theta0 *= DEG2RAD;
		  theta  *= DEG2RAD;

		  // Compute deviation from the reference angle
		  auto da  = theta - theta0
		  auto da2 = da*da
        
        auto fangle      = ka*da;
        energy          += half*ka*da2;
        
        auto = costh2 = gmx::square(costh);
        if (costh2 < 1)
        {
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            real nrkj_1, nrij_1;
            rvec f_i, f_j, f_k;

            st    = fangle*gmx::invsqrt(1 - costh2);   
            sth   = st*costh;                      

            nrij2 = iprod(r_ij, r_ij);                 
            nrkj2 = iprod(r_kj, r_kj);                 

            nrij_1 = gmx::invsqrt(nrij2);              
            nrkj_1 = gmx::invsqrt(nrkj2);              

            cik = st*nrij_1*nrkj_1;                    
            cii = sth*nrij_1*nrij_1;                   
            ckk = sth*nrkj_1*nrkj_1;                  

            for (auto m = 0; m < DIM; m++)
            {
                f_i[m]    = -(cik*r_kj[m] - cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m] - ckk*r_kj[m]);
                f_j[m]    = -f_i[m] - f_k[m];
                (*forces)[indices[0]][m] += f_i[m];
                (*forces)[indices[1]][m] += f_k[m];
                (*forces)[indices[2]][m] += f_j[m];
            }
        }                                          
    }
    return energy;
}


std::map<int, bondForceComputer> bondForceComputerMap = {
    { F_BONDS, computeBonds },
    { F_MORSE, computeMorse }
};

bondForceComputer getBondForceComputer(int gromacs_index)
{
    auto bfc = bondForceComputerMap.find(gromacs_index);
    if (bondForceComputerMap.end() != bfc)
    {
        return *bfc->second;
    }
    // GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such gromacs_index %d implemented to compute bonded forces for", gromacs_index).c_str()));
    // Keep the compiler happy
    return nullptr;
}

} // namespace
