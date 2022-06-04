#include <map>

#include "forcecomputerimpl.h"

#include "act/forces/forcecomputerutils.h"

#include "gromacs/gmxlib/nonbonded/nb_generic.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"

namespace alexandria
{

static double computeLJ(const ForceFieldParameterList      &ffpl,
                        const std::vector<TopologyEntry *> &pairs,
                        const std::vector<gmx::RVec>       *coordinates,
                        std::vector<gmx::RVec>             *forces)
{
    double ebond = 0;
    return ebond;
}

static double computeWBH(const ForceFieldParameterList      &ffpl,
                         const std::vector<TopologyEntry *> &pairs,
                         const std::vector<gmx::RVec>       *coordinates,
                         std::vector<gmx::RVec>             *forces)
{
    double ebond = 0;
    return ebond;
    auto   x     = *coordinates;
    auto  &f     = *forces;
    for (const auto b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        // For this to work, we need ffpl containing the parameters
        // with the combination rules applied. Something like that does
        // not exist in the input force field files.
        auto id      = b->id();
        auto sigma   = ffpl.findParameterTypeConst(id, "sigma").internalValue();
        auto epsilon = ffpl.findParameterTypeConst(id, "epsilon").internalValue();
        auto gamma   = ffpl.findParameterTypeConst(id, "gamma").internalValue();
        // Get the atom indices
        auto indices    = b->atomIndices();
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto rinv       = gmx::invsqrt(dr2);
        real ewbh, fwbh;
        wang_buckingham(sigma, epsilon, gamma, dr2, rinv, &ewbh, &fwbh);
        ebond      += ewbh;
        real fbond  = fwbh*rinv;
        for (int m = 0; (m < DIM); m++)
        {
            auto fij          = fbond*dx[m];
            f[indices[0]][m] += fij;
            f[indices[1]][m] -= fij;
        }
    }
    return ebond;
}

static double computeCoulomb(const ForceFieldParameterList      &ffpl,
                             const std::vector<TopologyEntry *> &pairs,
                             const std::vector<gmx::RVec>       *coordinates,
                             std::vector<gmx::RVec>             *forces)
{
    double ebond = 0;
    return ebond;
    auto   x     = *coordinates;
    auto  &f     = *forces;
    for (const auto b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        // For this to work, we need ffpl containing the parameters
        // with the combination rules applied. Something like that does
        // not exist in the input force field files.
        auto id    = b->id(); 
        auto qq    = ffpl.findParameterTypeConst(id, "qq").internalValue();
        auto izeta = ffpl.findParameterTypeConst(id, "izeta").internalValue();
        auto jzeta = ffpl.findParameterTypeConst(id, "jzeta").internalValue();
        // Get the atom indices
        auto indices    = b->atomIndices();
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto rinv       = gmx::invsqrt(dr2);
        real velec, felec;
        coulomb_gaussian(qq, izeta, jzeta, std::sqrt(dr2), &velec, &felec);
        ebond      += velec;
        real fbond  = felec*rinv;
        for (int m = 0; (m < DIM); m++)
        {
            auto fij          = fbond*dx[m];
            f[indices[0]][m] += fij;
            f[indices[1]][m] -= fij;
        }
    }
    return ebond;
}

static double computePartridge(const ForceFieldParameterList      &ffpl,
                               const std::vector<TopologyEntry *> &angles,
                               const std::vector<gmx::RVec>       *coordinates,
                               std::vector<gmx::RVec>             *forces)
{
    // Energy function according to 
    // The determination of an accurate isotope dependent potential energy 
    // surface for water from extensive ab initio calculations and 
    // experimental data
    // Harry Partridge and David W. Schwenke
    // J Chem Phys 106 (1997) p. 4618
    double  energy = 0, costh = 0;
    auto    x     = *coordinates;
    auto   &f     = *forces;
    const  real half = 0.5;
    for (const auto a : angles)
    {
        // Get the parameters. We have to know their names to do this.
        auto id         = a->id(); 
        auto rij0       = ffpl.findParameterTypeConst(id, "rij0").internalValue();
        auto rjk0       = ffpl.findParameterTypeConst(id, "rjk0").internalValue();
        auto betaij     = ffpl.findParameterTypeConst(id, "betaij").internalValue();
        auto betajk     = ffpl.findParameterTypeConst(id, "betajk").internalValue();
        auto theta0     = ffpl.findParameterTypeConst(id, "angle").internalValue();
        auto c000       = ffpl.findParameterTypeConst(id, "c000").internalValue();
        auto c111       = ffpl.findParameterTypeConst(id, "c111").internalValue();
        auto c222       = ffpl.findParameterTypeConst(id, "c222").internalValue();
        // Get the atom indices
        auto indices    = a->atomIndices();

        rvec r_ij, r_kj;
        auto theta = bond_angle(x[indices[0]], x[indices[1]], x[indices[2]], 
                                r_ij, r_kj, &costh);
        // Compute the components of the potential
        auto costh0 = std::cos(theta0);
        auto drij   = rij0-norm(r_ij);
        auto drkj   = rjk0-norm(r_kj);
        
        auto costh2 = gmx::square(costh);
        auto st    = theta0*gmx::invsqrt(1 - costh2);   
        auto sth   = st*costh;                      
        
        auto nrij2 = iprod(r_ij, r_ij);                 
        auto nrkj2 = iprod(r_kj, r_kj);                 
        
        auto nrij_1 = gmx::invsqrt(nrij2);              
        auto nrkj_1 = gmx::invsqrt(nrkj2);              
        
        auto cik = st*nrij_1*nrkj_1;                    
        auto cii = sth*nrij_1*nrij_1;                   
        auto ckk = sth*nrkj_1*nrkj_1;                  
        
        for (auto m = 0; m < DIM; m++)
        {
            auto f_im    = -(cik*r_kj[m] - cii*r_ij[m]);
            auto f_km    = -(cik*r_ij[m] - ckk*r_kj[m]);
            auto f_jm    = -f_im - f_km;
            f[indices[0]][m] += f_im;
            f[indices[1]][m] += f_km;
            f[indices[2]][m] += f_jm;
        }                                          
    }
    return energy;
}

static double computeBonds(const ForceFieldParameterList      &ffpl,
                           const std::vector<TopologyEntry *> &bonds,
                           const std::vector<gmx::RVec>       *coordinates,
                           std::vector<gmx::RVec>             *forces)
{
    double ebond = 0;
    auto   x     = *coordinates;
    auto  &f     = *forces;
    const  real half = 0.5;
    for (const auto b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto id         = b->id(); 
        auto bondlength = ffpl.findParameterTypeConst(id, "bondlength").internalValue();
        auto kb         = ffpl.findParameterTypeConst(id, "kb").internalValue();
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
            f[indices[0]][m] += fij;
            f[indices[1]][m] -= fij;
        }
    }
    return ebond;
}

static double computeMorse(const ForceFieldParameterList      &ffpl,
                           const std::vector<TopologyEntry *> &bonds,
                           const std::vector<gmx::RVec>       *coordinates,
                           std::vector<gmx::RVec>             *forces)
{
    double  ebond = 0;
    auto    x     = *coordinates;
    auto   &f     = *forces;
    for (const auto b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto id         = b->id(); 
        auto bondlength = ffpl.findParameterTypeConst(id, "bondlength").internalValue();
        auto beta       = ffpl.findParameterTypeConst(id, "beta").internalValue();
        auto De         = ffpl.findParameterTypeConst(id, "De").internalValue();
        auto D0         = ffpl.findParameterTypeConst(id, "D0").internalValue();
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

static double computeLinearAngles(const ForceFieldParameterList      &ffpl,
                                  const std::vector<TopologyEntry *> &angles,
                                  const std::vector<gmx::RVec>       *coordinates,
                                  std::vector<gmx::RVec>             *forces)
{
    double  ebond = 0;
    auto    x     = *coordinates;
    auto   &f     = *forces;
    for (const auto aaa : angles)
    {
        // Get the parameters. We have to know their names to do this.
        auto id      = aaa->id(); 
        auto a       = ffpl.findParameterTypeConst(id, "a").internalValue();
        auto klin    = ffpl.findParameterTypeConst(id, "klin").internalValue();
        // Get the atom indices
        auto indices = aaa->atomIndices();
        
        rvec r_ij, r_kj, r_ik;
        rvec_sub(x[indices[0]], x[indices[1]], r_ij);
        rvec_sub(x[indices[2]], x[indices[1]], r_kj);
        rvec_sub(r_ij, r_kj, r_ik);

        real dr2 = 0;
        real b   = 1-a;
        for (int m = 0; (m < DIM); m++)
        {
            real dr   = -a * r_ij[m] - b * r_kj[m];
            dr2      += dr*dr;
            real f_i  = a*klin*dr;
            real f_k  = b*klin*dr;
            real f_j  = -(f_i+f_k);
            f[indices[0]][m] += f_i;
            f[indices[1]][m] += f_j;
            f[indices[2]][m] += f_k;
        }
        ebond += 0.5*klin*dr2;
    }
    return ebond;
}

// It is not finished yes, needs to be double checked. 
static double computeAngles(const ForceFieldParameterList      &ffpl,
                            const std::vector<TopologyEntry *> &angles,
                            const std::vector<gmx::RVec>       *coordinates,
                            std::vector<gmx::RVec>             *forces)
{
    double  energy = 0, costh = 0;
    auto    x     = *coordinates;
    auto   &f     = *forces;
    const  real half = 0.5;
    for (const auto a : angles)
    {
        // Get the parameters. We have to know their names to do this.
        auto id         = a->id(); 
        auto theta0     = ffpl.findParameterTypeConst(id, "angle").internalValue();
        auto ka         = ffpl.findParameterTypeConst(id, "kt").internalValue();
        // Get the atom indices
        auto indices    = a->atomIndices();

		  // atom0 -> i
	     // atom1 -> k
		  // atom2 -> j
        rvec r_ij, r_kj;
        auto theta = bond_angle(x[indices[0]], x[indices[1]], x[indices[2]], 
                                r_ij, r_kj, &costh);

        // Compute deviation from the reference angle
        auto da  = theta - theta0;
        auto da2 = da*da;
        
        auto fangle      = ka*da;
        energy          += half*ka*da2;
        
        auto costh2 = gmx::square(costh);
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
                f[indices[0]][m] += f_i[m];
                f[indices[1]][m] += f_k[m];
                f[indices[2]][m] += f_j[m];
            }
        }                                          
    }
    return energy;
}

static double computePolarization(const ForceFieldParameterList      &ffpl,
                                  const std::vector<TopologyEntry *> &angles,
                                  const std::vector<gmx::RVec>       *coordinates,
                                  std::vector<gmx::RVec>             *forces)
{
    return 0;
}

static double computeImpropers(const ForceFieldParameterList      &ffpl,
                               const std::vector<TopologyEntry *> &angles,
                               const std::vector<gmx::RVec>       *coordinates,
                               std::vector<gmx::RVec>             *forces)
{
    return 0;
}

std::map<int, bondForceComputer> bondForceComputerMap = {
    { F_BONDS,         computeBonds        },
    { F_MORSE,         computeMorse        },
    { F_ANGLES,        computeAngles       },
    { F_LINEAR_ANGLES, computeLinearAngles },
    { F_LJ,            computeLJ           },
    { F_BHAM,          computeWBH          },
    { F_COUL_SR,       computeCoulomb      },
    { F_POLARIZATION,  computePolarization },
    { F_IDIHS,         computeImpropers    }
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
