#include <map>

#include "forcecomputerimpl.h"

#include "gromacs/topology/ifunc.h"
#include "gromacs/math/vec.h"

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
