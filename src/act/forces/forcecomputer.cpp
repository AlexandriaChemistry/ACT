/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */
#include "forcecomputer.h"

#include <algorithm>
#include <limits>
#include <set>
#include <utility>
#include <vector>

#include <cstdio>
#include <cstdlib>

#include "act/basics/chargemodel.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forces/conformationdata.h"
#include "act/forces/forces_simd.h"
#include "act/qgen/qtype.h"
#include "act/topology/actmol.h"
#include "act/topology/fragmenthandler.h"
#include "act/topology/topology.h"
#include "act/utility/timings.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

/*! \brief Compute the sum of the squared forces on the shells
 * Used to compute the mean square force on shell particles when
 * minimizing the energy with respect to the shell positions.
 * A shell at the maximum distance from its core with a force pointing
 * away from the core is held there, see ForceComputer::updateShells. For such
 * a shell only the force perpendicular to the core-shell vector counts,
 * since that is the only force that can still be relaxed, and only if it is
 * larger than shellWallAngleTolerance times the force.
 * \param[in] atoms             The particles
 * \param[in] isShell           Whether or not the particles are shells
 * \param[in] x                 The coordinates
 * \param[in] f                 The forces
 * \param[in] offset            Index of the first particle in x and f
 * \param[in] maxShellDistance2 Square of the maximum core-shell distance, 0 if there is none
 * \return The sum of the squared shell forces
 */
/*! \brief The squared length of the cross product of two vectors, in double precision
 * \param[in] a The first vector
 * \param[in] b The second vector
 * \return |a x b|^2
 */
static double perpendicularSquared(const gmx::RVec &a, const gmx::RVec &b)
{
    double c0 = static_cast<double>(a[YY])*b[ZZ] - static_cast<double>(a[ZZ])*b[YY];
    double c1 = static_cast<double>(a[ZZ])*b[XX] - static_cast<double>(a[XX])*b[ZZ];
    double c2 = static_cast<double>(a[XX])*b[YY] - static_cast<double>(a[YY])*b[XX];
    return c0*c0 + c1*c1 + c2*c2;
}

static double shellForce2(const std::vector<ActAtom>               &atoms,
                          const std::vector<bool>                  &isShell,
                          const CoordinateVectors<const Coordinate> &x,
                          const ParticleVectors<const real> &f,
                          size_t                               offset,
                          double                               maxShellDistance2)
{
    // Shells are put on the wall by scaling, allow for rounding
    const double atWall     = 0.9999*maxShellDistance2;
    const double wallAngle2 = shellWallAngleTolerance*shellWallAngleTolerance;
    double       f2sum      = 0;
    // Floating point operations, see Timings::addFlops
    double       flops  = 1;
    for(size_t i = 0; i < isShell.size(); i++)
    {
        if (!isShell[i])
        {
            continue;
        }
        gmx::RVec fi = f[i+offset].rvec();
        double    f2 = iprod(fi, fi);
        // The square and the sum
        flops += 6;
        if (maxShellDistance2 > 0 && atoms[i].parents().size() == 1)
        {
            gmx::RVec dx  = x[i+offset] - x[atoms[i].parents()[0]+offset];
            double    dx2 = iprod(dx, dx);
            double    fdx = iprod(fi, dx);
            flops += 13;
            if (dx2 >= atWall && fdx > 0)
            {
                // Only the force across the core-shell vector can be
                // relaxed, it converged if small compared to the force.
                // From the cross product, since f2 - fdx*fdx/dx2 is the
                // difference of two nearly equal numbers.
                double fperp2 = perpendicularSquared(fi, dx)/dx2;
                f2     = (fperp2 <= wallAngle2*f2) ? 0.0 : fperp2;
                flops += 16;
            }
        }
        f2sum += f2;
    }
    Timings::instance().addFlops(flops);
    return f2sum;
}

double shellMaxDistance(const ForceField *pd)
{
    double maxShellDistance = 0.02;
    auto   itype            = InteractionType::POLARIZATION;
    if (pd->interactionPresent(itype))
    {
        auto &fs = pd->findForcesConst(itype);
        if (fs.optionExists(shellMaxDistanceName))
        {
            maxShellDistance = std::atof(fs.optionValue(shellMaxDistanceName).c_str());
        }
    }
    return maxShellDistance;
}

void ForceComputer::init(double  msForce,
                         int     maxiter,
                         double  epsilonr)
{
    msForceToler_     = msForce;
    maxiter_          = maxiter;
    epsilonr_         = epsilonr;
}

void ForceComputer::constructVsiteCoordinates(const Topology         *top,
                                              std::vector<gmx::RVec> *coordinates) const
{
    // Construct virtual site coordinates
    vsiteHandler_.constructPositions(top, coordinates);
}

void ForceComputer::spreadVsiteForces(const Topology               *top,
                                      const std::vector<gmx::RVec> *coordinates,
                                      std::vector<gmx::RVec>       *forces) const
{
    // Spread virtual site forces
    vsiteHandler_.distributeForces(top, *coordinates, forces);
}

void ForceComputer::constructVsiteCoordinates(const Topology *top,
                                              Coordinates    *coordinates) const
{
    // As one slot of the arrays, relative to the parents
    vsiteHandler_.constructPositions(top, coordinates->view().data(), 1, 1, 1);
}

void ForceComputer::spreadVsiteForces(const Topology         *top,
                                      const Coordinates      *coordinates,
                                      std::vector<gmx::RVec> *forces) const
{
    GMX_RELEASE_ASSERT(coordinates->size() == forces->size(), "Coordinates and forces of different particles");
    if (forces->empty())
    {
        return;
    }
    vsiteHandler_.distributeForces(top, coordinates->view().data(), (*forces)[0].as_vec(), 1, 1, 1);
}

double calc_fscale(double k1, double k2, double F0)
{
    if (F0 == 0)
    {
        return 0;
    }
    if (k2 == 0)
    {
        // F0 = -k1 Delta x => Delta x = -F0/k1
        // However, F0 is taken out because we multiply by the
        // force vector down below.
        return -1/k1;
    }
    else
    {
        // Solution from Mathematica:
        // V[x_] := (1/2) k1 x^2 + k2 x^4
        // F[x_] := -D[V[x], x]
        // MyX = Assuming[{k1 > 0, k2 > 0}, Solve[F0 == F[x], x]]
        // three solutions are obtained, but only the first one is real.
        // Note that this solution assumes the force is independent of shell
        // position, which will never be the case.
        // Also, this assumes a 1D problem where the force on the shell is
        // parallel to the core-shell distance vector.
        // In double precision, since (k1 k2)^3 does not fit in a float
        static const double third = 1.0/3.0;
        static const double three = std::pow(3.0, third);
        static const double sqrt3 = std::sqrt(3.0);
        double k12    = k1*k2;
        double F0k22  = F0*k2*k2;
        double throot = std::pow(-9*F0k22 +
                               sqrt3*std::sqrt(k12*k12*k12 + 27*F0k22*F0k22), third);
        auto k = (-three*three*k12 + three*throot*throot)/(6*k2*throot);
        return k/F0;
    }
}

double shellForceConstant(double k1, double k2, double rh, double r)
{
    double fc = k1;
    if (k2 > 0 && r > rh)
    {
        double drh = r - rh;
        fc        += 4*k2*drh*drh*drh/r;
    }
    return fc;
}

/*! \brief Floating point operations of calc_fscale, see Timings::addFlops
 * \param[in] k2 The quartic force constant
 * \param[in] F  The norm of the force
 * \return the number of operations, with pow counted as log and exp
 */
static double calcFscaleFlops(real k2, real F)
{
    if (F == 0)
    {
        return 0;
    }
    return (k2 == 0) ? 6 : 93;
}

/*! \brief Move a shell on top of its core along the force
 * \param[inout] xs      Component XX of the coordinates of the shell, the others stride apart
 * \param[in]    fs      Likewise of the force on the shell
 * \param[in]    stride  The distance between the components
 * \param[in]    fcShell The harmonic force constant of the shell
 * \param[in]    fcHyper The quartic force constant of the shell
 * \param[in]    step    Fraction of the displacement to the minimum
 * \return the number of floating point operations, see Timings::addFlops
 */
static inline double update_shell_pos(Coordinate *xs,
                                      const real *fs,
                                      size_t      stride,
                                      double      fcShell,
                                      double      fcHyper,
                                      double      step)
{
    gmx::RVec fshell = { fs[0], fs[stride], fs[2*stride] };
    auto      Fnorm  = norm(fshell);
    auto      fscale = step*calc_fscale(fcShell, fcHyper, Fnorm);
    for(int m = 0; m < DIM; m++)
    {
        xs[m*stride] = moveCoordinate(xs[m*stride], -(fshell[m] * fscale));
    }
    // The norm with its square root, the step and the displacement
    return 10 + calcFscaleFlops(fcHyper, Fnorm) + 1 + 6;
}

/*! \brief The displacements that move shells towards the minimum of their potential
 * For consecutive slots, see ConformationArrays, in a loop without
 * branches that the compiler can vectorize. We want to determine the
 * displacement of the shell that cancels the total force, assuming the
 * other forces do not change, that is a Newton step for the spring:
 *
 *                 Ftotal
 *                   ^
 *                    ^
 *          Fspring    ^
 *  Core    <-------- Shell ---> r_core_shell
 *
 * For a spring potential V(r) the stiffness along the core-shell vector
 * is V''(r), perpendicular to it V'(r)/r. With u the unit vector along the
 * core-shell vector, the displacement is
 *
 *   step*(f/kperp + (f.u) u (1/kpar - 1/kperp))
 *
 * which takes three divisions. For a harmonic spring both stiffnesses
 * are the force constant, then the displacement is step*f/k without the
 * core-shell vector. For a quartic spring a shell on top of its core has
 * no core-shell vector, it is moved by update_shell_pos instead.
 * \param[in]  xs      Component XX of the coordinates of the shell in the first slot
 * \param[in]  xc      Likewise of the core of the shell
 * \param[in]  fs      Likewise of the force on the shell
 * \param[in]  stride  The distance between the components
 * \param[in]  ns      The number of slots
 * \param[in]  fcShell The harmonic force constant of the shell
 * \param[in]  fcHyper The quartic force constant of the shell
 * \param[in]  rHyper  The distance beyond which the quartic term acts
 * \param[in]  step    Per slot the fraction of the displacement that would
 *                     bring an isolated shell to its minimum
 * \param[out] disp    The displacements, component m of slot j at m*ns + j
 * \param[out] rcs     The core-shell distance of each slot, only with a
 *                     quartic spring
 * \return the number of floating point operations, see Timings::addFlops
 */
static double shellDisplacements(const Coordinate *xs,
                                 const Coordinate *xc,
                                 const real       *fs,
                                 size_t            stride,
                                 size_t            ns,
                                 double            fcShell,
                                 double            fcHyper,
                                 double            rHyper,
                                 const double     *step,
                                 double           *disp,
                                 double           *rcs)
{
    if (fcHyper == 0)
    {
        double invk = 1/fcShell;
        for(size_t j = 0; j < ns; j++)
        {
            double scale = step[j]*invk;
            for(int m = 0; m < DIM; m++)
            {
                disp[m*ns + j] = scale*fs[j + m*stride];
            }
        }
        return 5 + ns*(1 + DIM);
    }
    double flops = 0;
    for(size_t j = 0; j < ns; j++)
    {
        // Core shell vector, from the integers
        gmx::RVec r;
        gmx::RVec fshell;
        for(int m = 0; m < DIM; m++)
        {
            r[m]      = fromCoordinate(coordinateDifference(xs[j + m*stride], xc[j + m*stride]));
            fshell[m] = fs[j + m*stride];
        }
        double d      = norm(r);
        rcs[j]        = d;
        // On top of the core the displacement is not used
        double invd   = 1/((d == 0) ? 1 : d);
        double drh    = d - rHyper;
        bool   hyper  = d > rHyper;
        double kperp  = hyper ? fcShell + 4*fcHyper*drh*drh*drh*invd : fcShell;
        double kpar   = hyper ? fcShell + 12*fcHyper*drh*drh : fcShell;
        double ikperp = 1/kperp;
        double fpar   = iprod(fshell, r)*invd*invd*(1/kpar - ikperp);
        for(int m = 0; m < DIM; m++)
        {
            disp[m*ns + j] = step[j]*(fshell[m]*ikperp + fpar*r[m]);
        }
        // The vector and its norm with the square root, both force
        // constants, three divisions, the parallel force and the
        // displacement
        flops += 13 + 16 + 3*5 + 9 + DIM*4;
    }
    return flops;
}

void ForceComputer::compute(MsgHandler                        *msg_handler,
                            const ForceField                  *pd,
                            const Topology                    *top,
                            std::vector<gmx::RVec>            *coordinates,
                            std::vector<gmx::RVec>            *forces,
                            std::map<InteractionType, double> *energies,
                            const gmx::RVec                   &field,
                            bool                               resetShells,
                            std::set<int>                      relax,
                            bool                               shellsFromInput) const
{
    Coordinates x(*coordinates);
    compute(msg_handler, pd, top, &x, forces, energies, field, resetShells, relax, shellsFromInput);
    *coordinates = x.toRVec();
}

void ForceComputer::compute(MsgHandler                        *msg_handler,
                            const ForceField                  *pd,
                            const Topology                    *top,
                            Coordinates                       *coordinates,
                            std::vector<gmx::RVec>            *forces,
                            std::map<InteractionType, double> *energies,
                            const gmx::RVec                   &field,
                            bool                               resetShells,
                            std::set<int>                      relax,
                            bool                               shellsFromInput) const
{
    // One conformation, computed with the SIMD routines
    ConformationSet set;
    set.top = top;
    set.setCoordinates(std::vector<Coordinates>{ *coordinates });
    ShellOptions options;
    options.field       = field;
    options.resetShells     = resetShells;
    options.relaxShells     = relax;
    options.shellsFromInput = shellsFromInput;
    std::vector<std::map<InteractionType, double>> emaps;
    computeTotal(msg_handler, pd, &set, options, &emaps);
    *coordinates = Coordinates(set.data.x(0));
    *forces      = set.data.f(0).toRVec();
    *energies    = emaps[0];
}

/*! \brief The shells with one parent particle, with that parent
 * \param[in] atoms The particles
 * \return The shells and their parents
 */
static std::vector<std::pair<size_t, size_t>> shellParents(const std::vector<ActAtom> &atoms)
{
    std::vector<std::pair<size_t, size_t>> shells;
    for(size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].pType() == ActParticle::Shell &&
            atoms[i].parents().size() == 1)
        {
            shells.push_back({ i, atoms[i].parents()[0] });
        }
    }
    return shells;
}

/*! \brief Put shells on top of their parent particle
 * \param[in]    shells The shells and their parents, see shellParents
 * \param[inout] coords The coordinates
 */
static void resetShells(const std::vector<std::pair<size_t, size_t>> &shells,
                        const CoordinateVectors<Coordinate>          &coords)
{
    for(const auto &sp : shells)
    {
        for(int m = 0; m < DIM; m++)
        {
            coords[sp.first][m] = coords[sp.second][m];
        }
    }
}

/*! \brief Construct the virtual sites of one conformation
 * \param[in]    vsiteHandler The virtual site routines
 * \param[in]    top          The topology
 * \param[inout] x            The coordinates
 */
static void constructVsites(const VsiteHandler                  &vsiteHandler,
                            const Topology                      *top,
                            const CoordinateVectors<Coordinate> &x)
{
    vsiteHandler.constructPositions(top, x.data(), x.stride(), 1, 1);
}

/*! \brief Distribute the forces on the virtual sites of one conformation
 * \param[in]    vsiteHandler The virtual site routines
 * \param[in]    top          The topology
 * \param[in]    x            The coordinates
 * \param[inout] f            The forces, those on the virtual sites become zero
 */
static void spreadVsites(const VsiteHandler                        &vsiteHandler,
                         const Topology                            *top,
                         const CoordinateVectors<const Coordinate> &x,
                         const ParticleVectors<real>       &f)
{
    GMX_RELEASE_ASSERT(x.stride() == f.stride(), "Coordinates and forces have different layouts");
    vsiteHandler.distributeForces(top, x.data(), f.data(), x.stride(), 1, 1);
}

//! Smallest scaling of the step of the shells, see updateShells
static const double minShellStep = 0.1;
//! Largest scaling of the step of the shells, see updateShells
static const double maxShellStep = 1.0;

/*! \brief Pointer to a row of arrays in the layout of ConformationArrays
 * \param[in] p      The first slot of the arrays to use
 * \param[in] stride The distance between the rows
 * \param[in] a      The particle
 * \param[in] m      The dimension
 * \return the address of component m of particle a in the first slot
 */
template <typename T>
static inline T *slotRow(T *p, size_t stride, size_t a, int m)
{
    return p + (DIM*a + m)*stride;
}

/*! \brief The squared shell forces of consecutive slots, see shellForce2
 * The loop over the slots is the inner loop, such that it accesses the
 * arrays contiguously. Each slot is computed as by shellForce2.
 * \param[in]  atoms             The particles
 * \param[in]  isShell           Whether or not the particles are shells
 * \param[in]  x                 The coordinates in the first slot, see ConformationArrays
 * \param[in]  f                 The forces in the first slot
 * \param[in]  stride            The distance between the rows of x and f
 * \param[in]  ns                The number of slots
 * \param[in]  offset            Index of the first particle in the conformations
 * \param[in]  maxShellDistance2 Square of the maximum core-shell distance, 0 if there is none
 * \param[out] f2sum             The sum of the squared shell forces in each slot
 */
static void shellForces2(const std::vector<ActAtom> &atoms,
                         const std::vector<bool>    &isShell,
                         const Coordinate           *x,
                         const real                 *f,
                         size_t                      stride,
                         size_t                      ns,
                         size_t                      offset,
                         double                      maxShellDistance2,
                         std::vector<double>        *f2sum)
{
    // Shells are put on the wall by scaling, allow for rounding
    const double atWall     = 0.9999*maxShellDistance2;
    const double wallAngle2 = shellWallAngleTolerance*shellWallAngleTolerance;
    f2sum->assign(ns, 0.0);
    double      *sum    = f2sum->data();
    // Floating point operations, see Timings::addFlops
    double       flops  = ns;
    for(size_t i = 0; i < isShell.size(); i++)
    {
        if (!isShell[i])
        {
            continue;
        }
        bool        wall = maxShellDistance2 > 0 && atoms[i].parents().size() == 1;
        const real *fx   = slotRow(f, stride, i+offset, XX);
        const real *fy   = slotRow(f, stride, i+offset, YY);
        const real *fz   = slotRow(f, stride, i+offset, ZZ);
        // The square and the sum
        flops += 6*ns;
        if (!wall)
        {
            for(size_t j = 0; j < ns; j++)
            {
                real f2 = fx[j]*fx[j] + fy[j]*fy[j] + fz[j]*fz[j];
                sum[j] += f2;
            }
            continue;
        }
        size_t            parent = atoms[i].parents()[0] + offset;
        const Coordinate *xs     = slotRow(x, stride, i+offset, XX);
        const Coordinate *xp     = slotRow(x, stride, parent, XX);
        flops += 13*ns;
        for(size_t j = 0; j < ns; j++)
        {
            gmx::RVec fi = { fx[j], fy[j], fz[j] };
            gmx::RVec dx;
            for(int m = 0; m < DIM; m++)
            {
                dx[m] = fromCoordinate(coordinateDifference(xs[j + m*stride], xp[j + m*stride]));
            }
            double f2  = iprod(fi, fi);
            double dx2 = iprod(dx, dx);
            double fdx = iprod(fi, dx);
            if (dx2 >= atWall && fdx > 0)
            {
                // Only the force across the core-shell vector can be
                // relaxed, it converged if small compared to the force.
                // From the cross product, since f2 - fdx*fdx/dx2 is the
                // difference of two nearly equal numbers.
                double fperp2 = perpendicularSquared(fi, dx)/dx2;
                f2     = (fperp2 <= wallAngle2*f2) ? 0.0 : fperp2;
                flops += 16;
            }
            sum[j] += f2;
        }
    }
    Timings::instance().addFlops(flops);
}

size_t ForceComputer::updateShells(const Topology            *top,
                                   size_t                     offset,
                                   size_t                     nShell,
                                   const std::vector<bool>   &isShell,
                                   size_t                     fragment,
                                   double                     maxShellDistance,
                                   ConformationData          *cdata,
                                   const std::vector<size_t> *active) const
{
    auto  *myCoords   = &cdata->confs;
    auto   itype      = InteractionType::POLARIZATION;
    auto  &pols       = top->entry(itype);
    auto  &atoms      = top->atoms();
    size_t numberDone = 0;
    double maxShellDistance2 = maxShellDistance*maxShellDistance;
    // Floating point operations, see Timings::addFlops, those of shellForces2
    // are counted there
    double flops      = 1;
    auto  &ca         = cdata->arrays;
    auto  &state      = cdata->shells;
    // The slots from s0 up to s1 hold the conformations that are not done.
    // All slots in between are computed, those of conformations that are
    // done or not to update are masked, such that the loops over the slots
    // access the arrays contiguously, as on a vector unit or a GPU. Each
    // slot is computed independently of the others. The work vectors are
    // kept between calls, such that the memory is used again.
    auto  &msForce    = state.msForce;
    auto  &fmax       = state.fmax;
    auto  &after      = state.after;
    auto  &dxForce    = state.dxStep;
    auto  &disp       = state.disp;
    auto  &rcs        = state.rcs;
    auto  &relax      = state.relax;
    auto  &moved      = state.moved;
    size_t nconf      = active ? active->size() : myCoords->size();
    size_t s0         = ca.nconf;
    size_t s1         = 0;
    for(size_t i = 0; i < nconf; i++)
    {
        size_t c = active ? (*active)[i] : i;
        if ((*myCoords)[c].done)
        {
            numberDone += 1;
        }
        else
        {
            s0 = std::min(s0, ca.slot(c));
            s1 = std::max(s1, ca.slot(c) + 1);
        }
    }
    if (s1 <= s0)
    {
        return numberDone;
    }
    size_t ns     = s1 - s0;
    relax.assign(ns, 0);
    for(size_t i = 0; i < nconf; i++)
    {
        size_t c = active ? (*active)[i] : i;
        if (!(*myCoords)[c].done)
        {
            relax[ca.slot(c) - s0] = 1;
        }
    }
    GMX_RELEASE_ASSERT(state.stride == ca.stride && state.dx.size() == ca.x.size(),
                       "The shell relaxation was not started");
    size_t      stride  = ca.stride;
    Coordinate *x       = ca.x.data() + s0;
    const real *f       = ca.f.data() + s0;
    double     *step    = state.step.data() + fragment*stride + s0;
    double     *msPrev  = state.force.data() + fragment*stride + s0;
    double     *dxfPrev = state.dxForce.data() + fragment*stride + s0;
    double     *msRef   = state.reference.data() + fragment*stride + s0;
    double     *nstall  = state.stall.data() + fragment*stride + s0;
    real       *dx      = state.dx.data() + s0;
    // Mean square shell force
    shellForces2(atoms, isShell, x, f, stride, ns, offset, maxShellDistance2, &msForce);
    if (nShell > 0)
    {
        for(size_t j = 0; j < ns; j++)
        {
            msForce[j] /= nShell;
        }
        flops += 5*ns;
    }
    // The forces on the shells do not get below their rounding error,
    // which in single precision may be above the tolerance. When the
    // mean square force is within it and no longer decreases, moving
    // the shells does not reduce it.
    // That is only tested if it no longer decreases, then the largest
    // force is needed.
    auto stalled = [&](size_t j)
    {
        return msForce[j] > msForceToler_ && msForce[j] >= 0.5*msPrev[j];
    };
    bool needMax = false;
    for(size_t j = 0; j < ns; j++)
    {
        needMax = needMax || (relax[j] && stalled(j));
    }
    if (needMax)
    {
        fmax.assign(ns, 0.0);
        for(size_t i = 0; i < atoms.size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                const real *fm = slotRow(f, stride, i+offset, m);
                for(size_t j = 0; j < ns; j++)
                {
                    fmax[j] = std::max(fmax[j], std::abs(static_cast<double>(fm[j])));
                }
            }
        }
    }
    for(size_t j = 0; j < ns; j++)
    {
        if (!relax[j])
        {
            continue;
        }
        bool atFloor = false;
        if (stalled(j))
        {
            atFloor = msForce[j] <= gmx::square(32*GMX_REAL_EPS*fmax[j]);
        }
        if (!(msForce[j] <= msForceToler_ || atFloor))
        {
            // Without halving the mean square force for stallIterations_
            // steps, the shells stop without having converged
            if (msForce[j] <= 0.5*msRef[j])
            {
                msRef[j]  = msForce[j];
                nstall[j] = 0;
            }
            else
            {
                nstall[j] += 1;
            }
            if (stallIterations_ <= 0 || nstall[j] < stallIterations_)
            {
                continue;
            }
        }
        (*myCoords)[ca.conformation(s0 + j)].done = true;
        numberDone += 1;
        relax[j]    = 0;
    }
    // The step is scaled by the fraction of the previous step at which
    // the force along it would have been zero, from the forces before and
    // after it. For an isolated harmonic shell this is one, since the step
    // is a Newton step, coupled shells can overshoot or undershoot.
    after.assign(ns, 0.0);
    for(const auto &p : pols)
    {
        int shellLocal = p->atomIndex(1);
        if (!isShell[shellLocal])
        {
            continue;
        }
        size_t      shellGlobal = shellLocal + offset;
        const real *dxs         = slotRow(dx, stride, shellGlobal, XX);
        const real *fs          = slotRow(f, stride, shellGlobal, XX);
        for(size_t j = 0; j < ns; j++)
        {
            gmx::RVec d  = { dxs[j], dxs[j + stride], dxs[j + 2*stride] };
            gmx::RVec fj = { fs[j], fs[j + stride], fs[j + 2*stride] };
            after[j] += iprod(d, fj);
        }
        flops += 6*ns;
    }
    for(size_t j = 0; j < ns; j++)
    {
        if (!relax[j])
        {
            continue;
        }
        if (dxfPrev[j] > 0)
        {
            double before = dxfPrev[j];
            double scale  = (before > after[j]) ? before/(before - after[j]) : 1.0;
            step[j] = std::clamp(step[j]*scale, minShellStep, maxShellStep);
            flops  += (before > after[j]) ? 7 : 1;
        }
        // If the force increased anyway, the step is reduced, see compute
        if (msForce[j] > msPrev[j])
        {
            step[j] = std::min(step[j], 0.5);
        }
        msPrev[j] = msForce[j];
    }
    dxForce.assign(ns, 0.0);
    disp.resize(DIM*ns);
    rcs.resize(ns);
    moved.assign(ns, 0);
    // Loop over polarizabilities
    for(const auto &p : pols)
    {
        // Displace the shells according to the force
        // Since the potential is harmonic we use Hooke's law
        // F = k dx -> dx = F / k
        // Note that this only works for purely harmonic potentials,
        // see function calc_fscale for quartic potentials.
        // The atom indices stored in the fragment topology are local,
        // whereas the coordinates/forces of the dimer are global, so
        // we add the offset to map local -> global indices.
        int shellLocal   = p->atomIndex(1);
        int coreLocal    = p->atomIndex(0);
        // Only the shells in the mask are relaxed
        if (!isShell[shellLocal])
        {
            continue;
        }
        size_t shellGlobal = shellLocal + offset;
        size_t coreGlobal  = coreLocal + offset;
        // The force constant from the current polarizability and
        // shell charge, since both may be trained, as in compute.
        auto   &params = p->params();
        double  alpha  = params[polALPHA];
        double  qshell = atoms[shellLocal].charge();
        double  fc1    = 0;
        if (alpha > 0 && qshell != 0)
        {
            fc1    = (qshell*qshell*ONE_4PI_EPS0)/alpha;
            flops += 7*ns;
        }
        // The parent of the shell, to check the distance from it
        bool        wall        = maxShellDistance2 > 0 && atoms[shellLocal].parents().size() == 1;
        size_t      coreParent  = wall ? atoms[shellLocal].parents()[0] + offset : 0;
        Coordinate *xs          = slotRow(x, stride, shellGlobal, XX);
        Coordinate *xc          = slotRow(x, stride, coreGlobal, XX);
        Coordinate *xp          = slotRow(x, stride, coreParent, XX);
        const real *fs          = slotRow(f, stride, shellGlobal, XX);
        real       *dxs         = slotRow(dx, stride, shellGlobal, XX);
        double      fcHyper     = params[polFCHYPER];
        double      rHyper      = params[polRHYPER];
        bool        quartic     = fcHyper > 0;
        flops += shellDisplacements(xs, xc, fs, stride, ns, fc1, fcHyper, rHyper, step,
                                    disp.data(), rcs.data());
        for(size_t j = 0; j < ns; j++)
        {
            if (!relax[j])
            {
                continue;
            }
            // The displacement is computed from the integers
            Coordinate xold[DIM] = { xs[j], xs[j + stride], xs[j + 2*stride] };
            if (quartic && rcs[j] == 0)
            {
                // On top of the core the quartic term only contributes if
                // it starts at zero distance.
                flops += update_shell_pos(xs + j, fs + j, stride, fc1,
                                          rHyper > 0 ? 0 : fcHyper, step[j]);
            }
            else
            {
                for(int m = 0; m < DIM; m++)
                {
                    xs[j + m*stride] = moveCoordinate(xs[j + m*stride], disp[m*ns + j]);
                }
            }
            bool clipped = false;
            // Check distance from parent, if there is only one
            if (wall)
            {
                gmx::RVec d;
                for(int m = 0; m < DIM; m++)
                {
                    d[m] = fromCoordinate(coordinateDifference(xs[j + m*stride], xp[j + m*stride]));
                }
                real dx2 = iprod(d, d);
                flops   += 8;
                if (dx2 > maxShellDistance2)
                {
                    // Move back the shell/drude to the wall distance
                    real scale = std::sqrt(maxShellDistance2/dx2);
                    for (int m = 0; m < DIM; m++)
                    {
                        xs[j + m*stride] = toCoordinate(fromCoordinate(xp[j + m*stride]) + scale*d[m]);
                    }
                    clipped = true;
                    flops  += 16;
                }
            }
            for(int m = 0; m < DIM; m++)
            {
                moved[j] = moved[j] || xs[j + m*stride] != xold[m];
            }
            // Shells held at the maximum distance do not tell how far to step
            gmx::RVec dxj = { 0, 0, 0 };
            if (!clipped)
            {
                for(int m = 0; m < DIM; m++)
                {
                    dxj[m] = fromCoordinate(coordinateDifference(xs[j + m*stride], xold[m]));
                }
                gmx::RVec fj = { fs[j], fs[j + stride], fs[j + 2*stride] };
                dxForce[j] += iprod(dxj, fj);
                flops      += 9;
            }
            for(int m = 0; m < DIM; m++)
            {
                dxs[j + m*stride] = dxj[m];
            }
        }
    }
    for(size_t j = 0; j < ns; j++)
    {
        if (!relax[j])
        {
            continue;
        }
        dxfPrev[j] = dxForce[j];
        // The displacements of all shells were below the resolution of the
        // coordinates, see coordinateUnit, so the shells did not move. The
        // next steps would be the same, so the shells are as converged as
        // the coordinates allow. With very stiff springs this happens above
        // the tolerance on the force, for instance at 1e-4 kJ/mol nm with a
        // force constant of 1e13 kJ/mol nm^2.
        if (!moved[j])
        {
            (*myCoords)[ca.conformation(s0 + j)].done = true;
            numberDone += 1;
            relax[j]    = 0;
        }
    }
    Timings::instance().addFlops(flops);
    return numberDone;
}

//! \return Whether an interaction type changes upon shell relaxation
static inline bool isElecTerm(InteractionType itype)
{
    return (InteractionType::ELECTROSTATICS == itype ||
            InteractionType::POLARIZATION   == itype ||
            InteractionType::CHARGETRANSFER == itype);
}

/*! \brief Compute the induction, all electrostatics and exchange plus induction
 * The induction is the change of the electrostatic terms since their backup,
 * and these terms are restored to their backup.
 * \param[inout] cener The energies of all conformations
 */
static void calcInduction(ConformationEnergies *cener)
{
    // The terms present, in the same order as in an energy map.
    // Note that the INDUCTIONCORRECTION is treated in the calling routine.
    std::vector<std::pair<InteractionType, ConformationEnergies::Term>> terms;
    for(auto t : cener->types())
    {
        if (isElecTerm(t) || InteractionType::EXCHANGE == t)
        {
            terms.push_back({ t, cener->find(t) });
        }
    }
    auto induction = cener->find(InteractionType::INDUCTION);
    auto eAllElec  = cener->add(InteractionType::ALLELEC);
    auto eExchInd  = cener->add(InteractionType::EXCHIND);
    if (!cener->haveBackups())
    {
        // Nothing was relaxed, then the sums are those below without
        // induction from shells, added one term at a time without
        // branches, in the same order.
        size_t nconf = cener->nconf();
        for(auto eSum : { eAllElec, eExchInd })
        {
            for(size_t c = 0; c < nconf; c++)
            {
                eSum.before[c] = eSum.have[c] ? eSum.before[c] : 0.0;
                eSum.have[c]   = 1;
                eSum.act[c]    = 0.0;
            }
        }
        auto addTerm = [nconf](const ConformationEnergies::Term &et, double *sum)
        {
            for(size_t c = 0; c < nconf && et.have; c++)
            {
                sum[c] += et.have[c] ? et.act[c] : 0.0;
            }
        };
        for(const auto &tt : terms)
        {
            addTerm(tt.second, isElecTerm(tt.first) ? eAllElec.act : eExchInd.act);
        }
        if (induction.have)
        {
            addTerm(induction, eExchInd.act);
        }
        return;
    }
    for(size_t c = 0; c < cener->nconf(); c++)
    {
        bool   haveInduction = false;
        double eInduction    = 0;
        double allElec       = 0;
        double exchInd       = 0;
        for(auto &tt : terms)
        {
            auto &et = tt.second;
            if (!et.has(c))
            {
                continue;
            }
            if (isElecTerm(tt.first))
            {
                // Without backups nothing was relaxed
                if (cener->haveBackups() && et.before[c] != et.act[c])
                {
                    eInduction += et.act[c] - et.before[c];
                    et.restoreBackup(c);
                    haveInduction = true;
                }
                // Copy new value
                allElec += et.act[c];
            }
            else
            {
                exchInd += et.act[c];
            }
        }
        if (haveInduction)
        {
            // This is the total induction including the induction correction,
            // which may be split off later. The induction and induction
            // correction are not up to date yet.
            if (!induction.have)
            {
                induction = cener->add(InteractionType::INDUCTION);
            }
            induction.set(c, eInduction);
            allElec += eInduction;
            exchInd += eInduction;
        }
        else if (induction.has(c))
        {
            // Induction from a potential function rather than from shells
            exchInd += induction.act[c];
        }
        eAllElec.set(c, allElec);
        eExchInd.set(c, exchInd);
    }
}

/*! \brief Take the induction correction of a potential into account
 * A potential for the INDUCTIONCORRECTION, e.g. MACDANIEL_SCHMIDT, is part
 * of all electrostatics and of exchange plus induction, as in
 * ForceComputer::compute. If it is not treated separately, it is moved to
 * the induction.
 * \param[inout] cener      The energies of all conformations
 * \param[in]    correction The induction correction of the potential per conformation
 * \param[in]    toInduction Whether to move the correction to the induction
 */
static void addCorrection(ConformationEnergies      *cener,
                          const std::vector<double> &correction,
                          bool                       toInduction)
{
    auto eAllElec = cener->add(InteractionType::ALLELEC);
    auto eExchInd = cener->add(InteractionType::EXCHIND);
    for(size_t c = 0; c < cener->nconf() && c < correction.size(); c++)
    {
        if (eAllElec.has(c))
        {
            eAllElec.act[c] += correction[c];
        }
        if (eExchInd.has(c))
        {
            eExchInd.act[c] += correction[c];
        }
    }
    if (toInduction)
    {
        auto induction = cener->add(InteractionType::INDUCTION);
        auto icorr     = cener->find(InteractionType::INDUCTIONCORRECTION);
        for(size_t c = 0; c < cener->nconf() && c < correction.size(); c++)
        {
            if (icorr.has(c))
            {
                induction.set(c, (induction.has(c) ? induction.act[c] : 0.0) + correction[c]);
                icorr.have[c] = 0;
            }
        }
    }
}

/*! \brief The induction correction of a potential in all conformations
 * \param[in] cener The energies of all conformations
 * \return The correction of each conformation, empty if there is none
 */
static std::vector<double> potentialCorrection(ConformationEnergies *cener)
{
    std::vector<double> correction;
    auto icorr = cener->find(InteractionType::INDUCTIONCORRECTION);
    if (icorr.have)
    {
        correction.resize(cener->nconf(), 0.0);
        for(size_t c = 0; c < cener->nconf(); c++)
        {
            if (icorr.has(c))
            {
                correction[c] = icorr.act[c];
            }
        }
    }
    return correction;
}

/*! \brief Sum the ACT energy terms that change upon shell relaxation
 * \param[inout] cener The energies of all conformations
 * \return The sum of electrostatics, polarization and charge transfer for each conformation
 */
static std::vector<double> sumElecTerms(ConformationEnergies *cener)
{
    std::vector<double> esum(cener->nconf(), 0.0);
    for(const auto itype : { InteractionType::ELECTROSTATICS,
                             InteractionType::POLARIZATION,
                             InteractionType::CHARGETRANSFER })
    {
        auto ee = cener->find(itype);
        for(size_t c = 0; c < cener->nconf(); c++)
        {
            if (ee.has(c))
            {
                esum[c] += ee.act[c];
            }
        }
    }
    return esum;
}

/*! \brief Back up all ACT energies
 * \param[inout] cener The energies of all conformations
 */
static void backupEnergies(ConformationEnergies *cener)
{
    // The backups are only used for conformations with an energy, so
    // copying all of them is the same as copying only those.
    for(auto t : cener->types())
    {
        auto et = cener->find(t);
        std::copy(et.act, et.act + cener->nconf(), et.before);
    }
}

/*! \brief Subtract the energies of the monomers from those of the dimer
 * The backups are subtracted as well.
 * \param[inout] dimer    The energies of the dimer conformations
 * \param[inout] monomers The energies of the monomer conformations
 */
static void subtractMonomers(ConformationEnergies *dimer, ConformationEnergies *monomers[2])
{
    // Without backups, they are the same as the energies
    bool backups = dimer->haveBackups();
    for(auto t : dimer->types())
    {
        auto ed = dimer->find(t);
        ConformationEnergies::Term em[2] = { monomers[0]->find(t), monomers[1]->find(t) };
        for(size_t c = 0; c < dimer->nconf(); c++)
        {
            if (!ed.has(c))
            {
                continue;
            }
            double edim = ed.act[c];
            for(size_t j = 0; j < 2; j++)
            {
                if (c < monomers[j]->nconf() && em[j].has(c))
                {
                    edim -= em[j].act[c];
                }
            }
            ed.act[c] = edim;
            if (backups)
            {
                double edb4 = ed.before[c];
                for(size_t j = 0; j < 2; j++)
                {
                    if (c < monomers[j]->nconf() && em[j].has(c))
                    {
                        edb4 -= em[j].before[c];
                    }
                }
                ed.before[c] = edb4;
            }
        }
    }
}

/*! \brief Check whether the virtual sites have to be constructed
 * The atoms do not move, therefore the virtual sites only have to be
 * constructed the first time and when their parameters changed, unless
 * shells are used to construct them.
 * \param[in]    top   The topology
 * \param[inout] cdata The conformation data, updated if the virtual sites
 *                     have to be constructed
 * \return Whether the virtual sites have to be constructed
 */
static bool vsitesOutdated(const Topology *top, ConformationData *cdata)
{
    bool outdated = (!cdata->haveVsites || cdata->vsiteVersion != top->vsiteVersion());
    if (!outdated)
    {
        auto &atoms = top->atoms();
        for(const auto &entry : top->entries())
        {
            if (!isVsite(entry.first))
            {
                continue;
            }
            for(const auto &vs : entry.second)
            {
                for(auto a : vs->atomIndices())
                {
                    outdated = outdated || atoms[a].pType() == ActParticle::Shell;
                }
            }
        }
    }
    if (outdated)
    {
        cdata->haveVsites   = true;
        cdata->vsiteVersion = top->vsiteVersion();
        cdata->coordinates.update();
    }
    return outdated;
}

/*! \brief Whether the virtual sites can be constructed after resetting the shells
 * This is the case if no virtual site is constructed from a shell and no
 * shell is put on top of a virtual site, see resetShells. The energy
 * calculation then constructs them, see ConformationData::constructVsites.
 * \param[in] top The topology
 * \return Whether the virtual sites do not depend on the shells
 */
static bool vsitesIndependentOfShells(const Topology *top)
{
    auto &atoms = top->atoms();
    for(const auto &entry : top->entries())
    {
        if (!isVsite(entry.first))
        {
            continue;
        }
        for(const auto &vs : entry.second)
        {
            for(auto a : vs->atomIndices())
            {
                if (atoms[a].pType() == ActParticle::Shell)
                {
                    return false;
                }
            }
        }
    }
    for(const auto &atom : atoms)
    {
        if (atom.pType() == ActParticle::Shell)
        {
            for(auto p : atom.parents())
            {
                if (atoms[p].pType() != ActParticle::Atom)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

//! A component in units of coordinateUnit, to copy it exactly
struct CoordinateUnits
{
    //! The component
    Coordinate value;
    //! \return the component
    Coordinate units() const { return value; }
};

/*! \brief Call a function for each conformation of each topology of the molecules
 * \param[in] allmols The molecules
 * \param[in] what    0 or 1 for a monomer, 2 for the dimers, -2 for both monomers
 * \param[in] func    Called with the conformation data, its conformation and the topology
 */
template <typename F>
static void forEachConformation(const std::vector<ACTMol *> &allmols, int what, F func)
{
    for(auto molptr : allmols)
    {
        auto &tops = molptr->fragmentHandler()->topologies();
        for(int w = 0; w < 3; w++)
        {
            if (!(w == what || (what == -2 && w < 2)))
            {
                continue;
            }
            if ((w < 2 && static_cast<size_t>(w) >= tops.size()) || (w == 2 && !molptr->isDimer()))
            {
                continue;
            }
            const Topology *top   = (w < 2) ? &tops[w] : molptr->topology();
            auto           *cdata = molptr->conformationData(w);
            for(size_t c = 0; c < cdata->nconf(); c++)
            {
                func(cdata, c, top);
            }
        }
    }
}

/*! \brief Keep the shells of the conformations that converged, see ForceComputer::setKeepShells
 * \param[inout] allmols The molecules
 */
static void keepRelaxedShells(const std::vector<ACTMol *> &allmols)
{
    for(int what : { -2, 2 })
    {
        forEachConformation(allmols, what, [](ConformationData *cdata, size_t c, const Topology *top)
        {
            auto  x     = cdata->x(c);
            auto &kept  = cdata->confs[c].xkept;
            auto &atoms = top->atoms();
            kept.clear();
            // Shells that did not converge start from their parents next time
            if (!cdata->confs[c].done)
            {
                return;
            }
            for(size_t i = 0; i < atoms.size(); i++)
            {
                if (atoms[i].pType() == ActParticle::Shell)
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        kept.push_back(x[i][m].units());
                    }
                }
            }
        });
    }
}

/*! \brief Put back the shells that keepRelaxedShells kept
 * \param[inout] allmols The molecules
 * \param[in]    what    As for forEachConformation
 * \return whether shells were put back in any conformation
 */
static bool restoreRelaxedShells(const std::vector<ACTMol *> &allmols, int what)
{
    bool restored = false;
    forEachConformation(allmols, what, [&restored](ConformationData *cdata, size_t c, const Topology *top)
    {
        auto        x     = cdata->x(c);
        const auto &kept  = cdata->confs[c].xkept;
        auto       &atoms = top->atoms();
        size_t      k     = 0;
        for(size_t i = 0; i < atoms.size(); i++)
        {
            if (atoms[i].pType() == ActParticle::Shell)
            {
                if (k + DIM > kept.size())
                {
                    return;
                }
                for(int m = 0; m < DIM; m++)
                {
                    x[i][m] = CoordinateUnits{ kept[k++] };
                }
            }
        }
        restored = restored || k > 0;
    });
    return restored;
}

/*! \brief Copy the relaxed monomer shell positions into the dimer
 * \param[inout] allmols The molecules
 */
static void copyRelaxedMonomerShells(const std::vector<ACTMol *> &allmols)
{
    for(auto molptr : allmols)
    {
        auto &mol = *molptr;
        if (!mol.isDimer())
        {
            continue;
        }
        auto  fh     = mol.fragmentHandler();
        auto &tops   = fh->topologies();
        auto &astart = fh->atomStart();
        auto *dimers = mol.conformationData(2);
        for(size_t im = 0; im < tops.size(); im++)
        {
            auto &atoms = tops[im].atoms();
            auto &mons  = *mol.monomerCoords(im);
            for(size_t i = 0; i < dimers->nconf(); i++)
            {
                auto        xd = dimers->x(i);
                const auto &xr = mons[i].xrelax;
                for(size_t j = 0; j < atoms.size(); j++)
                {
                    if (atoms[j].pType() == ActParticle::Shell)
                    {
                        for(int m = 0; m < DIM; m++)
                        {
                            xd[astart[im]+j][m] = xr[j][m];
                        }
                    }
                }
            }
        }
    }
}

/*! \brief Start the shell relaxation with full steps
 * \param[inout] cdata The conformations
 */
static void resetShellSteps(ConformationData *cdata)
{
    cdata->shells.reset(cdata->arrays.stride, cdata->arrays.natom);
}

/*! \brief Whether the shells of a fragment of a conformation stopped without progress
 * \param[in] cdata    The conformations
 * \param[in] c        The conformation
 * \param[in] fragment The fragment
 * \param[in] stallIterations The steps without progress after which shells stop, 0 for never
 * \return Whether the relaxation stopped, see ForceComputer::setStallIterations
 */
static bool shellsStalled(const ConformationData *cdata, size_t c, size_t fragment, int stallIterations)
{
    const auto &st = cdata->shells;
    size_t      i  = fragment*st.stride + cdata->arrays.slot(c);
    return stallIterations > 0 && i < st.stall.size() && st.stall[i] >= stallIterations;
}

void ForceComputer::relaxMonomer(MsgHandler                  *msghandler,
                                 const ForceField            *pd,
                                 size_t                       im,
                                 const std::vector<ACTMol *> &allmols,
                                 bool                         allForces) const
{
    static const size_t timer = Timings::instance().index("Shell relaxation");
    ScopedTimer         scoped(timer);
    double              maxShellDistance = shellMaxDistance(pd);
    // The conformations that did not converge of each molecule, such that
    // the others are not visited in each iteration
    std::vector<std::vector<size_t>> active(allmols.size());
    for(size_t k = 0; k < allmols.size(); k++)
    {
        auto &mons = *allmols[k]->monomerCoords(im);
        resetShellSteps(allmols[k]->conformationData(im));
        for(size_t c = 0; c < mons.size(); c++)
        {
            if (!mons[c].done)
            {
                active[k].push_back(c);
            }
        }
    }
    // We have already done one energy calculation
    int iter     = 1;
    // Prevent shells from flying off.
    bool allDone = false;
    // On the GPU, if possible, for all molecules at once
    const bool relaxAll[2] = { true, true };
    size_t     nComputed   = 0;
    if (relaxShellsGpu(pd, allmols, im, relaxAll, maxiter_, msForceToler_, maxShellDistance,
                       minShellStep, maxShellStep, stallIterations_, allForces, &nComputed))
    {
        nShellIter_ += nComputed;
        allDone      = true;
    }
    while (iter < maxiter_ && !allDone)
    {
        // Set loop variable to true
        allDone = true;
        for(size_t k = 0; k < allmols.size(); k++)
        {
            auto &mol = *allmols[k];
            // Shortcuts...
            auto  fh    = mol.fragmentHandler();
            auto &tops  = fh->topologies();
            auto &pols  = fh->polInfo();
            // A compound with one fragment has no second monomer
            if (im >= tops.size())
            {
                continue;
            }
            auto &mons  = *mol.monomerCoords(im);
            auto &act   = active[k];
            // Relax monomers first.
            updateShells(&tops[im], 0, pols[im].nShell(), pols[im].isShell(),
                         0, maxShellDistance, mol.conformationData(im), &act);
            act.erase(std::remove_if(act.begin(), act.end(),
                                     [&mons](size_t c) { return mons[c].done; }),
                      act.end());
            // Not all done, recompute
            allDone = allDone && act.empty();
        }
        if (!allDone)
        {
            // Only for the conformations that did not converge
            for(const auto &act : active)
            {
                nShellIter_ += act.size();
            }
            computeEnergies(msghandler, pd, allmols, im, true, true, false, allForces);
        }
        iter++;
    }
    // Check the conformations computed last
    double maxShellDistance2 = maxShellDistance*maxShellDistance;
    for(auto molptr : allmols)
    {
        auto  fh   = molptr->fragmentHandler();
        auto &tops = fh->topologies();
        auto &pols = fh->polInfo();
        if (im >= tops.size() || pols[im].nShell() == 0)
        {
            continue;
        }
        const auto *cdata = molptr->conformationData(im);
        for(size_t c = 0; c < cdata->nconf(); c++)
        {
            if (cdata->confs[c].done ? shellsStalled(cdata, c, 0, stallIterations_) :
                (shellForce2(tops[im].atoms(), pols[im].isShell(), cdata->x(c), cdata->f(c),
                             0, maxShellDistance2)/pols[im].nShell() > msForceToler_))
            {
                nShellConvergeFailed_ += 1;
            }
        }
    }
}

void ForceComputer::relaxDimer(MsgHandler                  *msghandler,
                               const ForceField            *pd,
                               bool                         fixMonomer[2],
                               const std::vector<ACTMol *> &allmols,
                               bool                         allForces) const
{
    static const size_t timer = Timings::instance().index("Shell relaxation");
    ScopedTimer         scoped(timer);
    double              maxShellDistance = shellMaxDistance(pd);
    // The conformations that did not converge of each molecule, such that
    // the others are not visited in each iteration
    std::vector<std::vector<size_t>> active(allmols.size());
    for(size_t k = 0; k < allmols.size(); k++)
    {
        auto &dims = *allmols[k]->dimerCoords();
        resetShellSteps(allmols[k]->conformationData(2));
        for(size_t c = 0; c < dims.size(); c++)
        {
            if (!dims[c].done)
            {
                active[k].push_back(c);
            }
        }
    }
    // Whether each active conformation is done, kept to use the memory again
    std::vector<char> done;
    // We have already done one energy calculation
    int iter     = 1;
    // Prevent shells from flying off.
    bool allDone = false;
    // On the GPU, if possible, for all molecules at once
    const bool relaxMonomers[2] = { !fixMonomer[0], !fixMonomer[1] };
    size_t     nComputed        = 0;
    if (relaxShellsGpu(pd, allmols, 2, relaxMonomers, maxiter_, msForceToler_, maxShellDistance,
                       minShellStep, maxShellStep, stallIterations_, allForces, &nComputed))
    {
        nShellIter_ += nComputed;
        allDone      = true;
    }
    while (iter < maxiter_ && !allDone)
    {
        // Set loop variable to true
        allDone = true;
        for(size_t k = 0; k < allmols.size(); k++)
        {
            auto &mol = *allmols[k];
            // Shortcuts...
            auto  fh    = mol.fragmentHandler();
            auto &tops  = fh->topologies();
            auto &pols  = fh->polInfo();
            // A compound with one fragment has no dimer conformations
            if (tops.size() < 2)
            {
                continue;
            }
            // The done flag is per dimer conformation, but convergence
            // has to be checked for each monomer separately. A conformation
            // is done when the shells of all monomers that are relaxed converged.
            auto &dims = *mol.dimerCoords();
            auto &act  = active[k];
            done.assign(act.size(), 1);
            for(size_t im = 0; im < 2; im++)
            {
                if (!fixMonomer[im])
                {
                    // The active conformations are not done
                    for(auto c : act)
                    {
                        dims[c].done = false;
                    }
                    updateShells(&tops[im], fh->atomStart()[im],
                                 pols[im].nShell(), pols[im].isShell(),
                                 im, maxShellDistance, mol.conformationData(2), &act);
                    for(size_t i = 0; i < act.size(); i++)
                    {
                        done[i] = done[i] && dims[act[i]].done;
                    }
                }
            }
            size_t nleft = 0;
            for(size_t i = 0; i < act.size(); i++)
            {
                dims[act[i]].done = done[i];
                if (!done[i])
                {
                    act[nleft++] = act[i];
                }
            }
            act.resize(nleft);
            allDone = allDone && act.empty();
        }
        if (!allDone)
        {
            // Not all done, recompute the conformations that did not converge
            for(const auto &act : active)
            {
                nShellIter_ += act.size();
            }
            computeEnergies(msghandler, pd, allmols, 2, true, true, false, allForces);
        }
        iter++;
    }
    // Check the conformations computed last, the shells of the monomers that
    // are relaxed have to be converged
    double maxShellDistance2 = maxShellDistance*maxShellDistance;
    for(auto molptr : allmols)
    {
        auto  fh   = molptr->fragmentHandler();
        auto &tops = fh->topologies();
        auto &pols = fh->polInfo();
        const auto *cdata = molptr->conformationData(2);
        for(size_t c = 0; c < cdata->nconf(); c++)
        {
            bool done   = cdata->confs[c].done;
            bool failed = false;
            for(size_t im = 0; im < tops.size() && im < 2; im++)
            {
                if (fixMonomer[im] || pols[im].nShell() == 0)
                {
                    continue;
                }
                if (done ? shellsStalled(cdata, c, im, stallIterations_) :
                    (shellForce2(tops[im].atoms(), pols[im].isShell(), cdata->x(c), cdata->f(c),
                                 fh->atomStart()[im], maxShellDistance2)/pols[im].nShell() > msForceToler_))
                {
                    failed = true;
                }
            }
            nShellConvergeFailed_ += failed ? 1 : 0;
        }
    }
}

/*! \brief Extract the induction and electrostatics from the energies of relaxed shells
 * The last step of ForceComputer::computeRelaxedEnergies, also used when
 * the GPU computes the energies.
 * \param[in] allmols        The molecules
 * \param[in] intermolecular Whether only the intermolecular interactions were computed
 * \param[in] doInduc2       Whether the second order induction was computed separately
 * \param[in] separateInductionCorrection Whether the induction correction is separate
 * \param[in] eInduc2        The second order induction of each conformation of each molecule
 */
static void postProcessRelaxedEnergies(const std::vector<ACTMol *>            &allmols,
                                       bool                                    intermolecular,
                                       bool                                    doInduc2,
                                       bool                                    separateInductionCorrection,
                                       const std::vector<std::vector<double>> &eInduc2)
{
    // Now post-process the energies to extract induction and allelec
    for(const auto &molptr : allmols)
    {
        auto &mol = *molptr;
        auto *dimer = &mol.conformationData(2)->energies;
        ConformationEnergies *monomers[2] = { &mol.conformationData(0)->energies,
                                              &mol.conformationData(1)->energies };
        bool inter = intermolecular && mol.isDimer() && mol.topology()->hasIntermolecular();
        if (mol.isDimer())
        {
            // First, subtract monomer energies from dimer energies,
            // unless only the intermolecular interactions were computed
            if (!inter)
            {
                subtractMonomers(dimer, monomers);
            }
            calcInduction(dimer);
            // The induction correction from a potential, before the
            // induction correction from the shells is added
            auto correction = potentialCorrection(dimer);
            if (doInduc2)
            {
                // Split the total induction into second order induction and
                // the remaining higher order induction correction.
                size_t k     = &molptr - allmols.data();
                auto   eit   = dimer->find(InteractionType::INDUCTION);
                std::vector<double> eInd(dimer->nconf(), 0.0);
                for(size_t i = 0; i < dimer->nconf(); i++)
                {
                    if (eit.has(i))
                    {
                        eInd[i] = eit.act[i];
                    }
                }
                auto induction   = dimer->add(InteractionType::INDUCTION);
                auto icorrection = dimer->add(InteractionType::INDUCTIONCORRECTION);
                for(size_t i = 0; i < dimer->nconf() && i < eInduc2[k].size(); i++)
                {
                    double ecorr = i < correction.size() ? correction[i] : 0.0;
                    induction.set(i, eInduc2[k][i]);
                    icorrection.set(i, ecorr + eInd[i] - eInduc2[k][i]);
                }
            }
            if (!correction.empty())
            {
                addCorrection(dimer, correction, !separateInductionCorrection);
            }
        }
        auto &tops = mol.fragmentHandler()->topologies();
        for(size_t i = 0; i < tops.size() && !inter; i++)
        {
            calcInduction(monomers[i]);
            auto correction = potentialCorrection(monomers[i]);
            if (!correction.empty())
            {
                addCorrection(monomers[i], correction, false);
            }
        }
    }
}

//! Whether ACT_GPU_RESIDENT_CHECK is set, see ForceComputer::computeRelaxedEnergies
static bool gpuResidentCheck()
{
    static const bool check = std::getenv("ACT_GPU_RESIDENT_CHECK") != nullptr;
    return check;
}

bool ForceComputer::useGpuResident() const
{
    return (gpuResident_ || gpuResidentCheck()) && !gpuResidentUnavailable_;
}

bool ForceComputer::prepareGpuResident(const ForceField *pd, const std::vector<ACTMol *> &allmols) const
{
    if (!useGpuResident() || !pd->polarizable())
    {
        if (useGpuResident() && gpu::useGpu())
        {
            gpuResidentWhyNot_ = "the force field is not polarizable";
        }
        return false;
    }
    static const size_t timer = Timings::instance().index("GPU resident setup");
    ScopedTimer         scoped(timer);
    if (!gpuSystem_)
    {
        gpuSystem_ = std::make_shared<GpuTrainingSystem>();
        // The blocks of the structures of the multipoles, see setGpuMultipoles
        gpuSystem_->setMultipoleBlocks(gpuMultipoles_);
    }
    if (!gpuSystem_->prepare(pd, allmols, shellMaxDistance(pd), gpuCharges(pd), this, gpuReplicas_))
    {
        // The processor computes everything from now on
        gpuResidentUnavailable_ = true;
        gpuResidentWhyNot_      = gpuSystem_->whyNot();
        return false;
    }
    return true;
}

bool ForceComputer::gpuCharges(const ForceField *pd) const
{
    auto qdist = pd->chargeDistribution();
    return (useGpuResident() && !GMX_DOUBLE && gpu::useGpu() && pd->polarizable() &&
            pd->interactionPresent(InteractionType::BONDCORRECTIONS) &&
            pd->chargeGenerationAlgorithm() == ChargeGenerationAlgorithm::SQE &&
            (qdist == ChargeDistributionType::Gaussian || qdist == ChargeDistributionType::Point));
}

bool ForceComputer::gpuChargesOnly(const ForceField *pd) const
{
    return !gpuResidentCheck() && gpuCharges(pd);
}

void ForceComputer::computeRelaxedEnergiesGpu(MsgHandler                          *,
                                              const ForceField                    *pd,
                                              const std::vector<ACTMol *>         &allmols,
                                              bool                                 resetTheShells,
                                              bool                                 separateInductionCorrection,
                                              bool                                 checkCharges,
                                              const std::function<void(size_t)> *stage,
                                              bool                                 deviationsOnly) const
{
    static const size_t timer = Timings::instance().index("GPU resident evaluation");
    ScopedTimer         scoped(timer);
    auto &gs = *gpuSystem_;
    gs.setKeepEnergies(keepGpuEnergies_);
    // Only the deviations, the energies stay on the GPU, which does the
    // post-processing below as well, see GpuTrainingSystem::computeDeviations
    deviationsOnly = deviationsOnly && gs.deviationsReady() && !separateInductionCorrection;
    gs.setDeviationsOnly(deviationsOnly);
    // With several replicas, the work on the processor is done for each,
    // with its energies in the molecules, see GpuTrainingSystem::selectReplica
    const size_t nrep = gs.nreplica();
    auto forEachReplica = [&gs, nrep](const std::function<void(size_t)> &work)
    {
        for (size_t r = 0; r < nrep; r++)
        {
            gs.selectReplica(r);
            work(r);
        }
    };
    // Start the relaxation from the shells of the previous evaluation, see
    // setKeepShells. The energies before relaxing are still those with the
    // shells on their parents.
    const bool keepShells = keepShells_ && resetTheShells;
    gpu::RelaxSettings settings = { msForceToler_, gmx::square(shellMaxDistance(pd)), minShellStep, maxShellStep, maxiter_,
                                    stallIterations_, shellWallAngleTolerance };
    // For developers, the relaxation without asking the GPU whether it is
    // done, as in the steps of the Monte Carlo chain on the GPU
    static const bool noPoll = std::getenv("ACT_GPU_RELAX_NOPOLL") != nullptr;
    settings.poll = !noPoll;
    for (auto molptr : allmols)
    {
        auto &tops = molptr->fragmentHandler()->topologies();
        nEval_ += (stage ? nrep : 1)*molptr->conformationData(0)->nconf()*(tops.size() + (tops.size() > 1 ? 1 : 0));
    }
    static const size_t upload  = Timings::instance().index("GPU resident parameter upload");
    static const size_t stageT  = Timings::instance().index("GPU resident parameter staging");
    static const size_t relaxT  = Timings::instance().index("GPU resident relaxation");
    static const size_t energyT = Timings::instance().index("GPU resident energies");
    // With one set of parameters only the first replica is computed
    if (!stage)
    {
        for (size_t r = 1; r < nrep; r++)
        {
            gs.setActiveMolecules(r, std::vector<char>(allmols.size(), 0));
        }
    }
    // As computeRelaxedEnergies, with the coordinates on the GPU
    {
        if (stage)
        {
            {
                // The parameters of each replica, on the processor
                ScopedTimer st(stageT);
                for (size_t r = 0; r < nrep; r++)
                {
                    (*stage)(r);
                    gs.stageParameters(r, pd, this);
                }
            }
            ScopedTimer up(upload);
            gs.uploadParameters(pd);
        }
        else
        {
            {
                ScopedTimer st(stageT);
                for (size_t r = 0; r < nrep; r++)
                {
                    gs.stageParameters(r, pd, this);
                }
            }
            // As GpuTrainingSystem::updateParameters
            ScopedTimer up(upload);
            gs.uploadParameters(pd);
        }
    }
    if (resetTheShells)
    {
        gs.resetShells();
    }
    gs.constructVsites();
    // With the shells of the previous evaluation, the forces with the shells
    // on their parents are not needed, only the energies, see setKeepShells
    const bool restoreKept = keepShells && gs.haveKeptShells();
    // The energies before relaxing are backed up, to extract the induction
    { ScopedTimer et(energyT); gs.compute(-1, !restoreKept, true); }
    // The energies of computeForces are not needed, they are computed
    // after relaxing
    if (restoreKept && gs.restoreShells(-2))
    {
        ScopedTimer et(energyT);
        gs.computeForces(-2);
    }
    const bool all[2] = { true, true };
    { ScopedTimer rt(relaxT); size_t nf = 0; nShellIter_ += gs.relax(-1, all, settings, &nf); nShellConvergeFailed_ += nf; }
    // Both monomers at once
    { ScopedTimer et(energyT); gs.compute(-2, false); }
    // For each replica and molecule
    std::vector<std::vector<std::vector<double>>> eInduc2(nrep, std::vector<std::vector<double>>(allmols.size()));
    bool doInduc2 = separateInductionCorrection;
    if (doInduc2)
    {
        for (size_t r = 0; r < nrep; r++)
        {
            for (size_t k = 0; k < allmols.size(); k++)
            {
                if (allmols[k]->isDimer())
                {
                    eInduc2[r][k].resize(allmols[k]->dimerCoords()->size(), 0.0);
                }
            }
        }
        gs.saveDimers();
        gs.copyRelaxedMonomerShells();
        { ScopedTimer et(energyT); gs.compute(2); }
        std::vector<std::vector<std::vector<double>>> eRef(nrep, std::vector<std::vector<double>>(allmols.size()));
        forEachReplica([&allmols, &eRef](size_t r)
        {
            for (size_t k = 0; k < allmols.size(); k++)
            {
                eRef[r][k] = sumElecTerms(&allmols[k]->conformationData(2)->energies);
            }
        });
        for (size_t im = 0; im < 2; im++)
        {
            if (im > 0)
            {
                gs.copyRelaxedMonomerShells();
                { ScopedTimer et(energyT); gs.computeForces(2); }
            }
            const bool relaxMonomer[2] = { im == 0, im == 1 };
            { ScopedTimer rt(relaxT); size_t nf = 0; nShellIter_ += gs.relax(2, relaxMonomer, settings, &nf); nShellConvergeFailed_ += nf; }
            { ScopedTimer et(energyT); gs.compute(2, false); }
            forEachReplica([&allmols, &eRef, &eInduc2](size_t r)
            {
                for (size_t k = 0; k < allmols.size(); k++)
                {
                    auto esum = sumElecTerms(&allmols[k]->conformationData(2)->energies);
                    for (size_t i = 0; i < eInduc2[r][k].size() && i < esum.size(); i++)
                    {
                        eInduc2[r][k][i] += esum[i] - eRef[r][k][i];
                    }
                }
            });
        }
        gs.restoreDimers();
        if (restoreKept)
        {
            gs.restoreShells(2);
        }
        { ScopedTimer et(energyT); gs.computeForces(2); }
    }
    else if (restoreKept && gs.restoreShells(2))
    {
        ScopedTimer et(energyT);
        gs.computeForces(2);
    }
    { ScopedTimer rt(relaxT); size_t nf = 0; nShellIter_ += gs.relax(2, all, settings, &nf); nShellConvergeFailed_ += nf; }
    { ScopedTimer et(energyT); gs.compute(2, false); }
    if (!settings.poll)
    {
        // Whether to keep the shells depends on them, see keepShells
        size_t nf = 0;
        nShellIter_ += gs.relaxStatistics(&nf);
        nShellConvergeFailed_ += nf;
    }
    if (keepShells)
    {
        gs.keepShells();
    }
    // The deviations of the energies, if requested, see setGpuDeviations
    gs.computeDeviations();
    gs.setDeviationsOnly(false);
    // All molecules again in the next evaluation, unless told otherwise
    gs.setActiveMolecules(-1, std::vector<char>(allmols.size(), 1));
    if (!deviationsOnly)
    {
        forEachReplica([&allmols, doInduc2, separateInductionCorrection, &eInduc2](size_t r)
        {
            postProcessRelaxedEnergies(allmols, false, doInduc2, separateInductionCorrection, eInduc2[r]);
        });
    }
    if (checkCharges)
    {
        auto dq = gs.chargeDifference();
        fprintf(stderr, "GPU resident check: %-22s largest difference %.2e (of %.6g)\n",
                "charges", dq.first, dq.second);
    }
    // The energies and the charges of the first replica, the one computed
    // with one set of parameters
    if (!stage)
    {
        gs.selectReplica(0);
    }
    // The processor needs the charges too, for instance for the bounds
    gs.copyChargesToMolecules();
}

void ForceComputer::chainEvaluateGpu(const ForceField *pd, const std::vector<ACTMol *> &allmols) const
{
    // As computeRelaxedEnergiesGpu with the deviations only, for the steps
    // of the chain, see GpuTrainingSystem::chainPropose: the parameters are
    // on the GPU already, the masks of the blocks are those of the chain,
    // and nothing waits for the GPU
    auto &gs = *gpuSystem_;
    gs.setKeepEnergies(keepGpuEnergies_);
    gs.setDeviationsOnly(true);
    gpu::RelaxSettings settings = { msForceToler_, gmx::square(shellMaxDistance(pd)), minShellStep, maxShellStep, maxiter_,
                                    stallIterations_, shellWallAngleTolerance };
    // Asking the GPU every few iterations whether the relaxation is done
    // waits for it, but is faster than submitting all iterations (BFdb 6P
    // 57.6 s against 68.2 s), and a displacement that is not finite keeps
    // the shells from being kept, as in computeRelaxedEnergiesGpu. For
    // developers ACT_GPU_CHAIN_NOPOLL submits all.
    static const bool noPoll = std::getenv("ACT_GPU_CHAIN_NOPOLL") != nullptr;
    settings.poll = !noPoll;
    // Each replica computes the molecules, see GpuTrainingSystem::setChain
    for (auto molptr : allmols)
    {
        auto &tops = molptr->fragmentHandler()->topologies();
        nEval_ += gs.nreplica()*molptr->conformationData(0)->nconf()*(tops.size() + (tops.size() > 1 ? 1 : 0));
    }
    const bool keepShells = keepShells_;
    gs.resetShells();
    gs.constructVsites();
    const bool restoreKept = keepShells && gs.haveKeptShells();
    gs.compute(-1, !restoreKept, true);
    if (restoreKept && gs.restoreShells(-2))
    {
        gs.computeForces(-2);
    }
    const bool all[2] = { true, true };
    size_t nf = 0;
    nShellIter_ += gs.relax(-1, all, settings, &nf);
    nShellConvergeFailed_ += nf;
    gs.compute(-2, false);
    if (restoreKept && gs.restoreShells(2))
    {
        gs.computeForces(2);
    }
    nShellIter_ += gs.relax(2, all, settings, &nf);
    nShellConvergeFailed_ += nf;
    gs.compute(2, false);
    if (keepShells)
    {
        gs.keepShells();
    }
    gs.computeDeviations();
    gs.setDeviationsOnly(false);
}

void ForceComputer::chainStatisticsGpu() const
{
    size_t nf = 0;
    nShellIter_           += gpuSystem_->relaxStatistics(&nf);
    nShellConvergeFailed_ += nf;
}

bool ForceComputer::computeRelaxedEnergiesReplicas(MsgHandler                          *msghandler,
                                                   const ForceField                    *pd,
                                                   const std::vector<ACTMol *>         &allmols,
                                                   const std::function<void(size_t)> &stage,
                                                   bool                                 separateInductionCorrection,
                                                   bool                                 deviationsOnly) const
{
    if (!prepareGpuResident(pd, allmols))
    {
        return false;
    }
    computeRelaxedEnergiesGpu(msghandler, pd, allmols, true, separateInductionCorrection, false, &stage,
                              deviationsOnly);
    return true;
}

bool ForceComputer::computeGpuDeviations(MsgHandler                  *msghandler,
                                         const ForceField            *pd,
                                         const std::vector<ACTMol *> &allmols) const
{
    if (!prepareGpuResident(pd, allmols) || !gpuSystem_->deviationsReady())
    {
        return false;
    }
    computeRelaxedEnergiesGpu(msghandler, pd, allmols, true, false, false, nullptr, true);
    return true;
}

void ForceComputer::selectGpuReplica(size_t replica) const
{
    if (gpuSystem_)
    {
        gpuSystem_->selectReplica(replica);
    }
}

void ForceComputer::computeRelaxedEnergies(MsgHandler                  *msghandler,
                                           const ForceField            *pd,
                                           const std::vector<ACTMol *> &allmols,
                                           bool                         resetTheShells,
                                           bool                         separateInductionCorrection,
                                           bool                         computeForces) const
{
    // All on the GPU, if requested and possible, for now without forces,
    // see setGpuResident. For developers, ACT_GPU_RESIDENT_CHECK compares
    // the energies of the first evaluations with those of the processor.
    static int        checks   = gpuResidentCheck() ? 5 : 0;
    static thread_local bool onProcessor = false;
    std::vector<ConformationEnergies> gpuEnergies;
    if (useGpuResident() && !onProcessor && !computeForces && pd->polarizable())
    {
        bool chargesOnGpu = gpuChargesOnly(pd);
        if (prepareGpuResident(pd, allmols))
        {
            computeRelaxedEnergiesGpu(msghandler, pd, allmols, resetTheShells, separateInductionCorrection,
                                      checks > 0 && gpuCharges(pd));
            if (checks == 0)
            {
                return;
            }
            checks -= 1;
            for (auto molptr : allmols)
            {
                for (size_t what = 0; what < 3; what++)
                {
                    gpuEnergies.push_back(molptr->conformationData(what)->energies);
                }
            }
        }
        else if (chargesOnGpu)
        {
            // The charges were left to the GPU, but ACMFitnessComputer::prepareDeviation
            // calls prepareGpuResident before, so this should not happen
            msghandler->fatal("The GPU cannot compute these molecules, but was to generate their charges. Use train_ff -nogpu_resident.");
        }
    }
    if (!gpuEnergies.empty())
    {
        // Compute on the processor and compare at the end
        onProcessor = true;
        computeRelaxedEnergies(msghandler, pd, allmols, resetTheShells, separateInductionCorrection, computeForces);
        onProcessor = false;
        std::map<InteractionType, std::pair<double, double>> worst;
        size_t g = 0;
        for (auto molptr : allmols)
        {
            for (size_t what = 0; what < 3; what++, g++)
            {
                auto       &cpu = molptr->conformationData(what)->energies;
                auto       &gpu = gpuEnergies[g];
                for (auto itype : cpu.types())
                {
                    auto ec = cpu.find(itype);
                    auto eg = gpu.find(itype);
                    if (!ec.have || !eg.have)
                    {
                        continue;
                    }
                    for (size_t c = 0; c < cpu.nconf(); c++)
                    {
                        if (ec.has(c) && eg.has(c))
                        {
                            double rel = std::abs(ec.act[c] - eg.act[c])/(1 + std::abs(ec.act[c]));
                            auto  &w   = worst[itype];
                            if (rel > w.first)
                            {
                                w = { rel, ec.act[c] };
                            }
                        }
                    }
                }
            }
        }
        for (const auto &w : worst)
        {
            fprintf(stderr, "GPU resident check: %-22s largest relative difference %.2e (of %.6g)\n",
                    interactionTypeToString(w.first).c_str(), w.second.first, w.second.second);
        }
        return;
    }
    static const size_t timer = Timings::instance().index("Energy calculation overhead");
    ScopedTimer         scoped(timer);
    // Shells are relaxed using the forces
    bool forces = computeForces || pd->polarizable();
    // Without polarization, the interaction energies of dimers are those of
    // the intermolecular interactions, the monomers are not needed then.
    // The forces would lack the intramolecular contributions.
    bool intermolecular = !forces;
    // Loop over molecules
    for(auto molptr : allmols)
    {
        auto &mol = *molptr;
        bool  inter = intermolecular && mol.isDimer() && mol.topology()->hasIntermolecular();
        if (mol.isDimer())
        {
            // We have a dimer!
            auto top    = mol.topology();
            auto cdata  = mol.conformationData(2);
            bool vsites = vsitesOutdated(top, cdata);
            // The energy calculation constructs them if possible, and
            // spreads their forces
            cdata->constructVsites = vsites && vsitesIndependentOfShells(top);
            cdata->spreadVsites    = computeForces;
            vsites                 = vsites && !cdata->constructVsites;
            // Only the shells, if any, are reset
            std::vector<std::pair<size_t, size_t>> shells;
            if (resetTheShells)
            {
                shells = shellParents(top->atoms());
            }
            // Virtual sites for dimer
            for(size_t c = 0; c < cdata->nconf(); c++)
            {
                nEval_ += 1;
                if (vsites)
                {
                    constructVsites(vsiteHandler_, top, cdata->x(c));
                }
                // Reset shells if needed and if possible
                resetShells(shells, cdata->x(c));
            }
        }
        // Shortcuts...
        auto &tops  = mol.fragmentHandler()->topologies();
        // Virtual sites for monomers
        for(size_t index = 0; index < tops.size() && !inter; index++)
        {
            auto &top    = tops[index];
            auto  cdata  = mol.conformationData(index);
            bool  vsites = vsitesOutdated(&top, cdata);
            cdata->constructVsites = vsites && vsitesIndependentOfShells(&top);
            cdata->spreadVsites    = computeForces;
            vsites                 = vsites && !cdata->constructVsites;
            std::vector<std::pair<size_t, size_t>> shells;
            if (resetTheShells)
            {
                shells = shellParents(top.atoms());
            }
            for(size_t c = 0; c < cdata->nconf(); c++)
            {
                nEval_ += 1;
                if (vsites)
                {
                    constructVsites(vsiteHandler_, &top, cdata->x(c));
                }
                // Reset shells if needed and if possible
                resetShells(shells, cdata->x(c));
            }
        }
    }
    // Without forces in the result, only those on the shells are needed
    // for relaxing them
    bool allForces = computeForces;
    // Do first calculation every time.
    computeEnergies(msghandler, pd, allmols, 0, forces, false, intermolecular, allForces); // Monomer 0
    computeEnergies(msghandler, pd, allmols, 1, forces, false, intermolecular, allForces); // Monomer 1
    computeEnergies(msghandler, pd, allmols, 2, forces, false, intermolecular, allForces); // Complete topology
    // For developers, compare the first calculation with that of the GPU
    static int gpuChecks = std::getenv("ACT_GPU_SYSTEM_CHECK") ? 3 : 0;
    if (gpuChecks > 0)
    {
        gpuChecks -= 1;
        (void) checkGpuSystem(pd, allmols);
    }

    // Store the original (unrelaxed) energies so that we can extract the
    // induction later. Only terms that actually have an ACT energy are backed up.
    // Only if shells are relaxed the energies change.
    bool relax = pd->polarizable();
    for(auto molptr : allmols)
    {
        auto &mol = *molptr;
        for(size_t what = 0; what < 3; what++)
        {
            auto &energies = mol.conformationData(what)->energies;
            energies.setBackups(relax);
            if (relax)
            {
                backupEnergies(&energies);
            }
        }
    }
    // Now let's have a look whether we are polarizable
    if (pd->polarizable())
    {
        // Start from the shells of the previous relaxation, see setKeepShells
        if (keepShells_ && resetTheShells && restoreRelaxedShells(allmols, -2))
        {
            computeEnergies(msghandler, pd, allmols, 0, forces, false, intermolecular, allForces);
            computeEnergies(msghandler, pd, allmols, 1, forces, false, intermolecular, allForces);
        }
        for(size_t im = 0; im < 2; im++)
        {
            relaxMonomer(msghandler, pd, im, allmols, allForces);
        }
        // Backup the now-relaxed monomer coordinates so that we can return to
        // them later when computing the second order induction.
        for(auto molptr : allmols)
        {
            auto &mol = *molptr;
            // Shortcuts...
            auto &tops  = mol.fragmentHandler()->topologies();
            for(size_t index = 0; index < tops.size(); index++)
            {
                auto *cdata = mol.conformationData(index);
                for(size_t c = 0; c < cdata->nconf(); c++)
                {
                    cdata->confs[c].xrelax = cdata->x(c).toRVec();
                }
            }
        }
    }
    // Second order induction energy for each dimer conformation of each molecule
    std::vector<std::vector<double>> eInduc2(allmols.size());
    bool doInduc2 = pd->polarizable() && separateInductionCorrection;
    // From now on, only polarizable models relax the dimer further, but the
    // post-processing below runs in all cases.
    if (pd->polarizable())
    {
        if (doInduc2)
        {
            // Store the original dimer coordinates, such that the final
            // relaxation below does not depend on what we do here.
            std::vector<std::vector<std::vector<gmx::RVec>>> xdimer(allmols.size());
            for(size_t k = 0; k < allmols.size(); k++)
            {
                auto &mol = (*allmols[k]);
                if (mol.isDimer())
                {
                    auto *cdata = mol.conformationData(2);
                    for(size_t c = 0; c < cdata->nconf(); c++)
                    {
                        xdimer[k].push_back(cdata->x(c).toRVec());
                    }
                    eInduc2[k].resize(mol.dimerCoords()->size(), 0.0);
                }
            }
            // We have to compute second order induction first, starting
            // from the shells in the relaxed monomer positions.
            copyRelaxedMonomerShells(allmols);
            computeEnergies(msghandler, pd, allmols, 2, true, false, false, allForces);
            // Reference energy without mutual induction
            std::vector<std::vector<double>> eRef(allmols.size());
            for(size_t k = 0; k < allmols.size(); k++)
            {
                eRef[k] = sumElecTerms(&(*allmols[k]).conformationData(2)->energies);
            }
            for(size_t im = 0; im < 2; im++)
            {
                if (im > 0)
                {
                    // Reset coordinates to the relaxed monomer shells
                    copyRelaxedMonomerShells(allmols);
                    computeEnergies(msghandler, pd, allmols, 2, true, false, false, allForces);
                }
                // Relax monomer im in the presence of the other monomer
                bool fixMonomer[2] = { im != 0, im != 1 };
                relaxDimer(msghandler, pd, fixMonomer, allmols, allForces);
                // Extract energy
                for(size_t k = 0; k < allmols.size(); k++)
                {
                    auto esum = sumElecTerms(&(*allmols[k]).conformationData(2)->energies);
                    for(size_t i = 0; i < eInduc2[k].size() && i < esum.size(); i++)
                    {
                        eInduc2[k][i] += esum[i] - eRef[k][i];
                    }
                }
            }
            // Restore the original dimer coordinates and energies
            for(size_t k = 0; k < allmols.size(); k++)
            {
                auto *cdata = (*allmols[k]).conformationData(2);
                for(size_t i = 0; i < xdimer[k].size(); i++)
                {
                    cdata->x(i).assign(xdimer[k][i]);
                }
            }
            if (keepShells_ && resetTheShells)
            {
                (void) restoreRelaxedShells(allmols, 2);
            }
            computeEnergies(msghandler, pd, allmols, 2, true, false, false, allForces);
        }
        else if (keepShells_ && resetTheShells && restoreRelaxedShells(allmols, 2))
        {
            computeEnergies(msghandler, pd, allmols, 2, true, false, false, allForces);
        }
        // Final, total energy
        bool ff[2] = { false, false };
        relaxDimer(msghandler, pd, ff, allmols, allForces);
        if (keepShells_ && resetTheShells)
        {
            keepRelaxedShells(allmols);
        }
    }
    
    postProcessRelaxedEnergies(allmols, intermolecular, doInduc2, separateInductionCorrection, eInduc2);
    if (!allForces)
    {
        return;
    }
    // Spread forces to atoms, unless the energy calculations did
    // Loop over molecules
    for(auto molptr : allmols)
    {
        auto &mol = *molptr;
        // Shortcuts...
        auto &tops  = mol.fragmentHandler()->topologies();
        if (tops.size() == 2)
        {
            // We have a dimer!
            auto top    = mol.topology();
            auto cdata  = mol.conformationData(2);
            // Virtual sites for dimer
            for(size_t c = 0; c < cdata->nconf(); c++)
            {
                if (!cdata->vsiteForcesSpread)
                {
                    spreadVsites(vsiteHandler_, top, cdata->x(c), cdata->f(c));
                }
            }
            cdata->spreadVsites = false;
        }
        // Virtual sites for monomers
        for(size_t index = 0; index < tops.size(); index++)
        {
            auto &top   = tops[index];
            auto  cdata = mol.conformationData(index);
            for(size_t c = 0; c < cdata->nconf(); c++)
            {
                if (!cdata->vsiteForcesSpread)
                {
                    spreadVsites(vsiteHandler_, &top, cdata->x(c), cdata->f(c));
                }
            }
            cdata->spreadVsites = false;
        }
    }
}

/*! \brief Add the forces of an electric field on the charges
 * \param[in]    atoms  The particles
 * \param[in]    field  The electric field
 * \param[inout] cdata  The conformations, only those of which the last
 *                      calculation computed the forces are changed, see
 *                      ConformationArrays::computed
 */
static void addFieldForces(const std::vector<ActAtom>  &atoms,
                           const gmx::RVec             &field,
                           ConformationData            *cdata)
{
    if (0 == field[XX] && 0 == field[YY] && 0 == field[ZZ])
    {
        return;
    }
    for(size_t c = 0; c < cdata->nconf(); c++)
    {
        if (!cdata->arrays.computed(c))
        {
            continue;
        }
        auto f = cdata->f(c);
        for(size_t i = 0; i < atoms.size() && i < f.size(); i++)
        {
            real fac = FIELDFAC*atoms[i].charge();
            for(int m = 0; m < DIM; m++)
            {
                f[i][m] += fac*field[m];
            }
        }
    }
}

void ForceComputer::computeTotal(MsgHandler                                     *msghandler,
                                 const ForceField                               *pd,
                                 ConformationSet                                *set,
                                 const ShellOptions                             &options,
                                 std::vector<std::map<InteractionType, double>> *energies) const
{
    double maxShellDistance = shellMaxDistance(pd);
    auto  top   = set->top;
    auto &confs = set->data.confs;
    auto  cdata = &set->data;
    auto &atoms = top->atoms();
    // Construct the virtual sites and reset the shells. The energy
    // calculation constructs the virtual sites if possible, and spreads
    // their forces unless forces of the field are added afterwards.
    bool vsites = options.virtualSites && vsitesOutdated(top, cdata);
    bool noField = (0 == options.field[XX] && 0 == options.field[YY] && 0 == options.field[ZZ]);
    cdata->constructVsites = vsites && vsitesIndependentOfShells(top);
    cdata->spreadVsites    = options.virtualSites && noField;
    vsites                 = vsites && !cdata->constructVsites;
    std::vector<std::pair<size_t, size_t>> shells;
    if (options.resetShells)
    {
        shells = shellParents(atoms);
    }
    // The input positions of the shells, to start the relaxation from
    std::vector<std::vector<gmx::RVec>> inputShells;
    for(size_t c = 0; c < cdata->nconf(); c++)
    {
        auto x  = cdata->x(c);
        nEval_ += 1;
        if (vsites)
        {
            constructVsites(vsiteHandler_, top, x);
        }
        if (options.shellsFromInput && !shells.empty())
        {
            inputShells.emplace_back();
            for(const auto &sp : shells)
            {
                inputShells.back().push_back(x[sp.first].rvec());
            }
        }
        resetShells(shells, x);
    }
    computeConformations(msghandler, pd, top, cdata, true, false);
    addFieldForces(atoms, options.field, cdata);
    // The energies before relaxing the shells give the induction energy
    auto &cener = cdata->energies;
    auto  itype = InteractionType::POLARIZATION;
    bool  relax = options.relax && pd->polarizable() && top->hasEntry(itype);
    cener.setBackups(relax);
    if (relax)
    {
        backupEnergies(&cener);
        // The shells to relax
        std::vector<bool> isShell(atoms.size(), false);
        size_t            nShell = 0;
        for(size_t i = 0; i < atoms.size(); i++)
        {
            if (atoms[i].pType() == ActParticle::Shell &&
                (options.relaxShells.empty() ||
                 options.relaxShells.end() != options.relaxShells.find(i)))
            {
                isShell[i] = true;
                nShell    += 1;
            }
        }
        resetShellSteps(cdata);
        if (!inputShells.empty())
        {
            // Start from the input positions of the shells to relax, which
            // takes an energy calculation like an iteration
            for(size_t c = 0; c < confs.size(); c++)
            {
                for(size_t k = 0; k < shells.size(); k++)
                {
                    if (isShell[shells[k].first])
                    {
                        for(int m = 0; m < DIM; m++)
                        {
                            cdata->x(c)[shells[k].first][m] = inputShells[c][k][m];
                        }
                    }
                }
                nShellIter_ += 1;
            }
            computeConformations(msghandler, pd, top, cdata, true, true);
            addFieldForces(atoms, options.field, cdata);
        }
        // We have already done one energy calculation
        int  iter    = 1;
        bool allDone = false;
        while (iter < maxiter_ && !allDone)
        {
            allDone = (updateShells(top, 0, nShell, isShell, 0, maxShellDistance, cdata) == confs.size());
            if (!allDone)
            {
                // Only for the conformations that did not converge
                for(const auto &conf : confs)
                {
                    nShellIter_ += conf.done ? 0 : 1;
                }
                computeConformations(msghandler, pd, top, cdata, true, true);
                addFieldForces(atoms, options.field, cdata);
            }
            iter++;
        }
        // Check the conformations computed last
        double maxShellDistance2 = maxShellDistance*maxShellDistance;
        for(size_t c = 0; c < confs.size(); c++)
        {
            auto &conf = confs[c];
            if (conf.done)
            {
                nShellConvergeFailed_ += shellsStalled(cdata, c, 0, stallIterations_) ? 1 : 0;
                continue;
            }
            double msForce = shellForce2(atoms, isShell, cdata->x(c), cdata->f(c), 0, maxShellDistance2);
            if (nShell > 0)
            {
                msForce /= nShell;
            }
            if (msForce > msForceToler_)
            {
                nShellConvergeFailed_ += 1;
            }
            conf.done = true;
        }
    }
    // The energy terms, as in compute
    energies->clear();
    for(size_t c = 0; c < confs.size(); c++)
    {
        std::map<InteractionType, double> emap;
        for(auto t : cener.types())
        {
            auto et = cener.find(t);
            if (et.has(c))
            {
                emap.insert({ t, et.act[c] });
            }
        }
        if (!options.energyTerms)
        {
            energies->push_back(emap);
            continue;
        }
        // The change of the electrostatic terms due to the relaxation
        // of the shells is the induction energy.
        bool   haveInduction = false;
        double eInduction    = 0;
        for(auto et : { InteractionType::ELECTROSTATICS,
                        InteractionType::POLARIZATION,
                        InteractionType::CHARGETRANSFER })
        {
            auto tt = emap.find(et);
            if (emap.end() != tt)
            {
                double before = relax ? cener.find(et).before[c] : tt->second;
                eInduction   += tt->second - before;
                tt->second    = before;
                haveInduction = true;
            }
        }
        if (haveInduction)
        {
            emap.insert_or_assign(InteractionType::INDUCTION, eInduction);
        }
        auto sum = [&emap](std::initializer_list<InteractionType> itypes)
        {
            double esum = 0;
            for(auto it : itypes)
            {
                auto ee = emap.find(it);
                if (emap.end() != ee)
                {
                    esum += ee->second;
                }
            }
            return esum;
        };
        emap.insert_or_assign(InteractionType::ALLELEC,
                              sum({ InteractionType::ELECTROSTATICS,
                                    InteractionType::POLARIZATION,
                                    InteractionType::INDUCTION,
                                    InteractionType::INDUCTIONCORRECTION }));
        emap.insert_or_assign(InteractionType::EXCHIND,
                              sum({ InteractionType::EXCHANGE,
                                    InteractionType::INDUCTION,
                                    InteractionType::INDUCTIONCORRECTION }));
        energies->push_back(emap);
    }
    // Spread forces to atoms, unless the energy calculations did
    for(size_t c = 0; c < confs.size(); c++)
    {
        if (options.virtualSites && !(cdata->spreadVsites && cdata->vsiteForcesSpread))
        {
            spreadVsites(vsiteHandler_, top, cdata->x(c), cdata->f(c));
        }
    }
    cdata->spreadVsites = false;
}

void ForceComputer::computeDimerConformations(MsgHandler                                *msghandler,
                                              const ForceField                          *pd,
                                              ACTMol                                    *mol,
                                              const std::vector<std::vector<gmx::RVec>> &coords,
                                              bool                                       separateInductionCorrection) const
{
    // Replace the conformations of the dimer and its monomers
    for(size_t what = 0; what < 3; what++)
    {
        mol->conformationData(what)->invalidate();
    }
    size_t natom1 = mol->fragmentHandler()->atomStart()[1];
    size_t nconf  = coords.size();
    // The conformations of the previous call are used again, such that their
    // memory does not have to be allocated anew, but they are reset as new ones.
    auto reset = [](ConformationData                       *cdata,
                    size_t                                  c,
                    std::vector<gmx::RVec>::const_iterator  begin,
                    std::vector<gmx::RVec>::const_iterator  end)
    {
        auto *ce = &cdata->confs[c];
        ce->done = false;
        ce->dataFile.clear();
        cdata->setX(c, begin, end);
        ce->xrelax.clear();
        ce->fobs.clear();
        ce->eMap.clear();
        ce->eMin = 0;
    };
    for(size_t what = 0; what < 3; what++)
    {
        mol->conformationData(what)->confs.resize(nconf);
    }
    for(size_t c = 0; c < nconf; c++)
    {
        const auto &x = coords[c];
        reset(mol->conformationData(2), c, x.begin(), x.end());
        reset(mol->conformationData(0), c, x.begin(), x.begin()+natom1);
        reset(mol->conformationData(1), c, x.begin()+natom1, x.end());
        for(size_t k = 0; k < 2; k++)
        {
            auto *cdata = mol->conformationData(k);
            cdata->confs[c].xrelax = cdata->x(c).toRVec();
        }
    }
    computeRelaxedEnergies(msghandler, pd, { mol }, true, separateInductionCorrection, true);
}

void ForceComputer::computeInteraction(MsgHandler                                     *msghandler,
                                       const ForceField                               *pd,
                                       ACTMol                                         *mol,
                                       std::vector<std::vector<gmx::RVec>>            *coords,
                                       bool                                            separateInductionCorrection,
                                       std::vector<std::map<InteractionType, double>> *einter,
                                       std::vector<double>                            *epot,
                                       std::vector<std::vector<gmx::RVec>>            *forces) const
{
    if (!mol->isDimer())
    {
        einter->clear();
        epot->clear();
        forces->clear();
        return;
    }
    computeDimerConformations(msghandler, pd, mol, *coords, separateInductionCorrection);
    const ConformationData *cdatas[3] = { mol->conformationData(0), mol->conformationData(1),
                                          mol->conformationData(2) };
    size_t natom1 = mol->fragmentHandler()->atomStart()[1];
    size_t nconf  = coords->size();
    // The energies are taken from the conformation data, see ConformationData::energies
    auto &dimerEnergies = mol->conformationData(2)->energies;
    auto  dimerEpot     = dimerEnergies.find(InteractionType::EPOT);
    ConformationEnergies::Term monomerEpot[2] = {
        mol->conformationData(0)->energies.find(InteractionType::EPOT),
        mol->conformationData(1)->energies.find(InteractionType::EPOT)
    };
    size_t nmonomer[2] = { mol->conformationData(0)->energies.nconf(),
                           mol->conformationData(1)->energies.nconf() };
    // The polarization is part of the induction
    std::vector<std::pair<InteractionType, ConformationEnergies::Term>> terms;
    for(auto t : dimerEnergies.types())
    {
        if (InteractionType::POLARIZATION != t)
        {
            terms.push_back({ t, dimerEnergies.find(t) });
        }
    }
    // The output of the previous call is used again as well
    einter->resize(nconf);
    epot->resize(nconf);
    forces->resize(nconf);
    for(size_t c = 0; c < nconf; c++)
    {
        bool  inDimer = c < dimerEnergies.nconf();
        auto &emap    = (*einter)[c];
        // Remove the terms that are not there anymore, without
        // removing those that are only updated
        for(auto it = emap.begin(); it != emap.end(); )
        {
            auto t = dimerEnergies.find(it->first);
            if (InteractionType::POLARIZATION == it->first || !inDimer || !t.has(c))
            {
                it = emap.erase(it);
            }
            else
            {
                ++it;
            }
        }
        for(const auto &term : terms)
        {
            if (inDimer && term.second.has(c))
            {
                emap.insert_or_assign(term.first, term.second.act[c]);
            }
        }
        // The dimer energy is the interaction energy plus the energies of
        // the relaxed monomers, and the interaction forces are the dimer
        // forces minus the forces in the monomers.
        double e = (inDimer && dimerEpot.has(c)) ? dimerEpot.act[c] : 0;
        auto  &f = (*forces)[c];
        f = cdatas[2]->f(c).toRVec();
        for(size_t k = 0; k < 2; k++)
        {
            auto fmon = cdatas[k]->f(c);
            if (c < nmonomer[k] && monomerEpot[k].has(c))
            {
                e += monomerEpot[k].act[c];
            }
            size_t offset = (k == 0) ? 0 : natom1;
            for(size_t i = 0; i < fmon.size(); i++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    f[offset+i][m] -= fmon[i][m];
                }
            }
        }
        (*coords)[c] = cdatas[2]->x(c).toRVec();
        (*epot)[c]   = e;
    }
    // The monomer energies and the monomer forces subtracted, see Timings::addFlops
    Timings::instance().addFlops(nconf*(2 + DIM*(cdatas[0]->nconf() == 0 ? 0 : cdatas[0]->f(0).size() +
                                                  cdatas[1]->f(0).size())));
}

void ForceComputer::computeInteraction(MsgHandler                                *msghandler,
                                       const ForceField                          *pd,
                                       ACTMol                                    *mol,
                                       const std::vector<std::vector<gmx::RVec>> &coords,
                                       bool                                       separateInductionCorrection,
                                       DimerInteractions                         *result) const
{
    size_t nconf = mol->isDimer() ? coords.size() : 0;
    result->nconf  = nconf;
    result->stride = nconf;
    result->einter.resize(nconf);
    result->epot.resize(nconf);
    if (0 == nconf)
    {
        result->x.clear();
        result->f.clear();
        return;
    }
    computeDimerConformations(msghandler, pd, mol, coords, separateInductionCorrection);
    // The energies are taken from the conformation data, see the other computeInteraction
    auto  dimerEnergies = &mol->conformationData(2)->energies;
    auto  dimerEpot     = dimerEnergies->find(InteractionType::EPOT);
    ConformationEnergies::Term monomerEpot[2] = {
        mol->conformationData(0)->energies.find(InteractionType::EPOT),
        mol->conformationData(1)->energies.find(InteractionType::EPOT)
    };
    size_t nmonomer[2] = { mol->conformationData(0)->energies.nconf(),
                           mol->conformationData(1)->energies.nconf() };
    for(size_t c = 0; c < nconf; c++)
    {
        double e = (c < dimerEnergies->nconf() && dimerEpot.has(c)) ? dimerEpot.act[c] : 0;
        result->einter[c] = e;
        for(size_t k = 0; k < 2; k++)
        {
            if (c < nmonomer[k] && monomerEpot[k].has(c))
            {
                e += monomerEpot[k].act[c];
            }
        }
        result->epot[c] = e;
    }
    // The interaction forces are the dimer forces minus the forces in the monomers
    const ConformationData *cdatas[3] = { mol->conformationData(0), mol->conformationData(1),
                                          mol->conformationData(2) };
    auto  &astart  = mol->fragmentHandler()->atomStart();
    size_t natom   = cdatas[2]->x(0).size();
    size_t nmon[2] = { cdatas[0]->x(0).size(), cdatas[1]->x(0).size() };
    const ConformationArrays<real> *arrays[3] = { &mol->conformationData(0)->arrays,
                                                  &mol->conformationData(1)->arrays,
                                                  &mol->conformationData(2)->arrays };
    bool inArrays = true;
    for(size_t what = 0; what < 3; what++)
    {
        inArrays = (inArrays && arrays[what]->nconf == nconf && arrays[what]->stride == arrays[2]->stride &&
                    arrays[what]->inOrder());
    }
    if (inArrays)
    {
        // The rows of the arrays of the monomers are those of their
        // particles in the dimer
        size_t stride  = arrays[2]->stride;
        result->stride = stride;
        // The integer coordinates in nm, see coordinate.h
        result->x.resize(DIM*natom*stride);
        std::transform(arrays[2]->x.begin(), arrays[2]->x.begin() + DIM*natom*stride, result->x.begin(),
                       [](Coordinate c) { return fromCoordinate(c); });
        result->f.assign(arrays[2]->f.begin(), arrays[2]->f.begin() + DIM*natom*stride);
        for(size_t k = 0; k < 2; k++)
        {
            const real   *fm = arrays[k]->f.data();
            double       *fd = result->f.data() + DIM*astart[k]*stride;
            for(size_t i = 0; i < DIM*nmon[k]*stride; i++)
            {
                fd[i] -= fm[i];
            }
        }
    }
    else
    {
        result->x.resize(DIM*natom*nconf);
        result->f.resize(DIM*natom*nconf);
        for(size_t c = 0; c < nconf; c++)
        {
            auto xdim = cdatas[2]->x(c);
            auto fdim = cdatas[2]->f(c);
            for(size_t a = 0; a < natom; a++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    result->x[(DIM*a+m)*nconf + c] = xdim[a][m];
                    result->f[(DIM*a+m)*nconf + c] = fdim[a][m];
                }
            }
            for(size_t k = 0; k < 2; k++)
            {
                auto fm = cdatas[k]->f(c);
                for(size_t a = 0; a < fm.size(); a++)
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        result->f[(DIM*(astart[k]+a)+m)*nconf + c] -= fm[a][m];
                    }
                }
            }
        }
    }
    // The monomer energies and the monomer forces subtracted, see Timings::addFlops
    Timings::instance().addFlops(nconf*(2 + DIM*(nmon[0] + nmon[1])));
}

void ForceComputer::computeOnce(MsgHandler                        *msg_handler,
                                const ForceField                  *pd,
                                const Topology                    *top,
                                std::vector<gmx::RVec>            *coordinates,
                                std::vector<gmx::RVec>            *forces,
                                std::map<InteractionType, double> *energies,
                                const gmx::RVec                   &field) const
{
    Coordinates x(*coordinates);
    computeOnce(msg_handler, pd, top, &x, forces, energies, field);
}

void ForceComputer::computeOnce(MsgHandler                        *msg_handler,
                                const ForceField                  *pd,
                                const Topology                    *top,
                                const Coordinates                 *coordinates,
                                std::vector<gmx::RVec>            *forces,
                                std::map<InteractionType, double> *energies,
                                const gmx::RVec                   &field) const
{
    // One conformation, computed with the SIMD routines, the shells
    // and the virtual sites stay where they are
    ConformationSet set;
    set.top = top;
    set.setCoordinates(std::vector<Coordinates>{ *coordinates });
    ShellOptions options;
    options.relax        = false;
    options.resetShells  = false;
    options.field        = field;
    options.energyTerms  = false;
    options.virtualSites = false;
    std::vector<std::map<InteractionType, double>> emaps;
    computeTotal(msg_handler, pd, &set, options, &emaps);
    // The evaluations are counted by the callers of computeOnce
    nEval_      -= 1;
    *forces      = set.data.f(0).toRVec();
    *energies    = emaps[0];
}

Potential ForceComputer::ftype(const ForceField *pd,
                               InteractionType   itype) const
{
    return pd->potential(itype);
}

/*! \brief Compute the energy of one potential for one structure
 * \param[in] msghandler For warnings and errors
 * \param[in] pot        The potential
 * \param[in] entries    The interactions
 * \param[in] atoms      The atoms
 * \param[in] x          The coordinates
 * \param[in] epsilonr   Relative dielectric constant
 * \return The energy
 */
static double potentialEnergy(MsgHandler                   *msghandler,
                              Potential                     pot,
                              const TopologyEntryVector    &entries,
                              const std::vector<ActAtom>   &atoms,
                              const std::vector<gmx::RVec> &x,
                              double                        epsilonr)
{
    ConformationData cdata;
    cdata.confs.resize(1);
    cdata.setX(0, x);
    computePotential(msghandler, pot, entries, atoms, &cdata, epsilonr, false);
    auto epot = cdata.energies.find(InteractionType::EPOT);
    return epot.has(0) ? epot.act[0] : 0;
}

void ForceComputer::plot(MsgHandler        *msghandler,
                         const ForceField  *pd,
                         InteractionType    itype) const
{
    // The widths of the charge distributions are those of the particle types
    bool coulomb = (InteractionType::ELECTROSTATICS == itype);
    if (!coulomb && !pd->interactionPresent(itype))
    {
        fprintf(stderr, "No such interaction %s in the force field.\n",
                interactionTypeToString(itype).c_str());
        return;
    }
    std::string btype("bondtype");

    std::map<InteractionType, const std::string> i2s = {
        { InteractionType::BONDS,              btype },
        { InteractionType::ANGLES,             btype },
        { InteractionType::LINEAR_ANGLES,      btype },
        { InteractionType::PROPER_DIHEDRALS,   btype },
        { InteractionType::IMPROPER_DIHEDRALS, btype },
        { InteractionType::VDW,                "vdwtype" },
        { InteractionType::ELECTROSTATICS,     "acmtype" }
    };
    // The interactions to plot, of trained parameters, or for COULOMB of
    // all pairs of particle types with a charge distribution
    std::vector<Identifier> ids;
    if (coulomb)
    {
        std::vector<std::string> ptypes;
        for (const auto &pt : pd->particleTypesConst())
        {
            if (pt.second.apType() == ActParticle::Atom && pt.second.hasParameter(zetaName))
            {
                ptypes.push_back(pt.first.id());
            }
        }
        for (size_t i = 0; i < ptypes.size(); i++)
        {
            for (size_t j = i; j < ptypes.size(); j++)
            {
                ids.push_back(Identifier({ ptypes[i], ptypes[j] }, { 1 }, CanSwap::Yes));
            }
        }
    }
    else
    {
        for(const auto &f : pd->findForcesConst(itype).parametersConst())
        {
            unsigned int ntrain  = 0;
            for(const auto &pp : f.second)
            {
                ntrain = std::max(ntrain, pp.second.ntrain());
            }
            if (ntrain > 0)
            {
                ids.push_back(f.first);
            }
        }
    }
    auto potential = pd->potential(itype);
    {
        std::vector<std::pair<int, int>> linear_bonds = {
            { 0, 1}, {1, 2}, {2, 3}
        };
        std::vector<std::pair<int, int>> idih_bonds = {
            { 0, 1}, {0, 2}, {0, 3}
        };
        for(const auto &id : ids)
        {
            auto bos   = id.bondOrders();
            std::vector<Bond> bbb;
            for(size_t i = 0; i < bos.size(); i++)
            {
                if (InteractionType::IMPROPER_DIHEDRALS == itype)
                {
                    bbb.push_back(Bond(idih_bonds[i].first, idih_bonds[i].second, bos[i]));
                }
                else if (InteractionType::VDW == itype ||
                         InteractionType::ELECTROSTATICS == itype)
                {
                    ;
                }
                else
                {
                    bbb.push_back(Bond(linear_bonds[i].first, linear_bonds[i].second, bos[i]));
                }
            }
            Topology top;
            top.init(bbb);

            auto subtype = i2s.find(itype);
            if (i2s.end() != subtype)
            {
                auto ptypes = pd->particleTypesConst();
                // Loop over subtypes (depending on the InteractionType)
                for(const auto &atomsub : id.atoms())
                {
                    for (const auto &pt : ptypes)
                    {
                        if ((coulomb && pt.first.id() == atomsub) ||
                            (pt.second.hasOption(subtype->second) &&
                             pt.second.optionValue(subtype->second) == atomsub))
                        {
                            ActAtom a(pt.second);
                            a.setCharge(1);
                            top.addAtom(a);
                            break;
                        }
                    }
                }
            }
            if (top.nAtoms() < 2)
            {
                continue;
            }
            // Open outfile
            std::string filename = gmx::formatString("%s_%s.xvg",
                                                     interactionTypeToString(itype).c_str(),
                                                     id.id().c_str());
            std::vector<gmx::RVec> forces;
            gmx::RVec rvnul = { 0, 0, 0 };
            FILE *fp = gmx_ffopen(filename.c_str(), "w");
            std::map<InteractionType, double> energies;
            // Not every parameter makes a topology, e.g. when bonds are
            // missing for the atoms of an angle, these are skipped.
            bool built = false;
            auto buildTopology = [&](std::vector<gmx::RVec> *coordinates)
            {
                MsgHandler buildHandler;
                buildHandler.setPrintLevel(ACTStatus::Fatal);
                top.build(&buildHandler, pd, coordinates, 175.0, 5.0, missingParameters::Error);
                built = buildHandler.ok() && top.hasInteractions(itype);
                forces.resize(top.nAtoms(), rvnul);
                return built;
            };
            switch (itype)
            {
            case InteractionType::BONDS:
            case InteractionType::VDW:
            case InteractionType::ELECTROSTATICS:
                {
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 } };
                    // top.build will resize coordinates array.
                    if (!buildTopology(&coordinates))
                    {
                        break;
                    }
                    // First atom is zero, second must be the other particle
                    size_t jatom = 1+top.atoms()[0].shells().size()+top.atoms()[0].vsites().size();
                    if (jatom == 0)
                    {
                        GMX_THROW(gmx::InternalError(gmx::formatString("Could not find a second atom to make a plot, there are %zu atoms, interactionType %s, id %s", top.nAtoms(), interactionTypeToString(itype).c_str(), id.id().c_str()).c_str()));
                    }
                    {
                        // Now remove everything except what we want to compute.
                        for (const auto &ik : i2s)
                        {
                            if (itype != ik.first)
                            {
                                top.removeInteractionType(ik.first);
                            }
                        }
                    }
                    std::vector<double> rr, vv, ff;
                    // Now do the calculations and store the energy
                    double r0 = 0.05, r1 = 1.0, delta = 0.001;
                    int    nsteps = (r1-r0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double x = r0+i*delta;
                        coordinates[jatom][0] = x;
                        // Move shell position
                        for(int j : top.atoms()[jatom].shells())
                        {
                            coordinates[j][0] = x;
                        }
                        // Move vsite position
                        for(int j : top.atoms()[jatom].vsites())
                        {
                            coordinates[j][0] = x;
                        }
                        rr.push_back(x);
                        energies.clear();
                        for(size_t k = 0; k < top.nAtoms(); k++)
                        {
                            copy_rvec(rvnul, forces[k]);
                        }
                        // We need to compute the sum of all relevant forces and run shell minimization.
                        compute(msghandler, pd, &top, &coordinates, &forces, &energies);

                        auto ener = energies[InteractionType::EPOT];

                        if (ener == 0 && InteractionType::VDW == itype)
                        {
                            auto irep = InteractionType::EXCHANGE;
                            auto idsp = InteractionType::DISPERSION;
                            if (energies.find(irep) != energies.end() &&
                                energies.find(idsp) != energies.end())
                            {
                                ener = energies[irep] + energies[idsp];
                            }
                            auto iec = InteractionType::VDWCORRECTION;
                            if (energies.find(iec) != energies.end())
                            {
                                ener += energies[iec];
                            }
                        }
                        fprintf(fp, "%10g  %10g\n", x, ener);
                        vv.push_back(ener);
                        ff.push_back(forces[0][0]);
                    }
                    // Check whether the force is the derivative of the energy.
                    // The five point difference is accurate enough also where
                    // the potential is steep, and the absolute tolerance
                    // keeps the check meaningful where the force is zero.
                    for(size_t i = 2; i+2 < vv.size(); i++)
                    {
                        double fnumeric = (8*(vv[i+1]-vv[i-1]) - (vv[i+2]-vv[i-2]))/(12*delta);
                        if (std::abs(fnumeric-ff[i]) > 1e-2*std::abs(ff[i]) + 1e-6)
                        {
                            msghandler->msg(ACTStatus::Warning,
                                            gmx::formatString("%s: Force %g at r = %g, the derivative of the energy is %g",
                                                              filename.c_str(), ff[i], rr[i], fnumeric));
                        }
                    }
                }
                break;
            case InteractionType::ANGLES:
                {
                    std::vector<gmx::RVec> acoords     = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 } };
                    std::vector<gmx::RVec> coordinates = acoords;
                    // top.build will resize coordinates array.
                    if (!buildTopology(&coordinates))
                    {
                        break;
                    }
                    size_t ak = 0;
                    for(const auto &entry : top.entry(itype))
                    {
                        const auto &ind = entry->atomIndices();
                        GMX_RELEASE_ASSERT(ind.size() == acoords.size(), "Number of atoms in angle");
                        for(size_t m = 0; m < ind.size(); m++)
                        {
                            coordinates[ind[m]] = acoords[m];
                            ak = ind[m];
                        }
                        break;
                    }
                    double th0 = 45, th1 = 145, delta = 1;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[ak][0] = 1+std::cos(theta*DEG2RAD);
                        coordinates[ak][1] = std::sin(theta*DEG2RAD);
                        double ener = potentialEnergy(msghandler, potential, top.entry(itype),
                                                      top.atoms(), coordinates, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", theta, ener);
                    }
                }
                break;
            case InteractionType::LINEAR_ANGLES:
                {
                    //! \todo take a into account
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 }, { 2, 0, 0 } };
                    if (!buildTopology(&coordinates))
                    {
                        break;
                    }
                    double r0 = 0.0, r1 = 0.1, delta = 0.001;
                    int    nsteps = (r1-r0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double xx = (r0+i*delta);
                        coordinates[1][1] = xx;
                        double ener = potentialEnergy(msghandler, potential, top.entry(itype),
                                                      top.atoms(), coordinates, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", xx, ener);
                    }

                }
                break;
            case InteractionType::PROPER_DIHEDRALS:
                {
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 1 } };
                    if (!buildTopology(&coordinates))
                    {
                        break;
                    }
                    double th0 = 0, th1 = 360, delta = 2;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[3][0] = 1+std::cos(theta*DEG2RAD);
                        coordinates[3][1] = 1+std::sin(theta*DEG2RAD);
                        double ener = potentialEnergy(msghandler, potential, top.entry(itype),
                                                      top.atoms(), coordinates, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", theta, ener);
                    }
                }
                break;
            case InteractionType::IMPROPER_DIHEDRALS:
                {
                    std::vector<gmx::RVec> coordinates = { { 1, 0.5, 0 }, { 0, 0, 0 }, { 2, 0, 0 }, { 1, 1.5, 0 } };
                    if (!buildTopology(&coordinates))
                    {
                        break;
                    }
                    double th0 = -0.02, th1 = 0.02, delta = 0.001;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[3][2] = theta;
                        double ener = potentialEnergy(msghandler, potential, top.entry(itype),
                                                      top.atoms(), coordinates, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", theta, ener);
                    }
                }
                break;
            default: // will not plot this interaction
                break;
            }
            gmx_ffclose(fp);
            if (!built)
            {
                std::remove(filename.c_str());
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("Cannot plot %s %s, there is no topology for it",
                                                  interactionTypeToString(itype).c_str(),
                                                  id.id().c_str()));
            }
        }
    }
}

} // namespace alexandria
