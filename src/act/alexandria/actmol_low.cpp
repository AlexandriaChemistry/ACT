/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actmol_low.h"

#include <cstdio>
#include <cstring>
#include <cmath>

#include "act/alexandria/topology.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forces/combinationrules.h"
#include "act/forces/forcecomputerutils.h"
#include "act/qgen/qgen_acm.h"
#include "act/utility/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/strconvert.h"

namespace alexandria
{

bool is_planar(const rvec xi, const rvec xj, const rvec xk,
               const rvec xl, real phi_toler)
{
    rvec r_ij, r_kj, r_kl, m, n;
    real phi;

    phi = RAD2DEG*dih_angle(xi, xj, xk, xl, r_ij, r_kj, r_kl, m, n);

    return (fabs(phi) < phi_toler);
}

bool is_linear(const rvec xi, const rvec xj,
               const rvec xk, real th_toler)
{
    rvec r_ij, r_kj;
    real costh, th;

    th = fabs(RAD2DEG*bond_angle(xi, xj, xk, r_ij, r_kj, &costh));
    if ((th > th_toler) || (th < 180-th_toler))
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Angle is %g, th_toler is %g\n", th, th_toler);
        }
        return true;
    }
    return false;
}

real calc_relposition(const ForceField               *pd,
                      const std::vector<std::string> &atoms,
                      const std::vector<double>      &bondOrders)
{
    if (atoms.size() != 3)
    {
        GMX_THROW(gmx::InternalError("calc_relposition called with incorrect number of atoms"));
    }
    if (bondOrders.size() != 2)
    {
        GMX_THROW(gmx::InternalError("calc_relposition called with incorrect number of bondOrders"));
    }
    double b0 = 0, b1 = 0, relative_position = 0;

    Identifier aij({atoms[0], atoms[1]}, { bondOrders[0] }, CanSwap::Yes);
    Identifier ajk({atoms[1], atoms[2]}, { bondOrders[1] }, CanSwap::Yes);

    std::string type("bondlength");
    auto &fs = pd->findForcesConst(InteractionType::BONDS);
    auto bij = fs.findParameterTypeConst(aij, type);
    auto bjk = fs.findParameterTypeConst(ajk, type);

    b0 = convertToGromacs(bij.value(), bij.unit());
    b1 = convertToGromacs(bjk.value(), bjk.unit());

    relative_position = (b1/(b0+b1));

    return relative_position;
}

std::vector<double> getDoubles(const std::string &s)
{
    std::vector<double> d;

    for (auto &ss : gmx::splitString(s))
    {
        d.push_back(gmx::doubleFromString(ss.c_str()));
    }
    return d;
}

void put_in_box(int natom, matrix box, rvec x[], real dbox)
{
    int  i, m;
    rvec xmin, xmax, xcom;

    clear_rvec(xcom);
    copy_rvec(x[0], xmin);
    copy_rvec(x[0], xmax);
    for (i = 0; (i < natom); i++)
    {
        rvec_inc(xcom, x[i]);
        for (m = 0; (m < DIM); m++)
        {
            if (xmin[m] > x[i][m])
            {
                xmin[m] = x[i][m];
            }
            else if (xmax[m] < x[i][m])
            {
                xmax[m] = x[i][m];
            }
        }
    }
    for (m = 0; (m < DIM); m++)
    {
        xcom[m]  /= natom;
        box[m][m] = (dbox+xmax[m]-xmin[m]);
    }
}

void calc_rotmatrix(rvec target_vec, rvec ref_vec, matrix rotmatrix)
{
    rvec au = {0, 0, 0};
    rvec bu = {0, 0, 0};

    svmul((1.0/norm(target_vec)), target_vec, au);
    svmul((1.0/norm(ref_vec)), ref_vec, bu);

    rotmatrix[0][0] = bu[0]*au[0];
    rotmatrix[0][1] = bu[0]*au[1];
    rotmatrix[0][2] = bu[0]*au[2];
    rotmatrix[1][0] = bu[1]*au[0];
    rotmatrix[1][1] = bu[1]*au[1];
    rotmatrix[1][2] = bu[1]*au[2];
    rotmatrix[2][0] = bu[2]*au[0];
    rotmatrix[2][1] = bu[2]*au[1];
    rotmatrix[2][2] = bu[2]*au[2];
}

double computeAtomizationEnergy(const std::vector<ActAtom> &atoms,
                                const AtomizationEnergy    &atomenergy,
                                double                      temperature)
{
    double atomizationEnergy = 0;
    std::string H0("H(0)-H(T)");
    std::string DHf("DHf(T)");
    for (auto &a : atoms)
    {
        if (a.pType() == ActParticle::Atom)
        {
            std::string h0unit;
            double h0_hT = atomenergy.term(a.element(), 0, "exp",
                                           H0, temperature, &h0unit, nullptr);
            std::string dhFunit;
            double dhF   = atomenergy.term(a.element(), 0, "exp",
                                           DHf, temperature, &dhFunit, nullptr);
            if (debug)
            {
                fprintf(debug, "Found atomization energy terms for %2s at T = %g K. %s: %g (%s) %s: %g (%s)\n",
                        a.element().c_str(), temperature,
                        H0.c_str(), h0_hT, h0unit.c_str(),
                        DHf.c_str(), dhF, dhFunit.c_str());
            }
            atomizationEnergy += convertToGromacs(h0_hT, h0unit) + convertToGromacs(dhF, dhFunit);
        }
    }
    if (debug)
    {
        fprintf(debug, "Total atomization energy %g\n", atomizationEnergy);
    }
    return atomizationEnergy;
}

} // namespace alexandria
