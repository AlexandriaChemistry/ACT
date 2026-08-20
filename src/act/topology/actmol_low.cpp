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

#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forces/forcecomputer.h"
#include "act/forces/forcecomputerutils.h"
#include "act/molprop/multipole_names.h"
#include "act/qgen/qgen_acm.h"
#include "act/topology/actmol.h"
#include "act/topology/topology.h"
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

bool is_linear(MsgHandler *msghandler,
               const rvec xi, const rvec xj,
               const rvec xk, real th_toler)
{
    rvec r_ij, r_kj;
    real costh, th;

    th = fabs(RAD2DEG*bond_angle(xi, xj, xk, r_ij, r_kj, &costh));
    if ((th > th_toler) || (th < 180-th_toler))
    {
        if (msghandler && msghandler->debug())
        {
            msghandler->writeDebug(gmx::formatString("Angle is %g, th_toler is %g\n", th, th_toler));
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

void put_in_box(int natom, matrix box, const rvec x[], real dbox)
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

void calc_rotmatrix(const rvec target_vec, const rvec ref_vec, matrix rotmatrix)
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

double computeAtomizationEnergy(MsgHandler                 *msghandler,
                                const std::vector<ActAtom> &atoms,
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
            if (msghandler && msghandler->debug())
            {
                msghandler->writeDebug(gmx::formatString("Found atomization energy terms for %2s at T = %g K. %s: %g (%s) %s: %g (%s)\n",
                                                         a.element().c_str(), temperature,
                                                         H0.c_str(), h0_hT, h0unit.c_str(),
                                                         DHf.c_str(), dhF, dhFunit.c_str()));
            }
            atomizationEnergy += convertToGromacs(h0_hT, h0unit) + convertToGromacs(dhF, dhFunit);
        }
    }
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("Total atomization energy %g\n", atomizationEnergy));
    }
    return atomizationEnergy;
}

void analyse_multipoles(MsgHandler                                *msg_handler,
                        const ACTMol                              *mol,
                        const std::map<MolPropObservable, double> &toler,
                        const ForceField                          *pd,
                        const ForceComputer                                         *forceComputer,
                        std::map<MolPropObservable, std::map<iMolSelect, qtStats> > *lsq)
{
    gmx::TextWriter *tw = nullptr;
    if (msg_handler)
    {
        tw = msg_handler->tw();
    }
    auto topology = mol->topology();
    bool doForce  = pd->polarizable() || topology->hasVsites();
    auto qprops = mol->qPropsConst();
    for(auto qp = qprops.begin(); qp < qprops.end(); ++qp)
    {
        auto qelec = qp->qPqmConst();
        auto qcalc = qp->qPact();
        qcalc->initializeMoments();
        if (doForce)
        {
            std::vector<gmx::RVec>            forces(topology->nAtoms());
            std::map<InteractionType, double> energies;
            auto                              myx = qcalc->x();
            forceComputer->compute(msg_handler, pd, topology, &myx, &forces, &energies);
            qcalc->setX(myx);
        }
        qcalc->calcMoments(msg_handler);

        for(auto &mpo : mpoMultiPoles)
        {
            const char *name   = mpo_name(mpo);
            const char *unit   = mpo_unit2(mpo);
            double      factor = convertFromGromacs(1, unit);
            std::vector<double> Telec;
            if (qelec.hasMultipole(mpo))
            {
                tw->writeStringFormatted("Electronic %s (%s):\n", name, unit);
                Telec = qelec.getMultipole(msg_handler, mpo);
                for(const auto &pm : formatMultipole(mpo, Telec))
                {
                    tw->writeLine(pm);
                }
            }
            real delta = 0;
            auto Tcalc = qcalc->getMultipole(msg_handler, mpo);
            tw->writeStringFormatted("Calc %s (%s):\n", name, unit);
            for(const auto &pm : formatMultipole(mpo, Tcalc))
            {
                tw->writeLine(pm);
            }

            std::vector<double> diff;
            auto qt = qcalc->qtype();
            if (Telec.size() == Tcalc.size())
            {
                for(size_t i = 0; i < Tcalc.size(); i++)
                {
                    auto tc = Tcalc[i];
                    auto te = Telec[i];
                    diff.push_back(tc-te);
                    delta += gmx::square(tc-te);
                    if (lsq)
                    {
                        (*lsq)[mpo][mol->datasetType()][qt].add_point(factor*te, factor*tc, 0, 0);
                    }
                }
                double rms = std::sqrt(delta/Tcalc.size());
                std::string flag("");
                if (!toler.empty() && rms > toler.at(mpo))
                {
                    flag = " MULTI";
                }
                tw->writeStringFormatted("%s-Electronic Norm %g RMS = %g (%s)%s:\n",
                                         qPropertyTypeName(qt).c_str(),
                                         factor*std::sqrt(delta), factor*rms, unit, flag.c_str());
                for(const auto &pm : formatMultipole(mpo, diff))
                {
                    tw->writeLine(pm);
                }
            }
        }
    }
}
} // namespace alexandria
