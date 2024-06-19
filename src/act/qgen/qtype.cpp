/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2023
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

#include "qtype.h"

#include <map>
#include <random>

#include "act/alexandria/topology.h"
#include "act/forces/forcecomputer.h"
#include "act/molprop/multipole_names.h"
#include "act/utility/units.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

#include "qgen_resp.h"

namespace alexandria
{

static std::map<qType, std::string> qTypeNames = {
    { qType::ESP,       "qESP"        },
    { qType::Mulliken,  "qMulliken"   },
    { qType::Hirshfeld, "qHirshfeld"  },
    { qType::CM5,       "qCM5"        },
    { qType::Calc,      "Alexandria" },
    { qType::Gasteiger, "Gasteiger"  },
    { qType::Elec,      "Electronic" },
    { qType::ACM,       "qACM"        }
};

const std::string &qTypeName(qType qt)
{
    return qTypeNames[qt];
}

const std::map<qType, std::string> &qTypes()
{
    return qTypeNames;
}

qType stringToQtype(const std::string &type)
{
    for (const auto &qn : qTypeNames)
    {
        if (qn.second.compare(type) == 0)
        {
            return qn.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown charge type %s", type.c_str()).c_str()));
    // To make the compiler happy
    return qType::ESP;
}

QtypeProps::QtypeProps(qType                         qtype,
                       const std::vector<ActAtom>   &atoms,
                       const std::vector<gmx::RVec> &coords) : qtype_(qtype)
{
    // Copy coordinates
    x_ = coords;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        atomNumber_.push_back(atoms[i].atomicNumber());
        q_.push_back(atoms[i].charge());
    }
    // Compute CoC
    computeCoC();
}

void QtypeProps::computeCoC()
{
    if (atomNumber_.size() != x_.size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Mismatch in array sizes, atomNumber %zu coords %zu elements.",
                                                       atomNumber_.size(), x_.size()).c_str()));
    }
    int atntot = 0;
    clear_rvec(coc_);
    for (size_t i = 0; i < atomNumber_.size(); i++)
    {
        auto atn = atomNumber_[i];
        atntot  += atn;
        for (auto m = 0; m < DIM; m++)
        {
            coc_[m] += x_[i][m]*atn;
        }
    }
    /* Center of charge */
    if (atntot > 0)
    {
        svmul((1.0/atntot), coc_, coc_);
    }
}

void QtypeProps::initializeMoments()
{
    for(auto &mpo : mpoMultiPoles)
    {
        std::vector<double> dummy;
        dummy.resize(multipoleNames(mpo).size(), 0.0);
        multipoles_.insert({mpo, std::move(dummy)});
    }
}

void QtypeProps::resetMoments()
{
    for(auto &m : multipoles_)
    {
        std::fill(m.second.begin(), m.second.end(), 0.0);
    }
}

void QtypeProps::setX(const std::vector<gmx::RVec> &x)
{
    // Check that arrays are equally long
    GMX_RELEASE_ASSERT(q_.size() == x.size(),
                       gmx::formatString("Charge array (%d) should be the same size as coordinates (%d) for %s",
                                         static_cast<int>(q_.size()), static_cast<int>(x.size()), 
                                         qTypeName(qtype_).c_str()).c_str());
    x_ = x;
}

void QtypeProps::setQ(const std::vector<double> &q)
{
    q_.resize(q.size());
    std::copy(q.begin(), q.end(), q_.begin());
}

void QtypeProps::setQ(const std::vector<ActAtom> &atoms)
{
    q_.resize(atoms.size(), 0.0);
    atomNumber_.resize(atoms.size(), 0);
    for(size_t i = 0; i < atoms.size(); i++)
    {
        atomNumber_[i] = atoms[i].atomicNumber();
        q_[i] = atoms[i].charge();
    }
    computeCoC();
}

void QtypeProps::setQandX(const std::vector<double>    &q,
                          const std::vector<gmx::RVec> &x)
{
    q_.resize(q.size());
    setX(x);
    setQ(q);
    computeCoC();
}

void QtypeProps::setQandX(const std::vector<ActAtom>   &atoms,
                          const std::vector<gmx::RVec> &x)
{
    q_.resize(atoms.size());
    setX(x);
    setQ(atoms);
}

void QtypeProps::setPolarizabilityTensor(const tensor &alpha)
{
    copy_mat(alpha, alpha_);
    auto a = gmx::square(alpha_[XX][XX] - alpha_[YY][YY]);
    auto b = gmx::square(alpha_[XX][XX] - alpha_[ZZ][ZZ]);
    auto c = gmx::square(alpha_[ZZ][ZZ] - alpha_[YY][YY]);
    auto d = 6 * (gmx::square(alpha_[XX][YY]) + gmx::square(alpha_[XX][ZZ]) + gmx::square(alpha_[ZZ][YY]));

    anisotropy_ = sqrt(1/2.0) * sqrt(a + b + c + d);
    isotropy_   = (alpha_[XX][XX]+alpha_[YY][YY]+alpha_[ZZ][ZZ])/3;
    hasAlpha_   = true;
}

double QtypeProps::dipole() const
{
    auto mm = multipoles_.find(MolPropObservable::DIPOLE);
    if (multipoles_.end() == mm)
    {
        GMX_THROW(gmx::InternalError("No dipole!"));
    }
    double d2 = 0;
    for (int m = 0; m < DIM; m++)
    {
        d2 += gmx::square(mm->second[m]);
    }
    return std::sqrt(d2);
}

bool QtypeProps::hasMultipole(MolPropObservable mpo) const
{
    return (multipoles_.find(mpo) != multipoles_.end());
}

const std::vector<double> QtypeProps::getMultipole(MolPropObservable mpo) const
{
    auto mf = multipoles_.find(mpo);
    if (mf == multipoles_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No such observable %s in multipoles",
                                                       mpo_name(mpo)).c_str()));
    }
    return mf->second;
}

void QtypeProps::setMultipole(MolPropObservable mpo, const std::vector<double> &mult)
{
    if (multipoles_.find(mpo) == multipoles_.end())
    {
        multipoles_.insert({mpo, mult});
    }
    else
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Multipole %s already set before", mpo_name(mpo)).c_str()));
        multipoles_[mpo] = mult;
    }
}

void QtypeProps::calcPolarizability(const ForceField    *pd,
                                    const Topology      *top,
                                    const ForceComputer *forceComp)
{
    GMX_RELEASE_ASSERT(qType::Calc == qtype_, "Will only compute polarizability for Alexandria models");
    
    std::map<InteractionType, double> energies;
    gmx::RVec                         field  = { 0, 0, 0 };
    std::vector<gmx::RVec>            coords = x_;
    std::vector<gmx::RVec>            forces(x_.size());
    
    forceComp->compute(pd, top, &coords, &forces, &energies, field);
    setQ(top->atoms());
    setX(coords);
    calcMoments();
    auto mpo = MolPropObservable::DIPOLE;
    if (!hasMultipole(mpo))
    {
        GMX_THROW(gmx::InternalError("No dipole present to compute polarizablity."));
    }
    auto mu_ref = getMultipole(mpo);
    // Convert from e nm2/V to cubic nm
    double enm2_V = E_CHARGE*1e6*1e-18/(4*M_PI*EPSILON0_SI)*1e21;

    tensor alpha;
    clear_mat(alpha);
    // Units are not relevant since they drop out!
    // However, a field that is too large may lead to strange shell displacements.
    double efield = 0.1;
    for (int m = 0; m < DIM; m++)
    {
        field[m] = efield;
        forceComp->compute(pd, top, &coords, &forces, &energies, field);
        setX(coords);
        field[m] = 0;
        calcMoments();
        auto qmu = getMultipole(mpo);
        for (int n = 0; n < DIM; n++)
        {
            alpha[n][m] = enm2_V*((qmu[n]-mu_ref[n])/efield);
        }
    }
    setPolarizabilityTensor(alpha);
}

void QtypeProps::calcMoments()
{
    GMX_RELEASE_ASSERT(q_.size() > 0, gmx::formatString("No charges for %s", qTypeName(qtype_).c_str()).c_str());
    GMX_RELEASE_ASSERT(x_.size() > 0, gmx::formatString("No coordinates for %s", qTypeName(qtype_).c_str()).c_str());
    bool AllZero = true;
    for(size_t i = 0; AllZero && i < atomNumber_.size(); i++)
    {
        AllZero = AllZero && (atomNumber_[i] == 0);
    }
    if (AllZero)
    {
        fprintf(stderr, "All atomnumbers are zero when computing multipoles. Center of Charge not reliable.\n");
    }
    // distance of atoms to center of charge
    rvec   r;
    resetMoments();
    computeCoC();
    auto dip  = MolPropObservable::DIPOLE;
    auto quad = MolPropObservable::QUADRUPOLE;
    auto oct  = MolPropObservable::OCTUPOLE;
    auto hex  = MolPropObservable::HEXADECAPOLE;
    for (size_t i = 0; i < q_.size(); i++)
    {
        rvec_sub(x_[i], coc_, r);
        int qindex = 0;
        int oindex = 0;
        int hindex = 0;
        for (int m = 0; m < DIM; m++)
        {
            if (hasMultipole(dip))
            {
                multipoles_[dip][m] += r[m]*q_[i];
            }
            for (int n = m; n < DIM; n++)
            {
                if (hasMultipole(quad))
                {
                    multipoles_[quad][qindex++] += q_[i]*(r[m]*r[n]);
                }
                for (int o = n; o < DIM; o++)
                {
                    if (hasMultipole(oct))
                    {    
                        multipoles_[oct][oindex++] += q_[i]*(r[m]*r[n]*r[o]);
                    }
                    for (int p = o; p < DIM; p++)
                    {
                        if (hasMultipole(hex))
                        {
                            multipoles_[hex][hindex++] += q_[i]*(r[m]*r[n]*r[o]*r[p]);
                        }
                    }    
                }
            }
        }
    }
    if (hasMultipole(quad))
    {
        if (debug)
        {
            fprintf(debug, "Quadrupole: %7.3f %7.3f %7.3f\n", 
                    multipoles_[quad][multipoleIndex({XX, XX})],
                    multipoles_[quad][multipoleIndex({YY, YY})],
                    multipoles_[quad][multipoleIndex({ZZ, ZZ})]);
        }
        // Compute trace divided by 3
        double tr = 0.0;
        for (int m = 0; m < DIM; m++)
        {
            tr += multipoles_[quad][multipoleIndex({m, m})]/3.0;
        }
        // Subtract trace/3 to make the quadrupole traceless
        for (int m = 0; m < DIM; m++)
        {
            multipoles_[quad][multipoleIndex({m, m})] -= tr;
        }
        if (debug)
        {
            fprintf(debug, "Traceless:  %7.3f %7.3f %7.3f\n", 
                    multipoles_[quad][multipoleIndex({XX, XX})],
                    multipoles_[quad][multipoleIndex({YY, YY})],
                    multipoles_[quad][multipoleIndex({ZZ, ZZ})]);
        }
    }
}

} // namespace
