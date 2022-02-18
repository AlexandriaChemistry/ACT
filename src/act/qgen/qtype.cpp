/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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

#include "act/molprop/multipole_names.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

#include "qgen_resp.h"

static const double A2CM = E_CHARGE*1.0e-10;        /* e Angstrom to Coulomb meter */

static const double CM2D = SPEED_OF_LIGHT*1.0e+24;  /* Coulomb meter to Debye */

static inline double e2d(double a) {return a*ENM2DEBYE; }

// static inline int delta(int a, int b) { return ( a == b ) ? 1 : 0; }

namespace alexandria
{

static std::map<qType, std::string> qTypeNames = {
    { qType::ESP,       "ESP"        },
    { qType::Mulliken,  "Mulliken"   },
    { qType::Hirshfeld, "Hirshfeld"  },
    { qType::CM5,       "CM5"        },
    { qType::Calc,      "Alexandria" },
    { qType::Gasteiger, "Gasteiger"  },
    { qType::Elec,      "Electronic" }
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

QtypeProps::QtypeProps(qType qtype) : qtype_(qtype)
{
    for(auto &mpo : mpoMultiPoles)
    {
        std::vector<double> dummy;
        dummy.resize(multipoleNames(mpo).size(), 0.0);
        multipoles_.insert({mpo, std::move(dummy)});
    }
    resetMoments();
}

void QtypeProps::resetMoments()
{
    for(auto &m : multipoles_)
    {
        std::fill(m.second.begin(), m.second.end(), 0.0);
    }
}

void QtypeProps::setX(const gmx::HostVector<gmx::RVec> &x)
{
    // Check that arrays are equally long
    GMX_RELEASE_ASSERT(q_.size() - x.size() == 0,
                       gmx::formatString("Charge array (%d) should be the same size as coordinates (%d) for %s",
                                         static_cast<int>(q_.size()), static_cast<int>(x.size()), 
                                         qTypeName(qtype_).c_str()).c_str());
    x_.resizeWithPadding(x.size());
    for(int i = 0; i < x.size(); i++)
    {
        copy_rvec(x[i], x_[i]);
    }
    if (QgenResp_)
    {
        QgenResp_->updateAtomCoords(x);
    }
}

void QtypeProps::setQ(const std::vector<double> &q)
{
    q_.clear();
    for(auto &qq : q)
    {
        q_.push_back(qq);
    }
    // If QgenResp_ has not been allocated we likely don't need it.
    if (QgenResp_)
    {
        QgenResp_->updateAtomCharges(q);
    }
}

void QtypeProps::setQ(const t_atoms *atoms)
{
    q_.clear();
    for(int i = 0; i < atoms->nr; i++)
    {
        q_.push_back(atoms->atom[i].q);
    }
    // If QgenResp_ has not been allocated we likely don't need it.
    if (QgenResp_)
    {
        QgenResp_->updateAtomCharges(atoms);
    }
}

void QtypeProps::setQandX(const std::vector<double>        &q,
                          const gmx::HostVector<gmx::RVec> &x)
{
    // Check that arrays are equally long
    GMX_RELEASE_ASSERT(q.size() - x.size() == 0, 
                       gmx::formatString("Charge array should the same size as coordinates for %s",
                                         qTypeName(qtype_).c_str()).c_str());
    setQ(q);
    setX(x);
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

const std::vector<double> &QtypeProps::getMultipole(MolPropObservable mpo) const
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
        multipoles_[mpo] = mult;
    }
}

void QtypeProps::calcMoments()
{
    GMX_RELEASE_ASSERT(q_.size() > 0, gmx::formatString("No charges for %s", qTypeName(qtype_).c_str()).c_str());
    GMX_RELEASE_ASSERT(x_.size() > 0, gmx::formatString("No coordinates for %s", qTypeName(qtype_).c_str()).c_str());
    // distance of atoms to center of charge
    rvec   r;
    resetMoments(); 
    auto &dipole       = multipoles_[MolPropObservable::DIPOLE];
    auto &quadrupole   = multipoles_[MolPropObservable::QUADRUPOLE];
    auto &octupole     = multipoles_[MolPropObservable::OCTUPOLE];
    auto &hexadecapole = multipoles_[MolPropObservable::HEXADECAPOLE];
    for (size_t i = 0; i < q_.size(); i++)
    {
        rvec_sub(x_[i], coc_, r);
        int qindex = 0;
        int oindex = 0;
        int hindex = 0;
        for (int m = 0; m < DIM; m++)
        {
            dipole[m] += e2d(r[m]*q_[i]);
            for (int n = m; n < DIM; n++)
            {
                quadrupole[qindex++] += q_[i]*(r[m]*r[n])*NM2A*A2CM*CM2D*10;
                for (int o = n; o < DIM; o++)
                {    
                    octupole[oindex++] += q_[i]*(r[m]*r[n]*r[o])*NM2A*1e-12*A2CM*CM2D*1e10*10;
                    for (int p = o; p < DIM; p++)
                    {  
                        hexadecapole[hindex++] += q_[i]*(r[m]*r[n]*r[o]*r[p])*NM2A*1e-12*1e-12*A2CM*CM2D*1e10*1e10*10;
                    }    
                }
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "Quadrupole: %7.3f %7.3f %7.3f\n", 
                quadrupole[multipoleIndex({XX, XX})],
                quadrupole[multipoleIndex({YY, YY})],
                quadrupole[multipoleIndex({ZZ, ZZ})]);
    }
    // Compute trace divided by 3
    double tr = 0.0;
    for (int m = 0; m < DIM; m++)
    {
        tr += quadrupole[multipoleIndex({m, m})]/3.0;
    }
    // Subtract trace/3 to make the quadrupole traceless
    for (int m = 0; m < DIM; m++)
    {
        quadrupole[multipoleIndex({m, m})] -= tr;
    }
    if (debug)
    {
        fprintf(debug, "Traceless:  %7.3f %7.3f %7.3f\n", 
                quadrupole[multipoleIndex({XX, XX})],
                quadrupole[multipoleIndex({YY, YY})],
                quadrupole[multipoleIndex({ZZ, ZZ})]);
    }
}

QgenResp *QtypeProps::qgenResp()
{
    if (QgenResp_ == nullptr)
    {
        QgenResp_ = new QgenResp;
    }
    return QgenResp_;
}
    
const QgenResp *QtypeProps::qgenRespConst()
{
    if (QgenResp_ == nullptr)
    {
        QgenResp_ = new QgenResp;
    }
    return QgenResp_;
}

void QtypeProps::copyRespQ()
{
    GMX_RELEASE_ASSERT(QgenResp_ != nullptr, "QgenResp_ has not been initialized");
    int natom = q_.size();
    q_.clear();
    for(int i = 0; i < natom; i++)
    {
        q_.push_back(QgenResp_->getCharge(i));
    }
}

} // namespace
