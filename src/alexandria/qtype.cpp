/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021
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

#include "gromacs/math/units.h"

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

void QtypeProps::setMu(const rvec mu)
{
    copy_rvec(mu, mu_);
    dipole_ = norm(mu_);
}

void QtypeProps::setQuadrupole(const tensor quad)
{
    copy_mat(quad, quadrupole_);
}

void QtypeProps::calcMoments()
{
    GMX_RELEASE_ASSERT(q_.size() > 0, gmx::formatString("No charges for %s", qTypeName(qtype_).c_str()).c_str());
    GMX_RELEASE_ASSERT(x_.size() > 0, gmx::formatString("No coordinates for %s", qTypeName(qtype_).c_str()).c_str());
    // distance of atoms to center of charge
    rvec   r; 
    real   r2;
    clear_mat(quadrupole_);
    clear_rvec(mu_);
    for (size_t i = 0; i < q_.size(); i++)
    {
        rvec_sub(x_[i], coc_, r);
        for (int m = 0; m < DIM; m++)
        {
            mu_[m] += e2d(r[m]*q_[i]);
        }
        r2   = iprod(r, r);
        for (int m = 0; m < DIM; m++)
        {
            for (int n = m; n < DIM; n++)
            {
                quadrupole_[m][n] += q_[i]*(r[m]*r[n])*NM2A*A2CM*CM2D*10;
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "Quadrupole: %7.3f %7.3f %7.3f\n", quadrupole_[XX][XX],
                quadrupole_[YY][YY], quadrupole_[ZZ][ZZ]);
    }
    // Compute trace divided by 3
    double tr = trace(quadrupole_)/3.0;
    // Subtract trace/3 to make the quadrupole traceless
    for (int m = 0; m < DIM; m++)
    {
        quadrupole_[m][m] -= tr;
    }
    if (debug)
    {
        fprintf(debug, "Traceless:  %7.3f %7.3f %7.3f\n", quadrupole_[XX][XX],
                quadrupole_[YY][YY], quadrupole_[ZZ][ZZ]);
    }
    dipole_ = norm(mu_);
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
