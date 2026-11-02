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

#include "sensitivity.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "act/utility/regression.h"
#include "act/utility/stringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace alexandria
{

//! \brief Map from enum to string
std::map<CalcDev, const char *> cdMap =
    {
        { CalcDev::Compute, "Compute" },
        { CalcDev::ComputeAll, "ComputeAll" },
        { CalcDev::Parameters, "Parameter" },
        { CalcDev::Stop, "Stop" }
    };

//! \return a string corresponding to a CalcDev.
const char *calcDevName(CalcDev cd)
{
    return cdMap[cd];
}

//! \brief Map SensitivityMode to string
std::map<SensitivityMode, std::string> sens2str =
    {
        { SensitivityMode::No,   "No" },
        { SensitivityMode::OneD, "1D" },
        { SensitivityMode::TwoD, "2D" }
    };

SensitivityMode stringToSensitivityMode(const std::string &str)
{
    for(auto &s2s : sens2str)
    {
        if (str == s2s.second)
        {
            return s2s.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such sensitivity analysis %s",
                                                       str.c_str()).c_str()));
    // To satisfy the compiler
    return SensitivityMode::OneD;
}

const std::string &sensitivityModeToString(SensitivityMode mode)
{
    return sens2str[mode];
}

void Sensitivity::computeForceConstants(gmx::TextWriter *tw)
{
    if (p_.size() < 3)
    {
        if (tw)
        {
            tw->writeStringFormatted("Not enough parameters %zu to do sensitivty analysis\n",
                                     p_.size());
        }
        return;
    }
    // The centre of the reduced coordinates is the original parameter
    // value, or the mean of the values if there is none.
    if (original_ < p_.size())
    {
        p0_ = p_[original_];
    }
    else
    {
        p0_ = 0;
        for (auto p : p_)
        {
            p0_ += p;
        }
        p0_ /= p_.size();
    }
    // The scale is the largest distance to the centre, such that u is in [-1, 1]
    scale_ = 0;
    for (auto p : p_)
    {
        scale_ = std::max(scale_, std::abs(p - p0_));
    }
    if (scale_ == 0)
    {
        if (tw)
        {
            tw->writeStringFormatted("All %zu parameter values are the same, cannot do sensitivity analysis\n",
                                     p_.size());
        }
        scale_ = 1;
        return;
    }
    MatrixWrapper M(3, p_.size());
    std::vector<double> solution(3, 0.0);
    for (size_t i = 0; i < p_.size(); ++i)
    {
        double u = (p_[i] - p0_)/scale_;
        M.set(0, i, u*u);
        M.set(1, i, u);
        M.set(2, i, 1.0);
    }
    auto result = M.solve(chi2_, &solution);
    if (result == 0)
    {
        fitted_ = true;
        A_ = solution[0];
        B_ = solution[1];
        C_ = solution[2];
        // A curvature below the rounding errors of the chi2 values is
        // zero, such that a straight line has no minimum.
        double chi2max = 0;
        for (auto chi2 : chi2_)
        {
            chi2max = std::max(chi2max, std::abs(chi2));
        }
        if (std::abs(A_) <= 16*std::numeric_limits<double>::epsilon()*chi2max)
        {
            A_ = 0;
        }
        // chi2 = A ((p - p0)/s)^2 + B (p - p0)/s + C = a p^2 + b p + c
        a_ = A_/(scale_*scale_);
        b_ = B_/scale_ - 2*a_*p0_;
        c_ = a_*p0_*p0_ - B_*p0_/scale_ + C_;
    }
}

double Sensitivity::pMin() const
{
    if (A_ == 0)
    {
        return 0;
    }
    return p0_ - B_*scale_/(2*A_);
}

double Sensitivity::chi2Min() const
{
    if (A_ == 0)
    {
        return 0;
    }
    return C_ - B_*B_/(4*A_);
}

bool Sensitivity::checkConvergence(double lower, double upper, double tolerance)
{
    checked_   = true;
    converged_ = false;
    lower_     = lower;
    upper_     = upper;
    if (A_ > 0 && original_ < p_.size())
    {
        converged_ = std::abs(p_[original_] - pMinBounded()) <= tolerance;
    }
    return converged_;
}

double Sensitivity::pMinBounded() const
{
    return std::clamp(pMin(), lower_, upper_);
}

nlohmann::ordered_json Sensitivity::toJson(size_t             index,
                                           const std::string &label) const
{
    nlohmann::ordered_json js;
    js["index"] = index;
    js["name"]  = label;
    // The label contains the particle type and the parameter name
    auto pn = split(label, ' ');
    js["particle"]  = pn.size() > 0 ? pn[0] : "";
    js["parameter"] = pn.size() > 1 ? pn[1] : "";
    if (checked_)
    {
        js["lower"] = lower_;
        js["upper"] = upper_;
    }
    if (original_ < p_.size())
    {
        js["p_orig"]    = p_[original_];
        js["chi2_orig"] = chi2_[original_];
    }
    auto points = nlohmann::ordered_json::array();
    for (size_t i = 0; i < p_.size(); ++i)
    {
        points.push_back({ { "p", p_[i] }, { "chi2", chi2_[i] } });
    }
    js["points"] = points;
    js["fit_ok"] = fitted_;
    if (fitted_)
    {
        js["a"] = a_;
        js["b"] = b_;
        js["c"] = c_;
    }
    // The extremum is undefined for a straight line
    bool extremum = fitted_ && A_ != 0;
    js["p_min"]    = extremum ? nlohmann::ordered_json(pMin()) : nlohmann::ordered_json();
    js["chi2_min"] = extremum ? nlohmann::ordered_json(chi2Min()) : nlohmann::ordered_json();
    if (checked_)
    {
        js["p_min_bounded"] = extremum ? nlohmann::ordered_json(pMinBounded()) : nlohmann::ordered_json();
        js["converged"]     = converged_;
    }
    return js;
}


}
