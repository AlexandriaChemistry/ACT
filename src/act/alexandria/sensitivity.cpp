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

#include <Eigen/Eigenvalues>
#include <nlohmann/json.hpp>

#include "act/utility/regression.h"
#include "act/utility/stringutil.h"
#include "gromacs/math/functions.h"
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
        // Without any change of the chi2 there is no slope either
        if (flat())
        {
            A_ = 0;
            B_ = 0;
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

bool Sensitivity::flat() const
{
    if (chi2_.size() < 2)
    {
        return false;
    }
    const auto minmax  = std::minmax_element(chi2_.begin(), chi2_.end());
    double chi2max = std::max(std::abs(*minmax.first), std::abs(*minmax.second));
    return *minmax.second - *minmax.first <= 16*std::numeric_limits<double>::epsilon()*chi2max;
}

double Sensitivity::value(double p) const
{
    if (!fitted_)
    {
        return 0;
    }
    double u = (p - p0_)/scale_;
    return (A_*u + B_)*u + C_;
}

double Sensitivity::derivative(double p) const
{
    if (!fitted_)
    {
        return 0;
    }
    double u = (p - p0_)/scale_;
    return (2*A_*u + B_)/scale_;
}

nlohmann::ordered_json Sensitivity::toJson(size_t             index,
                                           const std::string &label) const
{
    nlohmann::ordered_json js;
    js["index"] = index;
    js["name"]  = label;
    // The label contains the particle type and the parameter name
    auto pn = split(label, ' ');
    js["particle"]  = !pn.empty() ? pn[0] : "";
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
    js["flat"]   = flat();
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


//! \brief Map from PairShape to string
static const std::map<PairShape, std::string> shapeNames =
    {
        { PairShape::Minimum,    "minimum"    },
        { PairShape::Saddle,     "saddle"     },
        { PairShape::Maximum,    "maximum"    },
        { PairShape::Degenerate, "degenerate" }
    };

const std::string &pairShapeName(PairShape shape)
{
    return shapeNames.at(shape);
}

/*! \brief Relative precision of the determinant of the Hessian
 * A determinant smaller than this times its terms is zero.
 */
static const double determinantPrecision = 1e-8;

void Sensitivity2D::computeCrossTerm()
{
    computed_  = false;
    aij_       = 0;
    cornerRms_ = 0;
    if (!si_.fitted() || !sj_.fitted() || corners_.empty() ||
        si_.original() >= si_.pValues().size() || sj_.original() >= sj_.pValues().size())
    {
        return;
    }
    double pi     = si_.pValues()[si_.original()];
    double pj     = sj_.pValues()[sj_.original()];
    double chi2_0 = si_.chi2Values()[si_.original()];
    // Least squares fit of r = a_ij dx dy, where r is what the
    // one-dimensional parabolas do not describe
    double num = 0, den = 0;
    std::vector<double> r, dxdy;
    for (const auto &corner : corners_)
    {
        double xy = (corner[0] - pi)*(corner[1] - pj);
        r.push_back(corner[2] - (si_.value(corner[0]) + sj_.value(corner[1]) - chi2_0));
        dxdy.push_back(xy);
        num += r.back()*xy;
        den += xy*xy;
    }
    if (den == 0)
    {
        return;
    }
    aij_ = num/den;
    double rms = 0;
    for (size_t k = 0; k < r.size(); k++)
    {
        rms += gmx::square(r[k] - aij_*dxdy[k]);
    }
    cornerRms_ = std::sqrt(rms/r.size());
    computed_  = true;
}

std::optional<double> Sensitivity2D::rho() const
{
    if (!computed_ || si_.a() <= 0 || sj_.a() <= 0)
    {
        return std::nullopt;
    }
    return aij_/(2*std::sqrt(si_.a()*sj_.a()));
}

std::array<std::array<double, 2>, 2> Sensitivity2D::hessian() const
{
    return { { { 2*si_.a(), aij_ }, { aij_, 2*sj_.a() } } };
}

std::array<double, 2> Sensitivity2D::gradient() const
{
    if (si_.original() >= si_.pValues().size() || sj_.original() >= sj_.pValues().size())
    {
        return { 0, 0 };
    }
    return { si_.derivative(si_.pValues()[si_.original()]),
             sj_.derivative(sj_.pValues()[sj_.original()]) };
}

PairShape Sensitivity2D::shape() const
{
    if (!computed_)
    {
        return PairShape::Degenerate;
    }
    auto   H     = hessian();
    double det   = H[0][0]*H[1][1] - H[0][1]*H[1][0];
    double scale = std::abs(H[0][0]*H[1][1]) + H[0][1]*H[1][0];
    if (H[0][0] == 0 || H[1][1] == 0 || std::abs(det) <= determinantPrecision*scale)
    {
        return PairShape::Degenerate;
    }
    else if (det < 0)
    {
        return PairShape::Saddle;
    }
    else if (H[0][0] > 0)
    {
        return PairShape::Minimum;
    }
    return PairShape::Maximum;
}

std::optional<std::array<double, 2>> Sensitivity2D::stationaryPoint() const
{
    if (shape() == PairShape::Degenerate)
    {
        return std::nullopt;
    }
    auto   H   = hessian();
    auto   g   = gradient();
    double det = H[0][0]*H[1][1] - H[0][1]*H[1][0];
    // Solve H d = -g
    double di  = -(H[1][1]*g[0] - H[0][1]*g[1])/det;
    double dj  = -(H[0][0]*g[1] - H[1][0]*g[0])/det;
    return std::array<double, 2>{ si_.pValues()[si_.original()] + di,
                                  sj_.pValues()[sj_.original()] + dj };
}

double Sensitivity2D::value(double x, double y) const
{
    double pi     = si_.pValues()[si_.original()];
    double pj     = sj_.pValues()[sj_.original()];
    double chi2_0 = si_.chi2Values()[si_.original()];
    return si_.value(x) + sj_.value(y) - chi2_0 + aij_*(x - pi)*(y - pj);
}

std::optional<std::array<double, 2>> Sensitivity2D::boxMinimum(const std::array<double, 4> &bounds) const
{
    if (shape() != PairShape::Minimum)
    {
        return std::nullopt;
    }
    auto sp = stationaryPoint();
    if (sp.has_value() &&
        (*sp)[0] >= bounds[0] && (*sp)[0] <= bounds[1] &&
        (*sp)[1] >= bounds[2] && (*sp)[1] <= bounds[3])
    {
        return sp;
    }
    // The minimum of a convex surface outside the bounds is on an edge.
    // Along an edge the surface is a parabola in the other parameter.
    double pi = si_.pValues()[si_.original()];
    double pj = sj_.pValues()[sj_.original()];
    auto   H  = hessian();
    auto   g  = gradient();
    std::vector<std::array<double, 2>> candidates;
    for (double x : { bounds[0], bounds[1] })
    {
        double dy = -(g[1] + H[0][1]*(x - pi))/H[1][1];
        candidates.push_back({ x, std::clamp(pj + dy, bounds[2], bounds[3]) });
    }
    for (double y : { bounds[2], bounds[3] })
    {
        double dx = -(g[0] + H[0][1]*(y - pj))/H[0][0];
        candidates.push_back({ std::clamp(pi + dx, bounds[0], bounds[1]), y });
    }
    auto best = candidates[0];
    for (const auto &cand : candidates)
    {
        if (value(cand[0], cand[1]) < value(best[0], best[1]))
        {
            best = cand;
        }
    }
    return best;
}

bool Sensitivity2D::checkConvergence(const std::array<double, 4> &bounds,
                                     const std::array<double, 2> &tolerances)
{
    checked_   = true;
    bounds_    = bounds;
    converged_ = false;
    auto pmin  = boxMinimum(bounds);
    if (pmin.has_value())
    {
        converged_ = (std::abs(si_.pValues()[si_.original()] - (*pmin)[0]) <= tolerances[0] &&
                      std::abs(sj_.pValues()[sj_.original()] - (*pmin)[1]) <= tolerances[1]);
    }
    return converged_;
}

nlohmann::ordered_json Sensitivity2D::toJson(size_t             i,
                                             size_t             j,
                                             const std::string &namei,
                                             const std::string &namej) const
{
    using json = nlohmann::ordered_json;
    json js;
    js["i"]     = i;
    js["j"]     = j;
    js["names"] = { namei, namej };
    // The grid of parameter values and the chi-squared at each point,
    // chi2[k][l] is at p_i[k], p_j[l]
    const auto &pvi = si_.pValues();
    const auto &pvj = sj_.pValues();
    auto chi2 = json::array();
    for (size_t k = 0; k < pvi.size(); k++)
    {
        auto row = json::array();
        for (size_t l = 0; l < pvj.size(); l++)
        {
            json value;
            if (k == si_.original())
            {
                value = sj_.chi2Values()[l];
            }
            else if (l == sj_.original())
            {
                value = si_.chi2Values()[k];
            }
            else
            {
                for (const auto &corner : corners_)
                {
                    if (corner[0] == pvi[k] && corner[1] == pvj[l])
                    {
                        value = corner[2];
                    }
                }
            }
            row.push_back(value);
        }
        chi2.push_back(row);
    }
    js["grid"] = { { "p_i", pvi }, { "p_j", pvj }, { "chi2", chi2 } };
    js["a_ij"]       = computed_ ? json(aij_) : json();
    auto r           = rho();
    js["rho"]        = r.has_value() ? json(*r) : json();
    js["corner_rms"] = computed_ ? json(cornerRms_) : json();
    auto H           = hessian();
    js["hessian"]    = { { H[0][0], H[0][1] }, { H[1][0], H[1][1] } };
    auto g           = gradient();
    js["gradient"]   = { g[0], g[1] };
    auto sh          = shape();
    js["shape"]      = pairShapeName(sh);
    auto sp          = stationaryPoint();
    bool minimum     = sh == PairShape::Minimum && sp.has_value();
    js["p_min"]      = minimum ? json({ (*sp)[0], (*sp)[1] }) : json();
    js["chi2_min"]   = minimum ? json(value((*sp)[0], (*sp)[1])) : json();
    if (checked_)
    {
        auto bm             = boxMinimum(bounds_);
        js["p_min_bounded"] = bm.has_value() ? json({ (*bm)[0], (*bm)[1] }) : json();
        js["converged"]     = converged_;
    }
    return js;
}

nlohmann::ordered_json HessianDiagnostic::toJson() const
{
    using json = nlohmann::ordered_json;
    json js;
    js["positive_definite"]        = positiveDefinite;
    js["eigenvalue_min"]           = eigenvalueMin;
    js["eigenvalue_max"]           = eigenvalueMax;
    js["max_newton_step_fraction"] = maxNewtonStepFraction.has_value() ? json(*maxNewtonStepFraction) : json();
    return js;
}

HessianDiagnostic diagnoseHessian(const std::vector<std::vector<double>> &hessian,
                                  const std::vector<double>              &gradient,
                                  const std::vector<double>              &range)
{
    HessianDiagnostic hd;
    auto n = static_cast<Eigen::Index>(hessian.size());
    GMX_RELEASE_ASSERT(gradient.size() == hessian.size() && range.size() == hessian.size(),
                       "Hessian, gradient and range should have the same size");
    if (n == 0)
    {
        return hd;
    }
    // In coordinates scaled by the range, y_i = x_i/r_i
    Eigen::MatrixXd S(n, n);
    Eigen::VectorXd sg(n);
    for (Eigen::Index i = 0; i < n; i++)
    {
        GMX_RELEASE_ASSERT(hessian[i].size() == hessian.size(), "The Hessian should be square");
        for (Eigen::Index j = 0; j < n; j++)
        {
            S(i, j) = hessian[i][j]*range[i]*range[j];
        }
        sg(i) = gradient[i]*range[i];
    }
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
    const auto &eval = solver.eigenvalues();
    hd.eigenvalueMin = eval.minCoeff();
    hd.eigenvalueMax = eval.maxCoeff();
    double scale     = std::max(std::abs(hd.eigenvalueMin), std::abs(hd.eigenvalueMax));
    hd.positiveDefinite = hd.eigenvalueMin > determinantPrecision*scale;
    if (hd.positiveDefinite)
    {
        // The Newton step in scaled coordinates, -S^-1 g
        const auto     &evec = solver.eigenvectors();
        const Eigen::VectorXd step = -evec*((evec.transpose()*sg).cwiseQuotient(eval));
        hd.maxNewtonStepFraction = step.cwiseAbs().maxCoeff();
    }
    return hd;
}

}
