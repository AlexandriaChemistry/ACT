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
#include <map>
#include <string>
#include <vector>

#include "act/utility/jsontree.h"
#include "act/utility/regression.h"
#include "act/utility/stringutil.h"
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
        A_ = solution[0];
        B_ = solution[1];
        C_ = solution[2];
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
    if (A_ > 0 && original_ < p_.size())
    {
        double pmin = std::clamp(pMin(), lower, upper);
        converged_  = std::abs(p_[original_] - pmin) <= tolerance;
    }
    return converged_;
}

void Sensitivity::print(gmx::TextWriter      *tw,
                        alexandria::JsonTree *jtree,
                        const std::string    &index,
                        const std::string    &label)
{
    // The extremum of the parabola, also needed without a text writer
    double p_min    = pMin();
    double chi2_min = chi2Min();
    if (tw)
    {
        tw->writeStringFormatted("Sensitivity %s Fit to parabola: a %10g b %10g c %10g\n",
                                 label.c_str(), a_, b_, c_);
        for(size_t i = 0; i < p_.size(); ++i)
        {
            tw->writeStringFormatted("    p[%zu] %g chi2[%zu] %g\n", i, p_[i], i, chi2_[i]);
        }
        if (A_ != 0)
        {
            tw->writeStringFormatted("    pmin %g chi2min %g (estimate based on parabola)\n",
                                     p_min, chi2_min);
        }
        if (checked_)
        {
            tw->writeStringFormatted("    %s\n", converged_ ? "converged" : "not converged");
        }
    }
    if (jtree)
    {
        JsonTree abc(index);
        // A bit of a hack since the input label contains two pieces of information
        auto pn = split(label, ' ');
        abc.addObject("particle", pn[0]);
        abc.addObject("parameter", pn[1]);
        abc.addObject("a", a_);
        abc.addObject("b", b_);
        abc.addObject("c", c_);
        abc.addObject("p_orig", p_[original_]);
        abc.addObject("chi2_orig", chi2_[original_]);
        if (A_ != 0)
        {
            abc.addObject("p_min", p_min);
            abc.addObject("chi2_min", chi2_min);
        }
        if (checked_)
        {
            abc.addObject("converged", converged_ ? "yes" : "no");
        }
        jtree->addObject(abc);
    }
}


}
