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

#include "bayes.h"

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
    if (p_.size() >= 3)
    {
        MatrixWrapper M(3, p_.size());
        std::vector<double> solution;
        solution.resize(3, 0.0);
        for (size_t i = 0; i < p_.size(); ++i)
        {
            M.set(0, i, p_[i]*p_[i]);
            M.set(1, i, p_[i]);
            M.set(2, i, 1.0);
        }
        auto result = M.solve(chi2_, &solution);
        if (result == 0)
        {
            a_ = solution[0];
            b_ = solution[1];
            c_ = solution[2];
        }
    }
    else if (tw)
    {
        tw->writeStringFormatted("Not enough parameters %zu to do sensitivty analysis\n",
                                 p_.size());
    }
}

void Sensitivity::print(gmx::TextWriter      *tw,
                        alexandria::JsonTree *jtree,
                        const std::string    &index,
                        const std::string    &label)
{
    double p_min    = 0;
    double chi2_min = 0;
    if (tw)
    {
        tw->writeStringFormatted("Sensitivity %s Fit to parabola: a %10g b %10g c %10g\n",
                                 label.c_str(), a_, b_, c_);
        for(size_t i = 0; i < p_.size(); ++i)
        {
            tw->writeStringFormatted("    p[%zu] %g chi2[%zu] %g\n", i, p_[i], i, chi2_[i]);
        }
        if (a_ != 0.0)
        {
            p_min = -b_/(2.0*a_);
            chi2_min = a_*p_min*p_min + b_*p_min + c_;
            tw->writeStringFormatted("    pmin %g chi2min %g (estimate based on parabola)\n",
                                     p_min, chi2_min);
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
        abc.addObject("p_orig", p_[1]);
        abc.addObject("chi2_orig", chi2_[1]);
        if (a_ != 0)
        {
            abc.addObject("p_min", p_min);
            abc.addObject("chi2_min", chi2_min);
        }
        jtree->addObject(abc);
    }
}


}
