/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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

#include "bayes.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functional>
#include <string>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "memory_check.h"
#include "regression.h"


namespace alexandria
{

void Sensitivity::computeForceConstants(FILE *fp)
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
    else if (fp)
    {
        fprintf(fp, "Not enough parameters %d to do sensitivty analysis\n",
                static_cast<int>(p_.size()));
    }
}

void Sensitivity::print(FILE *fp, const std::string &label)
{
    if (fp)
    {
        fprintf(fp, "Sensitivity %s Fit to parabola: a %10g b %10g c %10g\n",
                label.c_str(), a_, b_, c_);
        for(int i = 0; i < static_cast<int>(p_.size()); ++i)
        {
            fprintf(fp, "    p[%d] %g chi2[%d] %g\n", i, p_[i], i, chi2_[i]);
        }
        if (a_ != 0.0)
        {
            double p_min = -b_/(2.0*a_);
            double chi2_min = a_*p_min*p_min + b_*p_min + c_;
            fprintf(fp, "    pmin %g chi2min %g (estimate based on parabola)\n",
                    p_min, chi2_min);
        }
    }
}


}
