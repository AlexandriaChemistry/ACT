/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019
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

#include "optparam.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functional>
#include <string>
#include <vector>

#include "gromacs/random.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

namespace alexandria
{

void OptParam::add_pargs(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] = {
        { "-maxiter", FALSE, etINT, {&maxiter_},
          "Max number of iterations for optimization" },
        { "-temp",    FALSE, etREAL, {&temperature_},
          "'Temperature' for the Monte Carlo simulation" },
        { "-anneal", FALSE, etBOOL, {&anneal_},
          "Use annealing in Monte Carlo simulation." },
        { "-seed",   FALSE, etINT,  {&seed_},
          "Random number seed. If zero, a seed will be generated." },
        { "-step",  FALSE, etREAL, {&step_},
          "Step size in parameter optimization. Is used as a fraction of the starting value, should be less than 10%. At each reinit step the step size is updated." }
    };
    for (size_t i = 0; i < asize(pa); i++)
    {
        pargs->push_back(pa[i]);
    }
}

void OptParam::Init(const char             *xvgconv,
                    const char             *xvgepot,
                    const gmx_output_env_t *oenv)
{
    xvgconv_     = xvgconv;
    xvgepot_     = xvgepot;
    oenv_        = oenv;
}

double OptParam::computeBeta(int iter)
{
    double temp = temperature_;
    if (anneal_)
    {
        if (iter >= maxiter_)
        {
            temp = 0;
        }
        else
        {
            temp = temperature_*(1.0 - iter/(maxiter_ + 1.0));
        }
    }
    return 1/(BOLTZ*temp);
}

double OptParam::computeBeta(int maxiter, int iter, int ncycle)
{
    double temp = temperature_;
    if (anneal_)
    {
        if (iter >= maxiter_)
        {
            temp = 0;
        }
        else
        {
            temp = (0.5*temperature_)*((exp(-iter/(0.2*(maxiter+1)))) * (1.1 + cos((ncycle*M_PI*iter)/(maxiter+1))));
        }
    }
    return 1/(BOLTZ*temp);
}

}
