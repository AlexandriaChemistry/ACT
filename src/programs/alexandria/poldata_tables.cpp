/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
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

#include "poldata_tables.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/cstringutil.h"

#include "categories.h"
#include "chargemodel.h"
#include "composition.h"
#include "latex_util.h"
#include "poldata_low.h"

namespace alexandria
{

static void eemprops_zeta_header(LongTable &lt,
                                 const Poldata  *pd)
{
    char             longbuf[STRLEN];
    CompositionSpecs cs;
    
    auto  iModel = pd->getChargeModel();

    lt.setColumns("lccc");

    snprintf(longbuf, STRLEN, "The optimized parameters for the Alexandria charge model, %s.The atomic electronegativity and absolute hardness are represented by $\\chi$ and $\\eta$, respectively, in eV. The exponent of the charge density function is represented by $\\beta$ in nm$^{-1}$. The atomic polarizability is represented by $\\alpha$ in {\\AA}$^3$. The uncertainty in each paramter is represented by $\\sigma$", getEemtypeName(iModel));
    lt.setCaption(longbuf);
    lt.setLabel("eemprop");
    snprintf(longbuf, STRLEN, "Alexandria Type & $\\chi$($\\sigma$) & $\\eta$($\\sigma$) & $\\beta$($\\sigma$) & $\\alpha$($\\sigma$)");
    lt.addHeadLine(longbuf);
    lt.printHeader();
}

void alexandria_eemprops_table(FILE           *fp,
                               const Poldata  *pd)
{
    double     alpha       = 0;
    double     alpha_sigma = 0;
            
    char       longbuf[STRLEN];
    LongTable  lt(fp, false, nullptr);

    eemprops_zeta_header(lt, pd);
    for (auto eem = pd->BeginEemprops(); eem < pd->EndEemprops(); eem++)
    {
        if (eem != pd->EndEemprops())
        {
            alpha       = 0;
            alpha_sigma = 0;
            
            auto nzeta = eem->getNzeta();
            auto atype = pd->ztype2atype(eem->getName());
            
            pd->getZtypePol(eem->getName(), &alpha, &alpha_sigma);
            
            snprintf(longbuf, STRLEN, "%s & %0.5f (%0.3f) & %0.5f (%0.3f) & %0.5f (%0.3f) & %0.5f (%0.3f)",
                     atype.c_str(),
                     eem->getChi0(),
                     eem->getChi0_sigma() + 0.005,
                     eem->getJ0(),
                     eem->getJ0_sigma() + 0.005,
                     eem->getZeta(nzeta-1),
                     my_atof(gmx::splitString(eem->getZeta_sigma()).back().c_str(), "zeta") + 0.005,
                     alpha,
                     alpha_sigma + 0.005);
            lt.printLine(longbuf);
        }
    }
    lt.printFooter();
    fflush(fp);
}

} //namespace
