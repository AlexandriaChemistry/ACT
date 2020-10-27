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
    
    auto qt    = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    auto iType = name2ChargeType(qt.optionValue("chargetype"));

    lt.setColumns("lcccc");

    snprintf(longbuf, STRLEN, "The optimized parameters for the Alexandria charge model, %s.The atomic electronegativity and absolute hardness are represented by $\\chi$ and $\\eta$, respectively, in eV. The exponent of the charge density function is represented by $\\beta$ in nm$^{-1}$. The atomic polarizability is represented by $\\alpha$ in {\\AA}$^3$. The uncertainty in each paramter is represented by $\\sigma$", chargeTypeName(iType).c_str());
    lt.setCaption(longbuf);
    lt.setLabel("eemprop");
    snprintf(longbuf, STRLEN, "Alexandria Type & $\\chi$($\\sigma$) & $\\eta$($\\sigma$) & $\\beta$($\\sigma$) & $\\alpha$($\\sigma$)");
    lt.addHeadLine(longbuf);
    lt.printHeader();
}

void alexandria_eemprops_table(FILE           *fp,
                               const Poldata  *pd)
{
    LongTable  lt(fp, false, nullptr);

    eemprops_zeta_header(lt, pd);
    auto eep = pd->findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
    for (auto &eem : eep.parametersConst())
    {
        auto chi  = &eem.second.find("chi")->second;
        auto jaa  = &eem.second.find("jaa")->second;
        auto zeta = &eem.second.find("zeta")->second;
        auto buf  = gmx::formatString("%s & %0.5f (%0.3f) & %0.5f (%0.3f) & %0.5f (%0.3f)",
                                     eem.first.id().c_str(),
                                     chi->value(), chi->uncertainty() + 0.005,
                                     jaa->value(), jaa->uncertainty() + 0.005,
                                     zeta->value(), zeta->uncertainty() + 0.005);
        lt.printLine(buf.c_str());
    }
    lt.printFooter();
    fflush(fp);
}

void alexandria_eemprops_corr(const Poldata  *pd,
                              FILE           *fp)
{

    gmx_stats_t  chi_eta    = gmx_stats_init();
    gmx_stats_t  chi_zeta   = gmx_stats_init();
    gmx_stats_t  eta_zeta   = gmx_stats_init();
    
    real ce = 0;
    real cz = 0;
    real ez = 0;
    
    auto eep = pd->findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION);
    for (auto &eem : eep.parametersConst())
    {
        auto chi  = eem.second.find("chi")->second.value();
        auto jaa  = eem.second.find("jaa")->second.value();
        auto zeta = eem.second.find("zeta")->second.value();
        
        gmx_stats_add_point(chi_eta,  chi, jaa, 0, 0);
        gmx_stats_add_point(chi_zeta, chi, zeta, 0, 0);
        gmx_stats_add_point(eta_zeta, jaa, zeta, 0, 0);
    }
    
    gmx_stats_get_corr_coeff(chi_eta,    &ce);
    gmx_stats_get_corr_coeff(chi_zeta,   &cz);
    gmx_stats_get_corr_coeff(eta_zeta,   &ez);
    
    fprintf(fp, "\nCorrelation coefficient between eemprop parameters:\n");
    fprintf(fp, "Absolute Eelectronegativity and Absolute Hardness: %0.3f \n", 100*ce);
    fprintf(fp, "Absolute Eelectronegativity and Exponent of Charge Density: %0.3f \n", 100*cz);
    fprintf(fp, "Absolute Hardness and Exponent of Charge Density: %0.3f \n", 100*ez);    
}

} //namespace
