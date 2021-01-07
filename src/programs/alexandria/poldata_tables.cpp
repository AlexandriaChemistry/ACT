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
    std::string      longbuf;
    CompositionSpecs cs;
    
    auto pol   = pd->polarizable();
    auto qt    = pd->findForcesConst(InteractionType::CHARGEDISTRIBUTION);

    std::string model = gmx::formatString("ACM-%s", pol ? "pg" : "g");
    std::string polstring;
    if (pol)
    {    
        lt.setColumns("lccccc");
        polstring.assign("Atomic polarizability ($\\alpha$) in {\\AA}$^3$. ");
    }
    else
    {
        lt.setColumns("lcccc");
    }
    longbuf = gmx::formatString("Optimized parameters for the Alexandria Charge Model, %s. Atomic electronegativity $\\chi$ (eV), atomic hardness $\\eta$ (eV), exponent of the charge density function $\\beta$ in nm$^{-1}$. %sThe uncertainty in each parameter is given by $\\sigma$.", model.c_str(), polstring.c_str());
    lt.setCaption(longbuf.c_str());
    lt.setLabel("eemprop");
    longbuf.assign("Type & $\\chi$($\\sigma$) & $\\eta$($\\sigma$) & $\\beta$($\\sigma$)"); 
    if (pol)
    {
        longbuf.append("& $\\alpha$($\\sigma$)");
    }
    lt.addHeadLine(longbuf.c_str());
    lt.printHeader();
}

void alexandria_eemprops_table(FILE           *fp,
                               const Poldata  *pd)
{
    LongTable  lt(fp, false, nullptr);

    eemprops_zeta_header(lt, pd);
    auto itEem   = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
    auto eep     = pd->findForcesConst(itEem);
    auto itCdist = InteractionType::CHARGEDISTRIBUTION;
    auto cdist   = pd->findForcesConst(itCdist);
    auto itPol   = InteractionType::POLARIZATION;
    auto polariz = pd->findForcesConst(itPol);
    auto pol     = pd->polarizable();
    for (const auto &pt : pd->particleTypesConst())
    {
        if (pt.mass() == 0)
        {
            continue;
        }
        std::string line = pt.id().id() + " ";

        if (pt.hasInteractionType(itEem))
        {
            auto id  = pt.interactionTypeToIdentifier(itEem);
            if (eep.parameterExists(id))
            {
                auto chi = eep.findParameterTypeConst(id, "chi");
                auto jaa = eep.findParameterTypeConst(id, "jaa");
                line.append(gmx::formatString("& %0.3f(%0.3f) & %0.3f(%0.3f) ",
                                              chi.value(), chi.uncertainty() + 0.005,
                                              jaa.value(), jaa.uncertainty() + 0.005));
            }
            else
            {
                line.append(" & &");
            }
        }
        else
        {
            line.append(" & &");
        }
        if (pt.hasInteractionType(itCdist))
        {
            auto id   = pt.interactionTypeToIdentifier(itCdist);
            
            if (cdist.parameterExists(id))
            {
                auto zeta = cdist.findParameterTypeConst(id, "zeta");
                line.append(gmx::formatString("& %0.4f(%0.3f) ",
                                              zeta.value(), zeta.uncertainty() + 0.005));
            }
            else
            {
                line.append(" &");
            }
        }
        else
        {
            line.append(" &");
        }
        if (pol)
        { 
            if (pt.hasInteractionType(itPol))
            {
                auto id    = pt.interactionTypeToIdentifier(itPol);
                if (polariz.parameterExists(id))
                {
                    auto alpha = polariz.findParameterTypeConst(id, "alpha");
                    line.append(gmx::formatString("& %0.5f(%0.4f) ",
                                                  alpha.value(), alpha.uncertainty() + 0.005));
                }
                else
                {
                    line.append(" &");
                }
            }
            else
            {
                line.append(" &");
            }
        }
        lt.printLine(line);
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
