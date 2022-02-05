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

#include <cinttypes>
#include <cstdio>
#include <cstdlib>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/cstringutil.h"

#include "categories.h"
#include "chargemodel.h"
#include "poldata_low.h"
#include "utility/latex_util.h"

namespace alexandria
{

void alexandria_subtype_table(FILE          *fp,
                              const Poldata *pd)
{
    LongTable  lt(fp, false, nullptr);

    lt.setColumns("lcccc");
    auto longbuf = gmx::formatString("Mapping from atom type to subtypes: Alexandria Charge Model (ACM), Charge distribution (Beta), Bond and Polarizability (Pol) types.");
    lt.setCaption(longbuf.c_str());
    lt.setLabel("subtype");
    auto header = gmx::formatString("Atom type & ACM & Beta & Bond & Pol");
    lt.addHeadLine(header.c_str());
    lt.printHeader();
    std::vector<std::string> subtypes = { "acmtype", "zetatype", "bondtype", "poltype" };
    for (const auto &pt : pd->particleTypesConst())
    {
        if (pt.mass() == 0)
        {
            continue;
        }
        auto line = pt.id().id();
        for(const auto &st : subtypes)
        {
            if (pt.hasOption(st))
            {
                line.append(gmx::formatString(" & %s",
                                              pt.optionValue(st).c_str()));
            }
            else
            {
                line.append(" & -");
            }
        }
        lt.printLine(line);
    }
    lt.printFooter();
}

static void alexandria_itype_table(FILE           *fp,
                                   const Poldata  *pd,
                                   InteractionType itype)
{
    LongTable  lt(fp, false, nullptr);
    lt.setCaption(gmx::formatString("Parameters for %s. Average value(s) is/are given, with number of data points N and standard deviation $\\sigma$.", 
                                    interactionTypeToDescription(itype).c_str()).c_str());
    lt.setLabel(interactionTypeToString(itype).c_str());
    auto        eep      = pd->findForcesConst(itype);
    bool        first    = true;
    std::string header;
    int         ncolumns = 1;
    for (const auto &ep : eep.parametersConst())
    {
        std::string line = ep.first.id();
        if (first)
        {
            header = "Type";
        }
        int ntrain = 0;
        for (const auto &fp : ep.second)
        {
            // Round upwards the sigma values.
            if (fp.second.ntrain() > 0)
            {
                line += gmx::formatString(" & %.2f(%d, %.2f)",
                                          fp.second.value(),
                                          fp.second.ntrain(),
                                          fp.second.uncertainty()+0.005);
            }
            else
            {
                line += " & - ";
            }
            ntrain += fp.second.ntrain();
            if (first)
            {
                header   += gmx::formatString(" & %s (N, $\\sigma$)", fp.first.c_str());
                ncolumns += 1;
            }
        }
        if (first)
        {
            lt.setColumns(ncolumns);
            lt.addHeadLine(header.c_str());
            lt.printHeader();
            first = false;
        }
        if (ntrain > 0)
        {
            lt.printLine(line);
        }
    }
    lt.printFooter();
}

void alexandria_eemprops_table(FILE           *fp,
                               const Poldata  *pd)
{
    std::vector<InteractionType> itypes = {
        InteractionType::ELECTRONEGATIVITYEQUALIZATION,
        InteractionType::CHARGEDISTRIBUTION,
        InteractionType::POLARIZATION,
        InteractionType::BONDCORRECTIONS
    };
    for (const auto &itype : itypes)
    { 
        if (pd->interactionPresent(itype))
        {
            auto fs = pd->findForcesConst(itype);
            if (fs.numberOfParameters() > 0)
            {
                alexandria_itype_table(fp, pd, itype);
            }
        }
    }
}

void alexandria_charge_table(FILE           *fp,
                             const Poldata  *pd)
{
    LongTable   lt(fp, false, nullptr);
    std::string qq("charge");
    lt.setCaption(gmx::formatString("Parameters for %s. Average value(s) is/are given, with number of data points N and standard deviation $\\sigma$.", qq.c_str()).c_str());
    lt.setLabel(qq.c_str());
    lt.setColumns(4);
    lt.addHeadLine("Particle type & q & $\\Delta$q & N");
    lt.printHeader();
    for (const auto &ptype : pd->particleTypesConst())
    {
        if (ptype.hasParameter(qq))
        {
            auto pp = ptype.parameterConst(qq);
            if ((pp.mutability() == Mutability::Free ||
                 pp.mutability() == Mutability::Bounded) &&
                pp.ntrain() > 0)
            {
                auto line = gmx::formatString("%s & %.4f & %.4f & %d",
                                              ptype.id().id().c_str(),
                                              pp.value(), pp.uncertainty(),
                                              pp.ntrain());
                lt.printLine(line);
            }
        }
    }
    lt.printFooter();
}

} //namespace
