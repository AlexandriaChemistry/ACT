/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

#include "forcefield_tables.h"

#include <cinttypes>
#include <cstdio>
#include <cstdlib>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/cstringutil.h"

#include "act/basics/chargemodel.h"
#include "symcharges.h"
#include "act/utility/latex_util.h"

namespace alexandria
{

void ForceFieldTable::subtype_table(const std::string &info)
{
    LongTable  lt(fp_, false, nullptr);

    lt.setColumns("lcccc");
    auto longbuf = gmx::formatString("Mapping from atom type to subtypes: Alexandria Charge Model (ACM), Charge distribution (Beta), Bond and Polarizability (Pol) types. %s", info.c_str());
    lt.setCaption(longbuf.c_str());
    lt.setLabel("subtype");
    auto header = gmx::formatString("Atom type & ACM & Beta & Bond & Pol");
    lt.addHeadLine(header.c_str());
    lt.printHeader();
    std::vector<std::string> subtypes = { "acmtype", "zetatype", "bondtype", "poltype" };
    for (const auto &pt : pd_->particleTypesConst())
    {
        if (pt.second.mass() == 0)
        {
            continue;
        }
        auto line = pt.second.id().id();
        for(const auto &st : subtypes)
        {
            if (pt.second.hasOption(st))
            {
                line.append(gmx::formatString(" & %s",
                                              pt.second.optionValue(st).c_str()));
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

void ForceFieldTable::itype_table(InteractionType    itype,
                                  const std::string &info)
{
    if (!pd_->interactionPresent(itype))
    {
        return;
    }
    auto eep = pd_->findForcesConst(itype);
    // Check whether there is anything to print.
    int nparam = 0;
    for (const auto &p : eep.parametersConst())
    {
        for (const auto &pp : p.second)
        {
            if (pp.second.ntrain() > ntrain_)
            {
                nparam += 1;
            }
        }
    }
    if (nparam == 0)
    {
        return;
    }
    LongTable  lt(fp_, false, nullptr);
    if (printSigma_)
    {
        lt.setCaption(gmx::formatString("Parameters for %s. Average value(s) is/are given, with number of data points N and standard deviation $\\sigma$. %s",
                                        interactionTypeToDescription(itype).c_str(),
                                        info.c_str()).c_str());
    }
    else
    {
        lt.setCaption(gmx::formatString("Parameters for %s. %s",
                                        interactionTypeToDescription(itype).c_str(),
                                        info.c_str()).c_str());
    }
    lt.setLabel(interactionTypeToString(itype).c_str());
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
        unsigned int ntrain = 0;
        for (const auto &ffp : ep.second)
        {
            if (ffp.second.mutability() == Mutability::Dependent)
            {
                continue;
            }
            // Round upwards the sigma values.
            if (ffp.second.ntrain() >= ntrain_)
            {
                if (printSigma_)
                {
                    line += gmx::formatString(" & %.3f(%d, %.3f)",
                                              ffp.second.value(),
                                              ffp.second.ntrain(),
                                              ffp.second.uncertainty()+0.005);
                }
                else
                {
                    line += gmx::formatString(" & %.3f", ffp.second.value());
                }
            }
            else
            {
                line += " & - ";
            }
            ntrain += ffp.second.ntrain();
            if (first)
            {
                header   += gmx::formatString(" & %s (N, $\\sigma$)", ffp.first.c_str());
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
        if (ntrain >= ntrain_)
        {
            lt.printLine(line);
        }
    }
    lt.printFooter();
}

void ForceFieldTable::eemprops_table(const std::string &info)
{
    std::vector<InteractionType> itypes = {
        InteractionType::ELECTRONEGATIVITYEQUALIZATION,
        InteractionType::ELECTROSTATICS,
        InteractionType::POLARIZATION,
        InteractionType::INDUCTIONCORRECTION,
        InteractionType::BONDCORRECTIONS
    };
    for (const auto &my_itype : itypes)
    { 
        if (pd_->interactionPresent(my_itype))
        {
            auto fs = pd_->findForcesConst(my_itype);
            if (fs.numberOfParameters() > 0)
            {
                itype_table(my_itype, info);
            }
        }
    }
}

} //namespace
