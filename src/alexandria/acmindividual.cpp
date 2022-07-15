/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#include "acmindividual.h"

#include "gromacs/fileio/xvgr.h"
#include "act/poldata/poldata_xml.h"

namespace alexandria
{

/* * * * * * * * * * * * * * * * * * * * * *
 * BEGIN: Output stuff                     *
 * * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::printParameters(FILE *fp) const
{
    if (nullptr == fp)
    {
        return;
    }
    for(size_t i = 0; i < genome_.nBase(); i++)
    {
        fprintf(fp, "  %s  %e,", sii_->paramNames()[i].c_str(), genome_.base(i));
    }
    fprintf(fp, "\n");
}

/* * * * * * * * * * * * * * * * * * * * * *
 * END: Output stuff                       *
 * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
 * BEGIN: Poldata stuff                    *
 * * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::copyGenome(const ga::Genome &genome)
{
    genome_ = genome;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Poldata stuff                       *
* * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::addParam(const real val)
{
    initialGenome_.addBase(val);
    genome_.addBase(val);
}

} //namespace alexandria
