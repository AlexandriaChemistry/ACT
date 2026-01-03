/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
    
#include "import_utils.h"

#include "act/basics/msg_handler.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

void updateFragmentFromInchi(MsgHandler           *msg_handler,
                             const AlexandriaMols &amols,
                             const std::string    &inchi,
                             Fragment             *fptr)
{
    auto amol  = amols.findInChi(inchi);
    if (nullptr != amol)
    {
        fptr->setInchi(inchi);
        fptr->setIupac(amol->iupac);
        fptr->setCharge(amol->charge);
        fptr->setMass(amol->mass);
        fptr->setFormula(amol->formula);
    }
    else
    {
        msg_handler->msg(ACTStatus::Warning,
                         gmx::formatString("Could not find '%s' in database", inchi.c_str()));
        fptr->setInchi(inchi);
    }
}

std::vector<AtomBondtypeEntry> getAtomBondtypeDB(gmx_unused const std::string &dbname)
{
    // Note that the order is important, the vector will be checked from top to bottom
    std::vector<AtomBondtypeEntry> abe = {
        { "sulfate",
          "[#16](=[#8])(=[#8])(-[#8-])-[#8-]",
          -2, 1,
          { "s3", "o2", "o2", "o2", "o2" },
          { { 0, 1, 1.5 }, { 0, 2, 1.5 }, { 0, 3, 1.5 }, { 0, 4, 1.5 } }
        },
        { "phosphate", // hypervalent form
          "[#15D4]([#8D1])([#8D1])([#8-,#8D1])([#8-,#8D1])",
          -3, 1,
          { "p3", "o3", "o3", "o3", "o3" },
          { { 0, 1, 1.5 }, { 0, 2, 1.5 }, { 0, 3, 1.5 }, { 0, 4, 1.5 } }
        },
        { "phosphate2", // ion form, PO4(3-)
          "[#8-]-[#15](-[#8-])(-[#8-])=[#8]",
          -3, 1,
          { "o3", "p3", "o3", "o3", "o3" },
          { { 0, 1, 1.5 }, { 1, 2, 1.5 }, { 1, 3, 1.5 }, { 1, 4, 1.5 } }
        },
        { "phosphate2", // ion form, XPO3(2-)
          "[#8-]-[#15](-[#8-])(-[#8-])",
          -2, 1,
          { "o3", "p3", "o3", "o3" },
          { { 0, 1, 1.5 }, { 1, 2, 1.5 }, { 1, 3, 1.5 } }
        },
        { "phosphate3", // ion form, X2PO2(1-)
          "[#8-]-[#15](-[#8-])",
          -1, 1,
          { "o3", "p3", "o3" },
          { { 0, 1, 1 }, { 1, 2, 1 } }
        },
        { "nitro1",
          "[#7+](-[#8-])=[#8]",
          0, 1,
          { "n2", "o2", "o2" },
          { { 0, 1, 1.5 }, { 0, 2, 1.5 } }
        },
        { "nitro2",
          "[#7+](-[#8])=[#8]",
          0, 1,
          { "n2", "o2", "o2" },
          { { 0, 1, 1.5 }, { 0, 2, 1.5 } }
        }
    };
    return abe;
}

void writeAtomBondtypeDB(gmx_unused const std::string                    &filenm,
                         gmx_unused const std::vector<AtomBondtypeEntry> &db)
{
}

} // namespace alexandria
