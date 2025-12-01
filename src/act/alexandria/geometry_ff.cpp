/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include <algorithm>

#include "act/alexandria/actmol.h"
#include "act/alexandria/alex_modules.h"
#include "act/alexandria/allbondeds.h"
#include "act/alexandria/dissociation_energy.h"
#include "act/alexandria/train_utility.h"
#include "act/basics/identifier.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/molprop/molprop_util.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/memory_check.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

//! \brief Generate bond charge corrections from bonds
static void generate_bcc(ForceField *pd,
                         double   delta_eta)
{
    // Bonds should be present, so no checking
    auto &bonds = pd->findForcesConst(InteractionType::BONDS);
    auto itype  = InteractionType::BONDCORRECTIONS;
    if (!pd->interactionPresent(itype))
    {
        auto canSwap = CanSwap::No;
        ForceFieldParameterList newparam("", canSwap);
        pd->addForces(itype, newparam);
    }
    auto bcc   = pd->findForces(itype);
    bcc->clearParameters();

    auto delta_etaParam = ForceFieldParameter("eV/e", delta_eta, 0, 0, -8, 20, Mutability::Bounded, true, false);
    auto enpBounded     = ForceFieldParameter("eV", 0, 0, 0, -8, 8, Mutability::Bounded, true, false);
    auto enpFixed       = ForceFieldParameter("eV", 0, 0, 0, 0, 0, Mutability::Fixed, true, true);
    auto ptypes = pd->particleTypesConst();
    auto itpbond = InteractionType::BONDS;
    for (auto &ai : ptypes)
    {
        if (!ai.second.hasInteractionType(itpbond))
        {
            continue;
        }
        auto bi = ai.second.interactionTypeToIdentifier(itpbond).id();
        for (auto &aj : ptypes)
        {
            if (!aj.second.hasInteractionType(itpbond))
            {
                continue;
            }
            auto bj = aj.second.interactionTypeToIdentifier(itpbond).id();
            const double bondorders[] = { 1, 1.5, 2, 3 };
            const size_t nBondorder   = std::extent<decltype(bondorders)>::value;
            for(size_t bb = 0; bb < nBondorder; bb++)
            {
                Identifier bondId({ bi, bj }, { bondorders[bb] }, bonds.canSwap());
                if (bonds.parameterExists(bondId))
                {
                    auto entype = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
                    // Only mkae bcc entries if needed.
                    if (!ai.second.hasInteractionType(entype) ||
                        !aj.second.hasInteractionType(entype))
                    {
                        continue;
                    }
                    auto zi = ai.second.interactionTypeToIdentifier(entype).id();
                    auto zj = aj.second.interactionTypeToIdentifier(entype).id();
                    if (!zi.empty() && !zj.empty())
                    {
                        Identifier bccId1({ zi, zj }, { bondorders[bb] }, bcc->canSwap());
                        Identifier bccId2({ zj, zi }, { bondorders[bb] }, bcc->canSwap());
                        if (!bcc->parameterExists(bccId1) && 
                            !bcc->parameterExists(bccId2))
                        {
                            if (zi == zj)
                            {
                                bcc->addParameter(bccId1, "delta_chi", enpFixed);
                            }
                            else
                            {
                                bcc->addParameter(bccId1, "delta_chi", enpBounded);
                            }
                            bcc->addParameter(bccId1, "delta_eta", delta_etaParam);
                        }
                    }
                }
            }
        }
    }
    printf("Have generated %zu entries for bond charge correction (SQE algorithm).\n", bcc->parameters()->size());
}

/*! \brief Tool to analyze molecular geometries and update a force field
 * \param[in] argc Number of arguments on the cmd line
 * \param[in] argv The actual arguments
 * \return 0 if ok, 1 otherwise
 */
int geometry_ff(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "geometry_ff read a series of molecules and extracts average geometries from",
        "those. First atomtypes are determined and then bond-lengths, bond-angles",
        "and dihedral angles are extracted. The results are stored in the updated force field file (aff_out.xml).[PAR]",
        "The program can also generate (quite) realistic dissociation energies from",
        "experimental or QM data when the [TT]-dissoc[tt] option is given." 
    };

    std::vector<t_filenm> fnm = {
        { efXML, "-mp",  "allmols",      ffRDMULT },
        { efXML, "-ff",  "aff",          ffOPTRD },
        { efXML, "-o",   "aff_out",      ffWRITE },
        { efDAT, "-sel", "molselect",    ffOPTRD },
        { efCSV, "-de",  "dissociation", ffOPTWR }
    };

    static real                      delta_eta    = 5;
    static int                       compress    = 0;
    static int                       maxwarn     = 0;
    static int                       nBootStrap  = 0;
    static gmx_bool                  bHisto      = false;
    static gmx_bool                  bBondOrder  = true;
    static gmx_bool                  genBCC      = true;
    static gmx_bool                  bDissoc     = false;
    static gmx_bool                  strict      = true;
    static gmx_bool                  bQM         = true;
    std::vector<t_pargs> pa = {
        { "-strict", FALSE, etBOOL, {&strict},
          "Whether or not to be pedantic about the level of theory" },
        { "-maxwarn", FALSE, etINT, {&maxwarn},
          "Will only write output if number of warnings is at most this." },
        { "-dissoc",  FALSE, etBOOL, {&bDissoc},
          "Derive dissociation energy from the enthalpy of formation. If not chosen, the dissociation energy will be read from the aff.xml file." },
        { "-qm", FALSE, etBOOL, {&bQM},
          "Usa data from quantum chemistry to determine the dissociation energies." },
        { "-bootstrap", FALSE, etINT, {&nBootStrap},
          "Use bootstrap analysis for determining the uncertainty in the dissocation energy. If the value is less than 2 no bootstrapping will be done." },
        { "-histo", FALSE, etBOOL, {&bHisto},
          "Print (hundreds of) xvg files containing histograms for bonds, angles and dihedrals" },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" },
        { "-bondorder", FALSE, etBOOL, {&bBondOrder},
          "Make separate bonds for different bond orders" },
        { "-gen_bcc", FALSE, etBOOL, {&genBCC},
          "Re-generate bond charge corrections based on the list of bonds" },
        { "-delta_eta", FALSE, etREAL, {&delta_eta},
          "Default bond delta_eta when generating bond charge corrections based on the list of bonds" },
    };

    time_t                           my_t;
    gmx_output_env_t                *oenv      = nullptr;
    ForceField                       pd;
    MolSelect                        gms;
    std::vector<alexandria::MolProp> mp;
    std::string                      method, basis;
    AllBondeds bonds;
    MsgHandler msghandler;
    msghandler.addOptions(&pa, &fnm, "geometry_ff");
    bonds.addOptions(&pa);    
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, fnm.size(), fnm.data(),
                           pa.size(), pa.data(), desc.size(), desc.data(),
                           0, nullptr, &oenv))
    {
        return 0;
    }
    CommunicationRecord cr;
    cr.init(cr.size());
    msghandler.optionsFinished(fnm, &cr);

    auto tw = msghandler.tw();
    print_memory_usage(debug);
    time(&my_t);
    tw->writeStringFormatted("# This file was created %s", ctime(&my_t));
    tw->writeStringFormatted("# The Alexandria Chemistry Toolkit.\n#\n");

    auto selfile = opt2fn_null("-sel", fnm.size(), fnm.data());
    if (selfile)
    {
        gms.read(selfile);
        print_memory_usage(debug);
        if (gms.nMol() == gms.countDimers())
        {
            msghandler.fatal("Selection file contains dimers instead of monomers.");
        }
        printf("There are %d molecules in the selection file %s.\n",
               (gms.count(iMolSelect::Train) + gms.count(iMolSelect::Test)), selfile);
        tw->writeStringFormatted("# There are %d molecules.\n#\n", (gms.count(iMolSelect::Train) + gms.count(iMolSelect::Test)));
    }
    else
    {
        printf("Will use all molecules in the molprop file(s) as the selection.\n");
    }
    /* Read standard atom properties */
    print_memory_usage(debug);

    /* Read ForceField file */
    auto myff = opt2fn_null("-ff", fnm.size(), fnm.data());
    if (!myff)
    {
        msghandler.fatal("Please pass me a force field file name");
    }
    try
    {
        readForceField(myff, &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    print_memory_usage(debug);

    // This a hack to prevent that no bonds will be found to shells.
    bool polar = pd.polarizable();
    pd.setPolarizable(false);

    /* Read Molprops */
    merge_xml(&msghandler, opt2fns("-mp", fnm.size(), fnm.data()), &mp);
    print_memory_usage(debug);
    if (!selfile)
    {
        // Extract molecule names from the molprops
        int index = 0;
        for (const auto &m : mp)
        {
            gms.addOne(m.getIupac(), index++, iMolSelect::Train);
        }
    }
    if (!msghandler.ok())
    {
        fprintf(stderr, "Too many warnings. Terminating.\n");
        return 0;
    }
    std::vector<ACTMol> actmols;
    bonds.extractGeometries(&msghandler, &mp, &actmols, &pd, gms);
    
    print_memory_usage(debug);
    if (bHisto)
    {
        bonds.writeHistogram(oenv);
    }
    bonds.updateForceField(tw, &pd);
    pd.setPolarizable(polar);
    if (genBCC)
    {
        generate_bcc(&pd, delta_eta);
    }
    print_memory_usage(debug);
    if (bDissoc)
    {
        iqmType iqm = iqmType::Exp;
        if (bQM)
        {
            iqm = iqmType::QM;
        }
        double rmsd = getDissociationEnergy(&msghandler, &pd, &actmols, iqm,
                                            opt2fn_null("-de",  fnm.size(), fnm.data()), 
                                            nBootStrap);
        tw->writeStringFormatted("Root mean square deviation %.1f kJ/mol\n", rmsd);
    }
    pd.updateTimeStamp();
    pd.updateCheckSum();
    writeForceField(opt2fn("-o", fnm.size(), fnm.data()), &pd, compress);
    bonds.writeSummary(tw);
    print_memory_usage(debug);
    return 0;
}

} // namespace alexandria
