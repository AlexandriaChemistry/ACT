/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/real.h"

#include "alex_modules.h"
#include "allbondeds.h"
#include "dissociation_energy.h"
#include "act/basics/identifier.h"
#include "act/utility/memory_check.h"
#include "act/molprop/molprop_util.h"
#include "mymol.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/stringutil.h"
#include "tuning_utility.h"

namespace alexandria
{

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
        pd->addForces(interactionTypeToString(itype), newparam);
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
        if (!ai.hasInteractionType(itpbond))
        {
            continue;
        }
        auto bi = ai.interactionTypeToIdentifier(itpbond).id();
        for (auto &aj : ptypes)
        {
            if (!aj.hasInteractionType(itpbond))
            {
                continue;
            }
            auto bj = aj.interactionTypeToIdentifier(itpbond).id();
            const double bondorders[] = { 1, 1.5, 2, 3 };
            const size_t nBondorder   = std::extent<decltype(bondorders)>::value;
            for(size_t bb = 0; bb < nBondorder; bb++)
            {
                Identifier bondId({ bi, bj }, { bondorders[bb] }, bonds.canSwap());
                if (bonds.parameterExists(bondId))
                {
                    auto entype = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
                    auto zi = ai.interactionTypeToIdentifier(entype).id();
                    auto zj = aj.interactionTypeToIdentifier(entype).id();
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
    printf("Have generated %zu entries for BCC\n", bcc->parameters()->size());
}

int geometry_ff(int argc, char *argv[])
{
    static const char               *desc[] = {
        "geometry_ff read a series of molecules and extracts average geometries from",
        "those. First atomtypes are determined and then bond-lengths, bond-angles",
        "and dihedral angles are extracted. The results are stored in the updated force field file (aff_out.xml).[PAR]",
        "The program can also generate (quite) realistic dissociation energies from",
        "experimental or QM data when the [TT]-dissoc[tt] option is given." 
    };

    t_filenm                         fnm[] = {
        { efXML, "-mp",  "allmols",      ffRDMULT },
        { efXML, "-ff",  "aff",          ffOPTRD },
        { efXML, "-o",   "aff_out",      ffWRITE },
        { efDAT, "-sel", "molselect",    ffREAD },
        { efLOG, "-g",   "geometry",     ffWRITE },
        { efCSV, "-de",  "dissociation", ffOPTWR }
    };

    const int                        NFILE       = asize(fnm);
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

    FILE                            *fp;
    time_t                           my_t;
    gmx_output_env_t                *oenv      = nullptr;
    ForceField                          pd;
    MolSelect                        gms;
    std::vector<alexandria::MolProp> mp;
    std::string                      method, basis;
    AllBondeds bonds;
    
    bonds.addOptions(&pa);    
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm,
                           pa.size(), pa.data(), asize(desc), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }

    fp                 = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");
    print_memory_usage(debug);
    time(&my_t);
    fprintf(fp, "# This file was created %s", ctime(&my_t));
    fprintf(fp, "# The Alexandria Chemistry Toolkit.\n#\n");

    auto selfile = opt2fn("-sel", NFILE, fnm);
    gms.read(selfile);
    print_memory_usage(debug);
    printf("There are %d molecules in the selection file %s.\n",
           (gms.count(iMolSelect::Train) + gms.count(iMolSelect::Test)), selfile);
    fprintf(fp, "# There are %d molecules.\n#\n", (gms.count(iMolSelect::Train) + gms.count(iMolSelect::Test)));

    /* Read standard atom properties */
    print_memory_usage(debug);

    /* Read ForceField file */
    try
    {
        readForceField(opt2fn_null("-ff", NFILE, fnm), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    print_memory_usage(debug);

    // This a hack to prevent that no bonds will be found to shells.
    bool polar = pd.polarizable();
    pd.setPolarizable(false);

    /* Read Molprops */
    auto nwarn = merge_xml(opt2fns("-mp", NFILE, fnm), &mp, nullptr, nullptr, nullptr, false);
    print_memory_usage(debug);

    if (nwarn > maxwarn)
    {
        printf("Too many warnings (%d). Terminating.\n", nwarn);
        return 0;
    }
    std::vector<MyMol> mymols;
    bonds.extractGeometries(fp, mp, &mymols, pd, gms);
    
    print_memory_usage(debug);
    if (bHisto)
    {
        bonds.writeHistogram(oenv);
    }
    bonds.updateForceField(fp, &pd);
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
        double rmsd = getDissociationEnergy(fp, &pd, &mymols, iqm,
                                            opt2fn_null("-de",  NFILE, fnm), 
                                            nBootStrap);
        fprintf(fp, "Root mean square deviation %.1f kJ/mol\n", rmsd);
    }
    pd.updateTimeStamp();
    pd.updateCheckSum();
    writeForceField(opt2fn("-o", NFILE, fnm), &pd, compress);
    bonds.writeSummary(fp);
    print_memory_usage(debug);
    gmx_ffclose(fp);
    return 0;
}

} // namespace alexandria
