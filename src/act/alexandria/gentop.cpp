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

#include <cctype>
#include <cstdlib>

#include <algorithm>

#include "act/alexandria/actmol.h"
#include "act/alexandria/alex_modules.h"
#include "act/alexandria/compound_reader.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/openmm_xml.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/utility/stringutil.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

int gentop(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "gentop generates a topology from molecular coordinates",
        "from a database or, optionally, from a coordinate file or",
        "a gaussian log file.",
        "The program assumes all hydrogens are present when defining",
        "the hybridization from the atom name and the number of bonds.[PAR]",
        "If the [TT]-oi[tt] option is set an [TT]itp[tt] file will be generated",
        "instead of a [TT]top[tt] file.",
        "A force field file should be specified using the [TT]-ff[tt] flag.",
        "Depending on the selected force field file, the force field",
        "uses polarizable or non-polarizable charges or virtual sites.[PAR]",
        "When the [TT]-openmm[tt] flag is passed, an XML file will be created",
        "that can be used to run a simulation of the system using the OpenMM",
        "software."
    };
    gmx_output_env_t *oenv;
    gmx_atomprop_t    aps;

    std::vector<t_filenm> fnm = {
        { efXML, "-ff",         "aff",       ffREAD  },
        { efTOP, "-p",          "out",       ffOPTWR },
        { efXML, "-openmm",     "out",       ffOPTWR },
        { efITP, "-oi",         "out",       ffOPTWR },
        { efPDB, "-c",          "out",       ffWRITE },
        { efNDX, "-n",          "renum",     ffOPTWR },
        { efDAT, "-q",          "qout",      ffOPTWR },
        { efCUB, "-pot",        "potential", ffOPTWR },
        { efCUB, "-ref",        "refpot",    ffOPTRD },
        { efCUB, "-diff",       "diffpot",   ffOPTWR },
        { efCUB, "-rho",        "density",   ffOPTWR },
        { efXVG, "-diffhist",   "diffpot",   ffOPTWR },
        { efXVG, "-his",        "pot-histo", ffOPTWR },
        { efXVG, "-pc",         "pot-comp",  ffOPTWR },
        { efPDB, "-pdbdiff",    "pdbdiff",   ffOPTWR }
    };

    static int                       nsymm          = 0;
    static int                       ntrain         = 1;
    static real                      spacing        = 0.01;
    static real                      border         = 0.2;
    static real                      mDrude         = 0.1;
    static char                     *symm_string    = (char *)"";
    static char                     *conf           = (char *)"minimum";
    static char                     *jobtype        = (char *)"Opt";
    //! Write shells to trajectory and coordinates
    bool                             writeShells    = false;
    static gmx_bool                  bITP           = false;
    static gmx_bool                  bAllowMissing  = false;
    static gmx_bool                  addNumbersToAtoms = true;
    static rvec                      mybox           = { 0, 0, 0 };

    std::vector<t_pargs> pa =
    {
        { "-conf",  FALSE, etSTR, {&conf},
          "Conformation of the molecule" },
        { "-ws",     FALSE, etBOOL, {&writeShells},
          "Write coordinates of shell particles to trajectory and final coordinates as well." },
        { "-box",    FALSE, etRVEC, {mybox},
          "If set will replace the simulation box in the output coordinate file." },
        { "-mDrude", FALSE, etREAL, {&mDrude},
          "Mass to use for the drude particle if any. Default is 0.1 Da" },
        { "-allowmissing", FALSE, etBOOL, {&bAllowMissing},
          "Make a topology even if there are no force field parameters for all interactions" },
        { "-nsymm", FALSE, etINT, {&nsymm},
          "Symmetry number of the molecule can be supplied here if you know there is an error in the input file" },
        { "-spacing", FALSE, etREAL, {&spacing},
          "Spacing of grid points (nm) for computing the potential (not used when a reference file is read)." },
        { "-border", FALSE, etREAL, {&border},
          "Spacing around the compound (nm) for computing the potential (not used when a reference file is read)." },
        { "-ntrain", FALSE, etINT, {&ntrain},
          "Minimum number of training data to store a parameter to OpenMM. Parameters that have Mutability::Fixed are counted as having at least one training data point. This is rather fragile, so beware of this when increasing the number. Zero means taht even untrained parameters will be exported to OpenMM, use with care!" },
        { "-symm",   FALSE, etSTR, {&symm_string},
          "Use the order given here for symmetrizing, e.g. when specifying [TT]-symm '0 1 0'[tt] for a water molecule (H-O-H) the hydrogens will have obtain the same charge. For simple groups, like methyl (or water) this is done automatically, but higher symmetry is not detected by the program. The numbers should correspond to atom numbers minus 1, and point to either the atom itself or to a previous atom." },
        { "-numberAtypes", FALSE, etBOOL, {&addNumbersToAtoms},
          "Add a number index to OpenMM atomtypes when generating output for that software." },
         { "-jobtype",  FALSE, etSTR, {&jobtype},
          "The job type used in the Gaussian calculation: Opt, Polar, SP, and etc." }
    };
    MsgHandler          msghandler;
    CompoundReader      compR;
    compR.addOptions(&pa, &fnm, &desc);
    msghandler.addOptions(&pa, &fnm, "gentop.log");
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 1;
    }
    CommunicationRecord cr;
    cr.init(cr.size());
    msghandler.optionsFinished(fnm, &cr);
    
    compR.optionsOK(&msghandler, fnm);
    if (!msghandler.ok())
    {
        return 1;
    }
    ForceField          pd;
    std::string         method, basis;

    /* Check the options */
    bITP = opt2bSet("-oi", fnm.size(), fnm.data());

    const char *gentop_fnm = opt2fn_null("-ff", fnm.size(), fnm.data());
    if (nullptr == gentop_fnm)
    {
        fprintf(stderr, "Please pass me a force field file name with the -ff option.\n");
        return 1;
    }

    /* Read standard atom properties */
    aps = gmx_atomprop_init();
    if (!msghandler.verbose())
    {
        gmx_atomprop_quiet(aps);
    }
    try
    {
        readForceField(gentop_fnm, &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    (void) pd.verifyCheckSum(stderr);

    auto &fs = pd.findForcesConst(InteractionType::ELECTROSTATICS);
    std::string my_pol;
    if (pd.polarizable())
    {
        my_pol.assign(" polarizable");
    }
    auto qType = potentialToChargeType(fs.potential());
    printf("Using%s force field file %s and charge distribution model %s\n",
           my_pol.c_str(), gentop_fnm, chargeTypeName(qType).c_str());
    msghandler.write(gmx::formatString("Read force field information. There are %zu atomtypes.\n",
                                       pd.getNatypes()));

    auto forceComp = new ForceComputer();
    std::vector<ACTMol> actmols = compR.read(&msghandler, pd, forceComp);

    int mp_index   = 1;
    std::map<std::string, std::pair<ACTMessage, std::vector<std::string>>> errors;
    for(auto actmol = actmols.begin(); actmol < actmols.end(); ++actmol)
    {
        std::vector<gmx::RVec> forces(actmol->atomsConst().size());
        std::vector<gmx::RVec> coords = actmol->xOriginal();
        forceComp->generateVsites(actmol->topology(), &coords);

        actmol->GenerateCube(&pd, coords, forceComp,
                             spacing, border,
                             opt2fn_null("-ref",      fnm.size(), fnm.data()),
                             opt2fn_null("-pc",       fnm.size(), fnm.data()),
                             opt2fn_null("-pdbdiff",  fnm.size(), fnm.data()),
                             opt2fn_null("-pot",      fnm.size(), fnm.data()),
                             opt2fn_null("-rho",      fnm.size(), fnm.data()),
                             opt2fn_null("-his",      fnm.size(), fnm.data()),
                             opt2fn_null("-diff",     fnm.size(), fnm.data()),
                             opt2fn_null("-diffhist", fnm.size(), fnm.data()),
                             oenv);

        std::string index;
        if (actmols.size() > 1)
        {
            index = gmx::formatString("%d_", mp_index);
        }
        if (!opt2bSet("-openmm", fnm.size(), fnm.data()))
        {
            std::string tfn = gmx::formatString("%s%s", index.c_str(),
                                                bITP ? ftp2fn(efITP, fnm.size(), fnm.data()) : ftp2fn(efTOP, fnm.size(), fnm.data()));
            actmol->PrintTopology(tfn.c_str(), msghandler.verbose(), &pd, forceComp,
                                  &cr, coords, method, basis, bITP);
        }
        if (opt2bSet("-c", fnm.size(), fnm.data()))
        {
            matrix box = { { 5, 0, 0 }, { 0, 5, 0 }, { 0, 0, 5 }};
            std::string cfn = gmx::formatString("%s%s", index.c_str(),
                                                opt2fn("-c", fnm.size(), fnm.data()));
            if (opt2parg_bSet("-box", pa.size(), pa.data()))
            {
                for(int m = 0; m < DIM; m++)
                {
                    box[m][m] = mybox[m];
                    }
            }
            actmol->PrintConformation(cfn.c_str(), coords, writeShells, box);
        }
        mp_index++;
    }
    if (!actmols.empty())
    {
        if (opt2bSet("-openmm", fnm.size(), fnm.data()))
        {
            writeOpenMM(opt2fn("-openmm", fnm.size(), fnm.data()),
                        &pd, actmols, mDrude, ntrain, addNumbersToAtoms);
        }
    }
    printf("\nPlease check the %s file for warnings and error messages.\n", msghandler.filename().c_str());
    if (!msghandler.ok())
    {
        status = 1;
    }
    return status;
}

} // namespace alexandria
