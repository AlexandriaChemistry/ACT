/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#include <ctype.h>
#include <stdlib.h>

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
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/smalloc.h"

#include "alex_modules.h"
#include "atype_mapping.h"
#include "babel_io.h"
#include "fill_inputrec.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "mymol.h"
#include "openmm_xml.h"
#include "act/poldata/poldata_xml.h"

namespace alexandria
{

static void print_errors(const char                     *fn,
                         const std::vector<std::string> &errors,
                         immStatus                       imm)
{
    FILE  *fp   = gmx_ffopen(fn, "w");
    fprintf(fp, "The detected error is \"%s\"\n", alexandria::immsg(imm));
    for (auto error : errors)
    {
        fprintf(fp, "Error: %s\n", error.c_str());
    }
    fclose(fp);
}

int gentop(int argc, char *argv[])
{
    static const char               *desc[] = {
        "gentop generates a topology from molecular coordinates",
        "either from a file, from a database, or from a gaussian log file.",
        "The program assumes all hydrogens are present when defining",
        "the hybridization from the atom name and the number of bonds.[PAR]",
        "If the [TT]-oi[tt] option is set an [TT]itp[tt] file will be generated",
        "instead of a [TT]top[tt] file.",
        //"When [TT]-param[tt] is set, equilibrium distances and angles",
        //"and force constants will be printed in the topology for all",
        //"interactions. The equilibrium distances and angles are taken",
        //"from the input coordinates, the force constant are set with",
        //"command line options.",
        "With the [TT]-db molecule[tt] option a file is extracted from the",
        "database from the QM calculations.",
        "An alternative to the system-wide database [TT]molprops.xml[tt]",
        "can be passed along using the [TT]-mp[tt] flag.[PAR]",
        //"If the flag [TT]-act/qgen/qgen[tt] is given, charges will be generated using the",
        //"specified algorithm. Without the flag the charges from the QM calculation",
        //"will be used.",
        "Using the [TT]-ff[tt] flag, a force field file can be specified.",
        "Depending on the selected force field file, the force field",
        "uses polarizable or non-polarizable Gaussian charges. [PAR]",
        "When the [TT]-openmm[tt] flag is passed, an XML file will be created",
        "that can be used to run a simulation of the system using the OpenMM",
        "software."
    };
    gmx_output_env_t *oenv;
    gmx_atomprop_t    aps;
    immStatus         imm;

    t_filenm                         fnm[] = {
        { efTOP, "-p",        "out",           ffOPTWR },
        { efXML, "-openmm",   "out",           ffOPTWR },
        { efITP, "-oi",       "out",           ffOPTWR },
        { efSTO, "-c",        "out",           ffWRITE },
        { efNDX, "-n",        "renum",         ffOPTWR },
        { efDAT, "-q",        "qout",          ffOPTWR },
        { efXML, "-mp",       "molprops",      ffOPTRD },
        { efXML, "-ff",       "gentop",        ffOPTRD },
        { efCUB, "-pot",      "potential",     ffOPTWR },
        { efCUB, "-ref",      "refpot",        ffOPTRD },
        { efCUB, "-diff",     "diffpot",       ffOPTWR },
        { efCUB, "-rho",      "density",       ffOPTWR },
        { efXVG, "-diffhist", "diffpot",       ffOPTWR },
        { efXVG, "-his",      "pot-histo",     ffOPTWR },
        { efXVG, "-pc",       "pot-comp",      ffOPTWR },
        { efPDB, "-pdbdiff",  "pdbdiff",       ffOPTWR },
        { efXVG, "-plot_esp", "ESPcorr",       ffOPTWR },
        { efLOG, "-g",        "gentop_errors", ffWRITE }
    };

    const  int                       NFILE          = asize(fnm);

    static int                       maxpot         = 100;
    static int                       nsymm          = 0;
    static real                      qtot           = 0;
    static real                      watoms         = 0;
    static real                      spacing        = 0.01;
    static real                      border         = 0.2;
    static char                     *molnm          = (char *)"";
    static char                     *iupac          = (char *)"";
    static char                     *dbname         = (char *)"";
    static char                     *symm_string    = (char *)"";
    static char                     *qqm            = (char *)"";
    static char                     *conf           = (char *)"minimum";
    static char                     *jobtype        = (char *)"Opt";
    static char                     *filename       = (char *)"";
    static gmx_bool                  bQsym          = false;
    static gmx_bool                  bITP           = false;
    static gmx_bool                  bUsePDBcharge  = false;
    static gmx_bool                  bVerbose       = false;
    static gmx_bool                  bAllowMissing  = false;
    static gmx_bool                  addHydrogens   = false;

    //static const char               *ff[]           = {nullptr, "ACM-g", "ACM-pg", "ACM-s", "ACM-ps", "Verstraelen", nullptr};
    static const char               *qcustom        = nullptr;

    t_pargs                          pa[]     = 
    {
        { "-f",      FALSE, etSTR,  {&filename},
           "Input file name to be turned into GROMACS input" },
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Generate verbose output in the top file and on terminal." },
        { "-db",     FALSE, etSTR,  {&dbname},
          "Read a molecule from the database rather than from a file" },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-iupac",   FALSE, etSTR,  {&iupac},
          "IUPAC Name of your molecule" },
        { "-conf",  FALSE, etSTR, {&conf},
          "Conformation of the molecule" },
        { "-maxpot", FALSE, etINT, {&maxpot},
          "Fraction of potential points to read from the gaussian file (percent). If 100 all points are registered, else a selection of points evenly spread over the range of values is taken" },
        { "-allowmissing", FALSE, etBOOL, {&bAllowMissing},
          "Make a topology even if there are no force field parameters for all interactions" },
        { "-nsymm", FALSE, etINT, {&nsymm},
          "Symmetry number of the molecule can be supplied here if you know there is an error in the input file" },
        { "-pdbq",  FALSE, etBOOL, {&bUsePDBcharge},
          "HIDDENUse the B-factor supplied in a pdb file for the atomic charges" },
        { "-spacing", FALSE, etREAL, {&spacing},
          "Spacing of grid points (nm) for computing the potential (not used when a reference file is read)." },
        { "-border", FALSE, etREAL, {&border},
          "Spacing around the compound (nm) for computing the potential (not used when a reference file is read)." },
        { "-watoms", FALSE, etREAL, {&watoms},
          "Weight for the atoms when fitting the charges to the electrostatic potential. The potential on atoms is usually two orders of magnitude larger than on other points (and negative). For point charges or single smeared charges use 0. For point+smeared charges 1 is recommended." },
        //{ "-ff",     FALSE, etENUM, {ff},
        //  "Force field model. Note that only ACM-xx will yield complete topologies but see help text ([TT]-h[tt])." },
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Total charge of the molecule. This will be taken from the input file by default, but that is reliable only if the input is a Gaussian log file." },
        { "-qqm",    FALSE, etSTR,  {&qqm},
          "Use a method from quantum mechanics that needs to be present in the input file. Either ESP, Hirshfeld, CM5 or Mulliken may be available." },
        { "-addh",   FALSE, etBOOL, {&addHydrogens},
          "Add hydrogen atoms to the compound - useful for PDB files." },
        { "-qsymm",  FALSE, etBOOL, {&bQsym},
          "Symmetrize the charges on symmetric groups, e.g. CH3, NH2." },
        { "-symm",   FALSE, etSTR, {&symm_string},
          "Use the order given here for symmetrizing, e.g. when specifying [TT]-symm '0 1 0'[tt] for a water molecule (H-O-H) the hydrogens will have obtain the same charge. For simple groups, like methyl (or water) this is done automatically, but higher symmetry is not detected by the program. The numbers should correspond to atom numbers minus 1, and point to either the atom itself or to a previous atom." },
        { "-qcustom", FALSE, etSTR, {&qcustom}, 
          "Here a quoted string of custom charges can be provided such that a third party source can be used. It is then possible to generate multipoles and compare the ESP to a quantum chemistry result. The number of charges provided must match the number of particles (including shells if present in the force field used)." },
         { "-jobtype",  FALSE, etSTR, {&jobtype},
          "The job type used in the Gaussian calculation: Opt, Polar, SP, and etc." }
    };

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    Poldata        pd;
    t_inputrec    *inputrec = new t_inputrec();
    CommunicationRecord cr;
    gmx::MDLogger  mdlog {};
    std::string    method, basis;

    /* Check the options */
    bITP = opt2bSet("-oi", NFILE, fnm);
    // Check whether there is something to read
    if (strlen(dbname) == 0 && strlen(filename) == 0)
    {
        gmx_fatal(FARGS, "Specify either the -db or the -f option. No output without input");
    }
    const char *gentop_fnm = opt2fn_null("-ff", NFILE, fnm);
    //if (opt2parg_bSet("-ff", asize(pa), pa) && nullptr == gentop_fnm)
    //{
    //    gentop_fnm = gmx::formatString("%s.dat", ff[0]).c_str();
    // }
    //if (nullptr == gentop_fnm)
    //{
    //    fprintf(stderr, "Please specify either a force field file name or use the -ff flag");
    //    return 0;
    //}

    /* Read standard atom properties */
    aps = gmx_atomprop_init();
    if (!bVerbose)
    {
        gmx_atomprop_quiet(aps);
    }
    try
    {
        readPoldata(gentop_fnm, &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    (void) pd.verifyCheckSum(stderr);

    auto fs = pd.findForcesConst(InteractionType::CHARGEDISTRIBUTION);
    const char *ct = "chargetype";
    if (!fs.optionExists(ct))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The option %s is not defined in the force field file %s", ct, gentop_fnm).c_str()));
    }
    auto qType = fs.optionValue(ct);
    printf("Using force field file %s and charge distribution model %s\n",
           gentop_fnm, qType.c_str());
    if (bVerbose)
    {
        printf("Reading force field information. There are %d atomtypes.\n",
               static_cast<int>(pd.getNatypes()));
    }

    MyMol                mymol;
    if (strlen(dbname) > 0)
    {
        const char *molpropDatabase = opt2fn_null("-mp", NFILE, fnm);
        if (!molpropDatabase || strlen(molpropDatabase) == 0)
        {
            gmx_fatal(FARGS, "Empty database file name");
        }
        if (bVerbose)
        {
            printf("Looking up %s in molecule database %s.\n",
                   dbname, molpropDatabase);
        }
        MolPropIterator      mpi;
        std::vector<MolProp> mps;
        MolPropRead(molpropDatabase, &mps);
        for (mpi = mps.begin(); (mpi < mps.end()); mpi++)
        {
            if (strcasecmp(dbname, mpi->getMolname().c_str()) == 0)
            {
                break;
            }
        }
        if (mpi == mps.end())
        {
            gmx_fatal(FARGS, "Molecule %s not found in database", dbname);
        }
        mymol.Merge(&(*mpi));
    }
    else
    {
        if (strlen(molnm) == 0)
        {
            molnm = (char *)"MOL";
        }
        MolProp  mp;
        double qtot_babel = qtot;
        if (readBabel(filename,
                      &mp,
                      molnm,
                      iupac,
                      conf,
                      &method, 
                      &basis,
                      maxpot,
                      nsymm,
                      jobtype,
                      &qtot_babel,
                      addHydrogens))
        {
            std::map<std::string, std::string> g2a;
            gaffToAlexandria("", &g2a);
            bool mappingOK = true;
            if (!g2a.empty())
            {
                mappingOK = renameAtomTypes(&mp, g2a);
            }
            if (mappingOK)
            {
                mymol.Merge(&mp);
            }
        }
        else
        {
            gmx_fatal(FARGS, "No input file has been specified.");
        }
    }
    fill_inputrec(inputrec);
    mymol.setInputrec(inputrec);
    imm = mymol.GenerateTopology(stdout, &pd,
                                 bAllowMissing ? missingParameters::Ignore : missingParameters::Error);

    gmx_omp_nthreads_init(mdlog, cr.commrec(), 1, 1, 1, 0, false, false);
    if (immStatus::OK == imm)
    {
        mymol.symmetrizeCharges(&pd, bQsym, symm_string);
        maxpot = 100; // Use 100 percent of the ESP read from Gaussian file.
        
        mymol.initQgenResp(&pd, 0.0, maxpot);

        ChargeGenerationAlgorithm alg = ChargeGenerationAlgorithm::NONE;
        std::vector<double> myq;
        if (qcustom)
        {
            auto mycharges = gmx::splitString(qcustom);
            for(auto &q : mycharges)
            {
                myq.push_back(my_atof(q.c_str(), "custom q"));
            }
            alg = ChargeGenerationAlgorithm::Custom;
        }
        else if (strlen(qqm) > 0)
        {
            alg = nameToChargeGenerationAlgorithm(qqm);
        }
        imm    = mymol.GenerateCharges(&pd, mdlog, &cr, alg, myq);
    }
    /* Generate output file for debugging if requested */
    if (immStatus::OK == imm)
    {
        mymol.plotEspCorrelation(opt2fn_null("-plot_esp", NFILE, fnm), oenv);
    }

    if (immStatus::OK == imm)
    {
        mymol.GenerateCube(&pd,
                           spacing, border,
                           opt2fn_null("-ref",      NFILE, fnm),
                           opt2fn_null("-pc",       NFILE, fnm),
                           opt2fn_null("-pdbdiff",  NFILE, fnm),
                           opt2fn_null("-pot",      NFILE, fnm),
                           opt2fn_null("-rho",      NFILE, fnm),
                           opt2fn_null("-his",      NFILE, fnm),
                           opt2fn_null("-diff",     NFILE, fnm),
                           opt2fn_null("-diffhist", NFILE, fnm),
                           oenv);
    }

    if (immStatus::OK == imm && mymol.errors().size() == 0)
    {
        if (opt2bSet("-openmm", NFILE, fnm))
        {
            writeOpenMM(opt2fn("-openmm", NFILE, fnm), &pd, &mymol, 0);
        }
        else
        {
            mymol.PrintConformation(opt2fn("-c", NFILE, fnm));
            mymol.PrintTopology(bITP ? ftp2fn(efITP, NFILE, fnm) : ftp2fn(efTOP, NFILE, fnm),
                                bVerbose,
                                &pd,
                                &cr,
                                method,
                                basis,
                                bITP);
        }
    }
    else
    {
        auto fn = opt2fn("-g", NFILE, fnm);
        fprintf(stderr, "\nFatal Error. Please check the %s file for error messages.\n", fn);
        print_errors(fn, mymol.errors(), imm);
    }
    return 0;
}

} // namespace alexandria
