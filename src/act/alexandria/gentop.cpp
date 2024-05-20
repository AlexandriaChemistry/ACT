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

#include <cctype>
#include <cstdlib>

#include <algorithm>

#include "act/alexandria/actmol.h"
#include "act/alexandria/alex_modules.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/fetch_charges.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/openmm_xml.h"
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
    static const char               *desc[] = {
        "gentop generates a topology from molecular coordinates",
        "from a database or, optionally, from a coordinate file or",
        " a gaussian log file.",
        "The program assumes all hydrogens are present when defining",
        "the hybridization from the atom name and the number of bonds.[PAR]",
        "Recommend usage is to supply a molprop database",
        "with charges precalculate and then use the [TT]-db 'molecule(s)'[tt]",
        "option to extract compounds from the molprop database.[PAR]",
        "If the [TT]-oi[tt] option is set an [TT]itp[tt] file will be generated",
        "instead of a [TT]top[tt] file.",
        "A force field file should be specified using the [TT]-ff[tt] flag.",
        "Depending on the selected force field file, the force field",
        "uses polarizable or non-polarizable charges or virtual sites.[PAR]",
        "When the [TT]-openmm[tt] flag is passed, an XML file will be created",
        "that can be used to run a simulation of the system using the OpenMM",
        "software. In addition, a selection of parameters for the OpenMM simulation",
        "will be written to a [TT].dat[tt] file the name of which can be specified",
        "with the [TT]-openmm_sim[tt] flag."
    };
    gmx_output_env_t *oenv;
    gmx_atomprop_t    aps;
    immStatus         imm;

    t_filenm          fnm[] = {
        { efXML, "-ff",         "aff",       ffREAD  },
        { efTOP, "-p",          "out",       ffOPTWR },
        { efXML, "-openmm",     "out",       ffOPTWR },
        { efDAT, "-openmm_sim", "sim",       ffOPTWR },
        { efITP, "-oi",         "out",       ffOPTWR },
        { efPDB, "-c",          "out",       ffWRITE },
        { efNDX, "-n",          "renum",     ffOPTWR },
        { efDAT, "-q",          "qout",      ffOPTWR },
        { efXML, "-charges",    "molprops",  ffOPTRD },
        { efCUB, "-pot",        "potential", ffOPTWR },
        { efCUB, "-ref",        "refpot",    ffOPTRD },
        { efCUB, "-diff",       "diffpot",   ffOPTWR },
        { efCUB, "-rho",        "density",   ffOPTWR },
        { efXVG, "-diffhist",   "diffpot",   ffOPTWR },
        { efXVG, "-his",        "pot-histo", ffOPTWR },
        { efXVG, "-pc",         "pot-comp",  ffOPTWR },
        { efPDB, "-pdbdiff",    "pdbdiff",   ffOPTWR },
        { efLOG, "-g",          "errors",    ffWRITE }
    };

    const  int                       NFILE          = asize(fnm);

    static int                       maxpot         = 100;
    static int                       nsymm          = 0;
    static real                      qtot           = 0;
    static real                      watoms         = 0;
    static real                      spacing        = 0.01;
    static real                      border         = 0.2;
    static real                      mDrude         = 0.1;
    static char                     *molnm          = (char *)"";
    static char                     *iupac          = (char *)"";
    static char                     *dbname         = (char *)"";
    static char                     *symm_string    = (char *)"";
    static char                     *qqm            = (char *)"";
    static char                     *conf           = (char *)"minimum";
    static char                     *jobtype        = (char *)"Opt";
    static char                     *molFile        = (char *)"";
    //! Write shells to trajectory and coordinates
    bool                             writeShells    = false;
    static gmx_bool                  bITP           = false;
    static gmx_bool                  bUsePDBcharge  = false;
    static gmx_bool                  bVerbose       = false;
    static gmx_bool                  bAllowMissing  = false;
    static gmx_bool                  addHydrogens   = false;
    static gmx_bool                  addNumbersToAtoms = true;
    static gmx_bool                  genCharges      = false;
    static rvec                      mybox           = { 0, 0, 0 };
    //static const char               *ff[]           = {nullptr, "ACM-g", "ACM-pg", "ACM-s", "ACM-ps", "Verstraelen", nullptr};
    static const char               *qcustom        = nullptr;

    t_pargs                          pa[]     =
    {
        { "-f",      FALSE, etSTR,  {&molFile},
           "Input file name" },
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Generate verbose output in the top file and on terminal." },
        { "-db",     FALSE, etSTR,  {&dbname},
          "Read one or more molecules from the database rather than from a file. To specify multiple molecules please use quotes, e.g. [TT]-db[tt] 'water methane ammonia'." },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-iupac",   FALSE, etSTR,  {&iupac},
          "IUPAC Name of your molecule" },
        { "-conf",  FALSE, etSTR, {&conf},
          "Conformation of the molecule" },
        { "-generateCharges", FALSE, etBOOL, {&genCharges},
          "Generate charges for your compound using the Alexandria Charge Model." },
        { "-ws",     FALSE, etBOOL, {&writeShells},
          "Write coordinates of shell particles to trajectory and final coordinates as well." },
        { "-box",    FALSE, etRVEC, {mybox},
          "If set will replace the simulation box in the output coordinate file." },
        { "-mDrude", FALSE, etREAL, {&mDrude},
          "Mass to use for the drude particle if any. Default is 0.1 Da" },
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
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Total charge of the molecule. This will be taken from the input file by default, but that is reliable only if the input is a Gaussian log file." },
        { "-qqm",    FALSE, etSTR,  {&qqm},
          "Use a method from quantum mechanics that needs to be present in the input file. Either ESP, Hirshfeld, CM5 or Mulliken may be available." },
        { "-addh",   FALSE, etBOOL, {&addHydrogens},
          "Add hydrogen atoms to the compound - useful for PDB files." },
        { "-symm",   FALSE, etSTR, {&symm_string},
          "Use the order given here for symmetrizing, e.g. when specifying [TT]-symm '0 1 0'[tt] for a water molecule (H-O-H) the hydrogens will have obtain the same charge. For simple groups, like methyl (or water) this is done automatically, but higher symmetry is not detected by the program. The numbers should correspond to atom numbers minus 1, and point to either the atom itself or to a previous atom." },
        { "-qcustom", FALSE, etSTR, {&qcustom},
          "Here a quoted string of custom charges can be provided such that a third party source can be used. It is then possible to generate multipoles and compare the ESP to a quantum chemistry result. The number of charges provided must match the number of particles (including shells if present in the force field used)." },
        { "-numberAtypes", FALSE, etBOOL, {&addNumbersToAtoms},
          "Add a number index to OpenMM atomtypes when generating output for that software." },
         { "-jobtype",  FALSE, etSTR, {&jobtype},
          "The job type used in the Gaussian calculation: Opt, Polar, SP, and etc." }
    };
    int status = 0;
    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        status = 1;
        return status;
    }

    ForceField        pd;
    CommunicationRecord cr;
    gmx::MDLogger  mdlog {};
    std::string    method, basis;

    /* Check the options */
    bITP = opt2bSet("-oi", NFILE, fnm);

    const char *gentop_fnm = opt2fn_null("-ff", NFILE, fnm);
    if (nullptr == gentop_fnm)
    {
        fprintf(stderr, "Please pass me a force field file name with the -ff option.\n");
        status = 1;
        return status;
    }

    /* Read standard atom properties */
    aps = gmx_atomprop_init();
    if (!bVerbose)
    {
        gmx_atomprop_quiet(aps);
    }
    try
    {
        readForceField(gentop_fnm, &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    (void) pd.verifyCheckSum(stderr);

    auto &fs = pd.findForcesConst(InteractionType::COULOMB);
    std::string my_pol;
    if (pd.polarizable())
    {
        my_pol.assign(" polarizable");
    }
    auto qType = potentialToChargeType(fs.potential());
    printf("Using%s force field file %s and charge distribution model %s\n",
           my_pol.c_str(), gentop_fnm, chargeTypeName(qType).c_str());
    if (bVerbose)
    {
        printf("Read force field information. There are %zu atomtypes.\n",
               pd.getNatypes());
    }

    std::vector<ACTMol> actmols;
    matrix              box       = {{ 0 }};
    auto                forceComp = new ForceComputer();
    chargeMap           qmap;
    if (strlen(molFile) > 0)
    {
        if (strlen(molnm) == 0)
        {
            molnm = (char *)"MOL";
        }
        std::vector<MolProp> mps;
        double               qtot_babel = qtot;
        bool                 userqtot   = opt2parg_bSet("-qtot", asize(pa), pa);
        if (readBabel(&pd, molFile, &mps, molnm, iupac, conf, &method, &basis,
                      maxpot, nsymm, jobtype, userqtot, &qtot_babel, addHydrogens, box, false))
        {
            for(auto &mp : mps)
            {
                ACTMol mm;
                mm.Merge(&mp);
                actmols.push_back(mm);
            }
        }
        else
        {
            fprintf(stderr, "Could not read molecule(s) from  input file %s.\n", molFile);
            return 1;
        }
        if (nullptr != opt2fn_null("-pdbdiff",  NFILE, fnm))
        {
            fprintf(stderr, "WARNING: Cannot generate a pdbdiff file from user defined coordinates. Use the -db flag instead.\n");
        }
    }
    const char *molpropDatabase = opt2fn_null("-charges", NFILE, fnm);
    if (molpropDatabase && strlen(molpropDatabase) > 0)
    {
        std::vector<MolProp> mps;
        MolPropRead(molpropDatabase, &mps);
        if (strlen(molFile) == 0)
        {
            std::vector<std::string> mymols;
            if (strlen(dbname) > 0)
            {
                if (bVerbose)
                {
                    printf("Looking up %s in a molecule database %s.\n",
                           dbname, molpropDatabase);
                }
                mymols = split(dbname, ' ');
            }
            if (mymols.empty())
            {
                for(const auto &mp : mps)
                {
                    ACTMol mm;
                    mm.Merge(&(mp));
                    actmols.push_back(mm);
                }
            }
            else
            {
                // Throw away those compounds that are not in the selection
                for(auto mp = mps.begin(); mp < mps.end(); )
                {
                    const auto molnm = std::find_if(mymols.begin(), mymols.end(),
                                                    [mp](const std::string &x)
                                                    { return x == mp->getMolname(); });
                    if (molnm == mymols.end())
                    {
                        mp = mps.erase(mp);
                    }
                    else
                    {
                        ++mp;
                    }
                }
                for(auto mp : mps)
                {
                    ACTMol mm;
                    mm.Merge(&mp);
                    actmols.push_back(mm);
                }
            }
            if (actmols.empty())
            {
                fprintf(stderr, "Couldn't find any of the selected molecules\n");
                return 0;
            }
        }
        qmap = fetchChargeMap(&pd, forceComp, mps);
    }
    gmx_omp_nthreads_init(mdlog, cr.commrec(), 1, 1, 1, 0, false, false);
    int mp_index   = 1;
    std::map<std::string, std::pair<immStatus, std::vector<std::string>>> errors;
    for(auto actmol = actmols.begin(); actmol < actmols.end(); )
    {
        imm = actmol->GenerateTopology(stdout, &pd,
                                       bAllowMissing ? missingParameters::Ignore : missingParameters::Error);

        std::vector<gmx::RVec> forces(actmol->atomsConst().size());
        std::vector<gmx::RVec> coords = actmol->xOriginal();
        ForceComputer fcomp;
        fcomp.generateVsites(actmol->topology(), &coords);
        if (immStatus::OK == imm)
        {
            maxpot = 100; // Use 100 percent of the ESP read from QM file.
            std::map<MolPropObservable, iqmType> iqm = {
                { MolPropObservable::POTENTIAL, iqmType::QM },
                { MolPropObservable::CHARGE,    iqmType::QM }
            };
            actmol->getExpProps(&pd, iqm, 0.0, maxpot);
            auto fragments  = actmol->fragmentHandler();
            auto topologies = fragments->topologies();
            if (fragments->setCharges(qmap))
            {
                // Copy charges to the high-level topology as well
                fragments->fetchCharges(actmol->atoms());
            }
            else
            {
                auto qtype = qType::Calc;
                std::vector<double> myq;

                if (qcustom)
                {
                    // Second, if there are charges on the command line
                    fprintf(stderr, "WARNING: you provided charges on the command line. There is no guarantee that those charges can be used together with other parts of the force field and provide reasonable values..\n");
                    auto mycharges = gmx::splitString(qcustom);
                    for(auto &q : mycharges)
                    {
                        myq.push_back(my_atof(q.c_str(), "custom q"));
                    }
                    imm = actmol->GenerateCharges(&pd, forceComp, ChargeGenerationAlgorithm::Custom,
                                                  qtype, myq, &coords, &forces, true);
                }
                else if (genCharges)
                {
                    // Finally generate charges
                    auto alg   = pd.chargeGenerationAlgorithm();
                    fprintf(stderr, "WARNING: Using %s to generate charges. It is recommended to use a charge database instead of this option.\n", chargeGenerationAlgorithmName(alg).c_str());
                    imm = actmol->GenerateCharges(&pd, forceComp, alg, qtype, myq, &coords, &forces, true);
                }
                else
                {
                    fprintf(stderr, "Skipping %s since there are no charges available, please provide a charge database or use the -generateCharges flag.\n", actmol->getMolname().c_str());
                    imm = immStatus::ChargeGeneration;
                }
            }
        }
        if (immStatus::OK == imm)
        {
            actmol->GenerateCube(&pd, coords, forceComp,
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

        if (immStatus::OK == imm && actmol->errors().size() == 0)
        {
            std::string index;
            if (actmols.size() > 1)
            {
                index = gmx::formatString("%d_", mp_index);
            }
            if (!opt2bSet("-openmm", NFILE, fnm))
            {
                std::string tfn = gmx::formatString("%s%s", index.c_str(),
                                                    bITP ? ftp2fn(efITP, NFILE, fnm) : ftp2fn(efTOP, NFILE, fnm));
                actmol->PrintTopology(tfn.c_str(), bVerbose, &pd, forceComp,
                                     &cr, coords, method, basis, bITP);
            }
            if (opt2bSet("-c", NFILE, fnm))
            {
                std::string cfn = gmx::formatString("%s%s", index.c_str(),
                                                    opt2fn("-c", NFILE, fnm));
                if (opt2parg_bSet("-box", asize(pa), pa))
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        box[m][m] = mybox[m];
                    }
                }
                actmol->PrintConformation(cfn.c_str(), coords, writeShells, box);
            }
            ++actmol;
        }
        else
        {
            errors.insert({ actmol->getMolname(), { imm, actmol->errors() } });
            actmol = actmols.erase(actmol);
        }
        mp_index++;
    }
    if (!actmols.empty())
    {
        if (opt2bSet("-openmm", NFILE, fnm) || opt2bSet("-openmm_sim", NFILE, fnm))
        {
            writeOpenMM(opt2fn("-openmm", NFILE, fnm),
                        opt2fn("-openmm_sim", NFILE, fnm), &pd, actmols, mDrude, addNumbersToAtoms);
        }
    }
    if (!errors.empty())
    {
        auto fn = opt2fn("-g", NFILE, fnm);
        auto fp = gmx_ffopen(fn, "w");
        fprintf(fp, "Errors encountered during processing:\n");
        for(const auto &mess : errors)
        {
            fprintf(fp, "%s: %s\n", mess.first.c_str(), immsg(mess.second.first));
            for(const auto &m : mess.second.second)
            {
                fprintf(fp, "    %s\n", m.c_str());
            }
        }
        gmx_ffclose(fp);
        fprintf(stderr, "\nPlease check the %s file for error messages.\n", fn);
        status = 1;
    }
    return status;
}

} // namespace alexandria
