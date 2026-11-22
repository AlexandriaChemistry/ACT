/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/openmm_xml.h"
#include "act/basics/chargemodel.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forces/forcecomputer.h"
#include "act/import/compound_reader.h"
#include "act/molprop/molprop_xml.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_json.h"
#include "act/utility/jsontree.h"
#include "act/utility/output_environment.h"
#include "gentop.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

namespace
{

//! \return the name of a file for the functions that take nullptr for none
const char *fileOrNull(const std::string &fileName)
{
    return fileName.empty() ? nullptr : fileName.c_str();
}

} // namespace

void GentopOptions::check() const
{
    if (forceFieldFile.empty())
    {
        GMX_THROW(gmx::InvalidInputError("There is no force field file (-ff)"));
    }
    if (std::string(compounds.fileName()).empty() && std::string(compounds.smiles()).empty() &&
        std::string(compounds.dbName()).empty() && compounds.molselect().nMol() == 0)
    {
        GMX_THROW(gmx::InvalidInputError("There is no compound: give a structure file (-f), a SMILES string (-smiles), or compounds of the charge file (-db or -sel)"));
    }
    checkSettings();
}

void GentopOptions::checkSettings() const
{
    if (!std::isfinite(drudeMass) || drudeMass < 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -mDrude should be zero or more", drudeMass)));
    }
    if (minTraining < 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %d of -ntrain should be zero or more", minTraining)));
    }
    if (!std::isfinite(spacing) || spacing <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -spacing should be more than zero", spacing)));
    }
    if (!std::isfinite(border) || border < 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -border should be zero or more", border)));
    }
    if (box)
    {
        for (double side : *box)
        {
            if (!std::isfinite(side) || side < 0)
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -box should be zero or more", side)));
            }
        }
    }
}

int runGentop(MsgHandler *msghandler, const GentopOptions &options, gmx_output_env_t *oenv)
{
    options.check();
    ForceField pd;
    readForceField(options.forceFieldFile, &pd, msghandler);
    (void) pd.verifyCheckSum(msghandler);
    // GROMACS and OpenMM do not support second charges, use virtual sites instead
    pd.expandSecondCharges(msghandler);

    std::string my_pol;
    if (pd.polarizable())
    {
        my_pol.assign(" polarizable");
    }
    auto qType = pd.chargeDistribution();
    printf("Using%s force field file %s and charge distribution model %s\n",
           my_pol.c_str(), options.forceFieldFile.c_str(), chargeDistributionTypeName(qType).c_str());
    msghandler->write(gmx::formatString("Read force field information. There are %zu atomtypes.\n",
                                        pd.getNatypes()));

    ForceComputer       forceComputer;
    auto               *forceComp = &forceComputer;
    std::vector<ACTMol> actmols;
    // Reading changes the reader
    auto                compR = options.compounds;
    compR.read(msghandler, pd, forceComp, &actmols);
    // Compounds without a topology are left out
    if (actmols.empty())
    {
        msghandler->msg(ACTStatus::Error, "No compound with a topology, nothing to write.");
        return 1;
    }
    // Write information about force field
    msghandler->write("");
    for (const auto &ii : pd.info())
    {
        msghandler->write(ii);
    }

    // The same residue names in the OpenMM and the PDB file
    const auto residueNames = nameResidues(msghandler, &actmols, options.residueNames);
    if (!msghandler->ok())
    {
        return 1;
    }
    // A program without a command line has no output environment for the
    // histograms of the potential
    const OutputEnvironment output(oenv);
    int mp_index   = 1;
    // The compounds in json format, in the file of -json if it is given
    JsonTree    jmolecules = JsonTree::array("molecules");
    for(auto actmol = actmols.begin(); actmol < actmols.end(); ++actmol)
    {
        std::vector<gmx::RVec> forces(actmol->atomsConst().size());
        std::vector<gmx::RVec> coords = actmol->xOriginal();
        forceComp->generateVsites(actmol->topology(), &coords);
        actmol->updateQprops(&pd, forceComp, &forces);
        if (!options.jsonFile.empty())
        {
            JsonTree jmol("compound");
            jmol.addObject(moleculeToJson(&(*actmol)));
            JsonTree jcoords("coordinates");
            jcoords.addObject("unit", "nm");
            jcoords.addObject(coordinatesToJson("positions", coords));
            jmol.addObject(jcoords);
            jmolecules.addObject(jmol);
        }
        actmol->GenerateCube(msghandler, &pd, coords, forceComp,
                             options.spacing, options.border,
                             fileOrNull(options.referencePotential),
                             fileOrNull(options.potentialComparison),
                             fileOrNull(options.pdbDifference),
                             fileOrNull(options.potential),
                             fileOrNull(options.density),
                             fileOrNull(options.potentialHistogram),
                             fileOrNull(options.potentialDifference),
                             fileOrNull(options.differenceHistogram),
                             output.get());

        std::string index;
        if (actmols.size() > 1)
        {
            index = gmx::formatString("%d_", mp_index);
        }
        if (!options.topologyFile.empty())
        {
            const std::string tfn = index + options.topologyFile;
            actmol->PrintTopology(msghandler, tfn.c_str(), &pd, forceComp,
                                  coords, options.includeTopology);
        }
        for(const auto &c : actmol->generateCommercials(msghandler, &pd, forceComp, coords))
        {
            msghandler->write(c);
        }

        if (!options.structureFile.empty())
        {
            matrix box = { { 5, 0, 0 }, { 0, 5, 0 }, { 0, 0, 5 }};
            const std::string cfn = index + options.structureFile;
            if (options.box)
            {
                for(int m = 0; m < DIM; m++)
                {
                    box[m][m] = (*options.box)[m];
                }
            }
            actmol->PrintConformation(cfn.c_str(), coords, options.writeShells, box);
        }
        mp_index++;
    }
    if (!options.openmmFile.empty())
    {
        writeOpenMM(msghandler, options.openmmFile, &pd, actmols, options.drudeMass,
                    options.minTraining, options.numberTypes, residueNames);
    }
    if (!options.jsonFile.empty())
    {
        JsonTree jtree("gentop");
        jtree.addObject(jmolecules);
        jtree.write(options.jsonFile, true);
    }
    return msghandler->ok() ? 0 : 1;
}

int gentop(int argc, char *argv[])
{
    std::vector<const char *> desc = {
        "gentop generates a topology from molecular coordinates",
        "from a database or, optionally, from a coordinate file or",
        "a gaussian log file.",
        "The program assumes all hydrogens are present when defining",
        "the hybridization from the atom name and the number of bonds.[PAR]",
        "A GROMACS topology is written with [TT]-p[tt], or an include topology",
        "([TT]itp[tt]) with [TT]-oi[tt], also together with [TT]-openmm[tt].",
        "A force field file should be specified using the [TT]-ff[tt] flag.",
        "Depending on the selected force field file, the force field",
        "uses polarizable or non-polarizable charges or virtual sites.[PAR]",
        "When the [TT]-openmm[tt] flag is passed, an XML file will be created",
        "that can be used to run a simulation of the system using the OpenMM",
        "software.[PAR]",
        "The residue names in the [TT]-openmm[tt] file and the [TT]-c[tt] file",
        "are the same, as are the atom names, such that OpenMM can find the bonds",
        "in a PDB file without CONECT records, for instance one with many molecules",
        "made by another program, from the residue and atom names. The residue names",
        "can be given with [TT]-resnames[tt], otherwise they are made from the",
        "IUPAC names of the compounds.[PAR]",
        "With [TT]-json[tt], the compounds are written in json format to the",
        "given file, gentop.json by default: for each compound the particles",
        "(atoms, shells and virtual sites) with their types, masses and",
        "charges, the bonds and the coordinates of all particles."
    };
    gmx_output_env_t *oenv;

    std::vector<t_filenm> fnm = {
        { efXML, "-ff",         "aff",       ffREAD  },
        { efTOP, "-p",          "out",       ffOPTWR },
        { efXML, "-openmm",     "out",       ffOPTWR },
        { efITP, "-oi",         "out",       ffOPTWR },
        { efPDB, "-c",          "out",       ffWRITE },
        { efCUB, "-pot",        "potential", ffOPTWR },
        { efCUB, "-ref",        "refpot",    ffOPTRD },
        { efCUB, "-diff",       "diffpot",   ffOPTWR },
        { efCUB, "-rho",        "density",   ffOPTWR },
        { efXVG, "-diffhist",   "diffpot",   ffOPTWR },
        { efXVG, "-his",        "pot-histo", ffOPTWR },
        { efXVG, "-pc",         "pot-comp",  ffOPTWR },
        { efPDB, "-pdbdiff",    "pdbdiff",   ffOPTWR },
        { efJSON, "-json",      "gentop",    ffOPTWR }
    };

    GentopOptions                    options;
    int                              ntrain         = options.minTraining;
    real                             spacing        = options.spacing;
    real                             border         = options.border;
    real                             mDrude         = options.drudeMass;
    char                            *resnames       = (char *)"";
    //! Write shells to trajectory and coordinates
    bool                             writeShells    = options.writeShells;
    gmx_bool                         addNumbersToAtoms = options.numberTypes;
    rvec                             mybox           = { 0, 0, 0 };

    std::vector<t_pargs> pa =
    {
        { "-ws",     FALSE, etBOOL, {&writeShells},
          "Write coordinates of shell particles to trajectory and final coordinates as well." },
        { "-box",    FALSE, etRVEC, {mybox},
          "If set will replace the simulation box in the output coordinate file." },
        { "-mDrude", FALSE, etREAL, {&mDrude},
          "Mass to use for the drude/shell particles in OpenMM output, if there are any. If set to zero, the masses from the input force field will be retained." },
        { "-spacing", FALSE, etREAL, {&spacing},
          "Spacing of grid points (nm) for computing the potential (not used when a reference file is read)." },
        { "-border", FALSE, etREAL, {&border},
          "Spacing around the compound (nm) for computing the potential (not used when a reference file is read)." },
        { "-ntrain", FALSE, etINT, {&ntrain},
          "Minimum number of training data to store a parameter to OpenMM. Parameters that have Mutability::Fixed are counted as having at least one training data point. This is rather fragile, so beware of this when increasing the number. Zero means taht even untrained parameters will be exported to OpenMM, use with care!" },
        { "-numberAtypes", FALSE, etBOOL, {&addNumbersToAtoms},
          "Add a number index to OpenMM atomtypes when generating output for that software." },
        { "-resnames", FALSE, etSTR, {&resnames},
          "Residue names for the compounds, at most four letters or digits each, separated by spaces or commas. An entry [TT]methanol=MOH[tt] names the compound methanol, which can be given as its name or its InChi, an entry [TT]MOH[tt] names the next compound without a name, in the order of the input. Residues known to OpenMM, such as HOH, cannot be used. Compounds without a name get one made from their IUPAC name." }
    };
    MsgHandler          msghandler;
    options.compounds.addOptions(&pa, &fnm, &desc);
    msghandler.addOptions(&pa, &fnm, "gentop.log");
    if (!parse_common_args(&argc, argv, 0, fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 1;
    }
    CommunicationRecord cr;
    cr.init(cr.size());
    msghandler.optionsFinished(fnm, &cr);
    try
    {
        options.compounds.optionsFinished(&msghandler, fnm);
        if (!msghandler.ok())
        {
            return 1;
        }
        const char *forceField = opt2fn_null("-ff", fnm.size(), fnm.data());
        if (nullptr == forceField)
        {
            fprintf(stderr, "Please pass me a force field file name with the -ff option.\n");
            return 1;
        }
        options.forceFieldFile  = forceField;
        const auto setFile      = [&fnm](const char *option)
        {
            return opt2bSet(option, fnm.size(), fnm.data()) ? std::string(opt2fn(option, fnm.size(), fnm.data()))
                                                            : std::string();
        };
        // A topology only if it is asked for
        options.includeTopology = opt2bSet("-oi", fnm.size(), fnm.data());
        options.topologyFile    = setFile(options.includeTopology ? "-oi" : "-p");
        options.openmmFile    = setFile("-openmm");
        options.structureFile = setFile("-c");
        const auto optional   = [&fnm](const char *option)
        {
            const char *name = opt2fn_null(option, fnm.size(), fnm.data());
            return std::string(name ? name : "");
        };
        options.jsonFile            = optional("-json");
        options.referencePotential  = optional("-ref");
        options.potentialComparison = optional("-pc");
        options.pdbDifference       = optional("-pdbdiff");
        options.potential           = optional("-pot");
        options.density             = optional("-rho");
        options.potentialHistogram  = optional("-his");
        options.potentialDifference = optional("-diff");
        options.differenceHistogram = optional("-diffhist");
        if (opt2parg_bSet("-box", pa.size(), pa.data()))
        {
            options.box = std::array<double, 3>{ mybox[XX], mybox[YY], mybox[ZZ] };
        }
        options.writeShells  = writeShells;
        options.drudeMass    = mDrude;
        options.minTraining  = ntrain;
        options.numberTypes  = addNumbersToAtoms;
        options.residueNames = resnames;
        options.spacing      = spacing;
        options.border       = border;
        return runGentop(&msghandler, options, oenv);
    }
    catch (const gmx::InvalidInputError &ex)
    {
        fprintf(stderr, "%s\n", ex.what());
        return 1;
    }
}

} // namespace alexandria
