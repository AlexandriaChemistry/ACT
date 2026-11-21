/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2026
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
#include <set>
#include <string>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/basics/act_particle.h"
#include "act/basics/atomprops.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/basics/libraryfile.h"
#include "act/basics/mutability.h"
#include "act/forcefield/act_checksum.h"
#include "act/forcefield/combruleutil.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_json.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield_validation.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/particletype.h"
#include "act/forcefield/potential.h"
#include "act/forces/forcecomputer.h"
#include "act/utility/jsontree.h"
#include "act/utility/stringutil.h"
#include "gen_ff.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

namespace alexandria
{

//! \brief Add symmetric charges to force field
static void add_symm_charges(ForceField *pd)
{
    bool useCWD = false;
    const auto filename = findLibrary("symmetric_charges.csv", useCWD);
    if (filename.empty())
    {
        return;
    }
    gmx::TextReader tr(filename);
    std::string     tmp;
    while (tr.readLine(&tmp))
    {
        auto ptr = split(tmp, ',');
        if (3 == ptr.size())
        {
            pd->addSymcharges(ptr[0], ptr[1],
                              my_atoi(ptr[2], "numattch"));

        }
    }
}

//! \brief Read atomtypes from a filename
static std::map<std::string, std::map<std::string, std::string>> read_atomtypes(const char *filename)
{
    std::map<std::string, std::map<std::string, std::string> > table;
    gmx::TextReader          tr(filename);
    std::string              tmp;
    std::vector<std::string> value;
    int                      lineno = 1;
    while (tr.readLine(&tmp))
    {
        auto ptr = split(tmp, '|');
        // There should be at least 30 columns, but could be more.
        if (30 <= ptr.size())
        {
            if (value.empty())
            {
                // Read first line
                for(size_t i = 1; i < ptr.size(); i++)
                {
                    value.push_back(ptr[i]);
                }
            }
            else if (lineno > 1 && ptr[0].find("#") == std::string::npos)
            {
                // Read other lines
                std::map<std::string, std::string> row;
                for(size_t i = 1; i < ptr.size(); i++)
                {
                    row[value[i-1]] = ptr[i];
                }
                table[ptr[0]] = row;
            }
        }
        else
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Found %zu items on line %d in %s\n",
                                                               ptr.size(), lineno, filename).c_str()));
        }
        lineno++;
    }

    return table;
}

/*! \brief The number of atoms that a virtual site of a type is constructed from
 * \param[in] itype The interaction type
 * \return the number, 0 if not a virtual site
 */
static size_t vsiteAtoms(InteractionType itype)
{
    switch (itype)
    {
    case InteractionType::VSITE1:
        return 1;
    case InteractionType::VSITE2:
    case InteractionType::VSITE2FD:
        return 2;
    case InteractionType::VSITE3:
    case InteractionType::VSITE3S:
    case InteractionType::VSITE3FD:
    case InteractionType::VSITE3FAD:
    case InteractionType::VSITE3OUT:
    case InteractionType::VSITE3OUTS:
        return 3;
    case InteractionType::VSITE4:
    case InteractionType::VSITE4S:
    case InteractionType::VSITE4S3:
        return 4;
    default:
        return 0;
    }
}

/*! \brief Add the virtual sites of JSON files to a force field
 *
 * Each file is an object with a list "virtual_sites", and optionally a
 * "description". Each virtual site has its particle "type", the
 * "interaction" type, the "atoms" it is constructed from and, except for
 * VSITE1, the "min" and "max" of each of the "parameters" of its potential,
 * by name. A parameter starts in the middle of its range.
 * \param[in]    vsfiles The files
 * \param[inout] pd      The force field
 * \throws gmx::InvalidInputError for an invalid file
 */
static void add_vsites(const std::vector<std::string> &vsfiles,
                       ForceField                     *pd)
{
    std::map<InteractionType, ForceFieldParameterList> i2f;
    // How the atoms in the identifiers of the virtual sites can be swapped
    std::map<InteractionType, CanSwap> vsiteSwap = {
        { InteractionType::VSITE1,     CanSwap::Yes },
        { InteractionType::VSITE2,     CanSwap::Vsite2 },
        { InteractionType::VSITE2FD,   CanSwap::Vsite2 },
        { InteractionType::VSITE3,     CanSwap::Vsite3 },
        { InteractionType::VSITE3S,    CanSwap::Vsite3 },
        { InteractionType::VSITE3FD,   CanSwap::Vsite3 },
        { InteractionType::VSITE3FAD,  CanSwap::No },
        { InteractionType::VSITE3OUT,  CanSwap::No },
        { InteractionType::VSITE3OUTS, CanSwap::No },
        { InteractionType::VSITE4,     CanSwap::No },
        { InteractionType::VSITE4S,    CanSwap::No },
        { InteractionType::VSITE4S3,   CanSwap::No }
    };
    // The potential of each type of virtual site
    std::map<InteractionType, Potential> vsitePotential = {
        { InteractionType::VSITE1,     Potential::VSITE1 },
        { InteractionType::VSITE2,     Potential::VSITE2 },
        { InteractionType::VSITE2FD,   Potential::VSITE2FD },
        { InteractionType::VSITE3,     Potential::VSITE3 },
        { InteractionType::VSITE3S,    Potential::VSITE3S },
        { InteractionType::VSITE3FD,   Potential::VSITE3FD },
        { InteractionType::VSITE3FAD,  Potential::VSITE3FAD },
        { InteractionType::VSITE3OUT,  Potential::VSITE3OUT },
        { InteractionType::VSITE3OUTS, Potential::VSITE3OUTS },
        { InteractionType::VSITE4,     Potential::VSITE4 },
        { InteractionType::VSITE4S,    Potential::VSITE4S },
        { InteractionType::VSITE4S3,   Potential::VSITE4S3 }
    };
    // How the parameter lists can be swapped, as for the identifiers except VSITE3FD
    std::map<InteractionType, CanSwap> listSwap = vsiteSwap;
    listSwap[InteractionType::VSITE3FD] = CanSwap::No;
    for (const auto &vsfile : vsfiles)
    {
        const auto error = [&vsfile](size_t index, const std::string &what)
              {
                  GMX_THROW(gmx::InvalidInputError(gmx::formatString("Virtual site %zu in file %s: %s",
                                                                     index + 1, vsfile.c_str(), what.c_str()).c_str()));
              };
        nlohmann::json js;
        {
            std::ifstream in(vsfile);
            if (!in)
            {
                GMX_THROW(gmx::FileIOError(gmx::formatString("Cannot open file %s", vsfile.c_str()).c_str()));
            }
            try
            {
                in >> js;
            }
            catch (const nlohmann::json::exception &ex)
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("File %s is not valid JSON: %s",
                                                                   vsfile.c_str(), ex.what()).c_str()));
            }
        }
        if (!js.is_object() || !js.contains("virtual_sites") || !js["virtual_sites"].is_array())
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("File %s has no list \"virtual_sites\"",
                                                               vsfile.c_str()).c_str()));
        }
        const auto &sites = js["virtual_sites"];
        for (size_t k = 0; k < sites.size(); k++)
        {
            const auto &site = sites[k];
            for (const char *key : { "type", "interaction", "atoms" })
            {
                if (!site.is_object() || !site.contains(key))
                {
                    error(k, gmx::formatString("\"%s\" is missing", key));
                }
            }
            if (!site["type"].is_string() || !site["interaction"].is_string() || !site["atoms"].is_array())
            {
                error(k, "\"type\" and \"interaction\" have to be strings and \"atoms\" a list");
            }
            const auto type  = site["type"].get<std::string>();
            const auto inter = site["interaction"].get<std::string>();
            const Identifier vsid(type);
            if (!pd->hasParticleType(vsid))
            {
                error(k, gmx::formatString("unknown particle type %s", type.c_str()));
            }
            InteractionType itype;
            if (!stringToInteractionType(inter, &itype) || vsiteSwap.find(itype) == vsiteSwap.end())
            {
                error(k, gmx::formatString("%s is not an interaction type of virtual sites", inter.c_str()));
            }
            std::string atoms;
            for (const auto &a : site["atoms"])
            {
                if (!a.is_string())
                {
                    error(k, "the atoms have to be strings");
                }
                atoms += (atoms.empty() ? "" : "~") + a.get<std::string>();
            }
            if (site["atoms"].size() != vsiteAtoms(itype))
            {
                error(k, gmx::formatString("%s has %zu atoms, not %zu", inter.c_str(),
                                           vsiteAtoms(itype), site["atoms"].size()));
            }
            if (i2f.find(itype) == i2f.end())
            {
                i2f.emplace(itype, ForceFieldParameterList(potentialToString(vsitePotential[itype]), listSwap[itype]));
            }
            auto       &fs    = i2f.find(itype)->second;
            const Identifier  vs(itype, atoms + "!" + type, vsiteSwap[itype]);
            auto        names = potentialToParameterName(fs.potential());
            if (InteractionType::VSITE1 == itype)
            {
                // Not a lot of information here, but we need something for downstream processing.
                fs.addParameter(vs, names[vsite1A],
                                ForceFieldParameter("", 1, 0, 0, 1, 1, Mutability::Fixed, false, false));
                continue;
            }
            if (!site.contains("parameters") || !site["parameters"].is_object())
            {
                error(k, "\"parameters\" is missing");
            }
            const auto &params = site["parameters"];
            for (const auto &p : params.items())
            {
                if (std::find(names.begin(), names.end(), p.key()) == names.end())
                {
                    error(k, gmx::formatString("%s has no parameter %s", inter.c_str(), p.key().c_str()));
                }
            }
            for (const std::string name : names)
            {
                if (!params.contains(name) || !params[name].is_object() ||
                    !params[name].contains("min") || !params[name].contains("max") ||
                    !params[name]["min"].is_number() || !params[name]["max"].is_number())
                {
                    error(k, gmx::formatString("parameter %s needs a \"min\" and a \"max\"", name.c_str()));
                }
                double pmin = params[name]["min"].get<double>();
                double pmax = params[name]["max"].get<double>();
                if (pmin > pmax)
                {
                    error(k, gmx::formatString("the minimum of %s is larger than its maximum", name.c_str()));
                }
                fs.addParameter(vs, name,
                                ForceFieldParameter("", (pmin+pmax)/2, 0, 0, pmin, pmax,
                                                    Mutability::Bounded, false, false));
            }
        }
    }
    for(const auto &iii: i2f)
    {
        if (!iii.second.empty())
        {
            pd->addForces(iii.first, iii.second);
        }
    }
}

/*! \brief Tool to generate a force field from scratch
 * \param[in] argc Number of arguments on the cmd line
 * \param[in] argv The actual arguments
 * \return 0 if ok, 1 otherwise
 */
namespace
{

/*! \brief The potentials of an option of gen_ff
 * \param[in] potentials The potentials
 * \return Their names, as the parser of an option with a choice wants them
 */
std::vector<const char *> choices(const std::vector<Potential> &potentials)
{
    std::vector<const char *> names = { nullptr };
    for (const auto &p : potentials)
    {
        names.push_back(potentialToString(p).c_str());
    }
    names.push_back(nullptr);
    return names;
}

/*! \brief Check that a potential is one of those an option offers
 * \param[in] p       The potential
 * \param[in] offered The potentials of the option
 * \param[in] option  The option
 * \throws gmx::InvalidInputError if it is not
 */
void checkPotential(Potential p, const std::vector<Potential> &offered, const char *option)
{
    if (std::find(offered.begin(), offered.end(), p) == offered.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The potential %s cannot be used with %s",
                                                           potentialToString(p).c_str(), option)));
    }
}

/*! \brief The potential of an option of gen_ff from its name
 * \param[in] name   The name the parser gave
 * \param[in] option The option
 * \return The potential
 * \throws gmx::InvalidInputError for an unknown name
 */
Potential potentialOf(const char *name, const char *option)
{
    Potential p;
    if (nullptr == name || !stringToPotential(name, &p))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown potential '%s' of %s", name ? name : "", option)));
    }
    return p;
}

//! The charge distributions as gen_ff -qdist and the Coulomb potentials name them
const std::map<ChargeDistributionType, std::string> &distributionNames()
{
    static const std::map<ChargeDistributionType, std::string> names = {
        { ChargeDistributionType::Gaussian, "GAUSSIAN" },
        { ChargeDistributionType::Point, "POINT" },
        { ChargeDistributionType::Slater, "SLATER" }
    };
    return names;
}

//! \return the name of a charge distribution as gen_ff -qdist has it
std::string distributionName(ChargeDistributionType distribution)
{
    return distributionNames().at(distribution);
}

/*! \brief The charge distribution of gen_ff -qdist from its name
 * \param[in] name The name the parser gave
 * \return The distribution
 * \throws gmx::InvalidInputError for an unknown name
 */
ChargeDistributionType distributionOf(const char *name)
{
    for (const auto &d : distributionNames())
    {
        if (name && d.second == name)
        {
            return d.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown charge distribution '%s' of -qdist", name ? name : "")));
}

} // namespace

const std::vector<Potential> &genFFBondPotentials()
{
    static const std::vector<Potential> potentials = { Potential::HARMONIC_BONDS, Potential::MORSE_BONDS,
                                                       Potential::CUBIC_BONDS, Potential::HUA_BONDS };
    return potentials;
}

const std::vector<Potential> &genFFAnglePotentials()
{
    static const std::vector<Potential> potentials = { Potential::HARMONIC_ANGLES, Potential::UREY_BRADLEY_ANGLES };
    return potentials;
}

const std::vector<Potential> &genFFDihedralPotentials()
{
    static const std::vector<Potential> potentials = { Potential::FOURIER_DIHEDRALS, Potential::PROPER_DIHEDRALS };
    return potentials;
}

const std::vector<Potential> &genFFVdwPotentials()
{
    static const std::vector<Potential> potentials = {
        Potential::LJ12_6, Potential::LJ14_7, Potential::LJ8_6, Potential::LJ12_6_4,
        Potential::GENERALIZED_BUCKINGHAM, Potential::WANG_BUCKINGHAM,
        Potential::BUCKINGHAM, Potential::TANG_TOENNIES, Potential::TT2b,
        Potential::SLATER_ISA_TT };
    return potentials;
}

const std::vector<Potential> &genFFInductionCorrectionPotentials()
{
    static const std::vector<Potential> potentials = { Potential::MACDANIEL_SCHMIDT, Potential::MORSE_BONDS };
    return potentials;
}

void GenFFOptions::check() const
{
    if (atomTypesFile.empty())
    {
        GMX_THROW(gmx::InvalidInputError("There is no table of particle types (-f)"));
    }
    if (outputFile.empty())
    {
        GMX_THROW(gmx::InvalidInputError("There is no file for the force field (-o)"));
    }
    if (nexcl < 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %d of -nexcl should be zero or more", nexcl)));
    }
    if (!std::isfinite(epsilonr) || epsilonr <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -epsilonr should be more than zero", epsilonr)));
    }
    if (!std::isfinite(shellMaxDistance) || shellMaxDistance < 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of -shellmaxdistance should be zero or more",
                                                           shellMaxDistance)));
    }
    checkPotential(bondPotential, genFFBondPotentials(), "-bondfn");
    checkPotential(anglePotential, genFFAnglePotentials(), "-anglefn");
    checkPotential(dihedralPotential, genFFDihedralPotentials(), "-dihfn");
    checkPotential(vdwPotential, genFFVdwPotentials(), "-vdwfn");
    checkPotential(inductionCorrection, genFFInductionCorrectionPotentials(), "-icfn");
}

void runGenFF(const GenFFOptions &options)
{
    options.check();
    auto table = read_atomtypes(options.atomTypesFile.c_str());
    printf("There are %zu atom types in the force field with %zu properties.\n",
           table.size(), table.begin()->second.size());

    auto aprops = readAtomProps();
    printf("There are %zu element properties\n", aprops.size());

    alexandria::ForceField pd;
    // More stuff
    const std::vector<std::string> typeOptions = {
              "acmtype", "bondtype", "element", "poltype", "row", "zetatype", rowVSName
          };

    // Check for Point charges
    const bool bPoint = ChargeDistributionType::Point == options.chargeDistribution;
    if (bPoint)
    {
        printf("You selected point charges. Oh dear...\n");
    }
    
    std::map<InteractionType, ForceFieldParameterList> ffpl;
    for(const auto &entry : table)
    {
        // Generate particle type
        auto apType = ActParticle::Atom;
        const std::string elem = table[entry.first]["element"];
        if (elem == "X")
        {
            apType = ActParticle::Shell;
        }
        else if (elem == "VS")
        {
            apType = ActParticle::Vsite;
        }
        auto ptp = ParticleType(Identifier(entry.first),
                                table[entry.first]["comment"], apType);

        // Now add the "options"
        for(const auto &opt : typeOptions)
        {
            if (!table[entry.first][opt].empty())
            {
                ptp.setOption(opt, table[entry.first][opt]);
            }
        }
        const auto apropsptr = aprops.find(elem);
        if (aprops.end() == apropsptr)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown element '%s' of atom type %s", elem.c_str(), entry.first.c_str()).c_str()));
        }
        ptp.setOption("atomnumber", gmx::formatString("%d", apropsptr->second.atomnumber()));
        ptp.setOption("vdwtype", entry.first);
        ptp.setOption("vdwcorrtype", entry.first);
        ptp.setOption("induccorrtype", entry.first);
        // Now "parameters"
        auto mass       = apropsptr->second.mass();
        ptp.addForceFieldParameter("mass",
                                   ForceFieldParameter("Da", mass, 0, 1, mass,
                                                       mass, Mutability::Fixed,
                                                       true, true));
        double vmin = my_atof(table[entry.first][":q:min"], ":q:min");
        double vmax = my_atof(table[entry.first][":q:max"], ":q:max");
        const std::string vunit = table[entry.first][":q:unit"];
        Mutability  vmut;
        const std::string qmut(":q:mut");
        if (!nameToMutability(table[entry.first][qmut], &vmut))
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Invalid Mutability '%s' for charge in atomtype %s", table[entry.first][qmut].c_str(), entry.first.c_str()).c_str()));
        }
        ptp.addForceFieldParameter("charge",
                                   ForceFieldParameter(vunit, (vmin+vmax)/2,
                                                       0, 0, vmin,
                                                       vmax, vmut,
                                                       true, false));
        // A second charge on top of the particle, with its distribution,
        // see ParticleType::chargeVS. As for a virtual site, the charge
        // starts in the middle of its range and the width at its minimum.
        if (!table[entry.first][std::string(":") + chargeVSName + ":min"].empty())
        {
            for(const char *pname : { chargeVSName, zetaVSName })
            {
                const std::string prefix = std::string(":") + pname + ":";
                auto       &row    = table[entry.first];
                double      pmin   = my_atof(row[prefix+"min"], prefix+"min");
                double      pmax   = my_atof(row[prefix+"max"], prefix+"max");
                Mutability  pmut;
                if (!nameToMutability(row[prefix+"mut"], &pmut))
                {
                    GMX_THROW(gmx::InvalidInputError(gmx::formatString("Invalid Mutability '%s' for %s in atomtype %s", row[prefix+"mut"].c_str(), pname, entry.first.c_str()).c_str()));
                }
                bool   charge = (std::string(pname) == chargeVSName);
                double value  = charge ? (pmin+pmax)/2 : pmin;
                ptp.addForceFieldParameter(pname,
                                           ForceFieldParameter(row[prefix+"unit"], value, 0, 0,
                                                               pmin, pmax, pmut, true, !charge));
            }
        }

        pd.addParticleType(ptp);

        auto myatype = table[entry.first];
        std::map<InteractionType, std::map<std::string, ForceFieldParameter>> atomparams;
        std::map<InteractionType, std::map<std::string, bool>>                nonNegative;
        for (const auto &key : myatype)
        {
            auto content = split(key.first, ':');
            if (content.size() != 3)
            {
                continue;
            }
            InteractionType iType;
            if (!stringToInteractionType(content[0], &iType))
            {
                if (!content[0].empty())
                {
                    fprintf(stderr, "Unknown interaction type '%s' ignored\n",
                            content[0].c_str());
                }
                continue;
            }
            const auto value = myatype[key.first];
            if (atomparams.find(iType) == atomparams.end())
            {
                atomparams[iType] = {};
            }
            if (nonNegative.find(iType) == nonNegative.end())
            {
                nonNegative[iType] = {};
            }
            const auto& pName = content[1];
            if (atomparams[iType].find(pName) == atomparams[iType].end())
            {
                atomparams[iType][pName] = {};
            }
            if (content[2] == "unit")
            {
                atomparams[iType][pName].setUnit(value);
            }
            else if (content[2] == "nonneg")
            {
                nonNegative[iType][pName] = value == "True";
            }
            else if (content[2] == "mut")
            {
                Mutability mut;
                if (nameToMutability(value, &mut))
                {
                    atomparams[iType][pName].setMutability(mut);
                }
                else if (!value.empty())
                {
                    fprintf(stderr, "Incomprehensible Mutability '%s' for atype '%s' parameter %s\n", value.c_str(), entry.first.c_str(), pName.c_str());
                }
            }
            else
            {
                double v = my_atof(value, content[2]);

                if (content[2] == "min")
                {
                    atomparams[iType][pName].setMinimum(v);
                }
                else if (content[2] == "max")
                {
                    atomparams[iType][pName].setMaximum(v);
                }
            }
        }
        for (auto ap=atomparams.begin(); ap != atomparams.end(); ++ap)
        {
            std::string key = entry.first;
            if (ap->first == InteractionType::ELECTRONEGATIVITYEQUALIZATION)
            {
                key = table[entry.first]["acmtype"];
            }
            else if (ap->first == InteractionType::ELECTROSTATICS)
            {
                key = table[entry.first]["zetatype"];
            }
            if (!key.empty())
            {
                for(const auto &pname : ap->second)
                {
                    // The parameter is made anew from what was read, such
                    // that it starts at its minimum whatever its mutability
                    // and is not marked as corrected, see
                    // ForceFieldParameter::wasReset.
                    const bool nonneg = nonNegative[ap->first][pname.first];
                    double     pmin   = pname.second.minimum();
                    double     pmax   = pname.second.maximum();
                    Mutability pmut   = pname.second.mutability();
                    if (pmin > pmax)
                    {
                        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The minimum %g of %s %s is larger than the maximum %g for atomtype %s",
                                                                           pmin, interactionTypeToString(ap->first).c_str(),
                                                                           pname.first.c_str(), pmax, entry.first.c_str()).c_str()));
                    }
                    if (nonneg)
                    {
                        pmin = std::max(0.0, pmin);
                        pmax = std::max(pmin, pmax);
                    }
                    if (bPoint && ap->first == InteractionType::ELECTROSTATICS && pname.first == "zeta")
                    {
                        // Point charges have no width
                        pmin = 0;
                        pmax = 0;
                        pmut = Mutability::Fixed;
                    }
                    ffpl[ap->first].addParameter(key, pname.first,
                                                 ForceFieldParameter(pname.second.unit(), pmin, 0, 0,
                                                                     pmin, pmax, pmut, true, nonneg));
                }
            }
        }
    }
    // Now add the necessary terms to the final force field.
    std::set<InteractionType> nonbond = { InteractionType::ELECTROSTATICS,
        InteractionType::VDW, InteractionType::VDWCORRECTION,
        InteractionType::INDUCTIONCORRECTION };
    for(auto &fp : ffpl)
    {
        if (nonbond.find(fp.first) != nonbond.end())
        {
            fp.second.addOption("nexcl", gmx_itoa(options.nexcl));
            fp.second.setCanSwap(CanSwap::Yes);
        }
        if (fp.first == InteractionType::ELECTROSTATICS)
        {
            fp.second.addOption("epsilonr", gmx_ftoa(options.epsilonr));
            // COULOMB_POINT, COULOMB_GAUSSIAN or COULOMB_SLATER
            fp.second.setFunction("COULOMB_" + distributionName(options.chargeDistribution));
        }
        else if (fp.first == InteractionType::VDW)
        {
            fp.second.setFunction(potentialToString(options.vdwPotential));
        } 
        else if (fp.first == InteractionType::VDWCORRECTION && options.vdwCorrection)
        {
            fp.second.setFunction(potentialToString(Potential::BORN_MAYER));
        } 
        else if (fp.first == InteractionType::INDUCTIONCORRECTION)
        {
            fp.second.setFunction(potentialToString(options.inductionCorrection));
        }
        else if (fp.first == InteractionType::ELECTRONEGATIVITYEQUALIZATION)
        {
            fp.second.setCanSwap(CanSwap::Yes);
        }
        pd.addForces(fp.first, fp.second);
    }
    // Polarization function
    auto itPol = InteractionType::POLARIZATION;
    if (pd.interactionPresent(itPol))
    {
        auto fsp = pd.findForces(itPol);
        fsp->setFunction(potentialToString(Potential::POLARIZATION));
        fsp->setCanSwap(CanSwap::Yes);
        fsp->addOption(shellMaxDistanceName, gmx_ftoa(options.shellMaxDistance));
    }
    // Bonded functions need to be added
    const ForceFieldParameterList bonds(potentialToString(options.bondPotential), CanSwap::Yes);
    pd.addForces(InteractionType::BONDS, bonds);
    const ForceFieldParameterList angles(potentialToString(options.anglePotential), CanSwap::Yes);
    pd.addForces(InteractionType::ANGLES, angles);
    const ForceFieldParameterList linang("LINEAR_ANGLES", CanSwap::Linear);
    pd.addForces(InteractionType::LINEAR_ANGLES, linang);
    const ForceFieldParameterList idihs("IDIHS", CanSwap::Idih);
    pd.addForces(InteractionType::IMPROPER_DIHEDRALS, idihs);
    const ForceFieldParameterList pdihs(potentialToString(options.dihedralPotential), CanSwap::Yes);
    pd.addForces(InteractionType::PROPER_DIHEDRALS, pdihs);
    
    // Combination rules after adding the forces and parameters only, since it will look them up.
    CombRuleUtil crule;
    crule.setRules(options.combinationRules);
    crule.extract(&pd);
    // Otherwise the force field cannot be used
    crule.addDefaults(&pd);
    // Virtual sites
    if (!options.virtualSiteFiles.empty())
    {
        add_vsites(options.virtualSiteFiles, &pd);
    }
    // The charge distributions and the options of COULOMB are those of the
    // particle types, and virtual sites on top of atoms that only carry a
    // charge are stored as second charges, as when reading a force field.
    pd.convertCoulomb(nullptr);
    pd.convertVsite1(nullptr);
    if (options.symmetricCharges)
    {
        add_symm_charges(&pd);
    }
    // As when a force field is read, for the description in JSON
    pd.checkForPolarizability();
    pd.updateTimeStamp();
    pd.updateCheckSum();
    writeForceField(options.outputFile, &pd, false);
    if (!options.jsonFile.empty())
    {
        // The force field as it was written, for programs that show it
        JsonTree jtree("gen_ff");
        jtree.addObject(forceFieldToJson(pd));
        jtree.addObject(forceFieldIssuesToJson(validateForceField(pd)));
        jtree.write(options.jsonFile, true);
    }

}

int gen_ff(int argc, char*argv[])
{
    std::vector<const char *> desc =
    {
        "gen_ff generates a new force field file based on a specification",
        "provided by the user in one or more csv files and command line options.[PAR]",
        "Most important is the atomtypes file which is generated by the ACT GUI.",
        "In addition one or more JSON files with virtual sites can be provided",
        "with the [TT]-vs[tt] option. Default is not to add virtual sites.",
        "Each file has a list [TT]virtual_sites[tt] and optionally a [TT]description[tt].",
        "Each virtual site has its particle [TT]type[tt], the [TT]interaction[tt]",
        "type, such as VSITE2, the [TT]atoms[tt] it is constructed from, and, except",
        "for VSITE1, the [TT]min[tt] and [TT]max[tt] of each of the [TT]parameters[tt]",
        "of its potential, by name, for instance[BR]",
        "[TT]{ \"type\": \"v2f1\", \"interaction\": \"VSITE2\", \"atoms\": [\"f_b\", \"hf_b\"],",
        "\"parameters\": { \"vs2a\": { \"min\": -1, \"max\": 0.5 } } }[tt][BR]",
        "A parameter starts in the middle of its range. The data directory of ACT has",
        "examples, such as vsite2_halogen_sigma_holes.json.[PAR]",
        "A distributed charge on top of an atom, for instance in point+Gaussian",
        "models, can be specified as a second charge of the atom type with the",
        "columns [TT]:chargeVS:[tt], [TT]:zetaVS:[tt] and [TT]rowVS[tt]. Virtual",
        "sites of type VSITE1 that only carry a charge are replaced by such second",
        "charges as well.[PAR]"
        "Some input files are taken from the data directory of ACT, or from",
        "the directory in the ACTDATA environment variable if that is set, namely:[BR]",
        "atomprops.csv - containing properties for atoms.[BR]",
        "symmetric_charges.csv - containing groups of atoms that should have symmetric charges (optional).[PAR]",
        "With [TT]-json[tt] the force field that was made is written in json format",
        "to the given file, gen_ff.json by default, in the same form as",
        "edit_ff writes it, together with the result of the checks that",
        "edit_ff -validate does. This is meant for programs that show the",
        "force field, such as the user interface of ACT."
    };
    CombRuleUtil crule;
    crule.addInfo(&desc);
    gmx_output_env_t *oenv;
    GenFFOptions      options;
    int               nexcl    = options.nexcl;
    int               nexclqq  = -1;
    int               nexclvdw = -1;
    real              epsilonr = options.epsilonr;
    real              shellMaxDistance = options.shellMaxDistance;
    bool              qsymm    = options.symmetricCharges;
    std::vector<const char *>qdn2    = { nullptr, "GAUSSIAN", "POINT", "SLATER", nullptr };
    auto bondfn  = choices(genFFBondPotentials());
    auto anglefn = choices(genFFAnglePotentials());
    auto dihfn   = choices(genFFDihedralPotentials());
    auto vdwfn   = choices(genFFVdwPotentials());
    auto icfn    = choices(genFFInductionCorrectionPotentials());

    std::vector<t_filenm> fnm = {
        { efCSV, "-f",   "atomtypes", ffREAD  },
        { efJSON, "-vs", "vsites",    ffOPTRDMULT },
        { efXML, "-o" ,  "ffout",     ffWRITE },
        { efJSON, "-json", "gen_ff",  ffOPTWR }
    };
    gmx_bool useVdwCorr = options.vdwCorrection;
    std::vector<t_pargs> pa =
    {
        { "-nexcl", FALSE, etINT,  {&nexcl},
          "Number of bonds within which nonbonded interactions are excluded, the same for all of them." },
        { "-nexclqq", FALSE, etINT,  {&nexclqq},
          "HIDDENObsolete, use -nexcl. Accepted if it is the same as -nexclvdw." },
        { "-nexclvdw", FALSE, etINT,  {&nexclvdw},
          "HIDDENObsolete, use -nexcl. Accepted if it is the same as -nexclqq." },
        { "-epsilonr", FALSE, etREAL, {&epsilonr},
          "Relative dielectric constant. 1 is recommended for polarizable force fields, but maybe 1.7 might work for non-polarizable force fields instead of charge scaling." },
        { "-shellmaxdistance", FALSE, etREAL, {&shellMaxDistance},
          "Maximum distance (nm) between a shell and its core, 0 means no maximum. The shells of hyperpolarizable models may not need one." },
        { "-qsymm", FALSE, etBOOL, {&qsymm},
          "Implement symmetrization of charges in the force field" },
        { "-qdist", FALSE, etENUM, {qdn2.data()},
         "Charge distribution type, can be either" },
        { "-bondfn", FALSE, etENUM, {bondfn.data()},
          "Function to use for covalent bonds, can be either" },
        { "-anglefn", FALSE, etENUM, {anglefn.data()},
          "Function to use for covalent angles, can be either" },
        { "-dihfn", FALSE, etENUM, {dihfn.data()},
          "Function to use for proper dihedrals, can be either" },
        { "-vdwfn", FALSE, etENUM, {vdwfn.data()},
          "Function to use for Van der Waals interactions, can be either" },
        { "-vdwcorr", FALSE, etBOOL, {&useVdwCorr},
          "Add exchange correction to vsites" },
        { "-icfn", FALSE, etENUM, {icfn.data()},
          "Function to use for induction correction interactions, can be either" }
    };
    crule.addPargs(&pa);

    if (!parse_common_args(&argc, argv, 0, fnm.size(), fnm.data(),
                           pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 0;
    }
    // The obsolete options for the exclusions of Coulomb and van der Waals
    if (nexclqq >= 0 || nexclvdw >= 0)
    {
        if (nexclqq >= 0 && nexclvdw >= 0 && nexclqq != nexclvdw)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("The number of exclusions has to be the same for Coulomb (%d) and van der Waals (%d), use -nexcl",
                                                               nexclqq, nexclvdw).c_str()));
        }
        nexcl = std::max(nexclqq, nexclvdw);
    }

    options.atomTypesFile = opt2fn("-f", fnm.size(), fnm.data());
    options.outputFile    = opt2fn("-o", fnm.size(), fnm.data());
    if (opt2bSet("-vs", fnm.size(), fnm.data()))
    {
        const auto fns = opt2fns("-vs", fnm.size(), fnm.data());
        options.virtualSiteFiles.assign(fns.begin(), fns.end());
    }
    const char *jsonFile = opt2fn_null("-json", fnm.size(), fnm.data());
    options.jsonFile            = jsonFile ? jsonFile : "";
    options.nexcl               = nexcl;
    options.epsilonr            = epsilonr;
    options.shellMaxDistance    = shellMaxDistance;
    options.symmetricCharges    = qsymm;
    options.chargeDistribution  = distributionOf(qdn2[0]);
    options.bondPotential       = potentialOf(bondfn[0], "-bondfn");
    options.anglePotential      = potentialOf(anglefn[0], "-anglefn");
    options.dihedralPotential   = potentialOf(dihfn[0], "-dihfn");
    options.vdwPotential        = potentialOf(vdwfn[0], "-vdwfn");
    options.vdwCorrection       = useVdwCorr;
    options.inductionCorrection = potentialOf(icfn[0], "-icfn");
    options.combinationRules    = crule.rules();
    runGenFF(options);
    return 0;
}

} // namespace alexandria
