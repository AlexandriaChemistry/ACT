/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023,2024
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

#include <cstdlib>

#include <map>
#include <set>

#include "act/alexandria/alex_modules.h"
#include "act/basics/atomprops.h"
#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/forces/combinationrules.h"
#include "act/forces/forcecomputer.h"
#include "act/forcefield/act_checksum.h"
#include "act/forcefield/combruleutil.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/textreader.h"

namespace alexandria
{

static void add_symm_charges(ForceField *pd)
{
    auto            actdata  = getenv("ACTDATA");
    auto            filename = gmx::formatString("%s/symmetric_charges.csv", actdata);
    gmx::TextReader tr(filename);
    std::string     tmp;
    while (tr.readLine(&tmp))
    {
        auto ptr = split(tmp, ',');
        if (3 == ptr.size())
        {
            pd->addSymcharges(ptr[0], ptr[1],
                              my_atoi(ptr[2].c_str(), "numattch"));

        }
    }
}

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
        // There should be at least 37 columns, but could be more.
        if (37 <= ptr.size())
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

static bool minmaxmut(const std::string                  &atomtype,
                      std::map<std::string, std::string> &myatype,
                      const std::string                  &key,
                      double                             *min,
                      double                             *max,
                      Mutability                         *mut)
{
    auto minkey = key + "_min";
    auto maxkey = key + "_max";
    auto mutkey = key + "_mutability";

    *min = my_atof(myatype[minkey].c_str(), minkey.c_str());
    *max = my_atof(myatype[maxkey].c_str(), maxkey.c_str());

    if (!nameToMutability(myatype[mutkey], mut))
    {
        return false;
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Invalid Mutability '%s' for %s in atomtype %s", myatype[mutkey].c_str(), key.c_str(), atomtype.c_str()).c_str()));
    }
    return true;
}

static void add_vsites(const char *vsfile,
                       ForceField *pd)
{
    if (nullptr == vsfile)
    {
        return;
    }
    gmx::TextReader          tr(vsfile);
    std::string              tmp;
    std::vector<std::string> value;
    int                      lineno = 1;
    std::map<InteractionType, ForceFieldParameterList> i2f = {
        { InteractionType::VSITE2,     ForceFieldParameterList("vsite2",     CanSwap::Vsite2) },
        { InteractionType::VSITE2FD,   ForceFieldParameterList("vsite2fd",   CanSwap::Vsite2) },
        { InteractionType::VSITE3,     ForceFieldParameterList("vsite3",     CanSwap::Vsite3) },
        { InteractionType::VSITE3FD,   ForceFieldParameterList("vsite3fd",   CanSwap::No)     },
        { InteractionType::VSITE3FAD,  ForceFieldParameterList("vsite3fad",  CanSwap::No)     },
        { InteractionType::VSITE3OUT,  ForceFieldParameterList("vsite3out",  CanSwap::No)     },
        { InteractionType::VSITE3OUTS, ForceFieldParameterList("vsite3outs", CanSwap::No)     }
    };
    while (tr.readLine(&tmp))
    {
        auto ptr = split(tmp, '|');
        Identifier vsid(ptr[0]);
        if (!pd->hasParticleType(vsid))
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown vsite type %s on line %d in file %s",
                                                               ptr[0].c_str(), lineno, vsfile).c_str()));
        }
        auto itype = stringToInteractionType(ptr[1]);

        switch (itype)
        {
        case InteractionType::VSITE2:
        case InteractionType::VSITE2FD:
            {
                std::string myId = ptr[2] + "!" + ptr[0];
                Identifier vs(itype, myId, CanSwap::Vsite2);
                double amin = my_atof(ptr[3], "vsite_parameter_min");
                double amax = my_atof(ptr[4], "vsite_parameter_max");
                ForceFieldParameter vs2param("", (amin+amax)/2, 0, 0, amin, amax,
                                             Mutability::Bounded, false, false);
                if (InteractionType::VSITE2 == itype)
                {
                    i2f[itype].addParameter(vs, vsite2_name[vsite2A], vs2param);
                }
                else
                {
                    i2f[itype].addParameter(vs, vsite2fd_name[vsite2fdA], vs2param);
                }
            }
            break;
        case InteractionType::VSITE3:
        case InteractionType::VSITE3FD:
            {
                std::string myId = ptr[2] + "!" + ptr[0];
                
                Identifier vs(itype, myId, CanSwap::Vsite3);
                double amin = my_atof(ptr[3], "vsite_parameter_min");
                double amax = my_atof(ptr[4], "vsite_parameter_max");
                double bmin = my_atof(ptr[5], "vsite_parameter_min");
                double bmax = my_atof(ptr[6], "vsite_parameter_max");
                ForceFieldParameter vs3param_a("", (amin+amax)/2, 0, 0, amin, amax,
                                               Mutability::Bounded, false, false);
                ForceFieldParameter vs3param_b("", (bmin+bmax)/2, 0, 0, bmin, bmax,
                                               Mutability::Bounded, false, false);
                i2f[itype].addParameter(vs, vsite3_name[vsite3A], vs3param_a);
                i2f[itype].addParameter(vs, vsite3fd_name[vsite3fdB], vs3param_b);
            }
            break;
        case InteractionType::VSITE3OUT:
        case InteractionType::VSITE3OUTS:
            {
                std::string myId = ptr[2] + "!" + ptr[0];
                Identifier vs(itype, myId, CanSwap::No);

                double amin = my_atof(ptr[3], "vsite_parameter_min");
                double amax = my_atof(ptr[4], "vsite_parameter_max");
                double bmin = my_atof(ptr[5], "vsite_parameter_min");
                double bmax = my_atof(ptr[6], "vsite_parameter_max");
                ForceFieldParameter vs3outparam_a("", (amin+amax)/2, 0, 0, amin, amax,
                                                  Mutability::Bounded, false, false);
                if (InteractionType::VSITE3OUT == itype)
                {
                    double cmin = my_atof(ptr[7], "vsite_parameter_min");
                    double cmax = my_atof(ptr[8], "vsite_parameter_max");
                    ForceFieldParameter vs3outparam_b("", (bmin+bmax)/2, 0, 0, bmin, bmax,
                                                      Mutability::Bounded, false, false);
                    ForceFieldParameter vs3outparam_c("", (cmin+cmax)/2, 0, 0, cmin, cmax,
                                                      Mutability::Bounded, false, false);
                    i2f[itype].addParameter(vs, vsite3out_name[vsite3outA], vs3outparam_a);
                    i2f[itype].addParameter(vs, vsite3out_name[vsite3outB], vs3outparam_b);
                    i2f[itype].addParameter(vs, vsite3out_name[vsite3outC], vs3outparam_c);
                }
                else
                {
                    double cmin = my_atof(ptr[5], "vsite_parameter_min");
                    double cmax = my_atof(ptr[6], "vsite_parameter_max");
                    ForceFieldParameter vs3outparam_c("", (cmin+cmax)/2, 0, 0, cmin, cmax,
                                                      Mutability::Bounded, false, false);
                    i2f[itype].addParameter(vs, vsite3outs_name[vsite3outsA], vs3outparam_a);
                    i2f[itype].addParameter(vs, vsite3outs_name[vsite3outsC], vs3outparam_c);
                }
            }
            break;
        default:
            break;
        }

        lineno += 1;
    }
    for(const auto &iii: i2f)
    {
        if (!iii.second.empty())
        {
            pd->addForces(interactionTypeToString(iii.first), iii.second);
        }
    }
}

int gen_ff(int argc, char*argv[])
{
    std::vector<const char *> desc =
    {
        "gen_ff generates a new force field file based on a specification",
        "provided by the user in one or more csv files and command line options.[PAR]",
        "Most important is the atomtypes file for which there are some",
        "examples. In addition a file containing virtual sites can be provided.",
        "with the [TT]-vs[tt] option. Default is not to add virtual sites.[PAR]"
        "In order to work the ACTDATA environment variable should point to",
        "the directory where some input files are located, namely:[BR]",
        "atomprops.csv - containing properties for atoms.[BR]",
        "symmetric_charges.csv - containing groups of atoms that should have symmetric charges."
    };
    CombRuleUtil crule;
    crule.addInfo(&desc);
    gmx_output_env_t *oenv;
    int               nexclqq  = 0;
    int               nexclvdw = 0;
    double            epsilonr = 1;
    bool              qsymm    = true;
    const char *qdn2[]    = { nullptr, "Gaussian", "Point", "Slater", nullptr };
    const char *bondfn[]  = { nullptr, "CUBICBONDS", "BONDS", "MORSE", nullptr };
    const char *anglefn[] = { nullptr, "ANGLES", "UREYBRADLEY", nullptr };
    const char *dihfn[]   = { nullptr, "FOURDIHS", "PDIHS", nullptr };
    const char *vdwfn[]   = { nullptr, "WBHAM", "GBHAM", "LJ12_6", "LJ8_6", "LJ14_7", nullptr };

    std::vector<t_filenm> fnm = {
        { efCSV, "-f",   "atomtypes", ffREAD  },
        { efCSV, "-vs",  "vsites",    ffOPTRD },
        { efXML, "-o" ,  "ffout",     ffWRITE }
    };
    std::vector<t_pargs> pa =
    {
        { "-nexclqq", FALSE, etINT,  {&nexclqq},
          "Number of exclusions for Coulomb interactions, zero recommend for polarizable force fields." },
        { "-nexclvdw", FALSE, etINT,  {&nexclvdw},
          "Number of exclusions for Van der Waals interactions." },
        { "-epsilonr", FALSE, etREAL, {&epsilonr},
          "Relative dielectric constant. 1 is recommended for polarizable force fields, but maybe 1.7 might work for non-polarizable force fields instead of charge scaling." },
        { "-qsymm", FALSE, etBOOL, {&qsymm},
          "Implement symmetrization of charges in the force field" },
        { "-qdist", FALSE, etENUM, {qdn2},
          "Charge distribution type, can be either" },
        { "-bondfn", FALSE, etENUM, {bondfn},
          "Function to use for covalent bonds, can be either" },
        { "-anglefn", FALSE, etENUM, {anglefn},
          "Function to use for covalent angles, can be either" },
        { "-dihfn", FALSE, etENUM, {dihfn},
          "Function to use for proper dihedrals, can be either" },
        { "-vdwfn", FALSE, etENUM, {vdwfn},
          "Function to use for Van der Waals interactions, can be either" }
    };
    crule.addPargs(&pa);

    if (!parse_common_args(&argc, argv, 0, fnm.size(), fnm.data(),
                           pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 0;
    }

    auto table = read_atomtypes(opt2fn("-f", fnm.size(), fnm.data()));
    printf("There are %zu atom types in the force field with %zu properties.\n",
           table.size(), table.begin()->second.size());

    auto aprops = readAtomProps();
    printf("There are %zu element properties\n", aprops.size());

    alexandria::ForceField pd;
    // More stuff
    std::vector<std::string> options = {
        "acmtype", "bondtype", "element", "poltype", "row", "zetatype"
    };
    ForceFieldParameterList pols("Polarization", CanSwap::No);
    auto qdist = name2ChargeType(qdn2[0]);
    auto qpot  = chargeTypeToPotential(qdist);
    ForceFieldParameterList coulomb(potentialToString(qpot), CanSwap::Yes);
    coulomb.addOption("epsilonr", gmx_ftoa(epsilonr));
    coulomb.addOption("nexcl", gmx_itoa(nexclqq));
    ForceFieldParameterList vdw(vdwfn[0], CanSwap::Yes);
    // Combination rules
    crule.extract(&vdw);
    vdw.addOption("nexcl", gmx_itoa(nexclvdw));
    ForceFieldParameterList eem("", CanSwap::No);
    // Check for Point charges
    std::string ppp("Point");
    bool bPoint = ppp.compare(qdn2[0]) == 0;
    if (bPoint)
    {
        printf("You selected point charges. Oh dear...\n");
    }
    for(const auto &entry : table)
    {
        // Generate particle type
        int gmxtype = eptAtom;
        std::string elem = table[entry.first]["element"];
        if (elem == "X")
        {
            gmxtype = eptShell;
        }
        else if (elem == "VS")
        {
            gmxtype = eptVSite;
        }
        auto ptp = ParticleType(Identifier(entry.first),
                                table[entry.first]["comment"], gmxtype);

        // Now add the "options"
        for(const auto &opt : options)
        {
            if (!table[entry.first][opt].empty())
            {
                ptp.setOption(opt, table[entry.first][opt]);
            }
        }
        int atomnumber = 0;
        auto apropsptr = aprops.find(elem);
        if (aprops.end() != apropsptr)
        {
            atomnumber = apropsptr->second.atomnumber();
        }
        ptp.setOption("atomnumber", gmx::formatString("%d", atomnumber));
        ptp.setOption("vdwtype", entry.first);

        // Now "parameters"
        auto mass       = aprops.find(elem)->second.mass();
        ptp.addForceFieldParameter("mass",
                                   ForceFieldParameter("Da", mass, 0, 1, mass,
                                                       mass, Mutability::Fixed,
                                                       true, true));
        double vmin = my_atof(table[entry.first]["q_min"].c_str(), "q_min");
        double vmax = my_atof(table[entry.first]["q_max"].c_str(), "q_max");
        Mutability vmut;
        std::string qmut("charge_mutability");
        if (!nameToMutability(table[entry.first][qmut], &vmut))
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Invalid Mutability '%s' for charge in atomtype %s", table[entry.first][qmut].c_str(), entry.first.c_str()).c_str()));
        }
        ptp.addForceFieldParameter("charge",
                                   ForceFieldParameter("e", (vmin+vmax)/2,
                                                       0, 0, vmin,
                                                       vmax, vmut,
                                                       true, false));

        pd.addParticleType(ptp);

        auto myatype = table[entry.first];
        // Now add force field parameters
        if (eptShell == gmxtype)
        {
            if (minmaxmut(entry.first, myatype, "alpha", &vmin, &vmax, &vmut))
            {
                pols.addParameter(entry.first, "alpha",
                                  ForceFieldParameter("Angstrom3", (vmin+vmax)/2, 0, 0, vmin, vmax, vmut, true, true));
            }
        }
        // Charge distribution
        if (!myatype["zetatype"].empty() &&
            minmaxmut(entry.first, myatype, "zeta", &vmin, &vmax, &vmut))
        {
            if (bPoint)
            {
                vmin = vmax = 0;
                vmut = Mutability::Fixed;
            }
            coulomb.addParameter(myatype["zetatype"], "zeta",
                                 ForceFieldParameter("1/nm", (vmin+vmax)/2, 0, 0, vmin, vmax, vmut, true, true));
        }
        // Van der Waals
        {
            std::map<std::string, std::string> vdwlist;
            std::map<std::string, std::string> rename;
            switch (vdw.potential())
            {
            case Potential::LJ12_6:
                vdwlist = { { "sigma", "nm" }, { "epsilon", "kJ/mol" } };
                break;
            case Potential::LJ8_6:
                vdwlist = { { "sigma", "nm" }, { "epsilon", "kJ/mol" } };
                break;
            case Potential::LJ14_7:
                vdwlist = { { "sigma", "nm" }, { "epsilon", "kJ/mol" }, { "gamma", "" }, { "delta", "" }, { "aqt", "kJ/mol" }, { "bqt", "1/nm"} };
                break;
            case Potential::WANG_BUCKINGHAM:
                vdwlist = { { "sigma", "nm" }, { "epsilon", "kJ/mol" }, { "gamma", "" } };
                break;
            case Potential::GENERALIZED_BUCKINGHAM:
                vdwlist = { { "sigma", "nm" }, { "epsilon", "kJ/mol" }, { "gamma", "" }, { "delta", "" } };
                rename.insert({"sigma", "rmin"});
                break;
            default:
                GMX_THROW(gmx::InvalidInputError("Unknown function for Van der Waals interactions"));
            }
            for(const auto &vl : vdwlist)
            {
                if (minmaxmut(entry.first, myatype, vl.first, &vmin, &vmax, &vmut))
                {
                    std::string key = vl.first;
                    if (rename.find(key) != rename.end())
                    {
                        key = rename[key];
                    }
                    vdw.addParameter(entry.first, key,
                                     ForceFieldParameter(vl.second, (vmin+vmax)/2, 0, 0, vmin, vmax, vmut, true, true));
                }
            }
        }
        // EEM parameters
        {
            std::map<std::string, std::string> eemlist = {
                { "eta", "eV" },
                { "chi", "eV" }
            };
            for(const auto &vl : eemlist)
            {
                if (minmaxmut(entry.first, myatype, vl.first, &vmin, &vmax, &vmut))
                {
                    eem.addParameter(myatype["zetatype"], vl.first,
                                     ForceFieldParameter(vl.second, (vmin+vmax)/2, 0, 1, vmin, vmax, vmut, true, true));
                }
            }
        }
    }
    pd.addForces(interactionTypeToString(InteractionType::POLARIZATION), pols);
    pd.addForces(interactionTypeToString(InteractionType::COULOMB), coulomb);
    ForceFieldParameterList bonds(bondfn[0], CanSwap::Yes);
    pd.addForces(interactionTypeToString(InteractionType::BONDS), bonds);
    ForceFieldParameterList angles(anglefn[0], CanSwap::Yes);
    pd.addForces(interactionTypeToString(InteractionType::ANGLES), angles);
    ForceFieldParameterList linang("LINEAR_ANGLES", CanSwap::Linear);
    pd.addForces(interactionTypeToString(InteractionType::LINEAR_ANGLES), linang);
    ForceFieldParameterList idihs("IDIHS", CanSwap::Idih);
    pd.addForces(interactionTypeToString(InteractionType::IMPROPER_DIHEDRALS), idihs);
    ForceFieldParameterList pdihs(dihfn[0], CanSwap::Yes);
    pd.addForces(interactionTypeToString(InteractionType::PROPER_DIHEDRALS), pdihs);
    pd.addForces(interactionTypeToString(InteractionType::VDW), vdw);
    pd.addForces(interactionTypeToString(InteractionType::ELECTRONEGATIVITYEQUALIZATION), eem);
    // Virtual sites
    add_vsites(opt2fn_null("-vs", fnm.size(), fnm.data()), &pd);
    if (qsymm)
    {
        add_symm_charges(&pd);
    }
    pd.updateTimeStamp();
    pd.updateCheckSum();
    writeForceField(opt2fn("-o", fnm.size(), fnm.data()), &pd, 0);

    return 0;
}

} // namespace alexandria
