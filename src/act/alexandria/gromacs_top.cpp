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

#include "gromacs_top.h"

#include <map>

#include "act/alexandria/topology.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parametername.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

static int get_subtype(directive d, int ftype)
{
    int i;
    for (i = 1; i < 20; i++)
    {
        if (!(d == d_angles && i == 7))
        {
            if (ifunc_index(d, i) == ftype)
            {
                return i;
            }
        }
    }
    return 1;
}

static void print_bondeds(MsgHandler                *msg_handler,
                          FILE                      *out,
                          directive                  d,
                          int                        ftype,
                          const TopologyEntryVector &entries)
{
    if (entries.empty())
    {
        return;
    }
    fprintf(out, "[ %s ]\n", dir2str(d));
    fprintf(out, ";atom i");
    for (int j = 1; (j < NRAL(ftype)); j++)
    {
        fprintf(out, "  %5c", j+'i');
    }
    fprintf(out, "   type  parameters\n");
    int subtype = get_subtype(d, ftype);
    for (auto &entry : entries)
    {
        for (auto &j : entry->atomIndices())
        {
            fprintf(out, "  %5d", 1+j);
        }
        fprintf(out, "  %5d", subtype);
        auto params = entry->params();
        switch (ftype)
        {
        case F_BONDS:
            fprintf(out, "  %10g  %10g ; %s\n",
                    params[bondLENGTH],
                    params[bondKB],
                    entry->id().id().c_str());
            break;
        case F_MORSE:
            fprintf(out, "  %10g  %10g  %10g ; %s\n",
                    params[morseLENGTH],
                    params[morseDE],
                    params[morseBETA],
                    entry->id().id().c_str());
            break;
        case F_ANGLES:
            fprintf(out, "  %10g  %10g ; %s\n",
                    params[angleKT],
                    params[angleANGLE],
                    entry->id().id().c_str());
            break;
        case F_IDIHS:
            fprintf(out, "  0.0 %10g ; %s\n",
                    params[idihKPHI],
                    entry->id().id().c_str());
            break;
        default: // writes warning
            msg_handler->msg(ACTStatus::Warning,
                             gmx::formatString("Unsupported function type %s for writing GROMACS topology",
                                               interaction_function[ftype].name).c_str());
        }
    }
    fprintf(out, "\n");
}

static void print_excl(FILE *out, const Topology *top)
{
    // For GROMACS type output files we only look at VANDERWAALS
    auto itype = InteractionType::VDW;
    if (!top->hasExclusions(itype))
    {
        return;
    }
    auto excls = top->exclusions(itype);
    if (excls.empty())
    {
        return;
    }
    fprintf (out, "[ %s ]\n", dir2str(d_exclusions));
    fprintf (out, "; %4s    %s\n", "i", "excluded from i");
    for (size_t i = 0; i < excls.size(); i++)
    {
        if (!excls[i].empty())
        {
            fprintf (out, "%6lu", i+1);
            for (size_t j = 0; j < excls[i].size(); j++)
            {
                fprintf (out, " %5d", excls[i][j]+1);
            }
            fprintf (out, "\n");
        }
    }
    fprintf (out, "\n");
}

static void print_top_system(FILE *out, const char *title)
{
    fprintf(out, "[ %s ]\n", dir2str(d_system));
    fprintf(out, "; Name\n");
    fprintf(out, "%s\n\n", title[0] ? title : "Protein");
}

void print_top_mols(FILE *out, const char *title,
                    int nmol, t_mols *mols)
{
    print_top_system(out, title);

    if (nmol)
    {
        fprintf(out, "[ %s ]\n", dir2str(d_molecules));
        fprintf(out, "; %-15s %5s\n", "Compound", "#mols");
        for (int i = 0; i < nmol; i++)
        {
            fprintf(out, "%-15s %5d\n", mols[i].name, mols[i].nr);
        }
    }
}

static void print_atoms(FILE                           *out,
                        const std::vector<ActAtom>     &atoms,
                        const std::vector<std::string> &residueNames)
{
    fprintf(out, "[ %s ]\n", dir2str(d_atoms));
    fprintf(out, "; %4s %10s %6s %7s%6s %6s %10s %10s %6s %10s %10s\n",
            "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass", "typeB", "chargeB", "massB");

    double qtot  = 0;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        int atomnr = 1+i;
        /* This is true by construction, but static analysers don't know */
        fprintf(out, "%6d %10s %6d%c %5s %6s %6d %10g %10g",
                atomnr,
                atoms[i].ffType().c_str(),
                atoms[i].residueNumber()+1, ' ',
                residueNames[atoms[i].residueNumber()].c_str(),
                atoms[i].name().c_str(),
                atomnr,
                atoms[i].charge(), 
                atoms[i].mass());
        // Accumulate the total charge to help troubleshoot issues.
        qtot += static_cast<double>(atoms[i].charge());
        // Round it to zero if it is close to zero, because
        // printing -9.34e-5 confuses users.
        if (fabs(qtot) < 0.0001)
        {
            qtot = 0;
        }
        // Write the total charge for the last atom of the system
        // and/or residue, because generally that's where it is
        // expected to be an integer.
        if (i == atoms.size()-1)
        {
            fprintf(out, "   ; qtot %.4g\n", qtot);
        }
        else
        {
            fprintf(out, "   ; %.4g\n", qtot);
        }
    }
    fprintf(out, "\n");
}

void write_top(MsgHandler       *msg_handler,
               FILE             *out,
               char             *molname,
               const Topology   *topology,
               const ForceField *pd)
{
    std::map<int, directive> toPrint = {
        { F_POLARIZATION, d_polarization },
        { F_VSITEN,       d_vsitesn },
        { F_VSITE2,       d_vsites2 },
        { F_VSITE3,       d_vsites3 },
        { F_VSITE3FD,     d_vsites3 },
        { F_VSITE3FAD,    d_vsites3 },
        { F_VSITE3OUT,    d_vsites3 },
        { F_VSITE4FD,     d_vsites4 },
        { F_VSITE4FDN,    d_vsites4 }
    };
    auto myAtoms = topology->atoms();
    if (!myAtoms.empty())
    {
        std::vector<int> cgnr;
        for(size_t i = 0; i < myAtoms.size(); i++)
        {
            cgnr.push_back(i+1);
        }
        fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        int nexcl;
        if (!ffOption(*pd, InteractionType::VDW, "nexcl", &nexcl))
        {
            nexcl = 0;
        }
        fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", nexcl);
        print_atoms(out, myAtoms, topology->residueNames());
        for (auto &fs : pd->forcesConst())
        {
            auto iType = fs.first;
            if (!topology->hasEntry(iType))
            {
                continue;
            }
            auto fType = potentialToGromacsType(fs.second.potential());
            if (InteractionType::BONDS == fs.first)
            {
                print_bondeds(msg_handler, out, d_bonds, fType, topology->entry(iType));
            }
            else if (InteractionType::ANGLES == iType || InteractionType::LINEAR_ANGLES == iType)
            {
                print_bondeds(msg_handler, out, d_angles, fType, topology->entry(iType));
            }
            else if (InteractionType::PROPER_DIHEDRALS == iType || InteractionType::IMPROPER_DIHEDRALS == iType)
            {
                print_bondeds(msg_handler, out, d_dihedrals, fType, topology->entry(iType));
            }
            else if (toPrint.end() != toPrint.find(fType))
            {
                print_bondeds(msg_handler, out, toPrint.find(fType)->second, fType, topology->entry(iType));
            }
        }
        print_excl(out, topology);
    }
}

void print_top_header(FILE                    *fp,
                      const ForceField        *pd,
                      bool                     bPol,
                      const std::vector<std::string> &commercials,
                      bool                     bItp)
{
    std::string   gt_old, gt_type;
    auto qt          = pd->findForcesConst(InteractionType::ELECTROSTATICS);
    auto iChargeType = potentialToChargeType(qt.potential());

    fprintf(fp, ";\n");
    fprintf(fp, "; Topology generated by alexandria gentop.\n");
    fprintf(fp, "; WARNING: do not use for simulations!\n");
    fprintf(fp, "; Watch this space for information & commercials.\n");
    for (auto i = commercials.begin(); (i < commercials.end()); ++i)
    {
        fprintf(fp, "; %s\n", i->c_str());
    }
    fprintf(fp, ";\n");
    if (!bItp)
    {
        fprintf(fp, "[ defaults ]\n");
        fprintf(fp, "; nbfunc         comb-rule       gen-pairs       fudgeLJ     fudgeQQ\n");
        auto ftype = pd->findForcesConst(InteractionType::VDW).potential();
        std::string ff("BHAM");
        if (ftype == Potential::LJ12_6)
        {
            ff.assign("LJ");
        }
        auto myfs = pd->findForcesConst(InteractionType::VDW);
        const char *crule = "combination_rule";
        std::string combRule("geometric");
        if (myfs.optionExists(crule))
        {
            combRule = myfs.optionValue(crule);
        }
        fprintf(fp, "%-15s  %-15s no           %10g  %10g\n\n",
                ff.c_str(),
                combRule.c_str(),
                1.0, 1.0);

        fprintf(fp, "[ atomtypes ]\n");
        fprintf(fp, "%-7s%-6s  %6s  %11s  %10s  %5s %-s\n",
                ";atype ", "btype", "at.num", "mass", "charge", "ptype",
                "Van_der_Waals");

        gt_old = "";

        auto vdw = pd->findForcesConst(InteractionType::VDW);
        for (const auto &aType : pd->particleTypesConst())
        {
            gt_type    = aType.second.id().id();
            std::string bType;
            if (aType.second.hasInteractionType(InteractionType::BONDS))
            {
                bType = aType.second.interactionTypeToIdentifier(InteractionType::BONDS).id();
            }
            if ((0 ==  gt_old.size()) || (gt_old.compare(gt_type) != 0))
            {
                auto vdwtype = aType.second.interactionTypeToIdentifier(InteractionType::VDW);
                double sigma = 0, epsilon = 0, gamma = 0;
                if (!vdwtype.id().empty() && vdw.parameterExists(vdwtype))
                {
                    auto myvdw = vdw.findParametersConst(vdwtype);
                    sigma      = myvdw["sigma"].value();
                    epsilon    = myvdw["epsilon"].value();
                    gamma      = myvdw["gamma"].value();
                }
                fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f %s %g %g %g\n",
                        gt_type.c_str(), 
                        !bType.empty() ? bType.c_str() : gt_type.c_str(), 
                        aType.second.atomnumber(), 
                        aType.second.mass(), 0.0,
                        actParticleToString(aType.second.apType()).c_str(),
                        sigma, epsilon, gamma);
                if (false && bPol)
                {
                    auto sgt_type= aType.second.interactionTypeToIdentifier(InteractionType::POLARIZATION);
                    if (strcasecmp(ff.c_str(), "LJ") == 0)
                    {
                        fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  S     0  0\n",
                                sgt_type.id().c_str(),
                                sgt_type.id().c_str(),
                                0, 0.0, 0.0);
                    }
                    else
                    {
                        fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  S     0  0  0\n",
                                sgt_type.id().c_str(),
                                sgt_type.id().c_str(),
                                0, 0.0, 0.0);
                    }
                }
            }
            gt_old = gt_type;
        }
        fprintf(fp, "\n");
        bool printZeta = false;
        if (iChargeType != ChargeType::Point && printZeta)
        {
            auto eem = pd->findForcesConst(InteractionType::ELECTROSTATICS);
            fprintf(fp, "[ distributed_charges ]\n");
            for (const auto &atype : pd->particleTypesConst())
            {
                auto ztype     = atype.second.interactionTypeToIdentifier(InteractionType::ELECTROSTATICS);
                auto eep       = eem.findParametersConst(ztype);
                if (ChargeType::Slater == iChargeType)
                {
                    fprintf(fp, "%-7s  2  %d  %g\n", atype.second.id().id().c_str(),
                            atype.second.row(), eep["zeta"].value());
                }
                else if (ChargeType::Gaussian == iChargeType)
                {
                    fprintf(fp, "%-7s  1  %g\n", atype.second.id().id().c_str(),
                            eep["zeta"].value());
                }
            }
            fprintf(fp, "\n");
        }
    }
}

} // namespace alexandria
