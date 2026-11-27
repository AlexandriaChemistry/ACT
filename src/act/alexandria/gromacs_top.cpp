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

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "act/alexandria/pdbwriter.h"
#include "act/basics/act_particle.h"
#include "act/basics/chargemodel.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/potential.h"
#include "act/molprop/topologyentry.h"
#include "act/topology/topology.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs_top.h"

namespace alexandria
{

/*! \brief Look up the function number of GROMACS for an interaction
 * \param[in] d     The directive, bonds, angles or dihedrals
 * \param[in] ftype The function type
 * \return the number, or 1 if GROMACS has no such function for this
 * directive, as for the Hua and the quartic bonds of ACT
 */
static int get_subtype(directive d, int ftype)
{
    // The numbers that GROMACS has for each directive: ifunc_index
    // stops the program for others
    static const std::map<directive, std::vector<int>> numbers = {
        { d_bonds,     { 1, 2, 3, 4, 5, 6, 7, 8, 9 } },
        { d_angles,    { 1, 2, 3, 4, 5, 6, 8, 9, 10 } },
        { d_dihedrals, { 1, 2, 3, 4, 5, 8, 9, 10, 11 } }
    };
    const auto nn = numbers.find(d);
    if (numbers.end() != nn)
    {
        for (auto i : nn->second)
        {
            if (ifunc_index(d, i) == ftype)
            {
                return i;
            }
        }
    }
    return 1;
}

//! \brief Print bonded interactions to topology file
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
    // The Fourier dihedrals of ACT are the Ryckaert-Bellemans dihedrals of
    // GROMACS, sums of powers of the cosine of phi - 180 degrees
    int subtype      = (F_FOURDIHS == ftype) ? 3 : get_subtype(d, ftype);
    int nunsupported = 0;
    // The energy of the bonds at their minimum, which GROMACS does not have
    double offset    = 0;
    for (auto &entry : entries)
    {
        if (F_BONDS == ftype)
        {
            offset += entry->params()[bondENERGY];
        }
        else if (F_MORSE == ftype)
        {
            offset += entry->params()[morseD0] - entry->params()[morseDE];
        }
    }
    if (offset != 0)
    {
        fprintf(out, "; The energy of these bonds in ACT is that in GROMACS plus %.9g kJ/mol,\n", offset);
        fprintf(out, "; since GROMACS has an energy of zero for a bond at its minimum.\n");
    }
    for (auto &entry : entries)
    {
        auto params = entry->params();
        if (F_CROSS_BOND_ANGLES == ftype)
        {
            // One angle in ACT, up to three in GROMACS: the bond-angle
            // cross term, the bond-bond cross term and a Urey-Bradley
            // angle without angle part for the harmonic potential on the
            // distance between the outer atoms. The cross term with the
            // angle is always there, the others if they have a force
            // constant.
            const auto printAtoms = [&](int gmxFtype)
            {
                for (auto &j : entry->atomIndices())
                {
                    fprintf(out, "  %5d", 1+j);
                }
                fprintf(out, "  %5d", get_subtype(d, gmxFtype));
            };
            printAtoms(F_CROSS_BOND_ANGLES);
            fprintf(out, "  %10g  %10g  %10g  %10g ; %s\n",
                    params[cbaR1E], params[cbaR2E], params[cbaR3E], params[cbaKRT],
                    entry->id().id().c_str());
            if (params[cbaKRR] != 0)
            {
                printAtoms(F_CROSS_BOND_BONDS);
                fprintf(out, "  %10g  %10g  %10g ; %s\n",
                        params[cbaR1E], params[cbaR2E], params[cbaKRR],
                        entry->id().id().c_str());
            }
            if (params[cbaKUB] != 0)
            {
                printAtoms(F_UREY_BRADLEY);
                fprintf(out, "  %10g  %10g  %10g  %10g ; %s\n",
                        0.0, 0.0, params[cbaR3E], params[cbaKUB],
                        entry->id().id().c_str());
            }
            continue;
        }
        for (auto &j : entry->atomIndices())
        {
            fprintf(out, "  %5d", 1+j);
        }
        fprintf(out, "  %5d", subtype);
        switch (ftype)
        {
        case F_BONDS:
            fprintf(out, "  %12.9g  %12.9g ; %s\n",
                    params[bondLENGTH],
                    params[bondKB],
                    entry->id().id().c_str());
            break;
        case F_MORSE:
            fprintf(out, "  %12.9g  %12.9g  %12.9g ; %s\n",
                    params[morseLENGTH],
                    params[morseDE],
                    params[morseBETA],
                    entry->id().id().c_str());
            break;
        case F_ANGLES:
            // The angles are in radians in ACT and in degrees in GROMACS
            fprintf(out, "  %12.9g  %12.9g ; %s\n",
                    params[angleANGLE]*RAD2DEG,
                    params[angleKT],
                    entry->id().id().c_str());
            break;
        case F_UREY_BRADLEY:
            fprintf(out, "  %12.9g  %12.9g  %12.9g  %12.9g ; %s\n",
                    params[ubANGLE]*RAD2DEG,
                    params[ubKT],
                    params[ubR13],
                    params[ubKUB],
                    entry->id().id().c_str());
            break;
        case F_LINEAR_ANGLES:
            fprintf(out, "  %12.9g  %12.9g ; %s\n",
                    params[linangA],
                    params[linangKLIN],
                    entry->id().id().c_str());
            break;
        case F_QUARTIC_ANGLES:
            // The angle in degrees and the five coefficients, as in ACT
            fprintf(out, "  %12.9g", params[qaANGLE]*RAD2DEG);
            for (int k = qaC0; k < qaNR; k++)
            {
                fprintf(out, "  %12.9g", params[k]);
            }
            fprintf(out, " ; %s\n", entry->id().id().c_str());
            break;
        case F_PDIHS:
            fprintf(out, "  %12.9g  %12.9g  %d ; %s\n",
                    params[pdihANGLE]*RAD2DEG,
                    params[pdihKP],
                    static_cast<int>(std::lround(params[pdihMULT])),
                    entry->id().id().c_str());
            break;
        case F_FOURDIHS:
            for (int k = fdihC0; k < fdihNR; k++)
            {
                fprintf(out, "  %12.9g", k < static_cast<int>(params.size()) ? params[k] : 0.0);
            }
            fprintf(out, " ; %s\n", entry->id().id().c_str());
            break;
        case F_IDIHS:
            fprintf(out, "  0.0 %12.9g ; %s\n",
                    params[idihKPHI],
                    entry->id().id().c_str());
            break;
        default:
            // Without parameters, such that GROMACS does not accept the
            // topology
            fprintf(out, " ; %s, unsupported function %s\n", entry->id().id().c_str(),
                    interaction_function[ftype].name);
            nunsupported += 1;
        }
    }
    if (nunsupported > 0)
    {
        msg_handler->msg(ACTStatus::Warning,
                         gmx::formatString("GROMACS does not support %s: %d interactions are in the topology without parameters and GROMACS will not accept it",
                                           interaction_function[ftype].name, nunsupported));
    }
    fprintf(out, "\n");
}

/*! \brief Print the polarization of the shells to a topology file
 *
 * The polarizability is in nm^3 as in GROMACS. With a force constant
 * for the anharmonic term the function is the anharmonic polarization of
 * GROMACS.
 * \param[in] out      The file
 * \param[in] entries  The pairs of a core and its shell
 */
static void print_polarization(FILE                      *out,
                               const TopologyEntryVector &entries)
{
    if (entries.empty())
    {
        return;
    }
    fprintf(out, "[ polarization ]\n");
    fprintf(out, ";atom i      j   type  parameters\n");
    for (auto &entry : entries)
    {
        auto &indices = entry->atomIndices();
        auto &params  = entry->params();
        fprintf(out, "  %5d  %5d", 1+indices[0], 1+indices[1]);
        if (params.size() > static_cast<size_t>(polFCHYPER) && params[polFCHYPER] > 0)
        {
            fprintf(out, "  %5d  %12.9g  %12.9g  %12.9g ; %s\n", 2,
                    params[polALPHA], params[polRHYPER], params[polFCHYPER],
                    entry->id().id().c_str());
        }
        else
        {
            fprintf(out, "  %5d  %12.9g ; %s\n", 1, params[polALPHA],
                    entry->id().id().c_str());
        }
    }
    fprintf(out, "\n");
}

/*! \brief Print the virtual sites to a topology file
 *
 * In ACT the virtual site is the last particle of an entry and the
 * central atom of a site on three atoms the second, in GROMACS the site
 * comes first, followed by the central atom. The parameters are
 * converted to the definitions of GROMACS, see VsiteHandler for those
 * of ACT:
 * - VSITE1 is a site on top of an atom (virtual_sites1),
 * - VSITE2 and VSITE2FD are the functions 1 and 2 of virtual_sites2,
 * - VSITE3 and VSITE3S are function 1 of virtual_sites3,
 * - VSITE3FD is function 2 of virtual_sites3, where ACT adds a times
 *   the vector between the outer atoms to the vector from the central
 *   atom to the first and GROMACS takes a point between the two outer
 *   atoms, such that a becomes a/(1+a),
 * - VSITE3OUT and VSITE3OUTS are function 4 of virtual_sites3, with the
 *   sign of the site in the last parameter,
 * - VSITE4, VSITE4S and VSITE4S3 are linear combinations of four
 *   atoms, which GROMACS has as the center of weights of virtual_sitesn
 *   (function 3) if none of the weights is negative.
 * \param[in] msg_handler For warnings
 * \param[in] out         The file
 * \param[in] topology    The molecule
 */
static void print_vsites(MsgHandler     *msg_handler,
                         FILE           *out,
                         const Topology *topology)
{
    // The lines of each directive, in the order of GROMACS
    const std::vector<std::string>                  directives = {
        "virtual_sites1", "virtual_sites2", "virtual_sites3", "virtual_sitesn"
    };
    std::map<std::string, std::vector<std::string>> lines;
    const auto unsupported = [&](InteractionType itype, const TopologyEntry &entry, const std::string &why)
    {
        msg_handler->msg(ACTStatus::Warning,
                         gmx::formatString("Virtual site %s of type %s cannot be written to a GROMACS topology: %s",
                                           entry.id().id().c_str(),
                                           interactionTypeToString(itype).c_str(), why.c_str()));
        // The site without atoms or parameters, such that GROMACS does
        // not accept the topology
        lines["virtual_sites3"].push_back(gmx::formatString("  %5d ; %s of type %s is not supported: %s",
                                                            1 + entry.atomIndices().back(),
                                                            entry.id().id().c_str(),
                                                            interactionTypeToString(itype).c_str(), why.c_str()));
    };
    for (const auto &ee : topology->entries())
    {
        auto itype = ee.first;
        if (!isVsite(itype))
        {
            continue;
        }
        for (const auto &entry : ee.second)
        {
            auto       &ind    = entry->atomIndices();
            auto       &params = entry->params();
            // The site is the last particle
            const int   vs     = 1 + ind.back();
            const auto  id     = entry->id().id();
            // A line with the site, atoms in the order of GROMACS, the
            // function and the parameters
            const auto  line   = [&](const std::vector<int> &atoms, int func, const std::vector<double> &pp)
            {
                std::string str = gmx::formatString("  %5d", vs);
                for (auto a : atoms)
                {
                    str += gmx::formatString("  %5d", 1 + ind[a]);
                }
                str += gmx::formatString("  %5d", func);
                for (auto p : pp)
                {
                    str += gmx::formatString("  %12.9g", p);
                }
                return str + " ; " + id;
            };
            // A linear combination of atoms with weights
            const auto  weighted = [&](const std::vector<double> &weights)
            {
                for (auto w : weights)
                {
                    if (w < 0)
                    {
                        unsupported(itype, *entry->self(),
                                    "GROMACS has no negative weights for a site on four atoms");
                        return;
                    }
                }
                std::string str = gmx::formatString("  %5d  %5d", vs, 3);
                for (size_t k = 0; k < weights.size(); k++)
                {
                    str += gmx::formatString("  %5d  %12.9g", 1 + ind[k], weights[k]);
                }
                lines["virtual_sitesn"].push_back(str + " ; " + id);
            };
            switch (itype)
            {
            case InteractionType::VSITE1:
                lines["virtual_sites1"].push_back(line({ 0 }, 1, {}));
                break;
            case InteractionType::VSITE2:
                lines["virtual_sites2"].push_back(line({ 0, 1 }, 1, { params[vsite2A] }));
                break;
            case InteractionType::VSITE2FD:
                lines["virtual_sites2"].push_back(line({ 0, 1 }, 2, { params[vsite2fdA] }));
                break;
            case InteractionType::VSITE3:
                lines["virtual_sites3"].push_back(line({ 1, 0, 2 }, 1, { params[vsite3A], params[vsite3B] }));
                break;
            case InteractionType::VSITE3S:
                lines["virtual_sites3"].push_back(line({ 1, 0, 2 }, 1, { params[vsite3sA], params[vsite3sA] }));
                break;
            case InteractionType::VSITE3FD:
                {
                    double denom = 1 + params[vsite3fdA];
                    if (denom == 0)
                    {
                        unsupported(itype, *entry->self(), "the parameter a is -1");
                    }
                    else
                    {
                        // The direction turns around if 1+a is negative
                        double b = denom > 0 ? params[vsite3fdB] : -params[vsite3fdB];
                        lines["virtual_sites3"].push_back(line({ 1, 0, 2 }, 2, { params[vsite3fdA]/denom, b }));
                    }
                }
                break;
            case InteractionType::VSITE3OUT:
                {
                    auto vs3out = static_cast<const Vsite3OUT *>(entry->self());
                    lines["virtual_sites3"].push_back(line({ 1, 0, 2 }, 4,
                                                           { params[vsite3outA], params[vsite3outB],
                                                             vs3out->sign()*params[vsite3outC] }));
                }
                break;
            case InteractionType::VSITE3OUTS:
                {
                    auto vs3out = static_cast<const Vsite3OUT *>(entry->self());
                    lines["virtual_sites3"].push_back(line({ 1, 0, 2 }, 4,
                                                           { params[vsite3outsA], params[vsite3outsA],
                                                             vs3out->sign()*params[vsite3outsC] }));
                }
                break;
            case InteractionType::VSITE4:
                weighted({ params[vsite4A], 1 - params[vsite4A] - params[vsite4B] - params[vsite4C],
                           params[vsite4B], params[vsite4C] });
                break;
            case InteractionType::VSITE4S:
                weighted({ params[vsite4sA], 1 - params[vsite4sA] - 2*params[vsite4sB],
                           params[vsite4sB], params[vsite4sB] });
                break;
            case InteractionType::VSITE4S3:
                // The central atom comes first here
                weighted({ 1 - 3*params[vsite4s3A], params[vsite4s3A],
                           params[vsite4s3A], params[vsite4s3A] });
                break;
            default:
                unsupported(itype, *entry->self(), "ACT does not construct this type either");
            }
        }
    }
    for (const auto &d : directives)
    {
        const auto ll = lines.find(d);
        if (lines.end() == ll)
        {
            continue;
        }
        fprintf(out, "[ %s ]\n", d.c_str());
        if (d == "virtual_sitesn")
        {
            fprintf(out, "; site   type  atoms with their weights\n");
        }
        else
        {
            fprintf(out, "; site  atoms  type  parameters\n");
        }
        for (const auto &l : ll->second)
        {
            fprintf(out, "%s\n", l.c_str());
        }
        fprintf(out, "\n");
    }
}

//! \brief Print exclusions to a topology file
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
        // Each pair once
        std::vector<int> later;
        for (const auto j : excls[i])
        {
            if (j > static_cast<int>(i))
            {
                later.push_back(j);
            }
        }
        if (!later.empty())
        {
            fprintf (out, "%6lu", i+1);
            for (const auto j : later)
            {
                fprintf (out, " %5d", j+1);
            }
            fprintf (out, "\n");
        }
    }
    fprintf (out, "\n");
}

//! \brief Print the system part of a gromacs topology
static void print_top_system(FILE *out, const char *title)
{
    fprintf(out, "[ %s ]\n", dir2str(d_system));
    fprintf(out, "; Name\n");
    fprintf(out, "%s\n\n", title[0] ? title : "Protein");
}

void print_top_mols(FILE *out, const char *title,
                    int nmol, const t_mols *mols)
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

//! \brief Print the atoms to a gromacs top
static void print_atoms(FILE                           *out,
                        const std::vector<ActAtom>     &atoms,
                        const std::vector<std::string> &residueNames)
{
    fprintf(out, "[ %s ]\n", dir2str(d_atoms));
    fprintf(out, "; %4s %10s %6s %7s%6s %6s %10s %10s %6s %10s %10s\n",
            "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass", "typeB", "chargeB", "massB");

    // The same names as in the PDB file, such that GROMACS can match
    // them: there they have four characters at most
    auto   names = particleNames(atoms);
    for (auto &name : names)
    {
        if (name.size() > 4)
        {
            name.resize(4);
        }
    }
    double qtot  = 0;
    // The residues are numbered from one in each molecule type
    int    res0  = atoms.empty() ? 0 : atoms[0].residueNumber();
    for (size_t i = 0; i < atoms.size(); i++)
    {
        int atomnr = 1+i;
        /* This is true by construction, but static analysers don't know */
        fprintf(out, "%6d %10s %6d%c %5s %6s %6d %12.9g %12.9g",
                atomnr,
                atoms[i].ffType().c_str(),
                atoms[i].residueNumber()-res0+1, ' ',
                residueNames[atoms[i].residueNumber()].c_str(),
                names[i].c_str(),
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

void write_top(MsgHandler                     *msg_handler,
               FILE                           *out,
               const char                     *molname,
               const Topology                 *topology,
               const ForceField               *pd,
               const std::vector<ActAtom>     *atoms,
               const std::vector<std::string> *residueNames)
{
    const auto &myAtoms = atoms ? *atoms : topology->atoms();
    if (myAtoms.empty())
    {
        return;
    }
    fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
    fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
    int nexcl = std::max(0, pd->nexcl());
    fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", nexcl);
    print_atoms(out, myAtoms, residueNames ? *residueNames : topology->residueNames());
    for (auto &fs : pd->forcesConst())
    {
        auto iType = fs.first;
        if (!topology->hasEntry(iType) || isVsite(iType))
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
        else if (InteractionType::POLARIZATION == iType)
        {
            print_polarization(out, topology->entry(iType));
        }
    }
    print_vsites(msg_handler, out, topology);
    print_excl(out, topology);
}

/*! \brief Write a molecule type to a string
 * \param[in] msg_handler  For warnings etc.
 * \param[in] molname      The molecule name
 * \param[in] topology     ACT structure
 * \param[in] pd           The force field
 * \param[in] atoms        The particles
 * \param[in] residueNames The residue names
 * \return the text, empty if there is no temporary file to write to
 */
static std::string moleculeTypeString(MsgHandler                     *msg_handler,
                                      const char                     *molname,
                                      const Topology                 *topology,
                                      const ForceField               *pd,
                                      const std::vector<ActAtom>     &atoms,
                                      const std::vector<std::string> &residueNames)
{
    std::string text;
    FILE       *tmp = std::tmpfile();
    if (nullptr == tmp)
    {
        return text;
    }
    write_top(msg_handler, tmp, molname, topology, pd, &atoms, &residueNames);
    long size = std::ftell(tmp);
    if (size > 0)
    {
        text.resize(size);
        std::rewind(tmp);
        if (std::fread(&text[0], 1, size, tmp) != static_cast<size_t>(size))
        {
            text.clear();
        }
    }
    std::fclose(tmp);
    return text;
}

/*! \brief A text that is the same for compounds that are the same molecule
 *
 * The particle types and charges, and per interaction type the sorted
 * interactions, each with its atoms in an order that does not depend on
 * the direction it was stored in, and its parameters. Of a virtual site
 * the atoms it is made from are sorted and the signs of the parameters
 * are left out, since two sites of the same type on the same atoms, as
 * the lone pairs of water, may be numbered either way. The exclusions
 * are included as well.
 * \param[in] topology The compound
 * \param[in] atoms    The particles with their charges
 * \return the text
 */
static std::string compoundKey(const Topology             *topology,
                               const std::vector<ActAtom> &atoms)
{
    std::string key;
    for (const auto &a : atoms)
    {
        key += gmx::formatString("%s %.9g %.9g|", a.ffType().c_str(), a.charge(), a.mass());
    }
    for (const auto &ee : topology->entries())
    {
        std::vector<std::string> lines;
        for (const auto &entry : ee.second)
        {
            auto ind = entry->atomIndices();
            auto par = entry->params();
            if (isVsite(ee.first))
            {
                // The site last, the atoms sorted
                std::sort(ind.begin(), ind.end() - 1);
                for (auto &p : par)
                {
                    p = std::abs(p);
                }
            }
            else if (ind.size() > 1 && ind.back() < ind.front() &&
                     InteractionType::POLARIZATION != ee.first &&
                     InteractionType::IMPROPER_DIHEDRALS != ee.first)
            {
                std::reverse(ind.begin(), ind.end());
            }
            std::string line;
            for (auto i : ind)
            {
                line += gmx::formatString("%d ", i);
            }
            for (auto p : par)
            {
                line += gmx::formatString("%.9g ", p);
            }
            lines.push_back(line);
        }
        std::sort(lines.begin(), lines.end());
        key += interactionTypeToString(ee.first) + ":";
        for (const auto &l : lines)
        {
            key += l + ";";
        }
        key += "|";
    }
    auto itype = InteractionType::VDW;
    if (topology->hasExclusions(itype))
    {
        for (const auto &ex : topology->exclusions(itype))
        {
            auto sorted = ex;
            std::sort(sorted.begin(), sorted.end());
            for (auto j : sorted)
            {
                key += gmx::formatString("%d ", j);
            }
            key += ";";
        }
    }
    return key;
}

bool write_top_compounds(MsgHandler                     *msg_handler,
                         FILE                           *out,
                         const char                     *title,
                         const std::vector<std::string> &names,
                         const std::vector<Topology>    &compounds,
                         const Topology                 *topology,
                         const ForceField               *pd,
                         bool                            printMolecules)
{
    const auto &allAtoms = topology->atoms();
    size_t      natom    = 0;
    for (const auto &c : compounds)
    {
        natom += c.atoms().size();
    }
    if (compounds.empty() || names.size() != compounds.size() || natom != allAtoms.size())
    {
        return false;
    }
    // Compounds with the same key are the same molecule type: the same
    // particles, charges, interactions and parameters. The molecule
    // type is written for the first of these, with a place holder for
    // its name.
    const std::string        placeHolder("ACT_NAME_OF_THE_MOLECULE_TYPE");
    std::vector<std::string> text(compounds.size());
    std::vector<std::string> key;
    size_t                   start = 0;
    for (size_t f = 0; f < compounds.size(); f++)
    {
        const auto &c = compounds[f];
        // The particles as they are in the system, which has the charges
        const std::vector<ActAtom> atoms(allAtoms.begin() + start, allAtoms.begin() + start + c.atoms().size());
        for (size_t i = 0; i < atoms.size(); i++)
        {
            if (atoms[i].ffType() != c.atoms()[i].ffType())
            {
                return false;
            }
        }
        key.push_back(compoundKey(&c, atoms));
        bool isNew = true;
        for (size_t g = 0; g < f && isNew; g++)
        {
            isNew = key[g] != key[f];
        }
        if (isNew)
        {
            text[f] = moleculeTypeString(msg_handler, placeHolder.c_str(), &c, pd, atoms,
                                         topology->residueNames());
            if (text[f].find(placeHolder) == std::string::npos)
            {
                return false;
            }
        }
        start += atoms.size();
    }
    // The name of the molecule type of each compound: that of an earlier
    // compound with the same key, or a new one
    std::vector<std::string> typeName(compounds.size());
    std::vector<size_t>      unique;
    for (size_t f = 0; f < compounds.size(); f++)
    {
        for (auto u : unique)
        {
            if (key[u] == key[f])
            {
                typeName[f] = typeName[u];
                break;
            }
        }
        if (!typeName[f].empty())
        {
            continue;
        }
        // GROMACS reads a name up to the first white space
        std::string base = names[f].empty() ? std::string("Compound") : names[f];
        for (auto &ch : base)
        {
            if (std::isspace(static_cast<unsigned char>(ch)))
            {
                ch = '_';
            }
        }
        std::string name  = base;
        int         count = 1;
        const auto  taken = [&](const std::string &nm)
        {
            for (auto u : unique)
            {
                if (typeName[u] == nm)
                {
                    return true;
                }
            }
            return false;
        };
        while (taken(name))
        {
            count += 1;
            name   = gmx::formatString("%s_%d", base.c_str(), count);
        }
        typeName[f] = name;
        unique.push_back(f);
    }
    for (auto u : unique)
    {
        auto mt  = text[u];
        auto pos = mt.find(placeHolder);
        mt.replace(pos, placeHolder.size(), typeName[u]);
        fprintf(out, "%s", mt.c_str());
    }
    if (printMolecules)
    {
        print_top_system(out, title);
        fprintf(out, "[ %s ]\n", dir2str(d_molecules));
        fprintf(out, "; %-15s %5s\n", "Compound", "#mols");
        // Compounds of the same type that follow each other on one line
        for (size_t f = 0; f < compounds.size(); )
        {
            size_t n = 1;
            while (f + n < compounds.size() && typeName[f + n] == typeName[f])
            {
                n += 1;
            }
            fprintf(out, "%-15s %5zu\n", typeName[f].c_str(), n);
            f += n;
        }
    }
    return true;
}

//! \return the particle type of GROMACS for that of ACT
static const char *gromacsParticleType(ActParticle apType)
{
    switch (apType)
    {
    case ActParticle::Shell:
        return "S";
    case ActParticle::Vsite:
        return "V";
    default:
        return "A";
    }
}

void print_top_header(MsgHandler                     *msg_handler,
                      FILE                           *fp,
                      const ForceField               *pd,
                      const Topology                 *topology,
                      const std::vector<std::string> &commercials,
                      bool                            bItp)
{
    // What GROMACS does not have of this force field
    std::vector<std::string> missing;
    if (pd->chargeDistribution() != ChargeDistributionType::Point)
    {
        missing.push_back(gmx::formatString("%s, the charges are point charges here",
                                            potentialToString(pd->potential(InteractionType::ELECTROSTATICS)).c_str()));
    }
    if (pd->interactionPresent(InteractionType::VDW) &&
        Potential::LJ12_6 != pd->findForcesConst(InteractionType::VDW).potential())
    {
        missing.push_back(gmx::formatString("%s, there are no van der Waals interactions here",
                                            potentialToString(pd->findForcesConst(InteractionType::VDW).potential()).c_str()));
    }
    for (const auto &fs : pd->forcesConst())
    {
        // Potentials without a counterpart, for interactions in this topology
        const auto pot = fs.second.potential();
        if (topology->hasEntry(fs.first) && !topology->entry(fs.first).empty() &&
            (Potential::HUA_BONDS == pot || Potential::CUBIC_BONDS == pot ||
             Potential::QUARTIC_BONDS == pot))
        {
            missing.push_back(gmx::formatString("%s, GROMACS will not accept this topology",
                                                potentialToString(pot).c_str()));
        }
    }
    for (auto itype : { InteractionType::INDUCTIONCORRECTION, InteractionType::VDWCORRECTION })
    {
        if (pd->interactionPresent(itype) &&
            !pd->findForcesConst(itype).parametersConst().empty())
        {
            missing.push_back(gmx::formatString("the %s", interactionTypeToDescription(itype).c_str()));
        }
    }
    fprintf(fp, ";\n");
    fprintf(fp, "; Topology generated by alexandria gentop.\n");
    if (missing.empty())
    {
        fprintf(fp, "; GROMACS has all the interactions of this force field.\n");
    }
    else
    {
        fprintf(fp, "; WARNING: do not use for simulations! GROMACS does not have\n");
        for (const auto &m : missing)
        {
            fprintf(fp, ";   %s\n", m.c_str());
        }
    }
    if (pd->epsilonr() != 1)
    {
        fprintf(fp, "; The force field has a relative dielectric constant: use epsilon-r = %g\n", pd->epsilonr());
        fprintf(fp, "; in the mdp file.\n");
    }
    fprintf(fp, "; Watch this space for information & commercials.\n");
    for (auto i = commercials.begin(); (i < commercials.end()); ++i)
    {
        fprintf(fp, "; %s\n", i->c_str());
    }
    fprintf(fp, ";\n");
    if (bItp)
    {
        return;
    }
    if (pd->chargeDistribution() != ChargeDistributionType::Point)
    {
        msg_handler->msg(ACTStatus::Warning,
                         gmx::formatString("GROMACS does not support %s, the charges are point charges in the topology",
                                           potentialToString(pd->potential(InteractionType::ELECTROSTATICS)).c_str()));
    }
    // Lennard-Jones with sigma and epsilon, the parameters of the pairs of
    // particle types are given explicitly below, such that any combination
    // rule of ACT is reproduced. The exclusions of ACT are those of the
    // moleculetype and there are no special pairs.
    auto &vdw = pd->findForcesConst(InteractionType::VDW);
    bool  lj  = (Potential::LJ12_6 == vdw.potential());
    if (!lj)
    {
        msg_handler->msg(ACTStatus::Warning,
                         gmx::formatString("GROMACS does not support %s, the topology has no van der Waals interactions",
                                           potentialToString(vdw.potential()).c_str()));
    }
    fprintf(fp, "[ defaults ]\n");
    fprintf(fp, "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n");
    fprintf(fp, "%-8d  %-9d  %-9s  %7g  %7g\n\n", 1, 3, "no", 1.0, 1.0);

    // The sigma and epsilon of a van der Waals type, or of a pair of them
    const auto sigmaEpsilon = [&vdw, lj](const Identifier &id, double *sigma, double *epsilon)
    {
        *sigma   = 0;
        *epsilon = 0;
        if (lj && !id.id().empty() && vdw.parameterExists(id))
        {
            const auto &params = vdw.findParameterMapConst(id);
            *sigma   = params.find("sigma")->second.internalValue();
            *epsilon = params.find("epsilon")->second.internalValue();
        }
    };
    fprintf(fp, "[ atomtypes ]\n");
    fprintf(fp, "; %-6s  %6s  %12s  %10s  %5s  %12s  %12s\n",
            "name", "at.num", "mass", "charge", "ptype", "sigma", "epsilon");
    for (const auto &aType : pd->particleTypesConst())
    {
        double sigma, epsilon;
        sigmaEpsilon(aType.second.interactionTypeToIdentifier(InteractionType::VDW), &sigma, &epsilon);
        fprintf(fp, "%-8s  %6d  %12.6f  %10.4f  %5s  %12.6g  %12.6g\n",
                aType.second.id().id().c_str(), aType.second.atomnumber(),
                aType.second.mass(), 0.0, gromacsParticleType(aType.second.apType()),
                sigma, epsilon);
    }
    fprintf(fp, "\n");
    if (!lj)
    {
        return;
    }
    // The pairs of the particle types in the molecule
    std::set<std::string> types;
    for (const auto &atom : topology->atoms())
    {
        types.insert(atom.ffType());
    }
    fprintf(fp, "[ nonbond_params ]\n");
    fprintf(fp, "; %-6s  %-8s  %4s  %12s  %12s\n", "i", "j", "func", "sigma", "epsilon");
    for (auto ti = types.begin(); ti != types.end(); ++ti)
    {
        for (auto tj = ti; tj != types.end(); ++tj)
        {
            const auto vi = pd->findParticleType(*ti)->interactionTypeToIdentifier(InteractionType::VDW);
            const auto vj = pd->findParticleType(*tj)->interactionTypeToIdentifier(InteractionType::VDW);
            double sigma, epsilon;
            sigmaEpsilon(Identifier({ vi.id(), vj.id() }, { 1 }, CanSwap::Yes), &sigma, &epsilon);
            fprintf(fp, "%-8s  %-8s  %4d  %12.6g  %12.6g\n", ti->c_str(), tj->c_str(), 1, sigma, epsilon);
        }
    }
    fprintf(fp, "\n");
}

} // namespace alexandria
