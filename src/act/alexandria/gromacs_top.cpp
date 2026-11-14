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

#include "gromacs_top.h"

#include <cmath>
#include <map>
#include <set>
#include <vector>

#include "act/alexandria/pdbwriter.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/topology/topology.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

//! \return topology function index for bonds, angles etc.
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
    int subtype = (F_FOURDIHS == ftype) ? 3 : get_subtype(d, ftype);
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
            // The angles are in radians in ACT and in degrees in GROMACS
            fprintf(out, "  %10g  %10g ; %s\n",
                    params[angleANGLE]*RAD2DEG,
                    params[angleKT],
                    entry->id().id().c_str());
            break;
        case F_UREY_BRADLEY:
            fprintf(out, "  %10g  %10g  %10g  %10g ; %s\n",
                    params[ubANGLE]*RAD2DEG,
                    params[ubKT],
                    params[ubR13],
                    params[ubKUB],
                    entry->id().id().c_str());
            break;
        case F_LINEAR_ANGLES:
            fprintf(out, "  %10g  %10g ; %s\n",
                    params[linangA],
                    params[linangKLIN],
                    entry->id().id().c_str());
            break;
        case F_PDIHS:
            fprintf(out, "  %10g  %10g  %d ; %s\n",
                    params[pdihANGLE]*RAD2DEG,
                    params[pdihKP],
                    static_cast<int>(std::lround(params[pdihMULT])),
                    entry->id().id().c_str());
            break;
        case F_FOURDIHS:
            for (int k = fdihC0; k < fdihNR; k++)
            {
                fprintf(out, "  %10g", k < static_cast<int>(params.size()) ? params[k] : 0.0);
            }
            fprintf(out, " ; %s\n", entry->id().id().c_str());
            break;
        case F_IDIHS:
            fprintf(out, "  0.0 %10g ; %s\n",
                    params[idihKPHI],
                    entry->id().id().c_str());
            break;
        default: // writes warning
            fprintf(out, " ; %s, unsupported function %s\n", entry->id().id().c_str(),
                    interaction_function[ftype].name);
            msg_handler->msg(ACTStatus::Warning,
                             gmx::formatString("Unsupported function type %s for writing GROMACS topology",
                                               interaction_function[ftype].name));
        }
    }
    fprintf(out, "\n");
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

    // The same names as in the PDB file, such that GROMACS can match them
    auto   names = particleNames(atoms);
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

void write_top(MsgHandler       *msg_handler,
               FILE             *out,
               const char       *molname,
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
    const auto &myAtoms = topology->atoms();
    if (!myAtoms.empty())
    {
        std::vector<int> cgnr;
        for(size_t i = 0; i < myAtoms.size(); i++)
        {
            cgnr.push_back(i+1);
        }
        fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        int nexcl = std::max(0, pd->nexcl());
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
    fprintf(fp, ";\n");
    fprintf(fp, "; Topology generated by alexandria gentop.\n");
    fprintf(fp, "; WARNING: do not use for simulations!\n");
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
