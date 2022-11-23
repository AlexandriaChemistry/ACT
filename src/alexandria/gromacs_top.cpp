#include "gromacs_top.h"

#include <map>

#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/path.h"

#include "act/poldata/poldata.h"
#include "topology.h"

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

static void print_bondeds(FILE                               *out,
                          directive                           d,
                          int                                 ftype,
                          const std::vector<TopologyEntry *> &entries)
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
        for (const double &j : entry->params())
        {
            fprintf(out, "  %10g", j);
        }
        fprintf(out, "; %s\n", entry->id().id().c_str());
    }
    fprintf(out, "\n");
}

static void print_excl(FILE *out, int natoms, t_excls excls[])
{
    int         i;
    bool        have_excl;
    int         j;

    have_excl = FALSE;
    for (i = 0; i < natoms && !have_excl; i++)
    {
        have_excl = (excls[i].nr > 0);
    }

    if (have_excl)
    {
        fprintf (out, "[ %s ]\n", dir2str(d_exclusions));
        fprintf (out, "; %4s    %s\n", "i", "excluded from i");
        for (i = 0; i < natoms; i++)
        {
            if (excls[i].nr > 0)
            {
                fprintf (out, "%6d ", i+1);
                for (j = 0; j < excls[i].nr; j++)
                {
                    fprintf (out, " %5d", excls[i].e[j]+1);
                }
                fprintf (out, "\n");
            }
        }
        fprintf (out, "\n");
    }
}

static void print_top_system(FILE *out, const char *title)
{
    fprintf(out, "[ %s ]\n", dir2str(d_system));
    fprintf(out, "; Name\n");
    fprintf(out, "%s\n\n", title[0] ? title : "Protein");
}

static void print_top_water(FILE *out, const char *ffdir, const char *water)
{
    const char *p;
 
    fprintf(out, "; Include water topology\n");

    p = strrchr(ffdir, '/');
    p = (ffdir[0] == '.' || p == nullptr) ? ffdir : p+1;
    fprintf(out, "#include \"%s/%s.itp\"\n", p, water);

    fprintf(out, "\n");
    fprintf(out, "#ifdef POSRES_WATER\n");
    fprintf(out, "; Position restraint for each water oxygen\n");
    fprintf(out, "[ position_restraints ]\n");
    fprintf(out, ";%3s %5s %9s %10s %10s\n", "i", "funct", "fcx", "fcy", "fcz");
    fprintf(out, "%4d %4d %10g %10g %10g\n", 1, 1, 1000.0, 1000.0, 1000.0);
    fprintf(out, "#endif\n");
    fprintf(out, "\n");

    std::string buf = gmx::formatString("%s/ions.itp", p);

    if (fflib_fexist(buf.c_str()))
    {
        fprintf(out, "; Include topology for ions\n");
        fprintf(out, "#include \"%s\"\n", buf.c_str());
        fprintf(out, "\n");
    }
}

void print_top_mols(FILE *out,
                    const char *title, const char *ffdir, const char *water,
                    int nincl, char **incls, int nmol, t_mols *mols)
{

    if (nincl > 0)
    {
        fprintf(out, "; Include chain topologies\n");
        for (int i = 0; i < nincl; i++)
        {
            fprintf(out, "#include \"%s\"\n", gmx::Path::getFilename(incls[i]).c_str());
        }
        fprintf(out, "\n");
    }

    if (water)
    {
        print_top_water(out, ffdir, water);
    }
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

static void print_atoms(FILE                       *out,
                        const std::vector<ActAtom> &atoms,
                        const std::string          &residueName)
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
                atoms[i].residueNumber(), ' ',
                residueName.c_str(),
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


void write_top(FILE            *out,
               char            *molname,
               t_atoms         *at,
               const Topology  *topology,
               t_excls          excls[],
               struct gpp_atomtype *atype,
               const Poldata   *pd)
{
    std::map<int, directive> toPrint = {
        { F_CONSTR,       d_constraints },
        { F_CONSTRNC,     d_constraints },
        { F_LJ14,         d_pairs },
        { F_CMAP,         d_cmap },
        { F_POLARIZATION, d_polarization },
        { F_THOLE_POL,    d_thole_polarization },
        { F_VSITE2,       d_vsites2 },
        { F_VSITE3,       d_vsites3 },
        { F_VSITE3FD,     d_vsites3 },
        { F_VSITE3FAD,    d_vsites3 },
        { F_VSITE3OUT,    d_vsites3 },
        { F_VSITE4FD,     d_vsites4 },
        { F_VSITE4FDN,    d_vsites4 }
    };
    auto myAtoms = topology->atoms();
    if (!myAtoms.empty() && atype)
    {
        std::vector<int> cgnr;
        for(size_t i = 0; i < myAtoms.size(); i++)
        {
            cgnr.push_back(i+1);
        }
        fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", pd->getNexcl());
        print_atoms(out, myAtoms, topology->name());
        for (auto &fs : pd->forcesConst())
        {
            auto iType = fs.first;
            if (!topology->hasEntry(iType))
            {
                continue;
            }
            auto fType = fs.second.fType();
            if (InteractionType::BONDS == fs.first)
            {
                print_bondeds(out, d_bonds, fType, topology->entry(iType));
            }
            else if (InteractionType::ANGLES == iType || InteractionType::LINEAR_ANGLES == iType)
            {
                print_bondeds(out, d_angles, fType, topology->entry(iType));
            }
            else if (InteractionType::PROPER_DIHEDRALS == iType || InteractionType::IMPROPER_DIHEDRALS == iType)
            {
                print_bondeds(out, d_dihedrals, fType, topology->entry(iType));
            }
            else if (toPrint.end() != toPrint.find(fType))
            {
                print_bondeds(out, toPrint.find(fType)->second, fType, topology->entry(iType));
            }
        }
        print_excl(out, at->nr, excls);
    }
}

void print_top_header(FILE                    *fp,
                      const Poldata           *pd,
                      bool                     bPol,
                      const std::vector<std::string> &commercials,
                      bool                     bItp)
{
    std::string   gt_old, gt_type;
    auto qt          = pd->findForcesConst(InteractionType::COULOMB);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));

    fprintf(fp, ";\n");
    fprintf(fp, "; Topology generated by alexandria gentop.\n");
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
        auto ftype = pd->findForcesConst(InteractionType::VDW).fType();
        std::string ff;
        if (ftype == F_LJ)
        {
            ff.assign("LJ");
        }
        auto combRule = pd->findForcesConst(InteractionType::VDW).optionValue("combination_rule");
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
            gt_type    = aType.id().id();
            std::string bType;
            if (aType.hasInteractionType(InteractionType::BONDS))
            {
                bType = aType.interactionTypeToIdentifier(InteractionType::BONDS).id();
            }
            if ((0 ==  gt_old.size()) || (gt_old.compare(gt_type) != 0))
            {
                auto vdwtype = aType.interactionTypeToIdentifier(InteractionType::VDW);
                double sigma = 0, epsilon = 0, gamma = 0;
                if (!vdwtype.id().empty())
                {
                    auto myvdw = vdw.findParametersConst(vdwtype);
                    sigma      = myvdw["sigma"].value();
                    epsilon    = myvdw["epsilon"].value();
                    gamma      = myvdw["gamma"].value();
                }
                fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f %s %g %g %g\n",
                        gt_type.c_str(), 
                        !bType.empty() ? bType.c_str() : gt_type.c_str(), 
                        aType.atomnumber(), 
                        aType.mass(), 0.0,
                        ptype_str[aType.gmxParticleType()],
                        sigma, epsilon, gamma);
                if (false && bPol)
                {
                    auto sgt_type= aType.interactionTypeToIdentifier(InteractionType::POLARIZATION);
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
        if (iChargeType != ChargeType::Point)
        {
            auto eem = pd->findForcesConst(InteractionType::COULOMB);
            fprintf(fp, "[ distributed_charges ]\n");
            for (const auto &atype : pd->particleTypesConst())
            {
                auto ztype     = atype.interactionTypeToIdentifier(InteractionType::COULOMB);
                auto eep       = eem.findParametersConst(ztype);
                if (ChargeType::Slater == iChargeType)
                {
                    fprintf(fp, "%-7s  2  %d  %g\n", atype.id().id().c_str(),
                            atype.row(), eep["zeta"].value());
                }
                else if (ChargeType::Gaussian == iChargeType)
                {
                    fprintf(fp, "%-7s  1  %g\n", atype.id().id().c_str(),
                            eep["zeta"].value());
                }
            }
            fprintf(fp, "\n");
        }
    }
}

void excls_to_blocka(int natom, t_excls excls_[], t_blocka *blocka)
{
    int i, j, k, nra;

    if (blocka->nr < natom)
    {
        srenew(blocka->index, natom+1);
        for (int i = blocka->nr; i < natom+1; i++)
        {
            blocka->index[i] = 0;
        }
    }
    nra = 0;
    for (i = 0; (i < natom); i++)
    {
        nra += excls_[i].nr;
    }
    snew(blocka->a, nra+1);
    nra = 0;
    for (i = j = 0; (i < natom); i++)
    {
        blocka->index[i] = nra;
        for (k = 0; (k < excls_[i].nr); k++)
        {
            blocka->a[j++] = excls_[i].e[k];
        }
        nra += excls_[i].nr;
    }
    blocka->index[natom] = nra;
    blocka->nr           = natom;
    blocka->nra          = nra;
}

} // namespace alexandria
