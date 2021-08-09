/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
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

#include "molselect.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#include <strings.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/textreader.h"

#include "alex_modules.h"
#include "composition.h"
#include "molgen.h"
#include "molprop.h"
#include "molprop_xml.h"
#include "mymol.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "stringutil.h"

std::map<iMolSelect, const char *> MolSelect_Names = {
    { iMolSelect::Train,   "Train"   },
    { iMolSelect::Test,    "Test"    },
    { iMolSelect::Unknown, "Unknown" }
};

const char *iMolSelectName(iMolSelect ims)
{
    return MolSelect_Names[ims];
}

bool name2molselect(const std::string &name, iMolSelect *ims)
{
    for (auto iter = MolSelect_Names.begin(); iter != MolSelect_Names.end(); ++iter)
    {
        if (name.compare(iter->second) == 0)
        {
            *ims = iter->first;
            return true;
        }
    }
    return false;
}

namespace alexandria
{

static void sample_molecules(FILE                           *log,
                             FILE                           *fp,
                             std::vector<alexandria::MyMol> &mols,
                             alexandria::Poldata            *pd,
                             int                             mindata,
                             int                             maxatempt)
{

    int nmol      = 0;
    int atempt    = 0;

    CompositionSpecs                      cs;
    std::random_device                    rd;
    std::mt19937                          gen(rd());
    std::uniform_int_distribution<>       int_uniform(0, mols.size()-1);
    std::vector<alexandria::MyMol>        sample;

    const char  *alexandria = cs.searchCS(alexandria::iCalexandria)->name();
    for (const auto &atp : pd->particleTypesConst())
    {
        if (atp.element().compare("H") == 0)
        {
            nmol   = 0;
            atempt = 0;
            do
            {               
                auto  mol   = mols[int_uniform(gen)];
                auto  mci   = mol.SearchMolecularComposition(alexandria);
                if (mci != mol.EndMolecularComposition())
                {
                    for (auto &ani : mci->atomNumConst())
                    {
                        if (atp.id().id() == ani.getAtom())
                        {
                            if ((std::find(sample.begin(), sample.end(), mol) == sample.end()))
                            {
                                sample.push_back(mol);
                                nmol++;
                                break;
                            }
                        }
                    }
                }
                atempt++;
            }
            while (nmol < mindata && atempt < maxatempt);
            if (atempt >= maxatempt)
            {
                fprintf(log, "Picked only %d out of %d molecules required for %s atom type after %d attempts\n", 
                        nmol, mindata, atp.id().id().c_str(), atempt);
            }
        }
    }
    
    std::uniform_real_distribution<> real_uniform(0, 1);
    for (const auto &mol : sample)
    {
        if (real_uniform(gen) >= 0.5)
        {
            fprintf(fp, "%s|Train\n", mol.getMolname().c_str());
        }
        else
        {
            fprintf(fp, "%s|Test\n", mol.getMolname().c_str());
        }
    }
}

void MolSelect::read(const char *fn)
{
    gmx::TextReader tr(fn);
    std::string     tmp;
    int             index = 0;

    while (tr.readLine(&tmp))
    {
        while (!tmp.empty() && tmp[tmp.length()-1] == '\n')
        {
            tmp.erase(tmp.length()-1);
        }
        auto ptr = split(tmp, '|');
        if ((ptr.size() == 2) && (ptr[0].length() > 1)
            && (ptr[1].length() > 1))
        {
            iMolSelect status = iMolSelect::Unknown;
            for (auto iter = MolSelect_Names.begin(); iter != MolSelect_Names.end(); ++iter)
            {
                if (strcasecmp(iter->second, ptr[1].c_str()) == 0)
                {
                    status = iter->first;
                    break;
                }
            }
            if (status == iMolSelect::Unknown)
            {
                fprintf(stderr, "Unknown status '%s' for molecule %s on line %d in file %s\n",
                        ptr[1].c_str(), ptr[0].c_str(), index, fn);
            }
            ims_.push_back(IMolSelect(ptr[0], status, index++));
        }
        else
        {
            fprintf(stderr, "Invalid line '%s' in selection file\n",
                    tmp.c_str());
        }
    }
}

iMolSelect MolSelect::status(const std::string &iupac) const
{
    auto imi = std::find_if(ims_.begin(), ims_.end(),
                            [iupac](IMolSelect const &i)
                            {
                                return i.iupac().compare(iupac) == 0;
                            });

    if (imi != ims_.end())
    {
        return imi->status();
    }

    return iMolSelect::Unknown;
}

int MolSelect::index(const std::string &iupac) const
{
    auto imi = std::find_if(ims_.begin(), ims_.end(),
                            [iupac](IMolSelect const &i)
                            {
                                return i.iupac().compare(iupac) == 0;
                            });

    if (imi != ims_.end())
    {
        return imi->index();
    }

    return -1;
}

static void printAtomtypeStatistics(FILE                                 *fp,
                                    const alexandria::Poldata            *pd,
                                    const std::vector<alexandria::MyMol> &mymol)
{
    struct NN
    {
        std::string name;
        int         count;
    };
    std::vector<NN> nn;
    for (const auto &atype : pd->particleTypesConst())
    {
        struct NN n;
        n.name   = atype.id().id();
        n.count  = 0;
        nn.push_back(n);
    }
    for (const auto &mol : mymol)
    {
        auto atype = mol.gromppAtomtype();
        int ntypes = get_atomtype_ntypes(*atype);
        for (int i = 0; i < ntypes; i++)
        {
            char *tp = get_atomtype_name(i, *atype);
            for (auto &n : nn)
            {
                if (n.name.compare(tp) == 0)
                {
                    n.count += 1;
                    break;
                }
            }
        }
    }
    fprintf(fp, "Atomtype     Count\n");
    for (const auto &n : nn)
    {
        fprintf(fp, "%-8s  %8d\n", n.name.c_str(), n.count);
    }
}

}

int alex_molselect(int argc, char *argv[])
{
    const char *desc[] = {
        "molselect generates random samples from molprop database"
    };

    t_filenm    fnm[] =
    {
        { efDAT, "-f",    "allmols",   ffREAD },
        { efDAT, "-d",    "gentop",    ffREAD },
        { efDAT, "-o",    "selection", ffWRITE },
        { efLOG, "-g",    "molselect", ffWRITE },
        { efDAT, "-sel",  "molselect", ffREAD  },
    };

    const  int                  NFILE     = asize(fnm);

    static int                  nsample   = 1;
    static int                  maxatempt = 5000;
    static char                *opt_elem  = nullptr;
    static gmx_bool             bZero     = false;
    
    static const char          *select_types[]   = {nullptr, "Train", "Test", "Ignore", "Unknown", nullptr};
    
    t_pargs                     pa[]      =
    {
        { "-nsample",   FALSE, etINT, {&nsample},
          "Number of replicas." },
        { "-zero_dipole",    FALSE, etBOOL, {&bZero},
          "Take into account molecules with zero dipoles." },
        { "-maxatempt", FALSE, etINT, {&maxatempt},
          "Maximum number of atempts to sample mindata molecules per atom types." },
        { "-select", FALSE, etENUM, {select_types},
          "Select type for making the dataset for training or testing." },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of atom types to select molecules. If this variable is not set, all elements will be used." }
    };

    gmx_output_env_t       *oenv;
    alexandria::MolGen      mgn;
    alexandria::MolSelect   gms;
    time_t                  my_t;
    FILE                   *fp;

    std::vector<t_pargs>    pargs;
    for (size_t i = 0; i < sizeof(pa)/sizeof(pa[0]); i++)
    {
        pargs.push_back(pa[i]);
    }
    mgn.addOptions(&pargs, etuneNone);
    if (!parse_common_args(&argc, 
                           argv, 
                           PCA_CAN_VIEW, 
                           NFILE, fnm,
                           pargs.size(), 
                           pargs.data(),
                           asize(desc), 
                           desc, 
                           0, 
                           nullptr, 
                           &oenv))
    {
        return 0;
    }
    mgn.optionsFinished();

    fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

    time(&my_t);
    fprintf(fp, "# This file was created %s", ctime(&my_t));
    fprintf(fp, "# alexandria is part of GROMACS:\n#\n");
    fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());

    gms.read(opt2fn_null("-sel", NFILE, fnm));

    iMolSelect select_type;
    if (!name2molselect(select_types[0], &select_type))
    {
        gmx_fatal(FARGS, "No such selection type %s", select_types[0]);
    }
    mgn.Read(fp ? fp : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             bZero,
             gms, 
             false, 
             false, 
             nullptr,
             select_type);

    printAtomtypeStatistics(fp, mgn.poldata(), mgn.molset());

    for (int i = 0; i < nsample; i++)
    {
        char  buf[STRLEN];
        
        sprintf(buf, "%s_%d.dat", fnm[2].fn, i);
        
        fprintf(fp, "Sample file: %s\n\n", buf);
        
        FILE *dat = gmx_ffopen(buf, "w");
        
        sample_molecules(fp,
                         dat, 
                         mgn.molset(), 
                         mgn.poldata(),
                         mgn.mindata(), 
                         maxatempt);
        gmx_ffclose(dat);
    }
    gmx_ffclose(fp);
    //done_filenms(NFILE, fnm);

    return 0;
}
