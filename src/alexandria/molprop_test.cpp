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

#include <stdio.h>
#include <stdlib.h>

#include "gromacs/commandline/pargs.h"

#include "alex_modules.h"
#include "molprop.h"
#include "molprop_xml.h"
#include "mymol.h"
#include "poldata_xml.h"

static void dump_molecules(const char *logfilename,
                           const char *gentop,
                           const std::vector<alexandria::MolProp> &mpt)
{
    const std::string method("B3LYP");
    const std::string basis("aug-cc-pVTZ");
    alexandria::Poldata pd;
    try
    {
        alexandria::readPoldata(gentop, &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    std::map<const std::string, int> atomTypeCount;
    std::map<const std::string, int> bccTypeCount;
    FILE *fp = fopen(logfilename, "w");
    for (const auto &mp : mpt)
    {
        alexandria::MyMol mymol;
        mymol.Merge(&mp);
        auto imm = mymol.GenerateTopology(&pd,
                                          method,
                                          basis,
                                          nullptr,
                                          false,
                                          false,
                                          false,
                                          missingParameters::Error,
                                          nullptr);
        if (immStatus::OK != imm)
        {
            fprintf(fp, "Failed to generate topology for %s. Outcome: %s\n",
                    mymol.getMolname().c_str(), immsg(imm));
        }
        else
        {
            fprintf(fp, "Molecule: %s Formula: %s q: %d Mult: %d\n", mymol.getMolname().c_str(),
                    mymol.formula().c_str(), mymol.totalCharge(), mymol.getMultiplicity());
            // Atoms!
            auto &atoms = mymol.atomsConst();
            std::vector<Identifier> atomId;
            auto ztype = InteractionType::CHARGEDISTRIBUTION;
            for (int i = 0; i < atoms.nr; i++)
            {
                const char *atype = *atoms.atomtype[i];
                fprintf(fp, "atom: %2d  %5s  %5s", i+1, *atoms.atomname[i], atype);
                Identifier pid({ atype }, CanSwap::Yes);
                atomId.push_back(pid);
                if (pd.hasParticleType(pid))
                {
                    auto pIter = pd.findParticleType(pid);
                    if (pIter->hasInteractionType(ztype))
                    {
                        auto zid = pIter->interactionTypeToIdentifier(ztype);
                        fprintf(fp, "  %s", zid.id().c_str());
                        auto atypeMap = atomTypeCount.find(zid.id());
                        if (atypeMap == atomTypeCount.end())
                        {
                            atomTypeCount.insert(std::pair<const std::string, int>(zid.id(), 1));
                        }
                        else
                        {
                            atypeMap->second += 1;
                        }
                    }
                }
                fprintf(fp, "\n");
            }
            // Bonds!
            auto bctype = InteractionType::BONDCORRECTIONS;
            for (const auto &b : mp.bondsConst())
            {
                int ai = b.getAi()-1;
                int aj = b.getAj()-1;
                fprintf(fp, "bcc: %3d  %3d  %5g", ai+1, aj+1, b.getBondOrder());
                if (pd.hasParticleType(atomId[ai]) && pd.hasParticleType(atomId[aj]))
                {
                    auto pidI = pd.findParticleType(atomId[ai]);
                    auto pidJ = pd.findParticleType(atomId[aj]);
                    if (pidI->hasInteractionType(ztype) && pidJ->hasInteractionType(ztype))
                    {
                        auto zidI = pidI->interactionTypeToIdentifier(ztype);
                        auto zidJ = pidJ->interactionTypeToIdentifier(ztype);
                        Identifier mybond({ zidI.id(), zidJ.id()}, b.getBondOrder(), CanSwap::No);
                        auto btypeMap   = bccTypeCount.find(mybond.id());
                        bool bondExists = false;
                        auto fs         = pd.findForcesConst(bctype);
                        if (fs.parameterExists(mybond))
                        {
                            fprintf(fp, "  %s", mybond.id().c_str());
                            bondExists = true;
                        }
                        else
                        {
                            Identifier mybond2({ zidJ.id(), zidI.id()}, b.getBondOrder(), CanSwap::No);
                            mybond = mybond2;
                            btypeMap   = bccTypeCount.find(mybond.id());
                            auto fs = pd.findForcesConst(bctype);
                            if (fs.parameterExists(mybond))
                            {
                                fprintf(fp, "  %s", mybond.id().c_str());
                                bondExists = true;
                            }
                        }
                        if (bondExists)
                        {
                            if (btypeMap == bccTypeCount.end())
                            {
                                bccTypeCount.insert(std::pair<const std::string, int>(mybond.id(), 1));
                            }
                            else
                            {
                                btypeMap->second += 1;
                            }
                        }
                    }
                }
                fprintf(fp, "\n");
            }
        }
    }
    fprintf(fp, "Statistics\n");
    for(auto &atc : atomTypeCount)
    {
        fprintf(fp, "atom: %-6s  %5d\n", atc.first.c_str(), atc.second);
    }
    for(auto &bcc : bccTypeCount)
    {
        fprintf(fp, "bcc: %-12s  %5d\n", bcc.first.c_str(), bcc.second);
    }
    fclose(fp);
}

int alex_molprop_test(int argc, char*argv[])
{
    static const char               *desc[] = {
        "molprop_test reads a molprop file and writes a new one.",
        "If both a force field file and a log file name are gived",
        "a dump of atom types and bond types is done."
    };
    gmx_output_env_t                *oenv;
    std::vector<alexandria::MolProp> mpt;
    t_filenm                         fnm[] = {
        { efDAT, "-f", "molin",     ffREAD },
        { efDAT, "-o", "molout",    ffWRITE },
        { efDAT, "-d", "gentop",    ffOPTRD },
        { efLOG, "-g", "molecules", ffOPTWR }
    };
#define NFILE sizeof(fnm)/sizeof(fnm[0])

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, 0, nullptr,
                           1, desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    MolPropRead(opt2fn("-f", NFILE, fnm), &mpt);
    printf("Read %d molecules from %s\n", (int)mpt.size(), opt2fn("-f", NFILE, fnm));
    MolPropWrite(opt2fn("-o", NFILE, fnm), mpt, 1);
    if (opt2bSet("-g", NFILE, fnm) && opt2bSet("-d", NFILE, fnm))
    {
        dump_molecules(opt2fn("-g", NFILE, fnm),
                       opt2fn("-d", NFILE, fnm),
                       mpt);
    }

    return 0;
}
