/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#include <cstdio>
#include <cstdlib>

#include <map>
#include <string>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/strconvert.h"

#include "alexandria/alex_modules.h"
#include "alexandria/fill_inputrec.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "alexandria/mymol.h"
#include "act/poldata/poldata_xml.h"
#include "act/utility/units.h"

namespace alexandria
{

typedef std::map<const std::string, int> stringCount;

static void dump_molecule(FILE              *fp,
                          const std::string &method,
                          const std::string &basis,
                          stringCount       *atomTypeCount,
                          stringCount       *bccTypeCount,
                          const Poldata     &pd,
                          const MolProp     &mp,
                          t_inputrec        *inputrec)
{
    alexandria::MyMol mymol;
    mymol.Merge(&mp);
    mymol.setInputrec(inputrec);
    auto imm = mymol.GenerateTopology(fp,
                                      &pd,
                                      method,
                                      basis,
                                      missingParameters::Error,
                                      nullptr,
                                      false);
    if (immStatus::OK != imm)
    {
        fprintf(fp, "Failed to generate topology for %s. Outcome: %s\n",
                mymol.getMolname().c_str(), immsg(imm));
    }
    else
    {
        std::map<MolPropObservable, iqmType> iqm = {
            { MolPropObservable::DIPOLE, iqmType::Both },
            { MolPropObservable::QUADRUPOLE, iqmType::Both },
            { MolPropObservable::POLARIZABILITY, iqmType::Both },
        };
       
        fprintf(fp, "Molecule: %s Formula: %s q: %d Mult: %d\n", mymol.getMolname().c_str(),
                mymol.formula().c_str(), mymol.totalCharge(), mymol.getMultiplicity());
        mymol.getExpProps(iqm, true, "", "", &pd);
        mymol.Dump(fp);
        // Atoms!
        auto &atoms = mymol.atomsConst();
        std::vector<Identifier> atomId;
        auto ztype = InteractionType::CHARGEDISTRIBUTION;
        for (int i = 0; i < atoms.nr; i++)
        {
            const char *atype = *atoms.atomtype[i];
            fprintf(fp, "atom: %2d  %5s  %5s", i+1, *atoms.atomname[i], atype);
            Identifier pid(atype);
            atomId.push_back(pid);
            if (pd.hasParticleType(pid))
            {
                auto pIter = pd.findParticleType(pid);
                if (pIter->hasInteractionType(ztype))
                {
                    auto zid = pIter->interactionTypeToIdentifier(ztype);
                    fprintf(fp, "  %s", zid.id().c_str());
                    auto atypeMap = atomTypeCount->find(zid.id());
                    if (atypeMap == atomTypeCount->end())
                    {
                        atomTypeCount->insert(std::pair<const std::string, int>(zid.id(), 1));
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
            int ai = b.aI();
            int aj = b.aJ();
            fprintf(fp, "bcc: %3d  %3d  %5g", ai+1, aj+1, b.bondOrder());
            if (pd.hasParticleType(atomId[ai]) && pd.hasParticleType(atomId[aj]))
            {
                auto pidI = pd.findParticleType(atomId[ai]);
                auto pidJ = pd.findParticleType(atomId[aj]);
                if (pidI->hasInteractionType(ztype) && pidJ->hasInteractionType(ztype))
                {
                    auto zidI = pidI->interactionTypeToIdentifier(ztype);
                    auto zidJ = pidJ->interactionTypeToIdentifier(ztype);
                    Identifier mybond({ zidI.id(), zidJ.id()}, { b.bondOrder() }, CanSwap::No);
                    auto btypeMap   = bccTypeCount->find(mybond.id());
                    bool bondExists = false;
                    auto fs         = pd.findForcesConst(bctype);
                    if (fs.parameterExists(mybond))
                    {
                        fprintf(fp, "  %s", mybond.id().c_str());
                        bondExists = true;
                    }
                    else
                    {
                        Identifier mybond2({ zidJ.id(), zidI.id()}, { b.bondOrder() }, CanSwap::No);
                        mybond = mybond2;
                        btypeMap   = bccTypeCount->find(mybond.id());
                        auto fs = pd.findForcesConst(bctype);
                        if (fs.parameterExists(mybond))
                        {
                            fprintf(fp, "  %s", mybond.id().c_str());
                            bondExists = true;
                        }
                    }
                    if (bondExists)
                    {
                        if (btypeMap == bccTypeCount->end())
                        {
                            bccTypeCount->insert(std::pair<const std::string, int>(mybond.id(), 1));
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

int molprop_check(int argc, char*argv[])
{
    static const char               *desc[] = {
        "molprop_check checks calculations for missing hydrogens",
        "and inconsistent dipoles. It also tries to make a topology",
        "and reports errors doing this. Output is to a file."
    };
    t_filenm                         fnm[] =
    {
        { efXML, "-d",  "gentop",  ffREAD },
        { efXML, "-f",  "allmols",  ffREAD },
        { efLOG, "-g",  "molprop_check", ffWRITE }
    };
    int                              NFILE   = (sizeof(fnm)/sizeof(fnm[0]));
    const char *lot                          = "B3LYP/aug-cc-pVTZ";

    t_pargs pa[] = {
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges. Multiple levels can be specified which will be used in the order given, e.g.  B3LYP/aug-cc-pVTZ:HF/6-311G**" }
    };
std::vector<alexandria::MolProp> mp;
    gmx_output_env_t                *oenv;
    
    int  npa   = (sizeof(pa)/sizeof(pa[0]));
    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           npa, pa,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }
    MolPropRead(opt2fn("-f", NFILE, fnm), &mp);
    
    alexandria::Poldata pd;
    try
    {
        alexandria::readPoldata(opt2fn("-d", NFILE, fnm), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    std::string method, basis;
    splitLot(lot, &method, &basis);
    auto inputrec  = new t_inputrec();
    fill_inputrec(inputrec);

    stringCount atomTypeCount;
    stringCount bccTypeCount;

    FILE *mylog = gmx_fio_fopen(opt2fn("-g", NFILE, fnm), "w");
    fprintf(mylog, "Force field file %s\n", opt2fn("-d", NFILE, fnm));
    fprintf(mylog, "Molprop file     %s\n", opt2fn("-f", NFILE, fnm));
    for (auto &m : mp)
    {
        typedef struct
        {
            std::string name;
            rvec        mu;
        } name_mu;
        
        std::vector<name_mu> mus;
        for (auto &ci : m.experimentConst())
        {
            int nH = 0, nC = 0;
            for (auto &cai : ci.calcAtomConst())
            {
                std::string name = cai.getName();
                if (name.compare("H") == 0)
                {
                    nH++;
                }
                else if (name.compare("C") == 0)
                {
                    nC++;
                }
            }
            if (nC > 0 && nH == 0)
            {
                fprintf(mylog, "%s #C %d #H %d\n",
                        ci.getDatafile().c_str(), 
                        nC, nH);
            }
            double T = 0;
            auto gp = m.findProperty(MolPropObservable::DIPOLE, iqmType::QM, T, "", "", "");
            if (gp)
            {
                std::vector<double> mu = gp->getVector();
                name_mu nmu = { ci.getDatafile(), { mu[XX], mu[YY], mu[ZZ] } };
                mus.push_back(nmu);
            }
        
            auto Xcalc = ci.getCoordinates();
            auto Esp   = ci.electrostaticPotentialConst();
            if (Esp.size() >= Xcalc.size() && Xcalc.size() > 1)
            {
                double msd = 0;
                auto xunit = Esp[0].getXYZunit();
                double fac = convertToGromacs(1.0, xunit);
                for(size_t i = 0; i < Xcalc.size(); i++)
                {
                    msd += (gmx::square(Xcalc[i][XX]-fac*Esp[i].getX())+
                            gmx::square(Xcalc[i][YY]-fac*Esp[i].getY())+
                            gmx::square(Xcalc[i][ZZ]-fac*Esp[i].getZ()));
                }
                double rmsd = std::sqrt(msd/Xcalc.size());
                if (rmsd != 0)
                {
                    fprintf(mylog, "%s RMSD coordinates between ESP and QM %g\n",
                            m.getMolname().c_str(), rmsd);
                }
                if (rmsd > 1e-3)
                {
                    for(size_t i = 0; i < Xcalc.size(); i++)
                    {
                        fprintf(mylog, "%2d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                                static_cast<int>(i+1),
                                Xcalc[i][XX], Xcalc[i][YY], Xcalc[i][ZZ],
                                fac*Esp[i].getX(), fac*Esp[i].getY(),
                                fac*Esp[i].getZ());
                    }
                }
            }
        }
        // Check dipoles
        if (debug)
        {
            for(const auto &mi : mus)
            {
                fprintf(debug, "%s %s %.2f %.2f %.2f\n", m.getMolname().c_str(),
                        mi.name.c_str(), mi.mu[XX], mi.mu[YY], mi.mu[ZZ]);
            }
        }

        dump_molecule(mylog, method, basis, &atomTypeCount, &bccTypeCount, pd, m, inputrec);
    }
    fprintf(mylog, "Statistics\n");
    for(auto &atc : atomTypeCount)
    {
        fprintf(mylog, "atom: %-6s  %5d\n", atc.first.c_str(), atc.second);
    }
    for(auto &bcc : bccTypeCount)
    {
        fprintf(mylog, "bcc: %-12s  %5d\n", bcc.first.c_str(), bcc.second);
    }
    fclose(mylog);
    return 0;
}

} // namespace alexandria
