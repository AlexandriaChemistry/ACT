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

#include "alex_modules.h"
#include "fill_inputrec.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "mymol.h"
#include "poldata_xml.h"
#include "units.h"

int alex_molprop_check(int argc, char*argv[])
{
    static const char               *desc[] = {
        "molprop_check checks calculations for missing hydrogens",
        "and inconsistent dipoles. It also tries to make a topology",
        "and reports errors doing this. Output is to a file."
    };
    t_filenm                         fnm[] =
    {
        { efDAT, "-d",  "gentop",  ffREAD },
        { efDAT, "-f",  "allmols",  ffREAD },
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
            rvec   mu;
            tensor Q;
            double value, error, T = 0;
            std::string type;
            if (ci.getVal(type, MPO_DIPOLE, &value, &error,
                          &T, mu, Q))
            {
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

        alexandria::MyMol mymol;
        mymol.Merge(&m);
        mymol.setInputrec(inputrec);
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
            fprintf(mylog, "%s. Failed to generate topology: %s\n",
                    mymol.getMolname().c_str(), immsg(imm));
            for (const auto &error : mymol.errors())
            {
                fprintf(mylog, "Error: %s\n", error.c_str());
            }
        }
    }
    return 0;
}
