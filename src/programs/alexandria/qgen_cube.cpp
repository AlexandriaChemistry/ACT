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

#include "qgen_resp.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>

#include <map>

#include "gromacs/commandline/filenm.h"
#include "gromacs/coulombintegrals/coulombintegrals.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/textreader.h"

#include "poldata.h"
#include "regression.h"

namespace alexandria
{

void QgenResp::potcomp(const char             *potcomp,
                       const char             *pdbdiff,
                       const gmx_output_env_t *oenv)
{
    double  pp, exp, eem;
    FILE   *fp;
    int     unit = eg2cHartree_e;

    if (potcomp)
    {
        const char *pcleg[2] = { "Atoms", "ESP points" };
        fp = xvgropen(potcomp, "Electrostatic potential", unit2string(unit), unit2string(unit), oenv);
        xvgr_legend(fp, 2, pcleg, oenv);
        fprintf(fp, "@type xy\n");
        for (size_t i = 0; (i < nEsp()); i++)
        {
            /* Conversion may or may not be in vain depending on unit */
            exp = gmx2convert(ep_[i].v(), unit);
            eem = gmx2convert(ep_[i].vCalc(), unit);
            if (i == nRespAtom())
            {
                fprintf(fp, "&\n");
                fprintf(fp, "@type xy\n");
            }
            fprintf(fp, "%10.5e  %10.5e\n", exp, eem);
        }
        fprintf(fp, "&\n");
        fclose(fp);
    }
    if (pdbdiff)
    {
        fp = fopen(pdbdiff, "w");
        fprintf(fp, "REMARK All distances are scaled by a factor of two.\n");
        for (size_t i = 0; (i < nEsp()); i++)
        {
            exp = gmx2convert(ep_[i].v(), eg2cHartree_e);
            eem = gmx2convert(ep_[i].vCalc(), eg2cHartree_e);
            pp  = ep_[i].v()-ep_[i].vCalc();
            const gmx::RVec esp = ep_[i].esp();
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "ATOM", 1, "HE", "HE", ' ', static_cast<int>(i+1),
                    ' ', 20*esp[XX], 20*esp[YY], 20*esp[ZZ], 0.0, pp);
        }
        fclose(fp);
    }
}

void QgenResp::writeDiffCube(QgenResp               &src,
                             const std::string      &cubeFn,
                             const std::string      &histFn,
                             const std::string      &title,
                             const gmx_output_env_t *oenv,
                             int                     rho)
{
    FILE       *fp;
    int         i, m, ix, iy, iz;
    real        pp, q, r, rmin;
    gmx_stats_t gst = nullptr, ppcorr = nullptr;

    if (0 != histFn.size())
    {
        gst    = gmx_stats_init();
        ppcorr = gmx_stats_init();
    }
    if (0 != cubeFn.size())
    {
        fp = gmx_ffopen(cubeFn.c_str(), "w");
        fprintf(fp, "%s\n", title.c_str());
        fprintf(fp, "POTENTIAL\n");
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n",
                static_cast<int>(nRespAtom()),
                gmx2convert(origin_[XX], eg2cBohr),
                gmx2convert(origin_[YY], eg2cBohr),
                gmx2convert(origin_[ZZ], eg2cBohr));
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", nxyz_[XX],
                gmx2convert(space_[XX], eg2cBohr), 0.0, 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", nxyz_[YY],
                0.0, gmx2convert(space_[YY], eg2cBohr), 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", nxyz_[ZZ],
                0.0, 0.0, gmx2convert(space_[ZZ], eg2cBohr));

        for (size_t m = 0; (m < nRespAtom()); m++)
        {
            q = ra_[m].q();
            fprintf(fp, "%5d%12.6f%12.6f%12.6f%12.6f\n",
                    ra_[m].atomnumber(), q,
                    gmx2convert(ra_[m].x()[XX], eg2cBohr),
                    gmx2convert(ra_[m].x()[YY], eg2cBohr),
                    gmx2convert(ra_[m].x()[ZZ], eg2cBohr));
        }

        for (ix = m = 0; ix < nxyz_[XX]; ix++)
        {
            for (iy = 0; iy < nxyz_[YY]; iy++)
            {
                for (iz = 0; iz < nxyz_[ZZ]; iz++, m++)
                {
                    if (src.nEsp() > 0 && nullptr != ppcorr)
                    {
                        gmx_stats_add_point(ppcorr,
                                            gmx2convert(src.ep_[m].v(), eg2cHartree_e),
                                            gmx2convert(ep_[m].vCalc(), eg2cHartree_e), 0, 0);
                    }
                    pp = ep_[m].vCalc();
                    if (!src.ep_.empty())
                    {
                        pp -= src.ep_[m].v();
                    }
                    if (rho == 0)
                    {
                        pp = gmx2convert(pp, eg2cHartree_e);
                    }
                    else
                    {
                        pp = ep_[m].rho()*pow(BOHR2NM, 3);
                    }
                    fprintf(fp, "%13.5e", pp);
                    if (iz % 6 == 5)
                    {
                        fprintf(fp, "\n");
                    }
                    if (nullptr != gst)
                    {
                        rmin = 1000;
                        /* Add point to histogram! */
                        for (auto i = ra_.begin(); i < ra_.end(); ++i)
                        {
                            gmx::RVec dx;
                            rvec_sub(i->x(), ep_[m].esp(), dx);
                            r = norm(dx);
                            if (r < rmin)
                            {
                                rmin = r;
                            }
                        }
                        gmx_stats_add_point(gst, rmin, pp, 0, 0);
                    }
                }
                if ((iz % 6) != 0)
                {
                    fprintf(fp, "\n");
                }
            }
        }
        fclose(fp);
    }
    if (nullptr != gst)
    {
        int   nb = 0;
        real *x  = nullptr, *y = nullptr;

        fp = xvgropen(histFn.c_str(), "Absolute deviation from QM", "Distance (nm)",
                      "Potential", oenv);
        gmx_stats_dump_xy(gst, fp);
        if (0)
        {
            gmx_stats_make_histogram(gst, 0.01, &nb, ehistoX, 0, &x, &y);
            gmx_stats_free(gst);
            for (i = 0; (i < nb); i++)
            {
                fprintf(fp, "%10g  %10g\n", x[i], y[i]);
            }
            free(x);
            free(y);
        }
        fclose(fp);
        fp = xvgropen("diff-pot.xvg", "Correlation between QM and Calc", "Pot (QM)",
                      "Pot (Calc)", oenv);
        gmx_stats_dump_xy(ppcorr, fp);
        fclose(fp);
    }
}

void QgenResp::writeCube(const std::string      &fn, 
                         const std::string      &title,
                         const gmx_output_env_t *oenv)
{
    QgenResp    dummy;
    std::string empty;
    writeDiffCube(dummy,  fn, empty, title, oenv, 0);
}

void QgenResp::writeRho(const std::string      &fn, 
                        const std::string      &title,
                        const gmx_output_env_t *oenv)
{
    QgenResp    dummy;
    std::string empty;
    writeDiffCube(dummy,  fn, empty, title, oenv, 1);
}

void QgenResp::readCube(const std::string &fn, bool bESPonly)
{
    int                 natom, nxyz[DIM] = { 0, 0, 0 };
    double              space[DIM] = { 0, 0, 0 };
    std::vector<double> pot;

    gmx::TextReader     tr(fn);
    std::string         tmp;
    int                 line = 0;
    bool                bOK  = true;
    while (bOK && tr.readLine(&tmp))
    {
        while (!tmp.empty() && tmp[tmp.length()-1] == '\n')
        {
            tmp.erase(tmp.length()-1);
        }
        if (0 == line)
        {
            printf("%s\n", tmp.c_str());
        }
        else if (1 == line && tmp.compare("POTENTIAL") != 0)
        {
            bOK = false;
        }
        else if (2 == line)
        {
            double origin[DIM];
            bOK = (4 == sscanf(tmp.c_str(), "%d%lf%lf%lf",
                               &natom, &origin[XX], &origin[YY], &origin[ZZ]));
            if (bOK && !bESPonly)
            {
                origin_[XX] = origin[XX];
                origin_[YY] = origin[YY];
                origin_[ZZ] = origin[ZZ];
            }
        }
        else if (3 == line)
        {
            bOK = (2 == sscanf(tmp.c_str(), "%d%lf",
                               &nxyz[XX], &space[XX]));
        }
        else if (4 == line)
        {
            bOK = (2 == sscanf(tmp.c_str(), "%d%*s%lf",
                               &nxyz[YY], &space[YY]));
        }
        else if (5 == line)
        {
            bOK = (2 == sscanf(tmp.c_str(), "%d%*s%*s%lf",
                               &nxyz[ZZ], &space[ZZ]));
            if (bOK)
            {
                for (int m = 0; (m < DIM); m++)
                {
                    nxyz_[m]  = nxyz[m];
                    space_[m] = space[m];
                }
                for (int m = 0; (m < DIM); m++)
                {
                    origin_[m] = convert2gmx(origin_[m], eg2cBohr);
                    space_[m]  = convert2gmx(space_[m], eg2cBohr);
                }
            }
            pot.clear();
        }
        else if (line >= 6 && line < 6+natom)
        {
            double lx, ly, lz, qq;
            int    anr, m = line - 6;
            bOK = (5 == sscanf(tmp.c_str(), "%d%lf%lf%lf%lf",
                               &anr, &qq, &lx, &ly, &lz));
            if (bOK)
            {
                if (!bESPonly)
                {
                    ra_[m].setAtomnumber(anr);
                    ra_[m].setQ(qq);
                }
                gmx::RVec xx;
                xx[XX] = convert2gmx(lx, eg2cBohr);
                xx[YY] = convert2gmx(ly, eg2cBohr);
                xx[ZZ] = convert2gmx(lz, eg2cBohr);
                ra_[m].setX(xx);
            }
        }
        else if (line >= 6+natom)
        {
            std::vector<std::string> ss = gmx::splitString(tmp);
            for (const auto &s : ss)
            {
                pot.push_back(convert2gmx(my_atof(s.c_str(), "energy"), eg2cHartree_e));
            }
        }

        line++;
    }
    if (bOK)
    {
        ep_.clear();
        int m = 0;
        for (int ix = 0; ix < nxyz_[XX]; ix++)
        {
            for (int iy = 0; iy < nxyz_[YY]; iy++)
            {
                for (int iz = 0; iz < nxyz_[ZZ]; iz++, m++)
                {
                    gmx::RVec e;
                    e[XX] = origin_[XX] + ix*space_[XX];
                    e[YY] = origin_[YY] + iy*space_[YY];
                    e[ZZ] = origin_[ZZ] + iz*space_[ZZ];

                    ep_.push_back(EspPoint(e, pot[m]));
                }
            }
        }
    }
    if (!bOK)
    {
        gmx_fatal(FARGS, "Error reading %s. Found %d potential values, %d coordinates and %d atoms",
                  fn.c_str(), static_cast<int>(pot.size()), static_cast<int>(ep_.size()),
                  static_cast<int>(ra_.size()));
    }
}

} // namespace
