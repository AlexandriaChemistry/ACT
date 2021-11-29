/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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

#include "gromacs/fileio/xvgr.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/textreader.h"

#include "coulombintegrals/coulombintegrals.h"
#include "poldata.h"
#include "regression.h"
#include "units.h"

namespace alexandria
{

void QgenResp::potcomp(const char             *potcomp,
                       const t_atoms          *atoms,
                       const rvec             *x,
                       const char             *pdbdiff,
                       const gmx_output_env_t *oenv)
{
    FILE       *fp;
    std::string unit("Hartree/e");

    if (potcomp)
    {
        const char *pcleg[2] = { "Atoms", "ESP points" };
        fp = xvgropen(potcomp, "Electrostatic potential", unit, unit, oenv);
        xvgr_legend(fp, 2, pcleg, oenv);
        fprintf(fp, "@type xy\n");
        for (size_t i = 0; (i < nEsp()); i++)
        {
            /* Conversion may or may not be in vain depending on unit */
            auto myexp = convertFromGromacs(ep_[i].v(), unit);
            auto myeem = convertFromGromacs(ep_[i].vCalc(), unit);
            if (i == static_cast<size_t>(natoms()))
            {
                fprintf(fp, "&\n");
                fprintf(fp, "@type xy\n");
            }
            fprintf(fp, "%10.5e  %10.5e\n", myexp, myeem);
        }
        fprintf(fp, "&\n");
        fclose(fp);
    }
    std::string potUnit("Hartree/e");
    if (pdbdiff)
    {
        fp = fopen(pdbdiff, "w");
        for (int i = 0; i < atoms->nr; i++)
        {
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "ATOM", 1, *atoms->atomname[i], "MOL", 'A', i+1,
                    ' ', 10*x[i][XX], 10*x[i][YY], 10*x[i][ZZ], 0.0, 0.0);
        }
        double ymin, ymax;
        ymin = ymax = ep_[0].esp()[YY];
        for (size_t i = 0; (i < nEsp()); i++)
        {
            const gmx::RVec esp = ep_[i].esp();
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "HETATM", 1, "HE", "QM", 'B', static_cast<int>(i+1),
                    ' ', 10*esp[XX], 10*esp[YY], 10*esp[ZZ], 0.0, ep_[i].v());
            ymin = std::min(ymin, esp[YY]);
            ymax = std::max(ymax, esp[YY]);
        }
        double dy = 1.2*(ymax-ymin);
        for (size_t i = 0; (i < nEsp()); i++)
        {
            const gmx::RVec esp = ep_[i].esp();
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "HETATM", 1, "HE", "AX", 'C', static_cast<int>(i+1),
                    ' ', 10*esp[XX], 10*(esp[YY]-dy), 10*esp[ZZ], 0.0, ep_[i].vCalc());
        }
        for (size_t i = 0; (i < nEsp()); i++)
        {
            const gmx::RVec esp = ep_[i].esp();
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "HETATM", 1, "HE", "DIF", 'D', static_cast<int>(i+1),
                    ' ', 10*esp[XX], 10*(esp[YY]-2*dy), 10*esp[ZZ], 0.0, 
                                         ep_[i].vCalc() - ep_[i].v());
        }
        fclose(fp);
    }
}

void QgenResp::writeDiffCube(QgenResp               *src,
                             const std::string      &cubeFn,
                             const std::string      &histFn,
                             const std::string      &title,
                             const gmx_output_env_t *oenv,
                             int                     rho)
{
    FILE       *fp;
    int         i, m, ix, iy, iz;
    real        pp, r, rmin;
    gmx_stats_t gst = nullptr, ppcorr = nullptr;

    if (0 != histFn.size())
    {
        gst    = gmx_stats_init();
        ppcorr = gmx_stats_init();
    }
    std::string lengthUnit("Bohr");
    std::string potUnit("Hartree/e");
    if (0 != cubeFn.size())
    {
        fp = gmx_ffopen(cubeFn.c_str(), "w");
        fprintf(fp, "%s\n", title.c_str());
        fprintf(fp, "POTENTIAL\n");
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n",
                static_cast<int>(natoms()),
                convertFromGromacs(origin_[XX], lengthUnit),
                convertFromGromacs(origin_[YY], lengthUnit),
                convertFromGromacs(origin_[ZZ], lengthUnit));
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", nxyz_[XX],
                convertFromGromacs(space_[XX], lengthUnit), 0.0, 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", nxyz_[YY],
                0.0, convertFromGromacs(space_[YY], lengthUnit), 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", nxyz_[ZZ],
                0.0, 0.0, convertFromGromacs(space_[ZZ], lengthUnit));

        for (int m = 0; (m < natoms()); m++)
        {
            fprintf(fp, "%5d%12.6f%12.6f%12.6f%12.6f\n",
                    atomnumber_[m], q_[m],
                    convertFromGromacs(x_[m][XX], lengthUnit),
                    convertFromGromacs(x_[m][YY], lengthUnit),
                    convertFromGromacs(x_[m][ZZ], lengthUnit));
        }

        for (ix = m = 0; ix < nxyz_[XX]; ix++)
        {
            for (iy = 0; iy < nxyz_[YY]; iy++)
            {
                for (iz = 0; iz < nxyz_[ZZ]; iz++, m++)
                {
                    if (src->nEsp() > 0 && nullptr != ppcorr)
                    {
                        gmx_stats_add_point(ppcorr,
                                            convertFromGromacs(src->ep_[m].v(), potUnit),
                                            convertFromGromacs(ep_[m].vCalc(), potUnit), 0, 0);
                    }
                    pp = ep_[m].vCalc();
                    if (!src->ep_.empty())
                    {
                        pp -= src->ep_[m].v();
                    }
                    if (rho == 0)
                    {
                        pp = convertFromGromacs(pp, potUnit);
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
                        for (auto i = 0; i < nAtom_; ++i)
                        {
                            gmx::RVec dx;
                            rvec_sub(x_[i], ep_[m].esp(), dx);
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
    writeDiffCube(&dummy,  fn, empty, title, oenv, 0);
}

void QgenResp::writeRho(const std::string      &fn, 
                        const std::string      &title,
                        const gmx_output_env_t *oenv)
{
    QgenResp    dummy;
    std::string empty;
    writeDiffCube(&dummy,  fn, empty, title, oenv, 1);
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
    std::string         Bohr("Bohr");
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
                    origin_[m] = convertToGromacs(origin_[m], Bohr);
                    space_[m]  = convertToGromacs(space_[m], Bohr);
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
                    atomnumber_[m] = anr;
                    q_[m] = qq;
                }
                x_[m][XX] = convertToGromacs(lx, Bohr);
                x_[m][YY] = convertToGromacs(ly, Bohr);
                x_[m][ZZ] = convertToGromacs(lz, Bohr);
            }
        }
        else if (line >= 6+natom)
        {
            std::vector<std::string> ss = gmx::splitString(tmp);
            for (const auto &s : ss)
            {
                pot.push_back(convertToGromacs(my_atof(s.c_str(), "energy"), "Hartree/e"));
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
                  nAtom_);
    }
}

void QgenResp::makeGrid(real                              spacing, 
                        real                              border,
                        const gmx::HostVector<gmx::RVec> &x)
{
    if (0 != nEsp())
    {
        fprintf(stderr, "Overwriting existing ESP grid\n");
    }
    if (spacing <= 0)
    {
        spacing = 0.005;
        fprintf(stderr, "Spacing too small, setting it to %g nm\n", spacing);
    }
    /* Determine extent of compound */
    rvec xmin = {  100,  100,  100 };
    rvec xmax = { -100, -100, -100 };
    x_ = x;
    for (int i = 0; (i < nAtom_); i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            xmin[m]  = std::min(xmin[m], x[i][m]);
            xmax[m]  = std::max(xmax[m], x[i][m]);
        }
    }
    
    if (border <= 0)
    {
        border = 0.25; // nanometer
        fprintf(stderr, "Border too small, setting it to %g nm\n", border);
    }
    rvec box;
    for (int m = 0; (m < DIM); m++)
    {
        xmin[m]    -= border;
        xmax[m]    += border;
        origin_[m]  = xmin[m];
        box[m]      = xmax[m]-xmin[m];
    }
    for (int m = 0; (m < DIM); m++)
    {
        nxyz_[m]  = 1+(int) (box[m]/spacing);
        space_[m] = box[m]/nxyz_[m];
    }
    ep_.clear();
    for (int i = 0; (i < nxyz_[XX]); i++)
    {
        gmx::RVec xyz;
        xyz[XX] = xmin[XX] + i*space_[XX];
        for (int j = 0; (j < nxyz_[YY]); j++)
        {
            xyz[YY] = xmin[YY] + j*space_[YY];
            for (int k = 0; (k < nxyz_[ZZ]); k++)
            {
                xyz[ZZ] = xmin[ZZ] + k*space_[ZZ];
                ep_.push_back(EspPoint(xyz, 0));
            }
        }
    }
}

void QgenResp::calcRho()
{
    for (size_t i = 0; (i < nEsp()); i++)
    {
        double V = 0;
        for (int j = 0; j < nAtom_; j++)
        {
            double               vv = 0;
            gmx::RVec            dx;
            rvec_sub(ep_[i].esp(), x_[j], dx);
            double               r     = norm(dx);
            if (ChargeType::Gaussian == ChargeType_)
            {
                vv = q_[j]*Nuclear_GG(r, zeta_[j]);
            }
            else if (ChargeType::Slater == ChargeType_)
            {
                vv = q_[j]*Nuclear_SS(r, row_[j], zeta_[j]);
            }
            else
            {
                gmx_fatal(FARGS, "Unsupported charge model %s",
                          chargeTypeName(ChargeType_).c_str());
            }
            V  += vv;
        }
        ep_[i].setRho(V);
    }
}

} // namespace
