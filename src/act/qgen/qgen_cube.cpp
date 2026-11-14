/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021,2025,2026
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

#include "qgen_resp.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <map>

#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/statistics/statistics.h"
#include "act/topology/topology.h"
#include "act/utility/units.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/textreader.h"

namespace alexandria
{

void QgenResp::potcomp(const char             *potcomp,
                       const gmx_output_env_t *oenv)
{
    FILE       *fp;
    const std::string unit("Hartree/e");

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
}

void QgenResp::writePdbComparison(const std::vector<ActAtom>   &atoms,
                                  const std::string            &pdbdiff)
{
    const std::string potUnit("Hartree/e");
    double xfac = convertFromGromacs(1, "Angstrom");
    if (!pdbdiff.empty())
    {
        FILE *fp = gmx_ffopen(pdbdiff, "w");
        for (size_t i = 0; i < atoms.size(); i++)
        {
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4lu%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "ATOM", 1, atoms[i].name().c_str(), "MOL", 'A', i+1,
                    ' ', xfac*x_[i][XX], xfac*x_[i][YY], xfac*x_[i][ZZ], 0.0, 0.0);
        }
        double ymin, ymax;
        ymin = ymax = ep_[0].esp()[YY];
        for (size_t i = 0; (i < nEsp()); i++)
        {
            const gmx::RVec esp = ep_[i].esp();
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "HETATM", 1, "HE", "QM", 'B', static_cast<int>(i+1),
                    ' ', xfac*esp[XX], xfac*esp[YY], xfac*esp[ZZ], 0.0, ep_[i].v());
            ymin = std::min<double>(ymin, esp[YY]);
            ymax = std::max<double>(ymax, esp[YY]);
        }
        double dy = 1.2*(ymax-ymin);
        for (size_t i = 0; (i < nEsp()); i++)
        {
            const gmx::RVec esp = ep_[i].esp();
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "HETATM", 1, "HE", "AX", 'C', static_cast<int>(i+1),
                    ' ', xfac*esp[XX], xfac*(esp[YY]-dy), xfac*esp[ZZ], 0.0, ep_[i].vCalc());
        }
        for (size_t i = 0; (i < nEsp()); i++)
        {
            const gmx::RVec esp = ep_[i].esp();
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "HETATM", 1, "HE", "DIF", 'D', static_cast<int>(i+1),
                    ' ', xfac*esp[XX], xfac*(esp[YY]-2*dy), xfac*esp[ZZ], 0.0, 
                                         ep_[i].vCalc() - ep_[i].v());
        }
        gmx_ffclose(fp);
    }
}

void QgenResp::writeDiffCube(const QgenResp         *src,
                             const std::string      &cubeFn,
                             const std::string      &histFn,
                             const std::string      &title,
                             const gmx_output_env_t *oenv,
                             int                     rho)
{
    FILE       *fp;
    int         m, ix, iy, iz;
    real        pp, r, rmin;
    gmx_stats   gst, ppcorr;

    bool        histo = (!histFn.empty());
    const std::string lengthUnit("Bohr");
    const std::string potUnit("Hartree/e");
    if (!cubeFn.empty())
    {
        fp = gmx_ffopen(cubeFn, "w");
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

        for (size_t m = 0; (m < natoms()); m++)
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
                    if (src->nEsp() > 0 && histo)
                    {
                        ppcorr.add_point(convertFromGromacs(src->ep_[m].v(), potUnit),
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
                    if (histo)
                    {
                        rmin = 1000;
                        /* Add point to histogram! */
                        for (size_t i = 0; i < nAtom_; ++i)
                        {
                            gmx::RVec dx;
                            rvec_sub(x_[i], ep_[m].esp(), dx);
                            r = norm(dx);
                            if (r < rmin)
                            {
                                rmin = r;
                            }
                        }
                        gst.add_point(rmin, pp, 0, 0);
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
    if (histo)
    {
        fp = xvgropen(histFn.c_str(), "Absolute deviation from QM", "Distance (nm)",
                      "Potential", oenv);
        gst.dump_xy(fp);

        fclose(fp);
        // The correlation is written next to the histogram
        const auto corrFn = histFn.substr(0, histFn.rfind('.')) + "-corr.xvg";
        fp = xvgropen(corrFn.c_str(), "Correlation between QM and Calc", "Pot (QM)",
                      "Pot (Calc)", oenv);
        ppcorr.dump_xy(fp);
        fclose(fp);
    }
}

void QgenResp::writeCube(const std::string      &fn, 
                         const std::string      &title,
                         const gmx_output_env_t *oenv)
{
    const QgenResp dummy;
    const std::string empty;
    writeDiffCube(&dummy,  fn, empty, title, oenv, 0);
}

void QgenResp::writeRho(const std::string      &fn, 
                        const std::string      &title,
                        const gmx_output_env_t *oenv)
{
    const QgenResp dummy;
    const std::string empty;
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
    const std::string   Bohr("Bohr");
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
            bOK = (4 == sscanf(tmp.c_str(), "%d%lf%lf%lf",
                               &natom, &origin_[XX], &origin_[YY], &origin_[ZZ]));
            if (bOK && static_cast<size_t>(natom) != nAtom_)
            {
                gmx_fatal(FARGS, "%s has %d atoms but the molecule %zu",
                          fn.c_str(), natom, nAtom_);
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
            const std::vector<std::string> ss = gmx::splitString(tmp);
            for (const auto &s : ss)
            {
                pot.push_back(convertToGromacs(my_atof(s, "energy"), "Hartree/e"));
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

                    ep_.emplace_back(e, pot[m]);
                }
            }
        }
    }
    if (!bOK)
    {
        gmx_fatal(FARGS, "Error reading %s. Found %lu potential values, %lu coordinates and %lu atoms",
                  fn.c_str(), pot.size(), ep_.size(), nAtom_);
    }
}

void QgenResp::makeGrid(real                          spacing, 
                        real                          border,
                        const std::vector<gmx::RVec> &x)
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
    for (size_t i = 0; (i < nAtom_); i++)
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
                ep_.emplace_back(xyz, 0);
            }
        }
    }
}

/*! \brief The charge density of a unit charge
 * \param[in] cdt  The charge distribution
 * \param[in] r    Distance from the center in nm
 * \param[in] zeta Width in 1/nm
 * \param[in] row  Row in the periodic table, for Slater distributions
 * \return the density in 1/nm^3, of which Nuclear_GG or Nuclear_SS is the
 *         potential
 */
static double chargeDensity(ChargeDistributionType cdt, double r, double zeta, int row)
{
    if (zeta == 0)
    {
        // A point charge has no density outside the center
        return 0;
    }
    switch (cdt)
    {
    case ChargeDistributionType::Gaussian:
        return std::pow(zeta, 3)/std::pow(M_PI, 1.5)*std::exp(-gmx::square(zeta*r));
    case ChargeDistributionType::Slater:
        {
            int n = std::max(1, row);
            return std::pow(2*zeta, 2*n+1)/(4*M_PI*std::tgamma(2*n+1))*
                std::pow(r, 2*n-2)*std::exp(-2*zeta*r);
        }
    default:
        gmx_fatal(FARGS, "Unsupported charge model %s",
                  chargeDistributionTypeName(cdt).c_str());
    }
    return 0;
}

void QgenResp::calcRho()
{
    for (size_t i = 0; (i < nEsp()); i++)
    {
        double rho = 0;
        for (size_t j = 0; j < nAtom_; j++)
        {
            gmx::RVec dx;
            rvec_sub(ep_[i].esp(), x_[j], dx);
            double    r = norm(dx);
            rho += q_[j]*chargeDensity(ChargeDistributionType_, r, zeta_[j], row_[j]);
            if (qVS_[j] != 0)
            {
                rho += qVS_[j]*chargeDensity(ChargeDistributionType_, r, zetaVS_[j], rowVS_[j]);
            }
        }
        ep_[i].setRho(rho);
    }
}

} // namespace
