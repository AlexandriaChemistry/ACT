/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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
 
 
#include <cmath>
#include <cstdio>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/coulombintegrals/gaussian_integrals.h"
#include "gromacs/coulombintegrals/slater_integrals.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "alex_modules.h"
#include "act/forcefield.h"
#include "act/forcefield_xml.h"

using namespace alexandria;

static void filter_atypes(std::vector<Ffatype>  *FFatypes,
                          char                  *atypes)
{
    if (nullptr != atypes)
    {
        std::vector<std::string> ptr = gmx::splitString(atypes);
        for (auto atp = FFatypes->begin(); atp != FFatypes->end();)
        {
            if (std::find(ptr.begin(), ptr.end(), atp->getType()) == ptr.end())
            {
                atp = FFatypes->erase(atp);
            }
            else
            {
                ++atp;
            }
        }
    }
}

static int faculty(int n)
{
    if (n <= 0)
    {
        return 1;
    }
    else
    {
        return n*faculty(n-1);
    }
}

/*
  ref: Hogervorst, Physica, Volume: 51, Page: 77, Year: 1971.
*/
static void gmx_unused hogervorstCombination(double  si,  double  ei,  
                                             double  gi,
                                             double  sj,  double  ej,  
                                             double  gj,
                                             double *sij, double *eij, 
                                             double *gij)
{
    auto si6 = pow(si, 6);
    auto sj6 = pow(sj, 6);
    
    auto Eij = 2.0*((ei*ej)/(ei + ej));
    auto Gij = 0.5*(gi + gj);
    auto A   = sqrt(((ei*gi*si6)/(gi - 6))*((ej*gj*sj6)/(gj - 6)));    
    auto B   = (A*(Gij - 6))/(Eij*Gij);
    
    *eij = Eij;
    *gij = Gij;    
    *sij = pow(B , (1/6.0));
}

static double Coulomb_PP(double r)
{
    if (r == 0)
    {
        return 0;
    }
    return 1/r;
}

static double DCoulomb_PP(double r)
{
    if (r == 0)
    {
        return 0;
    }
    return 1/gmx::square(r);
}

static void lj(double r,
               double *vd, double *fd,
               double *vr, double *fr)
{
    double r2 = 0, r_6 = 0, r_12 = 0;

    r2    = r*r;
    r_6   = 1.0/(r2*r2*r2);
    r_12  = r_6*r_6;

    *vd   = -r_6;             /*  g(c)     Dispersion */
    *fd   =  6.0*(*vd)/r;     /* -g'(x)               */

    *vr   = r_12;             /*  h(x)     Repulsion  */
    *fr   = 12.0*(*vr)/r;     /* -h'(x)               */
}

static void gmx_unused wbk(double r,   double eps, double sig,
                           double gam, double *vd, double *fd,
                           double *vr, double *fr)
{
    auto r2    = r*r;
    auto r5    = r2*r2*r;
    auto r6    = r5*r;
    auto sig2  = sig*sig;
    auto sig6  = sig2*sig2*sig2;
    auto sig5  = sig2*sig2*sig;
    auto A     = std::exp(gam*(1-(r/sig)));
    auto B     = sig6 + r6;
    auto C     = gam + 3;
    
    *vd        = -2*eps*(1.0/(1-(3.0/C))*(sig6/B));                          /*  g(c)     Dispersion */
    *fd        = -2*eps*((6*C*r5*sig6)/(gam*(B*B)));                         /* -g'(x)               */
    
    *vr        = (2*eps*(1.0/(1-(3.0/C))*(sig6/B)))*((3.0/C)*A);             /*  h(x)     Repulsion  */            
    *fr        = (6*A*eps*sig5*(gam*r6 + 6*r5*sig + gam*sig6))/(gam*(B*B));  /* -h'(x)               */    
}


static void gmx_unused gen_alexandria_rho(ForceField                 &pd,
                                          const char             *fn,
                                          ChargeType              iType,
                                          real                    rcut,
                                          real                    spacing,
                                          const gmx_output_env_t *oenv)
{
    FILE                         *fp;
    int                           j     = 0;
    int                           n     = 0;
    int                           nmax  = 0;
    int                           nzeta = 0;
    int                          *row;
    std::string                   name;
    double                        rho  = 0;
    double                        rr   = 0;
    double                        qtot = 0;
    double                       *A;
    double                       *zeta;
    double                       *q;
    
    char                          buf[STRLEN];

    nmax = 1+(int)(rcut/spacing);
    for (auto &eep : pd.eempropsConst())
    {
        if (pd.chargeType() == iType)
        {
            name  = eep.first;
            nzeta = eep.second.getNzeta();
            snew(zeta, nzeta);
            snew(q, nzeta);
            snew(row, nzeta);
            snew(A, nzeta);
            qtot = 0;
            for (j = 0; j < nzeta; j++)
            {
                zeta[j] = eep.second.getZeta(j);
                q[j]    = eep.second.getQ(j);
                qtot   += q[j];
                row[j]  = eep.second.getRow(j);
                switch (iType)
                {
                case eqtGaussian:
                    A[j] = pow(zeta[j]*zeta[j]/M_PI, 1.5);
                    break;
                case eqtSlater:
                    A[j] = pow(2*zeta[j], 2*row[j]+1)/(4*M_PI*faculty(2*row[j]));
                    break;
                default: // gmx_fatal called
                    gmx_fatal(FARGS, "Don't know how to handle model %s",
                              chargeTypeName(pd.chargeType()).c_str());
                }
            }
            if (q[nzeta-1] == 0)
            {
                q[nzeta-1] = -qtot;
            }
            sprintf(buf, "%s_%s", name.c_str(), fn);
            fp = xvgropen(buf, "Rho", "r (nm)", "rho(r)", oenv);
            for (n = 0; n <= nmax; n++)
            {
                rr  = n*spacing;
                rho = 0;
                for (j = 0; j < nzeta; j++)
                {
                    if (zeta[j] > 0)
                    {
                        switch (iType)
                        {
                        case eqtGaussian:
                            rho += A[j]*exp(-gmx::square(rr*zeta[j]));
                            break;
                        case eqtSlater:
                            rho += A[j]*pow(rr, 2*row[j]-2)*exp(-2*zeta[j]*rr);
                            break;
                        default: // gmx_fatal called
                            gmx_fatal(FARGS, "Don't know how to handle model %s",
                                      chargeTypeName(iType).c_str());
                        }
                    }
                }
                fprintf(fp, "%10.5e  %10.5e\n", rr, rho);
            }
            fclose(fp);
            sfree(q);
            sfree(zeta);
            sfree(row);
            sfree(A);
        }
    }
}

static void gen_alexandria_tables(ForceField                 &pd,
                                  const char              *fn,
                                  real                     rcut,
                                  real                     spacing,
                                  const gmx_output_env_t  *oenv,
                                  char                    *atypes)
{
    ChargeType   iType = pd.chargeType();
                                  
    double       cv = 0;
    double       cf = 0;
    double       rr = 0;
    double       vd = 0;
    double       fd = 0;
    double       vr = 0;
    double       fr = 0;
    
    // TODO Fix this
    const char  *ns[2]  = {"", "_s"};
    
    FILE        *fp1 = nullptr;
    FILE        *fp2 = nullptr;

    auto nmax     = 1+(int)(rcut/spacing);
    auto FFatypes = pd.getAtypes();
    
    filter_atypes(&FFatypes, atypes);        
    for (auto atpi = FFatypes.begin(); atpi != FFatypes.end(); atpi++)
    {
        auto eei    = pd.ztype2Eem(atpi->getZtype());
        auto nzetaI = eei->getNzeta();
        for (auto atpj = atpi; atpj != FFatypes.end(); atpj++)
        {
            auto eej    = pd.ztype2Eem(atpj->getZtype());
            auto nzetaJ =  eej->getNzeta();
            for (auto i = 0; i < nzetaI; i++)
            {
                auto zetaI = eei->getZeta(i);
                auto rowI  = eei->getRow(i);
                for (auto j = 0; j < nzetaJ; j++)
                {
                    auto zetaJ = eej->getZeta(j);
                    auto rowJ  = eej->getRow(j);

                    std::string fnbuf(fn, fn+strlen(fn)-4);
                    std::string buf1 = 
                        gmx::formatString("%s_%s%s_%s%s.xvg",
                                          fnbuf.c_str(),
                                          atpi->getType().c_str(), 
                                          ns[i], atpj->getType().c_str(), 
                                          ns[j]);
                    if (atpi->getType() != atpj->getType())
                    {
                        std::string buf2 =
                            gmx::formatString("%s_%s%s_%s%s.xvg", 
                                              fnbuf.c_str(), 
                                              atpj->getType().c_str(), 
                                              ns[j], atpi->getType().c_str(), 
                                              ns[i]);
                        fp2 = xvgropen(buf2.c_str(), buf2.c_str(),
                                       "r (nm)", "V (kJ/mol e)", oenv);
                    }
                    fp1 = xvgropen(buf1.c_str(), buf1.c_str(),
                                   "r (nm)", "V (kJ/mol e)", oenv);
                    
                    for (auto n = 0; n <= nmax; n++)
                    {
                        rr = n*spacing;
                        switch (iType)
                        {
                        case eqtPoint:
                            cv = Coulomb_PP(rr);
                            cf = DCoulomb_PP(rr);
                            break;
                        case eqtGaussian:
                            cv = Coulomb_GG(rr, zetaI, zetaJ);
                            if (zetaI == 0 && zetaJ == 0 && rr != 0)
                            {
                                cf = -DCoulomb_GG(rr, zetaI, zetaJ);
                            }
                            else
                            {
                                cf = DCoulomb_GG(rr, zetaI, zetaJ);
                            }
                            break;
                        case eqtSlater:
                            cv = Coulomb_SS(rr, rowI, rowJ, zetaI, zetaJ);
                            if (rowI == 0 && rowJ == 0 && rr !=0)
                            {
                                cf = -DCoulomb_SS(rr, rowI, rowJ, zetaI, zetaJ);
                            }
                            else
                            {
                                cf = DCoulomb_SS(rr, rowI, rowJ, zetaI, zetaJ);
                            }
                            break;
                        default: // gmx_fatal called
                            gmx_fatal(FARGS, "Don't know how to handle model %s",
                                      chargeTypeName(iType).c_str());
                        }
                        if (rr > 0)
                        {
                            lj(rr, &vd, &fd, &vr, &fr);
                        }
                        else
                        {
                            vd = fd = vr = fr = 0;
                        }
                        fprintf(fp1, "%10.5e  %10.5e  %10.5e %10.5e %10.5e %10.5e %10.5e\n", rr, cv, cf, vd, fd, vr, fr);
                        if (atpi->getType() != atpj->getType())
                        {
                            fprintf(fp2, "%10.5e  %10.5e  %10.5e %10.5e %10.5e %10.5e %10.5e\n", rr, cv, cf, vd, fd, vr, fr);
                        }
                    }
                    fclose(fp1);
                    if (atpi->getType() != atpj->getType() && fp2)
                    {
                        fclose(fp2);
                    }
                }
            }
        }
    }
}

int alex_gen_table(int argc, char *argv[])
{
    static const char            *desc[] = {
        "gen_table generates tables for mdrun for use with the USER defined",
        "potentials. Note that the format has been update for higher",
        "accuracy in the forces starting with version 4.0. Using older",
        "tables with 4.0 will silently crash your simulations, as will",
        "using new tables with an older GROMACS version. This is because in the",
        "old version the second derevative of the potential was specified",
        "whereas in the new version the first derivative of the potential",
        "is used instead.[PAR]",
        "The program can read the [TT]ff.xml[tt] file (or otherwise as",
        "specified with the [TT]-di[tt] option) and generate tables for all",
        "possible interactions for the charge model specified in the file.[PAR]",
        "Broken: For Slater interactions four parameters must be passed: the 1/Width",
        "and the row number of the element. The interactions are computed analytically",
        "which may be slow due to the fact that arbitraray precision arithmetic is",
        "needed. If the width of one of the Slater is zero a Nucleus-Slater interaction",
        "will be generated."
    };

    t_filenm                      fnm[] = {
        { efXVG, "-o", "table",  ffWRITE },
        { efXML, "-ff", "aff", ffOPTRD }
    };
    
    const  int                    NFILE      = asize(fnm);
    
    static int                    nrow1      = 1;
    static int                    nrow2      = 1;
    static int                    nrep       = 12;
    static int                    ndisp      = 6;
    static int                    pts_nm     = 500;
    static double                 delta      = 0;
    static double                 efac       = 500;
    static double                 rc         = 0.9;
    static double                 rtol       = 1e-05;
    static double                 xi         = 0.15;
    static double                 xir        = 0.0615;
    static double                 w1         = 20;
    static double                 w2         = 20;
    static char                  *atypes     = nullptr;
    static const char            *opt[]      = {nullptr, "cut", "rf", "pme", nullptr};
    static const char            *vdw[]      = {nullptr, "lj", "bk", "wbk", nullptr};
    
    t_pargs                       pa[]       = {
        { "-el",     FALSE, etENUM, {opt},
          "Electrostatics type: cut, rf or pme" },
        { "-rc",     FALSE, etREAL, {&rc},
          "Cut-off required for rf or pme" },
        { "-rtol",   FALSE, etREAL, {&rtol},
          "Ewald tolerance required for pme" },
        { "-xi",   FALSE, etREAL, {&xi},
          "Width of the Gaussian diffuse charge of the G&G model" },
        { "-xir",   FALSE, etREAL, {&xir},
          "Width of erfc(z)/z repulsion of the G&G model (z=0.5 rOO/xir)" },
        { "-z1",   FALSE, etREAL, {&w1},
          "1/Width of the first Slater charge (unit 1/nm)" },
        { "-z2",   FALSE, etREAL, {&w2},
          "1/Width of the second Slater charge (unit 1/nm)" },
        { "-nrow1",   FALSE, etINT, {&nrow1},
          "Row number for the first Slater charge" },
        { "-nrow2",   FALSE, etINT, {&nrow2},
          "Row number for the first Slater charge" },
        { "-vdw",      FALSE, etENUM, {vdw},
          "Model for the tables" },
        { "-resol",  FALSE, etINT,  {&pts_nm},
          "Resolution of the table (points per nm)" },
        { "-delta",  FALSE, etREAL, {&delta},
          "Displacement in the Coulomb functions (nm), used as 1/(r+delta). Only for hard wall potential." },
        { "-efac",   FALSE, etREAL, {&efac},
          "Number indicating the steepness of the hardwall potential." },
        { "-nrep",   FALSE, etINT,  {&nrep},
          "Power for the repulsion potential (with model AB1 or maaren)" },
        { "-ndisp",   FALSE, etINT,  {&ndisp},
          "Power for the dispersion potential (with model AB1 or maaren)" },
        { "-atypes",  FALSE, etSTR, {&atypes},
          "List of atom types for which you want to generate table." }
    };
    
    gmx_output_env_t             *oenv;
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
     
    ForceField                  pd;
    gmx_atomprop_t           aps = gmx_atomprop_init();
    const char              *aff = opt2fn_null("-ff", NFILE, fnm);
    try
    {
        readForceField(aff, pd, aps);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    
    gen_alexandria_tables(pd, opt2fn("-o", NFILE, fnm), rc, 1.0/pts_nm, oenv, atypes);
    
    return 0;
}
