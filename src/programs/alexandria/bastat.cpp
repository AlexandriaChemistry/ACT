/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2019 
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
#include <string.h>

#include <cmath>

#include <algorithm>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/real.h"

#include "alex_modules.h"
#include "molprop_util.h"
#include "mymol.h"
#include "poldata_xml.h"
#include "stringutil.h"

typedef struct {
    std::string      a1, a2;
    int              order;
    std::vector<int> histo;
    gmx_stats_t      lsq;
} t_bond;

typedef struct {
    std::string      a1, a2, a3;
    std::vector<int> histo;
    gmx_stats_t      lsq;
} t_angle;

typedef struct {
    std::string      a1, a2, a3, a4;
    std::vector<int> histo;
    gmx_stats_t      lsq;
} t_dih;

typedef struct {
    std::vector<t_bond>  bond;
    std::vector<t_angle> angle;
    std::vector<t_angle> linangle;
    std::vector<t_dih>   dih;
    std::vector<t_dih>   imp;
} t_bonds;

static void sort_dihs(std::vector<t_dih> &dih)
{
    std::sort(dih.begin(), dih.end(),
              [](const t_dih &a, const t_dih &b)
              {
                  if (a.a1 == b.a1)
                  {
                      if (a.a2 == b.a2)
                      {
                          if (a.a3 == b.a3)
                          {
                              return (a.a4 < b.a4);
                          }
                          return (a.a3 < b.a3);
                      }
                      return (a.a2 < b.a2);
                  }
                  return (a.a1 < b.a1);
              });
}

static void sort_bonds(t_bonds *b)
{
    std::sort(b->bond.begin(), b->bond.end(),
              [](const t_bond &a, const t_bond &b)
              {
                  if (a.a1 == b.a1)
                  {
                      return (a.a2 < b.a2);
                  }
                  return (a.a1 < b.a1);
              });
    std::sort(b->angle.begin(), b->angle.end(),
              [](const t_angle &a, const t_angle &b)
              {
                  if (a.a1 == b.a1)
                  {
                      if (a.a2 == b.a2)
                      {
                          return (a.a3 < b.a3);
                      }
                      return (a.a2 < b.a2);
                  }
                  return (a.a1 < b.a1);
              });
    sort_dihs(b->dih);
    sort_dihs(b->imp);
}

static void add_bond(FILE *fplog, const char *molname, t_bonds *bonds,
                     const std::string a1, const std::string a2,
                     double blen, double spacing, int order)
{
    GMX_RELEASE_ASSERT(a1.size() > 0, "atom name a1 is empty");
    GMX_RELEASE_ASSERT(a2.size() > 0, "atom name a2 is empty");
    size_t index   = std::lround(blen/spacing);
    auto   b       = std::find_if(bonds->bond.begin(), bonds->bond.end(),
                                  [a1, a2, order](t_bond b)
                                  {
                                      return ((((b.a1.compare(a1) == 0) && (b.a2.compare(a2) == 0)) ||
                                               ((b.a1.compare(a2) == 0) && (b.a2.compare(a1) == 0))) &&
                                              (b.order == order));
                                  });
    if (b == bonds->bond.end())
    {
        t_bond bb;
        if (a1 < a2)
        {
            bb.a1.assign(a1);
            bb.a2.assign(a2);
        }
        else
        {
            bb.a1.assign(a2);
            bb.a2.assign(a1);
        }
        bb.order  = order;
        bb.histo.resize(2*index+1, 0);
        bb.lsq    = gmx_stats_init();
        bonds->bond.push_back(std::move(bb));
        b = bonds->bond.end()-1;
    }
    if (index >= b->histo.size())
    {
        b->histo.resize(index+100, 0);
    }
    gmx_stats_add_point(b->lsq, 0, blen, 0, 0);
    b->histo[index]++;
    if (nullptr != fplog)
    {
        fprintf(fplog, "%s bond-%s-%s-%d %g\n", 
                molname,
                a1.c_str(), 
                a2.c_str(), 
                order, 
                blen);
    }
}

static void lo_add_angle(FILE *fplog, const char *molname, std::vector<t_angle> &angle,
                         const std::string a1, const std::string a2, const std::string a3,
                         double refValue, double spacing)
{
    GMX_RELEASE_ASSERT(a1.size() > 0, "atom name a1 is empty");
    GMX_RELEASE_ASSERT(a2.size() > 0, "atom name a2 is empty");
    GMX_RELEASE_ASSERT(a3.size() > 0, "atom name a3 is empty");

    size_t index = std::lround(refValue/spacing);
    auto   a     = std::find_if(angle.begin(), angle.end(),
                                [a1, a2, a3](const t_angle &a)
                                {
                                    int d = a.a2.compare(a2);
                                    if (0 == d)
                                    {
                                        return ((a.a1.compare(a1) == 0 && a.a3.compare(a3) == 0) ||
                                                (a.a1.compare(a3) == 0 && a.a3.compare(a1) == 0));
                                    }
                                    return false;
                                });

    if (a == angle.end())
    {
        t_angle aa;
        aa.a2.assign(a2);
        if (a1 < a3)
        {
            aa.a1.assign(a1);
            aa.a3.assign(a3);
        }
        else
        {
            aa.a1.assign(a3);
            aa.a3.assign(a1);
        }
        aa.histo.resize((int) (180/spacing) + 1, 0);
        aa.lsq = gmx_stats_init();
        angle.push_back(std::move(aa));
        a = angle.end()-1;
    }
    gmx_stats_add_point(a->lsq, 0, refValue, 0, 0);
    a->histo[index]++;
    if (nullptr != fplog)
    {
        fprintf(fplog, "%s angle-%s-%s-%s %g\n", molname,
                a1.c_str(), a2.c_str(), a3.c_str(), refValue);
    }
}

static void add_angle(FILE *fplog, const char *molname, t_bonds *b,
                      const std::string a1, const std::string a2, const std::string a3,
                      double refValue, double spacing, InteractionType iType)
{
    lo_add_angle(fplog, molname,
                 (eitANGLES == iType) ? b->angle : b->linangle,
                 a1, a2, a3, refValue, spacing);
}

static void lo_add_dih(FILE *fplog, const char *molname,
                       std::vector<t_dih> &dih,
                       const std::string a1, const std::string a2, 
                       const std::string a3, const std::string a4,
                       double angle, double spacing, InteractionType iType)

{
    GMX_RELEASE_ASSERT(a1.size() > 0, "atom name a1 is empty");
    GMX_RELEASE_ASSERT(a2.size() > 0, "atom name a2 is empty");
    GMX_RELEASE_ASSERT(a3.size() > 0, "atom name a3 is empty");
    GMX_RELEASE_ASSERT(a4.size() > 0, "atom name a3 is empty");
    
    if (angle < 0)
    {
        angle += 360;
    }
    if (iType == eitIMPROPER_DIHEDRALS)
    {
        while (angle > 176)
        {
            angle -= 180;
        }
    }

    int index = std::lround(angle/spacing);
    if (index < 0)
    {
        index = 0;
    }
    auto d = std::find_if(dih.begin(), dih.end(),
                          [a1, a2, a3, a4](const t_dih &d)
                          {
                              return ((d.a1.compare(a1) == 0 && d.a2.compare(a2) == 0 &&
                                       d.a3.compare(a3) == 0 && d.a4.compare(a4) == 0) ||
                                      (d.a1.compare(a4) == 0 && d.a2.compare(a3) == 0 &&
                                       d.a3.compare(a2) == 0 && d.a4.compare(a1) == 0));
                          });

    if (dih.end() == d)
    {
        t_dih ddd;
        if (a1 < a4)
        {
            ddd.a1.assign(a1);
            ddd.a2.assign(a2);
            ddd.a3.assign(a3);
            ddd.a4.assign(a4);
        }
        else
        {
            ddd.a4.assign(a1);
            ddd.a3.assign(a2);
            ddd.a2.assign(a3);
            ddd.a1.assign(a4);
        }
        if (nullptr != debug)
        {
            fprintf(debug, "NEWDIH  %5s  %5s  %5s  %5s\n",
                    a1.c_str(), a2.c_str(), a3.c_str(), a4.c_str());
        }
        ddd.histo.resize((int) (360/spacing) + 1, 0);
        ddd.lsq = gmx_stats_init();
        dih.push_back(std::move(ddd));
        d = dih.end()-1;
    }
    gmx_stats_add_point(d->lsq, 0, angle, 0, 0);
    d->histo[index]++;
    if (nullptr != fplog)
    {
        fprintf(fplog, "%s dihedral-%s-%s-%s-%s %g\n", molname,
                d->a1.c_str(), d->a2.c_str(), d->a3.c_str(), d->a4.c_str(), angle);
    }
}

static void add_dih(FILE *fplog, const char *molname, t_bonds *b,
                    const std::string a1, const std::string a2, 
                    const std::string a3, const std::string a4,
                    double angle, double spacing, InteractionType iType)
{
    lo_add_dih(fplog, molname,
               (eitPROPER_DIHEDRALS == iType) ? b->dih : b->imp,
               a1, a2, a3, a4, angle, spacing, iType);
}

static void lo_dump_histo(char *fn, char *xaxis, const gmx_output_env_t *oenv, int Nsample,
                          int n, const int histo[], double spacing)
{
    FILE  *fp;
    int    j, j0, j1;
    double sum;

    for (j0 = 0; (j0 < n) && (histo[j0] == 0); j0++)
    {
        ;
    }
    j0 = std::max(j0-1, 0);
    for (j1 = n-1; (j1 > 0) && (histo[j1] == 0); j1--)
    {
        ;
    }
    j1  = std::min(j1+1, n-1);
    sum = 0;
    for (j = j0; (j <= j1); j++)
    {
        sum += histo[j];
    }
    if (sum > 0)
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "%s N = %d", fn, Nsample);
        fp = xvgropen(fn, buf, xaxis, "P (a.u.)", oenv);
        for (j = j0; (j <= j1); j++)
        {
            fprintf(fp, "%g  %g\n", spacing*j, histo[j]/sum);
        }
        fclose(fp);
    }
}

static void dump_histo(t_bonds *b, double bspacing,
                       double aspacing,
                       const gmx_output_env_t *oenv)
{
    int  N;
    char buf[256];

    for (const auto &i : b->bond)
    {
        if ((gmx_stats_get_npoints(i.lsq, &N) == 0) && (i.histo.size() > 0))
        {
            snprintf(buf, sizeof(buf), "bond-%s-%s-%d.xvg", i.a1.c_str(), i.a2.c_str(), i.order);
            lo_dump_histo(buf, (char *)"Distance (pm)", oenv, N,
                          i.histo.size(), i.histo.data(), bspacing);
        }
    }
    for (const auto &i : b->angle)
    {
        if ((gmx_stats_get_npoints(i.lsq, &N) == 0) && (i.histo.size() > 0))
        {
            snprintf(buf, sizeof(buf), "angle-%s-%s-%s.xvg",
                     i.a1.c_str(), i.a2.c_str(), i.a3.c_str());
            lo_dump_histo(buf, (char *)"Angle (deg.)", oenv, N,
                          i.histo.size(), i.histo.data(), aspacing);
        }
    }
    for (const auto &i : b->linangle)
    {
        if ((gmx_stats_get_npoints(i.lsq, &N) == 0) && (i.histo.size() > 0))
        {
            snprintf(buf, sizeof(buf), "linangle-%s-%s-%s.xvg",
                     i.a1.c_str(), i.a2.c_str(), i.a3.c_str());
            lo_dump_histo(buf, (char *)"Linear Angle (deg.)", oenv, N,
                          i.histo.size(), i.histo.data(), aspacing);
        }
    }
    for (const auto &i : b->dih)
    {
        if ((gmx_stats_get_npoints(i.lsq, &N) == 0)  && (i.histo.size() > 0))
        {
            snprintf(buf, sizeof(buf), "dih-%s-%s-%s-%s.xvg",
                     i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), i.a4.c_str());
            lo_dump_histo(buf, (char *)"Dihedral angle (deg.)", oenv, N,
                          i.histo.size(), i.histo.data(), aspacing);
        }
    }
    for (const auto &i : b->imp)
    {
        if ((gmx_stats_get_npoints(i.lsq, &N) == 0)  && (i.histo.size() > 0))
        {
            snprintf(buf, sizeof(buf), "imp-%s-%s-%s-%s.xvg",
                     i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), i.a4.c_str());
            lo_dump_histo(buf, (char *)"Improper angle (deg.)", oenv, N,
                          i.histo.size(), i.histo.data(), aspacing);
        }
    }
}

static void round_numbers(real *av, real *sig)
{
    *av  = ((int)(*av*100))/100.0;
    *sig = ((int)(*sig*100+50))/100.0;
}

static void update_pd(FILE          *fp,      
                      const t_bonds *b,        
                      Poldata       *pd,
                      real           Dm,       
                      real           beta,
                      real           kt, 
                      real           klin,
                      real           kp,
                      real           kimp, 
                      real           kub,
                      real           bond_tol,
                      real           angle_tol)
{
    int                      N;
    real                     av, sig;
    char                     pbuf[256];
    std::vector<std::string> atoms;

    auto                     morse             = pd->findForces(eitBONDS);
    auto                     angle             = pd->findForces(eitANGLES);
    auto                     linear_angle      = pd->findForces(eitLINEAR_ANGLES);
    auto                     proper_dihedral   = pd->findForces(eitPROPER_DIHEDRALS);
    auto                     improper_dihedral = pd->findForces(eitIMPROPER_DIHEDRALS);

    morse->eraseListedForce();
    angle->eraseListedForce();
    linear_angle->eraseListedForce();
    proper_dihedral->eraseListedForce();
    improper_dihedral->eraseListedForce();

    for (auto &i : b->bond)
    {
        gmx_stats_get_average(i.lsq, &av);
        gmx_stats_get_sigma(i.lsq, &sig);
        gmx_stats_get_npoints(i.lsq, &N);
        sprintf(pbuf, "%g  %g", Dm, beta);
        round_numbers(&av, &sig); // Rounding the numbers to 1/10 pm and 1/10 degree
        atoms = {i.a1, i.a2};
        morse->addForce(atoms, pbuf, false, av, sig, N, i.order);

        fprintf(fp, "bond-%s-%s-%d len %g sigma %g (pm) N = %d%s\n",
                i.a1.c_str(), i.a2.c_str(), i.order, av, sig, N,
                (sig > bond_tol) ? " WARNING" : "");
    }
    for (auto &i : b->angle)
    {
        gmx_stats_get_average(i.lsq, &av);
        gmx_stats_get_sigma(i.lsq, &sig);
        gmx_stats_get_npoints(i.lsq, &N);
        round_numbers(&av, &sig);
        atoms = {i.a1, i.a2, i.a3};
        if (angle->fType() == F_ANGLES)
        {   
            sprintf(pbuf, "%g", kt);
        }
        else
        {   
            sprintf(pbuf, "%g  %g", kt, kub);
        }
        angle->addForce(atoms, pbuf, false, av, sig, N);

        fprintf(fp, "harmonic_angle-%s-%s-%s angle %g sigma %g (deg) N = %d%s\n",
                i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), av, sig, N,
                (sig > angle_tol) ? " WARNING" : "");
    }
    for (auto &i : b->linangle)
    {
        gmx_stats_get_average(i.lsq, &av);
        gmx_stats_get_sigma(i.lsq, &sig);
        gmx_stats_get_npoints(i.lsq, &N);
        round_numbers(&av, &sig);
        atoms = {i.a1, i.a2, i.a3};
        sprintf(pbuf, "%g  %g",  klin, kub);
        linear_angle->addForce(atoms, pbuf, false, av, sig, N);

        fprintf(fp, "linear_angle-%s-%s-%s angle %g sigma %g (deg) N = %d%s\n",
                i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), av, sig, N,
                (sig > angle_tol) ? " WARNING" : "");
    }
    for (auto &i : b->dih)
    {
        gmx_stats_get_average(i.lsq, &av);
        gmx_stats_get_sigma(i.lsq, &sig);
        gmx_stats_get_npoints(i.lsq, &N);
        bool support = true;
        switch (proper_dihedral->fType())
        {
        case F_FOURDIHS:
            {
                sprintf(pbuf, "%g -1 1 0", kp);
            }
            break;
        case F_PDIHS:
            {
                sprintf(pbuf, "%g 3", kp);
            }
            break;
        default:
            fprintf(fp, "Unsupported dihedral type %s\n",
                    interaction_function[proper_dihedral->fType()].name);
            support = false;
        }
        if (support)
        {
            round_numbers(&av, &sig);
            atoms = {i.a1, i.a2, i.a3, i.a4};
            proper_dihedral->addForce(atoms, pbuf, false, av, sig, N);
            
            fprintf(fp, "dihedral-%s-%s-%s-%s angle %g sigma %g (deg)\n",
                    i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), i.a4.c_str(), av, sig);
        }
    }
    for (auto &i : b->imp)
    {
        gmx_stats_get_average(i.lsq, &av);
        gmx_stats_get_sigma(i.lsq, &sig);
        gmx_stats_get_npoints(i.lsq, &N);
        sprintf(pbuf, "%g", kimp);
        round_numbers(&av, &sig);
        atoms = {i.a1, i.a2, i.a3, i.a4};
        improper_dihedral->addForce(atoms, pbuf, false, av, sig, N);

        fprintf(fp, "improper-%s-%s-%s-%s angle %g sigma %g (deg)\n",
                i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), i.a4.c_str(), av, sig);
    }
}

int alex_bastat(int argc, char *argv[])
{
    static const char               *desc[] = {
        "bastat read a series of molecules and extracts average geometries from",
        "those. First atomtypes are determined and then bond-lengths, bond-angles",
        "and dihedral angles are extracted. The results are stored in a gentop.dat file."
    };

    t_filenm                         fnm[] = {
        { efDAT, "-f",   "allmols",    ffRDMULT },
        { efDAT, "-d",   "gentop",     ffOPTRD },
        { efDAT, "-o",   "bastat",     ffWRITE },
        { efDAT, "-sel", "molselect",  ffREAD },
        { efLOG, "-g",   "bastat",     ffWRITE }
    };

    const int                        NFILE       = asize(fnm);
    
    static int                       compress    = 0;
    static int                       maxwarn     = 0;
    static real                      Dm          = 0;
    static real                      kt          = 0;
    static real                      kp          = 0;
    static real                      kimp        = 0;
    static real                      beta        = 0;
    static real                      klin        = 0;
    static real                      kub         = 0;
    static real                      bond_tol    = 5;
    static real                      angle_tol   = 5;
    static char                     *lot         = (char *)"B3LYP/aug-cc-pVTZ";
    static gmx_bool                  bHisto      = false;
    static gmx_bool                  bDih        = false;
    static gmx_bool                  bBondOrder  = true;
    t_pargs                          pa[]        = {
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges" },
        { "-maxwarn", FALSE, etINT, {&maxwarn},
          "Will only write output if number of warnings is at most this." },
        { "-Dm",    FALSE, etREAL, {&Dm},
          "Dissociation energy (kJ/mol)" },
        { "-beta",    FALSE, etREAL, {&beta},
          "Steepness of the Morse potential (1/nm)" },
        { "-kt",    FALSE, etREAL, {&kt},
          "Angle force constant (kJ/mol/rad^2)" },
        { "-klin",  FALSE, etREAL, {&klin},
          "Linear angle force constant (kJ/mol/nm^2)" },
        { "-kp",    FALSE, etREAL, {&kp},
          "Dihedral angle force constant (kJ/mol/rad^2)" },
        { "-kimp",    FALSE, etREAL, {&kimp},
          "Improper dihedral angle force constant (kJ/mol/rad^2)" },
        { "-kub",   FALSE, etREAL, {&kub},
          "Urey_Bradley force constant" },
        { "-bond_tol",   FALSE, etREAL, {&bond_tol},
          "Tolerance for bond length" },
        { "-angle_tol",   FALSE, etREAL, {&angle_tol},
          "Tolerance for harmonic and linear angles" },  
        { "-dih",   FALSE, etBOOL, {&bDih},
          "Generate proper dihedral terms" },
        { "-histo", FALSE, etBOOL, {&bHisto},
          "Print (hundreds of) xvg files containing histograms for bonds, angles and dihedrals" },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" },
        { "-bondorder", FALSE, etBOOL, {&bBondOrder},
          "Make separate bonds for different bond orders" }
    };

    FILE                            *fp;
    ChargeModel          iModel;
    time_t                           my_t;
    t_bonds                         *bonds = new(t_bonds);
    rvec                             dx, dx2, r_ij, r_kj, r_kl, mm, nn;
    t_pbc                            pbc;
    int                              t1, t2, t3;
    matrix                           box;
    double                           bspacing = 1;   /* pm */
    double                           aspacing = 0.5; /* degree */
    double                           dspacing = 1;   /* degree */
    gmx_output_env_t                *oenv     = nullptr;
    Poldata                          pd;
    gmx_atomprop_t                   aps;
    MolSelect                        gms;
    std::vector<alexandria::MolProp> mp;
    std::string                      cai, caj, cak, cal;
    
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm, 
                           asize(pa), pa, asize(desc), desc, 
                           0, nullptr, &oenv))
    {
        return 0;
    }
    
    fp                 = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");
    
    time(&my_t);
    fprintf(fp, "# This file was created %s", ctime(&my_t));
    fprintf(fp, "# The Alexandria Chemistry Toolkit.\n#\n");
    
    gms.read(opt2fn("-sel", NFILE, fnm));
    printf("There are %d molecules.\n", (gms.count(imsTrain) + gms.count(imsTest)));
    fprintf(fp, "# There are %d molecules.\n#\n", (gms.count(imsTrain) + gms.count(imsTest)));

    /* Read standard atom properties */
    aps = gmx_atomprop_init();

    /* Read PolData */
    try
    {
        readPoldata(opt2fn_null("-d", NFILE, fnm), pd, aps);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    // This a hack to prevent that no bonds will be found to shells.
    iModel = pd.getChargeModel();
    pd.setChargeModel(eqdESP_p);

    /* Read Molprops */
    auto nwarn = merge_xml(opt2fns("-f", NFILE, fnm), &mp, nullptr, nullptr, nullptr, aps, pd, true);
    
    if (nwarn > maxwarn)
    {
        printf("Too many warnings (%d). Terminating.\n", nwarn);
        return 0;
    }
    for (auto mpi = mp.begin(); mpi < mp.end(); mpi++)
    {
        auto imol = gms.status(mpi->getIupac());
        if (imol == imsTrain || imol == imsTest)
        {
            alexandria::MyMol mmi;
            int               i;
            mmi.molProp()->Merge(mpi);
            if (mmi.molProp()->getMolname().size() == 0)
            {
                printf("Empty molname for molecule with formula %s\n",
                       mmi.molProp()->formula().c_str());
                continue;
            }      
                  
            auto imm = mmi.GenerateTopology(aps, 
                                            &pd, 
                                            lot, 
                                            false, 
                                            false, 
                                            bDih, 
                                            true, 
                                            nullptr);
            if (immOK != imm)
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Could not make topology for %s, reason %s\n",
                            mmi.molProp()->getMolname().c_str(),
                            immsg(imm) );
                }
                continue;
            }
            
#define ATP(ii) (*mmi.atoms_->atomtype[ii])
            
            for (i = 0; i < mmi.atoms_->nr; i++)
            {
                std::string btpi;
                if (!pd.atypeToBtype(*mmi.atoms_->atomtype[i], btpi))
                {
                    if (nullptr != debug)
                    {
                        fprintf(debug, "No bond-type support for atom %s in %s\n",
                                *mmi.atoms_->atomtype[i], mmi.molProp()->getMolname().c_str());
                    }
                    break;
                }
            }
            if ((mmi.atoms_->nr <= 0) || (i < mmi.atoms_->nr))
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "You may need to check the number of atoms for %s\n",
                            mmi.molProp()->getMolname().c_str());
                }
                continue;
            }
            auto x        = mmi.x();
            for (auto fs = pd.forcesBegin(); fs != pd.forcesEnd(); fs++)
            {
                if (eitBONDS == fs->iType())
                {
                    auto funcType = fs->fType();
                    for (auto j = 0; j < mmi.ltop_->idef.il[funcType].nr;
                         j += interaction_function[funcType].nratoms+1)
                    {
                        auto ai = mmi.ltop_->idef.il[funcType].iatoms[j+1];
                        auto aj = mmi.ltop_->idef.il[funcType].iatoms[j+2];
                        rvec_sub(x[ai], x[aj], dx);
                        if (pd.atypeToBtype(*mmi.atoms_->atomtype[ai], cai) &&
                            pd.atypeToBtype(*mmi.atoms_->atomtype[aj], caj))
                        {
                            for (auto bi = mmi.molProp()->BeginBond(); bi < mmi.molProp()->EndBond(); bi++)
                            {
                                auto xi = 0, xj = 0, xb = 0;
                                bi->get(&xi, &xj, &xb);
                                xi--;
                                xj--;
                                if (!bBondOrder)
                                {
                                    xb = 1;
                                }
                                if (((xi == ai) && (xj == aj)) || ((xj == ai) && (xi == aj)))
                                {
                                    add_bond(fp, mmi.molProp()->getMolname().c_str(), bonds,
                                             cai, caj, 1000*norm(dx),
                                             bspacing, xb);
                                    break;
                                }
                            }
                        }
                        else
                        {
                            fprintf(stderr, "No bond_atom type for either %s or %s\n",
                                    ATP(ai), ATP(aj));
                        }
                    }
                }
                else if (eitANGLES == fs->iType() ||
                         eitLINEAR_ANGLES == fs->iType())
                {
                    auto funcType = fs->fType();
                    for (auto j = 0; j < mmi.ltop_->idef.il[funcType].nr;
                         j += interaction_function[funcType].nratoms+1)
                    {
                        auto linear   = false;
                        auto refValue = 0.0;
                        auto ai       = mmi.ltop_->idef.il[funcType].iatoms[j+1];
                        auto aj       = mmi.ltop_->idef.il[funcType].iatoms[j+2];
                        auto ak       = mmi.ltop_->idef.il[funcType].iatoms[j+3];
                        rvec_sub(x[ai], x[aj], dx);
                        rvec_sub(x[ak], x[aj], dx2);
                        refValue = RAD2DEG*gmx_angle(dx, dx2);
                        if ( (refValue > 175) || (refValue < 5))
                        {
                            linear = true;
                        }
                        if (pd.atypeToBtype(*mmi.atoms_->atomtype[ai], cai) &&
                            pd.atypeToBtype(*mmi.atoms_->atomtype[aj], caj) &&
                            pd.atypeToBtype(*mmi.atoms_->atomtype[ak], cak))
                        {
                            add_angle(fp, mmi.molProp()->getMolname().c_str(), bonds,
                                      cai, caj, cak, refValue, aspacing, 
                                      (linear) ? eitLINEAR_ANGLES : eitANGLES);

                            if (nullptr != debug)
                            {
                                fprintf(debug, "Molname: %s  btype1: %s  btype2: %s  btype3: %s  angle: %0.2f\n",
                                        mmi.molProp()->getMolname().c_str(), cai.c_str(), caj.c_str(), cak.c_str(), refValue);
                            }
                        }
                        else
                        {
                            fprintf(stderr, "No bond_atom type for either %s, %s or %s in molecule %s\n",
                                    ATP(ai), ATP(aj), ATP(ak), mmi.molProp()->getMolname().c_str());
                        }
                    }
                }
                else if (eitPROPER_DIHEDRALS   == fs->iType() ||
                         eitIMPROPER_DIHEDRALS == fs->iType())
                {
                    auto funcType = fs->fType();
                    auto angle    = 0.0;
                    for (auto j = 0; j < mmi.ltop_->idef.il[funcType].nr;
                         j += interaction_function[funcType].nratoms+1)
                    {
                        auto ai  = mmi.ltop_->idef.il[funcType].iatoms[j+1];
                        auto aj  = mmi.ltop_->idef.il[funcType].iatoms[j+2];
                        auto ak  = mmi.ltop_->idef.il[funcType].iatoms[j+3];
                        auto al  = mmi.ltop_->idef.il[funcType].iatoms[j+4];
                        angle    = RAD2DEG*dih_angle(x[ai], x[aj], x[ak], x[al],
                                                     &pbc, r_ij, r_kj, r_kl, mm, nn,
                                                     &t1, &t2, &t3);
                        if (pd.atypeToBtype(*mmi.atoms_->atomtype[ai], cai) &&
                            pd.atypeToBtype(*mmi.atoms_->atomtype[aj], caj) &&
                            pd.atypeToBtype(*mmi.atoms_->atomtype[ak], cak) &&
                            pd.atypeToBtype(*mmi.atoms_->atomtype[al], cal))
                        {
                            add_dih(fp, mmi.molProp()->getMolname().c_str(), bonds,
                                    cai, caj, cak, cal, angle, dspacing, fs->iType());
                        }
                        else
                        {
                            fprintf(stderr, "No bond_atom type for either %s, %s, %s or %s\n",
                                    ATP(ai), ATP(aj), ATP(ak), ATP(al));
                        }
                    }
                }
                else
                {
                    fprintf(stderr, "Alexandria does not support the interaction type of %s\n",
                            iType2string(fs->iType()));
                }
                clear_mat(box);
                set_pbc(&pbc, epbcNONE, box);
            }
        }
    }
    sort_bonds(bonds);
    if (bHisto)
    {
        dump_histo(bonds, bspacing, aspacing, oenv);
    }
    update_pd(fp, bonds, &pd,
              Dm, beta, kt, klin, kp, kimp, kub, 
              bond_tol, angle_tol);
    pd.setChargeModel(iModel);
    writePoldata(opt2fn("-o", NFILE, fnm), &pd, compress);
    printf("Extracted %zu bondtypes, %zu angletypes, %zu linear-angletypes, %zu dihedraltypes and %zu impropertypes.\n",
           bonds->bond.size(), bonds->angle.size(), bonds->linangle.size(),
           bonds->dih.size(), bonds->imp.size());
    gmx_ffclose(fp);
    return 0;
}
