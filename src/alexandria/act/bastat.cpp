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
#include "memory_check.h"
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
    std::vector<t_bond>     bond;
    std::vector<t_angle>    angle;
    std::vector<t_angle> linangle;
    std::vector<t_dih>      dih;
    std::vector<t_dih>      imp;
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
                     double blen, double spacing, int order, InteractionType iType)
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
        fprintf(fplog, "%s %s-%s-%s-%d %g\n",
                molname, interactionTypeToString(iType).c_str(),
                a1.c_str(),
                a2.c_str(),
                order,
                blen);
    }
}

static void lo_add_angle(FILE *fplog, const char *molname, std::vector<t_angle> &angle,
                         const std::string a1, const std::string a2, const std::string a3,
                         double refValue, double spacing, InteractionType iType)
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
        fprintf(fplog, "%s %s-%s-%s-%s %g\n", molname, interactionTypeToString(iType).c_str(),
                a1.c_str(), a2.c_str(), a3.c_str(), refValue);
    }
}

static void add_angle(FILE *fplog, const char *molname, t_bonds *b,
                      const std::string a1, const std::string a2, const std::string a3,
                      double refValue, double spacing, InteractionType iType)
{
    lo_add_angle(fplog, molname, iType == InteractionType::ANGLES ? b->angle : b->linangle,
                 a1, a2, a3, refValue, spacing, iType);
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
    if (iType == InteractionType::IMPROPER_DIHEDRALS)
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
        fprintf(fplog, "%s %s-%s-%s-%s-%s %g\n", molname, interactionTypeToString(iType).c_str(),
                d->a1.c_str(), d->a2.c_str(), d->a3.c_str(), d->a4.c_str(), angle);
    }
}

static void add_dih(FILE *fplog, const char *molname, t_bonds *b,
                    const std::string a1, const std::string a2,
                    const std::string a3, const std::string a4,
                    double angle, double spacing, InteractionType iType)
{
    lo_add_dih(fplog, molname,
               (InteractionType::PROPER_DIHEDRALS == iType) ? b->dih : b->imp,
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

static void round_numbers(real *av, real *sig, int power10)
{
    *av  = ((int)(*av*power10))/(1.0*power10);
    *sig = ((int)(*sig*1.5*power10))/(1.0*power10);
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
                      real           angle_tol,
                      real           factor)
{
    std::vector<InteractionType> myIt =
        { 
            InteractionType::BONDS,
            InteractionType::ANGLES,
            InteractionType::LINEAR_ANGLES,
            InteractionType::PROPER_DIHEDRALS,
            InteractionType::IMPROPER_DIHEDRALS 
        };

    for(auto &iType : myIt)
    {
        auto fs    = pd->findForces(iType);
        auto fType = fs->fType();
        fs->eraseParameter();
    
        switch (iType)
        {
        case InteractionType::BONDS:
            // Note that the order of parameters is important!
            for (auto &i : b->bond)
            {
                int  N;
                real av, sig;
                gmx_stats_get_average(i.lsq, &av);
                gmx_stats_get_sigma(i.lsq, &sig);
                gmx_stats_get_npoints(i.lsq, &N);
                round_numbers(&av, &sig, 10); // Rounding the numbers to 1/10 pm and 1/10 degree
                Identifier bondId({i.a1, i.a2}, i.order, CanSwap::Yes);
                fs->addParameter(bondId, "bondlength",
                                 ForceFieldParameter("pm", av, sig, N, av*factor, av/factor, Mutability::Bounded, false));
                fs->addParameter(bondId, "Dm",
                                 ForceFieldParameter("kJ/mol", Dm, 0, 1, Dm*factor, Dm/factor, Mutability::Bounded, false));
                fs->addParameter(bondId, "beta",
                                 ForceFieldParameter("1/nm", beta, 0, 1, beta*factor, beta/factor, Mutability::Bounded, false));
        
                fprintf(fp, "bond-%s len %g sigma %g (pm) N = %d%s\n",
                        bondId.id().c_str(), av, sig, N, (sig > bond_tol) ? " WARNING" : "");
            }
        break;
        case InteractionType::ANGLES:
            for (auto &i : b->angle)
            {
                int  N;
                real av, sig;
                gmx_stats_get_average(i.lsq, &av);
                gmx_stats_get_sigma(i.lsq, &sig);
                gmx_stats_get_npoints(i.lsq, &N);
                round_numbers(&av, &sig, 10);
                Identifier bondId({i.a1, i.a2, i.a3}, CanSwap::Yes);
                fs->addParameter(bondId, "angle",
                                 ForceFieldParameter("degree", av, sig, N, av*factor, av/factor, Mutability::Bounded, false));
                fs->addParameter(bondId, "kt",
                                 ForceFieldParameter("kJ/mol/rad2", kt, 0, 1, kt*factor, kt/factor, Mutability::Bounded, false));
                if (fType == F_UREY_BRADLEY)
                {
                    fs->addParameter(bondId, "r13", 
                                     ForceFieldParameter("nm", 0, 0, 1, 0, 0, Mutability::Dependent, false));
                    fs->addParameter(bondId, "kub", 
                                     ForceFieldParameter("kJ/mol/nm2", kub, 0, 1, kub*factor, kub/factor, Mutability::Bounded, false));
                }
                
                fprintf(fp, "angle-%s angle %g sigma %g (deg) N = %d%s\n",
                        bondId.id().c_str(), av, sig, N, (sig > angle_tol) ? " WARNING" : "");
            }
            break;
        case InteractionType::LINEAR_ANGLES:
            for (auto &i : b->linangle)
            {
                int  N;
                real av, sig;
                gmx_stats_get_average(i.lsq, &av);
                gmx_stats_get_sigma(i.lsq, &sig);
                gmx_stats_get_npoints(i.lsq, &N);
                round_numbers(&av, &sig, 1000000);
                Identifier bondId({i.a1, i.a2, i.a3}, CanSwap::No);
                // TODO Fix the parameters to be correct!
                double myfactor = 0.99;
                fs->addParameter(bondId, "a",
                                 ForceFieldParameter("", av, sig, N, av*myfactor, av/myfactor, Mutability::Bounded, false));
                fs->addParameter(bondId, "klin", 
                                 ForceFieldParameter("kJ/mol/nm2", klin, 0, 1, av*factor, av/factor, Mutability::Bounded, false));
                
                fprintf(fp, "linear_angle-%s angle %g sigma %g N = %d%s\n",
                        bondId.id().c_str(), av, sig, N, (sig > angle_tol) ? " WARNING" : "");
            }
            break;
        case InteractionType::PROPER_DIHEDRALS:
            for (auto &i : b->dih)
            {
                int  N;
                real av, sig;
                Identifier bondId({i.a1, i.a2, i.a3, i.a4}, CanSwap::Yes);
                 
                switch (fType)
                {
                case F_FOURDIHS:
                    {
                        double val = 1;
                        std::vector<std::string> cname = { "c0", "c1", "c2", "c3" };
                        for(auto &c : cname)
                        {
                            fs->addParameter(bondId, c,
                                             ForceFieldParameter("kJ/mol", val, 0, 1, val*factor, val/factor, Mutability::Bounded, false));
                        }
                    }
                    break;
                case F_PDIHS:
                    {
                        gmx_stats_get_average(i.lsq, &av);
                        gmx_stats_get_sigma(i.lsq, &sig);
                        gmx_stats_get_npoints(i.lsq, &N);
                        round_numbers(&av, &sig, 10);
                        fs->addParameter(bondId, "angle", 
                                         ForceFieldParameter("degree", av, sig, N, av*factor, av/factor, Mutability::Bounded, false));
                        fs->addParameter(bondId, "kp", 
                                         ForceFieldParameter("kJ/mol", kp, 0, 1, kp*factor, kp/factor, Mutability::Bounded, false));
                        fs->addParameter(bondId, "mult", 
                                         ForceFieldParameter("", 3, 0, 1, 3, 3, Mutability::Fixed, true));
                    }
                    break;
                default:
                    GMX_THROW(gmx::InternalError(gmx::formatString("Unsupported dihedral type %s",
                                                                   interaction_function[fType].name).c_str()));
                }
                fprintf(fp, "dihedral-%s angle %g sigma %g (deg)\n",
                        bondId.id().c_str(), av, sig);
            }
            break;
        case InteractionType::IMPROPER_DIHEDRALS:
            for (auto &i : b->imp)
            {
                int  N;
                real av, sig;
                gmx_stats_get_average(i.lsq, &av);
                gmx_stats_get_sigma(i.lsq, &sig);
                gmx_stats_get_npoints(i.lsq, &N);
                round_numbers(&av, &sig, 10);
                Identifier bondId({i.a1, i.a2, i.a3, i.a4}, CanSwap::No);

                fs->addParameter(bondId, "phi", 
                                 ForceFieldParameter("degree", av, sig, N, av*factor, av/factor, Mutability::Bounded, false));
                fs->addParameter(bondId, "kimp", 
                                 ForceFieldParameter("kJ/mol", kimp, 0, 1, kimp*factor, kimp/factor, Mutability::Bounded, false));
                
                fprintf(fp, "improper-%s angle %g sigma %g (deg)\n",
                        bondId.id().c_str(), av, sig);
            }
            break;
        default:
            break;
        }
    }
}

static void generate_bcc(Poldata *pd,
                         double   hardness)
{
    // Bonds should be present, so no checking
    auto bonds = pd->findForcesConst(InteractionType::BONDS);
    auto itype = InteractionType::BONDCORRECTIONS;
    if (!pd->interactionPresent(itype))
    {
        auto canSwap = CanSwap::No;
        ForceFieldParameterList newparam("", canSwap);
        pd->addForces(interactionTypeToString(itype), newparam);
    }
    auto bcc   = pd->findForces(itype);
    bcc->clearParameters();

    auto hardnessParam = ForceFieldParameter("eV/e", hardness, 0, 0, 0, 2, Mutability::Bounded, true);
    auto enpBounded    = ForceFieldParameter("eV", 0, 0, 0, -2, 2, Mutability::Bounded, true);
    auto enpFixed      = ForceFieldParameter("eV", 0, 0, 0, 0, 0, Mutability::Fixed, true);
    auto ptypes = pd->particleTypesConst();
    for (auto &ai : ptypes)
    {
        for (auto &aj : ptypes)
        {
            for(int bondOrder = 1; bondOrder < 4; bondOrder++)
            {
                auto bi = ai.interactionTypeToIdentifier(InteractionType::BONDS).id();
                auto bj = aj.interactionTypeToIdentifier(InteractionType::BONDS).id();
                Identifier bondId({ bi, bj }, bondOrder, bonds.canSwap());
                if (bonds.parameterExists(bondId))
                {
                    auto entype = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
                    auto zi = ai.interactionTypeToIdentifier(entype).id();
                    auto zj = aj.interactionTypeToIdentifier(entype).id();
                    if (!zi.empty() && !zj.empty())
                    {
                        Identifier bccId1({ zi, zj }, bondOrder, bcc->canSwap());
                        Identifier bccId2({ zj, zi }, bondOrder, bcc->canSwap());
                        if (!bcc->parameterExists(bccId1) && 
                            !bcc->parameterExists(bccId2))
                        {
                            if (zi == zj)
                            {
                                bcc->addParameter(bccId1, "electronegativity", enpFixed);
                            }
                            else
                            {
                                bcc->addParameter(bccId1, "electronegativity", enpBounded);
                            }
                            bcc->addParameter(bccId1, "hardness", hardnessParam);
                        }
                    }
                }
            }
        }
    }
    printf("Have generated %zu entries for BCC\n", bcc->parameters()->size());
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
    static real                      hardness    = 1;
    static real                      kub         = 0;
    static real                      bond_tol    = 5;
    static real                      angle_tol   = 5;
    static real                      factor      = 0.8;
    static char                     *lot         = (char *)"B3LYP/aug-cc-pVTZ";
    static gmx_bool                  bHisto      = false;
    static gmx_bool                  bDih        = false;
    static gmx_bool                  bBondOrder  = true;
    static gmx_bool                  genBCC      = true;
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
          "Make separate bonds for different bond orders" },
        { "-genBCC", FALSE, etBOOL, {&genBCC},
          "Re-generate bond charge corrections based on the list of bonds" },
        { "-hardness", FALSE, etREAL, {&hardness},
          "Default bond hardness when generating bond charge corrections based on the list of bonds" },
        { "-factor", FALSE, etBOOL, {&factor},
          "Scale factor to set minimum and maximum values of parameters" }
    };

    FILE                            *fp;
    time_t                           my_t;
    t_bonds                         *bonds = new(t_bonds);
    rvec                             dx, dx2, r_ij, r_kj, r_kl, mm, nn;
    t_pbc                            pbc;
    int                              t1, t2, t3;
    matrix                           box;
    double                           bspacing  = 1;   /* pm */
    double                           aspacing  = 0.5; /* degree */
    double                           laspacing = 0.000001; /* relative number for linear angles */
    double                           dspacing  = 1;   /* degree */
    gmx_output_env_t                *oenv      = nullptr;
    Poldata                          pd;
    gmx_atomprop_t                   aps;
    MolSelect                        gms;
    std::vector<alexandria::MolProp> mp;
    std::string                      cai, caj, cak, cal;
    std::string                      method, basis;
    splitLot(lot, &method, &basis);
    
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm,
                           asize(pa), pa, asize(desc), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }

    fp                 = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");
    print_memory_usage(fp);
    time(&my_t);
    fprintf(fp, "# This file was created %s", ctime(&my_t));
    fprintf(fp, "# The Alexandria Chemistry Toolkit.\n#\n");

    auto selfile = opt2fn("-sel", NFILE, fnm);
    gms.read(selfile);
    print_memory_usage(fp);
    printf("There are %d molecules in the selection file %s.\n",
           (gms.count(imsTrain) + gms.count(imsTest)), selfile);
    fprintf(fp, "# There are %d molecules.\n#\n", (gms.count(imsTrain) + gms.count(imsTest)));

    /* Read standard atom properties */
    aps = gmx_atomprop_init();
    print_memory_usage(fp);

    /* Read PolData */
    try
    {
        readPoldata(opt2fn_null("-d", NFILE, fnm), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    print_memory_usage(fp);

    // This a hack to prevent that no bonds will be found to shells.
    bool polar = pd.polarizable();
    pd.setPolarizable(false);

    /* Read Molprops */
    auto nwarn = merge_xml(opt2fns("-f", NFILE, fnm), &mp, nullptr, nullptr, nullptr, aps, true);
    print_memory_usage(fp);

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
            mmi.Merge(&(*mpi));
            if (mmi.getMolname().size() == 0)
            {
                printf("Empty molname for molecule with formula %s\n",
                       mmi.formula().c_str());
                continue;
            }
            std::string mylot;
            auto        imm = mmi.GenerateTopology(&pd,
                                                   method,
                                                   basis,
                                                   &mylot,
                                                   false,
                                                   false,
                                                   bDih,
                                                   missingParameters::Generate,
                                                   nullptr);
            if (immStatus::OK != imm)
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Could not make topology for %s, reason %s\n",
                            mmi.getMolname().c_str(),
                            immsg(imm) );
                }
                continue;
            }

            auto myatoms = mmi.atomsConst();
#define ATP(ii) (*myatoms.atomtype[ii])
            for (i = 0; i < myatoms.nr; i++)
            {
                std::string btpi;
                if (!pd.atypeToBtype(*myatoms.atomtype[i], &btpi))
                {
                    if (nullptr != debug)
                    {
                        fprintf(debug, "No bond-type support for atom %s in %s\n",
                                *myatoms.atomtype[i], mmi.getMolname().c_str());
                    }
                    break;
                }
            }
            if ((myatoms.nr <= 0) || (i < myatoms.nr))
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "You may need to check the number of atoms for %s\n",
                            mmi.getMolname().c_str());
                }
                continue;
            }
            auto x        = mmi.x();
            for (auto &fs : pd.forcesConst())
            {
                auto iType    = fs.first;
                auto funcType = fs.second.fType();
                switch (iType)
                {
                case InteractionType::BONDS:
                    {
                        for (auto j = 0; j < mmi.ltop_->idef.il[funcType].nr;
                             j += interaction_function[funcType].nratoms+1)
                        {
                            auto ai = mmi.ltop_->idef.il[funcType].iatoms[j+1];
                            auto aj = mmi.ltop_->idef.il[funcType].iatoms[j+2];
                            rvec_sub(x[ai], x[aj], dx);
                            if (pd.atypeToBtype(*myatoms.atomtype[ai], &cai) &&
                                pd.atypeToBtype(*myatoms.atomtype[aj], &caj))
                            {
                                for (auto &bi : mmi.bondsConst())
                                {
                                    auto xi = 0, xj = 0, xb = 0;
                                    bi.get(&xi, &xj, &xb);
                                    xi--;
                                    xj--;
                                    if (!bBondOrder)
                                    {
                                        xb = 1;
                                    }
                                    if (((xi == ai) && (xj == aj)) || ((xj == ai) && (xi == aj)))
                                    {
                                        add_bond(fp, mmi.getMolname().c_str(), bonds,
                                                 cai, caj, 1000*norm(dx),
                                                 bspacing, xb, iType);
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
                    break;
                case InteractionType::ANGLES:
                case InteractionType::LINEAR_ANGLES:
                    {
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
                            auto spacing  = aspacing;
                            if ( (refValue > 175) || (refValue < 5))
                            {
                                linear   = true;
                                refValue = norm(dx2)/(norm(dx)+norm(dx2));
                                spacing  = laspacing;
                            }
                            if (pd.atypeToBtype(*myatoms.atomtype[ai], &cai) &&
                                pd.atypeToBtype(*myatoms.atomtype[aj], &caj) &&
                                pd.atypeToBtype(*myatoms.atomtype[ak], &cak))
                            {
                                add_angle(fp, mmi.getMolname().c_str(), bonds,
                                          cai, caj, cak, refValue, spacing,
                                          (linear) ? InteractionType::LINEAR_ANGLES : InteractionType::ANGLES);
                                
                                if (nullptr != debug)
                                {
                                    fprintf(debug, "Molname: %s  btype1: %s  btype2: %s  btype3: %s  angle: %0.2f\n",
                                            mmi.getMolname().c_str(), cai.c_str(), caj.c_str(), cak.c_str(), refValue);
                                }
                            }
                            else
                            {
                                fprintf(stderr, "No bond_atom type for either %s, %s or %s in molecule %s\n",
                                        ATP(ai), ATP(aj), ATP(ak), mmi.getMolname().c_str());
                            }
                        }
                    }
                    break;
                case InteractionType::PROPER_DIHEDRALS:
                case InteractionType::IMPROPER_DIHEDRALS:
                    {
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
                            if (pd.atypeToBtype(*myatoms.atomtype[ai], &cai) &&
                                pd.atypeToBtype(*myatoms.atomtype[aj], &caj) &&
                                pd.atypeToBtype(*myatoms.atomtype[ak], &cak) &&
                                pd.atypeToBtype(*myatoms.atomtype[al], &cal))
                            {
                                add_dih(fp, mmi.getMolname().c_str(), bonds,
                                        cai, caj, cak, cal, angle, dspacing, iType);
                            }
                            else
                            {
                                fprintf(stderr, "No bond_atom type for either %s, %s, %s or %s\n",
                                        ATP(ai), ATP(aj), ATP(ak), ATP(al));
                            }
                        }
                    }
                    break;
                default:
                    if (nullptr != debug)
                    {
                        fprintf(debug, "Alexandria does not support the interaction type of %s\n",
                                interactionTypeToString(iType).c_str());
                    }
                }
                clear_mat(box);
                set_pbc(&pbc, epbcNONE, box);
            }
        }
    }
    sort_bonds(bonds);
    print_memory_usage(fp);
    if (bHisto)
    {
        dump_histo(bonds, bspacing, aspacing, oenv);
    }
    update_pd(fp, bonds, &pd,
              Dm, beta, kt, klin, kp, kimp, kub,
              bond_tol, angle_tol, factor);
    pd.setPolarizable(polar);
    if (genBCC)
    {
        generate_bcc(&pd, hardness);
    }
    print_memory_usage(fp);
    writePoldata(opt2fn("-o", NFILE, fnm), &pd, compress);
    printf("Extracted %zu bondtypes, %zu angletypes, %zu linear-angletypes, %zu dihedraltypes and %zu impropertypes.\n",
           bonds->bond.size(), bonds->angle.size(), bonds->linangle.size(),
           bonds->dih.size(), bonds->imp.size());
    print_memory_usage(fp);
    gmx_ffclose(fp);
    return 0;
}
