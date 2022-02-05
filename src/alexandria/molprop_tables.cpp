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

#include "molprop_tables.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/strconvert.h"

#include "alex_modules.h"
#include "categories.h"
#include "utility/latex_util.h"

namespace alexandria
{

class ExpData
{
    public:
        double      val_, err_, temp_;
        std::string ref_, conf_;

        ExpData(double val, double err, double temp, 
                std::string ref, std::string conf) : 
            val_(val),
            err_(err), 
            temp_(temp), 
            ref_(ref), 
            conf_(conf) {}
};

class CalcData
{
    public:
        double val_, err_, temp_;
        int    found_;
        CalcData(double val, double err, double temp, int found) :
            val_(val), 
            err_(err), 
            temp_(temp), 
            found_(found) {};
};

typedef struct {
    char       *ptype;
    char       *miller;
    char       *bosque;
    gmx_stats   lsq;
    int         nexp;
    int         nqm;
} t_sm_lsq;

static void stats_header(LongTable         &lt,
                         MolPropObservable  mpo,
                         const QmCount     &qmc,
                         iMolSelect         ims)
{
    lt.setColumns(1+qmc.nCalc());

    for (int k = 0; (k < 2); k++)
    {
        if (0 == k)
        {
            // Caption
            char caption[STRLEN];
            snprintf(caption, sizeof(caption), "Performance of the different methods for predicting the molecular %s for molecules containing different chemical groups, given as the RMSD from experimental values (%s), and in brackets the number of molecules in this particular subset. {\\bf Data set: %s.} At the bottom the correlation coefficient R, the regression coefficient a and the intercept b are given as well as the normalized quality of the fit $\\chi^2$, the mean signed error (MSE) and the mean absolute error (MSA).",
                     mpo_name(mpo), mpo_unit(mpo), iMolSelectName(ims));
            lt.setCaption(caption);
            // Label
            char label[STRLEN];
            snprintf(label, sizeof(label), "%s_rmsd", mpo_name(mpo));
            lt.setLabel(label);
        }
        else
        {
            std::string hline("Method ");
            char        koko[STRLEN];
            for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q)
            {
                snprintf(koko, STRLEN, " & %s", q->method().c_str());
                hline.append(koko);
            }
            lt.addHeadLine(hline);

            hline.assign(" ");
            for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q)
            {
                snprintf(koko, STRLEN, "& %s ", q->basis().c_str());
                hline.append(koko);
            }
            lt.addHeadLine(hline);

            hline.assign(" ");
            for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q)
            {
                snprintf(koko, STRLEN, "& %s ", q->type().c_str());
                hline.append(koko);
            }
            lt.addHeadLine(hline);
        }
    }
    lt.printHeader();
}

void alexandria_molprop_stats_table(FILE                 *fp,
                                    MolPropObservable     mpo,
                                    std::vector<MolProp> &mp,
                                    const QmCount        &qmc,
                                    double                outlier,
                                    CategoryList          cList,
                                    const MolSelect      &gms,
                                    iMolSelect            ims)
{
    std::vector<MolProp>::iterator     mpi;
    std::vector<std::string>::iterator si;
    real                               rms, R, a, da, b, db, chi2;
    char                               buf[256];
    gmx_stats                          lsq;
    std::vector<gmx_stats>             lsqtot;
    LongTable                          lt(fp, true, nullptr);

    if (0 == cList.nCategories())
    {
        fprintf(stderr, "No categories. cList not initialized? Not doing category statistics.\n");
        return;
    }

    stats_header(lt, mpo, qmc, ims);

    for (auto &i : cList.elementsConst())
    {
        std::string catbuf;
        int         nqmres  = 0;
        int         nexpres = 0;
        std::string mybuf = i.getName();
        snprintf(buf, sizeof(buf), "%s", mybuf.c_str());
        catbuf.append(buf);
        for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q)
        {
            for (auto &mpi : mp)
            {
                if ((i.hasMolecule(mpi.getIupac())) &&
                    (mpi.SearchCategory(i.getName()) == 1))
                {
                    double Texp  = -1;
                    auto   gpexp = mpi.findProperty(mpo, iqmType::Exp, Texp, "", "", "");
                    if (gpexp)
                    {
                        double exp_val = gpexp->getValue();
                        double exp_err = gpexp->getError();
                        double Tqm     = -1;
                        auto   gpqm    = mpi.findProperty(mpo, iqmType::QM, Tqm, q->method(), q->basis(), "");
                        if (gpqm)
                        {
                            double qm_val = gpqm->getValue();
                            double qm_err = gpqm->getValue();
                            if (debug)
                            {
                                fprintf(debug, "%s %s - TAB4\n",
                                        mpi.getMolname().c_str(),
                                        i.getName().c_str());
                            }
                            lsq.add_point(exp_val, qm_val, exp_err, qm_err);
                            nexpres = 1;
                        }
                    }
                }
            }
            if (outlier > 0)
            {
                lsq.remove_outliers(outlier);
            }
            if (lsq.get_rmsd(&rms) == eStats::OK)
            {
                int N = lsq.get_npoints();
                snprintf(buf, sizeof(buf)-1, "& %8.1f(%d)", rms, N);
                catbuf.append(buf);
                nqmres++;
            }
            else
            {
                catbuf.append("& -");
            }
        }
        if ((nqmres > 0) && (nexpres > 0))
        {
            lt.printLine(catbuf);
        }
    }
    std::string catbuf;
    catbuf.assign("All");
    lsqtot.resize(qmc.nCalc());
    int k = 0;
    for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q, ++k)
    {
        for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
        {
            iMolSelect myIms;
            if (gms.status(mpi->getIupac(), &myIms) &&
                ims == myIms &&
                mpi->hasAllAtomTypes())
            {
                double Texp = -1;
                auto gpexp = mpi->findProperty(mpo, iqmType::Exp, Texp, "", "", "");
                if (gpexp)
                {
                    double Tqm     = -1;
                    double exp_err = gpexp->getError();
                    double exp_val = gpexp->getValue();
                    auto   gpqm = mpi->findProperty(mpo, iqmType::QM, Tqm, "", "", "");
                    if (gpqm)
                    {
                        double qm_err = gpqm->getError();
                        double qm_val = gpqm->getValue();
                        lsqtot[k].add_point(exp_val, qm_val, exp_err, qm_err);
                    }
                }
            }
        }
        if (lsqtot[k].get_rmsd(&rms) == eStats::OK)
        {
            int N = lsqtot[k].get_npoints();
            snprintf(buf, sizeof(buf),  "& %8.1f(%d)", rms, N);
            catbuf.append(buf);
        }
        else
        {
            catbuf.append("& -");
        }
    }
    lt.printLine(catbuf);
    lt.printHLine();

    catbuf.assign("a");
    for (auto &k : lsqtot)
    {
        if (k.get_ab(elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &R) ==
            eStats::OK)
        {
            snprintf(buf, sizeof(buf), "& %8.2f(%4.2f)", a, da);
            catbuf.append(buf);
        }
        else
        {
            catbuf.append("& -");
        }
    }
    lt.printLine(catbuf);

    catbuf.assign("b");
    for (auto &k : lsqtot)
    {
        if (k.get_ab(elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &R) ==
            eStats::OK)
        {
            snprintf(buf, sizeof(buf), "& %8.2f(%4.2f)", b, db);
            catbuf.append(buf);
        }
        else
        {
            catbuf.append("& -");
        }
    }
    lt.printLine(catbuf);

    catbuf.assign("R$^2$ (\\%)");
    for (auto &k : lsqtot)
    {
        if (k.get_corr_coeff(&R) == eStats::OK)
        {
            snprintf(buf, sizeof(buf), "& %8.2f", 100*R*R);
            catbuf.append(buf);
        }
        else
        {
            catbuf.append("& -");
        }
    }
    lt.printLine(catbuf);

    catbuf.assign("$\\chi^2$");
    for (auto &k : lsqtot)
    {
        if (k.get_ab(elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &R) ==
            eStats::OK)
        {
            snprintf(buf, sizeof(buf), "& %8.2f", chi2);
            catbuf.append(buf);
        }
        else
        {
            catbuf.append("& -");
        }
    }
    lt.printLine(catbuf);

    catbuf.assign("MSE");
    for (auto &k : lsqtot)
    {
        real mse;
        if (k.get_mse_mae(&mse, nullptr) ==
            eStats::OK)
        {
            snprintf(buf, sizeof(buf), "& %8.2f", mse);
            catbuf.append(buf);
        }
        else
        {
            catbuf.append("& -");
        }
    }
    lt.printLine(catbuf);

    catbuf.assign("MAE");
    for (auto &k : lsqtot)
    {
        real mae;
        if (k.get_mse_mae(nullptr, &mae) ==
            eStats::OK)
        {
            snprintf(buf, sizeof(buf), "& %8.2f", mae);
            catbuf.append(buf);
        }
        else
        {
            catbuf.append("& -");
        }
    }
    lt.printLine(catbuf);
    lt.printFooter();
}

#ifdef OLD
static void composition_header(LongTable             &lt,
                               iMolSelect             ims)
{
    char caption[STRLEN];

    snprintf(caption, STRLEN, "Decomposition of molecules into Alexandria atom types. {\\bf Data set: %s.} Charge is given when not zero, multiplicity is given when not 1.",
             iMolSelectName(ims));
    lt.setCaption(caption);
    lt.setLabel("frag_defs");
    lt.setColumns("p{75mm}ll");
    lt.addHeadLine("Molecule & Formula  & Types");
    lt.printHeader();
}

void alexandria_molprop_composition_table(FILE                 *fp,
                                          std::vector<MolProp>  mp,
                                          const MolSelect      &gms,
                                          iMolSelect            ims)
{
    std::vector<MolProp>::iterator mpi;
    MolecularCompositionIterator   mci;
    int                            q, m, iline, nprint;
    LongTable                      lt(fp, true, "small");

    nprint = 0;
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        iMolSelect myIms;
        if (gms.status(mpi->getIupac(), &myIms) &&
            myIms == ims &&
            mpi->hasAtomTypes())
        {
            nprint++;
        }
    }
    if (nprint <= 0)
    {
        return;
    }

    composition_header(lt, ims);
    iline = 0;
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        iMolSelect myIms;
        if (gms.status(mpi->getIupac(), &myIms) &&
            myIms == ims &&
            mpi->HasComposition(alex))
        {
            char qbuf[32];
            q = mpi->totalCharge();
            m = mpi->getMultiplicity();
            if ((q != 0) || (m != 1))
            {
                if ((q != 0) && (m != 1))
                {
                    snprintf(qbuf, sizeof(qbuf), " (q=%c%d, mult=%d)", (q < 0) ? '-' : '+', abs(q), m);
                }
                else if (q != 0)
                {
                    snprintf(qbuf, sizeof(qbuf), " (q=%c%d)", (q < 0) ? '-' : '+', abs(q));
                }
                else
                {
                    snprintf(qbuf, sizeof(qbuf), " (mult=%d)", m);
                }
            }
            else
            {
                qbuf[0] = '\0';
            }
            char buf[STRLEN];
            snprintf(buf, STRLEN, "%3d. %s%s & %s & ",
                     ++iline,
                     mpi->getIupac().c_str(),
                     qbuf,
                     mpi->getTexFormula().c_str());
            std::string longbuf(buf);
            mci = mpi->SearchMolecularComposition(cs.searchCS(iCalexandria)->name());
            if (mci != mpi->EndMolecularComposition())
            {
                AtomNumIterator ani;

                for (auto ani : mci->atomNumConst())
                {
                    snprintf(buf, 256, " %d %s\t",
                             ani.getNumber(), ani.getAtom().c_str());
                    longbuf.append(buf);
                }
            }
            lt.printLine(longbuf);
        }
    }
    lt.printFooter();
}
#endif

static void category_header(LongTable &lt)
{
    char longbuf[STRLEN];

    lt.setColumns("lcp{150mm}");
    snprintf(longbuf, sizeof(longbuf),
             "Molecules that are part of each category used for statistics.");
    lt.setCaption(longbuf);
    lt.setLabel("stats");
    lt.addHeadLine("Category & N & Molecule(s)");
    lt.printHeader();
}

void alexandria_molprop_category_table(FILE            *fp,
                                       int              catmin,
                                       CategoryList     cList)
{
    if (cList.nCategories() > 0)
    {
        LongTable lt(fp, true, "small");
        category_header(lt);
        for (const auto &i : cList.elementsConst())
        {
            std::string longbuf;
            char        buf[256];
            int         nMol = i.nMolecule();

            if ((nMol >= catmin) && (catmin > 1))
            {
                int n = 0;
                const std::vector<std::string> &mols = i.molecules();
                snprintf(buf, sizeof(buf), "%s & %d &", i.getName().c_str(), nMol);
                longbuf.append(buf);
                for (size_t j = 0; j < mols.size()-1; j++)
                {
                    snprintf(buf, sizeof(buf), "%s, ", mols[j].c_str());
                    longbuf.append(buf);
                    n++;
                    if (0 == (n % 50))
                    {
                        lt.printLine(longbuf);
                        longbuf.assign(" & &");
                    }
                }
                snprintf(buf, sizeof(buf), "%s", mols.back().c_str());
                longbuf.append(buf);
                lt.printLine(longbuf);
            }
        }
        lt.printFooter();
    }
}

#ifdef IGNORE
static void atomtype_tab_header(LongTable &lt)
{
    char             longbuf[STRLEN];
    CompositionSpecs cs;

    lt.setColumns("lcccccc");

    snprintf(longbuf, STRLEN, "Atomic polarizability (\\AA$^3$) obtained from the decomposition of the experimental isotropic molecular polarizability. $N$ is the number of experimental datapoints used. The columns Ahc and Ahp contain atomic hybrid components~\\protect\\cite{Miller1979a} and atomic hybrid polarizabilites~\\protect\\cite{Miller1990a, Kang1982a}, respectively. The column BS contains the polarizabilities of Bosque and Sales~\\protect\\cite{Bosque2002a}. The uncertainty, $\\sigma$, in the Alexandria polarizability values are computed by Bootstrapping with 1000 interations.");
    lt.setCaption(longbuf);
    lt.setLabel("fragments");
    snprintf(longbuf, STRLEN, "Polarizability Type  & $N$ & \\multicolumn{4}{c}{Polarizability}");
    lt.addHeadLine(longbuf);
    snprintf(longbuf, STRLEN, "& & %s ($\\sigma$) & Ahc & Ahp & %s ",
             cs.searchCS(iCalexandria)->name(),
             cs.searchCS(iCbosque)->abbreviation());
    lt.addHeadLine(longbuf);
    lt.printHeader();
}

static void alexandria_molprop_atomtype_polar_table(FILE                       *fp,
                                                    const Poldata              *pd,
                                                    std::vector<MolProp>        mp,
                                                    const char                 *lot)
{
    std::vector<MolProp>::iterator  mpi;
    double                          ahc, ahp, bos_pol;
    char                            longbuf[STRLEN];
    MolPropObservable               mpo = MolPropObservable::POLARIZABILITY;
    LongTable                       lt(fp, false, nullptr);
    CompositionSpecs                cs;
    const char                     *alexandria = cs.searchCS(iCalexandria)->name();

    /* Prepare printing it! */
    atomtype_tab_header(lt);

    /* First gather statistics from different input files.
     * Note that the input files
     * do not need to have the same sets of types,
     * as we check for the type name.
     */
    for (auto pType : pd->particleTypesConst())
    {
        if (pType.getPolarizability() > 0)
        {
            int atomnumber;
            int nexp   = 0;
            int nqm    = 0;

            /* Now loop over the poltypes and count number of occurrences
             * in molecules
             */
            for (auto &mpi : mp)
            {
                auto mci = mpi.SearchMolecularComposition(alexandria);

                if (mci != mpi.EndMolecularComposition())
                {
                    bool bFound = false;
                    for (auto ani = mci->BeginAtomNum(); !bFound && (ani < mci->EndAtomNum()); ++ani)
                    {
                        std::string pt;
                        if (pd->atypeToPtype(ani->getAtom(), &pt))
                        {
                            if (pt == pType->getType())
                            {
                                bFound = true;
                            }
                        }
                    }
                    if (bFound)
                    {
                        std::string method, basis;
                        splitLot(lot, &method, &basis);
                        double      val, T = -1;
                        if (mpi.getProp(mpo, iqmType::Exp, method, basis, "", &val, nullptr, &T))
                        {
                            nexp++;
                        }
                        else if (mpi.getProp(mpo, iqmType::QM, method, basis, "", (char *)"electronic", &val, nullptr, &T))
                        {
                            nqm++;
                        }
                    }
                }
            }

            /* Determine Miller and Bosque polarizabilities for this Alexandria element */
            size_t      pos   = pType->getType().find("p_");
            std::string ptype = pType->getType();
            if (pos != std::string::npos)
            {
                ptype = ptype.substr(pos+2);
            }
            snprintf(longbuf, STRLEN, "%s & %s & %s (%s) & %s & %s & %s",
                     ptype.c_str(),
                     (nexp > 0)     ? gmx_itoa(nexp).c_str()     : "-",
                     (pType->getPolarizability() > 0)  ? gmx_ftoa(pType->getPolarizability()).c_str()  : "-",
                     (pType->getSigPol() > 0)  ? gmx_ftoa(pType->getSigPol()).c_str() : "-",
                     (bos_pol > 0)     ? gmx_ftoa(bos_pol).c_str()     : "-");
            lt.printLine(longbuf);
        }
    }
    lt.printFooter();
    fflush(fp);
}

static void alexandria_molprop_atomtype_dip_table(FILE                       *fp,
                                                  const std::vector<Poldata> &pd)
{
    int         cur        = 0;
    std::string gt_type[2] = { "", "" };

#define prev (1-cur)
    std::string longbuf;
    LongTable   lt(fp, true, nullptr);

    lt.setCaption("Electronegativity equalization parameters for Alexandria models. $J_0$ and $\\chi_0$ in eV, $\\zeta$ in 1/nm.");
    lt.setLabel("eemparams");
    lt.setColumns(5);

    longbuf = gmx::formatString("Model & Type & $J_0$ & $\\chi_0$ & $\\zeta$");
    lt.addHeadLine(longbuf);
    lt.printHeader();

    for (auto ipd : pd)
    {
        for (auto aType = ipd.getAtypeBegin(); aType != ipd.getAtypeEnd(); aType++)
        {
            gt_type[cur] = aType->getType();
            if (((0 == gt_type[prev].size()) || (gt_type[cur] != gt_type[prev])))
            {
                longbuf = gmx::formatString("%s & %s",
                                            chargeTypeName(ipd.chargeType()).c_str(),
                                            gt_type[cur].c_str());
                if (ipd.haveEemSupport(gt_type[cur], false))
                {
                    longbuf.append(gmx::formatString(" & %.3f", ipd.getJ00(gt_type[cur])));
                    longbuf.append(gmx::formatString(" & %.3f", ipd.getChi0(gt_type[cur])));
                    int nzeta = ipd.getNzeta(gt_type[cur]);
                    int i     = 0;
                    for (; i < nzeta; i++)
                    {
                        longbuf.append(gmx::formatString(" & %.3f", ipd.getZeta(gt_type[cur], i+1)));
                    }
                    for (; i < 2; i++)
                    {
                        longbuf.append(" &");
                    }
                }
                else
                {
                    longbuf.append(" & & & &");
                }
                lt.printLine(longbuf);
            }
        }
        cur = prev;
    }
    lt.printFooter();
}

void alexandria_molprop_atomtype_table(FILE                       *fp,
                                       bool                        bPolar,
                                       const std::vector<Poldata> &pd,
                                       const std::vector<MolProp> &mp,
                                       const char                 *lot)
{
    if (bPolar)
    {
        alexandria_molprop_atomtype_polar_table(fp, &pd[0], mp, lot);
    }
    else
    {
        alexandria_molprop_atomtype_dip_table(fp, pd);
    }
}
#endif

static void prop_header(LongTable     &lt,
                        const char    *property,
                        const char    *unit,
                        real           rel_toler,
                        real           abs_toler,
                        const QmCount  qmc,
                        iMolSelect     ims,
                        bool           bPrintConf,
                        bool           bPrintBasis,
                        bool           bPrintMultQ)
{
    int  i, k, nc;
    char columns[256];
    char longbuf[STRLEN];
    char buf[256];

    snprintf(columns, 256, "p{75mm}");
    nc = 2 + qmc.nCalc();
    if (bPrintMultQ)
    {
        nc += 2;
    }
    if (bPrintConf)
    {
        nc++;
    }
    for (i = 0; (i < nc); i++)
    {
        strncat(columns, "c", 256-strlen(columns)-1);
    }
    lt.setColumns(columns);

    for (k = 0; (k < 2); k++)
    {
        if (0 == k)
        {
            snprintf(longbuf, STRLEN, "Comparison of experimental %s to calculated values. {\\bf Data set: %s}. Calculated numbers that are more than %.0f%s off the experimental values are printed in bold, more than %.0f%s off in bold red.",
                     property, iMolSelectName(ims),
                     (abs_toler > 0) ? abs_toler   : 100*rel_toler,
                     (abs_toler > 0) ? unit : "\\%",
                     (abs_toler > 0) ? 2*abs_toler : 200*rel_toler,
                     (abs_toler > 0) ? unit : "\\%");
            lt.setCaption(longbuf);
            snprintf(longbuf, STRLEN, "%s", iMolSelectName(ims));
            lt.setLabel(longbuf);
        }
        else
        {
            snprintf(longbuf, STRLEN, "Molecule & Form. %s %s & Exper. ",
                     bPrintMultQ ? "& q & mult" : "",
                     bPrintConf  ? "& Conf." : "");
            for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q)
            {
                snprintf(buf, 256, "& %s", q->method().c_str());
                strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
            }
            lt.addHeadLine(longbuf);

            if (bPrintBasis)
            {
                snprintf(longbuf, STRLEN, " & & %s %s",
                         bPrintMultQ ? "& &" : "",
                         bPrintConf ? "&" : "");
                for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q)
                {
                    snprintf(buf, 256, "& %s", q->basis().c_str());
                    strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
                }
                lt.addHeadLine(longbuf);
            }
            snprintf(longbuf, STRLEN, "Type & &%s %s",
                     bPrintMultQ ? "& &" : "",
                     bPrintConf ? "&" : "");
            for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q)
            {
                snprintf(buf, 256, "& %s", q->type().c_str());
                strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
            }
            lt.addHeadLine(longbuf);
        }
    }
    lt.printHeader();
}

static int outside(real vexp, real vcalc, real rel_toler, real abs_toler)
{
    real rdv, adv = fabs(vexp-vcalc);
    if (abs_toler > 0)
    {
        if (adv > 2*abs_toler)
        {
            return 2;
        }
        else if (adv > abs_toler)
        {
            return 1;
        }
        return 0;
    }
    else
    {
        if (vexp == 0)
        {
            return 0;
        }
        rdv = adv/vexp;

        if (rdv > 2*rel_toler)
        {
            return 2;
        }
        else if (rdv > rel_toler)
        {
            return 1;
        }
        return 0;
    }
}

void alexandria_molprop_prop_table(FILE                 *fp,
                                   MolPropObservable     mpo,
                                   real                  rel_toler,
                                   real                  abs_toler,
                                   std::vector<MolProp> &mp,
                                   const QmCount        &qmc,
                                   bool                  bPrintAll,
                                   bool                  bPrintBasis,
                                   bool                  bPrintMultQ,
                                   const MolSelect      &gms,
                                   iMolSelect            ims)
{
#define BLEN 1024
    std::vector<double>         vec;

    LongTable                   lt(fp, true, "small");

    int nprint = std::count_if(mp.begin(), mp.end(),
                               [ims, gms](const MolProp &mpi)
                               { iMolSelect myIms;
                                   return (gms.status(mpi.getIupac(), &myIms) &&
                                           ims == myIms &&
                                           mpi.hasAllAtomTypes()); });
    if (nprint <= 0)
    {
        return;
    }
    bool bPrintConf = false;
    prop_header(lt, mpo_name(mpo), mpo_unit(mpo),
                rel_toler, abs_toler, qmc,
                ims, bPrintConf, bPrintBasis, bPrintMultQ);
    for (auto &mpi : mp)
    {
        iMolSelect myIms;
        if (gms.status(mpi.getIupac(), &myIms) &&
            ims == myIms &&
            mpi.hasAllAtomTypes())
        {
            std::vector<ExpData>  ed;
            std::vector<CalcData> cd;
            for (auto ei : mpi.experimentConst())
            {
                for(const auto &prop : ei.propertyConst())
                {
                    for(const auto &pp : prop.second)
                    {
                        if (ei.dataSource() == dsExperiment)
                        {
                            ed.push_back(ExpData(pp->getValue(),
                                                 pp->getError(),
                                                 pp->getTemperature(),
                                                 ei.getReference(),
                                                 ei.getConformation()));
                        }
                        else
                        {
                            cd.push_back(CalcData(pp->getValue(),
                                                  pp->getError(),
                                                  pp->getTemperature(),
                                                  1));
                        }
                    }
                }
            }
#ifdef OLD
            int nqm = 0;
            for (int nexp = 0; (nexp < (int)ed.size()); nexp++)
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Found %d experiments and %d calculations\n",
                            (int)ed.size(), nqm);
                }
                if ((bPrintAll || (ed.size() > 0))  && (nqm > 0))
                {
                    if (0 == nexp)
                    {
                        iprint++;
                    }
                    std:: string myline;
                    if (nexp == 0)
                    {
                        if (bPrintMultQ)
                        {
                            snprintf(mylbuf, sizeof(mylbuf), "%d. %-15s & %s & %d & %d",
                                     iprint, mpi.getIupac().c_str(),
                                     mpi.getTexFormula().c_str(),
                                     mpi.totalCharge(),
                                     mpi.getMultiplicity());
                        }
                        else
                        {
                            snprintf(mylbuf, sizeof(mylbuf), "%d. %-15s & %s",
                                     iprint, mpi.getIupac().c_str(),
                                     mpi.getTexFormula().c_str());
                        }
                    }
                    else
                    {
                        snprintf(mylbuf, sizeof(mylbuf), " & ");
                    }
                    myline.append(mylbuf);
                    if (bPrintConf)
                    {
                        snprintf(mylbuf, BLEN, "      & %s ", ((ed[nexp].conf_.size() > 0) ?
                                                               ed[nexp].conf_.c_str() : "-"));
                        myline.append(mylbuf);
                    }
                    if (ed.size() > 0)
                    {
                        snprintf(mylbuf, sizeof(mylbuf), "& %8.3f", ed[nexp].val_);
                        myline.append(mylbuf);
                        if (ed[nexp].err_ > 0)
                        {
                            snprintf(mylbuf, sizeof(mylbuf), "(%.3f)", ed[nexp].err_);
                            myline.append(mylbuf);
                        }
                        if (strcmp(ed[nexp].ref_.c_str(), "Maaren2017a") == 0)
                        {
                            snprintf(mylbuf, sizeof(mylbuf), " (*)");
                        }
                        else
                        {
                            snprintf(mylbuf, sizeof(mylbuf), "~\\cite{%s} ", ed[nexp].ref_.c_str());
                        }
                        myline.append(mylbuf);
                    }
                    else
                    {
                        snprintf(mylbuf, sizeof(mylbuf), "& - ");
                        myline.append(mylbuf);
                    }
                    for (size_t j = 0; (j < qmc.nCalc()); j++)
                    {
                        if (cd[j].found_ > 0)
                        {
                            vc = cd[j].val_;
                            if (cd[j].err_ > 0)
                            {
                                snprintf(vbuf, sizeof(vbuf), "%8.2f(%.2f)", vc, cd[j].err_);
                            }
                            else
                            {
                                snprintf(vbuf, sizeof(vbuf), "%8.2f", vc);
                            }
                            if (ed.size() > 0)
                            {
                                int oo = outside(ed[nexp].val_, vc, rel_toler, abs_toler);
                                switch (oo)
                                {
                                    case 2:
                                        snprintf(mylbuf, sizeof(mylbuf), "& \\textcolor{Red}{\\bf %s} ", vbuf);
                                        break;
                                    case 1:
                                        snprintf(mylbuf, sizeof(mylbuf), "& {\\bf %s} ", vbuf);
                                        break;
                                    default:
                                        snprintf(mylbuf, sizeof(mylbuf), "& %s ", vbuf);
                                }
                            }
                            else
                            {
                                snprintf(mylbuf, sizeof(mylbuf), "& %s ", vbuf);
                            }
                            myline.append(mylbuf);
                        }
                        else
                        {
                            snprintf(mylbuf, sizeof(mylbuf), "& ");
                            myline.append(mylbuf);
                        }
                    }
                    lt.printLine(myline.c_str());
                }
            }
#endif
        }
    }
    lt.printFooter();
}

} // namespace
