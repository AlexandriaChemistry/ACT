/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include "tuning_utility.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "act/alexandria/actmol.h"
#include "act/alexandria/pdbwriter.h"
#include "act/alexandria/thermochemistry.h"
#include "act/basics/interactiontype.h"
#include "act/forces/forcecomputer.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/multipole_names.h"
#include "act/qgen/qtype.h"
#include "act/utility/stringutil.h"
#include "act/utility/units.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/gmxassert.h"

void doAddOptions(std::vector<t_pargs> *pargs, size_t npa, t_pargs pa[])
{
    for (size_t i = 0; i < npa; i++)
    {
        pargs->push_back(pa[i]);
    }
}

namespace alexandria
{

static void print_stats(FILE        *fp,
                        const char  *prop,
                        const char  *unit,
                        double       conversion,
                        gmx_stats   *lsq,
                        bool         bHeader,
                        const char  *xaxis,
                        const char  *yaxis,
                        bool         useOffset)
{
    real a    = 0, da  = 0, b    = 0, db   = 0;
    real mse  = 0, mae = 0, chi2 = 0, rmsd = 0;
    real Rfit = 0;
    int  n    = lsq->get_npoints();
    // TODO add error checking for lsq
    if (n == 0)
    {
        return;
    }
    std::string newprop = gmx::formatString("%s (%s)", prop, unit);
    if (useOffset)
    {
        if (bHeader)
        {
            fprintf(fp, "Fitting data to y = ax + b, where x = %s\n", xaxis);
            fprintf(fp, "%-32s %6s %13s %13s %7s %8s %8s %8s %10s\n",
                    "Property", "N", "a", "b", "R(%)", "RMSD", "MSE", "MAE", "Model");
            fprintf(fp, "------------------------------------------------------------------------------------------------\n");
        }
        eStats ok = lsq->get_ab(elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &Rfit);
        if (eStats::OK == ok)
        {
            ok = lsq->get_rmsd(&rmsd);
        }
        if (eStats::OK == ok)
        {
            ok = lsq->get_mse_mae(&mse, &mae);
        }
        if (eStats::OK == ok)
        {
            fprintf(fp, "%-32s %6d %6.3f(%5.3f) %6.3f(%5.3f) %7.2f %8.4f %8.4f %8.4f %10s\n",
                    newprop.c_str(), n, a, da, 
                    b, db, Rfit*100, 
                    conversion*rmsd, conversion*mse, 
                    conversion*mae, yaxis);
        }
        else
        {
            fprintf(fp, "Statistics problem for %s: %s\n", prop, gmx_stats_message(ok));
        }
    }
    else
    {
        if (bHeader)
        {
            fprintf(fp, "Fitting data to y = ax, where x = %s\n", xaxis);
            fprintf(fp, "%-32s %6s %13s %7s %8s %8s %8s %10s\n",
                    "Property", "N", "a", "R(%)", "RMSD", "MSE", "MAE", "Model");
            fprintf(fp, "----------------------------------------------------------------------------------------------\n");
        }
        eStats ok = lsq->get_a(elsqWEIGHT_NONE, &a, &da, &chi2, &Rfit);
        if (eStats::OK == ok)
        {
            ok = lsq->get_rmsd(&rmsd);
        }
        if (eStats::OK == ok)
        {
            ok = lsq->get_mse_mae(&mse, &mae);
        }
        if (eStats::OK == ok)
        {
            fprintf(fp, "%-32s %6d %6.3f(%5.3f) %7.2f %8.4f %8.4f %8.4f %10s\n",
                    newprop.c_str(), n, a, da, Rfit*100,
                    conversion*rmsd, conversion*mse, conversion*mae, yaxis);
        }
        else
        {
            fprintf(fp, "Statistics problem for %s: %s\n", prop, gmx_stats_message(ok));
        }
    }
}

static void print_lsq_set(FILE *fp, const gmx_stats &lsq)
{
    fprintf(fp, "@type nxy\n");
    auto x = lsq.getX();
    auto y = lsq.getY();
    for(size_t i = 0; i < x.size(); i++)
    {
        fprintf(fp, "%10g  %10g\n", x[i], y[i]);
    }
    fprintf(fp, "&\n");
}

static void xvgr_symbolize(FILE                           *xvgf,
                           const std::vector<std::string>  leg,
                           const gmx_output_env_t         *oenv)
{
    xvgrLegend(xvgf, leg, oenv);
    for (size_t i = 0; (i < leg.size()); i++)
    {
        xvgr_line_props(xvgf, i, elNone, ecBlack+i, oenv);
        fprintf(xvgf, "@ s%d symbol %d\n", static_cast<int>(i), static_cast<int>(i+1));
    }
}

static void print_one_alpha(FILE *fp, const std::string &label,
                            const tensor alpha)
{
    double fac = convertFromGromacs(1.0, "Angstrom3");
    fprintf(fp,
            "%s   (%6.2f %6.2f %6.2f)\n"
            "             (%6s %6.2f %6.2f)\n"
            "             (%6s %6s %6.2f)\n",
            label.c_str(),
            fac*alpha[XX][XX], fac*alpha[XX][YY], fac*alpha[XX][ZZ],
            "", fac*alpha[YY][YY], fac*alpha[YY][ZZ],
            "", "", fac*alpha[ZZ][ZZ]);
 }

static void print_polarizability(FILE              *fp,
                                 const ACTMol       *mol,
                                 const std::string &calc_name,
                                 real               alpha_toler,
                                 real               isopol_toler)
{
    tensor dalpha;
    real   delta      = 0;
    real   diso_pol   = 0;
    real   daniso_pol = 0;

    double fac = convertFromGromacs(1.0, "Angstrom3");
    if (!calc_name.empty())
    {
        auto qprops = mol->qPropsConst();
        for(auto qp = qprops.begin(); qp < qprops.end(); ++qp)
        {
            auto qcalc = qp->qPactConst();
            auto acalc = qcalc.polarizabilityTensor();
            auto qelec = qp->qPqmConst();
            if (!qelec.hasPolarizability())
            {
                print_one_alpha(fp, calc_name, acalc);
                continue;
            }
            tensor aelec;
            copy_mat(qelec.polarizabilityTensor(), aelec);
            print_one_alpha(fp, "Electronic", aelec);
        
            m_sub(aelec, acalc, dalpha);
            delta = fac*sqrt(gmx::square(dalpha[XX][XX])+gmx::square(dalpha[XX][YY])+
                             gmx::square(dalpha[XX][ZZ])+
                             gmx::square(dalpha[YY][YY])+gmx::square(dalpha[YY][ZZ]));
            diso_pol   = fac*std::abs(qcalc.isotropicPolarizability()-
                                      qelec.isotropicPolarizability());
            daniso_pol = fac*std::abs(qcalc.anisotropicPolarizability()-
                                      qelec.anisotropicPolarizability());
            fprintf(fp,
                    "%-12s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "             (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "             (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name.c_str(),
                    fac*acalc[XX][XX], fac*acalc[XX][YY], fac*acalc[XX][ZZ],
                    fac*dalpha[XX][XX], fac*dalpha[XX][YY], fac*dalpha[XX][ZZ], delta, (delta > alpha_toler) ? "ALPHA" : "",
                    "", fac*acalc[YY][YY], fac*acalc[YY][ZZ],
                    "", fac*dalpha[YY][YY], fac*dalpha[YY][ZZ],
                    "", "", fac*acalc[ZZ][ZZ],
                    "", "", fac*dalpha[ZZ][ZZ]);
            fprintf(fp,
                    "Isotropic polarizability:  %s Electronic: %6.2f  Calculated: %6.2f  Delta: %6.2f %s\n\n",
                    mol->getMolname().c_str(), 
                    fac*qelec.isotropicPolarizability(),
                    fac*qcalc.isotropicPolarizability(), 
                    diso_pol, (diso_pol > isopol_toler) ? "ISO" : "");
            fprintf(fp,
                    "Anisotropic polarizability:  %s Electronic: %6.2f  Calculated: %6.2f  Delta: %6.2f %s\n\n",
                    mol->getMolname().c_str(), 
                    fac*qelec.anisotropicPolarizability(),
                    fac*qcalc.anisotropicPolarizability(), 
                    daniso_pol, (daniso_pol > isopol_toler) ? "ANISO" : "");
        }
    }
}

static void analyse_multipoles(FILE                                                        *fp,
                               const std::vector<alexandria::ACTMol>::iterator             &mol,
                               std::map<MolPropObservable, double>                          toler,
                               std::map<MolPropObservable, std::map<iMolSelect, qtStats> > *lsq)
{
    auto qprops = mol->qPropsConst();
    for(auto qp = qprops.begin(); qp < qprops.end(); ++qp)
    {
        auto qelec = qp->qPqmConst();
        auto qcalc = qp->qPact();
        qcalc->initializeMoments();
        qcalc->calcMoments();

        for(auto &mpo : mpoMultiPoles)
        {
            const char *name   = mpo_name(mpo);
            const char *unit   = mpo_unit2(mpo);
            double      factor = convertFromGromacs(1, unit);
            std::vector<double> Telec;
            if (qelec.hasMultipole(mpo))
            {
                fprintf(fp, "Electronic %s (%s):\n", name, unit);
                Telec = qelec.getMultipole(mpo);
                printMultipole(fp, mpo, Telec);
            }
            real delta = 0;
            auto Tcalc = qcalc->getMultipole(mpo);
            fprintf(fp, "Calc %s (%s):\n", name, unit);
            printMultipole(fp, mpo, Tcalc);

            std::vector<double> diff;
            auto qt = qcalc->qtype();
            if (Telec.size() == Tcalc.size())
            {
                for(size_t i = 0; i < Tcalc.size(); i++)
                {
                    auto tc = Tcalc[i];
                    auto te = Telec[i];
                    diff.push_back(tc-te);
                        delta += gmx::square(tc-te);
                        (*lsq)[mpo][mol->datasetType()][qt].add_point(factor*te, factor*tc, 0, 0);
                }
                double rms = std::sqrt(delta/Tcalc.size());
                std::string flag("");
                if (rms > toler[mpo])
                {
                    flag = " MULTI";
                }
                fprintf(fp, "%s-Electronic Norm %g RMS = %g (%s)%s:\n",
                        qTypeName(qt).c_str(), factor*std::sqrt(delta), factor*rms, unit, flag.c_str());
                printMultipole(fp, mpo, diff);
            }
        }
    }
}

static void print_corr(const char                    *outfile,
                       const char                    *title,
                       const char                    *xaxis,
                       const char                    *yaxis, 
                       std::map<iMolSelect, qtStats> &stats,
                       const gmx_output_env_t        *oenv)
{
    std::vector<std::string> eprnm;
    std::vector<gmx_stats>   lsq;
    bool                     haveData = false;
    for (auto &ims : iMolSelectNames())
    {
        auto qs = stats.find(ims.first);
        if (qs != stats.end())
        {
            for (auto &i : qs->second)
            {
                int N = i.second.get_npoints();
                if (N > 0)
                {
                    eprnm.push_back(gmx::formatString("%s-%s", qTypeName(i.first).c_str(), ims.second));
                    lsq.push_back(i.second);
                    haveData = true;
                }
            }
        }
    }
    if (haveData && outfile)
    {
        FILE *muc = xvgropen(outfile, title, xaxis, yaxis, oenv);
        xvgr_symbolize(muc, eprnm, oenv);
        for(size_t i = 0; i < lsq.size(); i++)
        {
            print_lsq_set(muc, lsq[i]);
        }
        fclose(muc);
    }
}

static void print_corr(const char                      *outfile,
                       const char                      *title,
                       const char                      *xaxis,
                       const char                      *yaxis, 
                       std::map<iMolSelect, gmx_stats> &stats,
                       const gmx_output_env_t          *oenv)
{
    std::vector<std::string> eprnm;
    std::vector<gmx_stats>   lsq;
    bool                     haveData = false;
    for (auto &ims : iMolSelectNames())
    {
        auto qs = stats.find(ims.first);
        if (qs != stats.end())
        {
            int N = qs->second.get_npoints();
            if (N > 0)
            {
                eprnm.push_back(gmx::formatString("%s", ims.second));
                lsq.push_back(qs->second);
                haveData = true;
            }
        }
    }
    if (haveData && outfile)
    {
        FILE *muc = xvgropen(outfile, title, xaxis, yaxis, oenv);
        xvgr_symbolize(muc, eprnm, oenv);
        for(size_t i = 0; i < lsq.size(); i++)
        {
            print_lsq_set(muc, lsq[i]);
        }
        fclose(muc);
    }
}

static void write_q_histo(FILE                             *fplog,
                          const char                       *qhisto,
                          std::map<std::string, gmx_stats> *lsqt,
                          const gmx_output_env_t           *oenv,
                          qtStats                          *lsq_charge,
                          bool                              useOffset)
{
    std::vector<std::string> types;
    for (const auto &i: *lsqt)
    {
        int N = i.second.get_npoints();
        if (N > 0)
        {
            types.push_back(gmx::formatString("%s (N=%d)", i.first.c_str(), N));
        }
    }
    FILE *hh = nullptr;
    if (nullptr != qhisto)
    {
        hh = xvgropen(qhisto, "Histogram for charges", "q (e)", "a.u.", oenv);
        xvgrLegend(hh, types, oenv);
    }
    auto model = qTypeName(qType::Calc);
    auto gs    = lsq_charge->find(qType::Calc);
    if (gs != lsq_charge->end())
    {
        print_stats(fplog, "All Partial Charges", "e", 1.0, &gs->second, true,
                    qTypeName(qType::CM5).c_str(), model.c_str(), useOffset);
    }
    for (auto &q : *lsqt)
    {
        if (q.second.get_npoints() > 0)
        {
            int  nbins = 20;
            std::vector<double> x, y;
            if (q.second.make_histogram(0, &nbins, eHisto::Y, 1, &x, &y) == eStats::OK)
            {
                if (hh)
                {
                    fprintf(hh, "@type xy\n");
                    for (int i = 0; i < nbins; i++)
                    {
                        fprintf(hh, "%10g  %10g\n", x[i], y[i]);
                    }
                    fprintf(hh, "&\n");
                }
                print_stats(fplog, q.first.c_str(), "e", 1.0, &q.second, false, 
                            "CM5", model.c_str(), useOffset);
            }
        }
    }
    if (hh)
    {
        fclose(hh);
    }
    fprintf(fplog, "\n");
}

void TuneForceFieldPrinter::addOptions(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] =
    {
        { "-esp_toler", FALSE, etREAL, {&esp_toler_},
          "Tolerance (kJ/mol e) for marking ESP as an outlier in the log file" },
        { "-dip_toler", FALSE, etREAL, {&dip_toler_},
          "Tolerance (Debye) for marking dipole as an outlier in the log file" },
        { "-quad_toler", FALSE, etREAL, {&quad_toler_},
          "Tolerance (Buckingham) for marking quadrupole as an outlier in the log file" },
        { "-oct_toler", FALSE, etREAL, {&oct_toler_},
          "Tolerance (Buckingham) for marking octupole as an outlier in the log file" },
        { "-hex_toler", FALSE, etREAL, {&hex_toler_},
          "Tolerance (Buckingham) for marking hexapole as an outlier in the log file" },
        { "-alpha_toler", FALSE, etREAL, {&alpha_toler_},
          "Tolerance (A^3) for marking diagonal elements of the polarizability tensor as an outlier in the log file" },
        { "-isopol_toler", FALSE, etREAL, {&isopol_toler_},
          "Tolerance (A^3) for marking isotropic polarizability as an outlier in the log file" },
        { "-use_offset", FALSE, etBOOL,{&useOffset_},
          "Fit regression analysis of results to y = ax+b instead of y = ax" },
        { "-printSP", FALSE, etBOOL,{&printSP_},
          "Print information from all single point calculations in detail" },
        { "-calc_frequencies",  FALSE, etBOOL, {&calcFrequencies_},
          "Perform energy minimization and compute vibrational frequencies for each molecule (after optimizing the force field if -optimize is enabled). If turned on, this option will also generate thermochemistry values based on the force field." },
        { "-diatomics", FALSE, etBOOL, {&diatomic_},
          "Analyse diatomic compounds by plotting the reference energy minus the ACT Coulomb terms as a function of distance in one xvg file per compound." },
        { "-dump_outliers", FALSE, etREAL, {&dumpOutliers_},
          "Dump XYZ files for outliers with deviation larger than this number times the standard deviation of the energy (DeltaE0) or interaction energy. If this number is less than zero, nothing will be dumped." }
    };
    doAddOptions(pargs, sizeof(pa)/sizeof(pa[0]), pa);
}

void TuneForceFieldPrinter::addFileOptions(std::vector<t_filenm> *filenm)
{
    std::vector<t_filenm> fnm = {
        { efXVG, "-qhisto",    "q_histo",       ffOPTWR },
        { efXVG, "-epotcorr",  "epot_corr",     ffOPTWR },
        { efXVG, "-eintercorr","einter_corr",   ffOPTWR },
        { efXVG, "-forcecorr", "force_corr",    ffOPTWR },
        { efXVG, "-freqcorr",  "freq_corr",     ffOPTWR },
        { efXVG, "-espcorr",   "esp_corr",      ffOPTWR },
        { efXVG, "-alphacorr", "alpha_corr",    ffOPTWR },
        { efXVG, "-qcorr",     "q_corr",        ffOPTWR },
        { efXVG, "-isopol",    "isopol_corr",   ffOPTWR },
        { efXVG, "-anisopol",  "anisopol_corr", ffOPTWR }
    };
    for(size_t i = 0; i < fnm.size(); i++)
    {
        filenm->push_back(fnm[i]);
    }
    std::vector<t_filenm> multi;
    multi.resize(mpoMultiPoles.size());
    size_t ii = 0;
    for(auto &mpo : mpoMultiPoles)
    {
        std::string cmdFlag = gmx::formatString("-%scorr", mpo_name(mpo));
        std::string defFnm  = gmx::formatString("%s_corr", mpo_name(mpo));
        multi[ii].ftp  = efXVG;
        multi[ii].opt  = strdup(cmdFlag.c_str());
        multi[ii].fn   = strdup(defFnm.c_str());
        multi[ii].flag = ffOPTWR;
        filenm->push_back(multi[ii]);
    }
}

void TuneForceFieldPrinter::analysePolarisability(FILE                *fp,
                                                  const ForceField    *pd,
                                                  alexandria::ACTMol  *mol,
                                                  const ForceComputer *forceComp,
                                                  qtStats             *lsq_isoPol,
                                                  qtStats             *lsq_anisoPol,
                                                  qtStats             *lsq_alpha)
{
    for(auto qp = mol->qProps()->begin(); qp < mol->qProps()->end(); ++qp)
    {
        auto qcalc = qp->qPact();
        if (!qcalc->hasPolarizability())
        {
            continue;
        }
        qcalc->calcPolarizability(pd, mol->topology(), forceComp);
        auto acalc = qcalc->polarizabilityTensor();
    
        fprintf(fp, "Polarizability:\n");
        print_polarizability(fp, mol, qTypeName(qType::Calc), alpha_toler_, isopol_toler_);
        auto qelec = qp->qPqmConst();
        if (qelec.hasPolarizability())
        {
            auto aelec = qelec.polarizabilityTensor();
            (*lsq_isoPol)[qType::Calc].add_point(qelec.isotropicPolarizability(),
                                                 qcalc->isotropicPolarizability(),
                                                 0, 0);
            (*lsq_anisoPol)[qType::Calc].add_point(qelec.anisotropicPolarizability(),
                                                   qcalc->anisotropicPolarizability(),
                                                   0, 0);
            for (int mm = 0; mm < DIM; mm++)
            {
                (*lsq_alpha)[qType::Calc].add_point(aelec[mm][mm], acalc[mm][mm], 0, 0);
            }
        }
    }
}

void TuneForceFieldPrinter::printAtoms(FILE                         *fp,
                                       alexandria::ACTMol           *mol,
                                       const std::vector<gmx::RVec> &coords,
                                       const std::vector<gmx::RVec> &forces)
{
    std::map<qType, const std::vector<double> > qQM;
    std::vector<qType>                    typeQM = { 
        qType::CM5, qType::ESP, qType::Hirshfeld, qType::Mulliken
    };
    for(auto &qt : typeQM)
    {
        for(auto qp = mol->qProps()->begin(); qp < mol->qProps()->end(); ++qp)
        {
            auto qelec = qp->qPqmConst();
            if (qelec.qtype() == qt)
            {
                std::vector qq = qelec.charge();
                qQM.insert({qt, std::move(qq)});
            }
        }
    }
    fprintf(fp, "Atom   Type            ACM");
    for(auto &qt : qQM)
    {
        if (!qQM.find(qt.first)->second.empty())
        {
            fprintf(fp, "%10s", qTypeName(qt.first).c_str());
        }
    }
    fprintf(fp, "        x        y   z(pm)          fx         fy fz(kJ/mol nm)     qtot\n");
    int      i       = 0;
    double   qtot    = 0;
    auto    &myatoms = mol->atomsConst();
    for (size_t j = i = 0; j < myatoms.size(); j++)
    {
        if (myatoms[j].pType() == eptAtom ||
            myatoms[j].pType() == eptNucleus)
        {
            real qCalc = myatoms[j].charge();
            fprintf(fp, "%-2d%3lu  %-5s  %12g",
                    myatoms[j].atomicNumber(),
                    j+1,
                    myatoms[j].ffType().c_str(),
                    qCalc);
            qtot += qCalc;
            for(auto &qt : qQM)
            {
                if (!qQM.find(qt.first)->second.empty())
                {
                    fprintf(fp, "  %8.4f", qt.second[i]);
                }
            }
            fprintf(fp," %8.3f %8.3f %8.3f %10.3f %10.3f %10.3f  %10g\n", 
                    convertFromGromacs(coords[j][XX], "pm"),
                    convertFromGromacs(coords[j][YY], "pm"),
                    convertFromGromacs(coords[j][ZZ], "pm"),
                    forces[j][XX], forces[j][YY], forces[j][ZZ],
                    qtot);
            i++;
        }
        else
        {
            // Turned on printing of shells again
            fprintf(fp, "%-2d%3lu  %-5s  %12g",
                    0,
                    j+1,
                    myatoms[j].ffType().c_str(),
                    myatoms[j].charge());
            qtot += myatoms[j].charge();
            for(auto &qt : qQM)
            {
                if (!qQM.find(qt.first)->second.empty())
                {
                    fprintf(fp, "          ");
                }
            }
            fprintf(fp," %8.3f %8.3f %8.3f %10.3f %10.3f %10.3f  %10g\n", 
                    convertFromGromacs(coords[j][XX], "pm"),
                    convertFromGromacs(coords[j][YY], "pm"),
                    convertFromGromacs(coords[j][ZZ], "pm"),
                    forces[j][XX], forces[j][YY], forces[j][ZZ],
                    qtot);
        }
    }
}

/*! Store coordinates to a file, first original, then minimized.
 * \param[in] pdb   File name
 * \param[in] title Text string
 */
static void writeCoordinates(const std::vector<ActAtom>  &atoms,
                             const std::string           &pdb,
                             const std::string            &title,
                             const std::vector<gmx::RVec> &xorig,
                             const std::vector<gmx::RVec> &xmin)
{
    FILE *out = gmx_fio_fopen(pdb.c_str(), "w");
    matrix box = { { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 } };
    std::vector<std::string> resnames = { "A", "B" };
    pdbWriter(out, title.c_str(), atoms, xorig, resnames,
              epbcNONE, box, 'A', 1, {}, nullptr, false);

    pdbWriter(out, title.c_str(), atoms, xmin, resnames,
              epbcNONE, box, 'B', 1, {}, nullptr, false);

    gmx_fio_fclose(out);
}

static void normalize_vector(std::vector<double> *v)
{
    double total = 0;
    for(size_t i = 0; i < v->size(); i++)
    {
        total += (*v)[i];
    }
    if (total == 0)
    {
        return;
    }
    double tinv = 1.0/total;
    for(size_t i = 0; i < v->size(); i++)
    {
        (*v)[i] *= tinv;
    }
}

static void plot_spectrum(const char                *filenm,
                          const std::string         &molname,
                          double                     lineWidth,
                          const std::vector<double> &act_freq,
                          const std::vector<double> &act_intens,
                          const std::vector<double> &qm_freq,
                          const std::vector<double> &qm_intens,
                          gmx_output_env_t          *oenv,
                          bool                       normalize)
{
    if (nullptr == filenm)
    {
        return;
    }
    auto minfreq = *std::min_element(act_freq.begin(), act_freq.end());
    auto maxfreq = *std::max_element(act_freq.begin(), act_freq.end());
    if (!qm_freq.empty())
    {
        minfreq = std::min(minfreq, *std::min_element(qm_freq.begin(), qm_freq.end()));
        maxfreq = std::max(maxfreq, *std::max_element(qm_freq.begin(), qm_freq.end()));
    }
    int minf = convertFromGromacs(minfreq, "cm^-1")-5*lineWidth;
    int maxf = convertFromGromacs(maxfreq, "cm^-1")+5*lineWidth;
    std::vector<double> act_spec(maxf-minf+1, 0.0);
    std::vector<double> qm_spec;
    if (!qm_freq.empty())
    {
        qm_spec.resize(maxf-minf+1, 0.0);
    }
    for(size_t peak = 0; peak < act_freq.size(); peak++)
    {
        double act_ff, qm_ff = 0;
        act_ff = convertFromGromacs(act_freq[peak], "cm^-1");
        if (!qm_freq.empty())
        {
            qm_ff = convertFromGromacs(qm_freq[peak], "cm^-1");
        }
        for(int f = minf; f <= maxf; f++)
        {
            double act_ds = lineWidth/(2*M_PI*(gmx::square(f-act_ff)+gmx::square(lineWidth/2)));
            act_spec[f-minf] += act_ds*act_intens[peak];
            if (!qm_freq.empty())
            {
                double qm_ds = lineWidth/(2*M_PI*(gmx::square(f-qm_ff)+gmx::square(lineWidth/2)));
                qm_spec[f-minf] += qm_ds*qm_intens[peak];
            }
        }
    }
    if (normalize)
    {
        normalize_vector(&act_spec);
        normalize_vector(&qm_spec);
    }
    std::string title = gmx::formatString("Infrared spectrum for %s", molname.c_str());
    FILE *xvg = xvgropen(filenm, title.c_str(), "Frequency (1/cm)",
                         "Intensity (a.u.)", oenv);
    if (!qm_freq.empty())
    {
        std::vector<const char *> legend = { "QM", "ACT" };
        xvgr_legend(xvg, legend.size(), legend.data(), oenv);
    }
    for(int f = minf; f <= maxf; f++)
    {
        if (qm_freq.empty())
        {
            fprintf(xvg, "%6d  %10g\n", f, act_spec[f-minf]);
        }
        else
        {
            fprintf(xvg, "%6d  %10g  %10g\n", f, qm_spec[f-minf], act_spec[f-minf]);
        }
    }
    fprintf(xvg, "&\n");
    xvgrclose(xvg);
}

static void addThermo(JsonTree                  *jtree, 
                      alexandria::ACTMol         *mol,
                      std::vector<gmx::RVec>    *coords,
                      const AtomizationEnergy   &atomenergy,
                      const std::vector<double> &freq)
{
    real scale_factor = 1;
    real roomTemp     = 298.15;
    JsonTree jtc0("T=0");
    JsonTree jtcRT(gmx::formatString("T=%g", roomTemp));
    ThermoChemistry tc0(mol, *coords, atomenergy, freq, 0.0, 1, scale_factor);
    ThermoChemistry tcRT(mol, *coords, atomenergy, freq, roomTemp, 1, scale_factor);
    {
        std::string unit("kJ/mol");
        jtc0.addValueUnit("Zero point energy", gmx_ftoa(tc0.ZPE()), unit);
        jtcRT.addValueUnit("Zero point energy", gmx_ftoa(tcRT.ZPE()), unit);
    }
    {
        JsonTree jtS0("Standard entropy");
        JsonTree jtSRT("Standard entropy");
        for(const auto &tcc : tccmap())
        {
            std::string unit("J/mol K");
            jtS0.addValueUnit(tcc.second, gmx_ftoa(tc0.S0(tcc.first)), unit);
            jtSRT.addValueUnit(tcc.second, gmx_ftoa(tcRT.S0(tcc.first)), unit);
        }
        jtc0.addObject(jtS0);
        jtcRT.addObject(jtSRT);
    }
    {
        JsonTree jtcv0("Heat capacity cV");
        JsonTree jtcvRT("Heat capacity cV");
        for(const auto &tcc : tccmap())
        {
            std::string unit("J/mol K");
            jtcv0.addValueUnit(tcc.second, gmx_ftoa(tc0.cv(tcc.first)), unit);
            jtcvRT.addValueUnit(tcc.second, gmx_ftoa(tcRT.cv(tcc.first)), unit);
        }
        jtc0.addObject(jtcv0);
        jtcRT.addObject(jtcvRT);
    }
    {
        JsonTree jtcE0("Internal energy");
        JsonTree jtcERT("Internal energy");
        for(const auto &tcc : tccmap())
        {
            std::string unit("kJ/mol");
            jtcE0.addValueUnit(tcc.second, gmx_ftoa(tc0.Einternal(tcc.first)), unit);
            jtcERT.addValueUnit(tcc.second, gmx_ftoa(tcRT.Einternal(tcc.first)), unit);
        }
        jtc0.addObject(jtcE0);
        jtcRT.addObject(jtcERT);
    }
    {
        std::string unit("kJ/mol");
        jtc0.addValueUnit("Delta H formation", gmx_ftoa(tc0.DHform()), unit);
        jtcRT.addValueUnit("Delta H formation", gmx_ftoa(tcRT.DHform()), unit);
    }
    jtree->addObject(jtc0);
    jtree->addObject(jtcRT);
}

void doFrequencyAnalysis(const ForceField         *pd,
                         alexandria::ACTMol       *mol,
                         const MolHandler         &molhandler,
                         const ForceComputer      *forceComp,
                         std::vector<gmx::RVec>   *coords,
                         const AtomizationEnergy  &atomenergy,
                         gmx_stats                *lsq_freq_all,
                         JsonTree                 *jtree,
                         const char               *spectrumFileName,
                         double                    lineWidth,
                         gmx_output_env_t         *oenv,
                         bool                      useLapack,
                         bool                      debugNMA)
{
    std::vector<double> alex_freq, intensities;
    std::vector<std::string> output;
    molhandler.nma(pd, mol, forceComp, coords, &alex_freq, &intensities,
                   &output, useLapack, debugNMA);
    auto unit      = mpo_unit2(MolPropObservable::FREQUENCY);
    auto uniti     = mpo_unit2(MolPropObservable::INTENSITY);
    auto ref_freq  = mol->referenceFrequencies();
    auto ref_inten = mol->referenceIntensities();
    plot_spectrum(spectrumFileName, mol->getMolname(), lineWidth, alex_freq,
                  intensities, ref_freq, ref_inten, oenv, true);
    JsonTree ftree("Frequencies");
    JsonTree itree("Intensities");
    
    // Now loop over frequncies and intensities
    JsonTree aftree("Alexandria");
    JsonTree aitree("Alexandria");
    JsonTree rftree("Reference");
    JsonTree ritree("Reference");
    // Statistics
    gmx_stats lsq_freq, lsq_inten;
    for(size_t k = 0; k < alex_freq.size(); k++)
    {
        double fcalc = convertFromGromacs(alex_freq[k], unit);
        double icalc = convertFromGromacs(intensities[k], uniti);
        aftree.addValueUnit(gmx::formatString("%lu", k+1),
                            gmx_ftoa(fcalc), unit);
        aitree.addValueUnit(gmx::formatString("%lu", k+1),
                            gmx_ftoa(icalc), uniti);
        
        if (!ref_freq.empty())
        {
            double fref  = convertFromGromacs(ref_freq[k], unit);
            rftree.addValueUnit(gmx::formatString("%lu", k+1),
                                gmx_ftoa(fref), unit);
            
            if (nullptr != lsq_freq_all)
            {
                lsq_freq_all->add_point(fref, fcalc, 0, 0);
            }
            lsq_freq.add_point(fref, fcalc, 0, 0);
        }
        if (!ref_inten.empty())
        {
            double iref = convertFromGromacs(ref_inten[k], uniti);
            ritree.addValueUnit(gmx::formatString("%lu", k+1),
                                gmx_ftoa(iref), uniti);
            lsq_inten.add_point(iref, icalc, 0, 0);
        }
    }
    ftree.addObject(aftree);
    itree.addObject(aitree);
    if (!rftree.empty())
    {
        ftree.addObject(rftree);
        GMX_RELEASE_ASSERT(lsq_freq.get_npoints() > 0, "No statistics");
    
        real r, rmsd;
        if (eStats::OK == lsq_freq.get_corr_coeff(&r) &&
            eStats::OK == lsq_freq.get_rmsd(&rmsd))
        {
            ftree.addObject("Pearson r2", gmx_ftoa(100*r*r));
            ftree.addObject("RMSD", gmx_ftoa(rmsd));
        }
    }
    if (!ritree.empty())
    {
        itree.addObject(ritree);
        real r, rmsd;
        if (eStats::OK == lsq_inten.get_corr_coeff(&r) &&
            eStats::OK == lsq_inten.get_rmsd(&rmsd))
        {
            itree.addObject("Pearson r2", gmx_ftoa(100*r*r));
            itree.addObject("RMSD", gmx_ftoa(rmsd));
        }
    }
    jtree->addObject(ftree);
    jtree->addObject(itree);
    
    // Now time for thermochemistry output
    JsonTree tctree("Thermochemistry");
    JsonTree ajtc("Alexandria");
    addThermo(&ajtc, mol, coords, atomenergy, alex_freq);
    tctree.addObject(ajtc);
    if (!ref_freq.empty())
    {
        JsonTree rjtc("Reference");
        addThermo(&rjtc, mol, coords, atomenergy, ref_freq);
        tctree.addObject(rjtc);
    }
    jtree->addObject(tctree);
}

static void low_print_stats(std::vector<std::string> *tcout,
                            gmx_stats                *stats,
                            const char               *label)
{
    real mse, mae, rmsd, R = 0;
    stats->get_mse_mae(&mse, &mae);
    stats->get_rmsd(&rmsd);
    if (stats->get_npoints() > 2)
    {
        stats->get_corr_coeff(&R);
    }
    tcout->push_back(gmx::formatString("%s RMSD %8.2f MSE %8.2f (kJ/mol) R %4.1f%% #points = %zu",
                                       label, rmsd, mse, 100*R, stats->get_npoints()));
}

static void print_diatomics(const alexandria::ACTMol                                                  *mol,
                            const std::vector<std::pair<double, std::map<InteractionType, double> > > &energyComponentMap,
                            const gmx_output_env_t                                                    *oenv)
{
    std::vector<double> dist;
    for (const auto &ei : mol->experimentConst())
    {
        const std::vector<gmx::RVec> &xxx = ei.getCoordinates();
        gmx::RVec dx;
        rvec_sub(xxx[0], xxx[1], dx);
        dist.push_back(norm(dx));
    }
    auto filename = mol->getMolname() + ".xvg";
    FILE *fp = xvgropen(filename.c_str(), "QM - Coulomb Energy", "Distance (nm)", "Energy (kJ/mol)", oenv);
    size_t i = 0;
    for (const auto &ecm : energyComponentMap)
    {
        double eee = ecm.first;
        for(const auto itype : { InteractionType::COULOMB, InteractionType::POLARIZATION })
        {
            auto esm = ecm.second.find(itype);
            if (ecm.second.end() != esm)
            {
                eee -= esm->second;
            }
        }
        if (i < dist.size())
        {
            fprintf(fp, "%10g  %10g\n", dist[i], eee);
            i++;
        }
    }    
    xvgrclose(fp);
}
                                                
double TuneForceFieldPrinter::printEnergyForces(std::vector<std::string> *tcout,
                                                const ForceField         *pd,
                                                const ForceComputer      *forceComp,
                                                const AtomizationEnergy  &atomenergy,
                                                alexandria::ACTMol       *mol,
                                                gmx_stats                *lsq_rmsf,
                                                qtStats                  *lsq_epot,
                                                qtStats                  *lsq_eInter,
                                                gmx_stats                *lsq_freq,
                                                const gmx_output_env_t   *oenv)
{
    std::vector<ACTEnergy>                                  energyMap;
    std::vector<std::vector<std::pair<double, double> > >   forceMap;
    std::vector<std::pair<double, std::map<InteractionType, double> > > energyComponentMap, interactionEnergyMap;
    mol->forceEnergyMaps(pd, forceComp, &forceMap, &energyMap, &interactionEnergyMap,
                         &energyComponentMap);
    if (diatomic_ && mol->nRealAtoms() == 2)
    {
        print_diatomics(mol, energyComponentMap, oenv);
    }
    std::vector<std::string> dataFileNames;
    if (printSP_)
    {
        for (auto &ei : mol->experimentConst())
        {
            dataFileNames.push_back(ei.getDatafile());
        }
    }
    gmx_stats myepot;
    for(const auto &ff : energyMap)
    {
        (*lsq_epot)[qType::Calc].add_point(ff.eqm(), ff.eact(), 0, 0);
        myepot.add_point(ff.eqm(), ff.eact(), 0, 0);
    }
    gmx_stats myinter;
    for(const auto &ff : interactionEnergyMap)
    {
        auto eact = ff.second.find(InteractionType::EPOT)->second;
        (*lsq_eInter)[qType::Calc].add_point(ff.first, eact, 0, 0);
        myinter.add_point(ff.first, eact);
    }
    if (printSP_)
    {
        size_t ccc = 0;
        std::sort(energyComponentMap.begin(), energyComponentMap.end());
        for(auto eam = energyComponentMap.begin(); eam < energyComponentMap.end(); ++eam)
        {
            auto enerexp = eam->first;
            std::string ttt = gmx::formatString("%s Reference: %g ACT:", dataFileNames[ccc].c_str(),
                                                enerexp);
            for(const auto &terms : eam->second)
            {
                ttt += gmx::formatString(" %s: %g", interactionTypeToString(terms.first).c_str(), 
                                         terms.second);
            }
            ttt += gmx::formatString(" Diff: %g", eam->second[InteractionType::EPOT]-enerexp);
            tcout->push_back(ttt);
            ccc++;
        }
        std::sort(interactionEnergyMap.begin(), interactionEnergyMap.end());
        for(auto iem = interactionEnergyMap.begin(); iem < interactionEnergyMap.end(); ++iem)
        {
            auto eact = iem->second.find(InteractionType::EPOT)->second;
            std::string ttt = gmx::formatString("Reference Einteraction %g ACT %g Diff %g",
                                                iem->first, eact, eact-iem->first);
            std::map<InteractionType, const char *> terms = { 
                { InteractionType::COULOMB, "Coul." },
                { InteractionType::POLARIZATION, "Pol." },
                { InteractionType::DISPERSION, "Disp." },
                { InteractionType::REPULSION, "Rep." } };
            for(auto &term : terms)
            {
                auto tptr = iem->second.find(term.first);
                if (iem->second.end() != tptr)
                {
                    ttt += gmx::formatString(" %s %g", term.second, tptr->second);
                }
            }
            tcout->push_back(ttt);
        }
        std::sort(forceMap.begin(), forceMap.end());
        for(auto ff = forceMap.begin(); ff < forceMap.end(); ++ff)
        {
            for(const auto &fxyz : *ff)
            {
                if (fxyz.first != 0 || fxyz.second != 0)
                {
                    tcout->push_back(gmx::formatString("Force ref: %8.2f  act: %8.2f", fxyz.first, fxyz.second));
                    ccc++;
                }
            }
        }
    }
    // RMS energy
    if (energyMap.size() > 0)
    {
        low_print_stats(tcout, &myepot, "Energy");
    }
    if (interactionEnergyMap.size() > 0)
    {
        low_print_stats(tcout, &myinter, "Einteraction");
    }
    if (!forceMap.empty())
    {
        gmx_stats myforce;
        for(const auto &fstruct : forceMap)
        {
            for(auto &ff : fstruct)
            {
                if (ff.first != 0 || ff.second != 0)
                {
                    lsq_rmsf->add_point(ff.first, ff.second, 0, 0);
                    myforce.add_point(ff.first, ff.second, 0, 0);
                }
            }
        }
        low_print_stats(tcout, &myforce, "Force");
    }

    // Energy
    tcout->push_back(gmx::formatString("Energy terms (kJ/mol)"));
    std::map<InteractionType, double> eBefore;
    std::vector<gmx::RVec> coords = mol->xOriginal();
    std::vector<gmx::RVec> forces(coords.size());
    (void) forceComp->compute(pd, mol->topology(), &coords, &forces, &eBefore);

    double deltaE0 = 0;
    if (mol->energy(MolPropObservable::DELTAE0, &deltaE0))
    {
        //deltaE0 += mol->atomizationEnergy();
        tcout->push_back(gmx::formatString("   %-20s  %10.3f  Difference: %10.3f",
                                           "Reference EPOT", deltaE0, 
                                           eBefore[InteractionType::EPOT]-deltaE0));
    }
    if (mol->jobType() == JobType::OPT && calcFrequencies_)
    {
        // Now get the minimized structure RMSD and Energy
        SimulationConfigHandler simConfig;
        std::vector<gmx::RVec>  xmin   = coords;
        std::map<InteractionType, double> eAfter;
        molHandler_.minimizeCoordinates(pd, mol, forceComp, simConfig, 
                                        &xmin, &eAfter, nullptr, {});
        double rmsd = molHandler_.coordinateRmsd(mol, coords, &xmin);
        
        if (rmsd > 0.1) // nm
        {
            auto pdb   = gmx::formatString("inds/%s-original-minimized.pdb", mol->getMolname().c_str());
            auto title = gmx::formatString("%s RMSD %g Angstrom", mol->getMolname().c_str(), 10*rmsd);
            writeCoordinates(mol->atomsConst(), pdb, title, coords, xmin);
        }
        tcout->push_back(gmx::formatString("   %-20s  %10s  %10s  %10s minimization",
                                           "Term", "Before", "After", "Difference"));
        for(auto &ep : eBefore)
        {
            auto eb = ep.second;
            auto ea = eAfter[ep.first];
            tcout->push_back(gmx::formatString("   %-20s  %10.3f  %10.3f  %10.3f",
                                               interactionTypeToString(ep.first).c_str(),
                                               eb, ea, ea-eb));
        }
        tcout->push_back(gmx::formatString("Coordinate RMSD after minimization %10g pm", 1000*rmsd));
        
        // Do normal-mode analysis etc.
        JsonTree jtree("FrequencyAnalysis");
        doFrequencyAnalysis(pd, mol, molHandler_, forceComp, &coords,
                            atomenergy, lsq_freq, &jtree,
                            nullptr, 24, nullptr, false, false);
        int indent = 0;
        tcout->push_back(jtree.writeString(false, &indent));
    }
    else
    {
        for(auto &ep : eBefore)
        {
            tcout->push_back(gmx::formatString("   %-20s  %10.3f",
                                               interactionTypeToString(ep.first).c_str(),
                                               ep.second));
        }
    }
    return eBefore[InteractionType::EPOT];
}

static void dump_xyz(const ACTMol    *mol,
                     const ACTEnergy &actener)
{
    auto expconst = mol->experimentConst();
    if (static_cast<size_t>(actener.id()) < expconst.size())
    {
        auto filename = gmx::formatString("%s-%d.xyz", mol->getMolname().c_str(), actener.id());
        auto actatoms = mol->topology()->atoms();
        auto coords   = expconst[actener.id()].getCoordinates();
        FILE *fp = gmx_ffopen(filename.c_str(), "w");
        fprintf(fp, "%5d\n", mol->nRealAtoms());
        fprintf(fp, "Energy QM %g ACT %g QM datafile %s molname %s\n", actener.eqm(), actener.eact(),
                expconst[actener.id()].getDatafile().c_str(),
                mol->getMolname().c_str());
        int j = 0;
        for(int ra : mol->realAtoms())
        {
            // Convert coordinates to Angstrom
            fprintf(fp, "%3s  %12f  %12f  %12f\n", actatoms[ra].element().c_str(),
                    10*coords[j][XX], 10*coords[j][YY], 10*coords[j][ZZ]);
            j++;
        }
        gmx_ffclose(fp);
    }
}

static void printOutliers(FILE              *fp, 
                          iMolSelect         ims,
                          double             sigma,
                          const std::string &label,
                          real               dumpFactor,
                          const std::vector<alexandria::ACTMol>                *actmol,
                          const std::map<std::string, std::vector<ACTEnergy> > &allEpot)
{
    double epotMax = 1.5*sigma;
    if (allEpot.size() > 0)
    {
        fprintf(fp, "\nOverview of %s %s outliers for %s (Diff > %.3f)\n",
                iMolSelectName(ims), label.c_str(),
                qTypeName(qType::Calc).c_str(), epotMax);
        fprintf(fp, "----------------------------------\n");
        fprintf(fp, "%-40s  %12s  %12s  %12s\n", "Name",
                "Reference", qTypeName(qType::Calc).c_str(), "ACT-Ref.");
        int noutlier = 0;
        for (auto emm : allEpot)
        {
            std::string toFind(emm.first);
            auto actmolptr = std::find_if(actmol->begin(), actmol->end(),
                                          [&toFind](const ACTMol &x) { return x.getMolname() == toFind;});
            if (actmolptr->datasetType() == ims)
            {
                for (auto ener : emm.second)
                {
                    double deltaE = std::abs(ener.eqm()-ener.eact());
                    if (deltaE > epotMax)
                    {
                        fprintf(fp, "%-40s  %12g  %12g  %12g\n", emm.first.c_str(),
                                ener.eqm(), ener.eact(), ener.eact()-ener.eqm());
                        noutlier++;
                    }
                    if (dumpFactor >= 0 && deltaE >= dumpFactor*sigma)
                    {
                        if (actmol->end() != actmolptr)
                        {
                            dump_xyz(&(*actmolptr), ener);
                        }
                    }
                }
            }
        }
        if (noutlier)
        {
            printf("There were %d %s outliers for %s. Check the bottom of the log file\n", noutlier,
                   label.c_str(), iMolSelectName(ims));
        }
        else
        {
            printf("No %s outliers! Well done.\n", label.c_str());
        }
    }
}

void TuneForceFieldPrinter::print(FILE                            *fp,
                                  std::vector<alexandria::ACTMol> *actmol,
                                  const ForceField                *pd,
                                  const gmx_output_env_t           *oenv,
                                  const std::vector<t_filenm>      &filenm,
                                  const char                       *chargeMethod)
{
    int  n = 0;
    std::map<iMolSelect, qtStats>                               lsq_esp, lsq_alpha, lsq_isoPol,
        lsq_anisoPol, lsq_charge, lsq_epot, lsq_eInter;
    std::map<MolPropObservable, std::map<iMolSelect, qtStats> > lsq_multi;
    std::map<iMolSelect, std::map<std::string, gmx_stats> >     lsqt;
    std::map<iMolSelect, gmx_stats>                             lsq_rmsf, lsq_freq;
    for (auto &mpo : mpoMultiPoles)
    {
        std::map<iMolSelect, qtStats> mq_multi;
        for(auto &ims : iMolSelectNames())
        {
            qtStats qmpo;
            for (auto &i : qTypes())
            {
                gmx_stats gmult;
                qmpo.insert({ i.first, std::move(gmult) });
            }
            mq_multi.insert({ ims.first, qmpo });
        }
        lsq_multi.insert({ mpo, std::move(mq_multi) });
    }
    for(auto &ims : iMolSelectNames())
    {
        {
            gmx_stats rmsf;
            lsq_rmsf.insert({ ims.first, std::move(rmsf) });
            gmx_stats freq;
            lsq_freq.insert({ ims.first, std::move(freq) });
        }
        qtStats qesp, qepot, qinter;
        for (auto &i : qTypes())
        {
            //TODO Add checks for existence
            gmx_stats gesp;
            qesp.insert({ i.first, std::move(gesp) });
            gmx_stats gepot;
            qepot.insert({ i.first, std::move(gepot) });
            gmx_stats ginter;
            qinter.insert({ i.first, std::move(ginter) });
        }
        lsq_esp.insert({ ims.first, std::move(qesp) });
        lsq_epot.insert({ ims.first, std::move(qepot) });
        lsq_eInter.insert({ ims.first, std::move(qinter) });
        
        gmx_stats galpha;
        lsq_alpha.insert(std::pair<iMolSelect, qtStats>(ims.first, {{ qType::Calc, galpha }}));
        gmx_stats giso;
        lsq_isoPol.insert(std::pair<iMolSelect, qtStats>(ims.first, {{ qType::Calc, giso }}));
        gmx_stats ganiso;
        lsq_anisoPol.insert(std::pair<iMolSelect, qtStats>(ims.first, {{ qType::Calc, ganiso }}));
        gmx_stats gcharge;
        lsq_charge.insert(std::pair<iMolSelect, qtStats>(ims.first, {{ qType::Calc, gcharge }}));
    
        for (auto ai : pd->particleTypesConst())
        {
            auto qparam = ai.second.parameterConst("charge");
            if (Mutability::ACM == qparam.mutability())
            {
                auto lll = lsqt.find(ims.first);
                if (lll == lsqt.end())
                {
                    lsqt.insert({ ims.first, {} });
                    lll = lsqt.find(ims.first);
                }
                gmx_stats newz;
                lll->second.insert({ai.second.id().id(), std::move(newz)});
            }
        }
    }
    
    bool bPolar = pd->polarizable();
    std::map<MolPropObservable, double> multi_toler = {
        { MolPropObservable::DIPOLE,       dip_toler_  },
        { MolPropObservable::QUADRUPOLE,   quad_toler_ },
        { MolPropObservable::OCTUPOLE,     oct_toler_  },
        { MolPropObservable::HEXADECAPOLE, hex_toler_  },
    };
    
    auto forceComp = new ForceComputer();
    AtomizationEnergy atomenergy;
    std::map<std::string, double> molEpot;
    auto alg   = pd->chargeGenerationAlgorithm();
    auto qtype = qType::Calc;
    if (nullptr != chargeMethod && strlen(chargeMethod) > 0)
    {
        qtype = stringToQtype(chargeMethod);
        alg   = ChargeGenerationAlgorithm::Read;
    }
    std::map<std::string, std::vector<ACTEnergy > > allEpot; 
    std::map<std::string, std::vector<ACTEnergy > > allEinter; 
    for (auto mol = actmol->begin(); mol < actmol->end(); ++mol)
    {
        if (mol->support() != eSupport::No)
        {
            auto ims = mol->datasetType();
            fprintf(fp, "\nMolecule %d: Name: %s, Qtot: %d, Multiplicity: %d, MolWt: %g SymmetryNumber: %d Dataset: %s\n", n+1,
                    mol->getMolname().c_str(),
                    mol->totalCharge(),
                    mol->totalMultiplicity(),
                    mol->totalMass(),
                    mol->symmetryNumber(),
                    iMolSelectName(ims));

            // Recalculate the atomic charges using the optimized parameters.
            std::vector<double>    dummy;
            gmx::RVec vzero = { 0, 0, 0 };
            std::vector<gmx::RVec> forces(mol->atomsConst().size(), vzero);
            std::vector<gmx::RVec> coords = mol->xOriginal();
            mol->GenerateCharges(pd, forceComp, alg, qtype, dummy, &coords, &forces);
            // Now compute all the ESP RMSDs and multipoles and print it.
            fprintf(fp, "Electrostatic properties.\n");
            for (auto &i : qTypes())
            {
                for(auto qp = mol->qProps()->begin(); qp < mol->qProps()->end(); ++qp)
                {
                    auto qi    = i.first;
                    auto qresp = qp->qgenResp();
                    if (qresp->nEsp() > 0 && qp->qPactConst().qtype() == qi)
                    {
                        real rms, rrms, cosesp, mae, mse;
                        // Fetch coordinates and optimize shells if polarizable
                        auto myx = qresp->coords();
                        std::map<InteractionType, double> energies;
                        (void) forceComp->compute(pd, mol->topology(), &myx,
                                                  &forces, &energies);
                        qresp->updateAtomCoords(myx);
                        qresp->updateAtomCharges(mol->atomsConst());
                        qresp->calcPot(1.0);
                        rms = qresp->getStatistics(&rrms, &cosesp, &mae, &mse);
                        rms = convertToGromacs(rms, "Hartree/e");
                        std::string warning;
                        if (rms > esp_toler_ || cosesp < 0.5)
                        {
                            warning.assign(" EEE");
                        }
                        fprintf(fp, "ESP rms: %8.3f (kJ/mol e) rrms: %8.3f CosAngle: %6.3f - %s%s\n",
                                rms, rrms, cosesp, qTypeName(qi).c_str(), warning.c_str());   
                        if (mol->datasetType() == ims)
                        {
                            auto ep = qresp->espPoints();
                            for (size_t j = 0; j < ep.size(); j++)
                            {
                                lsq_esp[ims][qi].add_point(ep[j].v(), ep[j].vCalc(), 0, 0);
                            }
                        }
                    }
                }
            }
            // Charges
            auto atoms = mol->topology()->atoms();
            auto lll   = lsqt.find(ims);
            for(size_t ai = 0; ai < atoms.size(); ai++)
            {
                if (atoms[ai].pType() == eptAtom)
                {
                    auto llF = lll->second.find(atoms[ai].ffType());
                    if (lll->second.end() != llF)
                    {
                        auto qa = atoms[ai].charge();
                        if (ai < atoms.size() -1 && 
                            atoms[ai+1].pType() == eptShell)
                        {
                            qa += atoms[ai+1].charge();
                        }
                        llF->second.add_point(1, qa, 0, 0);
                    }
                }
            }
            // Multipoles
            analyse_multipoles(fp, mol, multi_toler, &lsq_multi);
            
            // Polarizability
            if (bPolar)
            {
                analysePolarisability(fp, pd, &(*mol), forceComp, &(lsq_isoPol[ims]),
                                      &(lsq_anisoPol[ims]), &(lsq_alpha[ims]));
            }

            // Atomic charges
            printAtoms(fp, &(*mol), coords, forces);
            // Energies
            std::vector<std::string> tcout;
            auto nepot   = lsq_epot[ims][qType::Calc].get_npoints();
            auto neInter = lsq_eInter[ims][qType::Calc].get_npoints();
            auto epot    = printEnergyForces(&tcout, pd, forceComp, atomenergy, &(*mol),
                                             &lsq_rmsf[ims], &lsq_epot[ims],
                                             &lsq_eInter[ims], &lsq_freq[ims], oenv);
            {
                // Copy the epot pairs
                auto x = lsq_epot[ims][qType::Calc].getX();
                auto y = lsq_epot[ims][qType::Calc].getY();
                if (x.size() > nepot)
                {
                    std::vector<ACTEnergy> myEpot;
                    for(size_t kk = nepot; kk < x.size(); kk++)
                    {
                        myEpot.push_back(ACTEnergy(kk-nepot, x[kk], y[kk]));
                    }
                    allEpot.insert({mol->getMolname(), myEpot});
                }
            }
            {
                // Copy the eInter pairs
                auto x = lsq_eInter[ims][qType::Calc].getX();
                auto y = lsq_eInter[ims][qType::Calc].getY();
                if (x.size() > neInter)
                {
                    std::vector<ACTEnergy> myEinter;
                    for(size_t kk = neInter; kk < x.size(); kk++)
                    {
                        myEinter.push_back(ACTEnergy(kk-neInter, x[kk], y[kk]));
                    }
                    allEinter.insert({ mol->getMolname(), myEinter });
                }
            }
            molEpot.insert({mol->getMolname(), epot});
            for(const auto &tout : tcout)
            {
                fprintf(fp, "%s\n", tout.c_str());
            }
            fprintf(fp, "\n");
            n++;
        }
    }

    for(auto &ims : iMolSelectNames())
    {
        bool header = true;
        
        fprintf(fp, "\n*** Results for %s data set ***\n", iMolSelectName(ims.first));
        for (auto &i : qTypes())
        {
            auto qt = i.first;
            if (qt == qType::Elec)
            {
                continue;
            }
            const char *name = qTypeName(qt).c_str();
            if (lsq_epot[ims.first][qt].get_npoints() > 0)
            {
                print_stats(fp, "Potential energy", "kJ/mol", 1.0, &lsq_epot[ims.first][qt],  header, "QM/DFT", name, useOffset_);
                header = false;
            }
            if (lsq_eInter[ims.first][qt].get_npoints() > 0)
            {
                print_stats(fp, "Interaction energy", "kJ/mol", 1.0, &lsq_eInter[ims.first][qt],  header, "QM/DFT", name, useOffset_);
                header = false;
            }
            if (lsq_rmsf[ims.first].get_npoints() > 0 && qt == qType::Calc)
            {
                print_stats(fp, "RMS Force", "kJ/mol nm", 1.0, &lsq_rmsf[ims.first], header, "QM/DFT", name, useOffset_);
                header = false;
            }
            if (lsq_freq[ims.first].get_npoints() > 0 && qt == qType::Calc)
            {
                print_stats(fp, "Frequencies", mpo_unit2(MolPropObservable::FREQUENCY),
                            1.0, &lsq_freq[ims.first], header, "QM/DFT", name, useOffset_);
                header = false;
            }
            print_stats(fp, "ESP", "kJ/mol e", 1.0, &lsq_esp[ims.first][qt],  header, "Electronic", name, useOffset_);
            header = false;
            for(auto &mpo : mpoMultiPoles)
            {
                print_stats(fp, mpo_name(mpo), mpo_unit2(mpo), 1.0, //convertFromGromacs(1.0, mpo_unit2(mpo)),
                            &lsq_multi[mpo][ims.first][qt],   header, "Electronic", name, useOffset_);
            }
            if (bPolar && qt == qType::Calc)
            {
                std::string polunit("Angstrom3");
                auto polfactor = convertFromGromacs(1.0, polunit);
                print_stats(fp, "Polariz. components", polunit.c_str(), polfactor,
                            &lsq_alpha[ims.first][qType::Calc],    header, "Electronic", name, useOffset_);
                print_stats(fp, "Isotropic Polariz.", polunit.c_str(), polfactor,
                            &lsq_isoPol[ims.first][qType::Calc],   header, "Electronic", name, useOffset_);
                print_stats(fp, "Anisotropic Polariz.", polunit.c_str(), polfactor,
                            &lsq_anisoPol[ims.first][qType::Calc], header, "Electronic", name, useOffset_);
            }
        }
    }
    write_q_histo(fp, opt2fn_null("-qhisto", filenm.size(), filenm.data()),
                  &lsqt[iMolSelect::Train], oenv,
                  &(lsq_charge[iMolSelect::Train]), useOffset_);

    for(auto &mpo : mpoMultiPoles)
    {
        std::string cmdFlag = gmx::formatString("-%scorr", mpo_name(mpo));
        std::string title   = gmx::formatString("%s components (%s)", mpo_name(mpo),
                                                mpo_unit2(mpo));
        print_corr(opt2fn_null(cmdFlag.c_str(), filenm.size(), filenm.data()),
                   title.c_str(), "Electronic", "Empirical", lsq_multi[mpo], oenv);
    }
    print_corr(opt2fn_null("-epotcorr", filenm.size(), filenm.data()),
               "Potential energy (kJ/mol)", "Reference", "Empirical",
               lsq_epot, oenv);
    print_corr(opt2fn_null("-eintercorr", filenm.size(), filenm.data()),
               "Interaction energy (kJ/mol)", "Reference", "Empirical",
               lsq_eInter, oenv);
    print_corr(opt2fn_null("-forcecorr", filenm.size(), filenm.data()),
               "Forces (kJ/mol nm)", "Reference", "Empirical",
               lsq_rmsf, oenv);
    print_corr(opt2fn_null("-freqcorr", filenm.size(), filenm.data()),
               "Frequencies (cm^-1)", "Reference", "Empirical",
               lsq_freq, oenv);
    print_corr(opt2fn_null("-espcorr", filenm.size(), filenm.data()),
               "Electrostatic Potential (Hartree/e)", "Electronic", "Calc",
               lsq_esp, oenv);
    print_corr(opt2fn_null("-qcorr", filenm.size(), filenm.data()),
               "Atomic Partial Charge", "q (e)", "a.u.", lsq_charge, oenv);
    
    if (bPolar)
    {
        print_corr(opt2fn_null("-alphacorr", filenm.size(), filenm.data()),
                   "Pricipal Components of Polarizability Tensor (A\\S3\\N)", "Electronic", "Calc", lsq_alpha, oenv);
        print_corr(opt2fn_null("-isopol", filenm.size(), filenm.data()),
                   "Isotropic Polarizability (A\\S3\\N)", "Electronic", "Calc", lsq_isoPol, oenv);
        print_corr(opt2fn_null("-anisopol", filenm.size(), filenm.data()),
                   "Anisotropic Polarizability (A\\S3\\N)", "Electronic", "Calc", lsq_anisoPol, oenv);
    }
    // List outliers based on the deviation in the Electrostatic Potential
    real espAver;
    if (eStats::OK == lsq_esp[iMolSelect::Train][qType::Calc].get_rmsd(&espAver))
    {
        int nout = 0;
        double espMax = 1.5*espAver;
        fprintf(fp, "\nOverview of ESP outliers for %s (RMSD > %.3f)\n",
                qTypeName(qType::Calc).c_str(), espMax);
        fprintf(fp, "----------------------------------\n");
        fprintf(fp, "%-40s  %12s  %12s\n", "Name",
                qTypeName(qType::Calc).c_str(), qTypeName(qType::ESP).c_str());
        for (auto mol = actmol->begin(); mol < actmol->end(); ++mol)
        {
            for(auto qp = mol->qProps()->begin(); qp < mol->qProps()->end(); ++qp)
            {
                real rms, rrms, cosesp, mae, mse;
                auto qresp = qp->qgenResp();
                rms        = convertToGromacs(qresp->getStatistics(&rrms, &cosesp, &mae, &mse), "Hartree/e");
                if ((mol->support() != eSupport::No) && (rms > espMax))
                {
                    fprintf(fp, "%-40s  %12.3f", mol->getMolname().c_str(), rms);
                    fprintf(fp, "  %s\n", iMolSelectName(mol->datasetType()));
                    nout++;
                }
            }
        }
        if (nout)
        {
            printf("There were %d ESP outliers. Check the bottom of the log file\n", nout);
        }
        else
        {
            printf("No ESP outliers! Well done.\n");
        }
    }
    real epotRmsd = 0;
    if (eStats::OK == lsq_epot[iMolSelect::Train][qType::Calc].get_rmsd(&epotRmsd))
    {
        for(const auto &ims : { iMolSelect::Train, iMolSelect::Test })
        {
            if (lsq_epot[ims][qType::Calc].get_npoints() > 0)
            {
                // List outliers based on the deviation in the Potential energy ...
                printOutliers(fp, ims, epotRmsd, "Epot", dumpOutliers_,
                              actmol, allEpot);
            }
        }
    }
    real einterRmsd = 0;
    if (eStats::OK == lsq_eInter[iMolSelect::Train][qType::Calc].get_rmsd(&einterRmsd))
    {
        for(const auto &ims : { iMolSelect::Train, iMolSelect::Test })
        {
            if (lsq_eInter[ims][qType::Calc].get_npoints() > 0)
            {
                // ... and the interaction energies.
                printOutliers(fp, ims, einterRmsd, "Einter", dumpOutliers_, 
                              actmol, allEinter);
            }
        }
    }
}

void print_header(FILE                        *fp, 
                  const std::vector<t_pargs>  &pargs,
                  const std::vector<t_filenm> &filenms)
{
    if (!fp)
    {
        return;
    }
    time_t my_t;
    time(&my_t);
    fprintf(fp, "# This file was created %s", ctime(&my_t));
    fprintf(fp, "# alexandria is the engine of the Alexandria Chemistry Toolkit\n#\n");
    fprintf(fp, "# https://github.com/dspoel/ACT\n#\n");
    for (auto &p: pargs)
    {
        std::string value, type;
        switch(p.type)
        {
        case etINT:
            type.assign("int");
            value = gmx::formatString("%d", *p.u.i);
            break;
        case etINT64:
            type.assign("int64");
            value = gmx::formatString("%" PRId64, *p.u.is);
            break;
        case etREAL:
            type.assign("real");
            value = gmx::formatString("%g", *p.u.r);
            break;
        case etSTR:
            type.assign("string");
            if (*p.u.c != nullptr)
            {
                value = *p.u.c;
            }
            break;
        case etENUM:
            type.assign("enum");
            if (*p.u.c != nullptr)
            {
                value = gmx::formatString("%s", *p.u.c);
            }
            break;
        case etBOOL:
            type.assign("bool");
            value = gmx::formatString("%s", *p.u.b ? "true" : "false");
            break;
        case etRVEC:
            type.assign("rvec");
            value = gmx::formatString("%g %g %g", *p.u.rv[XX],
                                      *p.u.rv[YY], *p.u.rv[ZZ]);
            break;
        case etTIME:
        case etNR:
        default:
            value.assign("help");
        }
        {
            std::string left = gmx::formatString("Option: %s (%s), Value: %s\n  Help: ", 
                                                 p.option, type.c_str(), value.c_str());
            gmx::TextLineWrapperSettings settings;
            settings.setLineLength(71);
            settings.setFirstLineIndent(0);
            settings.setIndent(8);
            gmx::TextLineWrapper wrapper(settings);
            std::string toWrap(p.desc);
            fprintf(fp, "%s%s\n", left.c_str(),
                    wrapper.wrapToString(toWrap).c_str());
        }
    }
    if (filenms.size() > 0)
    {
        fprintf(fp, "\nFiles used:\n");
    }
    std::map<unsigned int, const char *> fmap = {
        { ffREAD,  "R"   }, 
        { ffWRITE, "W"   },
        { ffSET,   "set" },
        { ffOPT,   "optional" }
    };
    for(const auto &f : filenms)
    {
        std::string flag;
        for(auto fm : fmap)
        {
            if ((fm.first & f.flag) == fm.first)
            {
                if (!flag.empty())
                {
                    flag += ",";
                }
                flag += fm.second;
            }
        }
        fprintf(fp, "Option: %s, Description: %s,", f.opt, ftp2desc(f.ftp));
        if (!flag.empty())
        {
            fprintf(fp, " Properties: %s, Filename(s)", flag.c_str());
        }
        if ((f.flag & ffMULT) == ffMULT)
        {
            for(auto &fnm : f.filenames)
            {
                fprintf(fp, " %s", fnm.c_str());
            }
        }
        else
        {
            if (nullptr != f.fn)
            {
                fprintf(fp, " %s", f.fn);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

} // namespace alexandria
