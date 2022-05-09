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

#include "tuning_utility.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/coolstuff.h"

#include "act/molprop/molprop_util.h"
#include "act/molprop/multipole_names.h"
#include "act/qgen/qtype.h"
#include "act/utility/units.h"
#include "alexandria/mymol.h"
#include "alexandria/thermochemistry.h"

void doAddOptions(std::vector<t_pargs> *pargs, size_t npa, t_pargs pa[])
{
    for (size_t i = 0; i < npa; i++)
    {
        pargs->push_back(pa[i]);
    }
}

namespace alexandria
{

class ZetaTypeLsq {
private:
    std::string ztype_;
public:
    gmx_stats lsq_;
    
    ZetaTypeLsq(const std::string &ztype) : ztype_(ztype) {}
    
    ZetaTypeLsq(const ZetaTypeLsq &zlsq)
    {
        lsq_   = zlsq.lsq_;
        ztype_ = zlsq.ztype_;
    }
    ZetaTypeLsq(ZetaTypeLsq &zlsq)
    {
        lsq_   = zlsq.lsq_;
        ztype_ = zlsq.ztype_;
    }
    ZetaTypeLsq& operator=(const ZetaTypeLsq &zlsq)
    {
        ZetaTypeLsq *zt = new ZetaTypeLsq(zlsq.name());
        zt->lsq_ = zlsq.lsq_;
        return *zt;
    }
    ZetaTypeLsq& operator=(ZetaTypeLsq &zlsq)
    {
        ZetaTypeLsq *zt = new ZetaTypeLsq(zlsq.name());
        zt->lsq_ = zlsq.lsq_;
        return *zt;
    }
        
    const std::string &name() const { return ztype_; }

    bool empty() const
    { 
        return 0 == lsq_.get_npoints();
    }
};

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

static void print_polarizability(FILE              *fp,
                                 const MyMol       *mol,
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
        auto qelec = mol->qTypeProps(qType::Elec);
        auto aelec = qelec->polarizabilityTensor();
        if (calc_name == qTypeName(qType::Calc))
        {
            auto qcalc = mol->qTypeProps(qType::Calc);
            auto acalc = qcalc->polarizabilityTensor();

            m_sub(aelec, acalc, dalpha);
            delta = fac*sqrt(gmx::square(dalpha[XX][XX])+gmx::square(dalpha[XX][YY])+
                             gmx::square(dalpha[XX][ZZ])+
                             gmx::square(dalpha[YY][YY])+gmx::square(dalpha[YY][ZZ]));
            diso_pol   = fac*std::abs(qcalc->isotropicPolarizability()-
                                      qelec->isotropicPolarizability());
            daniso_pol = fac*std::abs(qcalc->anisotropicPolarizability()-
                                      qelec->anisotropicPolarizability());
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
                    fac*qelec->isotropicPolarizability(),
                    fac*qcalc->isotropicPolarizability(), 
                    diso_pol, (diso_pol > isopol_toler) ? "ISO" : "");
            fprintf(fp,
                    "Anisotropic polarizability:  %s Electronic: %6.2f  Calculated: %6.2f  Delta: %6.2f %s\n\n",
                    mol->getMolname().c_str(), 
                    fac*qelec->anisotropicPolarizability(),
                    fac*qcalc->anisotropicPolarizability(), 
                    daniso_pol, (daniso_pol > isopol_toler) ? "ANISO" : "");

        }
        else
        {
            fprintf(fp, "Polarizability analysis\n");
            fprintf(fp,
                    "Electronic   (%6.2f %6.2f %6.2f)\n"
                    "             (%6s %6.2f %6.2f)\n"
                    "             (%6s %6s %6.2f)\n",
                    fac*aelec[XX][XX], fac*aelec[XX][YY], fac*aelec[XX][ZZ],
                    "", fac*aelec[YY][YY], fac*aelec[YY][ZZ],
                    "", "", fac*aelec[ZZ][ZZ]);
        }
    }
}

static void analyse_multipoles(FILE                                                        *fp,
                               const std::vector<alexandria::MyMol>::iterator              &mol,
                               std::map<MolPropObservable, double>                          toler,
                               std::map<MolPropObservable, std::map<iMolSelect, qtStats> > *lsq)
{
    auto qelec = mol->qTypeProps(qType::Elec);
    auto qcalc = mol->qTypeProps(qType::Calc);
    qcalc->calcMoments();
    for(auto &mpo : mpoMultiPoles)
    {
        const char *name   = mpo_name(mpo);
        const char *unit   = mpo_unit2(mpo);
        double      factor = convertFromGromacs(1, unit);
        fprintf(fp, "Electronic %s (%s):\n", name, unit);
        auto Telec = qelec->getMultipole(mpo);
        printMultipole(fp, mpo, Telec);

        for (auto &j : qTypes())
        {
            qType qt = j.first;
            if (qt == qType::Elec)
            {
                continue;
            }
            auto qcalc = mol->qTypeProps(qt);
            if (qcalc)
            {
                real delta = 0;
                auto Tcalc = qcalc->getMultipole(mpo);
                fprintf(fp, "%s %s (%s):\n",
                        qTypeName(qt).c_str(), name, unit);
                printMultipole(fp, mpo, Tcalc);

                std::vector<double> diff;
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
                fprintf(fp, "Difference %s Norm %g RMS = %g (%s)%s:\n",
                        name, factor*std::sqrt(delta), factor*rms, unit, flag.c_str());
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

static void write_q_histo(FILE                      *fplog,
                          const char                *qhisto,
                          std::vector<ZetaTypeLsq>  &lsqt,
                          const gmx_output_env_t    *oenv,
                          qtStats                   *lsq_charge,
                          bool                       useOffset)
{
    std::vector<std::string> types;
    for (auto i = lsqt.begin(); i < lsqt.end(); ++i)
    {
        types.push_back(i->name());
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
    for (auto &q : lsqt)
    {
        if (!q.empty())
        {
            int  nbins = 20;
            std::vector<double> x, y;
            if (q.lsq_.make_histogram(0, &nbins, eHisto::Y, 1, &x, &y) == eStats::OK)
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
                print_stats(fplog, q.name().c_str(), "e", 1.0, &q.lsq_, false,  "CM5", model.c_str(), useOffset);
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
          "Perform energy minimization and compute vibrational frequencies for each molecule (after optimizing the force field if -optimize is enabled). If turned on, this option will also generate thermochemistry values based on the force field." }
    };
    doAddOptions(pargs, sizeof(pa)/sizeof(pa[0]), pa);
}

void TuneForceFieldPrinter::addFileOptions(std::vector<t_filenm> *filenm)
{
    std::vector<t_filenm> fnm = {
        { efXVG, "-qhisto",    "q_histo",       ffOPTWR },
        { efXVG, "-epotcorr",  "epot_corr",     ffOPTWR },
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

void TuneForceFieldPrinter::analysePolarisability(FILE              *fp,
                                                  alexandria::MyMol *mol,
                                                  qtStats           *lsq_isoPol,
                                                  qtStats           *lsq_anisoPol,
                                                  qtStats           *lsq_alpha,
                                                  real               efield)
{
    auto qelec = mol->qTypeProps(qType::Elec);
    auto aelec = qelec->polarizabilityTensor();
    mol->CalcPolarizability(efield);
    auto qcalc = mol->qTypeProps(qType::Calc);
    auto acalc = qcalc->polarizabilityTensor();
    
    print_polarizability(fp, mol, qTypeName(qType::Elec), alpha_toler_, isopol_toler_);
    print_polarizability(fp, mol, qTypeName(qType::Calc), alpha_toler_, isopol_toler_);
    (*lsq_isoPol)[qType::Calc].add_point(qelec->isotropicPolarizability(),
                                         qcalc->isotropicPolarizability(),
                                         0, 0);
    (*lsq_anisoPol)[qType::Calc].add_point(qelec->anisotropicPolarizability(),
                                           qcalc->anisotropicPolarizability(),
                                           0, 0);
    for (int mm = 0; mm < DIM; mm++)
    {
        (*lsq_alpha)[qType::Calc].add_point(aelec[mm][mm], acalc[mm][mm], 0, 0);
    }
}

void TuneForceFieldPrinter::printAtoms(FILE              *fp,
                                       alexandria::MyMol *mol)
{
    std::map<qType, std::vector<double> > qQM;
    std::vector<qType>                    typeQM = { 
        qType::CM5, qType::ESP, qType::Hirshfeld, qType::Mulliken
    };
    for(auto &qt : typeQM)
    {
        auto qp = mol->qTypeProps(qt);
        if (qp)
        {
            std::vector<double> qqm;
            qqm = qp->charge();
            qQM.insert(std::pair<qType, std::vector<double>>(qt, qqm));
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
    auto x       = mol->x();
    int  i       = 0;
    double qtot  = 0;
    auto myatoms = mol->atomsConst();
    auto force   = mol->f();
    for (int j = i = 0; j < myatoms.nr; j++)
    {
        if (myatoms.atom[j].ptype == eptAtom ||
            myatoms.atom[j].ptype == eptNucleus)
        {
            real qCalc = myatoms.atom[j].q;
            fprintf(fp, "%-2d%3d  %-5s  %12g",
                    myatoms.atom[j].atomnumber,
                    j+1,
                    *(myatoms.atomtype[j]),
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
                    convertFromGromacs(x[j][XX], "pm"),
                    convertFromGromacs(x[j][YY], "pm"),
                    convertFromGromacs(x[j][ZZ], "pm"),
                    force[j][XX], force[j][YY], force[j][ZZ],
                    qtot);
            i++;
        }
        else
        {
            // Turned on printing of shells again
            fprintf(fp, "%-2d%3d  %-5s  %12g",
                    0,
                    j+1,
                    *(myatoms.atomtype[j]),
                    myatoms.atom[j].q);
            qtot += myatoms.atom[j].q;
            for(auto &qt : qQM)
            {
                if (!qQM.find(qt.first)->second.empty())
                {
                    fprintf(fp, "          ");
                }
            }
            fprintf(fp," %8.3f %8.3f %8.3f %10.3f %10.3f %10.3f  %10g\n", 
                    convertFromGromacs(x[j][XX], "pm"),
                    convertFromGromacs(x[j][YY], "pm"),
                    convertFromGromacs(x[j][ZZ], "pm"),
                    force[j][XX], force[j][YY], force[j][ZZ],
                    qtot);
        }
    }
}

void TuneForceFieldPrinter::printEnergyForces(std::vector<std::string> *tcout,
                                              alexandria::MyMol        *mol,
                                              const std::vector<int>   &ePlot,
                                              gmx_stats                *lsq_rmsf,
                                              qtStats                  *lsq_epot,
                                              gmx_stats                *lsq_freq)
{
    std::vector<std::pair<double, double> > eMap;
    std::vector<std::vector<std::pair<double, double> > > fMap;
    mol->forceEnergyMaps(&fMap, &eMap);
    std::vector<std::string> dataFileNames;
    if (printSP_)
    {
        for (auto &ei : mol->experimentConst())
        {
            dataFileNames.push_back(ei.getDatafile());
        }
    }
    double de2 = 0;
    size_t ccc = 0;
    for(const auto &ff : eMap)
    {
        auto enerexp = mol->atomizationEnergy() + ff.first;
        (*lsq_epot)[qType::Calc].add_point(enerexp, ff.second, 0, 0);
        de2 += gmx::square(enerexp - ff.second);
        if (printSP_)
        {
            tcout->push_back(gmx::formatString("%s Reference EPOT %g Calculated %g Difference %g",
                                               dataFileNames[ccc].c_str(),
                                               enerexp, ff.second, enerexp - ff.second));
        }
    }
    // RMS energy
    if (eMap.size() > 0)
    {
        tcout->push_back(gmx::formatString("RMS energy  %g (kJ/mol) #structures = %zu",
                                           std::sqrt(de2/eMap.size()), eMap.size()));
    }
    double df2 = 0;    
    for(const auto &fstruct : fMap)
    {
        for(const auto &ff : fstruct)
        {
            lsq_rmsf->add_point(ff.first, ff.second, 0, 0);
            df2 += gmx::square(ff.first - ff.second);
        }
    }
    // RMS force
    if (fMap.size() > 0)
    {
        tcout->push_back(gmx::formatString("RMS force   %g (kJ/mol nm) #structures = %zu",
                                           std::sqrt(df2/(mol->nRealAtoms())), fMap.size()));
    }
    // Energy
    tcout->push_back(gmx::formatString("Energy terms (kJ/mol, EPOT including atomization terms)"));
    std::vector<double> eBefore;
    auto terms = mol->energyTerms();
    for(int ii = 0; ii < F_NRE; ii++)
    {
        eBefore.push_back(terms[ii]);
    }
    double deltaE0 = 0;
    if (mol->energy(MolPropObservable::DELTAE0, &deltaE0))
    {
        deltaE0 += mol->atomizationEnergy();
        tcout->push_back(gmx::formatString("   %-20s  %10.3f  Difference: %10.3f",
                                           "Reference EPOT", deltaE0, eBefore[F_EPOT]-deltaE0));
    }
    if (mol->jobType() == JobType::OPT && calcFrequencies_)
    {
        double rmsd = 0;
        // Now get the minimized structure RMSD and Energy
        // TODO: Only do this for JobType::OPT
        molHandler_.minimizeCoordinates(mol, &rmsd);
        auto eAfter = mol->energyTerms();
        tcout->push_back(gmx::formatString("   %-20s  %10s  %10s  %10s minimization",
                                           "Term", "Before", "After", "Difference"));
        for(auto &ep : ePlot)
        {
            tcout->push_back(gmx::formatString("   %-20s  %10.3f  %10.3f  %10.3f",
                                               interaction_function[ep].name,
                                               eBefore[ep], eAfter[ep], eAfter[ep]-eBefore[ep]));
        }
        tcout->push_back(gmx::formatString("Coordinate RMSD after minimization %10g pm", 1000*rmsd));
        
        // Do normal-mode analysis
        std::vector<double> frequencies, intensities;
        molHandler_.nma(mol, &frequencies, &intensities, nullptr);
        auto unit = mpo_unit2(MolPropObservable::FREQUENCY);
        auto ref_freq = mol->referenceFrequencies();
        if (!ref_freq.empty())
        {
            tcout->push_back(gmx::formatString("Frequencies (%s)", mpo_unit2(MolPropObservable::FREQUENCY)));
            tcout->push_back(gmx::formatString("%10s  %10s", "Reference", "Alexandria"));
            for(size_t k = 0; k < frequencies.size(); k++)
            {
                double fref  = convertFromGromacs(ref_freq[k], unit);
                double fcalc = convertFromGromacs(frequencies[k], unit);
                lsq_freq->add_point(fref, fcalc, 0, 0);
                tcout->push_back(gmx::formatString("%10g  %10g", fref, fcalc));
            }
        }
        real scale_factor = 1;
        ThermoChemistry tc0(mol, frequencies, 0.0, 1, scale_factor);
        ThermoChemistry tcRT(mol, frequencies, 298.15, 1, scale_factor);
        tcout->push_back(gmx::formatString("%-30s  %10s  %10s", "Thermochemistry prediction", "0 K", "298.15 K"));
        
        tcout->push_back(gmx::formatString("%-30s  %10g  %10g (kJ/mol)", "Zero point energy", tc0.ZPE(), tcRT.ZPE()));
        tcout->push_back(gmx::formatString("%-30s  %10g  %10g (kJ/mol)", "Delta H formation", tc0.DHform(), tcRT.DHform()));
        for(const auto &tcc : tccmap())
        {
            tcout->push_back(gmx::formatString("%-30s  %10g  %10g (J/mol K)",
                                               gmx::formatString("Standard entropy - %s",
                                                                 tcc.second.c_str()).c_str(), 
                                               tc0.S0(tcc.first), tcRT.S0(tcc.first)));
        }
        for(const auto &tcc : tccmap())
        {
            tcout->push_back(gmx::formatString("%-30s  %10g  %10g (J/mol K)",
                                               gmx::formatString("Heat capacity cV - %s",
                                                                 tcc.second.c_str()).c_str(),
                                               tc0.cv(tcc.first), tcRT.cv(tcc.first)));
        }
        for(const auto &tcc : tccmap())
        {
            tcout->push_back(gmx::formatString("%-30s  %10g  %10g (kJ/mol)",
                                               gmx::formatString("Internal energy - %s",
                                                                 tcc.second.c_str()).c_str(),
                                               tc0.Einternal(tcc.first), tcRT.Einternal(tcc.first)));
        }
    }
    else
    {
        for(auto &ep : ePlot)
        {
            tcout->push_back(gmx::formatString("   %-20s  %10.3f",
                                               interaction_function[ep].name, eBefore[ep]));
        }
    }
}

void TuneForceFieldPrinter::print(FILE                           *fp,
                                  std::vector<alexandria::MyMol> *mymol,
                                  const Poldata                  *pd,
                                  const gmx::MDLogger            &fplog,
                                  const gmx_output_env_t         *oenv,
                                  const CommunicationRecord      *cr,
                                  real                            efield,
                                  const std::vector<t_filenm>    &filenm)
{
    int  n = 0;
    std::map<iMolSelect, qtStats>                               lsq_esp, lsq_alpha, lsq_isoPol,
                                                                lsq_anisoPol, lsq_charge, lsq_epot;
    std::map<MolPropObservable, std::map<iMolSelect, qtStats> > lsq_multi;
    std::map<iMolSelect, std::vector<ZetaTypeLsq> >             lsqt;
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
        qtStats qesp, qepot;
        for (auto &i : qTypes())
        {
            //TODO Add checks for existence
            gmx_stats gesp;
            qesp.insert({ i.first, std::move(gesp) });
            gmx_stats gepot;
            qepot.insert({ i.first, std::move(gepot) });
        }
        lsq_esp.insert({ ims.first, std::move(qesp) });
        lsq_epot.insert({ ims.first, std::move(qepot) });
        
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
            auto qparam = ai.parameterConst("charge");
            if (qparam.isMutable())
            {
                ZetaTypeLsq newz(ai.id().id());
                auto lll = lsqt.find(ims.first);
                if (lll == lsqt.end())
                {
                    std::vector<ZetaTypeLsq> zlsq = { std::move(newz) };
                    lsqt.insert(std::pair<iMolSelect, std::vector<ZetaTypeLsq>>(ims.first, zlsq));
                }
                else
                {
                    lll->second.push_back(std::move(newz));
                }
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
    // Extract terms to print for the enery terms
    std::vector<int> ePlot = { F_EPOT, F_COUL_SR, F_ATOMIZATION };
    {
        std::vector<InteractionType> includeTerms = { InteractionType::BONDS, InteractionType::ANGLES,
            InteractionType::LINEAR_ANGLES, InteractionType::PROPER_DIHEDRALS,
            InteractionType::IMPROPER_DIHEDRALS, InteractionType::VDW,
            InteractionType::POLARIZATION };
        for (const auto &ii : includeTerms)
        {
            if (pd->interactionPresent(ii))
            {
                ePlot.push_back(pd->findForcesConst(ii).fType());
            }
        }
    }
    for (auto mol = mymol->begin(); mol < mymol->end(); ++mol)
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
            std::vector<double> dummy;
            mol->GenerateCharges(pd, fplog, cr,
                                 ChargeGenerationAlgorithm::NONE, dummy);
            // Now compute all the ESP RMSDs and multipoles and print it.
            fprintf(fp, "Electrostatic properties.\n");
            mol->calcEspRms(pd);
            for (auto &i : qTypes())
            {
                auto qi  = i.first;
                if (qi == qType::Elec)
                {
                    continue;
                }
                auto qp = mol->qTypeProps(qi);
                if (nullptr == qp)
                {
                    continue;
                }
                real rms, rrms, cosesp, mae, mse;
                rms = qp->qgenResp()->getStatistics(&rrms, &cosesp, &mae, &mse);
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
                    auto ep = qp->qgenResp()->espPoint();
                    for (size_t j = 0; j < ep.size(); j++)
                    {
                        lsq_esp[ims][qi].add_point(ep[j].v(), ep[j].vCalc(), 0, 0);
                    }
                }
            }
            // Multipoles
            analyse_multipoles(fp, mol, multi_toler, &lsq_multi);
            
            // Polarizability
            if (bPolar)
            {
                analysePolarisability(fp, &(*mol), &(lsq_isoPol[ims]),
                                      &(lsq_anisoPol[ims]), &(lsq_alpha[ims]), efield);
            }

            // Atomic charges
            printAtoms(fp, &(*mol));
            // Energies
            std::vector<std::string> tcout;
            printEnergyForces(&tcout, &(*mol), ePlot,
                              &lsq_rmsf[ims], &lsq_epot[ims],
                              &lsq_freq[ims]);
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
                  lsqt[iMolSelect::Train], oenv,
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
        for (auto mol = mymol->begin(); mol < mymol->end(); ++mol)
        {
            real rms, rrms, cosesp, mae, mse;
            auto qcalc = mol->qTypeProps(qType::Calc);
            rms        = convertToGromacs(qcalc->qgenResp()->getStatistics(&rrms, &cosesp, &mae, &mse), "Hartree/e");
            if ((mol->support() != eSupport::No) && (rms > espMax))
            {
                fprintf(fp, "%-40s  %12.3f", mol->getMolname().c_str(), rms);
                auto qesp = mol->qTypeProps(qType::ESP);
                if (qesp)
                {
                    real rr, ce, mae, mse;
                    fprintf(fp, "  %12.3f",
                            convertToGromacs(qesp->qgenResp()->getStatistics(&rr, &ce, &mae, &mse), "Hartree/e"));
                }
                else
                {
                    fprintf(fp, "           N/A");
                }
                fprintf(fp, "  %s\n", iMolSelectName(mol->datasetType()));
                nout++;
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
    // List outliers based on the deviation in the Potential energy
    real epotAver;
    if (eStats::OK == lsq_epot[iMolSelect::Train][qType::Calc].get_rmsd(&epotAver))
    {
        int nout = 0;
        double epotMax = 1.5*epotAver;
        fprintf(fp, "\nOverview of Epot outliers for %s (Diff > %.3f)\n",
                qTypeName(qType::Calc).c_str(), epotMax);
        fprintf(fp, "----------------------------------\n");
        fprintf(fp, "%-40s  %12s  %12s\n", "Name",
                qTypeName(qType::Calc).c_str(), "QM/DFT");
        for (auto mol = mymol->begin(); mol < mymol->end(); ++mol)
        {
            double deltaE0;
            if (!mol->energy(MolPropObservable::DELTAE0, &deltaE0))
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("No molecular energy for %s",
                                                               mol->getMolname().c_str()).c_str()));
            }
            deltaE0 += mol->atomizationEnergy();
            auto diff = std::abs(mol->potentialEnergy()-deltaE0);
            if ((mol->support() != eSupport::No) && (diff > epotMax))
            {
                fprintf(fp, "%-40s  %12.3f  %12.3f  %s\n",
                        mol->getMolname().c_str(), mol->potentialEnergy(), deltaE0,
                        iMolSelectName(mol->datasetType()));
                nout++;
            }
        }
        if (nout)
        {
            printf("There were %d EPot outliers. Check the bottom of the log file\n", nout);
        }
        else
        {
            printf("No Epot outliers! Well done.\n");
        }
    }
}

void print_header(FILE                       *fp, 
                  const std::vector<t_pargs> &pargs)
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
}

} // namespace alexandria
