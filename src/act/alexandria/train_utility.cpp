/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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

#include "train_utility.h"

#include <algorithm>
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
#include "act/molprop/molprop_xml.h"
#include "act/molprop/multipole_names.h"
#include "act/qgen/qtype.h"
#include "act/utility/stringutil.h"
#include "act/utility/units.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/textwriter.h"

void doAddOptions(std::vector<t_pargs> *pargs, size_t npa, t_pargs pa[])
{
    for (size_t i = 0; i < npa; i++)
    {
        pargs->push_back(pa[i]);
    }
}

namespace alexandria
{

static void print_stats(gmx::TextWriter *tw,
                        const char      *prop,
                        const char      *unit,
                        double           conversion,
                        gmx_stats       *lsq,
                        bool             bHeader,
                        const char      *xaxis,
                        const char      *yaxis,
                        bool             useOffset)
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
            tw->writeStringFormatted("Fitting data to y = ax + b, where x = %s\n", xaxis);
            tw->writeStringFormatted("%-32s %6s %13s %13s %7s %8s %8s %8s %10s\n",
                    "Property", "N", "a", "b", "R(%)", "RMSD", "MSE", "MAE", "Dataset");
            tw->writeStringFormatted("------------------------------------------------------------------------------------------------\n");
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
            tw->writeStringFormatted("%-32s %6d %6.3f(%5.3f) %6.3f(%5.3f) %7.2f %8.4f %8.4f %8.4f %10s\n",
                    newprop.c_str(), n, a, da, 
                    b, db, Rfit*100, 
                    conversion*rmsd, conversion*mse, 
                    conversion*mae, yaxis);
        }
        else
        {
            tw->writeStringFormatted("Statistics problem for %s: %s\n", prop, gmx_stats_message(ok));
        }
    }
    else
    {
        if (bHeader)
        {
            tw->writeStringFormatted("Fitting data to y = ax, where x = %s\n", xaxis);
            tw->writeStringFormatted("%-32s %6s %13s %7s %8s %8s %8s %10s\n",
                    "Property", "N", "a", "R(%)", "RMSD", "MSE", "MAE", "Model");
            tw->writeStringFormatted("----------------------------------------------------------------------------------------------\n");
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
            tw->writeStringFormatted("%-32s %6d %6.3f(%5.3f) %7.2f %8.4f %8.4f %8.4f %10s\n",
                    newprop.c_str(), n, a, da, Rfit*100,
                    conversion*rmsd, conversion*mse, conversion*mae, yaxis);
        }
        else
        {
            tw->writeStringFormatted("Statistics problem for %s: %s\n", prop, gmx_stats_message(ok));
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
        fprintf(xvgf, "@ s%zu symbol %zu\n", i, i+1);
    }
}

static void print_one_alpha(gmx::TextWriter *tw,
                            const std::string &label,
                            const tensor alpha)
{
    double fac = convertFromGromacs(1.0, "Angstrom3");
    tw->writeStringFormatted("%s   (%6.2f %6.2f %6.2f)\n"
                             "             (%6s %6.2f %6.2f)\n"
                             "             (%6s %6s %6.2f)\n",
                             label.c_str(),
                             fac*alpha[XX][XX], fac*alpha[XX][YY], fac*alpha[XX][ZZ],
                             "", fac*alpha[YY][YY], fac*alpha[YY][ZZ],
                             "", "", fac*alpha[ZZ][ZZ]);
 }

static void print_polarizability(gmx::TextWriter   *tw,
                                 const std::string &molname,
                                 const QtypeProps  *qelec,
                                 const QtypeProps  *qcalc,
                                 real               alpha_toler,
                                 real               isopol_toler)
{
    tensor dalpha;
    real   delta      = 0;
    real   diso_pol   = 0;
    real   daniso_pol = 0;

    double fac = convertFromGromacs(1.0, "Angstrom3");
    auto acalc = qcalc->polarizabilityTensor();
    if (!qelec->hasPolarizability())
    {
        print_one_alpha(tw, "Alexandria", acalc);
    }
    else
    {
        tensor aelec;
        copy_mat(qelec->polarizabilityTensor(), aelec);
        print_one_alpha(tw, "Electronic", aelec);
        
        m_sub(aelec, acalc, dalpha);
        delta = fac*sqrt(gmx::square(dalpha[XX][XX])+gmx::square(dalpha[XX][YY])+
                         gmx::square(dalpha[XX][ZZ])+
                         gmx::square(dalpha[YY][YY])+gmx::square(dalpha[YY][ZZ]));
        diso_pol   = fac*std::abs(qcalc->isotropicPolarizability()-
                                  qelec->isotropicPolarizability());
        daniso_pol = fac*std::abs(qcalc->anisotropicPolarizability()-
                                  qelec->anisotropicPolarizability());
        tw->writeStringFormatted("%-12s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                                 "             (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                                 "             (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                                 "Alexandria",
                                 fac*acalc[XX][XX], fac*acalc[XX][YY], fac*acalc[XX][ZZ],
                                 fac*dalpha[XX][XX], fac*dalpha[XX][YY], fac*dalpha[XX][ZZ], delta, (delta > alpha_toler) ? "ALPHA" : "",
                                 "", fac*acalc[YY][YY], fac*acalc[YY][ZZ],
                                 "", fac*dalpha[YY][YY], fac*dalpha[YY][ZZ],
                                 "", "", fac*acalc[ZZ][ZZ],
                                 "", "", fac*dalpha[ZZ][ZZ]);
        tw->writeStringFormatted("Isotropic polarizability:  %s Electronic: %6.2f  Calculated: %6.2f  Delta: %6.2f %s\n\n",
                                 molname.c_str(),
                                 fac*qelec->isotropicPolarizability(),
                                 fac*qcalc->isotropicPolarizability(),
                                 diso_pol, (diso_pol > isopol_toler) ? "ISO" : "");
        tw->writeStringFormatted("Anisotropic polarizability:  %s Electronic: %6.2f  Calculated: %6.2f  Delta: %6.2f %s\n\n",
                                 molname.c_str(),
                                 fac*qelec->anisotropicPolarizability(),
                                 fac*qcalc->anisotropicPolarizability(),
                                 daniso_pol, (daniso_pol > isopol_toler) ? "ANISO" : "");
    }
}

TrainForceFieldPrinter::TrainForceFieldPrinter()
{
    terms_ = { InteractionType::EPOT,
               InteractionType::ELECTROSTATICS,
               InteractionType::INDUCTION, 
               InteractionType::INDUCTIONCORRECTION,
               InteractionType::DISPERSION,
               InteractionType::EXCHANGE,
               InteractionType::EXCHIND,
               InteractionType::ALLELEC };
}

void TrainForceFieldPrinter::analyse_multipoles(MsgHandler                                      *msg_handler,
                                                const std::vector<alexandria::ACTMol>::iterator &mol,
                                                std::map<MolPropObservable, double>              toler,
                                                const ForceField                                *pd,
                                                const ForceComputer                             *forceComputer)
{
    gmx::TextWriter *tw = nullptr;
    if (msg_handler)
    {
        tw = msg_handler->tw();
    }
    auto topology = mol->topology();
    bool doForce  = pd->polarizable() || topology->hasVsites();
    auto qprops = mol->qPropsConst();
    for(auto qp = qprops.begin(); qp < qprops.end(); ++qp)
    {
        auto qelec = qp->qPqmConst();
        auto qcalc = qp->qPact();
        qcalc->initializeMoments();
        if (doForce)
        {
            std::vector<gmx::RVec>            forces(topology->nAtoms());
            std::map<InteractionType, double> energies;
            auto                              myx = qcalc->x();
            forceComputer->compute(msg_handler, pd, topology, &myx, &forces, &energies);
            qcalc->setX(myx);
        }
        qcalc->calcMoments(msg_handler);

        for(auto &mpo : mpoMultiPoles)
        {
            const char *name   = mpo_name(mpo);
            const char *unit   = mpo_unit2(mpo);
            double      factor = convertFromGromacs(1, unit);
            std::vector<double> Telec;
            if (qelec.hasMultipole(mpo))
            {
                tw->writeStringFormatted("Electronic %s (%s):\n", name, unit);
                Telec = qelec.getMultipole(mpo);
                for(const auto &pm : formatMultipole(mpo, Telec))
                {
                    tw->writeLine(pm);
                }
            }
            real delta = 0;
            auto Tcalc = qcalc->getMultipole(mpo);
            tw->writeStringFormatted("Calc %s (%s):\n", name, unit);
            for(const auto &pm : formatMultipole(mpo, Tcalc))
            {
                tw->writeLine(pm);
            }

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
                        lsq_multi[mpo][mol->datasetType()][qt].add_point(factor*te, factor*tc, 0, 0);
                }
                double rms = std::sqrt(delta/Tcalc.size());
                std::string flag("");
                if (rms > toler[mpo])
                {
                    flag = " MULTI";
                }
                tw->writeStringFormatted("%s-Electronic Norm %g RMS = %g (%s)%s:\n",
                        qTypeName(qt).c_str(), factor*std::sqrt(delta), factor*rms, unit, flag.c_str());
                for(const auto &pm : formatMultipole(mpo, diff))
                {
                    tw->writeLine(pm);
                }
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
                    eprnm.push_back(gmx::formatString("%s", ims.second));
                    lsq.push_back(i.second);
                    haveData = true;
                }
            }
        }
    }
    if (haveData && outfile)
    {
        FILE *muc = xvgropen(outfile, title ? title : "", xaxis, yaxis, oenv);
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
        FILE *muc = xvgropen(outfile, title ? title : "", xaxis, yaxis, oenv);
        xvgr_symbolize(muc, eprnm, oenv);
        for(size_t i = 0; i < lsq.size(); i++)
        {
            print_lsq_set(muc, lsq[i]);
        }
        fclose(muc);
    }
}

static void write_q_histo(gmx::TextWriter                  *tw,
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
        print_stats(tw, "All Partial Charges", "e", 1.0, &gs->second, true,
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
                print_stats(tw, q.first.c_str(), "e", 1.0, &q.second, false, 
                            "CM5", model.c_str(), useOffset);
            }
        }
    }
    if (hh)
    {
        fclose(hh);
    }
}

void TrainForceFieldPrinter::addOptions(std::vector<t_pargs> *pargs)
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

void TrainForceFieldPrinter::addFileOptions(std::vector<t_filenm> *filenm)
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
        { efXVG, "-anisopol",  "anisopol_corr", ffOPTWR },
        { efXML, "-mpout",     "molprop_out",   ffOPTWR }
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

void TrainForceFieldPrinter::analysePolarisability(gmx::TextWriter     *tw,
                                                   const ForceField    *pd,
                                                   alexandria::ACTMol  *mol,
                                                   iMolSelect           ims,
                                                   const ForceComputer *forceComp)
{
    tw->writeStringFormatted("Polarizability:\n");
    for(auto qp = mol->qProps()->begin(); qp < mol->qProps()->end(); ++qp)
    {
        auto qcalc = qp->qPact();
        auto qelec = qp->qPqmConst();
        qcalc->calcPolarizability(pd, mol->topology(), forceComp);
        auto acalc = qcalc->polarizabilityTensor();
    
        print_polarizability(tw, mol->getMolname(), &qelec, qcalc, alpha_toler_, isopol_toler_);
        if (qelec.hasPolarizability())
        {
            auto aelec = qelec.polarizabilityTensor();
            lsq_isoPol_[ims][qType::Calc].add_point(qelec.isotropicPolarizability(),
                                               qcalc->isotropicPolarizability(),
                                               0, 0);
            lsq_anisoPol_[ims][qType::Calc].add_point(qelec.anisotropicPolarizability(),
                                                 qcalc->anisotropicPolarizability(),
                                                 0, 0);
            for (int mm = 0; mm < DIM; mm++)
            {
                lsq_alpha_[ims][qType::Calc].add_point(aelec[mm][mm], acalc[mm][mm], 0, 0);
            }
        }
    }
}

void TrainForceFieldPrinter::printAtoms(gmx::TextWriter              *tw,
                                        alexandria::ACTMol           *mol,
                                        const std::vector<gmx::RVec> &coords,
                                        const std::vector<gmx::RVec> &forces)
{
    std::map<qType, const std::vector<double> > qQM;
    std::vector<qType>                          typeQM = { 
        qType::CM5, qType::ESP, qType::Hirshfeld, qType::Mulliken,
        qType::RESP, qType::BCC
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
    tw->writeStringFormatted("Atom   Type            ACM");
    for(auto &qt : qQM)
    {
        if (!qQM.find(qt.first)->second.empty())
        {
            tw->writeStringFormatted("%10s", qTypeName(qt.first).c_str());
        }
    }
    tw->writeStringFormatted("        x        y   z(pm)          fx         fy fz(kJ/mol nm)     qtot\n");
    int      i       = 0;
    double   qtot    = 0;
    auto    &myatoms = mol->atomsConst();
    for (size_t j = i = 0; j < myatoms.size(); j++)
    {
        if (myatoms[j].pType() == ActParticle::Atom)
        {
            real qCalc = myatoms[j].charge();
            tw->writeStringFormatted("%-2d%3lu  %-5s  %12g",
                    myatoms[j].atomicNumber(),
                    j+1,
                    myatoms[j].ffType().c_str(),
                    qCalc);
            qtot += qCalc;
            for(auto &qt : qQM)
            {
                if (!qQM.find(qt.first)->second.empty())
                {
                    tw->writeStringFormatted("  %8.4f", qt.second[i]);
                }
            }
            tw->writeStringFormatted(" %8.3f %8.3f %8.3f %10.3f %10.3f %10.3f  %10g\n", 
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
            tw->writeStringFormatted("%-2d%3lu  %-5s  %12g",
                                     0,
                                     j+1,
                                     myatoms[j].ffType().c_str(),
                                     myatoms[j].charge());
            qtot += myatoms[j].charge();
            for(auto &qt : qQM)
            {
                if (!qQM.find(qt.first)->second.empty())
                {
                    tw->writeStringFormatted("          ");
                }
            }
            tw->writeStringFormatted(" %8.3f %8.3f %8.3f %10.3f %10.3f %10.3f  %10g\n", 
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
                      alexandria::ACTMol        *mol,
                      std::vector<gmx::RVec>    *coords,
                      double                     epot,
                      const AtomizationEnergy   &atomenergy,
                      const std::vector<double> &freq)
{
    real scale_factor = 1;
    real roomTemp     = 298.15;
    std::string eunit("kJ/mol");
    std::string sunit("J/mol K");
    jtree->addValueUnit("Epot", gmx_ftoa(epot), eunit);
    JsonTree jtc0("T=0");
    JsonTree jtcRT(gmx::formatString("T=%g", roomTemp));
    ThermoChemistry tc0(mol, *coords, atomenergy, freq, epot, 0.0, 1, scale_factor);
    ThermoChemistry tcRT(mol, *coords, atomenergy, freq, epot, roomTemp, 1, scale_factor);
    {
        jtc0.addValueUnit("Zero point energy", gmx_ftoa(tc0.ZPE()), eunit);
        jtcRT.addValueUnit("Zero point energy", gmx_ftoa(tcRT.ZPE()), eunit);
    }
    {
        JsonTree jtS0("Standard entropy");
        JsonTree jtSRT("Standard entropy");
        for(const auto &tcc : tccmap())
        {
            jtS0.addValueUnit(tcc.second, gmx_ftoa(tc0.S0(tcc.first)), sunit);
            jtSRT.addValueUnit(tcc.second, gmx_ftoa(tcRT.S0(tcc.first)), sunit);
        }
        jtc0.addObject(jtS0);
        jtcRT.addObject(jtSRT);
    }
    {
        JsonTree jtcv0("Heat capacity cV");
        JsonTree jtcvRT("Heat capacity cV");
        for(const auto &tcc : tccmap())
        {
            jtcv0.addValueUnit(tcc.second, gmx_ftoa(tc0.cv(tcc.first)), sunit);
            jtcvRT.addValueUnit(tcc.second, gmx_ftoa(tcRT.cv(tcc.first)), sunit);
        }
        jtc0.addObject(jtcv0);
        jtcRT.addObject(jtcvRT);
    }
    {
        JsonTree jtcE0("Internal energy");
        JsonTree jtcERT("Internal energy");
        for(const auto &tcc : tccmap())
        {
            jtcE0.addValueUnit(tcc.second, gmx_ftoa(tc0.Einternal(tcc.first)), eunit);
            jtcERT.addValueUnit(tcc.second, gmx_ftoa(tcRT.Einternal(tcc.first)), eunit);
        }
        jtc0.addObject(jtcE0);
        jtcRT.addObject(jtcERT);
    }
    {
        jtc0.addValueUnit("Delta H formation", gmx_ftoa(tc0.DHform()), eunit);
        jtcRT.addValueUnit("Delta H formation", gmx_ftoa(tcRT.DHform()), eunit);
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
                         bool                      debugNMA)
{
    std::vector<double> alex_freq, intensities;
    std::vector<std::string> output;
    molhandler.nma(pd, mol, forceComp, coords, &alex_freq, &intensities,
                   &output, debugNMA);
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
    // Energies
    std::map<InteractionType, double> energies;
    std::vector<gmx::RVec>    forces(coords->size());
    forceComp->compute(nullptr, pd, mol->topology(), coords, &forces, &energies);
    // Now time for thermochemistry output
    JsonTree tctree("Thermochemistry");
    JsonTree ajtc("Alexandria");
    auto     epot = energies[InteractionType::EPOT];
    addThermo(&ajtc, mol, coords, epot, atomenergy, alex_freq);
    tctree.addObject(ajtc);
    if (!ref_freq.empty())
    {
        JsonTree rjtc("Reference");
        addThermo(&rjtc, mol, coords, epot, atomenergy, ref_freq);
        tctree.addObject(rjtc);
    }
    jtree->addObject(tctree);
}

static std::string low_print_stats(gmx_stats  *stats,
                                   const char *label,
                                   const char *compound)
{
    real mse, mae, rmsd, R = 0;
    stats->get_mse_mae(&mse, &mae);
    stats->get_rmsd(&rmsd);
    if (stats->get_npoints() > 2)
    {
        stats->get_corr_coeff(&R);
    }
    return gmx::formatString("%18s RMSD %8.2f MSE %8.2f MAE %8.2f (kJ/mol) R %5.1f%% #points = %5zu %s",
                             label, rmsd, mse, mae, 100*R, stats->get_npoints(), compound);
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
        for(const auto itype : { InteractionType::ELECTROSTATICS, InteractionType::INDUCTION })
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

void TrainForceFieldPrinter::writeMolpropsEnergies(MsgHandler          *msghandler,
                                                   const char          *mpout,
                                                   const ForceField    *pd,
                                                   const ForceComputer *forceComp,
                                                   std::vector<ACTMol> *mols)
{
    std::vector<MolProp> mps;
    std::string version("ACT-");
    version += act_version();
    for (auto mol = mols->begin(); mol < mols->end(); ++mol)
    {
        std::vector<ACTEnergy>                                              energyMap;
        std::vector<std::vector<std::pair<double, double> > >               forceMap;
        std::vector<std::pair<double, std::map<InteractionType, double> > > energyComponentMap;
        ACTEnergyMapVector                                                  interactionEnergyMap;
        mol->forceEnergyMaps(msghandler, pd, forceComp, &forceMap,
                             &energyMap, &interactionEnergyMap, &energyComponentMap,
                             mol->hasMolPropObservable(MolPropObservable::INDUCTIONCORRECTION));
        // TODO: This is copied from actmol.cpp, store in one place instead?
        std::map<MolPropObservable, InteractionType> interE = {
            { MolPropObservable::INTERACTIONENERGY,   InteractionType::EPOT                 },
            { MolPropObservable::ELECTROSTATICS,      InteractionType::ELECTROSTATICS       },
            { MolPropObservable::DISPERSION,          InteractionType::DISPERSION           },
            { MolPropObservable::EXCHANGE,            InteractionType::EXCHANGE             },
            { MolPropObservable::VDWCORRECTION,       InteractionType::VDWCORRECTION        },
            { MolPropObservable::INDUCTIONCORRECTION, InteractionType::INDUCTIONCORRECTION  },
            { MolPropObservable::INDUCTION,           InteractionType::INDUCTION            },
            { MolPropObservable::CHARGETRANSFER,      InteractionType::CHARGETRANSFER       }
        };

        MolProp mpnew = *mol->molProp();

        size_t ii = 0;
        for (auto &exper : *mpnew.experiment())
        {
            exper.setReference("Spoel2025a");
            exper.setProgram(version);
            exper.setMethod("train_ff");
            exper.setBasisset(pd->filename());
            std::map<MolPropObservable, std::vector<GenericProperty *> > newprops;
            auto allprops = exper.properties();
            for(auto prop = allprops->begin(); prop != allprops->end(); prop++)
            {
                auto ie = interE.find(prop->first);
                if (interE.end() != ie)
                {
                    if (prop->second.size() != 1)
                    {
                        // Can only handle one energy at a time
                        fprintf(stderr, "More than one energy %s in experiment data structure for %s\n",
                                mpo_name(prop->first), mpnew.getMolname().c_str());
                    }
                    else
                    {
                        auto me    = static_cast<MolecularEnergy *>(prop->second[0]);
                        auto ae    = interactionEnergyMap[ii].find(ie->second);
                        if (ae != interactionEnergyMap[ii].end())
                        {
                            auto value = ae->second.eact();
                            me->Set(value, 0);
                            me->setInputUnit("kJ/mol");
                        }
                    }
                }
                else
                {
                    // Do nothing
                }
            }
            ii += 1;
        }
        mps.push_back(mpnew);
    }
    MolPropWrite(mpout, mps, false);
}

void TrainForceFieldPrinter::printEnergyForces(MsgHandler                          *msghandler,
                                               const ForceField                    *pd,
                                               const ForceComputer                 *forceComp,
                                               const std::map<eRMS, FittingTarget> &targets,
                                               const AtomizationEnergy             &atomenergy,
                                               alexandria::ACTMol                  *mol,
                                               iMolSelect                           ims,
                                               const gmx_output_env_t              *oenv,
                                               bool                                 printAll)
{
    std::vector<ACTEnergy>                                              energyMap;
    std::vector<std::vector<std::pair<double, double> > >               forceMap;
    std::vector<std::pair<double, std::map<InteractionType, double> > > energyComponentMap;
    ACTEnergyMapVector                                                  interactionEnergyMap;
    mol->forceEnergyMaps(msghandler, pd, forceComp, &forceMap,
                         &energyMap, &interactionEnergyMap,
                         &energyComponentMap,
                         mol->hasMolPropObservable(MolPropObservable::INDUCTIONCORRECTION));
    molEnergyMap_.insert({mol->getMolname(), energyMap});
    molInteractionEnergyMap_.insert({mol->getMolname(), interactionEnergyMap});
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
    if (targets.find(eRMS::EPOT) != targets.end())
    {
        for(const auto &ff : energyMap)
        {
            if (ff.haveQM() && ff.haveACT())
            {
                lsq_epot_[ims][qType::Calc].add_point(ff.eqm(), ff.eact(), 0, 0);
                myepot.add_point(ff.eqm(), ff.eact(), 0, 0);
            }
        }
    }
    std::map<InteractionType, gmx_stats> myeinter;
    std::map<eRMS, InteractionType> interE = {
        { eRMS::Interaction,    InteractionType::EPOT                },
        { eRMS::Electrostatics, InteractionType::ELECTROSTATICS      },
        { eRMS::Dispersion,     InteractionType::DISPERSION          },
        { eRMS::Exchange,       InteractionType::EXCHANGE            },
        { eRMS::Induction,      InteractionType::INDUCTION           },
        { eRMS::DeltaHF,        InteractionType::INDUCTIONCORRECTION },
        { eRMS::ExchInd,        InteractionType::EXCHIND             },
        { eRMS::AllElec,        InteractionType::ALLELEC             }
    };
    for(const auto &iem : interactionEnergyMap)
    {
        for(const auto &ie: iem)
        {
            auto val    = ie.first;
            auto result = std::find_if(interE.begin(), interE.end(),
                                       [val](const auto& mo) {return mo.second == val; });

            if (interE.end() == result)
            {
                continue;
            }
            auto tt = targets.find(result->first);
            if ((printAll || ( tt != targets.end() && tt->second.weight() > 0)))
            {
                if (ie.second.haveQM() && ie.second.haveACT())
                {
                    auto eqm  = ie.second.eqm();
                    auto eact = ie.second.eact();
                    lsq_einter_[ie.first][ims][qType::Calc].add_point(eqm, eact, 0, 0);
                    myeinter[ie.first].add_point(eqm, eact, 0, 0);
                }
            }
        }
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
            for(const auto &term : eam->second)
            {
                ttt += gmx::formatString(" %s: %g", interactionTypeToString(term.first).c_str(), 
                                         term.second);
            }
            ttt += gmx::formatString(" Diff: %g", eam->second[InteractionType::EPOT]-enerexp);
            msghandler->write(ttt);
            ccc++;
        }
        std::string myheader;
        std::string mysub;
        for(const auto &t : terms_)
        {
            myheader += gmx::formatString("%15s     ", interactionTypeToString(t).c_str());
            mysub    += gmx::formatString(" %9s %9s", "QM", "ACT");
        }
        msghandler->write(myheader);
        msghandler->write(mysub);
        
        std::sort(interactionEnergyMap.begin(), interactionEnergyMap.end(),
                  [](const ACTEnergyMap &a, const ACTEnergyMap &b)
        { 
            auto &aa = a.find(InteractionType::EPOT)->second;
            auto &bb = b.find(InteractionType::EPOT)->second;
            if (aa.haveQM() && bb.haveQM())
            {
                return (aa.eqm() < bb.eqm());
            }
            else
            {
                return false;
            }
        });
        auto expers = mol->experimentConst();
        for(const auto &iem : interactionEnergyMap)
        {
            std::string myline;
            for(const auto &t : terms_)
            {
                auto ttt = iem.find(t);
                
                if (iem.end() != ttt && ttt->second.haveQM())
                {
                    myline += gmx::formatString(" %9.2f", ttt->second.eqm());
                }
                else
                {
                    myline += ("         x");
                }
                if (iem.end() != ttt && ttt->second.haveACT())
                {
                    myline += gmx::formatString(" %9.2f", ttt->second.eact());
                }
                else
                {
                    myline += ("         x");
                }
            }
            auto ttt   = iem.find(InteractionType::EPOT);
            auto exper = std::find_if(expers.begin(), expers.end(), 
                                      [&ttt](const Experiment& x)
                                      { return x.id() == ttt->second.id(); });
            if (expers.end() != exper)
            {
                myline += " " + exper->getDatafile();
            }
            msghandler->write(myline);
        }
        std::sort(forceMap.begin(), forceMap.end());
        for(auto ff = forceMap.begin(); ff < forceMap.end(); ++ff)
        {
            for(const auto &fxyz : *ff)
            {
                if (fxyz.first != 0 || fxyz.second != 0)
                {
                    msghandler->write(gmx::formatString("Force ref: %8.2f  act: %8.2f", fxyz.first, fxyz.second));
                    ccc++;
                }
            }
        }
    }
    // RMS energy
    if (energyMap.size() > 0)
    {
        msghandler->write(low_print_stats(&myepot, "Energy", mol->getMolname().c_str()));
    }
    for(const auto &tt : terms_)
    {
        auto mystat = &myeinter[tt];
        if (mystat->get_npoints() > 0)
        {
            msghandler->write(low_print_stats(mystat, interactionTypeToString(tt).c_str(),
                                              mol->getMolname().c_str()));
        }
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
                    lsq_rmsf_[ims].add_point(ff.first, ff.second, 0, 0);
                    myforce.add_point(ff.first, ff.second, 0, 0);
                }
            }
        }
        msghandler->write(low_print_stats(&myforce, "Force", mol->getMolname().c_str()));
    }

    double deltaE0 = 0;
    if (mol->energy(MolPropObservable::DELTAE0, &deltaE0))
    {
        // Energy
        msghandler->write(gmx::formatString("Energy terms (kJ/mol)"));
        std::map<InteractionType, double> eBefore;
        std::vector<gmx::RVec> coords = mol->xOriginal();
        std::vector<gmx::RVec> forces(coords.size());
        forceComp->compute(msghandler, pd, mol->topology(), &coords, &forces, &eBefore);

        msghandler->write(gmx::formatString("   %-20s  %10.3f  Difference: %10.3f",
                                            "Reference EPOT", deltaE0, 
                                            eBefore[InteractionType::EPOT]-deltaE0));

        if (mol->jobType() == JobType::OPT && calcFrequencies_)
        {
            // Now get the minimized structure RMSD and Energy
            SimulationConfigHandler simConfig;
            std::vector<gmx::RVec>  xmin   = coords;
            std::map<InteractionType, double> eAfter;
            molHandler_.minimizeCoordinates(msghandler, pd, mol, forceComp, simConfig, 
                                            &xmin, &eAfter, {});
            double rmsd = molHandler_.coordinateRmsd(mol, coords, &xmin);

            if (rmsd > 0.1) // nm
            {
                auto pdb   = gmx::formatString("inds/%s-original-minimized.pdb", mol->getMolname().c_str());
                auto title = gmx::formatString("%s RMSD %g Angstrom", mol->getMolname().c_str(), 10*rmsd);
                writeCoordinates(mol->atomsConst(), pdb, title, coords, xmin);
            }
            msghandler->write(gmx::formatString("   %-20s  %10s  %10s  %10s minimization",
                                                "Term", "Before", "After", "Difference"));
            for(auto &ep : eBefore)
            {
                auto eb = ep.second;
                auto ea = eAfter[ep.first];
                msghandler->write(gmx::formatString("   %-20s  %10.3f  %10.3f  %10.3f",
                                                    interactionTypeToString(ep.first).c_str(),
                                                    eb, ea, ea-eb));
            }
            msghandler->write(gmx::formatString("Coordinate RMSD after minimization %10g pm", 1000*rmsd));

            // Do normal-mode analysis etc.
            JsonTree jtree("FrequencyAnalysis");
            doFrequencyAnalysis(pd, mol, molHandler_, forceComp, &coords,
                                atomenergy, &lsq_freq_, &jtree,
                                nullptr, 24, nullptr, false);
            int indent = 0;
            msghandler->write(jtree.writeString(false, &indent));
        }
        else
        {
            for(auto &ep : eBefore)
            {
                msghandler->write(gmx::formatString("   %-20s  %10.3f",
                                                    interactionTypeToString(ep.first).c_str(),
                                                    ep.second));
            }
        }
    }
}

static void dump_xyz(const std::string &label,
                     const ACTMol      *mol,
                     const ACTEnergy   &actener)
{
    auto expconst = mol->experimentConst();
    if (static_cast<size_t>(actener.id()) < expconst.size())
    {
        auto filename = gmx::formatString("%s-%d.xyz", label.c_str(), actener.id());
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
    else
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find experiment id %d for %s (there are %zu only)",
                                                       actener.id(),
                                                       mol->getMolname().c_str(),
                                                       expconst.size()).c_str()));
    }
}

static int printLow(gmx::TextWriter   *tw,
                    const ACTMol      *mol,
                    const ACTEnergy   &ener,
                    const std::string &label,
                    double             sigma,
                    double             dumpOutliers)
{
    int    noutlier = 0;
    double epotMax  = 1.5*sigma;
    if (ener.haveQM() && ener.haveACT() && mol)
    {
        double deltaE = std::abs(ener.eqm()-ener.eact());
        if (deltaE > epotMax)
        {
            tw->writeStringFormatted("%-40s  %12g  %12g  %12g\n", mol->getMolname().c_str(),
                    ener.eqm(), ener.eact(), ener.eact()-ener.eqm());
            noutlier++;
        }
        if (dumpOutliers >= 0 && deltaE >= dumpOutliers*sigma)
        {
            std::string myfnm = gmx::formatString("%s-%s.xyz", 
                                                  mol->getMolname().c_str(),
                                                  label.c_str());
            dump_xyz(myfnm, mol, ener);
        }
    }
    return noutlier;
}

void TrainForceFieldPrinter::printOutliers(gmx::TextWriter                       *tw,
                                           iMolSelect                             ims,
                                           double                                 sigma,
                                           bool                                   bIntermolecular,
                                           InteractionType                        itype,
                                           const std::vector<alexandria::ACTMol> *actmol)
{
    double epotMax = 1.5*sigma;
    auto   label   = gmx::formatString("%s-%s",
                                       interactionTypeToString(itype).c_str(), iMolSelectName(ims));
    tw->writeStringFormatted("\nOverview of %s outliers for %s (Diff > %.3f)\n",
            label.c_str(), qTypeName(qType::Calc).c_str(), epotMax);
    tw->writeStringFormatted("----------------------------------\n");
    tw->writeStringFormatted("%-40s  %12s  %12s  %12s\n", "Name",
            "Reference", qTypeName(qType::Calc).c_str(), "ACT-Ref.");
    int noutlier = 0;
    if (bIntermolecular)
    {
        for (const auto &miem : molInteractionEnergyMap_)
        {
            std::string toFind(miem.first);
            auto actmolptr = std::find_if(actmol->begin(), actmol->end(),
                                          [&toFind](const ACTMol &x) { return x.getMolname() == toFind;});
            if (actmolptr->datasetType() == ims)
            {
                for(const auto &emm : miem.second)
                {
                    auto emmptr = emm.find(itype);
                    if (emm.end() == emmptr)
                    {
                        continue;
                    }
                    noutlier += printLow(tw, &(*actmolptr), emmptr->second, label, sigma, dumpOutliers_);
                }
            }
        }
    }
    else
    {
        for (const auto &miem : molEnergyMap_)
        {
            std::string toFind(miem.first);
            auto actmolptr = std::find_if(actmol->begin(), actmol->end(),
                                          [&toFind](const ACTMol &x) { return x.getMolname() == toFind;});
            if (actmolptr->datasetType() == ims)
            {
                for(const auto &emm : miem.second)
                {
                    noutlier += printLow(tw, &(*actmolptr), emm, label, sigma, dumpOutliers_);
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

void TrainForceFieldPrinter::print(MsgHandler                  *msghandler,
                                   StaticIndividualInfo        *sii,
                                   std::vector<ACTMol>         *actmol,
                                   const gmx_output_env_t      *oenv,
                                   const std::vector<t_filenm> &filenm,
                                   bool                         printAll)
{
    const auto pd = sii->forcefield();
    int        n  = 0;
    std::map<iMolSelect, qtStats>                               lsq_charge, lsq_esp;
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
        lsq_epot_.insert({ ims.first, std::move(qepot) });
        
        gmx_stats galpha;
        lsq_alpha_.insert(std::pair<iMolSelect, qtStats>(ims.first, {{ qType::Calc, galpha }}));
        gmx_stats giso;
        lsq_isoPol_.insert(std::pair<iMolSelect, qtStats>(ims.first, {{ qType::Calc, giso }}));
        gmx_stats ganiso;
        lsq_anisoPol_.insert(std::pair<iMolSelect, qtStats>(ims.first, {{ qType::Calc, ganiso }}));
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
    for(const auto &tt : terms_)
    {
        std::map<iMolSelect, qtStats> qinter;
        for(auto &ims : iMolSelectNames())
        {
            for (auto &i : qTypes())
            {
                gmx_stats ginter;
                qinter[ims.first].insert({ i.first, std::move(ginter) });
            }
        }
        lsq_einter_[tt] = std::move(qinter);
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

    atomenergy.read(msghandler);
    auto tw = msghandler->tw();
    for (auto mol = actmol->begin(); mol < actmol->end(); ++mol)
    {
        auto topology = mol->topology();
        if (tw && mol->support() != eSupport::No)
        {
            auto ims = mol->datasetType();
            tw->writeStringFormatted("\nMolecule %d: Name: %s, Qtot: %d, Multiplicity: %d, MolWt: %g SymmetryNumber: %d Dataset: %s\n", n+1,
                    mol->getMolname().c_str(),
                    mol->totalCharge(),
                    mol->totalMultiplicity(),
                    mol->totalMass(),
                    mol->symmetryNumber(),
                    iMolSelectName(ims));

            gmx::RVec vzero = { 0, 0, 0 };
            std::vector<gmx::RVec> forces(mol->atomsConst().size(), vzero);

            // Now compute all the ESP RMSDs and multipoles and print it.
            tw->writeStringFormatted("Electrostatic properties.\n");
            for (auto &i : qTypes())
            {
                for(auto qp = mol->qProps()->begin(); qp < mol->qProps()->end(); ++qp)
                {
                    // For multipoles, further down.
                    auto qi    = i.first;
                    auto qresp = qp->qgenResp();
                    if (qresp->nEsp() > 0 && qp->qPactConst().qtype() == qi)
                    {
                        real rms, rrms, cosesp, mae, mse;
                        // Fetch coordinates and optimize shells if polarizable
                        auto myx = qresp->coords();
                        std::map<InteractionType, double> energies;
                        forceComp->compute(msghandler, pd, topology, &myx,
                                           &forces, &energies);
                        qresp->updateAtomCoords(myx);
                        qresp->updateAtomCharges(mol->atomsConst());
                        qresp->calcPot(msghandler, 1.0);
                        rms = qresp->getStatistics(msghandler, &rrms, &cosesp, &mae, &mse);
                        rms = convertToGromacs(rms, "Hartree/e");
                        std::string warning;
                        if (rms > esp_toler_ || cosesp < 0.5)
                        {
                            warning.assign(" EEE");
                        }
                        tw->writeStringFormatted("ESP rms: %8.3f (kJ/mol e) rrms: %8.3f CosAngle: %6.3f - %s%s\n",
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
            auto atoms = topology->atoms();
            auto lll   = lsqt.find(ims);
            for(size_t ai = 0; ai < atoms.size(); ai++)
            {
                if (atoms[ai].pType() == ActParticle::Atom)
                {
                    auto llF = lll->second.find(atoms[ai].ffType());
                    if (lll->second.end() != llF)
                    {
                        auto qa = atoms[ai].charge();
                        if (ai < atoms.size() -1 && 
                            atoms[ai+1].pType() == ActParticle::Shell)
                        {
                            qa += atoms[ai+1].charge();
                        }
                        llF->second.add_point(1, qa, 0, 0);
                    }
                }
            }
            // Multipoles
            analyse_multipoles(msghandler, mol, multi_toler, pd, forceComp);
            
            // Polarizability
            if (bPolar)
            {
                analysePolarisability(tw, pd, &(*mol), ims, forceComp);
            }

            // Atomic charges
            std::vector<gmx::RVec> coords = mol->xOriginal();
            {
                std::map<InteractionType, double> energies;
                forceComp->compute(msghandler, pd, topology, &coords,
                                   &forces, &energies);
            }
            printAtoms(tw, &(*mol), coords, forces);
            // Energies
            if (sii->haveFittingTargetSelection(ims))
            {
                std::vector<std::string> tcout;
                printEnergyForces(msghandler, pd, forceComp, sii->fittingTargetsConst(ims),
                                  atomenergy, &(*mol), ims, oenv, printAll);

                for(const auto &tout : tcout)
                {
                    tw->writeStringFormatted("%s\n", tout.c_str());
                }
                tw->writeStringFormatted("\n");
            }
            n++;
        }
    }

    for(auto &ims : iMolSelectNames())
    {
        bool header = true;
        
        tw->writeStringFormatted("\n*** Results for %s data set ***\n", iMolSelectName(ims.first));
        for (auto &i : qTypes())
        {
            auto qt = i.first;
            if (qt == qType::Elec)
            {
                continue;
            }
            std::string name = gmx::formatString("%s-%s", qTypeName(qt).c_str(), ims.second);
            if (lsq_epot_[ims.first][qt].get_npoints() > 0)
            {
                print_stats(tw, "Potential energy", "kJ/mol", 1.0, &lsq_epot_[ims.first][qt],
                            header, "QM/DFT", name.c_str(), useOffset_);
                header = false;
            }
            for(const auto &tt : terms_)
            {
                if (lsq_einter_[tt][ims.first][qt].get_npoints() > 0)
                {
                    print_stats(tw, interactionTypeToString(tt).c_str(),
                                "kJ/mol", 1.0, &lsq_einter_[tt][ims.first][qt],
                                header, "SAPT", name.c_str(), useOffset_);
                    header = false;
                }
            }
            if (lsq_rmsf_[ims.first].get_npoints() > 0 && qt == qType::Calc)
            {
                print_stats(tw, "RMS Force", "kJ/mol nm", 1.0, &lsq_rmsf_[ims.first],
                            header, "QM/DFT", name.c_str(), useOffset_);
                header = false;
            }
            if (lsq_freq_.get_npoints() > 0 && qt == qType::Calc)
            {
                print_stats(tw, "Frequencies", mpo_unit2(MolPropObservable::FREQUENCY),
                            1.0, &lsq_freq_, header, "QM/DFT", name.c_str(), useOffset_);
                header = false;
            }
            print_stats(tw, "ESP", "kJ/mol e", 1.0, &lsq_esp[ims.first][qt],
                        header, "Electronic", name.c_str(), useOffset_);
            header = false;
            for(auto &mpo : mpoMultiPoles)
            {
                print_stats(tw, mpo_name(mpo), mpo_unit2(mpo), 1.0,
                            &lsq_multi[mpo][ims.first][qt],   header, "Electronic", name.c_str(), useOffset_);
            }
            if (bPolar && qt == qType::Calc)
            {
                std::string polunit("Angstrom3");
                auto polfactor = convertFromGromacs(1.0, polunit);
                print_stats(tw, "Polariz. components", polunit.c_str(), polfactor,
                            &lsq_alpha_[ims.first][qType::Calc],    header, "Electronic", name.c_str(), useOffset_);
                print_stats(tw, "Isotropic Polariz.", polunit.c_str(), polfactor,
                            &lsq_isoPol_[ims.first][qType::Calc],   header, "Electronic", name.c_str(), useOffset_);
                print_stats(tw, "Anisotropic Polariz.", polunit.c_str(), polfactor,
                            &lsq_anisoPol_[ims.first][qType::Calc], header, "Electronic", name.c_str(), useOffset_);
            }
        }
    }
    write_q_histo(tw, opt2fn_null("-qhisto", filenm.size(), filenm.data()),
                  &lsqt[iMolSelect::Train], oenv,
                  &(lsq_charge[iMolSelect::Train]), useOffset_);

    const char *alex  = "Alexandria";
    for(auto &mpo : mpoMultiPoles)
    {
        std::string cmdFlag = gmx::formatString("-%scorr", mpo_name(mpo));
        std::string title   = gmx::formatString("%s components (%s)", mpo_name(mpo),
                                                mpo_unit2(mpo));
        print_corr(opt2fn_null(cmdFlag.c_str(), filenm.size(), filenm.data()),
                   nullptr, title.c_str(), alex, lsq_multi[mpo], oenv);
    }
    print_corr(opt2fn_null("-epotcorr", filenm.size(), filenm.data()),
               nullptr, "Potential energy (kJ/mol)", alex,
               lsq_epot_, oenv);
    for(const auto &tt : terms_)
    {
        std::string fnm   = gmx::formatString("%s.xvg", interactionTypeToString(tt).c_str());
        std::string title = gmx::formatString("%s (kJ/mol)", interactionTypeToString(tt).c_str());
        print_corr(fnm.c_str(), nullptr, title.c_str(), alex, lsq_einter_[tt], oenv);
    }
    print_corr(opt2fn_null("-forcecorr", filenm.size(), filenm.data()),
               nullptr, "Forces (kJ/mol nm)", alex,
               lsq_rmsf, oenv);
    print_corr(opt2fn_null("-freqcorr", filenm.size(), filenm.data()),
               nullptr, "Frequencies (cm^-1)", alex,
               lsq_freq, oenv);
    print_corr(opt2fn_null("-espcorr", filenm.size(), filenm.data()),
               nullptr, "Electrostatic Potential (kJ/mol e)", alex,
               lsq_esp, oenv);
    print_corr(opt2fn_null("-qcorr", filenm.size(), filenm.data()),
               nullptr, "Atomic Partial Charge (e)", "a.u.",
               lsq_charge, oenv);
    
    if (bPolar)
    {
        print_corr(opt2fn_null("-alphacorr", filenm.size(), filenm.data()),
                   nullptr, "Pricipal Components of Polarizability Tensor (A\\S3\\N)", alex, lsq_alpha_, oenv);
        print_corr(opt2fn_null("-isopol", filenm.size(), filenm.data()),
                   nullptr, "Isotropic Polarizability (A\\S3\\N)", alex, lsq_isoPol_, oenv);
        print_corr(opt2fn_null("-anisopol", filenm.size(), filenm.data()),
                   nullptr, "Anisotropic Polarizability (A\\S3\\N)", alex, lsq_anisoPol_, oenv);
    }
    // List outliers based on the deviation in the Electrostatic Potential
    real espAver;
    if (eStats::OK == lsq_esp[iMolSelect::Train][qType::Calc].get_rmsd(&espAver))
    {
        int nout = 0;
        double espMax = 1.5*espAver;
        tw->writeStringFormatted("\nOverview of ESP outliers for %s (RMSD > %.3f)\n",
                qTypeName(qType::Calc).c_str(), espMax);
        tw->writeStringFormatted("----------------------------------\n");
        tw->writeStringFormatted("%-40s  %12s  %12s\n", "Name",
                qTypeName(qType::Calc).c_str(), qTypeName(qType::ESP).c_str());
        for (auto mol = actmol->begin(); mol < actmol->end(); ++mol)
        {
            for(auto qp = mol->qProps()->begin(); qp < mol->qProps()->end(); ++qp)
            {
                real rms, rrms, cosesp, mae, mse;
                auto qresp = qp->qgenResp();
                rms        = convertToGromacs(qresp->getStatistics(msghandler, &rrms, &cosesp, &mae, &mse), "Hartree/e");
                if ((mol->support() != eSupport::No) && (rms > espMax))
                {
                    tw->writeStringFormatted("%-40s  %12.3f", mol->getMolname().c_str(), rms);
                    tw->writeStringFormatted("  %s\n", iMolSelectName(mol->datasetType()));
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
    if (
        eStats::OK == lsq_epot_[iMolSelect::Train][qType::Calc].get_rmsd(&epotRmsd))
    {
        for(const auto &ims : { iMolSelect::Train, iMolSelect::Test })
        {
            if (lsq_epot_[ims][qType::Calc].get_npoints() > 0)
            {
                // List outliers based on the deviation in the Potential energy ...
                printOutliers(tw, ims, epotRmsd, false, InteractionType::EPOT, actmol);
            }
        }
    }
    for(const auto &tt : terms_)
    {
        real einterRmsd = 0;
        if (eStats::OK == lsq_einter_[tt][iMolSelect::Train][qType::Calc].get_rmsd(&einterRmsd))
        {
            for(const auto &ims : { iMolSelect::Train, iMolSelect::Test })
            {
                if (lsq_einter_[tt][ims][qType::Calc].get_npoints() > 0)
                {
                    // ... and the interaction energies.
                    printOutliers(tw, ims, einterRmsd, true, tt, actmol);
                }
            }
        }
    }

    const char *mpout = opt2fn_null("-mpout", filenm.size(), filenm.data());
    if (mpout)
    {
        writeMolpropsEnergies(msghandler, mpout, pd, forceComp, actmol);
    }
}

void print_header(gmx::TextWriter             *tw,
                  const std::vector<t_pargs>  &pargs,
                  const std::vector<t_filenm> &filenms)
{
    if (!tw)
    {
        return;
    }
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
        default: // Output only, non-training related code.
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
            tw->writeStringFormatted("%s%s\n", left.c_str(),
                                     wrapper.wrapToString(toWrap).c_str());
        }
    }
    if (filenms.size() > 0)
    {
        tw->writeStringFormatted("\nFiles used:\n");
    }
    std::map<unsigned int, const char *> fmap = {
        { ffREAD,  "R"   }, 
        { ffWRITE, "W"   },
        { ffSET,   "set" },
        { ffOPT,   "optional" }
    };
    for(const auto &f : filenms)
    {
        if ((f.flag & ffOPT) == ffOPT && (f.flag & ffSET) == 0)
        {
            continue;
        }
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
        tw->writeStringFormatted("Option: %s, Description: %s,", f.opt, ftp2desc(f.ftp));
        if (!flag.empty())
        {
            tw->writeStringFormatted(" Properties: %s, Filename(s)", flag.c_str());
        }
        if ((f.flag & ffMULT) == ffMULT)
        {
            for(auto &fnm : f.filenames)
            {
                tw->writeStringFormatted(" %s", fnm.c_str());
            }
        }
        else
        {
            if (nullptr != f.fn)
            {
                tw->writeStringFormatted(" %s", f.fn);
            }
        }
        tw->writeStringFormatted("\n");
    }
    tw->writeStringFormatted("\n");
}

} // namespace alexandria
