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

#include "tuning_utility.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/coolstuff.h"

#include "molprop_util.h"
#include "mymol.h"
#include "qtype.h"
#include "units.h"

void doAddOptions(std::vector<t_pargs> *pargs, size_t npa, t_pargs pa[])
{
    for (size_t i = 0; i < npa; i++)
    {
        pargs->push_back(pa[i]);
    }
}

namespace alexandria
{

using qtStats = std::map<qType, gmx_stats>;

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
    if (useOffset)
    {
        if (bHeader)
        {
            fprintf(fp, "Fitting data to y = ax + b, where x = %s\n", xaxis);
            fprintf(fp, "%-26s %6s %13s %13s %7s %8s %8s %8s %10s\n",
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
            fprintf(fp, "%-26s %6d %6.3f(%5.3f) %6.3f(%5.3f) %7.2f %8.4f %8.4f %8.4f %10s\n",
                    prop, n, a, da, b, db, Rfit*100, rmsd, mse, mae, yaxis);
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
            fprintf(fp, "%-26s %6s %13s %7s %8s %8s %8s %10s\n",
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
            fprintf(fp, "%-26s %6d %6.3f(%5.3f) %7.2f %8.4f %8.4f %8.4f %10s\n",
                    prop, n, a, da, Rfit*100, rmsd, mse, mae, yaxis);
        }
        else
        {
            fprintf(fp, "Statistics problem for %s: %s\n", prop, gmx_stats_message(ok));
        }
    }
}

static void print_lsq_sets(FILE *fp, const std::vector<gmx_stats> &lsq)
{
    if (0 == lsq.size())
    {
        return;
    }
    std::vector<std::vector<double> > x, y;
    x.resize(lsq.size());
    y.resize(lsq.size());
    for (size_t s = 0; s < lsq.size(); s++)
    {
        auto xx = lsq[s].getX();
        auto yy = lsq[s].getY();
        for (size_t k=0; k < xx.size(); k++)
        {
                x[s].push_back(xx[k]);
                y[s].push_back(yy[k]);
        }
    }
    
    fprintf(fp, "@type xy\n");
    for(size_t i = 0; i < lsq.size(); i++)
    {
        fprintf(fp, "%10g", x[0][i]);
        for(size_t j = 0; j < lsq.size(); j++)
        {
            fprintf(fp, "  %10g", y[j][i]);
        }
        fprintf(fp, "\n");
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
                                 const std::vector<alexandria::MyMol>::iterator &mol,
                                 const std::string &calc_name,
                                 real               alpha_toler,
                                 real               isopol_toler)
{
    tensor dalpha;
    real   delta      = 0;
    real   diso_pol   = 0;
    real   daniso_pol = 0;

    if (!calc_name.empty())
    {
        if (calc_name == qTypeName(qType::Calc))
        {
            m_sub(mol->alpha_elec(), mol->alpha_calc(), dalpha);
            delta = sqrt(gmx::square(dalpha[XX][XX])+gmx::square(dalpha[XX][YY])+gmx::square(dalpha[XX][ZZ])+
                         gmx::square(dalpha[YY][YY])+gmx::square(dalpha[YY][ZZ]));
            diso_pol   = std::abs(mol->PolarizabilityDeviation());
            daniso_pol = std::abs(mol->AnisoPolarizabilityDeviation());
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name.c_str(),
                    mol->alpha_calc()[XX][XX], mol->alpha_calc()[XX][YY], mol->alpha_calc()[XX][ZZ],
                    dalpha[XX][XX], dalpha[XX][YY], dalpha[XX][ZZ], delta, (delta > alpha_toler) ? "ALPHA" : "",
                    "", mol->alpha_calc()[YY][YY], mol->alpha_calc()[YY][ZZ],
                    "", dalpha[YY][YY], dalpha[YY][ZZ],
                    "", "", mol->alpha_calc()[ZZ][ZZ],
                    "", "", dalpha[ZZ][ZZ]);
            fprintf(fp,
                    "Isotropic polarizability:  %s Electronic: %6.2f  Calculated: %6.2f  Delta: %6.2f %s\n\n",
                    mol->getMolname().c_str(), 
                    mol->ElectronicPolarizability(),
                    mol->CalculatedPolarizability(), 
                    diso_pol, (diso_pol > isopol_toler) ? "ISO" : "");
            fprintf(fp,
                    "Anisotropic polarizability:  %s Electronic: %6.2f  Calculated: %6.2f  Delta: %6.2f %s\n\n",
                    mol->getMolname().c_str(), 
                    mol->ElectronicAnisoPolarizability(),
                    mol->CalculatedAnisoPolarizability(), 
                    daniso_pol, (daniso_pol > isopol_toler) ? "ANISO" : "");

        }
        else
        {
            fprintf(fp, "Polarizability analysis\n");
            fprintf(fp,
                    "Electronic   (%6.2f %6.2f %6.2f)\n"
                    "             (%6s %6.2f %6.2f)\n"
                    "             (%6s %6s %6.2f)\n",
                    mol->alpha_elec()[XX][XX], mol->alpha_elec()[XX][YY], mol->alpha_elec()[XX][ZZ],
                    "", mol->alpha_elec()[YY][YY], mol->alpha_elec()[YY][ZZ],
                    "", "", mol->alpha_elec()[ZZ][ZZ]);
        }
    }
}

static void analyse_quadrapole(FILE                                           *fp,
                               const std::vector<alexandria::MyMol>::iterator &mol,
                               real                                            q_toler,
                               std::map<qType, gmx_stats>                      lsq_quad,
                               bool                                            bFullTensor)
{
    auto qelec = mol->qTypeProps(qType::Elec);
    fprintf(fp, "Quadrupole analysis (5 independent components only)\n");
    fprintf(fp,
            "Electronic   (%6.2f %6.2f %6.2f)\n"
            "             (%6s %6.2f %6.2f)\n"
            "             (%6s %6s %6.2f)\n",
            qelec->quad()[XX][XX], qelec->quad()[XX][YY], qelec->quad()[XX][ZZ],
            "", qelec->quad()[YY][YY], qelec->quad()[YY][ZZ],
            "", "", qelec->quad()[ZZ][ZZ]);
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
            tensor        dQ;
            real          delta = 0;
            
            m_sub(qelec->quad(), qcalc->quad(), dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-10s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "           (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "           (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    qTypeName(qt).c_str(),
                    qcalc->quad()[XX][XX], qcalc->quad()[XX][YY], qcalc->quad()[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "QUAD" : "",
                    "", qcalc->quad()[YY][YY], qcalc->quad()[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", qcalc->quad()[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
            for (int mm = 0; mm < DIM; mm++)
            {
                for (int nn = 0; nn < DIM; nn++)
                {
                    if (bFullTensor || (mm == nn))
                    {
                        lsq_quad[qt].add_point(qelec->quad()[mm][nn], qcalc->quad()[mm][nn], 0, 0);
                    }
                }
            }
        }
    }
}

static void print_dipole(FILE       *fp,
                         qType       qt,
                         const rvec  muelec,
                         const rvec  muqt,
                         real        toler)
{
    rvec dmu;
    real ndmu, cosa;
    char ebuf[32];
    
    rvec_sub(muelec, muqt, dmu);
    ndmu = norm(dmu);
    cosa = cos_angle(muelec, muqt);
    if (ndmu > toler || fabs(cosa) < 0.1)
    {
        sprintf(ebuf, "DIP");
    }
    else
    {
        ebuf[0] = '\0';
    }
    fprintf(fp, "%-10s (%6.3f,%6.3f,%6.3f) |Mu| = %6.3f",
            qTypeName(qt).c_str(), muqt[XX], muqt[YY], muqt[ZZ], norm(muqt));

    if (qt != qType::Elec)
    {
        fprintf(fp, " |Dev|: %6.3f CosAngle %6.3f %s", ndmu, cosa, ebuf);
    }
    fprintf(fp, "\n");
}

static void print_corr(const char                         *outfile,
                       const char                         *title,
                       const char                         *xaxis,
                       const char                         *yaxis, 
                       std::map<iMolSelect, std::map<qType, gmx_stats> > &stats,
                       const gmx_output_env_t             *oenv)
{
    std::vector<std::string> eprnm;
    std::vector<gmx_stats>   lsq;
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
                }
            }
        }
    }
    FILE *muc = xvgropen(outfile, title, xaxis, yaxis, oenv);
    xvgr_symbolize(muc, eprnm, oenv);
    print_lsq_sets(muc, lsq);
    /*
    for (auto &ims : iMolSelectNames())
    {
        auto qs = stats.find(ims.first);
        if (qs != stats.end())
        {
            for (auto &i : qs->second)
            {
                int N;
                if (eStats::OK == gmx_stats_get_npoints(i.second, &N) && N > 0)
                {
                    print_lsq_set(muc, i.second);
                }
            }
        }
        }*/
    fclose(muc);
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
    FILE *hh = xvgropen(qhisto, "Histogram for charges", "q (e)", "a.u.", oenv);
    xvgrLegend(hh, types, oenv);
    
    auto model = qTypeName(qType::Calc);
    auto gs    = lsq_charge->find(qType::Calc);
    if (gs != lsq_charge->end())
    {
        print_stats(fplog, "All Partial Charges (e)", &gs->second, true,
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
                fprintf(hh, "@type xy\n");
                for (int i = 0; i < nbins; i++)
                {
                    fprintf(hh, "%10g  %10g\n", x[i], y[i]);
                }
                fprintf(hh, "&\n");
                print_stats(fplog, q.name().c_str(),  &q.lsq_, false,  "CM5", model.c_str(), useOffset);
            }
        }
    }
    fclose(hh);
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
        { "-alpha_toler", FALSE, etREAL, {&alpha_toler_},
          "Tolerance (A^3) for marking diagonal elements of the polarizability tensor as an outlier in the log file" },
        { "-isopol_toler", FALSE, etREAL, {&isopol_toler_},
          "Tolerance (A^3) for marking isotropic polarizability as an outlier in the log file" },
        { "-use_offset", FALSE, etBOOL,{&useOffset_},
          "Fit regression analysis of results to y = ax+b instead of y = ax" }
    };
    doAddOptions(pargs, sizeof(pa)/sizeof(pa[0]), pa);
}

void TuneForceFieldPrinter::addFileOptions(std::vector<t_filenm> *filenm)
{
    t_filenm fnm[] =
        {
            { efXVG, "-qhisto",    "q_histo",       ffWRITE },
            { efXVG, "-dipcorr",   "dip_corr",      ffWRITE },
            { efXVG, "-mucorr",    "mu_corr",       ffWRITE },
            { efXVG, "-epotcorr",  "epot_corr",     ffWRITE },
            { efXVG, "-thetacorr", "theta_corr",    ffWRITE },
            { efXVG, "-espcorr",   "esp_corr",      ffWRITE },
            { efXVG, "-alphacorr", "alpha_corr",    ffWRITE },
            { efXVG, "-qcorr",     "q_corr",        ffWRITE },
            { efXVG, "-isopol",    "isopol_corr",   ffWRITE },
            { efXVG, "-anisopol",  "anisopol_corr", ffWRITE }
        };
    int nfm = sizeof(fnm)/sizeof(fnm[0]);
    for(int i = 0; i < nfm; i++)
    {
        filenm->push_back(fnm[i]);
    }
}

void TuneForceFieldPrinter::print(FILE                           *fp,
                                  std::vector<alexandria::MyMol> *mymol,
                                  const Poldata                  *pd,
                                  const gmx::MDLogger            &fplog,
                                  const char                     *lot,
                                  int                             qcycle,
                                  real                            qtol,
                                  const gmx_output_env_t         *oenv,
                                  bool                            bfullTensor,
                                  t_commrec                      *cr,
                                  real                            efield,
                                  const std::vector<t_filenm>    &filenm)
{
    int  n = 0;

    std::map<iMolSelect, qtStats>                   lsq_mu, lsq_dip, lsq_quad, lsq_esp, 
        lsq_alpha, lsq_isoPol, lsq_anisoPol, lsq_charge, lsq_epot;
    std::map<iMolSelect, std::vector<ZetaTypeLsq> > lsqt;

    for(auto &ims : iMolSelectNames())
    {
        qtStats qesp, qquad, qdip, qmu, qepot;
        for (auto &i : qTypes())
        {
            //TODO Add checks for existence
            gmx_stats gesp;
            qesp.insert(std::pair<qType, gmx_stats>(i.first, std::move(gesp)));
            gmx_stats gquad;
            qquad.insert(std::pair<qType, gmx_stats>(i.first, std::move(gquad)));
            gmx_stats gdip;
            qdip.insert(std::pair<qType, gmx_stats>(i.first, std::move(gdip)));
            gmx_stats gmu;
            qmu.insert(std::pair<qType, gmx_stats>(i.first, std::move(gmu)));
            gmx_stats gepot;
            qepot.insert(std::pair<qType, gmx_stats>(i.first, std::move(gepot)));
        }
        lsq_esp.insert(std::pair<iMolSelect, qtStats>(ims.first, std::move(qesp)));
        lsq_quad.insert(std::pair<iMolSelect, qtStats>(ims.first, std::move(qquad)));
        lsq_dip.insert(std::pair<iMolSelect, qtStats>(ims.first, std::move(qdip)));
        lsq_mu.insert(std::pair<iMolSelect, qtStats>(ims.first, std::move(qmu)));
        lsq_epot.insert(std::pair<iMolSelect, qtStats>(ims.first, std::move(qepot)));
        
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
    
    // Extract terms to print for the enery ters
    std::vector<int> ePlot = { F_EPOT, F_COUL_SR };
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
            fprintf(fp, "\nMolecule %d: Name: %s, Qtot: %d, Multiplicity: %d, Dataset: %s\n", n+1,
                    mol->getMolname().c_str(),
                    mol->totalCharge(),
                    mol->getMultiplicity(),
                    iMolSelectName(ims));

            // Recalculate the atomic charges using the optimized parameters.
            std::vector<double> dummy;
            mol->GenerateCharges(pd, fplog, cr, nullptr, qcycle, qtol,
                                 ChargeGenerationAlgorithm::NONE, dummy, lot);
            // Energy
            double emol = 0;
            if (mol->energy(MolPropObservable::EMOL, &emol))
            {
                lsq_epot[ims][qType::Calc].add_point(emol, mol->potentialEnergy(), 0, 0);
                fprintf(fp, "Reference potential energy %.2f\n", emol);
            }
            auto terms = mol->energyTerms();
            for(auto &ep : ePlot)
            {
                fprintf(fp, "%-20s  %.2f\n", interaction_function[ep].name, terms[ep]);
            }
            // Now compute all the ESP RMSDs.
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
            // Dipoles
            auto qelec = mol->qTypeProps(qType::Elec);
            print_dipole(fp, qType::Elec, qelec->mu(), qelec->mu(), dip_toler_);
            for (auto &j : qTypes())
            {
                qType qt = j.first;
                auto qprop = mol->qTypeProps(qt);
                if (qt != qType::Elec && qprop)
                {
                    if (qt == qType::Calc)
                    {
                        qprop->setX(mol->x());
                    }
                    qprop->setCenterOfCharge(mol->centerOfCharge());
                    qprop->calcMoments();
                    print_dipole(fp, qt, qelec->mu(), qprop->mu(), dip_toler_);
                    if (mol->datasetType() == ims)
                    {
                        lsq_dip[ims][qt].add_point(qelec->dipole(), qprop->dipole(), 0, 0);
                        for(int mm = 0; mm < DIM; mm++)
                        {
                            lsq_mu[ims][qt].add_point(qelec->mu()[mm], qprop->mu()[mm], 0, 0);
                        }
                    }
                }
            }

            // Quadrupoles
            analyse_quadrapole(fp, mol, quad_toler_, lsq_quad[ims], bfullTensor);
            
            // Polarizability
            if (bPolar)
            {
                mol->CalcPolarizability(efield, cr, nullptr);
                print_polarizability(fp, mol, qTypeName(qType::Elec), alpha_toler_, isopol_toler_);
                print_polarizability(fp, mol, qTypeName(qType::Calc), alpha_toler_, isopol_toler_);
                lsq_isoPol[ims][qType::Calc].add_point(mol->ElectronicPolarizability(),
                                                       mol->CalculatedPolarizability(),       0, 0);
                lsq_anisoPol[ims][qType::Calc].add_point(mol->anisoPolElec(), mol->anisoPolCalc(), 0, 0);
                for (int mm = 0; mm < DIM; mm++)
                {
                    lsq_alpha[ims][qType::Calc].add_point(mol->alpha_elec()[mm][mm], mol->alpha_calc()[mm][mm], 0, 0);
                }
            }

            // Atomic charges
            std::map<qType, std::vector<double> > qQM;
            std::vector<qType>                    typeQM = { qType::CM5, qType::ESP, qType::Hirshfeld, qType::Mulliken };
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
            fprintf(fp, "Atom   Type        ACM");
            for(auto &qt : qQM)
            {
                if (!qQM.find(qt.first)->second.empty())
                {
                    fprintf(fp, "%10s", qTypeName(qt.first).c_str());
                }
            }
            fprintf(fp, "        x        y        z (pm)\n");
            auto x       = mol->x();
            auto myatoms = mol->atomsConst();
            int  i       = 0;
            for (int j = i = 0; j < myatoms.nr; j++)
            {
                if (myatoms.atom[j].ptype == eptAtom ||
                    myatoms.atom[j].ptype == eptNucleus)
                {
                    real qCalc = myatoms.atom[j].q;
                    fprintf(fp, "%-2d%3d  %-5s  %8.4f",
                            myatoms.atom[j].atomnumber,
                            j+1,
                            *(myatoms.atomtype[j]),
                            qCalc);
                    for(auto &qt : qQM)
                    {
                        if (!qQM.find(qt.first)->second.empty())
                        {
                            fprintf(fp, "  %8.4f", qt.second[i]);
                        }
                    }
                    fprintf(fp," %8.3f %8.3f %8.3f\n", 
                            convertFromGromacs(x[j][XX], "pm"),
                            convertFromGromacs(x[j][YY], "pm"),
                            convertFromGromacs(x[j][ZZ], "pm"));
                    i++;
                }
                else
                {
                    // Turned on printing of shells again
                    fprintf(fp, "%-2d%3d  %-5s  %8.4f",
                            0,
                            j+1,
                            *(myatoms.atomtype[j]),
                            myatoms.atom[j].q);
                    for(auto &qt : qQM)
                    {
                        if (!qQM.find(qt.first)->second.empty())
                        {
                            fprintf(fp, "          ");
                        }
                    }
                    fprintf(fp," %8.3f %8.3f %8.3f\n", 
                            convertFromGromacs(x[j][XX], "pm"),
                            convertFromGromacs(x[j][YY], "pm"),
                            convertFromGromacs(x[j][ZZ], "pm"));
                }
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
            print_stats(fp, "ESP  (kJ/mol e)", &lsq_esp[ims.first][qt],  header, "Electronic", name, useOffset_);
            header = false;
            print_stats(fp, "Dipoles",         &lsq_mu[ims.first][qt],   header, "Electronic", name, useOffset_);
            print_stats(fp, "Dipole Moment",   &lsq_dip[ims.first][qt],  header, "Electronic", name, useOffset_);
            print_stats(fp, "Quadrupoles",     &lsq_quad[ims.first][qt], header, "Electronic", name, useOffset_);
            if (bPolar && qt == qType::Calc)
            {
                print_stats(fp, "Polariz. components (A^3)",  &lsq_alpha[ims.first][qType::Calc],    header, "Electronic", name, useOffset_);
                print_stats(fp, "Isotropic Polariz. (A^3)",   &lsq_isoPol[ims.first][qType::Calc],   header, "Electronic", name, useOffset_);
                print_stats(fp, "Anisotropic Polariz. (A^3)", &lsq_anisoPol[ims.first][qType::Calc], header, "Electronic", name, useOffset_);
            }
            fprintf(fp, "\n");
        }
    }
    write_q_histo(fp, opt2fn("-qhisto", filenm.size(), filenm.data()),
                  lsqt[iMolSelect::Train], oenv, &(lsq_charge[iMolSelect::Train]), useOffset_);
    
    print_corr(opt2fn("-dipcorr", filenm.size(), filenm.data()),
               "Dipole Moment (Debye)", "Electronic", "Empirical", lsq_dip, oenv);
    print_corr(opt2fn("-mucorr", filenm.size(), filenm.data()),
               "Dipole components (Debye)", "Electronic", "Empirical", lsq_mu, oenv);
    print_corr(opt2fn("-epotcorr", filenm.size(), filenm.data()),
               "Potential energy (kJ/mol)", "Reference", "Empirical", lsq_epot, oenv);
    print_corr(opt2fn("-thetacorr", filenm.size(), filenm.data()),
               "Quadrupole components (Buckingham)", "Electronic", "Empirical", lsq_quad, oenv);
    print_corr(opt2fn("-espcorr", filenm.size(), filenm.data()),
               "Electrostatic Potential (Hartree/e)", "Electronic", "Calc", lsq_esp, oenv);
    print_corr(opt2fn("-qcorr", filenm.size(), filenm.data()),
               "Atomic Partial Charge", "q (e)", "a.u.", lsq_charge, oenv);
    
    if (bPolar)
    {
        print_corr(opt2fn("-alphacorr", filenm.size(), filenm.data()),
                   "Pricipal Components of Polarizability Tensor (A\\S3\\N)", "Electronic", "Calc", lsq_alpha, oenv);
        print_corr(opt2fn("-isopol", filenm.size(), filenm.data()),
                   "Isotropic Polarizability (A\\S3\\N)", "Electronic", "Calc", lsq_isoPol, oenv);
        print_corr(opt2fn("-anisopol", filenm.size(), filenm.data()),
                   "Anisotropic Polarizability (A\\S3\\N)", "Electronic", "Calc", lsq_anisoPol, oenv);
    }
    // List outliers based on the deviation in the total dipole moment
    real sigma;
    if (eStats::OK == lsq_dip[iMolSelect::Train][qType::Calc].get_rmsd(&sigma))
    {
        fprintf(fp, "Overview of dipole moment outliers (> %.3f off)\n", 2*sigma);
        fprintf(fp, "----------------------------------\n");
        fprintf(fp, "%-40s  %12s  %12s  %12s\n", "Name", "Calc", "Electronic", "Deviation (Debye)");
        for (auto mol = mymol->begin(); mol < mymol->end(); ++mol)
        {
            auto qelec = mol->qTypeProps(qType::Elec);
            auto qcalc = mol->qTypeProps(qType::Calc);
            auto deviation = std::abs(qcalc->dipole() - qelec->dipole());
            if ((mol->support() != eSupport::No)  &&
                (qelec->dipole() > sigma) &&
                (deviation > 2*sigma))
            {
                fprintf(fp, "%-40s  %12.3f  %12.3f  %12.3f  %s\n",
                        mol->getMolname().c_str(),
                        qcalc->dipole(),qelec->dipole(), deviation,
                        iMolSelectName(mol->datasetType()));
            }
        }
    }
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
            printf("There were %d ESP outliers. See at the very bottom of the log file\n", nout);
        }
        else
        {
            printf("No outliers! Well done.\n");
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
    fprintf(fp, "%-15s  %s\n", "Option", "Value");
    for (auto &p: pargs)
    {
        std::string value;
        switch(p.type)
        {
        case etINT:
            value = gmx::formatString("%d", *p.u.i);
            break;
        case etINT64:
            value = gmx::formatString("%" PRId64, *p.u.is);
            break;
        case etREAL:
            value = gmx::formatString("%g", *p.u.r);
            break;
        case etSTR:
            value = *p.u.c;
            break;
        case etENUM:
            value = gmx::formatString("%s", *p.u.c);
            break;
        case etBOOL:
            value = gmx::formatString("%s", *p.u.b ? "true" : "false");
            break;
        case etRVEC:
            value = gmx::formatString("%g %g %g", *p.u.rv[XX],
                                      *p.u.rv[YY], *p.u.rv[ZZ]);
            break;
        case etTIME:
        case etNR:
        default:
            value.assign("help");
        }
        fprintf(fp, "%-15s  %s\n", p.option, value.c_str());
    }
}

} // namespace alexandria
