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

#include "tuning_utility.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "gromacs/utility/coolstuff.h"

#include "molprop_util.h"
#include "mymol.h"

namespace alexandria
{

static void print_stats(FILE        *fp,
                        const char  *prop,
                        gmx_stats_t  lsq,
                        bool         bHeader,
                        const char  *xaxis,
                        const char  *yaxis,
                        bool         useOffset)
{
    real a    = 0, da  = 0, b    = 0, db   = 0;
    real mse  = 0, mae = 0, chi2 = 0, rmsd = 0;
    real Rfit = 0;
    int  n;

    if (useOffset)
    {
        if (bHeader)
        {
            fprintf(fp, "Fitting data to y = ax + b, where x = %s and y = %s\n", xaxis, yaxis);
            fprintf(fp, "%-12s %5s %13s %13s %8s %8s %8s %8s\n",
                    "Property", "N", "a", "b", "R(%)", "RMSD", "MSE", "MAE");
            fprintf(fp, "---------------------------------------------------------------\n");
        }
        gmx_stats_get_ab(lsq, elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &Rfit);
        gmx_stats_get_rmsd(lsq,    &rmsd);
        gmx_stats_get_mse_mae(lsq, &mse, &mae);
        gmx_stats_get_npoints(lsq, &n);
        fprintf(fp, "%-12s %5d %6.3f(%5.3f) %6.3f(%5.3f) %7.2f %8.4f %8.4f %8.4f\n",
                prop, n, a, da, b, db, Rfit*100, rmsd, mse, mae);
    }
    else
    {
        if (bHeader)
        {
            fprintf(fp, "Fitting data to y = ax, where x = %s and y = %s\n", xaxis, yaxis);
            fprintf(fp, "%-12s %5s %13s %8s %8s %8s %8s\n",
                    "Property", "N", "a", "R(%)", "RMSD", "MSE", "MAE");
            fprintf(fp, "---------------------------------------------------------------\n");
        }
        gmx_stats_get_a(lsq, elsqWEIGHT_NONE, &a, &da, &chi2, &Rfit);
        gmx_stats_get_rmsd(lsq,    &rmsd);
        gmx_stats_get_mse_mae(lsq, &mse, &mae);
        gmx_stats_get_npoints(lsq, &n);
        fprintf(fp, "%-12s %5d %6.3f(%5.3f) %7.2f %8.4f %8.4f %8.4f\n",
                prop, n, a, da, Rfit*100, rmsd, mse, mae);
    }
}

void print_lsq_set(FILE *fp, gmx_stats_t lsq)
{
    real   x, y;

    fprintf(fp, "@type xy\n");
    while (gmx_stats_get_point(lsq, &x, &y, nullptr, nullptr, 0) == estatsOK)
    {
        fprintf(fp, "%10g  %10g\n", x, y);
    }
    fprintf(fp, "&\n");
}

void xvgr_symbolize(FILE                   *xvgf,
                    int                     nsym,
                    const char             *leg[],
                    const gmx_output_env_t *oenv)
{
    int i;

    xvgr_legend(xvgf, nsym, leg, oenv);
    for (i = 0; (i < nsym); i++)
    {
        xvgr_line_props(xvgf, i, elNone, ecBlack+i, oenv);
        fprintf(xvgf, "@ s%d symbol %d\n", i, i+1);
    }
}

void print_polarizability(FILE              *fp,
                          alexandria::MyMol *mol,
                          char              *calc_name,
                          real               alpha_toler,
                          real               isopol_toler)
{
    tensor dalpha;
    real   delta      = 0;
    real   diso_pol   = 0;
    real   daniso_pol = 0;

    if (nullptr != calc_name)
    {
        if (strcmp(calc_name, (char *)"Calculated") == 0)
        {
            m_sub(mol->alpha_elec_, mol->alpha_calc_, dalpha);
            delta = sqrt(gmx::square(dalpha[XX][XX])+gmx::square(dalpha[XX][YY])+gmx::square(dalpha[XX][ZZ])+
                         gmx::square(dalpha[YY][YY])+gmx::square(dalpha[YY][ZZ]));
            diso_pol   = std::abs(mol->PolarizabilityDeviation());
            daniso_pol = std::abs(mol->AnisoPolarizabilityDeviation());
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->alpha_calc_[XX][XX], mol->alpha_calc_[XX][YY], mol->alpha_calc_[XX][ZZ],
                    dalpha[XX][XX], dalpha[XX][YY], dalpha[XX][ZZ], delta, (delta > alpha_toler) ? "ALPHA" : "",
                    "", mol->alpha_calc_[YY][YY], mol->alpha_calc_[YY][ZZ],
                    "", dalpha[YY][YY], dalpha[YY][ZZ],
                    "", "", mol->alpha_calc_[ZZ][ZZ],
                    "", "", dalpha[ZZ][ZZ]);
            fprintf(fp,
                    "Isotropic polarizability:  %s Electronic: %6.2f  Calculated: %6.2f  Delta: %6.2f %s\n\n\n",
                    mol->getMolname().c_str(), 
                    mol->ElectronicPolarizability(),
                    mol->CalculatedPolarizability(), 
                    diso_pol, (diso_pol > isopol_toler) ? "ISO" : "");
            fprintf(fp,
                    "Anisotropic polarizability:  %s Electronic: %6.2f  Calculated: %6.2f  Delta: %6.2f %s\n\n\n",
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
                    mol->alpha_elec_[XX][XX], mol->alpha_elec_[XX][YY], mol->alpha_elec_[XX][ZZ],
                    "", mol->alpha_elec_[YY][YY], mol->alpha_elec_[YY][ZZ],
                    "", "", mol->alpha_elec_[ZZ][ZZ]);
        }
    }
}

static void print_quadrapole(FILE              *fp,
                             alexandria::MyMol *mol,
                             qType              qt,
                             real               q_toler)
{
    const tensor &qelec = mol->QQM(qtElec);

    if (qt != qtElec)
    {
        const tensor &qcalc = mol->QQM(qt);
        tensor        dQ;
        real          delta = 0;

        m_sub(qelec, qcalc, dQ);
        delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                     gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
        fprintf(fp,
                "%-10s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                "           (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                "           (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                qTypeName(qt),
                qcalc[XX][XX], qcalc[XX][YY], qcalc[XX][ZZ],
                dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "QUAD" : "",
                "", qcalc[YY][YY], qcalc[YY][ZZ],
                "", dQ[YY][YY], dQ[YY][ZZ],
                "", "", qcalc[ZZ][ZZ],
                "", "", dQ[ZZ][ZZ]);
    }
    else
    {
        fprintf(fp, "Quadrupole analysis (5 independent components only)\n");
        fprintf(fp,
                "Electronic   (%6.2f %6.2f %6.2f)\n"
                "             (%6s %6.2f %6.2f)\n"
                "             (%6s %6s %6.2f)\n",
                qelec[XX][XX], qelec[XX][YY], qelec[XX][ZZ],
                "", qelec[YY][YY], qelec[YY][ZZ],
                "", "", qelec[ZZ][ZZ]);
    }
}

static void print_dipole(FILE              *fp,
                         alexandria::MyMol *mol,
                         qType              qt,
                         real               toler)
{
    rvec dmu;
    real ndmu, cosa;
    char ebuf[32];

    rvec_sub(mol->muQM(qtElec), mol->muQM(qt), dmu);
    ndmu = norm(dmu);
    cosa = cos_angle(mol->muQM(qtElec), mol->muQM(qtCalc));
    if (ndmu > toler)
    {
        sprintf(ebuf, "DIP");
    }
    else if (fabs(cosa) < 0.1)
    {
        sprintf(ebuf, "YYY");
    }
    else
    {
        ebuf[0] = '\0';
    }
    fprintf(fp, "%-10s (%6.3f,%6.3f,%6.3f) |Mu| = %6.3f",
            qTypeName(qt), mol->muQM(qt)[XX], mol->muQM(qt)[YY], mol->muQM(qt)[ZZ],
            mol->dipQM(qt));
    if (qt != qtElec)
    {
        fprintf(fp, " Dev: (%6.3f,%6.3f,%6.3f) |%6.3f|%s",
                dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
    }
    fprintf(fp, "\n");
}

void print_electric_props(FILE                           *fp,
                          std::vector<alexandria::MyMol> &mymol,
                          const Poldata                  *pd,
                          const gmx::MDLogger            &fplog,
                          gmx_atomprop_t                  ap,
                          real                            hfac,
                          const char                     *lot,
                          const char                     *tabfn,
                          gmx_hw_info_t                  *hwinfo,
                          int                             qcycle,
                          real                            qtol,
                          const char                     *qhisto,
                          const char                     *DipCorr,
                          const char                     *MuCorr,
                          const char                     *Qcorr,
                          const char                     *EspCorr,
                          const char                     *alphaCorr,
                          const char                     *isopolCorr,
                          const char                     *anisopolCorr,
                          const char                     *qCorr,
                          real                            esp_toler,
                          real                            dip_toler,
                          real                            quad_toler,
                          real                            alpha_toler,
                          real                            isopol_toler,
                          const gmx_output_env_t         *oenv,
                          bool                            bPolar,
                          bool                            bDipole,
                          bool                            bQuadrupole,
                          bool                            bfullTensor,
                          IndexCount                     *indexCount,
                          t_commrec                      *cr,
                          real                            efield,
                          bool                            useOffset)
{
    int            i    = 0, j     = 0, n     = 0;
    int            nout = 0, mm    = 0, nn    = 0;
    real           sse  = 0, sigma = 0, qCalc = 0;

    FILE          *dipc, *muc,  *Qc;
    FILE          *hh,   *espc, *alphac, *isopolc, *anisopolc, *qc;

    struct ZetaTypeLsq {
        std::string ztype;
        gmx_stats_t lsq;
    };

    gmx_stats_t               lsq_mu[qtNR], lsq_dip[qtNR], lsq_quad[qtNR];
    gmx_stats_t               lsq_esp, lsq_alpha, lsq_isoPol, lsq_anisoPol, lsq_charge;
    const char               *eprnm[qtNR];
    std::vector<ZetaTypeLsq>  lsqt;

    for (int i = 0; i < qtNR; i++)
    {
        eprnm[i]    = qTypeName(static_cast<qType>(i));
        lsq_quad[i] = gmx_stats_init();
        lsq_dip[i]  = gmx_stats_init();
        lsq_mu[i]   = gmx_stats_init();
    }
    lsq_esp      = gmx_stats_init();
    lsq_alpha    = gmx_stats_init();
    lsq_isoPol   = gmx_stats_init();
    lsq_anisoPol = gmx_stats_init();
    lsq_charge   = gmx_stats_init();
    n            = 0;

    for (auto ai = indexCount->beginIndex(); ai < indexCount->endIndex(); ++ai)
    {
        // More than 1 data point is needed for doing statistical analysis for each zeta type.
        if (ai->count() > 1)
        {
            ZetaTypeLsq k;
            k.ztype.assign(ai->name());
            k.lsq = gmx_stats_init();
            lsqt.push_back(std::move(k));
        }
    }
    std::string method, basis;
    splitLot(lot, &method, &basis);
    for (auto &mol : mymol)
    {
        if (mol.eSupp_ != eSupportNo)
        {
            fprintf(fp, "Molecule %d: %s Qtot: %d, Multiplicity %d\n", n+1,
                    mol.getMolname().c_str(),
                    mol.getCharge(),
                    mol.getMultiplicity());

            // Recalculate the atomic charges using the optmized parameters.
            mol.GenerateCharges(pd, fplog, ap, hfac, false, nullptr,
                                cr, tabfn, hwinfo, qcycle, qtol);

            // Electrostatic potentials
            mol.Qgresp_.updateZeta(mol.atoms_, pd);
            mol.Qgresp_.updateAtomCharges(mol.atoms_);
            mol.Qgresp_.updateAtomCoords(mol.x());
            mol.Qgresp_.calcPot(pd->getEpsilonR());
            real rrms    = 0, wtot = 0, cosangle = 0;
            bool espZero = false;
            for (size_t i = 0; i < mol.Qgresp_.nEsp(); i++)
            {
                auto vqm   = mol.Qgresp_.espPoint()[i].v();
                auto valex = mol.Qgresp_.espPoint()[i].vCalc();
                if (std::abs(vqm) > 100 && std::abs(valex) < 10)
                {
                    espZero = true;
                }
            }
            auto rms    = mol.Qgresp_.getRms(&wtot, &rrms, &cosangle);
            auto espRms = convert2gmx(rms, eg2cHartree_e);
            fprintf(fp, "ESP rms: %g (kJ/mol e) %s\n", espRms, (espRms > esp_toler) ? " EEE" : "");
            fprintf(fp, "ESP cosangle: %.3f%s\n", cosangle, (cosangle < 0.5) ? " WWW" : "");
            if (espZero)
            {
                fprintf(fp, "Orthogonal ESP\n");
            }
            for (size_t i = 0; i < mol.Qgresp_.nEsp(); i++)
            {
                gmx_stats_add_point(lsq_esp, mol.Qgresp_.espPoint()[i].v(), mol.Qgresp_.espPoint()[i].vCalc(), 0, 0);
            }

            // Dipoles
            mol.CalcDipole();
            mol.rotateDipole(mol.muQM(qtCalc), mol.muQM(qtElec));
            print_dipole(fp, &mol, qtElec, dip_toler);
            for (int j = 0; j < qtElec; j++)
            {
                qType qt = static_cast<qType>(j);
                print_dipole(fp, &mol, qt,   dip_toler);
                gmx_stats_add_point(lsq_dip[j], mol.dipQM(qtElec), mol.dipQM(qt), 0, 0);
            }
            sse += gmx::square(mol.dipQM(qtElec) - mol.dipQM(qtCalc));

            // Quadrupoles
            mol.CalcQuadrupole();
            print_quadrapole(fp, &mol, qtElec, quad_toler);
            for (int j = 0; j < qtElec; j++)
            {
                qType qt = static_cast<qType>(j);
                print_quadrapole(fp, &mol, qt, quad_toler);
                for (mm = 0; mm < DIM; mm++)
                {
                    gmx_stats_add_point(lsq_mu[j], mol.muQM(qtElec)[mm], mol.muQM(qt)[mm], 0, 0);
                    for (nn = 0; nn < DIM; nn++)
                    {
                        if (bfullTensor || (mm == nn))
                        {
                            gmx_stats_add_point(lsq_quad[j], mol.QQM(qtElec)[mm][nn], mol.QQM(qt)[mm][nn], 0, 0);
                        }
                    }
                }
            }

            // Polarizability
            if (bPolar)
            {
                mol.CalcPolarizability(efield, cr, nullptr);
                print_polarizability(fp, &mol, (char *)"Electronic", alpha_toler, isopol_toler);
                print_polarizability(fp, &mol, (char *)"Calculated", alpha_toler, isopol_toler);
                gmx_stats_add_point(lsq_isoPol, mol.ElectronicPolarizability(),
                                    mol.CalculatedPolarizability(),       0, 0);
                gmx_stats_add_point(lsq_anisoPol, mol.anisoPol_elec_, mol.anisoPol_calc_, 0, 0);
                for (mm = 0; mm < DIM; mm++)
                {
                    gmx_stats_add_point(lsq_alpha, mol.alpha_elec_[mm][mm], mol.alpha_calc_[mm][mm], 0, 0);
                }
            }

            // Atomic charges
            fprintf(fp, "Atom   Type      q_Calc     q_ESP     q_CM5     q_HPA     q_MPA       x       y       z\n");
            auto qcm5  = mol.chargeQM(qtCM5);
            auto qESP  = mol.chargeQM(qtESP);
            auto qHir  = mol.chargeQM(qtHirshfeld);
            auto qMul  = mol.chargeQM(qtMulliken);
            auto x     = mol.x();
            auto qrmsd = 0.0;
            for (j = i = 0; j < mol.atoms_->nr; j++)
            {
                if (mol.atoms_->atom[j].ptype == eptAtom ||
                    mol.atoms_->atom[j].ptype == eptNucleus)
                {
                    auto  atp = pd->findAtype(*(mol.atoms_->atomtype[j]));
                    auto  ztp = atp->getZtype();
                    auto  k   = std::find_if(lsqt.begin(), lsqt.end(),
                                             [ztp](const ZetaTypeLsq &atlsq)
                                             {
                                                 return atlsq.ztype.compare(ztp) == 0;
                                             });
                    qCalc = mol.atoms_->atom[j].q;
                    if (nullptr != mol.shellfc_)
                    {
                        //qCalc += mol.atoms_->atom[j+1].q;
                    }
                    if (k != lsqt.end())
                    {
                        gmx_stats_add_point(k->lsq, qcm5[i], qCalc, 0, 0);
                    }
                    gmx_stats_add_point(lsq_charge, qcm5[i], qCalc, 0, 0);
                    qrmsd += gmx::square(qcm5[i]-qCalc);
                    fprintf(fp, "%-2d%3d  %-5s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f %8.3f%8.3f%8.3f\n",
                            mol.atoms_->atom[j].atomnumber,
                            j+1,
                            *(mol.atoms_->atomtype[j]),
                            qCalc,
                            qESP.size() > 0 ? qESP[i] : 0.0,
                            qcm5.size() > 0 ? qcm5[i] : 0.0,
                            qHir.size() > 0 ? qHir[i] : 0.0,
                            qMul.size() > 0 ? qMul[i] : 0.0,
                            x[j][XX],
                            x[j][YY],
                            x[j][ZZ]);
                    i++;
                }
                else
                {
                    fprintf(fp, "%-2d%3d  %-5s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f %8.3f%8.3f%8.3f\n",
                            0,
                            j+1,
                            *(mol.atoms_->atomtype[j]),
                            mol.atoms_->atom[j].q,
                            0.0, 0.0, 0.0, 0.0,
                            x[j][XX],
                            x[j][YY],
                            x[j][ZZ]);
                }
            }
            fprintf(fp, "\n");
            qrmsd = sqrt(qrmsd/mol.atoms_->nr);
            fprintf(fp, "q rmsd: %g (e) %s\n", qrmsd, (qrmsd > 5e-2) ? "XXX" : "");
            n++;
        }
    }

    fprintf(fp, "Dipoles are %s in Calc Parametrization.\n",     (bDipole ?     "used" : "not used"));
    fprintf(fp, "Quadrupoles are %s in Calc Parametrization.\n", (bQuadrupole ? "used" : "not used"));
    fprintf(fp, "\n");

    print_stats(fp, "ESP  (kJ/mol e)",  lsq_esp, true,  "Electronic", "Calculated", useOffset);
    fprintf(fp, "\n");

    for (int i = 0; i < qtElec; i++)
    {
        const char *name = qTypeName(static_cast<qType>(i));
        print_stats(fp, "Dipoles",       lsq_mu[i],   true,  "Electronic", name, useOffset);
        print_stats(fp, "Dipole Moment", lsq_dip[i],  false, "Electronic", name, useOffset);
        print_stats(fp, "Quadrupoles",   lsq_quad[i], false, "Electronic", name, useOffset);
        if (bPolar && i == qtCalc)
        {
            print_stats(fp, "Principal Components of Polarizability (A^3)",  lsq_alpha, false,  "Electronic", name, useOffset);
            print_stats(fp, "Isotropic Polarizability (A^3)",    lsq_isoPol,   false,  "Electronic", name, useOffset);
            print_stats(fp, "Anisotropic Polarizability (A^3)",  lsq_anisoPol, false,  "Electronic", name, useOffset);
        }
        fprintf(fp, "\n");
    }

    std::vector<const char*> types;
    for (const auto &k : lsqt)
    {
        types.push_back(k.ztype.c_str());
    }
    hh = xvgropen(qhisto, "Histogram for charges", "q (e)", "a.u.", oenv);
    xvgr_legend(hh, types.size(), types.data(), oenv);
    print_stats(fp, "All Partial Charges  (e)",  lsq_charge, true,  "CM5", "Calculated", useOffset);
    for (auto k = lsqt.begin(); k < lsqt.end(); ++k)
    {
        int   nbins;
        if (gmx_stats_get_npoints(k->lsq, &nbins) == estatsOK)
        {
            real *x, *y;
            if (gmx_stats_make_histogram(k->lsq, 0, &nbins, ehistoY, 1, &x, &y) == estatsOK)
            {
                fprintf(hh, "@type xy\n");
                for (int i = 0; i < nbins; i++)
                {
                    fprintf(hh, "%10g  %10g\n", x[i], y[i]);
                }
                fprintf(hh, "&\n");
                print_stats(fp, k->ztype.c_str(),  k->lsq, false,  "CM5", "Calculated", useOffset);
                free(x);
                free(y);
            }
        }
        gmx_stats_free(k->lsq);
    }
    fclose(hh);
    fprintf(fp, "\n");

    dipc = xvgropen(DipCorr, "Dipole Moment (Debye)", "Electronic", "Empirical", oenv);
    xvgr_symbolize(dipc, 5, eprnm, oenv);
    for (int i = 0; i < qtElec; i++)
    {
        print_lsq_set(dipc, lsq_dip[i]);
    }
    fclose(dipc);

    muc = xvgropen(MuCorr, "Dipoles (Debye)", "Electronic", "Empirical", oenv);
    xvgr_symbolize(muc, 5, eprnm, oenv);
    for (int i = 0; i < qtElec; i++)
    {
        print_lsq_set(muc, lsq_mu[i]);
    }
    fclose(muc);

    Qc = xvgropen(Qcorr, "Quadrupoles (Buckingham)", "Electronic", "Empirical", oenv);
    xvgr_symbolize(Qc, 5, eprnm, oenv);
    for (int i = 0; i < qtElec; i++)
    {
        print_lsq_set(Qc, lsq_quad[i]);
    }
    fclose(Qc);

    espc = xvgropen(EspCorr, "Electrostatic Potential (Hartree/e)", "Electronic", "Calc", oenv);
    xvgr_symbolize(espc, 1, eprnm, oenv);
    print_lsq_set(espc, lsq_esp);
    fclose(espc);

    qc = xvgropen(qCorr, "Atomic Partial Charge", "q (e)", "a.u.", oenv);
    xvgr_symbolize(qc, 1, eprnm, oenv);
    print_lsq_set(qc, lsq_charge);
    fclose(qc);

    if (bPolar)
    {
        alphac = xvgropen(alphaCorr, "Pricipal Components of Polarizability Tensor (A\\S3\\N)", "Electronic", "Calc", oenv);
        xvgr_symbolize(alphac, 1, eprnm, oenv);
        print_lsq_set(alphac, lsq_alpha);
        fclose(alphac);

        isopolc = xvgropen(isopolCorr, "Isotropic Polarizability (A\\S3\\N)", "Electronic", "Calc", oenv);
        xvgr_symbolize(isopolc, 1, eprnm, oenv);
        print_lsq_set(isopolc, lsq_isoPol);
        fclose(isopolc);

        anisopolc = xvgropen(anisopolCorr, "Anisotropic Polarizability (A\\S3\\N)", "Electronic", "Calc", oenv);
        xvgr_symbolize(anisopolc, 1, eprnm, oenv);
        print_lsq_set(anisopolc, lsq_anisoPol);
        fclose(anisopolc);
    }

    // List outliers based on the deviation in the total dipole moment
    sigma = sqrt(sse/n);
    fprintf(fp, "Overview of dipole moment outliers (> %.3f off)\n", 2*sigma);
    fprintf(fp, "----------------------------------\n");
    fprintf(fp, "%-20s  %12s  %12s  %12s\n", "Name", "Calc", "Electronic", "Deviation (Debye)");
    for (auto &mol : mymol)
    {
        auto deviation = std::abs(mol.dipQM(qtCalc) - mol.dipQM(qtElec));
        if ((mol.eSupp_ != eSupportNo)  &&
            (mol.dipQM(qtElec) > sigma) &&
            (deviation > 2*sigma))
        {
            fprintf(fp, "%-20s  %12.3f  %12.3f  %12.3f\n",
                    mol.getMolname().c_str(),
                    mol.dipQM(qtCalc), mol.dipQM(qtElec), deviation);
            nout++;
        }
    }
    if (nout)
    {
        printf("There were %d outliers. See at the very bottom of the log file\n", nout);
    }
    else
    {
        printf("No outliers! Well done.\n");
    }

    // Free allocated memory
    for (int i = 0; i < qtNR; i++)
    {
        gmx_stats_free(lsq_quad[i]);
        gmx_stats_free(lsq_mu[i]);
        gmx_stats_free(lsq_dip[i]);
    }
    gmx_stats_free(lsq_esp);
    gmx_stats_free(lsq_alpha);
    gmx_stats_free(lsq_isoPol);
    gmx_stats_free(lsq_anisoPol);
    gmx_stats_free(lsq_charge);
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
    fprintf(fp, "# alexandria is part of GROMACS:\n#\n");
    fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());
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
            value = gmx::formatString("%ld", *p.u.is);
            break;
        case etREAL:
            value = gmx::formatString("%g", *p.u.r);
            break;
        case etSTR:
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
