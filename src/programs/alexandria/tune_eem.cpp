/*
 * This source file is part of the Alexandria program.
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

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <random>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "alex_modules.h"
#include "gentop_core.h"
#include "gmx_simple_comm.h"
#include "molgen.h"
#include "mymol_low.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_tables.h"
#include "poldata_xml.h"
#include "tuning_utility.h"

namespace alexandria
{

class OptACM : public MolGen
{
    using param_type = std::vector<double>;

    private:
        gmx_bool       bFullTensor_;
        gmx_bool       bFitAlpha_;
        gmx_bool       bFitZeta_;
        gmx_bool       bSameZeta_;
        gmx_bool       bFitChi_;
        gmx_bool       bUseCM5_;

        Bayes <double> TuneACM_;
        param_type     param_, lower_, upper_, best_;
        param_type     orig_, psigma_, pmean_;

        real           penalty_;

    public:

        OptACM()
            :
              bFullTensor_(false),
              bFitAlpha_(false),
              bFitZeta_(false),
              bSameZeta_(false),
              bFitChi_(true),
              bUseCM5_(false),
              penalty_(0)
        {}

        ~OptACM() {}

        gmx_bool bESP() const { return weight(ermsESP); }

        gmx_bool dipole() const { return weight(ermsMU); }

        gmx_bool quadrupole() const { return weight(ermsQUAD); }

        gmx_bool fullTensor() const { return bFullTensor_; }

        gmx_bool fitZeta() const { return bFitZeta_; }
       
        gmx_bool sameZeta() const { return bSameZeta_; }
        
        gmx_bool fitChi() const { return bFitChi_; }

        gmx_bool useCM5() const {return bUseCM5_; }

        gmx_bool penalize() const {return true ? (penalty_ > 0) : false; }
        
        void add_pargs(std::vector<t_pargs> *pargs)
        {
            t_pargs pa[] =
            {
                { "-fullTensor", FALSE, etBOOL, {&bFullTensor_},
                  "Consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization" },
                { "-fitalpha", FALSE, etBOOL, {&bFitAlpha_},
                  "Calibrate atomic polarizability." },
                { "-fitzeta", FALSE, etBOOL, {&bFitZeta_},
                  "Calibrate orbital exponent." },
                { "-samezeta", FALSE, etBOOL, {&bSameZeta_},
                  "Use the same zeta for both the core and the shell of the Drude model." },
                { "-fitchi", FALSE, etBOOL, {&bFitChi_},
                  "Calibrate electronegativity and hardness." },
                { "-penalty", FALSE, etREAL, {&penalty_},
                  "penalty to keep the Chi0 and J0 in order." },
                { "-cm5", FALSE, etBOOL, {&bUseCM5_},
                  "Reproduce CM5 charges in fitting." },
            };
            for (size_t i = 0; i < asize(pa); i++)
            {
                pargs->push_back(pa[i]);
            }
            addOptions(pargs, etuneEEM);
            TuneACM_.add_pargs(pargs);
        }
        void optionsFinished()
        {
            MolGen::optionsFinished();
            TuneACM_.setBounds(weight(ermsBOUNDS) > 0);
        }

        double l2_regularizer (double x,
                               double min,
                               double max)
        {
            return (x < min) ? (0.5 * gmx::square(x-min)) : ((x > max) ? (0.5 * gmx::square(x-max)) : 0);
        }

        void initQgresp();

        void initQgacm();

        void polData2TuneACM();

        void TuneACM2PolData();

        void InitOpt(real factor);

        void calcDeviation();

        double objFunction(const double v[]);

        double calcPenalty(AtomIndexIterator ai);

        /*! \brief
         * Do the actual optimization.
         * \param[in] fp     FILE pointer for logging
         * \param[in] fplog  FILE pointer for logging
         * \param[in] oenv   Output environment for managing xvg files etc.
         * \param[in] xvgconv Output file monitoring parameters
         * \param[in] xvgepot Output file monitoring penalty function
         * \return true if better parameters were found.
         */
        bool optRun(FILE                   *fp,
                    FILE                   *fplog,
                    int                     nrun,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot);
};

void OptACM::initQgresp()
{
    for (auto &mymol : mymols())
    {
        if (mymol.eSupp_ != eSupportNo)
        {
            mymol.initQgresp(poldata(), 
                             lot(),
                             watoms(), 
                             maxPot());
        }
    }
}

void OptACM::initQgacm()
{
    bool bHaveShells = false;
    for (auto &mymol : mymols())
    {
        if (mymol.eSupp_ != eSupportNo)
        {
            bHaveShells = false;
            if (nullptr != mymol.shellfc_)
            {
                bHaveShells = true;
            }
            mymol.Qgacm_.setInfo(poldata(), 
                                 &(mymol.topology_->atoms),
                                 hfac(),
                                 mymol.molProp()->getCharge(),
                                 bHaveShells);
        }
    }
}

void OptACM::calcDeviation()
{
    int                 i         = 0;
    int                 j         = 0;
    int                 iter      = 0;
    double              qtot      = 0;
    double              EemRms    = 0;
    bool                converged = false;
    std::vector<double> qq;

    if (PAR(commrec()))
    {
        bool bFinal = final();
        gmx_bcast(sizeof(final()), &bFinal, commrec());
        if (bFinal)
        {
            setFinal();
        }
    }
    if (PAR(commrec()) && !final())
    {
        poldata()->broadcast(commrec());
    }
    resetEnergies();
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupportLocal) ||
            (final() && (mymol.eSupp_ == eSupportRemote)))
        {
            auto q     = mymol.Qgacm_.q();
            auto natom = mymol.Qgacm_.natom();

            qq.resize(natom + 1);
            for (auto i = 0; i < natom + 1; i++)
            {
                qq[i] = q[i][0];
            }

            converged = false;
            iter      = 0;
            do
            {

                if (nullptr != mymol.shellfc_)
                {
                    if (bFitAlpha_)
                    {
                        mymol.UpdateIdef(poldata(), eitPOLARIZATION);
                    }
                    mymol.computeForces(nullptr, commrec());
                }

                for (auto i = 0; i < mymol.mtop_->natoms; i++)
                {
                    mymol.mtop_->moltype[0].atoms.atom[i].q      =
                        mymol.mtop_->moltype[0].atoms.atom[i].qB = mymol.topology_->atoms.atom[i].q;
                }

                mymol.Qgacm_.generateCharges(debug,
                                             mymol.molProp()->getMolname().c_str(),
                                             poldata(),
                                             &(mymol.topology_->atoms),
                                             mymol.x());

                q       = mymol.Qgacm_.q();
                EemRms  = 0;
                for (auto i = 0; i < natom + 1; i++)
                {
                    EemRms   += gmx::square(qq[i] - q[i][0]);
                    qq[i]     = q[i][0];
                }
                EemRms   /= natom;
                converged = (EemRms < qtol()) || (nullptr == mymol.shellfc_);
                iter++;
            }
            while ((!converged) && (iter < qcycle()));
            for (auto i = 0; i < mymol.mtop_->natoms; i++)
            {
                mymol.mtop_->moltype[0].atoms.atom[i].q      =
                    mymol.mtop_->moltype[0].atoms.atom[i].qB = mymol.topology_->atoms.atom[i].q;
            }

            if (weight(ermsCHARGE))
            {
                int    nChargeResidual = 0; // number of charge residuals added per molecule
                double ChargeResidual  = 0;
                qtot = 0;
                for (j = i = 0; j < mymol.topology_->atoms.nr; j++)
                {
                    auto atomnr = mymol.topology_->atoms.atom[j].atomnumber;
                    auto qq     = mymol.topology_->atoms.atom[j].q;
                    qtot       += qq;
                    if (mymol.topology_->atoms.atom[j].ptype == eptAtom ||
                        mymol.topology_->atoms.atom[j].ptype == eptNucleus)
                    {
                        auto q_H        = 0 ? (nullptr != mymol.shellfc_) : 1;
                        auto q_OFSClBrI = 0 ? (nullptr != mymol.shellfc_) : 2;
                        if (((qq < q_H) && (atomnr == 1)) ||
                            ((qq > q_OFSClBrI) && ((atomnr == 8)  || (atomnr == 9) ||
                                                   (atomnr == 16) || (atomnr == 17) ||
                                                   (atomnr == 35) || (atomnr == 53))))
                        {
                            ChargeResidual += gmx::square(qq);
                            nChargeResidual++;
                        }
                        if (useCM5())
                        {
                            if (nullptr != mymol.shellfc_)
                            {
                                qq += mymol.topology_->atoms.atom[j+1].q;
                            }
                            ChargeResidual += gmx::square(qq - mymol.chargeQM(qtCM5)[i++]);
                            nChargeResidual++;
                        }
                    }
                }
                ChargeResidual += gmx::square(qtot - mymol.molProp()->getCharge());
                nChargeResidual++;
                increaseEnergy(ermsCHARGE, (ChargeResidual/nChargeResidual));
            }
            if (weight(ermsESP))
            {
                real rrms     = 0;
                real wtot     = 0;
                real cosangle = 0;
                if (nullptr != mymol.shellfc_)
                {
                    mymol.Qgresp_.updateAtomCoords(mymol.x());
                }
                if (bFitZeta_)
                {
                    mymol.Qgresp_.updateZeta(&mymol.topology_->atoms, poldata());
                }
                mymol.Qgresp_.updateAtomCharges(&mymol.topology_->atoms);
                mymol.Qgresp_.calcPot();
                increaseEnergy(ermsESP, convert2gmx(mymol.Qgresp_.getRms(&wtot, &rrms, &cosangle), eg2cHartree_e));
            }
            if (weight(ermsMU))
            {
                mymol.CalcDipole();
                mymol.rotateDipole(mymol.muQM(qtCalc), mymol.muQM(qtElec));
                if (bQM())
                {
                    rvec dmu;
                    rvec_sub(mymol.muQM(qtCalc), mymol.muQM(qtElec), dmu);
                    increaseEnergy(ermsMU, iprod(dmu, dmu));
                }
                else
                {
                    increaseEnergy(ermsMU, gmx::square(mymol.dipQM(qtCalc) - mymol.dipExper()));
                }
            }
            if (weight(ermsQUAD))
            {
                mymol.CalcQuadrupole();
                for (auto mm = 0; mm < DIM; mm++)
                {
                    for (auto nn = 0; nn < DIM; nn++)
                    {
                        if (bFullTensor_ || mm == nn)
                        {
                            increaseEnergy(ermsQUAD, gmx::square(mymol.QQM(qtCalc)[mm][nn] - mymol.QQM(qtElec)[mm][nn]));
                        }
                    }
                }
            }
        }
    }
    sumEnergies();
    normalizeEnergies();
    printEnergies(debug);
}

void OptACM::polData2TuneACM()
{
    param_.clear();
    auto *ic = indexCount();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto ei   = poldata()->findEem(ai->name());
            GMX_RELEASE_ASSERT(ei != poldata()->EndEemprops(), "Cannot find eemprops");
            
            if (bFitChi_)
            {
                auto J00  = ei->getJ0();
                param_.push_back(std::move(J00));
                
                if (ai->name().compare(fixchi()) != 0)
                {
                    auto Chi0 = ei->getChi0();
                    param_.push_back(std::move(Chi0));
                }
            }
            if (bFitZeta_)
            {
                auto nzeta = ei->getNzeta();
                auto zeta  = ei->getZeta(nzeta-1); // We only optimize zeta for shell.
                if (0 != zeta)
                {
                    param_.push_back(std::move(zeta));
                }
                else
                {
                    gmx_fatal(FARGS, "Zeta is zero for atom %s in model %s\n",
                              ai->name().c_str(), getEemtypeName(poldata()->getEqdModel()));
                }
            }
            if (bFitAlpha_)
            {
                auto alpha = 0.0;
                auto sigma = 0.0;
                if (poldata()->getAtypePol(ai->name(), &alpha, &sigma))
                {
                    if (0 != alpha)
                    {
                        param_.push_back(std::move(alpha));
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Polarizability is zero for atom %s\n", ai->name().c_str());
                    }
                }
            }
        }
    }
    if (optHfac())
    {
        param_.push_back(hfac());
    }
}

void OptACM::TuneACM2PolData()
{
    char     zstr[STRLEN];
    char     z_sig[STRLEN];
    char     buf[STRLEN];
    char     buf_sig[STRLEN];

    int      n     = 0;
    double   zeta  = 0;
    double   sigma = 0;

    auto     pd    = poldata();
    auto    *ic    = indexCount();

    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto ei = pd->findEem(ai->name());
            GMX_RELEASE_ASSERT(ei != pd->EndEemprops(), "Cannot find eemprops");

            if (bFitChi_)
            {
                ei->setJ0(param_[n]);
                ei->setJ0_sigma(psigma_[n++]);
                
                if (ai->name().compare(fixchi()) != 0)
                {
                    ei->setChi0(param_[n]);
                    ei->setChi0_sigma(psigma_[n++]);
                }
            }
            if (bFitZeta_)
            {
                zstr[0]            = '\0';
                z_sig[0]           = '\0';
                std::string qstr   = ei->getQstr();
                std::string rowstr = ei->getRowstr();
                auto iModel = poldata()->getEqdModel();
                if (getEemtypeDistributed(iModel))
                {
                    if (bSameZeta_)
                    {
                        zeta   = param_[n]; // Same zeta will be used for both core and shell
                        sigma  = psigma_[n++];

                        sprintf(buf, "%g %g ", zeta, zeta);
                        sprintf(buf_sig, "%g %g ", sigma, sigma);
                        strcat(zstr, buf);
                        strcat(z_sig, buf_sig);
                    }
                    else
                    {
                        for (auto i = 0; i < ei->getNzeta(); i++)
                        {   
                            zeta  = ei->getZeta(i); //zeta for core read from poladata
                            sigma = 0;
                            if (i > 0)
                            {
                                zeta   = param_[n]; //zeta for shell taken from param_ vector
                                sigma  = psigma_[n++];
                            }
                            sprintf(buf, "%g ", zeta);
                            sprintf(buf_sig, "%g ", sigma);
                            strcat(zstr, buf);
                            strcat(z_sig, buf_sig);
                        }
                    }
                }
                else
                {
                    for (auto i = 0; i < ei->getNzeta(); i++)
                    {
                        zeta   = param_[n];
                        sigma  = psigma_[n++];
                        
                        sprintf(buf, "%g ", zeta);
                        sprintf(buf_sig, "%g ", sigma);
                        strcat(zstr, buf);
                        strcat(z_sig, buf_sig);
                    }
                }
                ei->setRowZetaQ(rowstr, zstr, qstr);
                ei->setZetastr(zstr);
                ei->setZeta_sigma(z_sig);
            }
            if (bFitAlpha_)
            {
                std::string ptype;
                if (pd->atypeToPtype(ai->name(), ptype))
                {
                    pd->setPtypePolarizability(ptype, param_[n], psigma_[n]);
                    n++;
                }
                else
                {
                    gmx_fatal(FARGS, "No Ptype for atom type %s\n", ai->name().c_str());
                }
            }
        }
    }
    if (optHfac())
    {
        setHfac(param_[n++]);
    }
}

void OptACM::InitOpt(real  factor)
{
    polData2TuneACM();

    orig_.resize(param_.size(), 0);
    best_.resize(param_.size(), 0);
    lower_.resize(param_.size(), 0);
    upper_.resize(param_.size(), 0);
    psigma_.resize(param_.size(), 0);
    pmean_.resize(param_.size(), 0);

    if (factor < 1)
    {
        factor = 1/factor;
    }
    for (size_t i = 0; i < param_.size(); i++)
    {
        best_[i]  = orig_[i] = param_[i];
        lower_[i] = orig_[i]/factor;
        upper_[i] = orig_[i]*factor;
    }
}

double OptACM::calcPenalty(AtomIndexIterator ai)
{
    double         penalty = 0;
    const auto     pd      = poldata();

    auto           ei      = pd->findEem(ai->name());
    auto           ai_elem = pd->ztype2elem(ei->getName());
    auto           ai_chi  = ei->getChi0();
    auto           ai_J0   = ei->getJ0();
    auto           ai_atn  = gmx_atomprop_atomnumber(atomprop(), ai_elem.c_str());

    if (strlen(fixchi()) != 0)
    {
        const auto ref_eem  = pd->findEem(fixchi());
        if (ai_chi < ref_eem->getChi0())
        {
            penalty += penalty_;
        }
    }
    
    if (ai->name() == "z_c3" && (ai_chi < 5 or ai_chi > 8))
    {
        penalty += (ai_atn * penalty_);
    }

    if (ai->name() == "z_h1" && ai_chi > 2.5)
    {
       penalty += (6 * penalty_);
    }

    auto *ic = indexCount();
    for (auto aj = ic->beginIndex(); aj < ic->endIndex(); ++aj)
    {
        if (!aj->isConst())
        {
            const auto ej      = pd->findEem(aj->name());
            const auto aj_elem = pd->ztype2elem(ej->getName());
            auto       aj_atn  = gmx_atomprop_atomnumber(atomprop(), aj_elem.c_str());

            if (ai_atn != aj_atn)
            {
                //Penalize if HeavyAtoms_chi <= H_chi or HeavyAtoms_J0 <= H_J0                
                auto aj_chi = ej->getChi0();
                auto aj_J0  = ej->getJ0();
                if ((ai_atn == 1 && aj_atn > 1  && (aj_chi <= ai_chi || aj_J0 <= ai_J0)) ||
                    (ai_atn > 1  && aj_atn == 1 && (aj_chi <= ai_chi || aj_J0 <= ai_J0)))
                {
                    penalty += (std::abs((aj_atn - ai_atn)) * penalty_);
                }
            }
        }
    }
    return penalty;
}

double OptACM::objFunction(const double v[])
{
    double bound   = 0;
    double penalty = 0;
    int    n       = 0;

    auto   np = param_.size();
    for (size_t i = 0; i < np; i++)
    {
        param_[i] = v[i];
    }

    TuneACM2PolData();
    auto *ic = indexCount();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto name = ai->name();
            
            if (bFitChi_)
            {
                auto J00  = param_[n++];
                bound    += l2_regularizer(J00, J0Min(), J0Max());
                if (strcasecmp(name.c_str(), fixchi()) != 0)
                {
                    auto Chi0 = param_[n++];
                    bound    += l2_regularizer(Chi0, chi0Min(), chi0Max());
                }
                if (penalize())
                {
                    penalty += calcPenalty(ai);
                }
            }
            if (bFitZeta_)
            {
                auto nzeta = poldata()->getNzeta(ai->name());
                for (auto zz = 0; zz < (nzeta-1); zz++)
                {
                    auto zeta = param_[n++];
                    bound += l2_regularizer(zeta, zetaMin(), zetaMax());
                }
            }
        }
    }
    if (optHfac())
    {
        setHfac(param_[n++]);
        bound += 100*gmx::square(hfacDiff());
    }
    calcDeviation();
    increaseEnergy(ermsBOUNDS, bound);
    increaseEnergy(ermsTOT, bound);
    if (penalize())
    {
        increaseEnergy(ermsTOT, penalty);
    }
    return energy(ermsTOT);
}

bool OptACM::optRun(FILE                   *fp,
                    FILE                   *fplog,
                    int                     nrun,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot)
{
    std::vector<double> optb, opts, optm;
    double              chi2, chi2_min;
    bool                bMinimum = false;

    auto                func = [&] (const double v[]) {
            return objFunction(v);
        };

    if (MASTER(commrec()))
    {
        if (PAR(commrec()))
        {
            for (int dest = 1; dest < commrec()->nnodes; dest++)
            {
                gmx_send_int(commrec(), dest, (nrun*TuneACM_.maxIter()*param_.size()));
            }
        }
        chi2 = chi2_min = GMX_REAL_MAX;
        TuneACM_.setFunc(func, param_, lower_, upper_, &chi2);
        TuneACM_.Init(xvgconv, xvgepot, oenv);

        for (auto n = 0; n < nrun; n++)
        {
            if ((nullptr != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n, nrun);
            }
            TuneACM_.MCMC();
            TuneACM_.getBestParam(optb);
            TuneACM_.getPsigma(opts);
            TuneACM_.getPmean(optm);
            if (chi2 < chi2_min)
            {
                bMinimum = true;
                for (size_t k = 0; k < param_.size(); k++)
                {
                    best_[k]   = optb[k];
                    pmean_[k]  = optm[k];
                    psigma_[k] = opts[k];
                }
                chi2_min = chi2;
            }
            TuneACM_.setParam(best_);
        }
        if (bMinimum)
        {
            param_    = best_;
            auto emin = objFunction(best_.data());
            if (fplog)
            {
                fprintf(fplog, "\nMinimum rmsd value during optimization: %.3f.\n", sqrt(emin));
                fprintf(fplog, "Statistics of parameters after optimization\n");
                for (size_t k = 0; k < param_.size(); k++)
                {
                    fprintf(fplog, "Parameter %3zu  Best value:%10g  Mean value:%10g  Sigma:%10g\n",
                            k, best_[k], pmean_[k], psigma_[k]);
                }
            }
        }
    }
    else
    {
        /* S L A V E   N O D E S */
        auto niter = gmx_recv_int(commrec(), 0);
        for (auto n = 0; n < niter + 2; n++)
        {
            calcDeviation();
        }
    }
    setFinal();
    if (MASTER(commrec()))
    {
        chi2 = objFunction(best_.data());
        printEnergies(fp);
        printEnergies(fplog);
    }
    return bMinimum;
}
}

int alex_tune_eem(int argc, char *argv[])
{
    static const char          *desc[] = {
        "tune_eem read a series of molecules and corresponding experimental",
        "dipole moments from a file, and tunes parameters in an algorithm",
        "until the experimental dipole moments are reproduced by the",
        "charge generating algorithm AX as implemented in the gentop program.[PAR]",
        "Minima and maxima for the parameters can be set, these are however",
        "not strictly enforced, but rather they are penalized with a harmonic",
        "function, for which the force constant can be set explicitly.[PAR]",
        "At every reinit step parameters are changed by a random amount within",
        "the fraction set by step size, and within the boundaries given",
        "by the minima and maxima. If the [TT]-random[tt] flag is",
        "given a completely random set of parameters is generated at the start",
        "of each run. At reinit steps however, the parameters are only changed",
        "slightly, in order to speed-up local search but not global search."
        "In other words, complete random starts are done only at the beginning of each",
        "run, and only when explicitly requested.[PAR]",
        "The absolut dipole moment of a molecule remains unchanged if all the",
        "atoms swap the sign of the charge. To prevent this kind of mirror",
        "effects a penalty is added to the square deviation ",
        "if hydrogen atoms have a negative charge. Similarly a penalty is",
        "added if atoms from row VI or VII in the periodic table have a positive",
        "charge. The penalty is equal to the force constant given on the command line",
        "time the square of the charge.[PAR]",
        "One of the electronegativities (chi) is redundant in the optimization,",
        "only the relative values are meaningful.",
        "Therefore by default we fix the value for hydrogen to what is written",
        "in the eemprops.dat file (or whatever is given with the [tt]-d[TT] flag).",
        "A suitable value would be 2.3, the original, value due to Pauling,",
        "this can by overridden by setting the [tt]-fixchi[TT] flag to something else (e.g. a non-existing atom).[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };

    t_filenm                    fnm[] = {
        { efDAT, "-f",         "allmols",       ffREAD  },
        { efDAT, "-d",         "gentop",        ffOPTRD },
        { efDAT, "-o",         "tunedip",       ffWRITE },
        { efDAT, "-sel",       "molselect",     ffREAD  },
        { efXVG, "-table",     "table",         ffOPTRD },
        { efLOG, "-g",         "charges",       ffWRITE },
        { efXVG, "-qhisto",    "q_histo",       ffWRITE },
        { efXVG, "-dipcorr",   "dip_corr",      ffWRITE },
        { efXVG, "-mucorr",    "mu_corr",       ffWRITE },
        { efXVG, "-thetacorr", "theta_corr",    ffWRITE },
        { efXVG, "-espcorr",   "esp_corr",      ffWRITE },
        { efXVG, "-alphacorr", "alpha_corr",    ffWRITE },
        { efXVG, "-qcorr",     "q_corr",        ffWRITE },
        { efXVG, "-isopol",    "isopol_corr",   ffWRITE },
        { efXVG, "-anisopol",  "anisopol_corr", ffWRITE },
        { efXVG, "-conv",      "param-conv",    ffWRITE },
        { efXVG, "-epot",      "param-epot",    ffWRITE },
        { efTEX, "-latex",     "eemprop",       ffWRITE }
    };

    const int                   NFILE         = asize(fnm);

    int                         nrun          = 1;
    int                         reinit        = 0;
    real                        th_toler      = 170;
    real                        ph_toler      = 5;
    real                        dip_toler     = 0.5;
    real                        quad_toler    = 5;
    real                        alpha_toler   = 3;
    real                        isopol_toler  = 2;
    real                        factor        = 0.8;
    real                        efield        = 10;
    char                       *opt_elem      = nullptr;
    gmx_bool                    bRandom       = false;
    gmx_bool                    bcompress     = false;
    gmx_bool                    bPrintTable   = false;
    gmx_bool                    bZero         = true;
    gmx_bool                    bOptimize     = true;
    gmx_bool                    bForceOutput  = true;

    t_pargs                     pa[]         = {
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-nrun",   FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of atom types to optimize, e.g. \"H C Br\". The other available atom types in gentop.dat are left unmodified. If this variable is not set, all elements will be optimized." },
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
        { "-zero", FALSE, etBOOL, {&bZero},
          "Use molecules with zero dipole in the fit as well" },
        { "-dip_toler", FALSE, etREAL, {&dip_toler},
          "Tolerance (Debye) for marking dipole as an outlier in the log file" },
        { "-quad_toler", FALSE, etREAL, {&quad_toler},
          "Tolerance (Buckingham) for marking quadrupole as an outlier in the log file" },
        { "-alpha_toler", FALSE, etREAL, {&alpha_toler},
          "Tolerance (A^3) for marking diagonal elements of the polarizability tensor as an outlier in the log file" },
        { "-isopol_toler", FALSE, etREAL, {&isopol_toler},
          "Tolerance (A^3) for marking isotropic polarizability as an outlier in the log file" },
        { "-th_toler", FALSE, etREAL, {&th_toler},
          "Minimum angle to be considered a linear A-B-C bond" },
        { "-ph_toler", FALSE, etREAL, {&ph_toler},
          "Maximum angle to be considered a planar A-B-C/B-C-D torsion" },
        { "-compress", FALSE, etBOOL, {&bcompress},
          "Compress output XML file" },
        { "-btex", FALSE, etBOOL, {&bPrintTable},
          "[HIDDEN]Print the latex table for the Gaussian and Slater exponents" },
        { "-factor", FALSE, etREAL, {&factor},
          "Factor for generating random parameters. Parameters will be taken within the limit factor*x - x/factor" },
        { "-efield",  FALSE, etREAL, {&efield},
          "The magnitude of the external electeric field to calculate polarizability tensor." },
        { "-optimize",     FALSE, etBOOL, {&bOptimize},
          "Do parameter optimization when true, or a single calculation otherwise." },
        { "-force_output", FALSE, etBOOL, {&bForceOutput},
          "Write output even if no new minimum is found" }
          
    };

    FILE                       *fp;
    gmx_output_env_t           *oenv;
    time_t                      my_t;
    MolSelect                   gms;

    std::vector<t_pargs>        pargs;
    for (size_t i = 0; i < asize(pa); i++)
    {
        pargs.push_back(pa[i]);
    }
    alexandria::OptACM opt;
    opt.add_pargs(&pargs);

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm,
                           pargs.size(), pargs.data(),
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    opt.optionsFinished();

    if (MASTER(opt.commrec()))
    {
        fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

        time(&my_t);
        fprintf(fp, "# This file was created %s", ctime(&my_t));
        fprintf(fp, "# alexandria is part of GROMACS:\n#\n");
        fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());
    }
    else
    {
        fp = nullptr;
    }
    if (MASTER(opt.commrec()))
    {
        gms.read(opt2fn_null("-sel", NFILE, fnm));
    }

    const char *tabfn = opt2fn_null("-table", NFILE, fnm);

    opt.Read(fp ? fp : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             bZero,
             opt_elem,
             gms,
             true,
             false,
             false,
             false,
             opt.fitZeta(),
             tabfn);

    if (nullptr != fp)
    {
        fprintf(fp, "In the total data set of %zu molecules we have:\n",
                opt.mymols().size());
    }
    
    if (opt.iChargeGenerationAlgorithm() != eqgESP)
    {
        opt.initQgresp();
    }
    if (opt.iChargeGenerationAlgorithm() != eqgACM)
    {
        opt.initQgacm();
    }

    bool bMinimum = false;
    if (bOptimize)
    {
        if (MASTER(opt.commrec()))
        {
            opt.InitOpt(factor);
        }
        
        bMinimum = opt.optRun(MASTER(opt.commrec()) ? stderr : nullptr,
                              fp,
                              nrun,
                              oenv,
                              opt2fn("-conv", NFILE, fnm),
                              opt2fn("-epot", NFILE, fnm));                   
    }

    if (MASTER(opt.commrec()))
    {
        if (bMinimum || bForceOutput)
        {
            auto iModel = opt.poldata()->getEqdModel();
            gmx_bool bPolar = getEemtypePolarizable(iModel);
            
            auto *ic = opt.indexCount();
            print_electric_props(fp,
                                 opt.mymols(),
                                 opt.poldata(),
                                 opt.mdlog(),
                                 opt.atomprop(),
                                 eqgACM,
                                 opt.watoms(),
                                 opt.hfac(),
                                 opt.lot(),
                                 tabfn,
                                 opt.hwinfo(),
                                 500, //opt.qcycle(),
                                 opt.maxPot(),
                                 1e-6, //opt.qtol(),                          
                                 opt2fn("-qhisto",    NFILE, fnm),
                                 opt2fn("-dipcorr",   NFILE, fnm),
                                 opt2fn("-mucorr",    NFILE, fnm),
                                 opt2fn("-thetacorr", NFILE, fnm),
                                 opt2fn("-espcorr",   NFILE, fnm),
                                 opt2fn("-alphacorr", NFILE, fnm),
                                 opt2fn("-isopol",    NFILE, fnm),
                                 opt2fn("-anisopol",  NFILE, fnm),
                                 opt2fn("-qcorr",     NFILE, fnm),
                                 dip_toler,
                                 quad_toler,
                                 alpha_toler,
                                 isopol_toler,
                                 oenv,
                                 bPolar,
                                 opt.dipole(),
                                 opt.quadrupole(),
                                 opt.fullTensor(),
                                 ic,
                                 opt.commrec(),
                                 efield);
            writePoldata(opt2fn("-o", NFILE, fnm), opt.poldata(), bcompress);
            if (bPrintTable)
            {
                FILE        *tp;
                tp = gmx_ffopen(opt2fn("-latex", NFILE, fnm), "w");
                alexandria_poldata_eemprops_table(tp, opt.poldata());
                gmx_ffclose(tp);
            }
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
        gmx_ffclose(fp);
    }
    return 0;
}
