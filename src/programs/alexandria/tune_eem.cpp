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

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <random>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "alex_modules.h"
#include "gentop_core.h"
#include "gmx_simple_comm.h"
#include "molgen.h"
#include "molprop_util.h"
#include "mymol_low.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_tables.h"
#include "poldata_xml.h"
#include "tuning_utility.h"
#include "units.h"

namespace alexandria
{

static void my_fclose(FILE *fp)
{
    int myerrno = gmx_ffclose(fp);
    if (myerrno != 0)
    {
        fprintf(stderr, "Error %d closing file\n", myerrno);
    }
}

class OptACM : public MolGen, Bayes
{
    using param_type = std::vector<double>;

    private:
        bool       bFullTensor_                  = false;
        bool       bSameZeta_                    = true;
        bool       bUseCM5_                      = false;
        bool       bRemoveMol_                   = true;
        bool       bPointCore_                   = false;
        gmx::unique_cptr<FILE, my_fclose> fplog_ = nullptr;

    public:

        OptACM() {}

        ~OptACM() {}

        bool bESP() const { return weight(ermsESP); }

        bool dipole() const { return weight(ermsMU); }

        bool quadrupole() const { return weight(ermsQUAD); }

        bool fullTensor() const { return bFullTensor_; }

        bool sameZeta() const { return bSameZeta_; }

        bool useCM5() const {return bUseCM5_; }
        
        bool removeMol() const {return bRemoveMol_; }
        
        void set_pointCore(bool pointCore) {bPointCore_ =  pointCore;}

        void add_pargs(std::vector<t_pargs> *pargs)
        {
            t_pargs pa[] =
            {
                { "-fullTensor", FALSE, etBOOL, {&bFullTensor_},
                  "Consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization" },
                { "-samezeta", FALSE, etBOOL, {&bSameZeta_},
                  "Use the same zeta for both the core and the shell of the Drude model." },
                { "-cm5", FALSE, etBOOL, {&bUseCM5_},
                  "Reproduce CM5 charges in fitting." },
                { "-removemol", FALSE, etBOOL, {&bRemoveMol_},
                  "Remove a molecule from training set if shell minimzation does not converge." },
            };
            for (size_t i = 0; i < asize(pa); i++)
            {
                pargs->push_back(pa[i]);
            }
            addOptions(pargs, etuneEEM);
            Bayes::add_pargs(pargs);
        }

        void optionsFinished()
        {
            MolGen::optionsFinished();
            setBoxConstraint(bConstrain());
        }

        void openLogFile(const char *logfileName)
        {
            fplog_.reset(gmx_ffopen(logfileName, "w"));
        }
        
        FILE *logFile()
        {
            if (fplog_)
            { 
                return fplog_.get();
            }
            else
            {
                return nullptr;
            }
        }
        
        double l2_regularizer (double             x,
                               double             min,
                               double             max,
                               const std::string &label,
                               bool               verbose)
        {
            double p = 0;
            if (x < min)
            {
                p = (0.5 * gmx::square(x-min));
            }
            else if (x > max)
            {
                p = (0.5 * gmx::square(x-max));
            }
            if (verbose && p != 0.0)
            {
                fprintf(logFile(), "Variable %s is %g, should be within %g and %g\n", label.c_str(), x, min, max);
            }
            return p;
        }
        
        void initChargeGeneration();

        /*! \brief
         *
         * Fill parameter vector based on Poldata.
         * \param[in] bRandom Generate random initial values for parameters if true
         */
        void InitOpt(bool bRandom);

        /*! \brief
         * Copy the optimization parameters to the poldata structure
         * \param[in] changed List over the parameters that have changed.
         */
        virtual void toPoldata(const std::vector<bool> &changed);

        virtual double calcDeviation();

        /*! \brief
         * Do the actual optimization.
         * \param[in] fp     FILE pointer for logging
         * \param[in] oenv   Output environment for managing xvg files etc.
         * \param[in] xvgconv Output file monitoring parameters
         * \param[in] xvgepot Output file monitoring penalty function
         * \return true if better parameters were found.
         */
        bool optRun(FILE                   *fp,
                    int                     nrun,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot);
};

void OptACM::initChargeGeneration()
{
    std::string method, basis, conf, type, myref, mylot;
    splitLot(lot(), &method, &basis);
    tensor           polar      = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    rvec             vec;
    for (auto &mymol : mymols())
    {
        if (fit("alpha"))
        {
            // For fitting alpha we need a reference polarizability
            double ref_pol, error, T = 0;
            if (mymol.getPropRef(MPO_POLARIZABILITY, iqmQM,
                                 method, basis, "",
                                 (char *)"electronic",
                                 &ref_pol, &error, &T,
                                 &myref, &mylot, vec, polar))
            {
                mymol.SetElectronicPolarizability(ref_pol);
            }
            else
            {
                if (logFile())
                {
                    fprintf(logFile(), "Removing %s due to lacking reference polarizability at the %s/%s LoT.\n",
                            mymol.getMolname().c_str(),
                            method.c_str(), basis.c_str());
                }
                mymol.eSupp_ = eSupport::No;
            }
        }
        if (mymol.eSupp_ != eSupport::No)
        {
            mymol.QgenAcm_ = new QgenAcm(poldata(), mymol.atoms_, 
                                         mymol.getCharge());
        }
    }
}

static void dumpQX(FILE *fp, MyMol *mol, const std::string &info)
{
    if (false && fp)
    {
        std::string label = mol->getMolname() + "-" + info;
        fprintf(fp, "%s q:", label.c_str());
        t_mdatoms *md = mol->getMdatoms();
        for (int i = 0; i < mol->atoms_->nr; i++)
        {
            fprintf(fp, " %g (%g)", mol->atoms_->atom[i].q,
                    md->chargeA[i]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "%s alpha", label.c_str());
        int ft = F_POLARIZATION;
        for(int i = 0; i < mol->ltop_->idef.il[ft].nr; i += interaction_function[ft].nratoms+1)
        {
            auto tp = mol->ltop_->idef.il[ft].iatoms[i];
            fprintf(fp, " %g", mol->ltop_->idef.iparams[tp].polarize.alpha);
        }
        fprintf(fp, "\n");
        pr_rvecs(fp, 0, label.c_str(), mol->x().rvec_array(), mol->atoms_->nr);
    }
}

double OptACM::calcDeviation()
{
    bool  verbose = final();
    resetEnergies();
    if (MASTER(commrec()))
    {
        if (weight(ermsBOUNDS))
        {
            const param_type &param = Bayes::getParam();
            double            bound = 0;
            size_t            n     = 0;
            for (auto &optIndex : optIndex_)
            {
                auto p = poldata()->findForcesConst(optIndex.iType()).findParameterTypeConst(optIndex.id(), optIndex.type());

                bound  += l2_regularizer(param[n++], p.minimum(), p.maximum(),
                                         optIndex.name(), verbose);
            }
            increaseEnergy(ermsBOUNDS, bound);
            GMX_RELEASE_ASSERT(n == param.size(), 
                               gmx::formatString("Death horror error. n=%zu param.size()=%zu", n, param.size()).c_str());
        }
    }

    if (PAR(commrec()))
    {
        bool bFinal = final();
        gmx_bcast(sizeof(final()), &bFinal, commrec());
        if (bFinal)
        {
            setFinal();
        }
        else
        {
            // Communicate the force field data.
            // TODO: Optimize to only do relevant stuff
            poldata()->broadcast_eemprop(commrec());
        }
    }
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupport::Local) ||
            (final() && (mymol.eSupp_ == eSupport::Remote)))
        {
            dumpQX(logFile(), &mymol, "BEFORE");
            std::vector<double> qq;
            auto q = mymol.QgenAcm_->q();

            qq.resize(q.size(), 0);
            for (size_t i = 0; i < q.size(); i++)
            {
                qq[i] = q[i];
            }
            // Update the polarizabilities only once before the loop
            if (fit("alpha"))
            {
                mymol.UpdateIdef(poldata(), InteractionType::POLARIZATION);
            }
            bool converged = false;
            int  iter      = 0;
            do
            {
                // Update charges in mtop before doing
                // shell optimization.
                for (int i = 0; i < mymol.mtop_->natoms; i++)
                {
                    mymol.mtop_->moltype[0].atoms.atom[i].q      =
                        mymol.mtop_->moltype[0].atoms.atom[i].qB =
                            mymol.atoms_->atom[i].q;
                }
                mymol.zeta2atoms(poldata());
                dumpQX(logFile(), &mymol, "LOOP1");
                if (nullptr != mymol.shellfc_)
                {
                    double rmsf;
                    auto imm = mymol.computeForces(logFile(), commrec(), &rmsf);
                    if (immStatus::OK != imm)
                    {
                        if (logFile())
                        {
                            fprintf(logFile(), "Could not compute forces for %s.\n",
                                    mymol.getMolname().c_str());
                        }
                        if (removeMol())
                        {
                            mymol.eSupp_ = eSupport::No;
                            fprintf(logFile()," Removing it from the data set.\n");
                            break;
                        }
                    }
                }
                dumpQX(logFile(), &mymol, "LOOP2");
                auto qgen =  mymol.QgenAcm_->generateCharges(debug,
                                                             mymol.getMolname(),
                                                             poldata(),
                                                             mymol.atoms_,
                                                             mymol.x(),
                                                             mymol.bondsConst());
                if (qgen != eQgen::OK)
                {
                    gmx_fatal(FARGS, "Could not generate charges for %s: %s",
                              mymol.getMolname().c_str(),
                              mymol.QgenAcm_->message());
                }
                double EemRms = 0;
                for (size_t i = 0; i < q.size(); i++)
                {
                    EemRms   += gmx::square(qq[i] - q[i]);
                    qq[i]     = q[i];
                }
                EemRms   /= q.size();
                converged = (EemRms < qtol()) || (nullptr == mymol.shellfc_);
                iter++;
            }
            while ((!converged) && (iter < qcycle()));
            if (!converged)
            {
                if (logFile())
                {
                    fprintf(logFile(), "Could not generate charges for %s. Removing compound.\n",
                            mymol.getMolname().c_str());
                }
                if (removeMol())
                {
                    mymol.eSupp_ = eSupport::No;
                }
            }
            // Check whether we have disabled this compound
            if (mymol.eSupp_ == eSupport::No)
            {
                continue;
            }
            for (int i = 0; i < mymol.mtop_->natoms; i++)
            {
                mymol.mtop_->moltype[0].atoms.atom[i].q      =
                    mymol.mtop_->moltype[0].atoms.atom[i].qB =
                    mymol.atoms_->atom[i].q;
            }
            dumpQX(logFile(), &mymol, "AFTERLOOP");
            if (weight(ermsCHARGE))
            {
                double qtot = 0;
                int    i;
                for (int j = i = 0; j < mymol.atoms_->nr; j++)
                {
                    auto atype = poldata()->findParticleType(*mymol.atoms_->atomtype[j]);
                    auto qparm = atype->parameterConst("charge");
                    double qj  = mymol.atoms_->atom[j].q;
                    qtot += qj;
                    switch (qparm.mutability())
                    {
                    case Mutability::Fixed:
                        if (qparm.value() != qj)
                        {
                            GMX_THROW(gmx::InternalError(gmx::formatString("Fixed charge for atom %s in %s was changed from %g to %g",
                                                                           *mymol.atoms_->atomname[j], mymol.getMolname().c_str(), qparm.value(), qj).c_str()));
                        }
                        break;
                    case Mutability::Bounded:
                        {
                            real dq = 0;
                            if (qj < qparm.minimum())
                            {
                                dq = qparm.minimum() - qj;
                            }
                            else if (qj > qparm.maximum())
                            {
                                dq = qj - qparm.maximum();
                            } 
                            increaseEnergy(ermsCHARGE, dq*dq);
                        }
                        break;
                    default:
                        break;
                    }
                }
                increaseEnergy(ermsCHARGE, gmx::square(qtot - mymol.getCharge()));
            }
            if (weight(ermsESP))
            {
                real rrms     = 0;
                real cosangle = 0;
                if (nullptr != mymol.shellfc_)
                {
                    mymol.QgenResp_->updateAtomCoords(mymol.x());
                }
                if (fit("zeta"))
                {
                    mymol.QgenResp_->updateZeta(mymol.atoms_, poldata());
                }
                dumpQX(logFile(), &mymol, "ESP");
                mymol.QgenResp_->updateAtomCharges(mymol.atoms_);
                mymol.QgenResp_->calcPot(poldata()->getEpsilonR());
                auto myRms =
                    convertToGromacs(mymol.QgenResp_->getRms(&rrms, &cosangle),
                                     "Hartree/e");
                increaseEnergy(ermsESP, gmx::square(myRms));
                if (debug)
                {
                    fprintf(debug, "%s ESPrms = %g cosangle = %g\n",
                            mymol.getMolname().c_str(),
                            myRms, cosangle);
                }
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
                for (int mm = 0; mm < DIM; mm++)
                {
                    for (int nn = 0; nn < DIM; nn++)
                    {
                        if (bFullTensor_ || mm == nn)
                        {
                            increaseEnergy(ermsQUAD, gmx::square(mymol.QQM(qtCalc)[mm][nn] - mymol.QQM(qtElec)[mm][nn]));
                        }
                    }
                }
            }
            if (weight(ermsPolar))
            {
                double diff2 = 0;
                mymol.CalcPolarizability(10, commrec(), nullptr);
                if (bFullTensor_)
                {
                    // It is already squared
                    diff2 = mymol.PolarizabilityTensorDeviation();
                }
                else
                {
                    diff2 = gmx::square(mymol.PolarizabilityDeviation());
                }
                if (false && logFile())
                {
                    fprintf(logFile(), "DIFF %s %g\n", mymol.getMolname().c_str(), diff2);
                }
                increaseEnergy(ermsPolar, diff2);
            }
        }
    }
    sumEnergies();
    if (verbose && logFile())
    {
        printParameters(logFile());
        printEnergies(logFile());
    }
    return energy(ermsTOT)/nMolSupport();
}

void OptACM::InitOpt(bool bRandom)
{
    for(auto &optIndex : optIndex_)
    {
        auto param = poldata()->findForcesConst(optIndex.iType()).findParameterTypeConst(optIndex.id(), optIndex.type());
        if (param.ntrain() >= static_cast<uint64_t>(mindata()))
        {
            Bayes::addParam(optIndex.name(),
                            param.value(), param.minimum(), param.maximum(),
                            bRandom);
        }
    }
}

void OptACM::toPoldata(const std::vector<bool> &changed)
{
    size_t   n      = 0;
    auto     param  = Bayes::getParam();
    auto     psigma = Bayes::getPsigma();
    if (psigma.empty())
    {
        psigma.resize(param.size(), 0);
    }
    Bayes::printParameters(debug);
    for (const auto &optIndex : optIndex_)
    {
        if (changed[n])
        {
            auto p = poldata()->findForces(optIndex.iType())->findParameterType(optIndex.id(), optIndex.type());
            
            p->setValue(param[n]);
            p->setUncertainty(psigma[n]);
        }
        n++;
    }
    GMX_RELEASE_ASSERT(n == changed.size(),
                       gmx::formatString("n = %zu changed.size() = %zu",
                                         n, changed.size()).c_str());
}

bool OptACM::optRun(FILE                   *fp,
                    int                     nrun,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot)
{
    bool bMinimum = false;
    if (MASTER(commrec()))
    {
        if (PAR(commrec()))
        {
            int niter = 2 + nrun*Bayes::maxIter()*Bayes::nParam();
            for (int dest = 1; dest < commrec()->nnodes; dest++)
            {
                gmx_send_int(commrec(), dest, niter);
            }
        }
        std::vector<std::string> paramClass;
        for(const auto &fm : typesToFit())
        {
            paramClass.push_back(fm.first);
        }
        Bayes::setOutputFiles(xvgconv, paramClass, xvgepot, oenv);
        double     chi2_min = calcDeviation();
        fprintf(logFile(), "Initial chi2 value %g\n", chi2_min);
        printEnergies(logFile());
        for (auto n = 0; n < nrun; n++)
        {
            if ((nullptr != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n, nrun);
            }
            double chi2 = 0;
            if (Bayes::adaptive())
            {
                chi2 = Bayes::Adaptive_MCMC(logFile());
            }
            else
            {
                chi2 = Bayes::MCMC(logFile());
            }
            if (chi2 < chi2_min)
            {
                bMinimum = true;
                chi2_min = chi2;
            }
        }
        if (bMinimum)
        {
            auto i_param        = Bayes::getInitialParam();
            auto best           = Bayes::getBestParam();
            auto pmean          = Bayes::getPmean();
            auto psigma         = Bayes::getPsigma();
            auto attemptedMoves = Bayes::getAttemptedMoves();
            auto acceptedMoves  = Bayes::getAcceptedMoves();
            auto paramNames     = Bayes::getParamNames();
            if (logFile())
            {
                fprintf(logFile(), "\nMinimum RMSD value during optimization: %.3f.\n", sqrt(chi2_min));
                fprintf(logFile(), "Statistics of parameters after optimization\n");
                fprintf(logFile(), "#best %zu #mean %zu #sigma %zu #param %zu\n",
                        best.size(), pmean.size(), psigma.size(), paramNames.size());
                if (best.size() == Bayes::nParam())
                {
                    for (size_t k = 0; k < Bayes::nParam(); k++)
                    {
                        double acceptance_ratio = 100*(double(acceptedMoves[k])/attemptedMoves[k]);
                        fprintf(logFile(), "%-10s  Initial:%10g  Best:%10g  Mean:%10g  Sigma:%10g  Attempted moves:%3d  Acceptance ratio:%5g\n",
                                paramNames[k].c_str(), i_param[k], best[k], pmean[k], psigma[k], attemptedMoves[k], acceptance_ratio);
                    }
                }
            }
        }
    }
    else
    {
        /* S L A V E   N O D E S */
        auto niter = gmx_recv_int(commrec(), 0);
        for (auto n = 0; n < niter; n++)
        {
            (void) calcDeviation();
        }
        if (debug)
        {
            fprintf(debug, "Ready doing calculations on slave nodes\n");
        }
    }
    setFinal();
    if (MASTER(commrec()))
    {
        auto best = Bayes::getBestParam();
        if (!best.empty())
        {
            // Restore best parameter set
            Bayes::setParam(best);
            // Copy it to Poldata
            std::vector<bool> changed;
            changed.resize(best.size(), true);
            toPoldata(changed);
            // Compute the deviation once more
            double chi2 = calcDeviation();
            printEnergies(fp);
            printEnergies(logFile());
            fprintf(logFile(), "Minimum energy chi2 %g\n", chi2);
        }
        else
        {
            fprintf(logFile(), "Did not find a better parameter set\n");
        }
    }
    return bMinimum;
}

} // namespace alexandria

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
        { efDAT, "-o",         "tune_eem",      ffWRITE },
        { efDAT, "-sel",       "molselect",     ffREAD  },
        { efXVG, "-table",     "table",         ffOPTRD },
        { efLOG, "-g",         "tune_eem",      ffWRITE },
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
    real                        esp_toler     = 30;
    real                        dip_toler     = 0.5;
    real                        quad_toler    = 5;
    real                        alpha_toler   = 3;
    real                        isopol_toler  = 2;
    real                        efield        = 10;
    bool                        bRandom       = false;
    bool                        bcompress     = false;
    bool                        bZero         = true;
    bool                        bOptimize     = true;
    bool                        bForceOutput  = true;
    bool                        useOffset     = false;
    
    static const char          *select_types[]   = {nullptr, "Train", "Test", "Ignore", "Unknown", nullptr};

    t_pargs                     pa[]         = {
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-nrun",   FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
        { "-zero", FALSE, etBOOL, {&bZero},
          "Use molecules with zero dipole in the fit as well" },
        { "-esp_toler", FALSE, etREAL, {&esp_toler},
          "Tolerance (kJ/mol e) for marking ESP as an outlier in the log file" },
        { "-dip_toler", FALSE, etREAL, {&dip_toler},
          "Tolerance (Debye) for marking dipole as an outlier in the log file" },
        { "-quad_toler", FALSE, etREAL, {&quad_toler},
          "Tolerance (Buckingham) for marking quadrupole as an outlier in the log file" },
        { "-alpha_toler", FALSE, etREAL, {&alpha_toler},
          "Tolerance (A^3) for marking diagonal elements of the polarizability tensor as an outlier in the log file" },
        { "-isopol_toler", FALSE, etREAL, {&isopol_toler},
          "Tolerance (A^3) for marking isotropic polarizability as an outlier in the log file" },
        { "-use_offset", FALSE, etBOOL,{&useOffset},
          "Fit regression analysis of results to y = ax+b instead of y = ax" },
        { "-compress", FALSE, etBOOL, {&bcompress},
          "Compress output XML file" },
        { "-efield",  FALSE, etREAL, {&efield},
          "The magnitude of the external electeric field to calculate polarizability tensor." },
        { "-optimize",     FALSE, etBOOL, {&bOptimize},
          "Do parameter optimization when true, or a single calculation otherwise." },
        { "-force_output", FALSE, etBOOL, {&bForceOutput},
          "Write output even if no new minimum is found" },
        { "-select", FALSE, etENUM, {select_types},
          "Select type for making the dataset for training or testing." }

    };

    gmx_output_env_t           *oenv;
    MolSelect                   gms;

    std::vector<t_pargs>        pargs;
    for (size_t i = 0; i < asize(pa); i++)
    {
        pargs.push_back(pa[i]);
    }
    alexandria::OptACM opt;
    opt.add_pargs(&pargs);

    if (!parse_common_args(&argc, 
                           argv, 
                           PCA_CAN_VIEW, 
                           NFILE, 
                           fnm,
                           pargs.size(), 
                           pargs.data(),
                           asize(desc), 
                           desc, 
                           0, 
                           nullptr, 
                           &oenv))
    {
        return 0;
    }
    
    opt.optionsFinished();
    
    if (MASTER(opt.commrec()))
    {
        opt.openLogFile(opt2fn("-g", NFILE, fnm));
        print_header(opt.logFile(), pargs);
    }
    else if (false)
    {
        std::string logf = gmx::formatString("log%d.log",
                                             opt.commrec()->nodeid);
        opt.openLogFile(logf.c_str());
    }
    
    if (MASTER(opt.commrec()))
    {
        gms.read(opt2fn_null("-sel", NFILE, fnm));
    }

    const char *tabfn = opt2fn_null("-table", NFILE, fnm);   
    
    iMolSelect select_type = name2molselect(select_types[0]);
    opt.Read(opt.logFile() ? opt.logFile() : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             bZero,
             gms,
             false,
             false,
             tabfn,
             select_type);
             
    bool  pointCore = opt.poldata()->corePointCharge();
    
    opt.set_pointCore(pointCore);
    opt.initChargeGeneration();

    bool bMinimum = false;
    if (bOptimize)
    {
        if (MASTER(opt.commrec()))
        {
            opt.InitOpt(bRandom);
        }

        bMinimum = opt.optRun(MASTER(opt.commrec()) ? stderr : nullptr,
                              nrun,
                              oenv,
                              opt2fn("-conv", NFILE, fnm),
                              opt2fn("-epot", NFILE, fnm));
    }

    if (MASTER(opt.commrec()))
    {
        if (bMinimum || bForceOutput)
        {
            bool  bPolar = opt.poldata()->polarizable();
            print_electric_props(opt.logFile(),
                                 opt.mymols(),
                                 opt.poldata(),
                                 opt.mdlog(),
                                 opt.lot(),
                                 tabfn,
                                 opt.hwinfo(),
                                 opt.qcycle(),
                                 opt.qtol(),
                                 opt2fn("-qhisto",    NFILE, fnm),
                                 opt2fn("-dipcorr",   NFILE, fnm),
                                 opt2fn("-mucorr",    NFILE, fnm),
                                 opt2fn("-thetacorr", NFILE, fnm),
                                 opt2fn("-espcorr",   NFILE, fnm),
                                 opt2fn("-alphacorr", NFILE, fnm),
                                 opt2fn("-isopol",    NFILE, fnm),
                                 opt2fn("-anisopol",  NFILE, fnm),
                                 opt2fn("-qcorr",     NFILE, fnm),
                                 esp_toler,
                                 dip_toler,
                                 quad_toler,
                                 alpha_toler,
                                 isopol_toler,
                                 oenv,
                                 bPolar,
                                 opt.dipole(),
                                 opt.quadrupole(),
                                 opt.fullTensor(),
                                 opt.commrec(),
                                 efield,
                                 useOffset,
                                 opt.optIndex_);
            writePoldata(opt2fn("-o", NFILE, fnm), opt.poldata(), bcompress);            
            if (opt2bSet("-latex", NFILE, fnm))
            {
                FILE        *tp;
                tp = gmx_ffopen(opt2fn("-latex", NFILE, fnm), "w");
                alexandria_eemprops_table(tp, opt.poldata());
                alexandria_eemprops_corr(opt.poldata(), opt.logFile());
                gmx_ffclose(tp);
            }
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
    }
    return 0;
}
