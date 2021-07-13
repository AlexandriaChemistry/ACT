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
#include "gromacs/utility/unique_cptr.h"

#include "alex_modules.h"
#include "gentop_core.h"
#include "gmx_simple_comm.h"
#include "memory_check.h"
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
        int        numberCalculateDeviation_     = 0;
        int        numberCalcDevCalled_          = 0;
        gmx::unique_cptr<FILE, my_fclose> fplog_ = nullptr;
        std::string outputFile_;
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
        
        void saveState();

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
            for (int i = 0; i < asize(pa); i++)
            {
                pargs->push_back(pa[i]);
            }
            addOptions(pargs, etuneEEM);
            Bayes::add_pargs(pargs);
        }

        void optionsFinished(const std::string &outputFile)
        {
            MolGen::optionsFinished();
            setBoxConstraint(bConstrain());
            outputFile_ = outputFile;
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

        virtual double calcDeviation(bool verbose,
                                     bool calcAll);

        /*! \brief
         * Master sends number of calcdevs to slave
         * \param[in] bOptimize    Whether or not to optimize
         * \param[in] bSensitivity Whether or not to do sensitivity analysis
         */
        void communicateNumberOfCalculations(bool bOptimize,
                                             bool bSensitivity);
        /*! \brief
         * Do the actual optimization.
         * \param[in] fp          FILE pointer for logging
         * \param[in] oenv        Output environment for managing xvg files etc.
         * \param[in] xvgconv     Output file monitoring parameters
         * \param[in] xvgepot     Output file monitoring penalty function
         * \param[in] optimize    If true an optimization will be done
         * \param[in] sensitivity If true, a sensitivity analysis will be done
         * \return true if better parameters were found.
         */
        bool runMaster(FILE                   *fp,
                       const gmx_output_env_t *oenv,
                       const char             *xvgconv,
                       const char             *xvgepot,
                       bool                    optimize,
                       bool                    sensitivity);
        /*! \brief 
         * For the slave nodes.
         */
        void runSlave();
    
        /* \brief
         * Print the final results to the logFile.
         * \param[in] chi2_min Final chi2 in optimization
         */
        void printResults(double chi2_min);
};

void OptACM::saveState()
{
    writePoldata(outputFile_, poldata(), false);
}

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
    }
}

static void dumpQX(FILE *fp, MyMol *mol, const std::string &info)
{
    if (false && fp)
    {
        std::string label = mol->getMolname() + "-" + info;
        fprintf(fp, "%s q:", label.c_str());
        t_mdatoms *md = mol->getMdatoms();
        auto myatoms = mol->atomsConst();
        for (int i = 0; i < myatoms.nr; i++)
        {
            fprintf(fp, " %g (%g)", myatoms.atom[i].q,
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
        pr_rvecs(fp, 0, label.c_str(), mol->x().rvec_array(), myatoms.nr);
    }
}

double OptACM::calcDeviation(bool verbose,
                             bool calcAll)
{
    resetChiSquared();
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
            increaseChiSquared(ermsBOUNDS, 1, bound);
            GMX_RELEASE_ASSERT(n == param.size(), 
                               gmx::formatString("Death horror error. n=%zu param.size()=%zu", n, param.size()).c_str());
        }
    }

    if (PAR(commrec()))
    {
        if (! calcAll)
        {
            poldata()->broadcast_eemprop(commrec());
        }
    }
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupport::Local) ||
            (calcAll && (mymol.eSupp_ == eSupport::Remote)))
        {
            // Update the polarizabilities only once before the loop
            if (fit("alpha"))
            {
                mymol.UpdateIdef(poldata(), InteractionType::POLARIZATION);
            }
            // Update the electronegativity parameters
            mymol.zetaToAtoms(poldata(), mymol.atoms());
            // Run charge generation including shell minimization
            auto imm = mymol.GenerateAcmCharges(poldata(), commrec(),
                                                qcycle(), qtol());

            // Check whether we have to disable this compound
            if (immStatus::OK != imm && removeMol())
            {
                mymol.eSupp_ = eSupport::No;
                continue;
            }
            
            if (weight(ermsCHARGE))
            {
                double qtot = 0;
                int    i;
                auto myatoms = mymol.atomsConst();
                for (int j = i = 0; j < myatoms.nr; j++)
                {
                    auto atype = poldata()->findParticleType(*myatoms.atomtype[j]);
                    auto qparm = atype->parameterConst("charge");
                    double qj  = myatoms.atom[j].q;
                    qtot += qj;
                    switch (qparm.mutability())
                    {
                    case Mutability::Fixed:
                        if (qparm.value() != qj)
                        {
                            GMX_THROW(gmx::InternalError(gmx::formatString("Fixed charge for atom %s in %s was changed from %g to %g",
                                                                           *myatoms.atomname[j], mymol.getMolname().c_str(), qparm.value(), qj).c_str()));
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
                            increaseChiSquared(ermsCHARGE, 1, dq*dq);
                        }
                        break;
                    default:
                        break;
                    }
                }
                increaseChiSquared(ermsCHARGE, 1, gmx::square(qtot - mymol.totalCharge()));
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
                    mymol.QgenResp_->updateZeta(mymol.atoms(), poldata());
                }
                dumpQX(logFile(), &mymol, "ESP");
                mymol.QgenResp_->updateAtomCharges(mymol.atoms());
                mymol.QgenResp_->calcPot(poldata()->getEpsilonR());
                auto myRms =
                    convertToGromacs(mymol.QgenResp_->getRms(&rrms, &cosangle),
                                     "Hartree/e");
                auto nEsp = mymol.QgenResp_->nEsp();
                increaseChiSquared(ermsESP, nEsp, gmx::square(myRms)*nEsp);
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
                mymol.rotateDipole(mymol.muQM(qType::Calc), mymol.muQM(qType::Elec));
                if (bQM())
                {
                    rvec dmu;
                    rvec_sub(mymol.muQM(qType::Calc), mymol.muQM(qType::Elec), dmu);
                    increaseChiSquared(ermsMU, 1, iprod(dmu, dmu));
                }
                else
                {
                    increaseChiSquared(ermsMU, 1, gmx::square(mymol.dipQM(qType::Calc) - mymol.dipExper()));
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
                            increaseChiSquared(ermsQUAD, 1,
                                               gmx::square(mymol.QQM(qType::Calc)[mm][nn] - mymol.QQM(qType::Elec)[mm][nn]));
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
                increaseChiSquared(ermsPolar, 1, diff2);
            }
        }
    }
    sumChiSquared(!calcAll);
    if (verbose && logFile())
    {
        printParameters(logFile());
        printChiSquared(logFile());
    }
    numberCalcDevCalled_ += 1;
    return chiSquared(ermsTOT);
}

void OptACM::InitOpt(bool bRandom)
{
    for(auto &optIndex : optIndex_)
    {
        auto param = poldata()->findForcesConst(optIndex.iType()).findParameterTypeConst(optIndex.id(), optIndex.type());
        if (param.ntrain() >= mindata())
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

void OptACM::printResults(double chi2_min)
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
                fprintf(logFile(), "%-10s  Initial: %6.3f  Best: %6.3f  Mean: %6.3f  Sigma: %6.3f  Attempted moves: %4d  Acceptance ratio: %5.1f%%\n",
                        paramNames[k].c_str(), i_param[k], best[k], pmean[k], psigma[k], attemptedMoves[k], acceptance_ratio);
            }
        }
    }
}

void OptACM::communicateNumberOfCalculations(bool bOptimize,
                                             bool bSensitivity)
{
    if (MASTER(commrec()))
    {
        if (PAR(commrec()))
        {
            numberCalculateDeviation_ = 0;
            if (bOptimize)
            {
                numberCalculateDeviation_ += 2 + Bayes::maxIter()*Bayes::nParam();
            }
            if (bSensitivity)
            {
                numberCalculateDeviation_ += 1 + 4*Bayes::nParam();
            }
            fprintf(logFile(), "Will do %d force evaluations on %d processors\n",
                    numberCalculateDeviation_, commrec()->nnodes);
            for (int dest = 1; dest < commrec()->nnodes; dest++)
            {
                gmx_send_int(commrec(), dest, numberCalculateDeviation_);
            }
        }
    }
    else
    {
        numberCalculateDeviation_ = gmx_recv_int(commrec(), 0);
    }
}

bool OptACM::runMaster(FILE                   *fp,
                       const gmx_output_env_t *oenv,
                       const char             *xvgconv,
                       const char             *xvgepot,
                       bool                    optimize,
                       bool                    sensitivity)
{
    bool bMinimum = false;
    GMX_RELEASE_ASSERT(MASTER(commrec()), "WTF");
    
    print_memory_usage(logFile());
    std::vector<std::string> paramClass;
    for(const auto &fm : typesToFit())
    {
        paramClass.push_back(fm.first);
    }
    if (optimize)
    {
        Bayes::setOutputFiles(xvgconv, paramClass, xvgepot, oenv);
        double     chi2_min = calcDeviation(true, false);
        fprintf(logFile(), "Initial chi2 value %g\n", chi2_min);
        printChiSquared(logFile());
        
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
    if (sensitivity)
    {
        Bayes::SensitivityAnalysis(logFile());
    }
    
    if (bMinimum)
    {
        auto best = Bayes::getBestParam();
        if (best.empty())
        {
            GMX_THROW(gmx::InternalError("Minimum found but not best parameters"));
        }
        // Restore best parameter set
        Bayes::setParam(best);
        // Copy it to Poldata
        std::vector<bool> changed;
        changed.resize(best.size(), true);
        toPoldata(changed);
        // Compute the deviation once more
        double chi2 = calcDeviation(true, true);
        printResults(chi2);
        printChiSquared(fp);
        printChiSquared(logFile());
        fprintf(logFile(), "Minimum chi2 %g\n", chi2);
    }
    else
    {
        fprintf(logFile(), "Did not find a better parameter set\n");
    }
    return bMinimum;
}

void OptACM::runSlave()
{
    /* S L A V E   N O D E S */
    for (auto n = 0; n < numberCalculateDeviation_; n++)
    {
        (void) calcDeviation(false, false);
    }
    //setFinal();
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
        { efXVG, "-epot",      "param-epot",    ffWRITE }
    };

    const int                   NFILE         = asize(fnm);

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
    bool                        bSensitivity  = true;
    bool                        bForceOutput  = true;
    bool                        useOffset     = false;
    
    static const char          *select_types[]   = {nullptr, "Train", "Test", "Ignore", "Unknown", nullptr};

    t_pargs                     pa[]         = {
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
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
        { "-sensitivity",  FALSE, etBOOL, {&bSensitivity},
          "Do a sensitivity analysis." },
        { "-force_output", FALSE, etBOOL, {&bForceOutput},
          "Write output even if no new minimum is found" },
        { "-select", FALSE, etENUM, {select_types},
          "Select type for making the dataset for training or testing." }

    };

    gmx_output_env_t           *oenv;
    MolSelect                   gms;

    std::vector<t_pargs>        pargs;
    for (int i = 0; i < asize(pa); i++)
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
    
    opt.optionsFinished(opt2fn("-o", NFILE, fnm));
    
    if (MASTER(opt.commrec()))
    {
        opt.openLogFile(opt2fn("-g", NFILE, fnm));
        print_memory_usage(opt.logFile());
        print_header(opt.logFile(), pargs);
        gms.read(opt2fn_null("-sel", NFILE, fnm));
        print_memory_usage(opt.logFile());
    }

    const char *tabfn = opt2fn_null("-table", NFILE, fnm);   
    
    iMolSelect select_type;
    if (!name2molselect(select_types[0], &select_type))
    {
        gmx_fatal(FARGS, "No such selection type %s", select_types[0]);
    }
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
    if (MASTER(opt.commrec()))
    {
        if (bOptimize || bSensitivity)
        {
            opt.InitOpt(bRandom);
        }
        opt.communicateNumberOfCalculations(bOptimize, bSensitivity);
        bMinimum = opt.runMaster(MASTER(opt.commrec()) ? stderr : nullptr,
                                 oenv,
                                 opt2fn("-conv", NFILE, fnm),
                                 opt2fn("-epot", NFILE, fnm),
                                 bOptimize,
                                 bSensitivity);
    }
    else if (bOptimize || bSensitivity)
    {
        opt.communicateNumberOfCalculations(bOptimize, bSensitivity);
        opt.runSlave();
    }
   
    if (MASTER(opt.commrec()))
    {
        if (bMinimum || bForceOutput || !bOptimize)
        {
            bool bPolar = opt.poldata()->polarizable();
            auto mymols = opt.mymols();
            print_electric_props(opt.logFile(),
                                 &mymols,
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
                                 useOffset);
            print_memory_usage(opt.logFile());
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
    }
    return 0;
}
