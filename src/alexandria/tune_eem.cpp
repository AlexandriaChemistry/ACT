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
        bool       bRemoveMol_                   = true;
        bool       bPointCore_                   = false;
        int        numberCalcDevCalled_          = 0;
        gmx::unique_cptr<FILE, my_fclose> fplog_ = nullptr;
        std::string outputFile_;

    public:

        OptACM() {}

        ~OptACM() {}

        bool fullTensor() const { return bFullTensor_; }

        bool sameZeta() const { return bSameZeta_; }

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
        
        /*! \brief Initialize charge generation 
         * \param[in] ims The data set to do the work for
         */
        void initChargeGeneration(iMolSelect ims);

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

        /*! \brief
         * Computes deviation from target
         * \param[in] verbose Whether or not to prin a lot
         * \param[in] calcDev The type of calculation to do
         * \param[in] ims     The data set to do computations on
         * \return the squayre deviation
         */
        virtual double calcDeviation(bool       verbose,
                                     CalcDev    calcDev,
                                     iMolSelect ims);

        /*! \brief
         * Do the actual optimization.
         * \param[in] oenv        Output environment for managing xvg files etc.
         * \param[in] xvgconv     Output file monitoring parameters
         * \param[in] xvgepot     Output file monitoring penalty function
         * \param[in] optimize    If true an optimization will be done
         * \param[in] sensitivity If true, a sensitivity analysis will be done
         * \return true if better parameters were found.
         */
        bool runMaster(const gmx_output_env_t *oenv,
                       const char             *xvgconv,
                       const char             *xvgepot,
                       bool                    optimize,
                       bool                    sensitivity,
		       bool		       bEvaluate_testset);
        /*! \brief 
         * For the slave nodes.
         */
        void runSlave();
    
};

void OptACM::saveState()
{
    writePoldata(outputFile_, poldata(), false);
}

void OptACM::initChargeGeneration(iMolSelect ims)
{
    std::string method, basis, conf, type, myref, mylot;
    splitLot(lot(), &method, &basis);
    tensor           polar      = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    rvec             vec;
    for (auto &mymol : mymols())
    {
        if (mymol.datasetType() != ims)
        {
            continue;
        }
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
            mymol.QgenAcm_ = new QgenAcm(poldata(), mymol.atoms(), 
                                         mymol.totalCharge());
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

double OptACM::calcDeviation(bool       verbose,
                             CalcDev    calcDev,
                             iMolSelect ims)
{
    auto cr = commrec();
    if (PAR(cr))
    {
        if (MASTER(cr))
        {
            for(int i = 1; i < cr->nnodes; i++)
            {
                gmx_send_int(cr, i, static_cast<int>(calcDev));
                gmx_send_int(cr, i, static_cast<int>(ims));
            }
        }
        else
        {
            calcDev = static_cast<CalcDev>(gmx_recv_int(cr, 0));
            ims     = static_cast<iMolSelect>(gmx_recv_int(cr, 0));
        }
    }
    if (calcDev == CalcDev::Final)
    {
        return -1;
    }
    resetChiSquared(ims);
    auto targets = fittingTargets(ims);
    if (MASTER(commrec()))
    {
        if ((*targets).find(eRMS::BOUNDS)->second.weight() > 0)
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
            (*targets).find(eRMS::BOUNDS)->second.increase(1, bound);
            GMX_RELEASE_ASSERT(n == param.size(), 
                               gmx::formatString("Death horror error. n=%zu param.size()=%zu", n, param.size()).c_str());
        }
    }


    if (PAR(commrec()))
    {
        if (calcDev == CalcDev::Parallel)
        {
            poldata()->broadcast_eemprop(commrec());
        }
    }
    int nmolCalculated = 0;
    for (auto &mymol : mymols())
    {
        if (ims != mymol.datasetType())
        {
            continue;
        }
        if ((mymol.eSupp_ == eSupport::Local) ||
            (calcDev == CalcDev::Master && mymol.eSupp_ == eSupport::Remote))
        {
            nmolCalculated += 1;
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
            
            if ((*targets).find(eRMS::CHARGE)->second.weight() > 0 ||
                (*targets).find(eRMS::CM5)->second.weight() > 0)
            {
                double qtot    = 0;
                int    i = 0;
                auto   myatoms = mymol.atomsConst();
                std::vector<double> qcm5;
                auto qp = mymol.qTypeProps(qType::CM5);
                if (qp)
                {
                    qcm5 = qp->charge();
                    if (debug)
                    {
                        for (int j = 0; j < myatoms.nr; j++)
                        {
                            fprintf(debug, "Charge %d. CM5 = %g ACM = %g\n", j, qcm5[j], myatoms.atom[j].q);
                        }
                    }
                }
                for (int j = 0; j < myatoms.nr; j++)
                {
                    if (myatoms.atom[j].ptype == eptShell)
                    {
                        continue;
                    }
                    auto atype = poldata()->findParticleType(*myatoms.atomtype[j]);
                    auto qparm = atype->parameterConst("charge");
                    double qj  = myatoms.atom[j].q;
                    double qjj = qj;
                    // TODO: only count in real shells
                    if (nullptr != mymol.shellfc_ && 
                        j < myatoms.nr-1 && 
                        myatoms.atom[j+1].ptype == eptShell)
                    {
                        qjj += myatoms.atom[j+1].q;
                    }
                    qtot += qjj;
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
                            if ((*targets).find(eRMS::CHARGE)->second.weight() > 0)
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
                                (*targets).find(eRMS::CHARGE)->second.increase(1, dq*dq);
                            }
                        }
                        break;
                    default:
                        break;
                    }
                    if (qp &&
                        qparm.mutability() != Mutability::Fixed &&
                        (*targets).find(eRMS::CM5)->second.weight() > 0)
                    {
                        // TODO: Add charge of shell!
                        real dq2 = gmx::square(qjj - qcm5[i]);
                        (*targets).find(eRMS::CM5)->second.increase(1, dq2);
                    }
                    i += 1;
                }
                (*targets).find(eRMS::CHARGE)->second.increase(1, gmx::square(qtot - mymol.totalCharge()));
            }
            if ((*targets).find(eRMS::ESP)->second.weight() > 0)
            {
                real rrms     = 0;
                real cosangle = 0;
                auto qgr = mymol.qTypeProps(qType::Calc)->qgenResp();
                if (nullptr != mymol.shellfc_)
                {
                    qgr->updateAtomCoords(mymol.x());
                }
                if (fit("zeta"))
                {
                    qgr->updateZeta(mymol.atoms(), poldata());
                }
                dumpQX(logFile(), &mymol, "ESP");
                qgr->updateAtomCharges(mymol.atoms());
                qgr->calcPot(poldata()->getEpsilonR());
                auto myRms =
                    convertToGromacs(qgr->getRms(&rrms, &cosangle),
                                     "Hartree/e");
                auto nEsp = qgr->nEsp();
                (*targets).find(eRMS::ESP)->second.increase(nEsp, gmx::square(myRms)*nEsp);
                if (debug)
                {
                    fprintf(debug, "%s ESPrms = %g cosangle = %g\n",
                            mymol.getMolname().c_str(),
                            myRms, cosangle);
                }
            }
            // These two things need to be present, if not the code will crash
            auto qelec = mymol.qTypeProps(qType::Elec);
            auto qcalc = mymol.qTypeProps(qType::Calc);
            if ((*targets).find(eRMS::MU)->second.weight() > 0 ||
                (*targets).find(eRMS::QUAD)->second.weight() > 0)
            {
                qcalc->setQ(mymol.atoms());
                qcalc->setX(mymol.x());
                qcalc->calcMoments();
            }
            if ((*targets).find(eRMS::MU)->second.weight() > 0)
            {
                real delta = 0;
                if (bQM())
                {
                    rvec dmu;
                    rvec_sub(qcalc->mu(), qelec->mu(), dmu);
                    delta = iprod(dmu, dmu);
                }
                else
                {
                    delta = gmx::square(qcalc->dipole() - mymol.dipExper());
                }
                (*targets).find(eRMS::MU)->second.increase(1, delta);
            }
            if ((*targets).find(eRMS::QUAD)->second.weight() > 0)
            {
                double delta    = 0; 
                for (int mm = 0; mm < DIM; mm++)
                {
                    for (int nn = 0; nn < DIM; nn++)
                    {
                        if (bFullTensor_ || mm == nn)
                        {
                            delta += gmx::square(qcalc->quad()[mm][nn] - qelec->quad()[mm][nn]);
                        }
                    }
                }
                (*targets).find(eRMS::QUAD)->second.increase(1, delta);
            }
            if ((*targets).find(eRMS::Polar)->second.weight() > 0)
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
                (*targets).find(eRMS::Polar)->second.increase(1, diff2);
            }
        }
    }
    if (debug)
    {
        printParameters(debug);
    }
    sumChiSquared(calcDev == CalcDev::Parallel, ims);
    if (verbose && logFile())
    {
        printChiSquared(logFile(), ims);
    }
    numberCalcDevCalled_ += 1;
    return (*targets).find(eRMS::TOT)->second.chiSquared();
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
                            param.ntrain(), bRandom);
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

bool OptACM::runMaster(const gmx_output_env_t *oenv,
                       const char             *xvgconv,
                       const char             *xvgepot,
                       bool                    optimize,
                       bool                    sensitivity,
                       bool 		       bEvaluate_testset)
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
        double chi2     = 0;
        if (Bayes::adaptive())
        {
            bMinimum = Bayes::Adaptive_MCMC(logFile(), &chi2);
        }
        else
        {
            bMinimum = Bayes::MCMC(logFile(), bEvaluate_testset, &chi2);
        }
    }
    if (sensitivity)
    {
        // only on the training set
        Bayes::SensitivityAnalysis(logFile(), iMolSelect::Train);
    }
    // Finalize the calculations on the slaves
    GMX_RELEASE_ASSERT(calcDeviation(false, CalcDev::Final, iMolSelect::Train) < 0, "Result for final parallel calcDeviation should be less than zero");

    printMonteCarloStatistics(logFile());
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
        for (const auto &ims : iMolSelectNames())
        {
            double chi2 = calcDeviation(true, CalcDev::Master, ims.first);
            fprintf(logFile(), "Minimum chi2 for %s %g\n",
                    iMolSelectName(ims.first), chi2);
        }
    }
    else
    {
        fprintf(logFile(), "Did not find a better parameter set\n");
    }
    return bMinimum;
}

void OptACM::runSlave()
{
    // S L A V E   N O D E S
    // The second and third variable are set by the master, but
    // we have to pass something.
    // If the result is less than zero, we are done.
    while (calcDeviation(false, CalcDev::Parallel, iMolSelect::Train) >= 0)
    {
        ;
    }

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
        "slightly, in order to speed-up local search but not global search.",
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
    bool                        bEvaluate_testset = false;    

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
        { "-evaluate_testset", FALSE, etBOOL, {&bEvaluate_testset},
          "Evaluate the MCMC energy on the test set." }

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
    
    opt.Read(opt.logFile() ? opt.logFile() : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             bZero,
             gms,
             false,
             false,
             tabfn);
             
    bool  pointCore = opt.poldata()->corePointCharge();
    
    opt.set_pointCore(pointCore);
    // init charge generation for compounds in the 
    // training set
    opt.initChargeGeneration(iMolSelect::Train);
    if (bEvaluate_testset)
    {
        // init charge generation for compounds in the
        // test set
        opt.initChargeGeneration(iMolSelect::Test);
        opt.initChargeGeneration(iMolSelect::Ignore);
    }

    bool bMinimum = false;
    if (MASTER(opt.commrec()))
    {
        if (bOptimize || bSensitivity)
        {
            opt.InitOpt(bRandom);
        }
        bMinimum = opt.runMaster(oenv,
                                 opt2fn("-conv", NFILE, fnm),
                                 opt2fn("-epot", NFILE, fnm),
                                 bOptimize,
                                 bSensitivity,
                                 bEvaluate_testset);
    }
    else if (bOptimize || bSensitivity)
    {
        opt.runSlave();
    }
   
    if (MASTER(opt.commrec()))
    {
        if (bMinimum || bForceOutput || !bOptimize)
        {
            bool bPolar = opt.poldata()->polarizable();
            auto mymols = opt.mymols();
            alexandria::print_electric_props(opt.logFile(),
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
