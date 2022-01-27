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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "actpre.h"

#include "tune_ff.h"

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

#include "acm_ga.h"
#include "acthelper.h"
#include "actmiddleman.h"
#include "alex_modules.h"
#include "bayes.h"
#include "memory_check.h"
#include "mcmcmutator.h"
#include "molgen.h"
#include "molprop_util.h"
#include "mymol_low.h"
#include "npointcrossover.h"
#include "percentmutator.h"
#include "poldata.h"
#include "poldata_tables.h"
#include "poldata_xml.h"
#include "tuning_utility.h"
#include "units.h"

namespace alexandria
{

/*! \defgroup tune_ff Schematic flowchart for tune_ff
 *
 * This diagram shows how the program is parallelized.
 * Boxes in cyan run on all nodes, pink on the master node
 * and yellow on the helper nodes.
 * \dot
digraph tune_ff {
    compound = true;
    splines=true;
    node [shape=box,style=filled,color=pink] rif rm;
    rif [label="Read initial force field\nand molecules for train and test set"];
    mcmc [ label="Monte Carlo\nWrite chi-square and parameters"];

    subgraph cluster_3 {
        label = "Helpers 0 ... N-1";
        rm  [label="Calc deviation 0"];
        node [shape=box,style=filled,color=yellow];
        rh1 -> rh2 -> rh3  [style=invisible,rankdir=TB,dir=none];
        rh1  [label="Calc deviation 1"];
        rh2  [label="Calc deviation 2"];
        rh3  [label="Calc deviation N-1"];
        rm -> rh2 [style=invisible,rankdir=LR,dir=none,rank=same];
    }

    node [shape=box] [label="Start",color=cyan]; start;
    node [shape=box] [label="Parse command-line options",color=cyan]; parse;
    node [shape=diamond] [label="Master node?",color=cyan]; is_master;
    start -> parse -> is_master;
    rif -> mcmc;
    mcmc -> rm [dir=both];
    rm  -> rh1 [dir=both] [label="communication"];
    rm  -> rh2 [dir=both] [label="communication"];
    rm  -> rh3 [dir=both] [label="communication"];
    node [shape=box] [ label="Write optimized force field\nand statistics",color=pink ]; ready;
    node [shape=box][ label="Finish", color=cyan]; finished;
    is_master -> rh1 [ label="no" ];
    is_master -> rh2 [ label="no" ];
    is_master -> rh3 [ label="no" ];
    is_master -> rif [ label="yes" ];
    mcmc -> ready [ label="done" ];
    ready -> finished;
    rh1 -> finished;
    rh2 -> finished;
    rh3 -> finished;
}
 * \enddot
 */

void my_fclose(FILE *fp)
{
    const int myerrno = gmx_ffclose(fp);
    if (myerrno != 0)
    {
        fprintf(stderr, "Error %d closing file\n", myerrno);
    }
}

void OptACM::add_pargs(std::vector<t_pargs> *pargs) {
    t_pargs pa[] =
        {
            { "-fullQuadrupole", FALSE, etBOOL, {&bFullQuadrupole_},
              "Consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization"},
            { "-removemol",      FALSE, etBOOL, {&bRemoveMol_},
              "Remove a molecule from training set if shell minimization does not converge."},
            { "-v",              FALSE, etBOOL, {&verbose_},
              "Flush output immediately rather than letting the OS buffer it. Don't use for production simulations."}
        };
    for (int i = 0; i < asize(pa); i++) {
        pargs->push_back(pa[i]);
    }
    mg_.addOptions(pargs, sii_->fittingTargets(iMolSelect::Train));
    bch_.add_pargs(pargs);
    gach_.add_pargs(pargs);
}

void OptACM::check_pargs()
{
    bch_.check_pargs();
    gach_.check_pargs();
}

void OptACM::optionsFinished(const std::string &outputFile)
{
    mg_.optionsFinished();
    sii_->setOutputFile(outputFile);
    int nmiddlemen = 0;
    if (gach_.optimizer() != OptimizerAlg::MCMC)
    {
        nmiddlemen = gach_.popSize();
    }
    // Update the communication record and do necessary checks.
    commRec_.init(nmiddlemen);
}

void OptACM::openLogFile(const char *logfileName) {
    fplog_.reset(gmx_ffopen(logfileName, "w"));
}

FILE *OptACM::logFile() {
    if (fplog_) 
    {
        return fplog_.get();
    }
    else
    {
        return nullptr;
    }
}

void OptACM::initChargeGeneration(iMolSelect ims)
{
    std::string method, basis, conf, type, myref, mylot;
    splitLot(mg_.lot(), &method, &basis);
    std::vector<double> vec;
    for (MyMol &mymol : mg_.mymols())
    {
        if (mymol.datasetType() != ims)
        {
            continue;
        }
        if (mg_.fit("alpha"))
        {
            // For fitting alpha we need a reference polarizability
            double T = 0;
            auto gp = mymol.findProperty(MolPropObservable::POLARIZABILITY, iqmType::QM, T,
                                         method, basis, "");
            if (gp)
            {
                mymol.SetElectronicPolarizability(gp->getValue());
            }
            else
            {
                if (logFile())
                {
                    fprintf(logFile(), "Removing %s due to lacking reference polarizability at the %s/%s LoT.\n",
                            mymol.getMolname().c_str(),
                            method.c_str(), basis.c_str());
                }
                mymol.setSupport(eSupport::No);
            }
        }
        if (mymol.support() != eSupport::No)
        {
            mymol.setQgenAcm(new QgenAcm(sii_->poldata(), mymol.atoms(),
                                         mymol.totalCharge()));
        }
    }
}

void OptACM::initMaster(const std::string &outputFile)
{
    ga::ProbabilityComputer *probComputer = nullptr;
    // ProbabilityComputer
    switch (gach_.probabilityComputerAlg())
    {
    case ProbabilityComputerAlg::pcRANK:
        {
            probComputer = new ga::RankProbabilityComputer(gach_.popSize());
            break;
        }
    case ProbabilityComputerAlg::pcFITNESS:
        {
            probComputer = new ga::FitnessProbabilityComputer();
            break;
        }
    case ProbabilityComputerAlg::pcBOLTZMANN:
        {
            probComputer = new ga::BoltzmannProbabilityComputer(gach_.boltzTemp(), gach_.popSize());
            break;
        }
    }
    
    // Fitness computer FIXME: do we want to give the pointer to the logfile instead of nullptr?
    fitComp_ = new ACMFitnessComputer(nullptr, sii_, &mg_, 
                                      false, false, false);
    
    // Create and initialize the mutator
    ga::Mutator *mutator;
    if (gach_.optimizer() == OptimizerAlg::GA)
    {
        mutator = new alexandria::PercentMutator(sii_, gach_.percent());
    }
    else
    {
        // auto mut = new alexandria::MCMCMutator(nullptr, verbose(), &bch_, fitComp_, sii_);
        auto mut = new alexandria::MCMCMutator(logFile(), verbose(), &bch_, fitComp_, sii_);
        mut->openParamConvFiles(oenv());
        mut->openChi2ConvFile(oenv(), bch()->evaluateTestset());
        mutator = mut;
    }
    
    // Selector
    auto *selector = new ga::RouletteSelector();
    
    // Crossover
    GMX_RELEASE_ASSERT(gach_.nCrossovers() < static_cast<int>(sii_->nParam()),
                       gmx::formatString("The order of the crossover operator should be smaller than the amount of parameters. You chose -nCrossovers %i, but there are %lu parameters. Please adjust -nCrossovers.", gach_.nCrossovers(), sii_->nParam()).c_str() );
    
    auto *crossover = new alexandria::NPointCrossover(sii_->nParam(),
                                                      gach_.nCrossovers());
    
    // Terminator
    auto *terminator = new ga::GenerationTerminator(gach_.maxGenerations());
    
    if (gach_.optimizer() == OptimizerAlg::MCMC)
    {
        auto initializer = new ACMInitializer(sii_, gach_.randomInit(),
                                              outputFile, bch_.seed());
    
        ga_ = new ga::MCMC(logFile(), initializer,
                           fitComp_, probComputer,
                           selector, crossover, mutator, terminator, sii_,
                           &gach_, bch_.evaluateTestset());
    }
    else
    {
        ga_ = new ga::HybridGAMC(logFile(), nullptr,
                                 fitComp_, probComputer,
                                 selector, crossover, mutator, terminator, sii_,
                                 &gach_);
    }
}

bool OptACM::runMaster(bool        optimize,
                       bool        sensitivity)
{
    GMX_RELEASE_ASSERT(commRec_.nodeType() == NodeType::Master,
                       "I thought I was the master...");

    print_memory_usage(debug);
    bool bMinimum = false;
    if (optimize)
    {
        bMinimum = ga_->evolve(ga_->bestGenomePtr());
    }
    auto mut = static_cast<MCMCMutator *>(ga_->mutator());
    if (gach_.optimizer() != OptimizerAlg::GA && sensitivity)
    {
        // Do sensitivity analysis only on the training set
        mut->sensitivityAnalysis(ga_->bestGenomePtr(), iMolSelect::Train);
    }
    // Stop the middlemen ...
    if (commRec_.nmiddlemen() > 0)
    {
        for(auto &dest : commRec_.middlemen())
        {
            commRec_.send_done(dest);
        }
    }
    else
    {
        // ... or the helpers if there are no middlemen.
        for(auto &dest : commRec_.helpers())
        {
            commRec_.send_done(dest);
        }
    }
    if (gach_.optimizer() != OptimizerAlg::GA)
    {
        mut->printMonteCarloStatistics(logFile(), ga_->bestGenomePtr());
    }

    if (bMinimum)
    {
        auto best = ga_->bestGenome().bases();
        if (best.empty())
        {
            GMX_THROW(gmx::InternalError("Minimum found but no best parameters"));
        }
        // If MCMC was chosen as optimizer, restore best parameter set
        // TODO CHECK THIS
        if (gach_.optimizer() == OptimizerAlg::MCMC)
        {
            //bestGenome()->setParam(best);
        }
        // Copy it to Poldata
        std::vector<bool> changed;
        changed.resize(best.size(), true);
        sii_->updatePoldata(changed, ga_->bestGenomePtr());
        for (const auto &ims : iMolSelectNames())
        {
            // TODO printing
            double chi2 = fitComp_->calcDeviation(ga_->bestGenomePtr()->basesPtr(),
                                                  CalcDev::Master, ims.first);
            fprintf(logFile(), "Minimum chi2 for %s %g\n",
                    iMolSelectName(ims.first), chi2);
        }
        // Save force field of best individual
        sii_->saveState(true);
    }
    else if (optimize)
    {
        fprintf(logFile(), "Did not find a better parameter set\n");
    }
    return bMinimum;
}


int tune_ff(int argc, char *argv[])
{
    static const char          *desc[] = {
        "tune_ff read a series of molecules and corresponding experimental",
        "dipole moments from a file, and tunes parameters in an algorithm",
        "until the experimental dipole moments are reproduced by the",
        "charge generating algorithm AX as implemented in the gentop program.[PAR]",
        "Minima and maxima for the parameters can be set, these are however",
        "not strictly enforced, but rather they are penalized with a harmonic",
        "function, for which the force constant can be set explicitly.[PAR]",
        "At every reinit step parameters are changed by a random amount within",
        "the fraction set by step size, and within the boundaries given",
        "by the minima and maxima. If the [TT]-randomInit[tt] flag is",
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

    real                efield              = 10;
    bool                bcompress           = false;
    bool                bZero               = true;
    bool                bOptimize           = true;
    bool                bSensitivity        = true;
    bool                bForceOutput        = false;
    bool                bEvaluate_testset   = false;

    gmx_output_env_t           *oenv;
    MolSelect                   gms;
    TuneForceFieldPrinter       printer;

    std::vector<t_pargs>        pargs;
    {
        t_pargs                     pa[]         = {
            { "-zero", FALSE, etBOOL, {&bZero},
              "Use molecules with zero dipole in the fit as well" },
            { "-compress", FALSE, etBOOL, {&bcompress},
              "Compress output XML file" },
            { "-efield",  FALSE, etREAL, {&efield},
              "The magnitude of the external electric field to calculate polarizability tensor." },
            { "-optimize",     FALSE, etBOOL, {&bOptimize},
              "Do parameter optimization when true, or a single calculation otherwise." },
            { "-sensitivity",  FALSE, etBOOL, {&bSensitivity},
              "Do a sensitivity analysis." },
            { "-force_output", FALSE, etBOOL, {&bForceOutput},
              "Write output even if no new minimum is found" },
            { "-evaluate_testset", FALSE, etBOOL, {&bEvaluate_testset},
              "Evaluate the MCMC energy on the test set." }
        };

        for (int i = 0; i < asize(pa); i++)
        {
            pargs.push_back(pa[i]);
        }
    }
    alexandria::OptACM opt;
    opt.add_pargs(&pargs);
    printer.addOptions(&pargs);

    std::vector<t_filenm>       filenms =
    {
        { efXML, "-f",         "allmols",    ffREAD   },
        { efXML, "-d",         "gentop",     ffRDMULT },
        { efXML, "-o",         "tune_ff",    ffWRITE  },
        { efDAT, "-sel",       "molselect",  ffREAD   },
        { efLOG, "-g",         "tune_ff",    ffWRITE  },
        { efXVG, "-conv",      "param-conv", ffWRITE  },
        { efXVG, "-epot",      "param-epot", ffWRITE  }
    };

    printer.addFileOptions(&filenms);

    if (!parse_common_args(&argc,
                           argv,
                           PCA_CAN_VIEW,
                           filenms.size(),
                           filenms.data(),
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

    // Set output environment in the optimization driver
    opt.setOenv(oenv);

    // Check validity of arguments with check_pargs() in ConfigHandler(s)
    opt.check_pargs();

    // Finishing MolGen stuff and setting output file for FF in OptACM.
    // Calls optionsFinished() for MolGen instance.
    opt.optionsFinished(opt2fn("-o", filenms.size(), filenms.data()));

    // Propagate weights from training set to other sets
    // TODO: is this necessary if all processors parse the command line?
    opt.sii()->propagateWeightFittingTargets();

    if (opt.commRec()->isMaster())
    {
        opt.openLogFile(opt2fn("-g", filenms.size(), filenms.data()));
        print_memory_usage(debug);
        print_header(opt.logFile(), pargs);
        gms.read(opt2fn_null("-sel", filenms.size(), filenms.data()));
        fprintf(opt.logFile(), "Found %d Train and %d Test compounds in %s\n\n",
                gms.count(iMolSelect::Train), gms.count(iMolSelect::Test),
                opt2fn("-sel", filenms.size(), filenms.data()));
        print_memory_usage(debug);
    }

    // Figure out a logfile to pass down :)
    FILE *fp = opt.logFile() ? opt.logFile() : (debug ? debug : nullptr);

    // Read poldata in StaticIndividualInfo sii_
    {
        auto fns = opt2fns("-d", filenms.size(), filenms.data());
        GMX_RELEASE_ASSERT(fns.size() == 1 || fns.size() == opt.gach()->popSize(),
                           gmx::formatString("Please pass exactly one or %d (popSize) force field file names", opt.gach()->popSize()).c_str());
        
        int fnIndex = 0;
        if (opt.commRec()->isMiddleMan() && fns.size() > 1)
        {
            fnIndex = opt.commRec()->middleManOrdinal();
        }
        opt.sii()->fillPoldata(fp, fns[fnIndex].c_str());
        printf("On proc %d, found %d particle types\n",
               opt.sii()->commRec()->rank(),
               opt.sii()->poldata()->nParticleTypes());
    }
    // MolGen read being called here!
    if (0 == opt.mg()->Read(fp,
                            opt2fn("-f", filenms.size(), filenms.data()),
                            opt.sii()->poldata(),
                            bZero,
                            gms,
                            nullptr,
                            opt.verbose()))
    {
        if (opt.logFile())
        {
            fprintf(opt.logFile(), "Training set is empty, check your input. Rerun with -v option or -debug 1.\n");
        }
        return 0;
    }

    // StaticIndividualInfo things
    opt.sii()->generateOptimizationIndex(fp, opt.mg());
    opt.sii()->fillVectors(opt.mg()->mindata());
    opt.sii()->computeWeightedTemperature(opt.bch()->temperatureWeighting());
    // Set the output file names, has to be done before
    // creating a mutator.
    if (bOptimize)
    {
        std::vector<std::string> paramClass;
        for(const auto &fm : opt.mg()->typesToFit())
        {
            paramClass.push_back(fm.first);
        }
        opt.sii()->setOutputFiles(opt2fn("-conv", filenms.size(), filenms.data()), 
                                  paramClass,
                                  opt2fn("-epot", filenms.size(), filenms.data()));
        opt.sii()->assignParamClassIndex();  // paramClass needs to be defined when we call this!
    }

    // init charge generation for compounds in the
    // training set
    opt.initChargeGeneration(iMolSelect::Train);
    if (bEvaluate_testset)
    {
        // init charge generation for compounds in the
        // test and ignore sets
        opt.initChargeGeneration(iMolSelect::Test);
        opt.initChargeGeneration(iMolSelect::Ignore);
    }

    // Create ACMFitnessComputer and fill the DevComputers
    // This is needed on all nodes.

    if (NodeType::Master == opt.commRec()->nodeType())
    {
        opt.initMaster(opt2fn("-o", filenms.size(), filenms.data()));

        // Master only
        bool bMinimum = opt.runMaster(bOptimize, bSensitivity);
        
        if (bMinimum || bForceOutput || !bOptimize)
        {
            // if (bForceOutput)
            // {
            // FIXME: this is not true! The best parameters are fed back to params_ (and to pd_) at the end
            // of runMaster if a better parameter set was found. So no, the final step will not be the output
            fprintf(opt.logFile(), "Output based on last step of MC simulation per your specification.\nUse the -noforce_output flag to prevent this.\nThe force field output file %s is based on the last MC step as well.\n", opt2fn("-o", filenms.size(), filenms.data()));
            opt.sii()->saveState(true);
            // }
            MolGen *tmpMg = opt.mg();
            printer.print(opt.logFile(), &(tmpMg->mymols()),
                          opt.sii()->poldata(),
                          tmpMg->mdlog(), tmpMg->lot(),
                          tmpMg->qcycle(), tmpMg->qtol(),
                          oenv, opt.fullQuadrupole(),
                          opt.commRec(), efield, filenms);
            print_memory_usage(debug);
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
    }
    else if (NodeType::MiddleMan == opt.commRec()->nodeType())
    {
        // Master and Individuals (middle-men) need to initialize more,
        // so let's go.
        ACTMiddleMan middleman(opt.logFile(),
                               opt.mg(), opt.sii(), opt.gach(), opt.bch(),
                               opt2fn("-o", filenms.size(), filenms.data()),
                               opt.verbose(), opt.oenv());
        middleman.run();
    }
    else if (bOptimize || bSensitivity)
    {
        ACTHelper helper(opt.sii(), opt.mg());
        helper.run();
    }
    return 0;
}

} // namespace alexandria
