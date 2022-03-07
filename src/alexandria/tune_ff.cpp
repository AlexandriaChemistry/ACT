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
#include "act/utility/memory_check.h"
#include "mcmcmutator.h"
#include "molgen.h"
#include "act/molprop/molprop_util.h"
#include "mymol_low.h"
#include "npointcrossover.h"
#include "percentmutator.h"
#include "act/poldata/poldata.h"
#include "act/poldata/poldata_tables.h"
#include "act/poldata/poldata_xml.h"
#include "tuning_utility.h"
#include "act/utility/units.h"

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
    // COMMENTED TO MAKE MCMC WORK WITH ANY POPULATION SIZE
    // int nmiddlemen = 0;
    // if (gach_.optimizer() != OptimizerAlg::MCMC)
    // {
    //     nmiddlemen = gach_.popSize();
    // }
    const int nmiddlemen = gach_.popSize();  // MASTER now makes the work of a middleman too
    // Update the communication record and do necessary checks.
    commRec_.init(nmiddlemen);
    // Set prefix and id in sii_
    sii_->fillIdAndPrefix();
    // Set output file for FF parameters
    baseOutputFileName_ = outputFile;
    sii_->setOutputFile(outputFile);
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
            auto gp = reinterpret_cast<const MolecularPolarizability*>(mymol.findProperty(MolPropObservable::POLARIZABILITY, iqmType::QM,
                                                                                          T, method, basis, ""));
            if (gp)
            {
                auto qelec = mymol.qTypeProps(qType::Elec);
                qelec->setPolarizabilityTensor(gp->getTensor());
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

void OptACM::initMaster()
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

    // Fitness computer
    // FIXME: what about the flags? Here it is a bit more clear that they should be all false?
    fitComp_ = new ACMFitnessComputer(nullptr, sii_, &mg_, false, false);

    // Adjust the seed that gets passed around to components of the optimizer
    int seed = bch_.seed();
    // Create random number generator and feed it the global seed
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> dis(0); // Default constructor to cover all available (positive) range
    gen.seed(seed);
    seed = dis(gen);

    // Initializer
    auto *initializer = new ACMInitializer(sii_, gach_.randomInit(), seed);

    // Create and initialize the mutator
    ga::Mutator *mutator;
    if (gach_.optimizer() == OptimizerAlg::GA)
    {
        mutator = new alexandria::PercentMutator(sii_, seed, gach_.percent());
    }
    else
    {
        // auto mut = new alexandria::MCMCMutator(nullptr, verbose(), &bch_, fitComp_, sii_);
        auto mut = new alexandria::MCMCMutator(logFile(), verbose_, seed, &bch_, fitComp_, sii_, bch_.evaluateTestset());
        sii_->makeIndividualDir();
        mut->openParamConvFiles(oenv_);
        mut->openChi2ConvFile(oenv_);
        // if (sii_->commRec()->nmiddlemen() == 0)  // If we are running pure MCMC
        // {
        //     mut->openParamConvFiles(oenv());
        //     mut->openChi2ConvFile(oenv(), bch()->evaluateTestset());
        // }
        mutator = mut;
    }

    // Selector
    auto *selector = new ga::RouletteSelector(seed);

    // Crossover
    GMX_RELEASE_ASSERT(gach_.nCrossovers() < static_cast<int>(sii_->nParam()),
                       gmx::formatString("The order of the crossover operator should be smaller than the amount of parameters. You chose -nCrossovers %i, but there are %lu parameters. Please adjust -nCrossovers.", gach_.nCrossovers(), sii_->nParam()).c_str() );

    auto *crossover = new alexandria::NPointCrossover(sii_->nParam(),
                                                      gach_.nCrossovers(),
                                                      seed);

    // Terminator(s)
    std::vector<ga::Terminator*> terminators;
    terminators.push_back(new ga::GenerationTerminator(gach_.maxGenerations(), logFile()));  // maxGenerations will always be positive!
    if (gach_.maxTestGenerations() != -1)  // If maxTestGenerations is enabled...
    {
        fprintf(logFile(), "Appending a TestGenTerminator to the list...\n");
        terminators.push_back(new ga::TestGenTerminator(gach_.maxTestGenerations(), logFile()));
    }

    if (gach_.optimizer() == OptimizerAlg::MCMC)
    {
        // auto initializer = new ACMInitializer(sii_, gach_.randomInit(), bch_.seed());

        ga_ = new ga::MCMC(logFile(), initializer, fitComp_, mutator, sii_, &gach_);
    }
    else
    {
        // We pass the global seed to the optimizer
        ga_ = new ga::HybridGAMC(logFile(), initializer, fitComp_, probComputer, selector, crossover, mutator, terminators, sii_, &gach_,
                                 bch_.seed());
    }
}

void OptACM::printNumCalcDevEstimate()
{
    long nCalcDevTrain = 0;
    long nCalcDevTest = 0;
    long nCalcDevIgnore = 0;

    if (gach_.optimizer() == OptimizerAlg::MCMC)
    {
        nCalcDevTrain = bch_.maxIter() * sii_->nParam() + 1; // Initial one
        if (bch_.evaluateTestset())
        {
            nCalcDevTest = bch_.maxIter() + 1;  // Initial one
        }
    }
    else if (gach_.optimizer() == OptimizerAlg::GA)
    {
        nCalcDevTrain = gach_.maxGenerations() + 1;  // Extra initial generation
        nCalcDevTest  = gach_.maxGenerations() + 1;  // Extra initial generation
    }
    else if (gach_.optimizer() == OptimizerAlg::HYBRID)
    {
        nCalcDevTrain = (bch_.maxIter() * sii_->nParam() + 1) * gach_.maxGenerations() + 1;
        if (bch_.evaluateTestset())
        {
            // Extra one per generation at the end (on top of the MCMC one)
            nCalcDevTest = (bch_.maxIter() + 2) * gach_.maxGenerations();
        }
    }

    // Multiply by amount of individuals
    nCalcDevTrain *= gach_.popSize();
    nCalcDevTest *= gach_.popSize();

    nCalcDevTrain += 1;  // At the end, for the best
    nCalcDevTest += 1;  // At the end, for the best
    nCalcDevIgnore += 1;  // At the end, for the best

    fprintf(
        logFile(),
        "\nTotal number of fitness computations to be done:\n  - Train: %ld\n  - Test: %ld\n  - Ignore: %ld\n\n",
        nCalcDevTrain, nCalcDevTest, nCalcDevIgnore
    );
}

bool OptACM::runMaster(bool        optimize,
                       bool        sensitivity)
{
    GMX_RELEASE_ASSERT(commRec_.nodeType() == NodeType::Master,
                       "I thought I was the master...");

    print_memory_usage(debug);
    bool bMinimum = false;
    ga::Genome bestGenome;
    if (optimize)
    {
        // Estimate number of fitness computations per dataset
        if (logFile())
        {
            printNumCalcDevEstimate();
        }
        // Optimize!
        bMinimum = ga_->evolve(&bestGenome);
    }
    if (gach_.optimizer() != OptimizerAlg::GA && sensitivity)
    {
        // Do sensitivity analysis only on the training set
        // mut->sensitivityAnalysis(&bestGenome, iMolSelect::Train);
    }
    // Stop the middlemen ...
    if (commRec_.nmiddlemen() > 1)
    {
        for(auto &dest : commRec_.middlemen())
        {
            commRec_.send_done(dest);
        }
    }

    // else
    // {
    //     // ... or the helpers if there are no middlemen.
    //     for(auto &dest : commRec_.helpers())
    //     {
    //         commRec_.send_done(dest);
    //     }
    // }

    if (bMinimum)
    {
        if (bestGenome.bases().empty())
        {
            GMX_THROW(gmx::InternalError("Minimum found but no best parameters"));
        }
        
        for (const auto &ims : iMolSelectNames())
        {
            // TODO printing
            fitComp_->compute(&bestGenome, ims.first);
            double chi2 = bestGenome.fitness(ims.first);
            fprintf(logFile(), "Minimum chi2 for %s %g\n",
                    iMolSelectName(ims.first), chi2);
    }
        // Save force field of best individual
        sii_->setFinalOutputFile(baseOutputFileName_);
        sii_->saveState(true);
    }
    else if (optimize)
    {
        fprintf(logFile(), "Did not find a better parameter set\n");
    }

    // Stop MASTER's helpers
    std::vector<double> dummy;
    fitComp_->calcDeviation(&dummy, CalcDev::Final, iMolSelect::Train);

    return bMinimum;
}


int tune_ff(int argc, char *argv[])
{
    static const char          *desc[] = {
        "tune_ff read a series of molecules and corresponding physicochemical",
        "properties from a file. The properties can either originate from",
        "experimental data or from quantum chemistry calculations.",
        "The program then tunes empirical force field parameters using",
        "one of a series of algorithms",
        "until the experimental properties are reproduced reasonably.",
        "The properties include the electrostatic potential and multipole",
        "moments generated by molecules and the molecular energy.[PAR]"
        "Minima and maxima for some parameters can be set that are",
        "not strictly enforced, but rather they are penalized with a harmonic",
        "function. For other parameters (in particular atomic charges) strict",
        "enforcement is not possible, but unrealistic charges can be",
        "penalized with a harmonic force for  which the force constant",
        "can be set explicitly on the command line.[PAR]",
        "The algorithmes available for tuning are:[PAR]", 
        "1) a Monte Carlo algorithm",
        "that can be combined with simulated annealing.", 
        "Choose [TT]-optimizer MCMC[tt] to select this option,[PAR]", 
        "2) a Genetic Algorithm, that can be selected with [TT]-optimizer GA[tt]",
        "and[PAR]", 
        "3) a hybrid GA/MC algorithm that uses a Genetic Algorithm at the",
        "top level and uses MC as the mutator part of the algorithm.", 
        "Choose [TT]-optimizer HYBRID[tt]",
        "for this option.[PAR]", 
        "All algorithms work in parallel using MPI, under",
        "some constraints. The GA and HYBRID need an even number of processor cores.",
        "Each individual [TT]P[tt] (set with the [TT]-popSize[tt] flag) can use",
        "multiple ([TT]M[tt]) additional processor cores to compute the fitness. The total",
        "number of cores [TT]N[tt] should therefore be exactly equal to ",
        "[TT]PxM[tt].[PAR]",
        "If the [TT]-randomInit[tt] flag is",
        "given a completely random set of parameters is generated at the start",
        "of each run, within the bounds given in the input force field file.[PAR]"
        "The absolut dipole moment of a molecule remains unchanged if all the",
        "atoms swap the sign of the charge. To prevent this kind of mirror",
        "effects a penalty is added to the square deviation ",
        "if hydrogen atoms have a negative charge. Similarly a penalty is",
        "added if atoms from row VI or VII in the periodic table have a positive",
        "charge. The penalty is equal to the force constant given on the command line",
        "time the square of the charge.[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[PAR]",
        "iupac|Train[PAR]",
        "iupac|Test[PAR]",
        "iupac|Ignore[PAR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored. Selection files",
        " can be generated using the [TT]molselect[tt] script."
    };

    real                efield              = 1;
    bool                bcompress           = false;
    bool                bOptimize           = true;
    bool                bSensitivity        = true;
    bool                bForceOutput        = false;

    gmx_output_env_t           *oenv;
    MolSelect                   gms;
    TuneForceFieldPrinter       printer;

    std::vector<t_pargs>        pargs;
    {
        t_pargs                     pa[]         = {
            { "-compress", FALSE, etBOOL, {&bcompress},
              "Compress output XML file" },
            { "-efield",  FALSE, etREAL, {&efield},
              "The magnitude of the external electric field to calculate polarizability tensor." },
            { "-optimize",     FALSE, etBOOL, {&bOptimize},
              "Do parameter optimization when true, or a single calculation otherwise." },
            { "-sensitivity",  FALSE, etBOOL, {&bSensitivity},
              "Do a sensitivity analysis." }
            // { "-force_output", FALSE, etBOOL, {&bForceOutput},
            //   "Write output force field even if no new minimum is found beyond the initial set of candidate solutions." },
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
    // Initializes commRec_ in opt
    // Fills id and prefix in sii
    // Calls optionsFinished() for MolGen instance.
    opt.optionsFinished(opt2fn("-o", filenms.size(), filenms.data()));

    // Propagate weights from training set to other sets
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
                            gms,
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
    if (opt.bch()->evaluateTestset())
    {
        // init charge generation for compounds in the
        // test and ignore sets
        opt.initChargeGeneration(iMolSelect::Test);
        opt.initChargeGeneration(iMolSelect::Ignore);
    }

    // Create ACMFitnessComputer and fill the DevComputers
    // This is needed on all nodes.

    if (opt.commRec()->isMaster())
    {
        opt.initMaster();

        // Master only
        bool bMinimum = opt.runMaster(bOptimize, bSensitivity);

        if (bMinimum || bForceOutput || !bOptimize)
        {
            if (bForceOutput && !bMinimum)
            {
                fprintf(opt.logFile(), "No better minimum than the best initial candidate solution was found but -force_output was selected. This means that a global best force field file %s has been written written anyway.\n", opt2fn("-o", filenms.size(), filenms.data()));
                // fprintf(opt.logFile(), "Output based on last step of MC simulation per your specification.\nUse the -noforce_output flag to prevent this.\nThe force field output file %s is based on the last MC step as well.\n", opt2fn("-o", filenms.size(), filenms.data()));
                opt.sii()->saveState(true);
            }
            MolGen *tmpMg = opt.mg();
            printer.print(opt.logFile(), &(tmpMg->mymols()),
                          opt.sii()->poldata(),
                          tmpMg->mdlog(), tmpMg->lot(),
                          oenv, opt.commRec(), efield, filenms);
            print_memory_usage(debug);
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
    }
    else if (opt.commRec()->isMiddleMan())
    {
        // Master and Individuals (middle-men) need to initialize more,
        // so let's go.
        ACTMiddleMan middleman(opt.mg(), opt.sii(), opt.gach(), opt.bch(),
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
