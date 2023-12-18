/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "actpre.h"

#include "train_ff.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <random>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"

#include "acm_ga.h"
#include "acthelper.h"
#include "actmiddleman.h"
#include "alex_modules.h"
#include "bayes.h"
#include "act/utility/memory_check.h"
#include "mcmcmutator.h"
#include "molgen.h"
#include "act/molprop/molprop_util.h"
#include "actmol_low.h"
#include "act/ga/npointcrossover.h"
#include "act/ga/penalizer.h"
#include "percentmutator.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_tables.h"
#include "act/forcefield/forcefield_xml.h"
#include "train_utility.h"
#include "act/utility/units.h"

namespace alexandria
{

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
            { "-flush",              FALSE, etBOOL, {&flush_},
              "Flush output immediately rather than letting the OS buffer it. Don't use for production simulations."},
            { "-v",              FALSE, etBOOL, {&verbose_},
              "Print extra information to the log file during optimization. Also create convergence files for all parameters."}
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
    const int nmiddlemen = gach_.popSize();  // MASTER now makes the work of a middleman too
    if (debug)
    {
        fprintf(debug, "nmiddlemen = %d\n", nmiddlemen);
        fflush(debug);
    }
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

// TODO rename function and make a coupling between targets from the command line
// and properties to fetch from the training dataset.
void OptACM::initChargeGeneration(iMolSelect ims)
{
    std::string method, basis, conf, type, myref, mylot;
    std::vector<double> vec;
    auto mymols = mg_.actmolsPtr();
    for (auto actmol = mymols->begin(); actmol < mymols->end(); ++actmol)
    {
        if (actmol->datasetType() != ims)
        {
            continue;
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
            probComputer = new ga::FitnessProbabilityComputer(gach_.popSize());
            break;
        }
    case ProbabilityComputerAlg::pcBOLTZMANN:
        {
            probComputer = new ga::BoltzmannProbabilityComputer(
                gach_.boltzTemp(), gach_.maxGenerations(),
                gach_.boltzAnneal(), gach_.popSize()
            );
            break;
        }
    }
    // Force computer
    forceComp_ = new ForceComputer;
    // Fitness computer
    // FIXME: what about the flags? Here it is a bit more clear that they should be all false?
    fitComp_ = new ACMFitnessComputer(logFile(), verbose_, sii_, &mg_, false, forceComp_);

    // Adjust the seed that gets passed around to components of the optimizer
    int seed = bch_.seed();
    if (0 == seed)
    {
        // Create random number generator and feed it the global seed
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        seed = rd();
    }
    
    // Standard mersenne_twister_engine seeded with seed
    std::mt19937 gen(seed);
    // Default constructor to cover all available (positive) range
    std::uniform_int_distribution<int> dis(0);
    // Distribute different random number seeds to the middlemen
    for(const auto mman : sii_->commRec()->middlemen())
    {
        // Generate new seed for each of the middlemen
        // and send it over.
        sii_->commRec()->send_int(mman, dis(gen));
    }

    // Initializer
    auto *initializer = new ACMInitializer(sii_, gach_.randomInit(), dis(gen));

    // Do this only when explicitly requested.
    if (bch_.checkPoint())
    {
        sii_->makeIndividualDir();  // We need to call this before opening working files!
    }
    if (gach_.optimizer() == OptimizerAlg::GA)
    {
        mutator_ = new alexandria::PercentMutator(sii_, dis(gen), gach_.percent());
    }
    else
    {
        auto mut = new alexandria::MCMCMutator(logFile(), verbose_, flush_, dis(gen), &bch_, fitComp_, sii_, bch_.evaluateTestset());
        // TODO Only open these files when we are optimizing in verbose mode.
        if (verbose_)
        {
            printf("Will open param conv fle on the master!\n");
            mut->openParamConvFiles(oenv_);
            mut->openChi2ConvFile(oenv_);
        }
        mutator_ = mut;
    }

    // Selector
    auto *selector = new ga::RouletteSelector(dis(gen));

    // Crossover
    GMX_RELEASE_ASSERT(gach_.nCrossovers() < static_cast<int>(sii_->nParam()),
                       gmx::formatString("The order of the crossover operator should be smaller than the amount of parameters. You chose -n_crossovers %i, but there are %lu parameters. Please adjust -n_crossovers.", gach_.nCrossovers(), sii_->nParam()).c_str() );

    auto *crossover = new ga::NPointCrossover(sii_->nParam(),
                                              gach_.nCrossovers(),
                                              dis(gen));

    // Penalizer(s)
    std::vector<ga::Penalizer*> *penalizers = new std::vector<ga::Penalizer*>();
    // VolumeFractionPenalizer
    const double totalVolume = sii_->getParamSpaceVolume(gach_.logVolume());
    if (logFile())
    {
        fprintf(
            logFile(),
            "\nTotal %s(hyper)volume of the parameter space is %g.\n",
            gach_.logVolume() ? "log " : "",
            totalVolume
        );
    }
    if (gach_.vfpVolFracLimit() != -1)  // VolumeFractionPenalizer enabled
    {
        if (logFile())
        {
            fprintf(
                logFile(),
                "Appending a VolumeFractionPenalizer to the list of penalizers...\n"
            );
        }
        penalizers->push_back(
            new ga::VolumeFractionPenalizer(
                oenv_, gach_.logVolume(), logFile(), totalVolume,
                gach_.vfpVolFracLimit(), gach_.vfpPopFrac(), initializer
            )
        );
    }
    if (gach_.cpGenInterval() != -1)  // CatastrophePenalizer enabled
    {
        if (logFile())
        {
            fprintf(
                logFile(),
                "Appending a CatastrophePenalizer to the list of penalizers...\n"
            );
        }
        penalizers->push_back(
            new ga::CatastrophePenalizer(logFile(), dis(gen),
                                         gach_.cpGenInterval(),
                                         gach_.cpPopFrac(),
                                         initializer, gach_.popSize()
            )
        );
    }

    // Terminator(s)
    std::vector<ga::Terminator*> *terminators = new std::vector<ga::Terminator*>;
    if (logFile())
    {
        fprintf(
            logFile(),
            "Appending a GenerationTerminator to the list of terminators...\n"
        );
    }
    terminators->push_back(new ga::GenerationTerminator(gach_.maxGenerations(), logFile()));  // maxGenerations will always be positive!
    if (gach_.maxTestGenerations() != -1)  // If maxTestGenerations is enabled...
    {
        if (logFile())
        {
            fprintf(
                logFile(),
                "Appending a TestGenTerminator to the list of terminators...\n"
            );
        }
        terminators->push_back(new ga::TestGenTerminator(gach_.maxTestGenerations(), logFile()));
    }

    // Initialize the optimizer
    if (gach_.optimizer() == OptimizerAlg::MCMC)
    {
        ga_ = new ga::MCMC(logFile(), initializer, fitComp_, mutator_, sii_, &gach_);
    }
    else
    {
        // We pass the global seed to the optimizer
        ga_ = new ga::HybridGAMC(
            logFile(), initializer, fitComp_, probComputer, selector, crossover,
            mutator_, terminators, penalizers, sii_, &gach_, dis(gen)
        );
    }
    if (logFile())
    {
        fprintf(logFile(), "Done initializing master node\n");
        fflush(logFile());
    }
}

void OptACM::printNumCalcDevEstimate()
{
    long nCalcDevTrain = 0;
    long nCalcDevTest = 0;
    long nCalcDevIgnore = 0;

    if (gach_.optimizer() == OptimizerAlg::MCMC)
    {
        nCalcDevTrain += bch_.maxIter() * sii_->nParam() + 1 + 1; // Initial one by MCMC mutator and initial one by middlemen/master
        if (bch_.evaluateTestset())
        {
            nCalcDevTest += bch_.maxIter() + 1;  // Initial one in MCMC mutator
        }
    }
    else if (gach_.optimizer() == OptimizerAlg::GA)
    {
        nCalcDevTrain += gach_.maxGenerations() + 1;  // Extra initial generation
        if (gach_.cpGenInterval() > 0)  // Extra when catastrophe penalizer
        {
            nCalcDevTrain += gach_.maxGenerations() / gach_.cpGenInterval();
        }
        if (gach_.evaluateTestset())
        {
            nCalcDevTest += gach_.maxGenerations() + 1;  // Extra initial generation
        }
    }
    else if (gach_.optimizer() == OptimizerAlg::HYBRID)
    {
        nCalcDevTrain += (bch_.maxIter() * sii_->nParam() + 1) * gach_.maxGenerations() + 1;
        if (gach_.cpGenInterval() > 0)  // Extra when catastrophe penalizer
        {
            nCalcDevTrain += gach_.maxGenerations() / gach_.cpGenInterval();
        }
        if (bch_.evaluateTestset())
        {
            nCalcDevTest += (bch_.maxIter() + 1) * gach_.maxGenerations();
        }
        if (gach_.evaluateTestset())
        {
            nCalcDevTest += gach_.maxGenerations() + 1;  // Extra initial generation
        }
    }

    // Multiply by amount of individuals
    nCalcDevTrain  *= gach_.popSize();
    nCalcDevTest   *= gach_.popSize();
    nCalcDevIgnore *= gach_.popSize();

    // Add extra ones just in master
    nCalcDevTest += 1;  // Always done in master to print the components of the initial vector
    nCalcDevTrain += 1;  // At the end, for the best
    nCalcDevTest += 1;  // At the end, for the best
    nCalcDevIgnore += 1;  // At the end, for the best

    if (gach_.optimizer() == OptimizerAlg::MCMC)
    {
        fprintf(
            logFile(),
            "\nMaximum number of fitness computations to be done for MCMC:\n  - Train: %ld\n  - Test: %ld\n  - Ignore: %ld\n\n",
            nCalcDevTrain, nCalcDevTest, nCalcDevIgnore
        );
    }
    else  // GA/HYBRID
    {
        fprintf(
            logFile(),
            "\nMaximum number of fitness computations to be done for %d generations:\n  - Train: %ld\n  - Test: %ld\n  - Ignore: %ld\nConsider that in GA and HYBRID each penalty conveys an extra fitness computation per genome for the training set. The count above includes the Catastrophe penalizer but we cannot predict when the VolumePenalizer will kick in...\n\n",
            gach_.maxGenerations(), nCalcDevTrain, nCalcDevTest, nCalcDevIgnore
        );
    }
}

void OptACM::printGenomeTable(const std::map<iMolSelect, ga::Genome> &genome,
                              const ga::GenePool                     &pop)
{
    if (!logFile())
    {
        return;
    }
    // Get header names
    std::vector<std::string> headerNames{ "CLASS", "NAME" };
    for (const auto &pair : genome)
    {
        headerNames.push_back(gmx::formatString("BEST (%s)", iMolSelectName(pair.first)));
    }
    const std::vector<std::string> tmpHeaderNames{"MIN", "MAX", "MEAN", "STDEV", "MEDIAN"};
    headerNames.insert(headerNames.end(), tmpHeaderNames.begin(), tmpHeaderNames.end());
    // Get header sizes
    const int FLOAT_SIZE = 14;  // Adjusted for the %g formatting plus negative numbers
    std::vector<int> sizes{5, 4};
    for (size_t i = 0; i < genome.size(); i++)
    {
        sizes.push_back(FLOAT_SIZE);
    }
    std::vector<int> tmpSizes{FLOAT_SIZE, FLOAT_SIZE, FLOAT_SIZE, FLOAT_SIZE, FLOAT_SIZE};
    sizes.insert(sizes.end(), tmpSizes.begin(), tmpSizes.end());
    // Adjust size for class field
    const auto paramClass = sii_->paramClass();
    for (auto pClass : paramClass)
    {
        if (static_cast<int>(pClass.size()) > sizes[0])
        {
            sizes[0] = pClass.size();
        }
    }
    // Adjust size for name field
    const auto paramClassIndex = sii_->paramClassIndex();
    const auto paramNames = sii_->paramNamesWOClass();
    for (auto pName : paramNames)
    {
        if (static_cast<int>(pName.size()) > sizes[1])
        {
            sizes[1] = pName.size();
        }
    }
    const size_t TOTAL_WIDTH = static_cast<size_t>(std::accumulate(sizes.begin(), sizes.end(), 0)) + 3*(sizes.size()-1) + 4;
    const std::string HLINE(TOTAL_WIDTH, '-');
    // Print header
    fprintf(logFile(), "%s\n|", HLINE.c_str());
    for (size_t i = 0; i < headerNames.size(); i++)
    {
        fprintf(logFile(), " %-*s |", sizes[i], headerNames[i].c_str());
    }
    fprintf(logFile(), "\n%s\n", HLINE.c_str());
    // Gather statistics from the population
    const std::vector<double> min    = pop.min();
    const std::vector<double> max    = pop.max();
    const std::vector<double> mean   = pop.mean();
    const std::vector<double> stdev  = pop.stdev();
    const std::vector<double> median = pop.median();
    // Print the parameter information
    for (size_t i = 0; i < paramClass.size(); i++)
    {
        for (size_t j = 0; j < paramNames.size(); j++)
        {
            if (paramClassIndex[j] != i)
            {
                continue;
            }
            fprintf(
                logFile(),
                "| %-*s | %-*s |",
                sizes[0], paramClass[i].c_str(),
                sizes[1], paramNames[j].c_str()
            );
            size_t k = 2;
            for (const auto &pair : genome)
            {
                fprintf(
                    logFile(),
                    " %-*g |",
                    sizes[k], pair.second.base(j)
                );
                k++;
            }
            fprintf(
                logFile(),
                " %-*g | %-*g | %-*g | %-*g | %-*g |\n%s\n",
                sizes[k], min[j],
                sizes[k+1], max[j],
                sizes[k+2], mean[j],
                sizes[k+3], stdev[j],
                sizes[k+4], median[j],
                HLINE.c_str()
            );
        }
    }
}

bool OptACM::runMaster(bool optimize,
                       bool sensitivity)
{
    GMX_RELEASE_ASSERT(commRec_.nodeType() == NodeType::Master,
                       "I thought I was the master...");

    print_memory_usage(debug);
    bool bMinimum = false;
    std::map<iMolSelect, ga::Genome> bestGenome;
    if (optimize)
    {
        // Estimate number of fitness computations per dataset
        if (logFile())
        {
            printNumCalcDevEstimate();
            fflush(logFile());
        }
        // Train the force field!
        bMinimum = ga_->evolve(&bestGenome);
        if (logFile())
        {
            fprintf(logFile(), "\nHere are the best parameters I found, together with some summary statistics of the last population:\n");
        }
        printGenomeTable(bestGenome, ga_->getLastPop());
    }
    if (gach_.optimizer() != OptimizerAlg::GA && sensitivity)
    {
        // Do sensitivity analysis only on the training set
        mutator_->sensitivityAnalysis(&bestGenome[iMolSelect::Train], iMolSelect::Train);
    }
    // Stop the middlemen ...
    if (commRec_.nmiddlemen() > 1)
    {
        for(auto &dest : commRec_.middlemen())
        {
            commRec_.send_done(dest);
        }
    }

    if (bMinimum)
    {
        if (bestGenome.empty())
        {
            GMX_THROW(gmx::InternalError("Minimum found but the <Dataset, Genome> map is empty"));
        }
        if (logFile())
        {   
            for (auto it = bestGenome.begin(); it != bestGenome.end(); it++)
            {
                it->second.print("Final best genome", logFile());
                fprintf(logFile(), "\nChi2 components of the best parameter vector found (for %s):\n", iMolSelectName(it->first));
                for (const auto &ims : iMolSelectNames())
                {
                    fitComp_->compute(&(it->second), ims.first, true);
                    double chi2 = it->second.fitness(ims.first);
                    fprintf(logFile(), "Minimum chi2 for %s %g\n", ims.second, chi2);
                }
            }
        }
        // Save force field of best individual(s)
        for (const auto &pair : bestGenome)
        {
            std::set<int> changed;
            sii_->updateForceField(changed, pair.second.bases());
            auto myfilenm = gmx::formatString("%s-", iMolSelectName(pair.first)) + baseOutputFileName_;
            fprintf(logFile(), "Will save best force field to %s\n", myfilenm.c_str());
            sii_->saveState(true, myfilenm);
        }
        // FIXME: resetting the train parameters for the TrainFFPrinter. We may have to work on that if we want to show the best test parameters too
        std::set<int> changed;
        sii_->updateForceField(changed, bestGenome.find(iMolSelect::Train)->second.bases());
    }
    else if (optimize)
    {
        if (logFile())
        {
            fprintf(logFile(), "Did not find a better parameter set\n");
        }
    }

    // Stop MASTER's helpers
    fitComp_->distributeTasks(CalcDev::Stop);
    // Turn off further communication
    sii_->commRecPtr()->done();

    // Final energy calculation for all molecules
    // TODO: parallellize this. FIXME: there is no need to do this I believe, it's done above, and parallel!
    if (!bestGenome.empty())
    {
        std::set<int> changed;
        auto ims = iMolSelect::Train;
        auto bbb = bestGenome.find(ims);
        if (bestGenome.end() != bbb)
        {
            sii_->updateForceField(changed, bbb->second.bases());
            fitComp_->calcDeviation(CalcDev::ComputeAll, ims);
        }
        // Now compute the test compounds, with the best Train parameters.
        ims = iMolSelect::Test;
        fitComp_->calcDeviation(CalcDev::ComputeAll, ims);
    }
    // Delete the penalizers
    if (nullptr != ga_->penalizers())
    {
        for (auto pen : *ga_->penalizers())
        {
            delete pen;
        }
    }

    return bMinimum;
}


int train_ff(int argc, char *argv[])
{
    static const char          *desc[] = {
        "train_ff reads a series of molecules and corresponding physico-chemical",
        "properties from a file. The properties can either originate from",
        "experimental data or from quantum chemistry calculations.",
        "The program then trains empirical force field parameters using",
        "one of a series of algorithms",
        "until the experimental properties are reproduced reasonably.",
        "The properties include the electrostatic potential and multipole",
        "moments generated by molecules, the molecular energy, atomic"
        "forces or interaction energy between compounds.[PAR]"
        "Minima and maxima for some parameters can be set that are",
        "strictly enforced. For other parameters (in particular atomic charges) strict",
        "enforcement is not possible, but unrealistic charges can be",
        "penalized with a harmonic force for  which the force constant",
        "can be set explicitly on the command line.[PAR]",
        "The algorithmes available for training are:[PAR]", 
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
        "If the [TT]-random_init[tt] flag is",
        "given a completely random set of parameters is generated at the start",
        "of each run, within the bounds given in the input force field file.[PAR]",
        "The absolute dipole moment of a molecule remains unchanged if all the",
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
        "([TT]-mp[tt] option). Missing molecules will be ignored. Selection files",
        "can be generated using the [TT]molselect[tt] script."
    };

    bool                bcompress           = false;
    bool                bOptimize           = true;
    bool                bSensitivity        = true;
    bool                bForceOutput        = false;

    gmx_output_env_t           *oenv;
    MolSelect                   gms;
    TrainForceFieldPrinter      printer;  // TODO: pargs is a ConfigHandler, maybe we could inherit the superclass?

    std::vector<t_pargs>        pargs;
    {
        t_pargs                     pa[]         = {
            { "-compress", FALSE, etBOOL, {&bcompress},
              "Compress output XML file" },
            { "-optimize",     FALSE, etBOOL, {&bOptimize},
              "Do parameter optimization when true, or a single calculation otherwise." },
            { "-sensitivity",  FALSE, etBOOL, {&bSensitivity},
              "Do a sensitivity analysis." }
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
        { efXML, "-mp",   "allmols",     ffREAD   },
        { efXML, "-ff",   "aff",         ffRDMULT },
        { efXML, "-o",    "train_ff",    ffWRITE  },
        { efDAT, "-sel",  "molselect",   ffREAD   },
        { efLOG, "-g",    "train_ff",    ffWRITE  },
        { efXVG, "-conv", "param_conv" , ffWRITE  },
        { efXVG, "-chi2", "chi_squared", ffWRITE  }
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
        print_header(opt.logFile(), pargs, filenms);

        gms.read(opt2fn_null("-sel", filenms.size(), filenms.data()));
        fprintf(opt.logFile(), "Found %d Train and %d Test compounds in %s\n\n",
                gms.count(iMolSelect::Train), gms.count(iMolSelect::Test),
                opt2fn("-sel", filenms.size(), filenms.data()));
        print_memory_usage(debug);
    }

    // Figure out a logfile to pass down :)
    FILE *fp = opt.logFile() ? opt.logFile() : (debug ? debug : nullptr);

    // Read forcefield in StaticIndividualInfo sii_
    {
        auto fns = opt2fns("-ff", filenms.size(), filenms.data());
        GMX_RELEASE_ASSERT(fns.size() == 1 || fns.size() == opt.gach()->popSize(),
                           gmx::formatString("Please pass exactly one or %d (popSize) force field file names", opt.gach()->popSize()).c_str());

        int fnIndex = 0;
        if (opt.commRec()->isMiddleMan() && fns.size() > 1)
        {
            fnIndex = opt.commRec()->middleManOrdinal();
        }
        opt.sii()->fillForceField(fp, fns[fnIndex].c_str());
        if (fp)
        {
            fprintf(fp, "On proc %d, found %d particle types\n",
                    opt.sii()->commRec()->rank(),
                    opt.sii()->forcefield()->nParticleTypes());
        }
    }

    // MolGen read being called here!
    if (0 == opt.mg()->Read(fp,
                            opt2fn("-mp", filenms.size(), filenms.data()),
                            opt.sii()->forcefield(),
                            gms,
                            opt.sii()->fittingTargetsConst(iMolSelect::Train),
                            opt.verbose()))
    {
        GMX_THROW(gmx::InvalidInputError("Training set is empty, check your input. Rerun with -v option or -debug 1"));
    }


    // StaticIndividualInfo things
    opt.sii()->generateOptimizationIndex(fp, opt.mg(), opt.commRec());
    opt.sii()->fillVectors(opt.mg()->mindata());
    opt.sii()->computeWeightedTemperature(opt.bch()->temperatureWeighting());
    // Set the output file names, has to be done before
    // creating a mutator.
    if (bOptimize)
    {
        if (opt.sii()->nParam() == 0)
        {
            if (opt.sii()->commRec()->isMaster())
            {
                fprintf(stderr, "Nothing to optimize. Check your input.\n");
            }
            return 0;
        }
        std::vector<std::string> paramClass;
        for(const auto &fm : opt.mg()->typesToFit())
        {
            paramClass.push_back(fm.first);
        }
        opt.sii()->setOutputFiles(opt2fn("-conv", filenms.size(), filenms.data()),
                                  paramClass,
                                  opt2fn("-chi2", filenms.size(), filenms.data()));
        opt.sii()->assignParamClassIndex();  // paramClass needs to be defined when we call this!
    }

    // init charge generation for compounds in the
    // training set
    opt.initChargeGeneration(iMolSelect::Train);
    if (opt.bch()->evaluateTestset() || opt.gach()->evaluateTestset())
    {
        // init charge generation for compounds in the
        // test and ignore sets
        opt.initChargeGeneration(iMolSelect::Test);
        opt.initChargeGeneration(iMolSelect::Ignore);
    }

    if (opt.commRec()->isMaster())
    {
        bool bMinimum = false;
        if (opt.sii()->nParam() > 0)
        {
            opt.initMaster();

            // Master only
            bMinimum = opt.runMaster(bOptimize, bSensitivity);
            if (bOptimize)
            {
                printf("DONE WITH OPTIMIZATION\n");
            }
        }
        if (bMinimum || bForceOutput || !bOptimize)
        {
            if (bForceOutput && !bMinimum)
            {
                fprintf(opt.logFile(), "No better minimum than the best initial candidate solution was found but -force_output was selected. This means that a global best force field file %s has been written written anyway.\n", opt2fn("-o", filenms.size(), filenms.data()));
                opt.sii()->saveState(true);
            }
            MolGen *tmpMg = opt.mg();
            printer.print(opt.logFile(), tmpMg->actmolsPtr(),
                          opt.sii()->forcefield(),
                          oenv, filenms, tmpMg->chargeMethod());
            print_memory_usage(debug);
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
    }
    else if (opt.commRec()->isMiddleMan())
    {
        if (opt.sii()->nParam() > 0)
        {
            // Master and Individuals (middle-men) need to initialize more,
            // so let's go.
            ACTMiddleMan middleman(opt.mg(), opt.sii(), opt.gach(), opt.bch(),
                                   opt.verbose(), opt.oenv(), opt.verbose());
            middleman.run();
        }
    }
    else if (bOptimize || bSensitivity)
    {
        if (opt.sii()->nParam() > 0)
        {
            ACTHelper helper(opt.sii(), opt.mg());
            helper.run();
        }
    }
    return 0;
}

} // namespace alexandria
