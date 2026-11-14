/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */

#include "actpre.h"

#include "train_ff.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <random>

#include "act/alexandria/acthelper.h"
#include "act/alexandria/actmiddleman.h"
#include "act/alexandria/alex_modules.h"
#include "act/alexandria/hybrid.h"
#include "act/alexandria/mcmc.h"
#include "act/alexandria/mcmcmutator.h"
#include "act/alexandria/molgen.h"
#include "act/alexandria/percentmutator.h"
#include "act/alexandria/queuebalance.h"
#include "act/alexandria/sensitivity.h"
#include "act/alexandria/sensitivity_analysis.h"
#include "act/alexandria/train_utility.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_tables.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/gpu/gpu.h"
#include "act/ga/penalizer.h"
#include "act/molprop/molprop_util.h"
#include "act/topology/actmol_low.h"
#include "act/utility/jsontree.h"
#include "act/utility/memory_check.h"
#include "act/utility/portable_random.h"
#include "act/utility/units.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/gmxfio.h"
#include "act/forces/simd_information.h"
#include "act/utility/timings.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"

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

OptACM::~OptACM()
{
    delete ga_;
    delete forceComp_;
    delete fitComp_;
    delete mutator_;
    delete probComputer_;
    delete initializer_;
    delete crossover_;
    delete sii_;
    delete selector_;
    // Delete the penalizers
    for (auto x : penalizers_)
    {
        delete x;
    }
    // Delete the terminators
    for (auto x : terminators_)
    {
        delete x;
    }
    if (oenv_)
    {
        output_env_done(oenv_);
        oenv_ = nullptr;
    }
}

void OptACM::add_options(std::vector<t_pargs>  *pargs,
                         std::vector<t_filenm> *fnms)
{
    std::vector<t_pargs> pa =
    {
        { "-removemol",      FALSE, etBOOL, {&bRemoveMol_},
          "Remove a molecule from training set if shell minimization does not converge."}
    };
    for (size_t i = 0; i < pa.size(); i++) {
        pargs->push_back(pa[i]);
    }
    std::vector<t_filenm> filenames = {
        { efXML, "-o",     "ff_train", ffWRITE },
        { efXML, "-otest", "ff_test",  ffWRITE },
        { efCSV, "-gpin",  "genepool", ffOPTRD },
        { efCSV, "-gpout", "genepool", ffWRITE }
    };
    for (size_t i = 0; i < filenames.size(); ++i)
    {
        fnms->push_back(filenames[i]);
    }
    mg_.addOptions(pargs, sii_->fittingTargets(iMolSelect::Train));
    mg_.addFilenames(fnms);

    bch_.add_options(pargs, fnms);
    gach_.add_options(pargs, fnms);
    static const std::string defname("train_ff");
    msghandler_.addOptions(pargs, fnms, defname);
}

void OptACM::check_pargs(MsgHandler *)
{
    bch_.check_pargs(&msghandler_);
    gach_.check_pargs(&msghandler_);
}

void OptACM::optionsFinished(const std::vector<t_filenm> &filenames)
{
    // With -nogpu_master the master does not start the GPU, before the log
    // file reports it, leaving the GPU to the other processes
    if (commRec_.isMaster() && commRec_.isParallel() && !bch_.gpuMaster())
    {
        gpu::turnOffGpu("in the master process with train_ff -nogpu_master");
    }
    msghandler_.optionsFinished(filenames, &commRec_);
    mg_.optionsFinished(&msghandler_);
    const int nmiddlemen = gach_.popSize();  // MASTER now makes the work of a middleman too
    msghandler_.writeDebug(gmx::formatString("nmiddlemen = %d", nmiddlemen));
    // Update the communication record and do necessary checks.
    commRec_.init(nmiddlemen);
    if (commRec_.queue() && gach_.optimizer() == OptimizerAlg::MCMC)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The population (%d) is larger than the number of cores/threads (%d), which only the GA and HYBRID optimizers support",
                                                           nmiddlemen, commRec_.size()).c_str()));
    }
    // Set prefix and id in sii_
    sii_->fillIdAndPrefix();
    // Set output file for FF parameters, the middle men need it for
    // their checkpoints
    sii_->setOutputFile(opt2fn("-o", filenames.size(), filenames.data()));
    if (commRec_.isMaster())
    {
        outputFileName_.insert({iMolSelect::Train, opt2fn("-o", filenames.size(), filenames.data())});
        outputFileName_.insert({iMolSelect::Test, opt2fn("-otest", filenames.size(), filenames.data())});
    }
}

void OptACM::balanceQueue(gmx::TextWriter *tw)
{
    const auto note = [tw](const std::string &text)
    {
        if (tw)
        {
            tw->writeLine(text);
        }
        fprintf(stderr, "\n%s\n", text.c_str());
    };
    // Several individuals at once only help if the GPU computes them
    // together, see ACMFitnessComputer::calcDeviations. Otherwise a middleman
    // computes them one after another, and individuals handed out one at a
    // time are shared more evenly over the middlemen.
    if (gach_.lockstep() > 1 && !forceComp_->gpuCharges(sii_->forcefield()))
    {
        note(gmx::formatString("Note: -lockstep %d has no effect here, since only the GPU computes several individuals\n"
                               "at once, for polarizable force fields with SQE charges in single precision builds.\n"
                               "Each middleman gets one individual at a time instead.", gach_.lockstep()));
        gach_.setLockstep(1);
    }
    size_t nmiddlemen = sii_->commRec()->middlemen().size();
    const auto qb         = queueBalance(gach_.popSize(), gach_.lockstep(), nmiddlemen);
    if (qb.idle >= 0.05)
    {
        const std::string group = gach_.lockstep() > 1 ?
            gmx::formatString("%zu groups of %d individuals", qb.groups, gach_.lockstep()) :
            gmx::formatString("%zu individuals", qb.groups);
        std::string text = gmx::formatString("Note: the %zu middlemen share the %s of each generation unevenly: some\n"
                                             "compute %zu, while the others wait on average %.0f%% of the time. With %zu or %zu\n"
                                             "processes (one more than the middlemen) nobody waits.",
                                             nmiddlemen, group.c_str(), qb.rounds, 100*qb.idle,
                                             qb.more + 1, qb.fewer + 1);
        if (qb.fewer == 0)
        {
            text = gmx::formatString("Note: the %zu middlemen share the %s of each generation unevenly: some\n"
                                     "compute %zu, while the others wait on average %.0f%% of the time. With %zu\n"
                                     "processes (one more than the middlemen) nobody waits.",
                                     nmiddlemen, group.c_str(), qb.rounds, 100*qb.idle, qb.more + 1);
        }
        note(text);
    }
}

void OptACM::initForceComputer()
{
    // Force computer may be needed even with nooptimize
    forceComp_ = new ForceComputer();
    double epsilonr;
    epsilonr = sii_->forcefield()->epsilonr();
    forceComp_->init(bch_.shellToler(), bch_.shellMaxIter(), epsilonr);
    forceComp_->setKeepShells(bch_.shellKeep());
    forceComp_->setGpuResident(bch_.gpuResident());
    if (bch_.gpuResident())
    {
        // Start the GPU now in all processes at once, not at the first
        // evaluation. With many processes on one GPU this may take seconds,
        // which is then not counted as part of an evaluation.
        static const size_t timer = Timings::instance().index("GPU initialization");
        const ScopedTimer   scoped(timer);
        (void) gpu::useGpu();
    }
}

//! \todo rename function and make a coupling between targets from the command line
// and properties to fetch from the training dataset.
void OptACM::initChargeGeneration(iMolSelect ims)
{
    std::string method, basis, conf, type, myref, mylot;
    const std::vector<double> vec;
    auto mymols = mg_.actmolsPtr();
    for (auto actmol = mymols->begin(); actmol < mymols->end(); ++actmol)
    {
        if (actmol->datasetType() != ims)
        {
            continue;
        }
    }
}

int OptACM::initMaster(const std::vector<t_filenm> &fnm,
                       ChargeGenerationAlgorithm    algorithm)
{
    // ProbabilityComputer
    switch (gach_.probabilityComputerAlg())
    {
    case ProbabilityComputerAlg::pcRANK:
        {
            probComputer_ = new ga::RankProbabilityComputer(gach_.popSize());
            break;
        }
    case ProbabilityComputerAlg::pcFITNESS:
        {
            probComputer_ = new ga::FitnessProbabilityComputer(gach_.popSize());
            break;
        }
    case ProbabilityComputerAlg::pcBOLTZMANN:
        {
            probComputer_ = new ga::BoltzmannProbabilityComputer(
                gach_.boltzTemp(), gach_.maxGenerations(),
                gach_.boltzAnneal(), gach_.popSize()
            );
            break;
        }
    }
    // Do this only when explicitly requested.
    if (bch_.checkPoint() || msghandler_.verbose())
    {
        sii_->makeIndividualDir();  // We need to call this before opening working files!
    }
    // Fitness computer
    fitComp_ = new ACMFitnessComputer();
    fitComp_->init(&msghandler_, sii_, &mg_, bRemoveMol_, forceComp_, algorithm, gach_.lossFunction());
    // Check whether we have to do anything
    if (fitComp_->numDevComputers() == 0)
    {
        msghandler_.msg(ACTStatus::Error,
                        "Nothing to train! Check your input and your code.\n");
        return 0;
    }
    // Adjust the seed that gets passed around to components of the optimizer
    int seed = bch_.seed();
    if (0 == seed)
    {
        // Create random number generator and feed it the global seed
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        seed = rd();
    }
    
    // Standard mersenne_twister_engine seeded with seed, from which the
    // seeds of the components are drawn in the same way on all platforms
    std::mt19937 gen(seed);
    // Distribute different random number seeds to the middlemen. With a
    // queue they do not affect the results, see CommunicationRecord::queue,
    // and they come from a generator of their own, such that the seeds of
    // the master do not depend on the number of middlemen.
    std::mt19937 mmgen(seed + 1);
    for(const auto mman : sii_->commRec()->middlemen())
    {
        // Generate new seed for each of the middlemen
        // and send it over.
        sii_->commRec()->send(mman, portableSeed(sii_->commRec()->queue() ? &mmgen : &gen));
    }

    // Initializer
    initializer_ = new ACMInitializer(sii_, gach_.randomInit(), portableSeed(&gen));

    if (gach_.optimizer() == OptimizerAlg::GA)
    {
        mutator_ = new alexandria::PercentMutator(sii_, portableSeed(&gen), algorithm, gach_.percent());
    }
    else
    {
        auto mut = new alexandria::MCMCMutator(portableSeed(&gen), &bch_, fitComp_, sii_, bch_.evaluateTestset(),
                                               algorithm, gach_.maxGenerations());
        mut->setNstFitness(gach_.nstFitness());
        //! \todo Only open these files when we are optimizing in verbose mode.
        if (msghandler_.verbose())
        {
            msghandler_.msg(ACTStatus::Verbose,
                            "Will open parameter convergence files on the master!\n");
            mut->openParamConvFiles(oenv_);
            mut->openChi2ConvFile(oenv_);
        }
        mutator_ = std::move(mut);
    }

    // Selector
    selector_ = new ga::RouletteSelector(portableSeed(&gen));

    auto tw = msghandler_.tw();
    // VolumeFractionPenalizer
    const double totalVolume = sii_->getParamSpaceVolume(gach_.logVolume());
    if (tw)
    {
        tw->writeStringFormatted("\nTotal %s(hyper)volume of the parameter space is %g.\n",
                                 gach_.logVolume() ? "log " : "",
                                 totalVolume);
    }
    if (gach_.vfpVolFracLimit() != -1)  // VolumeFractionPenalizer enabled
    {
        if (tw)
        {
            tw->writeStringFormatted("Appending a VolumeFractionPenalizer to the list of penalizers...\n");
        }
        penalizers_.push_back(
            new ga::VolumeFractionPenalizer(
                oenv_, gach_.logVolume(), totalVolume,
                gach_.vfpVolFracLimit(), gach_.vfpPopFrac(), initializer_
            )
        );
    }
    if (gach_.cpGenInterval() != -1)  // CatastrophePenalizer enabled
    {
        if (tw)
        {
            tw->writeStringFormatted("Appending a CatastrophePenalizer to the list of penalizers...\n");
        }
        penalizers_.push_back(
            new ga::CatastrophePenalizer(portableSeed(&gen),
                                         gach_.cpGenInterval(),
                                         gach_.cpPopFrac(),
                                         initializer_, gach_.popSize()
            )
        );
    }

    // Terminator(s)
    if (tw)
    {
        tw->writeStringFormatted("Appending a GenerationTerminator to the list of terminators...\n");
    }
    terminators_.push_back(new ga::GenerationTerminator(gach_.maxGenerations()));
    // maxGenerations will always be positive!
    // If maxTestGenerations is enabled and there is something to test as well
    if (gach_.maxTestGenerations() != -1 && mg_.iMolSelectSize(iMolSelect::Test) > 0)
    {
        if (tw)
        {
            tw->writeString("Appending a TestGenTerminator to the list of terminators...\n");
        }
        terminators_.push_back(new ga::TestGenTerminator(gach_.maxTestGenerations()));
    }

    // If maxLocalMinimumGenerations is enabled add a terminator
    if (gach_.maxLocalMinimumGenerations() > 0)
    {
        if (tw)
        {
            tw->writeString("Appending a CurvatureTerminator to the list of terminators...\n");
        }
        terminators_.push_back(new ga::CurvatureTerminator(gach_.maxLocalMinimumGenerations()));
    }

    // Initialize the optimizer
    if (gach_.optimizer() == OptimizerAlg::MCMC)
    {
        ga_ = new ga::MCMC(initializer_, fitComp_, mutator_, sii_, &gach_);
    }
    else
    {
        // Crossover
        if (gach_.nCrossovers() >= static_cast<int>(sii_->nParam()))
        {
            tw->writeStringFormatted("The order of the crossover operator should be smaller than the amount of parameters. You chose -n_crossovers %i, but there are %lu parameters. Changing n_crossovers to 1.",
                                     gach_.nCrossovers(), sii_->nParam());
            gach_.setCrossovers(1);
        }

        if (sii_->commRec()->queue())
        {
            balanceQueue(tw);
        }
        crossover_ = new ga::NPointCrossover(sii_->nParam(),
                                             gach_.nCrossovers(),
                                             portableSeed(&gen));

        // We pass the global seed to the optimizer
        ga_ = new ga::HybridGAMC(initializer_, fitComp_, probComputer_,
                                 selector_, crossover_,
                                 mutator_, &terminators_, &penalizers_,
                                 sii_, &gach_,
                                 opt2fn("-fitness", fnm.size(), fnm.data()),
                                 opt2fn_null("-gpin", fnm.size(), fnm.data()),
                                 opt2fn("-gpout", fnm.size(), fnm.data()),
                                 portableSeed(&gen), bch_.adaptiveStep());
    }
    if (tw)
    {
        tw->writeStringFormatted("Done initializing master node\n");
    }
    return 1;
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
    if (gach_.optimizer() != OptimizerAlg::GA)
    {
        // The sensitivity analysis after training, two per parameter and
        // in 2D four per pair of parameters
        long np = sii_->nParam();
        switch (gach_.sensitivity())
        {
        case SensitivityMode::OneD:
            nCalcDevTrain += 2*np;
            break;
        case SensitivityMode::TwoD:
            nCalcDevTrain += 2*np*np;
            break;
        case SensitivityMode::No:
            break;
        }
    }
    if (gach_.optimizer() == OptimizerAlg::HYBRID &&
        (gach_.maxLocalMinimumGenerations() > 0 || bch_.adaptiveStep()))
    {
        // The curvature analysis of the best genome in each generation,
        // two-dimensional only if the analysis after training is
        long np = sii_->nParam();
        long perGeneration = (gach_.sensitivity() == SensitivityMode::TwoD) ? 2*np*np : 2*np;
        nCalcDevTrain += perGeneration*gach_.maxGenerations();
    }
    nCalcDevTest += 1;  // Always done in master to print the components of the initial vector
    nCalcDevTrain += 1;  // At the end, for the best
    nCalcDevTest += 1;  // At the end, for the best
    nCalcDevIgnore += 1;  // At the end, for the best

    auto tw = msghandler_.tw();
    if (gach_.optimizer() == OptimizerAlg::MCMC)
    {
        tw->writeStringFormatted("\nMaximum number of fitness computations to be done for MCMC:\n  - Train: %ld\n  - Test: %ld\n  - Ignore: %ld\n\n",
                                 nCalcDevTrain, nCalcDevTest, nCalcDevIgnore);
    }
    else  // GA/HYBRID
    {
        tw->writeStringFormatted("\nMaximum number of fitness computations to be done for %d generations:\n  - Train: %ld\n  - Test: %ld\n  - Ignore: %ld\nConsider that in GA and HYBRID each penalty conveys an extra fitness computation per genome for the training set. The count above includes the Catastrophe penalizer but we cannot predict when the VolumePenalizer will kick in...\n\n",
                                 gach_.maxGenerations(), nCalcDevTrain, nCalcDevTest, nCalcDevIgnore);
    }
}

void OptACM::printGenomeTable(const std::map<iMolSelect, ga::Genome> &genome,
                              const ga::GenePool                     &pop,
                              JsonTree                               *jtree)
{
    auto tw = msghandler_.tw();
    if (!tw && !jtree)
    {
        return;
    }
    // Get header names
    std::vector<std::string> headerNames{ "CLASS", "NAME" };
    for (const auto &pair : genome)
    {
        headerNames.push_back(gmx::formatString("BEST (%s)      ", iMolSelectName(pair.first).c_str()));
    }
    const std::vector<std::string> tmpHeaderNames{"MIN", "MAX", "MEAN", "STDEV", "MEDIAN" };
    headerNames.insert(headerNames.end(), tmpHeaderNames.begin(), tmpHeaderNames.end());
    // Get header sizes
    const int FLOAT_SIZE = 14;  // Adjusted for the %g formatting plus negative numbers
    std::vector<int> sizes{5, 4};
    for (size_t i = 0; i < genome.size(); i++)
    {
        sizes.push_back(FLOAT_SIZE);
    }
    std::vector<int> tmpSizes{FLOAT_SIZE, FLOAT_SIZE, FLOAT_SIZE, FLOAT_SIZE, FLOAT_SIZE };
    sizes.insert(sizes.end(), tmpSizes.begin(), tmpSizes.end());
    // Adjust size for class field
    const auto paramClass = sii_->paramClass();
    for (const auto &pClass : paramClass)
    {
        if (static_cast<int>(pClass.size()) > sizes[0])
        {
            sizes[0] = pClass.size();
        }
    }
    // Adjust size for name field
    const auto paramClassIndex = sii_->paramClassIndex();
    const auto paramNames = sii_->paramNamesWOClass();
    for (const auto &pName : paramNames)
    {
        if (static_cast<int>(pName.size()) > sizes[1])
        {
            sizes[1] = pName.size();
        }
    }
    const size_t TOTAL_WIDTH = static_cast<size_t>(std::accumulate(sizes.begin(), sizes.end(), 0)) + 3*(sizes.size()-1) + 4 + 4;
    const std::string HLINE(TOTAL_WIDTH, '-');
    // Print header
    if (tw)
    {
        tw->writeStringFormatted("%s\n|", HLINE.c_str());
        for (size_t i = 0; i < headerNames.size(); i++)
        {
            tw->writeStringFormatted(" %-*s |", sizes[i], headerNames[i].c_str());
        }
        tw->writeStringFormatted("\n%s\n", HLINE.c_str());
    }
    // Gather statistics from the population
    const std::vector<double> min    = pop.min();
    const std::vector<double> max    = pop.max();
    const std::vector<double> mean   = pop.mean();
    const std::vector<double> stdev  = pop.stdev();
    const std::vector<double> median = pop.median();
    auto optIndex = *sii_->optIndexPtr();

    // JsonTree with same information as log file
    JsonTree params("parameters");
    // Print the parameter information
    bool hitWall = false;
    for (size_t i = 0; i < paramClass.size(); i++)
    {
        for (size_t j = 0; j < paramNames.size(); j++)
        {
            if (paramClassIndex[j] != i)
            {
                continue;
            }
            if (tw)
            {
                tw->writeStringFormatted("| %-*s | %-*s |",
                                         sizes[0], paramClass[i].c_str(),
                                         sizes[1], paramNames[j].c_str());
            }
            JsonTree param(std::to_string(j));
            param.addObject(JsonTree("class", paramClass[i]));
            param.addObject(JsonTree("name", paramNames[j]));
            size_t k = 2;
            double value = 0;
            auto ffp = optIndex[j].forceFieldParameter();
            for (const auto &pair : genome)
            {
                value = pair.second.base(j);
                if (ffp->minimum() == value || ffp->maximum() == value)
                {
                    param.addObject("wall_hit", ffp->minimum() == value ? "lower" : "upper");
                    if (tw)
                    {
                        tw->writeStringFormatted(" %-*g (W) |", sizes[k], value);
                    }
                    hitWall = true;
                }
                else
                {
                    if (tw)
                    {
                        tw->writeStringFormatted(" %-*g     |", sizes[k], value);
                    }
                }
                param.addObject("best", std::to_string(value));
                k++;
            }
            if (tw)
            {
                param.addObject("min", std::to_string(min[j]));
                param.addObject("max", std::to_string(max[j]));
                param.addObject("mean", std::to_string(mean[j]));
                param.addObject("stdev", std::to_string(stdev[j]));
                param.addObject("median", std::to_string(median[j]));
                tw->writeStringFormatted(" %-*g | %-*g | %-*g | %-*g | %-*g |\n%s\n",
                                         sizes[k], min[j],
                                         sizes[k+1], max[j],
                                         sizes[k+2], mean[j],
                                         sizes[k+3], stdev[j],
                                         sizes[k+4], median[j],
                                         HLINE.c_str());
            }
            params.addObject(param);
        }
    }
    if (jtree)
    {
        jtree->addObject(params);
    }
    if (hitWall)
    {
        if (tw)
        {
            tw->writeStringFormatted("(W) indicates this parameter hit the wall.\n\n");
        }
    }
}

bool OptACM::runMaster(bool               optimize,
                       const std::string &sensitivityFile,
                       JsonTree          *jtree)
{
    GMX_RELEASE_ASSERT(commRec_.nodeType() == NodeType::Master,
                       "I thought I was the master...");

    msghandler_.writeDebug(memory_usage());
    auto tw = msghandler_.tw();
    bool bMinimum = false;
    std::map<iMolSelect, ga::Genome> bestGenome;
    if (optimize)
    {
        // Estimate number of fitness computations per dataset
        if (tw)
        {
            printNumCalcDevEstimate();
        }
        // Train the force field!
        bMinimum = ga_->evolve(&msghandler_, &bestGenome);
        if (tw)
        {
            tw->writeStringFormatted("\nHere are the best parameters found, together with some summary statistics of the last population:\n");
        }
        printGenomeTable(bestGenome, ga_->getLastPop(), jtree);
    }
    else
    {
        // This is to allow a sensitivity analysis with nooptimize. It is
        // done for the parameters of the force field, not random ones.
        if (gach_.randomInit())
        {
            msghandler_.msg(ACTStatus::Warning,
                            "Random initialization (-random_init) is turned off, since nothing is optimized. The parameters are those of the force field.");
        }
        ACMInitializer initializer(sii_, false, 0);
        auto ind = static_cast<alexandria::ACMIndividual *>(initializer.initialize());

        bestGenome[iMolSelect::Train] = ind->genome();
        delete ind;
        // The middlemen wait for a genome to start from, as in training,
        // before they can be given work or be stopped.
        for(auto &dest : commRec_.middlemen())
        {
            int read = 1;
            commRec_.send(dest, read);
            commRec_.send(dest, bestGenome[iMolSelect::Train].bases());
        }
        // They compute its fitness and send it back, see ACTMiddleMan::run.
        // Nothing is done with it, but it has to be received before the
        // middlemen are given work.
        for(auto &src : commRec_.middlemen())
        {
            ga::Genome genome;
            genome.Receive(&commRec_, src);
        }
    }
    auto sensitivity = gach_.sensitivity();
    if (gach_.optimizer() != OptimizerAlg::GA &&
        sensitivity != SensitivityMode::No && msghandler_.ok())
    {
        // Do sensitivity analysis only on the training set. The middlemen
        // wait for work and compute part of it.
        SensitivityAnalysis sens;
        (void) sens.run(&msghandler_, sii_, fitComp_,
                        &bestGenome[iMolSelect::Train], iMolSelect::Train, sensitivity,
                        false, true);
        const auto &res = sens.result();
        if (res.mode != SensitivityMode::No)
        {
            // The details go to their own file, train_ff.json gets a summary
            sens.write(sensitivityFile);
            msghandler_.write(gmx::formatString("Details of the sensitivity analysis are in %s",
                                                sensitivityFile.c_str()));
            if (jtree)
            {
                JsonTree summary("sensitivity");
                summary.addObject("mode", sensitivityModeToString(res.mode));
                summary.addObject("file", sensitivityFile);
                summary.addObject("minimum", res.minimum ? "yes" : "no");
                summary.addObject("n_parameters", std::to_string(res.nParameters));
                summary.addObject("n_parameters_converged", std::to_string(res.nParametersConverged));
                summary.addObject("n_parameters_flat", std::to_string(res.nParametersFlat));
                if (res.mode == SensitivityMode::TwoD)
                {
                    summary.addObject("n_pairs", std::to_string(res.nPairs));
                    summary.addObject("n_pairs_converged", std::to_string(res.nPairsConverged));
                }
                jtree->addObject(summary);
            }
        }
    }
    // Stop the middlemen ...
    if (commRec_.nmiddlemen() > 1)
    {
        for(auto &dest : commRec_.middlemen())
        {
            commRec_.send_done(dest);
        }
    }

    if (bMinimum && msghandler_.ok())
    {
        if (bestGenome.empty())
        {
            GMX_THROW(gmx::InternalError("Minimum found but the <Dataset, Genome> map is empty"));
        }
        if (tw)
        {
            for (auto it = bestGenome.begin(); it != bestGenome.end(); it++)
            {
                tw->writeString(it->second.print("Final best genome"));
                tw->writeStringFormatted("\nChi2 components of the best parameter vector found (for %s):\n", iMolSelectName(it->first).c_str());
                auto chi2before = it->second.fitness(it->first);
                // Print components
                fitComp_->compute(&msghandler_, &(it->second), it->first);
                // Check for consistency (read: bugs)
                double chi2 = it->second.fitness(it->first);
                // Shells relaxed from those of the previous evaluation, see
                // ForceComputer::setKeepShells, depend slightly on where they
                // started, in particular with poor parameters.
                double toler = bch_.shellKeep() ? 1e-3*std::max(1.0, std::abs(chi2before)) : 1e-3;
                if (std::abs(chi2before - chi2) > toler)
                {
                    msghandler_.msg(ACTStatus::Fatal,
                                    gmx::formatString("Chi2 after training was %.10g, after recomputing it is %.10g.", chi2before, chi2));
                }
                tw->writeStringFormatted("Minimum chi2 for %s %g\n",
                                         iMolSelectName(it->first).c_str(), chi2);
            }
        }
        // Save force field of best individual(s)
        {
            const std::set<int> changed;
            for (const auto &pair : bestGenome)
            {
                sii_->updateForceField(&msghandler_, changed, pair.second.bases());
                const auto myfilenm = outputFileName_[pair.first];
                tw->writeStringFormatted("Will save best %s force field to %s\n",
                                         iMolSelectName(pair.first).c_str(),
                                         myfilenm.c_str());
                sii_->saveState(true, myfilenm);
            }
            //! \todo resetting the train parameters for the TrainFFPrinter. We may have to work on that if we want to show the best test parameters too
            sii_->updateForceField(&msghandler_, changed, bestGenome.find(iMolSelect::Train)->second.bases());
        }
    }
    else if (optimize)
    {
        if (tw)
        {
            tw->writeString("Did not find a better parameter set\n");
        }
    }

    // Stop MASTER's helpers
    fitComp_->distributeTasks(CalcDev::Stop);
    // Turn off further communication
    sii_->commRecPtr()->done();

    // Final energy calculation for all molecules
    //! \todo parallellize this. There is no need to do this I believe, it's done above, and parallel!
    if (!bestGenome.empty() && msghandler_.ok())
    {
        const std::set<int> changed;
        auto ims = iMolSelect::Train;
        const auto bbb = bestGenome.find(ims);
        if (bestGenome.end() != bbb)
        {
            sii_->updateForceField(&msghandler_, changed, bbb->second.bases());
            fitComp_->calcDeviation(&msghandler_, CalcDev::ComputeAll, ims);
        }
        // Now compute the test compounds, with the best Train parameters.
        ims = iMolSelect::Test;
        fitComp_->calcDeviation(&msghandler_, CalcDev::ComputeAll, ims);
    }

    return bMinimum;
}

/*! \brief Tool to perform force field trainings
 * \param[in] argc Number of arguments on the cmd line
 * \param[in] argv The actual arguments
 * \return 0 if ok, 1 otherwise
 */
int train_ff(int argc, char *argv[])
{
    // Time the phases of the program, the parts of the training are
    // timed where they are done, see printTimings.
    PhaseTimer phase("Startup");
    std::vector<const char *> desc = {
        "train_ff reads a series of molecules and corresponding physico-chemical",
        "properties from a file. The properties can either originate from",
        "experimental data or from quantum chemistry calculations.",
        "The program then trains an empirical force field using",
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
        "[TT]P x M[tt].[PAR]",
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
        "Following training with MCMC or HYBRID a sensitivity analysis is done on",
        "the trained parameters, selected with [TT]-sensitivity[tt]: [TT]1D[tt]",
        "(default) changes each parameter on its own, [TT]2D[tt] also changes",
        "every pair of parameters, and [TT]No[tt] turns the analysis off. The",
        "results are written to the [TT]-sensout[tt] file, the log file gets a summary.",
        "It is also possible",
        "to do a sensitivity analysis on an existing force field. In that case",
        "please select either [TT]-fc_inter[tt] (Intermolecular energy) or [TT]-fc_epot[tt]",
        "(Intramolecular energy), set the parameters to be investigated with",
        "[TT]-fit 'a b c'[tt], and, finally, turn off optimization and random",
        "parameter initialization with the [TT]-nooptimize[tt] and",
        "[TT]-norandom_init[tt] flags.[PAR]"
    };

    bool                bcompress           = false;
    bool                bOptimize           = true;
    bool                bForceOutput        = false;

    gmx_output_env_t           *oenv;
    MolSelect                   gms;
    TrainForceFieldPrinter      printer;
    std::vector<t_pargs>        pargs;
    {
        t_pargs                     pa[]         = {
            { "-compress", FALSE, etBOOL, {&bcompress},
              "Compress output XML file" },
            { "-optimize",     FALSE, etBOOL, {&bOptimize},
              "Do parameter optimization when true, or a single calculation otherwise." }
        };

        for (int i = 0; i < asize(pa); i++)
        {
            pargs.push_back(pa[i]);
        }
    }
    std::vector<t_filenm>       filenms =
    {
        { efXML,  "-ff",      "aff",         ffRDMULT },
        { efXVG,  "-conv",    "param_conv" , ffWRITE  },
        { efXVG,  "-chi2",    "chi_squared", ffWRITE  },
        { efDAT,  "-fitness", "ga_fitness",  ffWRITE  },
        { efJSON, "-json",    "train_ff",    ffWRITE  },
        { efJSON, "-sensout", "sensitivity", ffWRITE  }
    };

    alexandria::OptACM opt;
    opt.add_options(&pargs, &filenms);
    CompoundReader compR;
    compR.addOptions(&pargs, &filenms, &desc);
    printer.addOptions(&pargs);
    printer.addFileOptions(&filenms);

    if (!parse_common_args(&argc,
                           argv,
                           PCA_CAN_VIEW,
                           filenms.size(), filenms.data(),
                           pargs.size(), pargs.data(),
                           desc.size(), desc.data(),
                           0,
                           nullptr,
                           &oenv))
    {
        return 1;
    }

    // Set output environment in the optimization driver
    opt.setOenv(oenv);

    // Check validity of arguments with check_pargs() in ConfigHandler(s)
    opt.check_pargs(nullptr);

    // Finishing MolGen stuff and setting output file for FF in OptACM.
    // Initializes commRec_ in opt
    // Fills id and prefix in sii
    // Calls optionsFinished() for MolGen instance.
    opt.optionsFinished(filenms);
    if (opt.commRec()->isIdle())
    {
        // Left over in the layout of the processes, see CommunicationRecord::init
        return 0;
    }

    // Check charge reader options
    compR.optionsFinished(opt.msgHandler(), filenms);
    if (!opt.msgHandler()->ok())
    {
        opt.msgHandler()->msg(ACTStatus::Error, "Please fix errors written above.");
        return 1;
    }
    auto tw = opt.msgHandler()->tw();
    if (opt.commRec()->isMaster())
    {
        if (opt.msgHandler())
        {
            opt.msgHandler()->writeDebug(memory_usage());
        }
        gms = compR.molselect();
        if (tw)
        {
            print_header(tw, pargs, filenms);
            tw->writeStringFormatted("Found %d Train and %d Test compounds in %s\n\n",
                                     gms.count(iMolSelect::Train), gms.count(iMolSelect::Test),
                                     opt2fn("-sel", filenms.size(), filenms.data()));
        }
        if (opt.msgHandler())
        {
            opt.msgHandler()->writeDebug(memory_usage());
        }
    }
    gms.bcast(opt.commRec());
    if (gms.count(iMolSelect::Test) > 0)
    {
        opt.sii()->fillFittingTargets(iMolSelect::Test);
    }
    else
    {
        opt.gach()->setEvaluateTestset(false);
        opt.bch()->setEvaluateTestset(false);
        if (opt.commRec()->isMaster() && tw)
        {
            tw->writeString("Turning off the evaluate of test set since it is empty.\n");
        }
    }
    // Propagate weights from training set to other sets
    opt.sii()->propagateWeightFittingTargets();

    // Read forcefield in StaticIndividualInfo sii_
    {
        const auto fns = opt2fns("-ff", filenms.size(), filenms.data());
        GMX_RELEASE_ASSERT(fns.size() == 1 || fns.size() == opt.gach()->popSize(),
                           gmx::formatString("Please pass exactly one or %d (popSize) force field file names", opt.gach()->popSize()).c_str());

        int fnIndex = 0;
        if (opt.commRec()->isMiddleMan() && fns.size() > 1)
        {
            fnIndex = opt.commRec()->middleManOrdinal();
        }
        opt.sii()->fillForceField(opt.msgHandler(), fns[fnIndex].c_str());
        if (tw)
        {
            tw->writeStringFormatted("On proc %d, found %d particle types\n",
                                     opt.sii()->commRec()->rank(),
                                     opt.sii()->forcefield()->nParticleTypes());
        }
    }
    if (opt.sii()->commRec()->isMaster())
    {
        for(const auto &ss : opt.sii()->forcefield()->info())
        {
            tw->writeLine(ss);
        }
    }
    // Now initiate the force computer
    opt.initForceComputer();
    // Check inputs
    bool optionsOk = true;
    if (opt.commRec()->isMaster())
    {
        opt.mg()->checkOptions(opt.msgHandler(), opt.sii()->forcefield(), compR.algorithm());
    }
    opt.commRec()->bcast(&optionsOk, opt.commRec()->comm_world());
    if (!optionsOk)
    {
        return 0;
    }
    // Read charges file if needed. No need to store the ACTmols being returned
    {
        std::vector<ACTMol> actmols;
        compR.read(opt.msgHandler(), *opt.sii()->forcefield(), opt.forceComp(), &actmols);
    }

    // MolGen read being called here!
    if (0 == opt.mg()->Read(opt.msgHandler(), filenms, opt.sii()->forcefield(), gms,
                            opt.sii()->fittingTargetsConst(iMolSelect::Train), &compR,
                            opt.forceComp()))
    {
        opt.msgHandler()->fatal("Training set is empty on one or more of the nodes, check your input and your command line flags. Rerun with -v 5 flag");
    }
    if (compR.algorithm() == ChargeGenerationAlgorithm::ESP)
    {
        // By now, we have initialized all the charges.
        // This means we can fix those, but also that no charge related parameters can be trained.
        //! \todo The latter needs to be implemented to prevent errors.
        compR.setAlgorithm(ChargeGenerationAlgorithm::Read);
    }

    // StaticIndividualInfo things
    {
        opt.sii()->generateOptimizationIndex(opt.msgHandler(), opt.mg(), opt.commRec());
        opt.sii()->fillVectors(opt.mg()->mindata(), opt.bch()->step());
        opt.sii()->computeWeightedTemperature(opt.bch()->temperatureWeighting());
    }
    // Set the output file names, has to be done before
    // creating a mutator.
    if (bOptimize)
    {
        if (opt.sii()->nParam() == 0)
        {
            if (opt.sii()->commRec()->isMaster())
            {
                opt.msgHandler()->msg(ACTStatus::Error,
                                      "No parameters to train. Did you run alexandria geometry_ff?");
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

    phase.next("Training (other)");
    int initOK   = 0;
    // Written at the end, after the timings have been added.
    JsonTree jtree("train_ff");
    if (opt.commRec()->isMaster())
    {
        if (opt.sii()->nParam() > 0)
        {
            initOK = opt.initMaster(filenms, compR.algorithm());
        }
        // Let the other nodes know whether all is well.
        if (opt.commRec()->isParallel())
        {
            opt.commRec()->bcast(&initOK, opt.commRec()->comm_world());
        }
        bool bMinimum = false;
        if (initOK)
        {
            bMinimum = opt.runMaster(bOptimize,
                                     opt2fn("-sensout", filenms.size(), filenms.data()),
                                     &jtree);
            if (bOptimize)
            {
                const auto msg = gmx::formatString("Training finished %s.", bMinimum ? "succesfully" : "without finding a better force field");
                opt.msgHandler()->msg(ACTStatus::Info, msg);
            }
        }
        phase.next("Output");
        if (bMinimum || bForceOutput || !bOptimize)
        {
            if (bForceOutput && !bMinimum)
            {
                opt.msgHandler()->msg(ACTStatus::Warning,
                                      gmx::formatString("No better minimum than the best initial candidate solution was found but -force_output was selected. This means that a global best force field file %s has been written written anyway.\n", opt2fn("-o", filenms.size(), filenms.data())));
                opt.sii()->saveState(true);
            }
            MolGen *tmpMg = opt.mg();
            printer.print(opt.msgHandler(), &jtree, opt.sii(), opt.forceComp(),
                          tmpMg->actmolsPtr(), oenv, filenms, !bOptimize);
            if (opt.msgHandler())
            {
                opt.msgHandler()->writeDebug(memory_usage());
            }
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
    }
    else
    {
        // Verify that all is well with the master.
        opt.commRec()->bcast(&initOK, opt.commRec()->comm_world());
        if (initOK > 0)
        {
            if (opt.commRec()->isMiddleMan())
            {
                if (opt.sii()->nParam() > 0)
                {
                    // Master and Individuals (middle-men) need to initialize more,
                    // so let's go.
                    ACTMiddleMan middleman(opt.msgHandler(), opt.mg(), opt.sii(), opt.gach(),
                                           opt.bch(), opt.oenv(), compR.algorithm(),
                                           opt.msgHandler()->verbose());
                    middleman.run(opt.msgHandler());
                }
            }
            else if (bOptimize || opt.gach()->sensitivity() != SensitivityMode::No)
            {
                if (opt.sii()->nParam() > 0)
                {
                    ACTHelper helper(opt.msgHandler(), opt.sii(), opt.mg(),
                                     opt.bch()->shellToler(),
                                     opt.bch()->shellMaxIter(),
                                     opt.bch()->shellKeep(),
                                     opt.bch()->gpuResident(),
                                     compR.algorithm(),
                                     opt.gach()->lossFunction());
                    helper.run(opt.msgHandler());
                }
            }
        }
    }
    phase.stop();
    auto json_file = opt.commRec()->isMaster() ? opt2fn_null("-json", filenms.size(), filenms.data()) : nullptr;
    const auto peak      = theoreticalPeakPerCore();
    printTimings(opt.commRec(), opt.commRec()->isMaster() ? opt.msgHandler()->tw() : nullptr,
                 phase.wallTime(), json_file ? &jtree : nullptr, &peak);
    if (json_file)
    {
        jtree.write(json_file, true);
    }
    return 0;
}

} // namespace alexandria
