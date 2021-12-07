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

#include "tune_eem.h"

#include "actpre.h"

#include "tune_eem.h"

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

#include "alex_modules.h"
#include "gentop_core.h"
#include "gmx_simple_comm.h"
#include "memory_check.h"
#include "molgen.h"
#include "molprop_util.h"
#include "mymol_low.h"
#include "bayes.h"
#include "poldata.h"
#include "poldata_tables.h"
#include "poldata_xml.h"
#include "tuning_utility.h"
#include "units.h"

namespace alexandria
{

/*! \defgroup tune_eem Schematic flowchart for tune_eem
 *
 * This diagram shows how the program is parallelized.
 * Boxes in cyan run on all nodes, pink on the master node
 * and yellow on the helper nodes.
 * \dot
digraph tune_eem {
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
                    {"-fullQuadrupole", FALSE, etBOOL, {&bFullQuadrupole_},
                            "Consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization"},
                    {"-removemol",      FALSE, etBOOL, {&bRemoveMol_},
                            "Remove a molecule from training set if shell minimization does not converge."},
            };
    for (int i = 0; i < asize(pa); i++) {
        pargs->push_back(pa[i]);
    }
    addOptions(pargs, eTune::EEM);
    configHandlerPtr()->add_pargs(pargs);
}

void OptACM::optionsFinished(const std::string &outputFile) {
    MolGen::optionsFinished();
    outputFile_ = outputFile;
}

void OptACM::openLogFile(const char *logfileName) {
    fplog_.reset(gmx_ffopen(logfileName, "w"));
}

FILE *OptACM::logFile() {
    if (fplog_) {
        return fplog_.get();
    } else {
        return nullptr;
    }
}

void OptACM::initChargeGeneration(iMolSelect ims)
{
    std::string method, basis, conf, type, myref, mylot;
    splitLot(lot(), &method, &basis);
    tensor              polar      = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<double> vec;
    for (MyMol &mymol : mymols())
    {
        if (mymol.datasetType() != ims)
        {
            continue;
        }
        if (fit("alpha"))
        {
            // For fitting alpha we need a reference polarizability
            double ref_pol, error, T = 0;
            if (mymol.getPropRef(MolPropObservable::POLARIZABILITY, iqmType::QM,
                                 method, basis, "",
                                 (char *)"electronic",
                                 &ref_pol, &error, &T,
                                 &myref, &mylot, &vec, polar))
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

void OptACM::initOpt(bool bRandom)
{
    for(auto &optIndex : optIndex_)
    {
        auto                iType = optIndex.iType();
        ForceFieldParameter p;
        if (iType == InteractionType::CHARGE)
        {
            if (poldata()->hasParticleType(optIndex.particleType()))
            {
                p = poldata()->findParticleType(optIndex.particleType())->parameterConst(optIndex.parameterType());
            }
        }
        else if (poldata()->interactionPresent(iType))
        {
            p = poldata()->findForcesConst(iType).findParameterTypeConst(optIndex.id(), optIndex.parameterType());
        }
        if (p.ntrain() >= mindata())
        {
            Bayes::addParam(optIndex.name(),
                            p.value(), p.mutability(),
                            p.minimum(), p.maximum(),
                            p.ntrain(), bRandom);
        }
    }
}

bool OptACM::runMaster(const gmx_output_env_t *oenv,
                       const char             *xvgconv,
                       const char             *xvgepot,
                       bool                    optimize,
                       bool                    sensitivity,
                       bool 		           bEvaluate_testset)
{
    bool bMinimum = false;
    GMX_RELEASE_ASSERT(MASTER(commrec()), "WTF");

    print_memory_usage(debug);
    std::vector<std::string> paramClass;
    for(const auto &fm : typesToFit())
    {
        paramClass.push_back(fm.first);
    }
    if (optimize)
    {
        configHandlerPtr()->setOutputFiles(xvgconv, paramClass, xvgepot, oenv);
        double chi2     = 0;
        bMinimum = Bayes::MCMC(logFile(), bEvaluate_testset, &chi2);
    }
    if (sensitivity)
    {
        // only on the training set
        Bayes::SensitivityAnalysis(logFile(), iMolSelect::Train);
    }
    // Finalize the calculations on the helpers
    GMX_RELEASE_ASSERT(calcDeviation(false, CalcDev::Final, iMolSelect::Train) < 0,
                       "Result for final parallel calcDeviation should be less than zero");

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
    else if (optimize)
    {
        fprintf(logFile(), "Did not find a better parameter set\n");
    }
    return bMinimum;
}

void OptACM::runHelper()
{
    // H E L P E R   N O D E S
    // The second and third variable are set by the master, but
    // we have to pass something.
    // If the result is less than zero (-1), we are done.
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
        { efXVG, "-conv",      "param-conv",    ffWRITE },
        { efXVG, "-epot",      "param-epot",    ffWRITE }
    };
    const int           NFILE               = asize(fnm);

    real                efield              = 10;
    bool                bRandom             = false;
    bool                bcompress           = false;
    bool                bZero               = true;
    bool                bOptimize           = true;
    bool                bSensitivity        = true;
    bool                bForceOutput        = false;
    bool                bEvaluate_testset   = false;

    t_pargs                     pa[]         = {
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
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

    gmx_output_env_t           *oenv;
    MolSelect                   gms;
    TuneForceFieldPrinter       printer;

    std::vector<t_pargs>        pargs;
    for (int i = 0; i < asize(pa); i++)
    {
        pargs.push_back(pa[i]);
    }
    alexandria::OptACM opt;
    opt.add_pargs(&pargs);
    printer.addOptions(&pargs);

    std::vector<t_filenm>       filenms;
    for(int i = 0; i < asize(fnm); i++)
    {
        filenms.push_back(fnm[i]);
    }
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

    // TODO: Check validity of arguments with check_pargs() in ConfigHandler(s)
    opt.configHandlerPtr()->check_pargs();

    opt.optionsFinished(opt2fn("-o", NFILE, fnm));

    opt.fillDevComputers();

    if (MASTER(opt.commrec()))
    {
        opt.openLogFile(opt2fn("-g", NFILE, fnm));
        print_memory_usage(debug);
        print_header(opt.logFile(), pargs);
        gms.read(opt2fn_null("-sel", NFILE, fnm));
        fprintf(opt.logFile(), "Found %d Train and %d Test compounds in %s\n\n",
                gms.count(iMolSelect::Train), gms.count(iMolSelect::Test),
                opt2fn("-sel", NFILE, fnm));
        print_memory_usage(debug);
    }

    // MolGen read being called here!
    if (0 == opt.Read(opt.logFile() ? opt.logFile() : (debug ? debug : nullptr),
                      opt2fn("-f", NFILE, fnm),
                      opt2fn_null("-d", NFILE, fnm),
                      bZero,
                      gms,
                      false,
                      false,
                      opt2fn_null("-table", NFILE, fnm),
                      opt.verbose()))
    {
        if (opt.logFile())
        {
            fprintf(opt.logFile(), "Training set is empty, check your input. Rerun with -v option or -debug 1.\n");
        }
        return 0;
    }
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

    if (MASTER(opt.commrec()))
    {
        if (bOptimize || bSensitivity)
        {
            opt.initOpt(bRandom);
        }
        bool bMinimum = opt.runMaster(oenv,
                                      opt2fn("-conv", NFILE, fnm),
                                      opt2fn("-epot", NFILE, fnm),
                                      bOptimize,
                                      bSensitivity,
                                      bEvaluate_testset);

        if (bMinimum || bForceOutput || !bOptimize)
        {
            if (bForceOutput)
            {
                fprintf(opt.logFile(), "Output based on last step of MC simulation per your specification.\nUse the -noforce_output flag to prevent this.\nThe force field output file %s is based on the last MC step as well.\n", opt2fn("-o", NFILE, fnm));
                opt.saveState();
            }
            printer.print(opt.logFile(),
                          &(opt.mymols()),
                          opt.poldata(),
                          opt.mdlog(),
                          opt.lot(),
                          opt.qcycle(),
                          opt.qtol(),
                          oenv,
                          opt.fullQuadrupole(),
                          opt.commrec(),
                          efield,
                          filenms);
            print_memory_usage(debug);
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
    }
    else if (bOptimize || bSensitivity)
    {
        opt.runHelper();
    }
    return 0;
}
