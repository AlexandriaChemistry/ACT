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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <map>
#include <nlohmann/json.hpp>
#include <string>

#include "act/basics/dataset.h"
#include "act/ga/genome.h"
#include "gromacs/utility/futil.h"
#include "train_progress.h"

namespace alexandria
{

namespace
{

/*! \brief Make the object for a chi squared
 * \param[in] chi2     The chi squared for the data sets for which it is known
 * \param[in] withTest Whether to include the one for the test set
 * \return the object, with a key for each of the data sets train and test
 */
nlohmann::ordered_json chi2Json(const std::map<iMolSelect, double> &chi2,
                                bool                                withTest)
{
    nlohmann::ordered_json js = nlohmann::ordered_json::object();
    for (const auto &ims : { std::make_pair(iMolSelect::Train, "train"),
                             std::make_pair(iMolSelect::Test, "test") })
    {
        const auto c = chi2.find(ims.first);
        if (c != chi2.end() && (withTest || ims.first != iMolSelect::Test))
        {
            js[ims.second] = c->second;
        }
    }
    return js;
}

/*! \brief Write an object on one line and make sure it is on disk
 * \param[in] fp The file
 * \param[in] js The object
 */
void writeLine(FILE *fp, const nlohmann::ordered_json &js)
{
    // Replacing what is not a valid character or number
    fprintf(fp, "%s\n", js.dump(-1, ' ', false, nlohmann::json::error_handler_t::replace).c_str());
    fflush(fp);
}

} // namespace

TrainProgress::TrainProgress(const std::string &fileName,
                             const std::string &optimizer,
                             int                maxGenerations,
                             int                maxIter,
                             int                popSize,
                             size_t             nParam,
                             bool               testSet)
    : fileName_(fileName), optimizer_(optimizer), maxGenerations_(maxGenerations),
      maxIter_(maxIter), popSize_(popSize), nParam_(nParam), testSet_(testSet)
{
}

TrainProgress::~TrainProgress()
{
    if (fp_)
    {
        gmx_ffclose(fp_);
    }
}

std::map<iMolSelect, double> TrainProgress::chi2(const ga::Genome &genome)
{
    std::map<iMolSelect, double> chi2;
    for (const auto ims : { iMolSelect::Train, iMolSelect::Test })
    {
        if (genome.hasFitness(ims))
        {
            chi2[ims] = genome.fitness(ims);
        }
    }
    return chi2;
}

double TrainProgress::elapsed() const
{
    return std::chrono::duration<double>(Clock::now() - start_).count();
}

void TrainProgress::start(const std::map<iMolSelect, double> &chi2)
{
    if (fp_)
    {
        return;
    }
    fp_    = gmx_ffopen(fileName_, "w");
    start_ = Clock::now();
    nlohmann::ordered_json js;
    js["event"]        = "start";
    js["optimizer"]    = optimizer_;
    js["ngenerations"] = maxGenerations_;
    js["niterations"]  = maxIter_;
    js["population"]   = popSize_;
    js["nparameters"]  = nParam_;
    js["chi2"]         = chi2Json(chi2, testSet_);
    writeLine(fp_, js);
}

void TrainProgress::step(double                              iteration,
                         const std::map<iMolSelect, double> &chi2,
                         const std::map<iMolSelect, double> &best)
{
    if (!fp_)
    {
        return;
    }
    const auto now = Clock::now();
    if (haveStep_ && std::chrono::duration<double>(now - lastStep_).count() < stepInterval_)
    {
        return;
    }
    haveStep_ = true;
    lastStep_ = now;
    nlohmann::ordered_json js;
    js["event"]      = "step";
    js["generation"] = generation_;
    js["iteration"]  = iteration;
    // The chain has the test set only if it is computed at each step, and
    // then not for the best parameters
    js["chi2"]       = chi2Json(chi2, true);
    js["best"]       = chi2Json(best, false);
    js["elapsed"]    = { { "value", elapsed() }, { "unit", "s" } };
    writeLine(fp_, js);
}

void TrainProgress::generation(int                                 generation,
                               const std::map<iMolSelect, double> &best,
                               bool                                improved)
{
    if (!fp_)
    {
        return;
    }
    nlohmann::ordered_json js;
    js["event"]        = "generation";
    js["generation"]   = generation;
    js["ngenerations"] = maxGenerations_;
    js["best"]         = chi2Json(best, testSet_);
    js["improved"]     = improved;
    js["elapsed"]      = { { "value", elapsed() }, { "unit", "s" } };
    writeLine(fp_, js);
}

} // namespace alexandria
