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
#include "sensitivity_analysis.h"

#include "act/basics/msg_handler.h"
#include "act/utility/communicationrecord.h"
#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/actmiddleman.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/ga/genome.h"
#include "act/utility/timings.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/textwriter.h"

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

namespace alexandria
{

/*! \brief Fraction of the range of a parameter within which it has to be
 * from the minimum of the parabola to be converged, see
 * Sensitivity::checkConvergence
 */
static const double convergenceFraction = 0.05;

namespace
{

//! \brief A set of changed parameters and the resulting chi-squared
struct EvalPoint
{
    //! Index and new value of each changed parameter
    std::vector<std::pair<size_t, double>> changes;
    //! The chi-squared with these changes
    double                                 chi2 = 0;
};

/*! \brief The two values around a parameter at which to compute the chi-squared
 * Half a percent of the range below and above the parameter, or twice that
 * on one side when the other one is beyond a bound. The values stay apart,
 * also at a bound, and are in increasing order.
 * \param[in] p     The parameter
 * \param[in] lower Its lower bound
 * \param[in] upper Its upper bound
 * \return the two values
 */
std::pair<double, double> perturbedValues(double p, double lower, double upper)
{
    double deltap = (upper - lower)/200;
    if (p - deltap < lower)
    {
        return { p + deltap, p + 2*deltap };
    }
    else if (p + deltap > upper)
    {
        return { p - 2*deltap, p - deltap };
    }
    return { p - deltap, p + deltap };
}

/*! \brief Compute the chi-squared for each set of changed parameters
 * The points are done in batches of one for the master and one for each
 * worker. The master computes the first point of a batch itself, the
 * workers get a copy of the genome with the changes applied and compute
 * the chi-squared of the training set, see ACTMiddleMan::run. On return
 * the genome and the force field are as on entry.
 * \param[in]    msghandler The message and status handler
 * \param[in]    sii        The force field and parameter information
 * \param[in]    fitComp    Computes the chi-squared on the master
 * \param[inout] genome     The parameters, restored after each point
 * \param[in]    ims        The dataset
 * \param[in]    workers    The middlemen that compute part of the points
 * \param[inout] points     The changes, on return with their chi-squared
 */
void evaluate(MsgHandler              *msghandler,
              StaticIndividualInfo    *sii,
              ACMFitnessComputer      *fitComp,
              ga::Genome              *genome,
              iMolSelect               ims,
              const std::vector<int>  &workers,
              std::vector<EvalPoint>  *points)
{
    auto                 cr    = sii->commRec();
    auto                 cdc   = CalcDev::Compute;
    std::vector<double> *param = genome->basesPtr();
    for (size_t first = 0; first < points->size(); first += workers.size() + 1)
    {
        // The middlemen get the points after the first one
        size_t nsent = std::min(workers.size(), points->size() - first - 1);
        for (size_t w = 0; w < nsent; w++)
        {
            ga::Genome g = *genome;
            for (const auto &change : (*points)[first + 1 + w].changes)
            {
                (*g.basesPtr())[change.first] = change.second;
            }
            if (cr->queue())
            {
                // The random numbers and the generation, not used for a fitness
                sendQueueWork(cr, workers[w], { &g }, TrainFFMiddlemanMode::FITNESS, { 0 }, 0);
            }
            else
            {
                cr->send_data(workers[w]);
                g.Send(cr, workers[w]);
                cr->send(workers[w], TrainFFMiddlemanMode::FITNESS);
            }
        }
        // The master computes the first one
        auto               &point = (*points)[first];
        std::set<int>       changed;
        std::vector<double> pstore;
        for (const auto &change : point.changes)
        {
            changed.insert(static_cast<int>(change.first));
            pstore.push_back((*param)[change.first]);
            (*param)[change.first] = change.second;
        }
        sii->updateForceField(msghandler, changed, *param);
        fitComp->distributeTasks(cdc);
        point.chi2 = fitComp->calcDeviation(msghandler, cdc, ims);
        // Restore in reverse order, should an index occur more than once
        for (size_t c = point.changes.size(); c > 0; c--)
        {
            (*param)[point.changes[c-1].first] = pstore[c-1];
        }
        sii->updateForceField(msghandler, changed, *param);
        for (size_t w = 0; w < nsent; w++)
        {
            cr->recv(workers[w], &(*points)[first + 1 + w].chi2);
        }
    }
}

} // namespace

bool SensitivityAnalysis::run(MsgHandler           *msghandler,
                              StaticIndividualInfo *sii,
                              ACMFitnessComputer   *fitComp,
                              ga::Genome           *genome,
                              iMolSelect            ims,
                              SensitivityMode       mode,
                              bool                  quiet,
                              bool                  useMiddlemen)
{
    static const size_t timer = Timings::instance().index("Sensitivity analysis");
    const ScopedTimer   scoped(timer);
    std::vector<double> *param = genome->basesPtr();
    const auto upperBound      = sii->upperBound();
    const auto lowerBound      = sii->lowerBound();
    const auto paramNames      = sii->paramNames();

    // Reset the output of a previous run
    forceConstant_.clear();
    result_ = SensitivityResult();
    json_   = nlohmann::ordered_json();
    if (param->empty() || mode == SensitivityMode::No)
    {
        return true;
    }
    bool   twoD   = mode == SensitivityMode::TwoD;
    size_t np     = param->size();
    size_t npairs = twoD ? np*(np - 1)/2 : 0;
    const std::set<int> changed;
    sii->updateForceField(msghandler, changed, genome->bases());
    auto cdc    = CalcDev::Compute;
    fitComp->distributeTasks(cdc);
    auto chi2_0 = fitComp->calcDeviation(msghandler, cdc, ims);
    auto tw = msghandler->tw();
    if (quiet)
    {
        tw = nullptr;
    }
    // The chi2 of each parameter is known for its value, and computed
    // for two other values. In 2D, the chi2 is computed for the four
    // combinations of these values of each pair of parameters as well.
    // Each value is computed on its own, by the master or by one of the
    // middlemen, which only compute the training set.
    auto                cr       = sii->commRec();
    std::vector<int>    workers;
    if (useMiddlemen && iMolSelect::Train == ims && cr->isMaster())
    {
        workers = cr->middlemen();
    }
    std::vector<double>    pvalue(2*np);
    std::vector<EvalPoint> points(2*np + 4*npairs);
    for (size_t i = 0; i < np; ++i)
    {
        auto pv         = perturbedValues((*param)[i], lowerBound[i], upperBound[i]);
        pvalue[2*i]     = pv.first;
        pvalue[2*i + 1] = pv.second;
        points[2*i].changes     = { { i, pv.first } };
        points[2*i + 1].changes = { { i, pv.second } };
    }
    // The corners of each pair follow the points of the parameters
    size_t k = 2*np;
    for (size_t i = 0; twoD && i < np; ++i)
    {
        for (size_t j = i + 1; j < np; ++j)
        {
            for (size_t a = 0; a < 2; a++)
            {
                for (size_t b = 0; b < 2; b++)
                {
                    points[k++].changes = { { i, pvalue[2*i + a] }, { j, pvalue[2*j + b] } };
                }
            }
        }
    }
    evaluate(msghandler, sii, fitComp, genome, ims, workers, &points);

    bool   minimum    = true;
    auto   parameters = nlohmann::ordered_json::array();
    size_t nconv      = 0;
    size_t nflat      = 0;
    std::vector<Sensitivity> sens1d;
    std::vector<double>      range(np);
    for (size_t i = 0; i < np; ++i)
    {
        Sensitivity s;
        double      pstore = (*param)[i];
        // The points in increasing order of the parameter
        if (pvalue[2*i] > pstore)
        {
            s.add(pstore, chi2_0, true);
        }
        s.add(pvalue[2*i], points[2*i].chi2);
        if (pvalue[2*i] < pstore && pvalue[2*i + 1] > pstore)
        {
            s.add(pstore, chi2_0, true);
        }
        s.add(pvalue[2*i + 1], points[2*i + 1].chi2);
        if (pvalue[2*i + 1] < pstore)
        {
            s.add(pstore, chi2_0, true);
        }
        s.computeForceConstants(tw);
        range[i]     = upperBound[i] - lowerBound[i];
        bool   conv  = s.checkConvergence(lowerBound[i], upperBound[i], convergenceFraction*range[i]);
        minimum      = minimum && conv;
        if (conv)
        {
            nconv += 1;
        }
        if (s.flat())
        {
            nflat += 1;
        }
        // Compute dimensionless force constant, or zero if undefined
        double fc = 0;
        if (pstore != 0)
        {
            fc = s.a() / (pstore*pstore);
        }
        forceConstant_.push_back(fc);
        auto js = s.toJson(i, paramNames[i]);
        js["delta"]          = range[i]/200;
        js["force_constant"] = fc;
        parameters.push_back(js);
        sens1d.push_back(s);
    }

    // The pairs, and the full Hessian from all of them
    auto   pairs       = nlohmann::ordered_json::array();
    size_t npairsConv  = 0;
    std::array<size_t, 4> nshape = { 0, 0, 0, 0 };
    std::vector<std::vector<double>> hessian(np, std::vector<double>(np, 0.0));
    std::vector<double>              gradient(np);
    for (size_t i = 0; twoD && i < np; ++i)
    {
        hessian[i][i] = 2*sens1d[i].a();
        gradient[i]   = sens1d[i].derivative((*param)[i]);
    }
    k = 2*np;
    for (size_t i = 0; twoD && i < np; ++i)
    {
        for (size_t j = i + 1; j < np; ++j)
        {
            Sensitivity2D s2(sens1d[i], sens1d[j]);
            for (size_t a = 0; a < 2; a++)
            {
                for (size_t b = 0; b < 2; b++)
                {
                    s2.addCorner(pvalue[2*i + a], pvalue[2*j + b], points[k++].chi2);
                }
            }
            s2.computeCrossTerm();
            bool conv = s2.checkConvergence({ lowerBound[i], upperBound[i], lowerBound[j], upperBound[j] },
                                            { convergenceFraction*range[i], convergenceFraction*range[j] });
            minimum = minimum && conv;
            if (conv)
            {
                npairsConv += 1;
            }
            nshape[static_cast<size_t>(s2.shape())] += 1;
            hessian[i][j] = hessian[j][i] = s2.aij();
            pairs.push_back(s2.toJson(i, j, paramNames[i], paramNames[j]));
        }
    }
    result_.mode                 = mode;
    result_.chi2_0               = chi2_0;
    result_.nParameters          = np;
    result_.nParametersConverged = nconv;
    result_.nParametersFlat      = nflat;
    result_.nPairs               = npairs;
    result_.nPairsConverged      = npairsConv;
    result_.nPairShape           = nshape;
    result_.nEvaluations         = points.size();
    result_.minimum              = minimum;
    if (twoD)
    {
        result_.hessian = diagnoseHessian(hessian, gradient, range);
    }

    nlohmann::ordered_json sens;
    sens["mode"]                   = sensitivityModeToString(mode);
    sens["dataset"]                = iMolSelectName(ims);
    sens["chi2_0"]                 = chi2_0;
    sens["n_parameters"]           = np;
    if (twoD)
    {
        sens["n_pairs"]            = npairs;
    }
    sens["n_evaluations"]          = points.size();
    sens["tolerance_fraction"]     = convergenceFraction;
    sens["minimum"]                = minimum;
    sens["n_parameters_converged"] = nconv;
    sens["n_parameters_flat"]      = nflat;
    if (twoD)
    {
        sens["n_pairs_converged"]  = npairsConv;
        nlohmann::ordered_json shapes;
        for (auto shape : { PairShape::Minimum, PairShape::Saddle, PairShape::Maximum, PairShape::Degenerate })
        {
            shapes[pairShapeName(shape)] = nshape[static_cast<size_t>(shape)];
        }
        sens["pair_shapes"]        = shapes;
    }
    sens["parameters"]             = parameters;
    if (twoD)
    {
        sens["pairs"]              = pairs;
        sens["hessian_diagnostic"] = result_.hessian->toJson();
    }
    json_["sensitivity"]           = sens;
    if (tw)
    {
        if (twoD)
        {
            tw->writeStringFormatted("\nSensitivity analysis (%s) of %zu parameters and %zu pairs on the %s set, %zu chi2 evaluations, chi2_0 = %g\n",
                                     sensitivityModeToString(mode).c_str(), np, npairs,
                                     iMolSelectName(ims).c_str(), points.size(), chi2_0);
            tw->writeStringFormatted("Converged: %zu/%zu parameters, %zu/%zu pairs (%zu minimum, %zu saddle, %zu maximum, %zu degenerate), %s.\n",
                                     nconv, np, npairsConv, npairs,
                                     nshape[static_cast<size_t>(PairShape::Minimum)],
                                     nshape[static_cast<size_t>(PairShape::Saddle)],
                                     nshape[static_cast<size_t>(PairShape::Maximum)],
                                     nshape[static_cast<size_t>(PairShape::Degenerate)],
                                     minimum ? "in a minimum" : "not in a minimum");
            tw->writeStringFormatted("The Hessian of all parameters is %spositive definite.\n",
                                     result_.hessian->positiveDefinite ? "" : "not ");
        }
        else
        {
            tw->writeStringFormatted("\nSensitivity analysis (%s) of %zu parameters on the %s set, %zu chi2 evaluations, chi2_0 = %g\n",
                                     sensitivityModeToString(mode).c_str(), np,
                                     iMolSelectName(ims).c_str(), points.size(), chi2_0);
            tw->writeStringFormatted("Converged: %zu/%zu parameters, %s.\n", nconv, np,
                                     minimum ? "in a minimum" : "not in a minimum");
        }
        if (nflat > 0)
        {
            tw->writeStringFormatted("%zu of the %zu parameters do not change the chi2 (flat), for instance since they are zero in a force field that is not trained yet.%s\n",
                                     nflat, np, twoD ? " Their pairs are degenerate." : "");
        }
    }
    return minimum;
}

void SensitivityAnalysis::write(const std::string &fileName) const
{
    FILE *fp = gmx_ffopen(fileName, "w");
    fprintf(fp, "%s\n", json_.dump(2).c_str());
    gmx_ffclose(fp);
}

}
