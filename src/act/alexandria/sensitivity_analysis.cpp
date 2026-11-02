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
            cr->send_data(workers[w]);
            g.Send(cr, workers[w]);
            cr->send(workers[w], TrainFFMiddlemanMode::FITNESS);
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
    ScopedTimer         scoped(timer);
    std::vector<double> *param = genome->basesPtr();
    const auto upperBound      = sii->upperBound();
    const auto lowerBound      = sii->lowerBound();
    const auto paramNames      = sii->paramNames();

    // Reset the output of a previous run
    forceConstant_.clear();
    result_ = SensitivityResult();
    json_   = nlohmann::ordered_json();
    if (param->size() == 0 || mode == SensitivityMode::No)
    {
        return true;
    }
    if (mode == SensitivityMode::TwoD)
    {
        //! \todo Implement the two-dimensional analysis
        msghandler->msg(ACTStatus::Warning,
                        "The two-dimensional sensitivity analysis is not implemented yet, doing the one-dimensional one.");
        mode = SensitivityMode::OneD;
    }
    std::set<int> changed;
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
    // for two other values, all of them first. Each value is computed on
    // its own, by the master or by one of the middlemen, which only
    // compute the training set.
    auto                cr       = sii->commRec();
    std::vector<int>    workers;
    if (useMiddlemen && iMolSelect::Train == ims && cr->isMaster())
    {
        workers = cr->middlemen();
    }
    size_t                 np = param->size();
    std::vector<double>    pvalue(2*np), chi2(2*np);
    std::vector<EvalPoint> points(2*np);
    for (size_t i = 0; i < np; ++i)
    {
        auto pv         = perturbedValues((*param)[i], lowerBound[i], upperBound[i]);
        pvalue[2*i]     = pv.first;
        pvalue[2*i + 1] = pv.second;
        points[2*i].changes     = { { i, pv.first } };
        points[2*i + 1].changes = { { i, pv.second } };
    }
    evaluate(msghandler, sii, fitComp, genome, ims, workers, &points);
    for (size_t k = 0; k < points.size(); ++k)
    {
        chi2[k] = points[k].chi2;
    }
    bool minimum   = true;
    auto parameters = nlohmann::ordered_json::array();
    size_t nconv    = 0;
    for (size_t i = 0; i < np; ++i)
    {
        Sensitivity s;
        double      pstore = (*param)[i];
        // The points in increasing order of the parameter
        if (pvalue[2*i] > pstore)
        {
            s.add(pstore, chi2_0, true);
        }
        s.add(pvalue[2*i], chi2[2*i]);
        if (pvalue[2*i] < pstore && pvalue[2*i + 1] > pstore)
        {
            s.add(pstore, chi2_0, true);
        }
        s.add(pvalue[2*i + 1], chi2[2*i + 1]);
        if (pvalue[2*i + 1] < pstore)
        {
            s.add(pstore, chi2_0, true);
        }
        s.computeForceConstants(tw);
        double range = upperBound[i] - lowerBound[i];
        bool   conv  = s.checkConvergence(lowerBound[i], upperBound[i], convergenceFraction*range);
        minimum      = minimum && conv;
        if (conv)
        {
            nconv += 1;
        }
        // Compute dimensionless force constant, or zero if undefined
        double fc = 0;
        if (pstore != 0)
        {
            fc = s.a() / (pstore*pstore);
        }
        forceConstant_.push_back(fc);
        auto js = s.toJson(i, paramNames[i]);
        js["delta"]          = range/200;
        js["force_constant"] = fc;
        parameters.push_back(js);
    }
    result_.mode                 = mode;
    result_.chi2_0               = chi2_0;
    result_.nParameters          = np;
    result_.nParametersConverged = nconv;
    result_.nEvaluations         = points.size();
    result_.minimum              = minimum;

    nlohmann::ordered_json sens;
    sens["mode"]                   = sensitivityModeToString(mode);
    sens["dataset"]                = iMolSelectName(ims);
    sens["chi2_0"]                 = chi2_0;
    sens["n_parameters"]           = np;
    sens["n_evaluations"]          = points.size();
    sens["tolerance_fraction"]     = convergenceFraction;
    sens["minimum"]                = minimum;
    sens["n_parameters_converged"] = nconv;
    sens["parameters"]             = parameters;
    json_["sensitivity"]           = sens;
    if (tw)
    {
        tw->writeStringFormatted("\nSensitivity analysis (%s) of %zu parameters on the %s set, %zu chi2 evaluations, chi2_0 = %g\n",
                                 sensitivityModeToString(mode).c_str(), np,
                                 iMolSelectName(ims).c_str(), points.size(), chi2_0);
        tw->writeStringFormatted("Converged: %zu/%zu parameters, %s.\n", nconv, np,
                                 minimum ? "in a minimum" : "not in a minimum");
    }
    return minimum;
}

void SensitivityAnalysis::write(const std::string &fileName) const
{
    FILE *fp = gmx_ffopen(fileName.c_str(), "w");
    fprintf(fp, "%s\n", json_.dump(2).c_str());
    gmx_ffclose(fp);
}

}
