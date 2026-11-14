/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
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

#ifndef ALEXANDRIA_SENSITIVITY_H
#define ALEXANDRIA_SENSITIVITY_H

#include <array>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "act/alexandria/sensitivity.h"
#include "act/basics/dataset.h"

namespace ga
{
class Genome;
}

namespace alexandria
{

class MsgHandler;
class StaticIndividualInfo;
class ACMFitnessComputer;

//! \brief Summary of a sensitivity analysis
struct SensitivityResult
{
    //! The analysis that was done, No if nothing was done
    SensitivityMode mode                 = SensitivityMode::No;
    //! The chi-squared of the parameters themselves
    double          chi2_0               = 0;
    //! Number of parameters
    size_t          nParameters          = 0;
    //! Number of converged parameters
    size_t          nParametersConverged = 0;
    //! Number of parameters that do not change the chi-squared, see Sensitivity::flat
    size_t          nParametersFlat      = 0;
    //! Number of pairs of parameters, zero in 1D
    size_t          nPairs               = 0;
    //! Number of converged pairs
    size_t          nPairsConverged      = 0;
    //! Number of pairs of each shape, indexed by PairShape
    std::array<size_t, 4> nPairShape     = { 0, 0, 0, 0 };
    //! Number of chi-squared computations, not counting chi2_0
    size_t          nEvaluations         = 0;
    //! Whether the parameters are in a minimum
    bool            minimum              = true;
    //! The full Hessian, only in 2D
    std::optional<HessianDiagnostic> hessian;
};

/*! \brief Sensitivity analysis of the chi-squared with respect to all
 * force field parameters, using a Sensitivity object for each parameter
 */
class SensitivityAnalysis
{
private:
    //! Dimensionless force constants
    std::vector<double>    forceConstant_;
    //! Summary of the last run
    SensitivityResult      result_;
    //! Machine readable output of the last run
    nlohmann::ordered_json json_;
public:
    //! Constructor
    SensitivityAnalysis() = default;
    /*!
     * \brief Perform a sensitivity analysis by systematically changing all parameters and
     * re-evaluating the \f$ \chi^2 \f$.
     * The results are stored for json() and write(), the log file gets
     * a summary only.
     * \param[in] msghandler The message and status handler
     * \param[in] sii        The force field and parameter information
     * \param[in] fitComp    Computes the chi-squared
     * \param[in] genome     Pointer to genome
     * \param[in] ims        Dataset to perform sensitivity analysis on
     * \param[in] mode       The analysis to do: 1D changes each parameter
     *                       on its own, 2D also each pair of parameters,
     *                       see Sensitivity2D. With No nothing is done.
     * \param[in] quiet      Do not print to log file
     * \param[in] useMiddlemen Whether the middlemen compute part of the
     *                       deviations for the training set. They have to
     *                       be waiting for work from the master, see
     *                       ACTMiddleMan::run.
     * \return true if all parameters are converged, that is the parabola
     *         through the chi2 of each of them curves upwards, and the
     *         parameter is within 5% of its range of the minimum of the
     *         parabola, limited to the bounds of the parameter. See
     *         Sensitivity::checkConvergence. In 2D all pairs have to be
     *         converged as well, see Sensitivity2D::checkConvergence.
     */
    bool run(MsgHandler           *msghandler,
             StaticIndividualInfo *sii,
             ACMFitnessComputer   *fitComp,
             ga::Genome           *genome,
             iMolSelect            ims,
             SensitivityMode       mode,
             bool                  quiet        = false,
             bool                  useMiddlemen = false);

    //! \return The force constants
    const std::vector<double> forceConstants() const { return forceConstant_; }

    //! \return The summary of the last run
    const SensitivityResult &result() const { return result_; }

    /*! \return The output of the last run, an object with the key
     * "sensitivity", or null if nothing was analysed
     */
    const nlohmann::ordered_json &json() const { return json_; }

    /*! \brief Write the output of the last run to a JSON file
     * \param[in] fileName The file name
     */
    void write(const std::string &fileName) const;
};

} // namespace

#endif
