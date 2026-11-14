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
 */

#ifndef ALEXANDRIA_BAYES_H
#define ALEXANDRIA_BAYES_H

#include <array>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json_fwd.hpp>

namespace gmx
{
class TextWriter;
}

namespace alexandria
{

//! \brief Which sensitivity analysis to do, see SensitivityAnalysis
enum class SensitivityMode {
    //! No sensitivity analysis
    No,
    //! Each parameter on its own
    OneD,
    //! Each parameter on its own and every pair of parameters
    TwoD
};

/*! \brief Convert string into SensitivityMode
 * \param[in] str The string, No, 1D or 2D
 * \return the sensitivity mode
 * \throws if there is no mode matching the string
 */
SensitivityMode stringToSensitivityMode(const std::string &str);

/*! \brief Convert SensitivityMode to string
 * \param[in] mode The sensitivity mode
 * \return The corresponding string, No, 1D or 2D
 */
const std::string &sensitivityModeToString(SensitivityMode mode);

//! How to perform the calculation of deviations (chi-squared)
enum class CalcDev {
    //! Do a calculation local
    Compute = 7,
    //! Compute all the deviations
    ComputeAll = 13,
    //! Distribute parameters
    Parameters = 17,
    //! Ready to stop
    Stop = 23,
};

const char *calcDevName(CalcDev cd);

/*! \brief Sensitivity of the chi-squared to one force field parameter,
 * determined by fitting a parabola chi2 = a p^2 + b p + c through
 * chi-squared values computed for different values p of the parameter
 */
class Sensitivity
{
private:
    //! \brief parameter values used
    std::vector<double> p_;
    //! \brief chi2 values obtains
    std::vector<double> chi2_;
    //! \brief Constants for the parabola fitting
    double a_ = 0, b_ = 0, c_ = 0;
    /*! \brief The parabola in reduced coordinates u = (p - p0_)/scale_
     * chi2 = A_ u^2 + B_ u + C_. The fit is done in these coordinates,
     * since the fit in p is ill-conditioned when the points are close
     * together compared to their value.
     */
    double A_ = 0, B_ = 0, C_ = 0;
    //! \brief Centre of the reduced coordinates
    double p0_    = 0;
    //! \brief Scale of the reduced coordinates
    double scale_ = 1;
    //! \brief Index of the point with the original parameter value
    size_t original_ = 1;
    //! \brief Whether the fit succeeded
    bool   fitted_    = false;
    //! \brief Whether checkConvergence was called
    bool   checked_   = false;
    //! \brief The bounds given to checkConvergence
    double lower_ = 0, upper_ = 0;
    //! \brief Whether the parameter is converged, see checkConvergence
    bool   converged_ = false;
public:
    //! \brief Constructor
    Sensitivity() = default;

    /*! \brief
     * Add a point
     * \param[in] p        The parameter value
     * \param[in] chi2     The chi-squared value
     * \param[in] original Whether this is the original parameter value,
     *                     otherwise it is the second point
     */
    void add(double p, double chi2, bool original = false)
    {
        if (original)
        {
            original_ = p_.size();
        }
        p_.push_back(p);
        chi2_.push_back(chi2);
    }
    /*! \brief
     * Compute the fit to the curve
     * The parabola is fitted in reduced coordinates around the original
     * parameter value, or the mean if there is none, and converted to
     * a p^2 + b p + c.
     * \param[in] tw A text writer
     */
    void computeForceConstants(gmx::TextWriter *tw);

    /*! \brief Determine whether the original parameter value is converged
     * This is the case if the parabola curves upwards, and the parameter
     * is within the tolerance of the minimum of the parabola. The minimum
     * is limited to the bounds, such that a parameter at a bound is
     * converged if the minimum is beyond that bound.
     * Call computeForceConstants first.
     * \param[in] lower     Lower bound of the parameter
     * \param[in] upper     Upper bound of the parameter
     * \param[in] tolerance Largest distance to the minimum
     * \return whether the parameter is converged
     */
    bool checkConvergence(double lower, double upper, double tolerance);

    //! \return the quadratic coefficient a of the parabola after computation
    double a() const { return a_; }
    //! \return the linear coefficient b of the parabola after computation
    double b() const { return b_; }
    //! \return the constant c of the parabola after computation
    double c() const { return c_; }

    /*! \brief The extremum of the parabola, computed in reduced coordinates
     * \return the parameter value at the extremum, or zero if the parabola
     *         is a straight line
     */
    double pMin() const;

    /*! \brief The chi-squared at the extremum of the parabola
     * \return the chi-squared at pMin(), or zero if the parabola is a
     *         straight line
     */
    double chi2Min() const;

    //! \return the parameter value of the minimum limited to the bounds of checkConvergence
    double pMinBounded() const;

    //! \return whether the parabola was fitted
    bool fitted() const { return fitted_; }

    /*! \return whether the chi-squared is the same at all points, within
     * its rounding errors, such that it does not depend on the parameter
     * here. This is the case for instance for a Lennard-Jones sigma when
     * epsilon is zero, as in a force field that is not trained yet.
     */
    bool flat() const;

    //! \return the parameter values, in the order they were added
    const std::vector<double> &pValues() const { return p_; }

    //! \return the chi-squared values, in the order they were added
    const std::vector<double> &chi2Values() const { return chi2_; }

    //! \return the index of the original parameter value
    size_t original() const { return original_; }

    /*! \brief The parabola at a parameter value, computed in reduced coordinates
     * \param[in] p The parameter value
     * \return the chi-squared of the fitted parabola, or zero without a fit
     */
    double value(double p) const;

    /*! \brief The derivative of the parabola, computed in reduced coordinates
     * \param[in] p The parameter value
     * \return the derivative of the fitted parabola, or zero without a fit
     */
    double derivative(double p) const;

    /*! \brief The results in machine readable form
     * Numbers are JSON numbers and flags are booleans. Values that are
     * undefined, like the minimum of a straight line, are null.
     * \param[in] index Parameter index
     * \param[in] label Label for identifying the parameter, the particle
     *                  type and the parameter name separated by a space
     * \return a JSON object
     */
    nlohmann::ordered_json toJson(size_t             index,
                                  const std::string &label) const;
};

//! \brief The shape of the chi-squared around a pair of parameters
enum class PairShape {
    //! Both curvatures and the determinant of the Hessian are positive
    Minimum,
    //! The determinant of the Hessian is negative
    Saddle,
    //! Both curvatures are negative and the determinant is positive
    Maximum,
    //! A curvature or the determinant is zero within the precision
    Degenerate
};

//! \return the name of a PairShape, minimum, saddle, maximum or degenerate
const std::string &pairShapeName(PairShape shape);

/*! \brief The chi-squared as a function of two parameters
 * The surface is the sum of the two one-dimensional parabolas plus a
 * coupling term,
 * chi2(x, y) = P_i(x) + P_j(y) - chi2_0 + a_ij (x - p_i) (y - p_j),
 * where p_i and p_j are the original parameter values and chi2_0 the
 * chi-squared there. The coupling a_ij is fitted to the corners of the
 * grid, the points where both parameters differ from the original values.
 * For a symmetric grid it is the mixed central difference, and it is
 * exact for a quadratic surface.
 */
class Sensitivity2D
{
private:
    //! \brief The one-dimensional analyses of the two parameters
    Sensitivity si_, sj_;
    //! \brief The corners: both parameter values and the chi-squared
    std::vector<std::array<double, 3>> corners_;
    //! \brief The coupling term
    double aij_       = 0;
    //! \brief Root mean square deviation of the corners from the surface
    double cornerRms_ = 0;
    //! \brief Whether the coupling could be computed
    bool   computed_  = false;
    //! \brief Whether checkConvergence was called
    bool   checked_   = false;
    //! \brief Whether the pair is converged, see checkConvergence
    bool   converged_ = false;
    //! \brief The bounds given to checkConvergence, lower and upper of i and of j
    std::array<double, 4> bounds_ = { 0, 0, 0, 0 };
public:
    /*! \brief Constructor
     * \param[in] si The fitted one-dimensional analysis of the first parameter
     * \param[in] sj The same for the second parameter
     */
    Sensitivity2D(const Sensitivity &si, const Sensitivity &sj) : si_(si), sj_(sj) {}

    /*! \brief Add a corner
     * \param[in] pi   Value of the first parameter
     * \param[in] pj   Value of the second parameter
     * \param[in] chi2 The chi-squared
     */
    void addCorner(double pi, double pj, double chi2) { corners_.push_back({ pi, pj, chi2 }); }

    //! \brief Compute the coupling term from the corners
    void computeCrossTerm();

    //! \return the coupling term a_ij
    double aij() const { return aij_; }

    //! \return the root mean square deviation of the corners from the surface
    double cornerRms() const { return cornerRms_; }

    /*! \return the dimensionless coupling a_ij/(2 sqrt(a_i a_j)), which is
     * in (-1, 1) for a minimum, or nothing if a curvature is not positive
     */
    std::optional<double> rho() const;

    //! \return the Hessian at the original parameter values
    std::array<std::array<double, 2>, 2> hessian() const;

    //! \return the gradient at the original parameter values
    std::array<double, 2> gradient() const;

    //! \return the shape of the surface
    PairShape shape() const;

    /*! \return the stationary point of the surface, or nothing if the
     * Hessian is singular
     */
    std::optional<std::array<double, 2>> stationaryPoint() const;

    /*! \brief The minimum of the surface within bounds
     * Only defined for a minimum. If the stationary point is outside the
     * bounds, the minimum is on one of the edges, where the surface is
     * a parabola in the other parameter.
     * \param[in] bounds Lower and upper bound of the first and of the second parameter
     * \return the minimum, or nothing if the surface is not a minimum
     */
    std::optional<std::array<double, 2>> boxMinimum(const std::array<double, 4> &bounds) const;

    /*! \brief The surface at parameter values
     * \param[in] x Value of the first parameter
     * \param[in] y Value of the second parameter
     * \return the chi-squared of the surface
     */
    double value(double x, double y) const;

    /*! \brief Determine whether the pair is converged
     * This is the case if the surface is a minimum, and both parameters
     * are within their tolerance of the minimum within the bounds.
     * \param[in] bounds     Lower and upper bound of the first and of the second parameter
     * \param[in] tolerances Largest distance to the minimum of the first and the second parameter
     * \return whether the pair is converged
     */
    bool checkConvergence(const std::array<double, 4> &bounds,
                          const std::array<double, 2> &tolerances);

    /*! \brief The results in machine readable form
     * \param[in] i     Index of the first parameter
     * \param[in] j     Index of the second parameter
     * \param[in] namei Name of the first parameter
     * \param[in] namej Name of the second parameter
     * \return a JSON object
     */
    nlohmann::ordered_json toJson(size_t             i,
                                  size_t             j,
                                  const std::string &namei,
                                  const std::string &namej) const;
};

/*! \brief Properties of the full Hessian of the chi-squared
 * The Hessian is scaled by the ranges of the parameters, H_ij r_i r_j,
 * such that its eigenvalues do not depend on the units of the parameters.
 * Positive definite pairs do not make a positive definite Hessian, so this
 * is a stricter test for a minimum than that of all pairs.
 */
struct HessianDiagnostic
{
    //! Whether the Hessian is positive definite
    bool                  positiveDefinite = false;
    //! Smallest eigenvalue of the scaled Hessian
    double                eigenvalueMin    = 0;
    //! Largest eigenvalue of the scaled Hessian
    double                eigenvalueMax    = 0;
    /*! The largest component of the Newton step to the minimum, relative
     * to the range of the parameter, only for a positive definite Hessian.
     * The bounds are not taken into account.
     */
    std::optional<double> maxNewtonStepFraction;

    //! \return the diagnostic in machine readable form
    nlohmann::ordered_json toJson() const;
};

/*! \brief Diagnose the full Hessian of the chi-squared
 * \param[in] hessian  The symmetric Hessian, N x N
 * \param[in] gradient The gradient, N
 * \param[in] range    The range of each parameter, N, all positive
 * \return the diagnostic
 */
HessianDiagnostic diagnoseHessian(const std::vector<std::vector<double>> &hessian,
                                  const std::vector<double>              &gradient,
                                  const std::vector<double>              &range);

}  //namespace alexandria

#endif //ALEXANDRIA_BAYES_H
