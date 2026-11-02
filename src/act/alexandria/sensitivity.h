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

#include <string>
#include <vector>

namespace gmx
{
class TextWriter;
}

namespace alexandria
{

class JsonTree;

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
    //! \brief Whether checkConvergence was called
    bool   checked_   = false;
    //! \brief Whether the parameter is converged, see checkConvergence
    bool   converged_ = false;
public:
    //! \brief Constructor
    Sensitivity() {}

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

    //! Return the constants after computation
    double a() const { return a_; }
    double b() const { return b_; }
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

    /*! \brief Print output
     * \param[in] tw    A text writer
     * \param[in] jtree For machine readable output
     * \param[in] index Parameter index, must be unique
     * \param[in] label Label for identifying the parameter
     */
    void print(gmx::TextWriter      *tw,
               alexandria::JsonTree *jtree,
               const std::string    &index,
               const std::string    &label);
};

}  //namespace alexandria

#endif //ALEXANDRIA_BAYES_H
