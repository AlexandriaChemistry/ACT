/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
 * Tests for the Sensitivity class and the SensitivityAnalysis routine.
 * The analysis is done for the electronegativities of the atom types in
 * a training set of butanol, acetone and water, with the electrostatic
 * potential as the fitting target.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <string>
#include <vector>

#include <Eigen/Eigenvalues>
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/acmindividual.h"
#include "act/alexandria/acminitializer.h"
#include "act/alexandria/molgen.h"
#include "act/alexandria/molselect.h"
#include "act/alexandria/sensitivity_analysis.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/basics/msg_handler.h"
#include "act/forces/forcecomputer.h"
#include "act/ga/genome.h"
#include "act/import/compound_reader.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/filetypes.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

// Test: With an empty genome there is nothing to analyse, the routine
// reports a minimum and adds nothing to the output.
TEST(SensitivityAnalysisTest, EmptyGenome)
{
    MsgHandler           msghandler;
    CommunicationRecord  cr;
    cr.init(cr.size());
    StaticIndividualInfo sii(&cr);
    ga::Genome           genome;
    SensitivityAnalysis  sens;
    for (auto mode : { SensitivityMode::No, SensitivityMode::OneD, SensitivityMode::TwoD })
    {
        EXPECT_TRUE(sens.run(&msghandler, &sii, nullptr, &genome, iMolSelect::Train, mode, true));
        EXPECT_TRUE(sens.json().is_null());
        EXPECT_TRUE(sens.forceConstants().empty());
        EXPECT_EQ(SensitivityMode::No, sens.result().mode);
    }
}

/*! \brief Sensitivity analysis of a small training set
 * The electronegativities of the atom types of butanol, acetone and water
 * are the parameters, the deviation from the electrostatic potentials of
 * these molecules is the chi-squared.
 */
class SensitivityAnalysisTrainTest : public gmx::test::CommandLineTestBase
{
protected:
    //! Messages
    MsgHandler                     msghandler_;
    //! Communication, a single process
    CommunicationRecord            cr_;
    //! The parameters and the force field
    StaticIndividualInfo          *sii_    = nullptr;
    //! The molecules
    MolGen                        *molgen_ = nullptr;
    //! Computes energies and forces
    ForceComputer                  forceComp_;
    //! Computes the chi-squared
    ACMFitnessComputer             fitComp_;
    //! The parameters of the force field
    ga::Genome                     genome_;

    //! Set up the training set
    SensitivityAnalysisTrainTest()
    {
        msghandler_.setPrintLevel(ACTStatus::Warning);
        cr_.init(cr_.size());
        sii_ = new StaticIndividualInfo(&cr_);
        const auto ffName = fileManager().getInputFilePath("forcefields/ACS-g1.xml");
        sii_->fillForceField(&msghandler_, ffName.c_str());
        // All molecules are in the training set
        const auto selName = fileManager().getTemporaryFilePath("selection.dat");
        FILE *fp      = fopen(selName.c_str(), "w");
        for(const char *mol : { "1-butanol", "acetone", "water" })
        {
            fprintf(fp, "%s|Train\n", mol);
        }
        fclose(fp);
        MolSelect gms;
        gms.read(selName.c_str());
        // Fit the electronegativities to the electrostatic potential
        sii_->target(iMolSelect::Train, eRMS::ESP)->setWeight(1);
        sii_->propagateWeightFittingTargets();
        molgen_ = new MolGen(&cr_);
        molgen_->addFitOption("chi");
        auto mpName = fileManager().getInputFilePath("molprop/molprop.xml");
        const std::vector<t_filenm> fnms = { { efXML, "-mp", nullptr, ffREAD, { mpName } } };
        CompoundReader        compR;
        EXPECT_EQ(3u, molgen_->Read(&msghandler_, fnms, sii_->forcefield(), gms,
                                    sii_->fittingTargetsConst(iMolSelect::Train),
                                    &compR, &forceComp_));
        sii_->generateOptimizationIndex(&msghandler_, molgen_, &cr_);
        sii_->fillVectors(molgen_->mindata(), 0.1);
        sii_->computeWeightedTemperature(true);
        fitComp_.init(&msghandler_, sii_, molgen_, false, &forceComp_,
                      sii_->forcefield()->chargeGenerationAlgorithm(), LossFunction::MSE);
        // Start from the parameters in the force field
        ACMInitializer initializer(sii_, false, 0);
        auto           ind = initializer.initialize();
        genome_ = ind->genome();
        delete ind;
        EXPECT_TRUE(msghandler_.ok());
    }

    //! Clean up
    ~SensitivityAnalysisTrainTest() override
    {
        delete molgen_;
        delete sii_;
    }

    /*! \brief Run the analysis and compare it to the parabolas through the chi-squared
     * The chi-squared is computed for each parameter and two values
     * half a percent of its range apart, below and above it, or on one
     * side of it if the other one is beyond a bound.
     * \param[in]  param     The parameters to start from
     * \param[out] converged Whether the analysis found all parameters
     *                       converged, may be nullptr
     * \param[out] pmin      The minima of the parabolas, limited to the
     *                       bounds, may be nullptr
     */
    void checkAnalysis(const std::vector<double> &param,
                       bool                      *converged = nullptr,
                       std::vector<double>       *pmin      = nullptr)
    {
        size_t np = sii_->nParam();
        ASSERT_EQ(np, param.size());
        *genome_.basesPtr() = param;
        double chi2_0 = chi2(param);
        EXPECT_GT(chi2_0, 0);

        SensitivityAnalysis sens;
        bool                minimum = sens.run(&msghandler_, sii_, &fitComp_, &genome_,
                                               iMolSelect::Train, SensitivityMode::OneD, true);
        EXPECT_TRUE(msghandler_.ok());
        EXPECT_EQ(param, genome_.bases());
        ASSERT_EQ(np, sens.forceConstants().size());
        // The summary
        const auto &res = sens.result();
        EXPECT_EQ(SensitivityMode::OneD, res.mode);
        EXPECT_DOUBLE_EQ(chi2_0, res.chi2_0);
        EXPECT_EQ(np, res.nParameters);
        EXPECT_EQ(2*np, res.nEvaluations);
        EXPECT_EQ(minimum, res.minimum);
        // The output has the summary and an entry per parameter
        ASSERT_TRUE(sens.json().contains("sensitivity"));
        const auto &sensJson = sens.json()["sensitivity"];
        EXPECT_EQ("1D", sensJson["mode"].get<std::string>());
        EXPECT_EQ("Train", sensJson["dataset"].get<std::string>());
        EXPECT_DOUBLE_EQ(chi2_0, sensJson["chi2_0"].get<double>());
        EXPECT_EQ(np, sensJson["n_parameters"].get<size_t>());
        EXPECT_EQ(2*np, sensJson["n_evaluations"].get<size_t>());
        EXPECT_EQ(0.05, sensJson["tolerance_fraction"].get<double>());
        EXPECT_EQ(minimum, sensJson["minimum"].get<bool>());
        EXPECT_FALSE(sensJson.contains("pairs"));
        const auto &entries = sensJson["parameters"];
        ASSERT_EQ(np, entries.size());
        size_t nconv = 0;

        bool allConverged = true;
        if (pmin)
        {
            pmin->clear();
        }
        for(size_t i = 0; i < np; i++)
        {
            double lower = sii_->lowerBound()[i];
            double upper = sii_->upperBound()[i];
            double delta = (upper - lower)/200;
            // The three parameter values in increasing order
            std::vector<double> x = { param[i] - delta, param[i], param[i] + delta };
            if (x[0] < lower)
            {
                x = { param[i], param[i] + delta, param[i] + 2*delta };
            }
            else if (x[2] > upper)
            {
                x = { param[i] - 2*delta, param[i] - delta, param[i] };
            }
            std::vector<double> y;
            for(double xx : x)
            {
                EXPECT_GE(xx, lower);
                EXPECT_LE(xx, upper);
                auto p = param;
                p[i]   = xx;
                y.push_back(chi2(p));
            }
            // The parabola chi2 = a x^2 + b x + c through the points
            double s01  = (y[1] - y[0])/(x[1] - x[0]);
            double s12  = (y[2] - y[1])/(x[2] - x[1]);
            double a    = (s12 - s01)/(x[2] - x[0]);
            double b    = s01 - a*(x[0] + x[1]);
            double pmin_i = -b/(2*a);
            // The force constant is a divided by the parameter squared
            double fc   = sens.forceConstants()[i];
            EXPECT_NEAR(a, fc*param[i]*param[i], 1e-5*std::abs(a)) << sii_->paramNames()[i];
            // Converged if the parabola curves upwards and the parameter
            // is within 5% of the range of the minimum, limited to the bounds
            bool conv = (a > 0 &&
                         std::abs(param[i] - std::clamp(pmin_i, lower, upper)) <= 0.05*(upper - lower));
            allConverged = allConverged && conv;
            if (conv)
            {
                nconv += 1;
            }
            if (pmin)
            {
                pmin->push_back(std::clamp(pmin_i, lower, upper));
            }

            const auto &entry = entries[i];
            const auto &name  = sii_->paramNames()[i];
            EXPECT_EQ(i, entry["index"].get<size_t>());
            EXPECT_EQ(name, entry["name"].get<std::string>());
            EXPECT_EQ(name, entry["particle"].get<std::string>() + " " + entry["parameter"].get<std::string>());
            EXPECT_EQ(lower, entry["lower"].get<double>());
            EXPECT_EQ(upper, entry["upper"].get<double>());
            EXPECT_EQ((upper - lower)/200, entry["delta"].get<double>());
            EXPECT_EQ(param[i], entry["p_orig"].get<double>()) << name;
            EXPECT_DOUBLE_EQ(chi2_0, entry["chi2_orig"].get<double>()) << name;
            // The three points, in increasing order of the parameter
            ASSERT_EQ(3u, entry["points"].size()) << name;
            for(size_t k = 0; k < 3; k++)
            {
                EXPECT_DOUBLE_EQ(x[k], entry["points"][k]["p"].get<double>()) << name;
                EXPECT_NEAR(y[k], entry["points"][k]["chi2"].get<double>(), 1e-12*y[k]) << name;
            }
            EXPECT_TRUE(entry["fit_ok"].get<bool>()) << name;
            EXPECT_NEAR(a, entry["a"].get<double>(), 1e-5*std::abs(a)) << name;
            EXPECT_EQ(fc, entry["force_constant"].get<double>()) << name;
            EXPECT_NEAR(pmin_i, entry["p_min"].get<double>(), 1e-5*std::abs(pmin_i) + 1e-6) << name;
            EXPECT_NEAR(std::clamp(pmin_i, lower, upper), entry["p_min_bounded"].get<double>(),
                        1e-5*std::abs(pmin_i) + 1e-6) << name;
            EXPECT_EQ(conv, entry["converged"].get<bool>()) << name;
            // The electronegativities all change the chi2
            EXPECT_FALSE(entry["flat"].get<bool>()) << name;
        }
        // The genome is a minimum if all parameters are converged
        EXPECT_EQ(allConverged, minimum);
        EXPECT_EQ(nconv, res.nParametersConverged);
        EXPECT_EQ(nconv, sensJson["n_parameters_converged"].get<size_t>());
        EXPECT_EQ(0u, res.nParametersFlat);
        EXPECT_EQ(0u, sensJson["n_parameters_flat"].get<size_t>());
        if (converged)
        {
            *converged = minimum;
        }
    }

    /*! \brief The three values of a parameter in the analysis
     * Half a percent of the range apart, on one side of the parameter
     * at a bound.
     * \param[in] param The parameters
     * \param[in] i     The parameter
     * \return the values in increasing order
     */
    std::vector<double> gridOf(const std::vector<double> &param, size_t i)
    {
        double lower = sii_->lowerBound()[i];
        double upper = sii_->upperBound()[i];
        double delta = (upper - lower)/200;
        if (param[i] - delta < lower)
        {
            return { param[i], param[i] + delta, param[i] + 2*delta };
        }
        else if (param[i] + delta > upper)
        {
            return { param[i] - 2*delta, param[i] - delta, param[i] };
        }
        return { param[i] - delta, param[i], param[i] + delta };
    }

    /*! \brief The parabola through three points
     * \param[in] x The abscissae
     * \param[in] y The ordinates
     * \param[in] t Where to evaluate it
     * \return the value of the parabola at t
     */
    static double lagrange(const std::vector<double> &x, const std::vector<double> &y, double t)
    {
        return (y[0]*(t - x[1])*(t - x[2])/((x[0] - x[1])*(x[0] - x[2])) +
                y[1]*(t - x[0])*(t - x[2])/((x[1] - x[0])*(x[1] - x[2])) +
                y[2]*(t - x[0])*(t - x[1])/((x[2] - x[0])*(x[2] - x[1])));
    }

    //! \return the curvature a of the parabola a x^2 + b x + c through three points
    static double curvature(const std::vector<double> &x, const std::vector<double> &y)
    {
        double s01 = (y[1] - y[0])/(x[1] - x[0]);
        double s12 = (y[2] - y[1])/(x[2] - x[1]);
        return (s12 - s01)/(x[2] - x[0]);
    }

    /*! \brief Run the two-dimensional analysis and compare it to the chi-squared
     * The chi-squared of each corner of each pair is computed here, the
     * coupling is fitted to what the parabolas through the points of each
     * parameter do not describe.
     * \param[in]  param     The parameters to start from
     * \param[out] converged Whether the analysis found the parameters and
     *                       all pairs converged, may be nullptr
     * \param[out] newton    The parameters after a projected Newton step,
     *                       may be nullptr. Parameters at a bound with the
     *                       gradient pointing out of the bounds stay there,
     *                       the step of the others uses the pseudo-inverse
     *                       of their Hessian, and is limited to the bounds.
     * \param[out] output    The output of the analysis, may be nullptr
     */
    void checkPairs(const std::vector<double> &param,
                    bool                      *converged = nullptr,
                    std::vector<double>       *newton    = nullptr,
                    nlohmann::ordered_json    *output    = nullptr)
    {
        size_t np     = sii_->nParam();
        size_t npairs = np*(np - 1)/2;
        ASSERT_EQ(np, param.size());
        *genome_.basesPtr() = param;
        double chi2_0 = chi2(param);

        SensitivityAnalysis sens;
        bool                minimum = sens.run(&msghandler_, sii_, &fitComp_, &genome_,
                                               iMolSelect::Train, SensitivityMode::TwoD, true);
        EXPECT_TRUE(msghandler_.ok());
        EXPECT_EQ(param, genome_.bases());
        const auto &res = sens.result();
        EXPECT_EQ(SensitivityMode::TwoD, res.mode);
        EXPECT_EQ(npairs, res.nPairs);
        EXPECT_EQ(2*np*np, res.nEvaluations);
        EXPECT_EQ(minimum, res.minimum);
        ASSERT_TRUE(res.hessian.has_value());

        const auto &sj = sens.json()["sensitivity"];
        EXPECT_EQ("2D", sj["mode"].get<std::string>());
        EXPECT_EQ(npairs, sj["n_pairs"].get<size_t>());
        EXPECT_EQ(2*np*np, sj["n_evaluations"].get<size_t>());
        EXPECT_EQ(minimum, sj["minimum"].get<bool>());
        ASSERT_EQ(np, sj["parameters"].size());
        ASSERT_EQ(npairs, sj["pairs"].size());
        ASSERT_TRUE(sj.contains("hessian_diagnostic"));
        EXPECT_EQ(res.hessian->positiveDefinite, sj["hessian_diagnostic"]["positive_definite"].get<bool>());

        // The points of each parameter
        std::vector<std::vector<double>> xs(np), ys(np);
        bool allConverged = true;
        for (size_t i = 0; i < np; i++)
        {
            xs[i] = gridOf(param, i);
            for (double x : xs[i])
            {
                auto p = param;
                p[i]   = x;
                ys[i].push_back(chi2(p));
            }
            allConverged = allConverged && sj["parameters"][i]["converged"].get<bool>();
        }
        // The Hessian and the gradient from the analysis, for the Newton step
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(np, np);
        Eigen::VectorXd g = Eigen::VectorXd::Zero(np);
        for (size_t i = 0; i < np; i++)
        {
            H(i, i) = 2*sj["parameters"][i]["a"].get<double>();
        }
        size_t k = 0, nconv = 0, nminimum = 0;
        for (size_t i = 0; i < np; i++)
        {
            for (size_t j = i + 1; j < np; j++)
            {
                const auto &pair  = sj["pairs"][k++];
                const std::string label = sii_->paramNames()[i] + " / " + sii_->paramNames()[j];
                EXPECT_EQ(i, pair["i"].get<size_t>());
                EXPECT_EQ(j, pair["j"].get<size_t>());
                EXPECT_EQ(sii_->paramNames()[i], pair["names"][0].get<std::string>());
                EXPECT_EQ(sii_->paramNames()[j], pair["names"][1].get<std::string>());
                double num = 0, den = 0;
                for (size_t a = 0; a < 3; a++)
                {
                    double x = xs[i][a];
                    EXPECT_GE(x, sii_->lowerBound()[i]);
                    EXPECT_LE(x, sii_->upperBound()[i]);
                    EXPECT_DOUBLE_EQ(x, pair["grid"]["p_i"][a].get<double>()) << label;
                    for (size_t b = 0; b < 3; b++)
                    {
                        double y = xs[j][b];
                        EXPECT_DOUBLE_EQ(y, pair["grid"]["p_j"][b].get<double>()) << label;
                        double f;
                        if (x == param[i])
                        {
                            f = ys[j][b];
                        }
                        else if (y == param[j])
                        {
                            f = ys[i][a];
                        }
                        else
                        {
                            auto p = param;
                            p[i]   = x;
                            p[j]   = y;
                            f      = chi2(p);
                            double xy = (x - param[i])*(y - param[j]);
                            double r  = f - (lagrange(xs[i], ys[i], x) + lagrange(xs[j], ys[j], y) - chi2_0);
                            num += r*xy;
                            den += xy*xy;
                        }
                        EXPECT_NEAR(f, pair["grid"]["chi2"][a][b].get<double>(), 1e-12*std::abs(f))
                            << label << " " << a << " " << b;
                    }
                }
                double aij   = num/den;
                double di    = (sii_->upperBound()[i] - sii_->lowerBound()[i])/200;
                double dj    = (sii_->upperBound()[j] - sii_->lowerBound()[j])/200;
                EXPECT_NEAR(aij, pair["a_ij"].get<double>(), 1e-6*std::abs(aij) + 1e-9*chi2_0/(di*dj)) << label;
                // The shape from the curvatures computed here, unless the
                // determinant is too close to zero to tell
                double hi    = 2*curvature(xs[i], ys[i]);
                double hj    = 2*curvature(xs[j], ys[j]);
                double det   = hi*hj - aij*aij;
                const auto shape = pair["shape"].get<std::string>();
                if (std::abs(det) > 1e-4*(std::abs(hi*hj) + aij*aij))
                {
                    const std::string expected = det < 0 ? "saddle" : (hi > 0 ? "minimum" : "maximum");
                    EXPECT_EQ(expected, shape) << label;
                }
                bool conv = pair["converged"].get<bool>();
                if (shape != "minimum")
                {
                    EXPECT_FALSE(conv) << label;
                    EXPECT_TRUE(pair["p_min"].is_null()) << label;
                }
                else
                {
                    nminimum += 1;
                }
                allConverged = allConverged && conv;
                nconv       += conv ? 1 : 0;
                H(i, j) = H(j, i) = pair["a_ij"].get<double>();
                g(i)    = pair["gradient"][0].get<double>();
                g(j)    = pair["gradient"][1].get<double>();
            }
        }
        // The genome is a minimum if all parameters and pairs are converged
        EXPECT_EQ(allConverged, minimum);
        EXPECT_EQ(nconv, res.nPairsConverged);
        EXPECT_EQ(nconv, sj["n_pairs_converged"].get<size_t>());
        EXPECT_EQ(nminimum, sj["pair_shapes"]["minimum"].get<size_t>());
        size_t nshapes = 0;
        for (const auto &s : sj["pair_shapes"].items())
        {
            nshapes += s.value().get<size_t>();
        }
        EXPECT_EQ(npairs, nshapes);
        if (converged)
        {
            *converged = minimum;
        }
        if (output)
        {
            *output = sens.json();
        }
        if (newton)
        {
            // The parameters that are free to move
            std::vector<Eigen::Index> free;
            for (size_t i = 0; i < np; i++)
            {
                bool outLower = param[i] <= sii_->lowerBound()[i] && g(i) > 0;
                bool outUpper = param[i] >= sii_->upperBound()[i] && g(i) < 0;
                if (!outLower && !outUpper)
                {
                    free.push_back(i);
                }
            }
            auto            nf = static_cast<Eigen::Index>(free.size());
            Eigen::MatrixXd Hf(nf, nf);
            Eigen::VectorXd gf(nf);
            for (Eigen::Index m = 0; m < nf; m++)
            {
                gf(m) = g(free[m]);
                for (Eigen::Index n = 0; n < nf; n++)
                {
                    Hf(m, n) = H(free[m], free[n]);
                }
            }
            // Pseudo-inverse, since a shift of all electronegativities
            // may not change the charges
            Eigen::VectorXd step = Eigen::VectorXd::Zero(nf);
            if (nf > 0)
            {
                const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Hf);
                const auto &eval = solver.eigenvalues();
                const auto &evec = solver.eigenvectors();
                double      emax = eval.cwiseAbs().maxCoeff();
                for (Eigen::Index m = 0; m < eval.size(); m++)
                {
                    if (std::abs(eval(m)) > 1e-6*emax)
                    {
                        step -= evec.col(m)*(evec.col(m).dot(gf)/eval(m));
                    }
                }
            }
            *newton = param;
            for (Eigen::Index m = 0; m < nf; m++)
            {
                size_t i     = free[m];
                (*newton)[i] = std::clamp(param[i] + step(m), sii_->lowerBound()[i], sii_->upperBound()[i]);
            }
        }
    }

    /*! \brief Describe the state of an analysis, for when a test fails
     * \param[in] js The output of the analysis
     * \return the chi-squared, and for each parameter its value relative
     *         to its bounds, its gradient and whether it is converged, and
     *         the pairs that are not converged
     */
    std::string describe(const nlohmann::ordered_json &js)
    {
        const auto &sj  = js["sensitivity"];
        std::string str = "chi2_0 " + std::to_string(sj["chi2_0"].get<double>()) +
            " pairs converged " + std::to_string(sj["n_pairs_converged"].get<size_t>()) +
            " Hessian " + sj["hessian_diagnostic"].dump() + "\n";
        for (const auto &p : sj["parameters"])
        {
            double lower = p["lower"].get<double>();
            double upper = p["upper"].get<double>();
            str += "  " + p["name"].get<std::string>() +
                " p " + std::to_string(p["p_orig"].get<double>()) +
                " in [" + std::to_string(lower) + ", " + std::to_string(upper) + "]" +
                " p_min_bounded " + p["p_min_bounded"].dump() +
                " converged " + p["converged"].dump() + "\n";
        }
        for (const auto &pair : sj["pairs"])
        {
            if (!pair["converged"].get<bool>())
            {
                str += "  not converged: " + pair["names"].dump() + " " + pair["shape"].get<std::string>() +
                    " rho " + pair["rho"].dump() + " gradient " + pair["gradient"].dump() +
                    " p_min_bounded " + pair["p_min_bounded"].dump() + "\n";
            }
        }
        return str;
    }

    /*! \brief Compute the chi-squared for some parameters
     * \param[in] param The parameters
     * \return the chi-squared of the training set
     */
    double chi2(const std::vector<double> &param)
    {
        sii_->updateForceField(&msghandler_, {}, param);
        auto cdc = CalcDev::Compute;
        fitComp_.distributeTasks(cdc);
        return fitComp_.calcDeviation(&msghandler_, cdc, iMolSelect::Train);
    }
};

// Test: The genome keeps the parameters exactly, also in single precision,
// such that parameters at their bounds stay within them
TEST_F(SensitivityAnalysisTrainTest, GenomeKeepsParameters)
{
    EXPECT_EQ(sii_->defaultParam(), genome_.bases());
    ACMIndividual ind(0, sii_);
    ind.addParam(0.1);
    ASSERT_EQ(1u, ind.genome().bases().size());
    EXPECT_EQ(0.1, ind.genome().bases()[0]);
}

// Test: The force constants and the output are those of the parabolas

// through the chi-squared of each parameter and the two values around it.
TEST_F(SensitivityAnalysisTrainTest, ParabolaThroughChi2)
{
    // The atom types c2, c3, h, o2 and o3
    ASSERT_EQ(5u, sii_->nParam());
    checkAnalysis(genome_.bases());
}

// Test: Parameters at or close to a bound are analysed with two values on
// the other side, and the output has the parameters themselves.
TEST_F(SensitivityAnalysisTrainTest, ParametersAtBounds)
{
    ASSERT_EQ(5u, sii_->nParam());
    auto param  = genome_.bases();
    const auto delta  = [this](size_t i) { return (sii_->upperBound()[i] - sii_->lowerBound()[i])/200; };
    param[0]    = sii_->lowerBound()[0];
    param[1]    = sii_->upperBound()[1];
    param[2]    = sii_->lowerBound()[2] + 0.5*delta(2);
    param[3]    = sii_->upperBound()[3] - 0.5*delta(3);
    checkAnalysis(param);
}

// Test: Moving the parameters to the minima of their parabolas, the
// analysis finds them converged after a few steps
TEST_F(SensitivityAnalysisTrainTest, ConvergedAtMinimum)
{
    auto param     = genome_.bases();
    bool converged = false;
    std::vector<double> pmin;
    checkAnalysis(param, &converged, &pmin);
    // At the parameters of the force field the parameters are not converged
    EXPECT_FALSE(converged);
    for(int step = 0; step < 20 && !converged; step++)
    {
        param = pmin;
        checkAnalysis(param, &converged, &pmin);
    }
    EXPECT_TRUE(converged);
}

// Test: After the analysis the genome and the force field are as before
TEST_F(SensitivityAnalysisTrainTest, RestoresParameters)
{
    const auto param  = genome_.bases();
    double chi2_0 = chi2(param);
    SensitivityAnalysis sens;
    (void) sens.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train, SensitivityMode::OneD, true);
    EXPECT_TRUE(msghandler_.ok());
    EXPECT_EQ(param, genome_.bases());
    // Without updating the force field, the chi-squared is the same
    auto cdc = CalcDev::Compute;
    fitComp_.distributeTasks(cdc);
    EXPECT_EQ(chi2_0, fitComp_.calcDeviation(&msghandler_, cdc, iMolSelect::Train));
}

// Test: With mode No nothing is computed and nothing is in the output
TEST_F(SensitivityAnalysisTrainTest, ModeNoDoesNothing)
{
    SensitivityAnalysis sens;
    EXPECT_TRUE(sens.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train,
                         SensitivityMode::No, true));
    EXPECT_TRUE(sens.json().is_null());
    EXPECT_TRUE(sens.forceConstants().empty());
    EXPECT_EQ(SensitivityMode::No, sens.result().mode);
    EXPECT_EQ(0u, sens.result().nEvaluations);
}

// Test: A second run replaces the output of the first
TEST_F(SensitivityAnalysisTrainTest, SecondRunReplacesOutput)
{
    SensitivityAnalysis sens;
    (void) sens.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train,
                    SensitivityMode::OneD, true);
    const auto first = sens.json();
    (void) sens.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train,
                    SensitivityMode::OneD, true);
    EXPECT_EQ(first, sens.json());
    EXPECT_EQ(sii_->nParam(), sens.forceConstants().size());
    (void) sens.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train,
                    SensitivityMode::No, true);
    EXPECT_TRUE(sens.json().is_null());
}

// Test: The file has the same content as the output in memory
TEST_F(SensitivityAnalysisTrainTest, WritesJsonFile)
{
    SensitivityAnalysis sens;
    (void) sens.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train,
                    SensitivityMode::OneD, true);
    const auto fileName = fileManager().getTemporaryFilePath("sensitivity.json");
    sens.write(fileName);
    std::ifstream in(fileName);
    ASSERT_TRUE(in.good());
    const auto fromFile = nlohmann::ordered_json::parse(in);
    EXPECT_EQ(sens.json(), fromFile);
}

// Test: The coupling of each pair of parameters is fitted to the chi-squared
// of the corners, the output has the grid, the shape and the convergence
TEST_F(SensitivityAnalysisTrainTest, TwoDCrossTerms)
{
    ASSERT_EQ(5u, sii_->nParam());
    checkPairs(genome_.bases());
}

// Test: The same for parameters at or close to bounds, where the grids are
// on one side of the parameters
TEST_F(SensitivityAnalysisTrainTest, TwoDParametersAtBounds)
{
    ASSERT_EQ(5u, sii_->nParam());
    auto param  = genome_.bases();
    const auto delta  = [this](size_t i) { return (sii_->upperBound()[i] - sii_->lowerBound()[i])/200; };
    param[0]    = sii_->lowerBound()[0];
    param[1]    = sii_->upperBound()[1];
    param[2]    = sii_->lowerBound()[2] + 0.5*delta(2);
    param[3]    = sii_->upperBound()[3] - 0.5*delta(3);
    checkPairs(param);
}

// Test: The one-dimensional part of the two-dimensional analysis is the
// one-dimensional analysis
TEST_F(SensitivityAnalysisTrainTest, TwoDContainsOneD)
{
    SensitivityAnalysis sens1, sens2;
    (void) sens1.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train,
                     SensitivityMode::OneD, true);
    (void) sens2.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train,
                     SensitivityMode::TwoD, true);
    EXPECT_EQ(sens1.json()["sensitivity"]["parameters"], sens2.json()["sensitivity"]["parameters"]);
    EXPECT_EQ(sens1.forceConstants(), sens2.forceConstants());
    EXPECT_EQ(sens1.result().nParametersConverged, sens2.result().nParametersConverged);
    // In 2D the pairs have to be converged as well
    if (!sens1.result().minimum)
    {
        EXPECT_FALSE(sens2.result().minimum);
    }
}

// Test: Stepping to the minimum of the surface given by the analysis,
// the analysis finds the parameters and all pairs converged
TEST_F(SensitivityAnalysisTrainTest, TwoDConvergedAfterNewtonSteps)
{
    auto                   param     = genome_.bases();
    bool                   converged = false;
    std::vector<double>    next;
    nlohmann::ordered_json js;
    checkPairs(param, &converged, &next, &js);
    // At the parameters of the force field they are not converged
    EXPECT_FALSE(converged);
    std::string history = "Step 0: " + describe(js);
    for (int step = 1; step <= 20 && !converged; step++)
    {
        param = next;
        checkPairs(param, &converged, &next, &js);
        history += "Step " + std::to_string(step) + ": " + describe(js);
    }
    EXPECT_TRUE(converged) << history;
}

// Test: After the two-dimensional analysis the genome and the force field
// are as before
TEST_F(SensitivityAnalysisTrainTest, TwoDRestoresParameters)
{
    const auto param  = genome_.bases();
    double chi2_0 = chi2(param);
    SensitivityAnalysis sens;
    (void) sens.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train,
                    SensitivityMode::TwoD, true);
    EXPECT_TRUE(msghandler_.ok());
    EXPECT_EQ(param, genome_.bases());
    auto cdc = CalcDev::Compute;
    fitComp_.distributeTasks(cdc);
    EXPECT_EQ(chi2_0, fitComp_.calcDeviation(&msghandler_, cdc, iMolSelect::Train));
}

} // namespace

} // namespace alexandria
