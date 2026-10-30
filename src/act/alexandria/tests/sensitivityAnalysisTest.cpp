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

#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

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
#include "act/utility/jsontree.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/filetypes.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/*! \brief Find an object in a JsonTree
 * \param[in] tree The tree
 * \param[in] key  The key of the object
 * \return the object or nullptr if not found
 */
const JsonTree *findObject(const JsonTree &tree, const std::string &key)
{
    for(const auto &obj : tree.objects())
    {
        if (obj.key() == key)
        {
            return &obj;
        }
    }
    return nullptr;
}

/*! \brief The value of a number in a JsonTree
 * \param[in] tree The tree
 * \param[in] key  The key of the number
 * \return the number
 */
double findNumber(const JsonTree &tree, const std::string &key)
{
    auto obj = findObject(tree, key);
    EXPECT_NE(nullptr, obj) << key;
    return obj ? std::stod(obj->value()) : 0;
}

// Test: With an empty genome there is nothing to analyse, the routine
// reports a minimum and adds nothing to the output.
TEST(SensitivityAnalysisTest, EmptyGenome)
{
    MsgHandler           msghandler;
    CommunicationRecord  cr;
    cr.init(cr.size());
    StaticIndividualInfo sii(&cr);
    ga::Genome           genome;
    JsonTree             jtree("SensitivityAnalysisTest");
    SensitivityAnalysis  sens;
    EXPECT_TRUE(sens.run(&msghandler, &sii, nullptr, &genome, iMolSelect::Train, &jtree, true));
    EXPECT_TRUE(jtree.objects().empty());
    EXPECT_TRUE(sens.forceConstants().empty());
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
        auto ffName = fileManager().getInputFilePath("ACS-g.xml");
        sii_->fillForceField(&msghandler_, ffName.c_str());
        // All molecules are in the training set
        auto  selName = fileManager().getTemporaryFilePath("selection.dat");
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
        auto mpName = fileManager().getInputFilePath("../../molprop/tests/molprop.xml");
        std::vector<t_filenm> fnms = { { efXML, "-mp", nullptr, ffREAD, { mpName } } };
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
    ~SensitivityAnalysisTrainTest()
    {
        delete molgen_;
        delete sii_;
    }

    /*! \brief Run the analysis and compare it to the parabolas through the chi-squared
     * The chi-squared is computed for each parameter and two values
     * half a percent of its range apart, below and above it, or on one
     * side of it if the other one is beyond a bound.
     * \param[in] param The parameters to start from
     */
    void checkAnalysis(const std::vector<double> &param)
    {
        size_t np = sii_->nParam();
        ASSERT_EQ(np, param.size());
        *genome_.basesPtr() = param;
        double chi2_0 = chi2(param);
        EXPECT_GT(chi2_0, 0);

        JsonTree            jtree("SensitivityAnalysisTest");
        SensitivityAnalysis sens;
        bool                minimum = sens.run(&msghandler_, sii_, &fitComp_, &genome_,
                                               iMolSelect::Train, &jtree, true);
        EXPECT_TRUE(msghandler_.ok());
        EXPECT_EQ(param, genome_.bases());
        ASSERT_EQ(np, sens.forceConstants().size());
        ASSERT_EQ(1u, jtree.objects().size());
        const auto &sensTree = jtree.objects()[0];
        EXPECT_EQ("sensitivity", sensTree.key());
        ASSERT_EQ(np, sensTree.objects().size());

        bool allConvex = true;
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
            double pmin = -b/(2*a);
            // The force constant is a divided by the parameter squared
            double fc   = sens.forceConstants()[i];
            EXPECT_NEAR(a, fc*param[i]*param[i], 1e-5*std::abs(a)) << sii_->paramNames()[i];
            allConvex   = allConvex && a >= 0;

            const auto &entry = sensTree.objects()[i];
            EXPECT_EQ(std::to_string(i), entry.key());
            EXPECT_EQ(sii_->paramNames()[i],
                      findObject(entry, "particle")->value() + " " + findObject(entry, "parameter")->value());
            EXPECT_NEAR(param[i], findNumber(entry, "p_orig"), 1e-6) << sii_->paramNames()[i];
            EXPECT_NEAR(chi2_0, findNumber(entry, "chi2_orig"), 1e-6*chi2_0 + 1e-6);
            EXPECT_NEAR(pmin, findNumber(entry, "p_min"), 1e-5*std::abs(pmin) + 1e-6)
                << sii_->paramNames()[i];
        }
        // The genome is a minimum if all parabolas curve upwards
        EXPECT_EQ(allConvex, minimum);
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
    auto delta  = [this](size_t i) { return (sii_->upperBound()[i] - sii_->lowerBound()[i])/200; };
    param[0]    = sii_->lowerBound()[0];
    param[1]    = sii_->upperBound()[1];
    param[2]    = sii_->lowerBound()[2] + 0.5*delta(2);
    param[3]    = sii_->upperBound()[3] - 0.5*delta(3);
    checkAnalysis(param);
}

// Test: After the analysis the genome and the force field are as before
TEST_F(SensitivityAnalysisTrainTest, RestoresParameters)
{
    auto   param  = genome_.bases();
    double chi2_0 = chi2(param);
    SensitivityAnalysis sens;
    (void) sens.run(&msghandler_, sii_, &fitComp_, &genome_, iMolSelect::Train, nullptr, true);
    EXPECT_TRUE(msghandler_.ok());
    EXPECT_EQ(param, genome_.bases());
    // Without updating the force field, the chi-squared is the same
    auto cdc = CalcDev::Compute;
    fitComp_.distributeTasks(cdc);
    EXPECT_EQ(chi2_0, fitComp_.calcDeviation(&msghandler_, cdc, iMolSelect::Train));
}

} // namespace

} // namespace alexandria
