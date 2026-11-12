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
/*! \internal \file
 * \brief
 * Tests that updating the force field one parameter at a time, as the
 * MCMC steps of training do, gives the same results as filling all
 * parameters from scratch.
 *
 * Training only updates what depends on the parameters that changed: the
 * topology entries of their interaction types, the parameters of the
 * non-bonded pairs, the charges, and on the GPU the rows of the parameter
 * tables. Each of these keeps track of versions to know what to do again.
 * A mistake there gives energies of a mix of old and new parameters, which
 * a test of a single evaluation does not find. Here, for each kind of
 * parameter of a force field, the parameter is changed in one set of
 * molecules that was evaluated before, as in an MCMC step, and the result
 * is compared with that of another set of molecules that is evaluated for
 * the first time, with the same parameters. Then the change is undone, as
 * for a rejected step, which has to give the result from before again.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_alexandria_tests
 */
#include "actpre.h"

#include <cmath>
#include <cstdlib>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/molgen.h"
#include "act/alexandria/optimizationindex.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forces/forcecomputer.h"
#include "act/ga/ga_test_helper.h"
#include "act/gpu/gpu.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/gputests.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! A force field with molecules to train on
struct UpdateCase
{
    //! The input files
    GaTestFiles       files;
    //! The fitting targets
    std::vector<eRMS> targets;
};

//! The result of an evaluation
struct Evaluation
{
    //! The value of the objective function
    double                 chi2 = 0;
    //! The contribution of each target with a weight
    std::map<eRMS, double> targets;
    //! The charges of the particles of the molecules for training
    std::vector<double>    charges;
};

/*! \brief The names of all parameters of a force field that can be trained
 * \param[in] ffFile The force field, relative to the test data
 * \return The names
 */
std::vector<std::string> trainableNames(const std::string &ffFile)
{
    ForceField pd;
    readForceField(gmx::test::TestFileManager::getInputFilePath(ffFile), &pd);
    std::set<std::string> names;
    for (const auto &fs : pd.forcesConst())
    {
        for (const auto &fpl : fs.second.parametersConst())
        {
            for (const auto &param : fpl.second)
            {
                if (param.second.isMutable())
                {
                    names.insert(param.first);
                }
            }
        }
    }
    for (const auto &pt : pd.particleTypesConst())
    {
        for (const auto &p : pt.second.parametersConst())
        {
            if (p.second.isMutable())
            {
                names.insert(p.first);
            }
        }
    }
    return { names.begin(), names.end() };
}

/*! \brief The parameter of the force field of an optimization index
 * \param[in] sii The information shared by the individuals
 * \param[in] k   The index
 * \return The parameter
 */
ForceFieldParameter *parameter(StaticIndividualInfo *sii, size_t k)
{
    auto &oi = (*sii->optIndexPtr())[k];
    if (!oi.forceFieldParameter())
    {
        oi.findForceFieldParameter(sii->forcefield());
    }
    return oi.forceFieldParameter();
}

//! \return the values of the optimized parameters of the force field of sii
std::vector<double> parameterValues(StaticIndividualInfo *sii)
{
    std::vector<double> values;
    for (size_t k = 0; k < sii->optIndex().size(); k++)
    {
        auto p = parameter(sii, k);
        values.push_back(p ? p->value() : 0.0);
    }
    return values;
}

/*! \brief Another value of a parameter, within its bounds
 * \param[in] p The parameter
 * \return The new value
 */
double otherValue(const ForceFieldParameter &p)
{
    double v     = p.value();
    double range = p.maximum() - p.minimum();
    double d     = range > 0 ? 0.05*range : 0.05*std::abs(v) + 1e-3;
    double nv    = v + d;
    if (nv > p.maximum())
    {
        nv = v - d;
    }
    if (nv < p.minimum())
    {
        nv = 0.5*(p.minimum() + p.maximum());
    }
    return nv;
}

/*! \brief Evaluate the molecules for training with these parameters
 * As MCMCMutator::mutate, the parameters are distributed, after which the
 * deviations are computed.
 * \param[in] helper  The training components
 * \param[in] values  The values of all optimized parameters
 * \param[in] changed The parameters that changed, all if empty
 * \param[in] charges Whether to keep the charges, which the processor only
 *                    has if it generates them
 * \return The result
 */
Evaluation evaluate(GaTestHelper              *helper,
                    const std::vector<double> &values,
                    const std::set<int>       &changed,
                    bool                       charges)
{
    auto *fc = helper->fitnessComputer;
    fc->distributeParameters(helper->msghandler, &values, changed);
    Evaluation ev;
    ev.chi2 = fc->calcDeviation(helper->msghandler, CalcDev::Compute, iMolSelect::Train);
    for (const auto &t : helper->sii->fittingTargetsConst(iMolSelect::Train))
    {
        if (t.second.weight() > 0)
        {
            ev.targets[t.first] = t.second.chiSquared();
        }
    }
    if (charges)
    {
        for (const auto &mol : helper->molgen->actmols())
        {
            if (mol.datasetType() == iMolSelect::Train)
            {
                for (const auto &a : mol.atomsConst())
                {
                    ev.charges.push_back(a.charge());
                }
            }
        }
    }
    return ev;
}

/*! \brief Compare two evaluations
 * \param[in] a         One
 * \param[in] b         The other
 * \param[in] tolerance Relative tolerance
 */
void compare(const Evaluation &a, const Evaluation &b, double tolerance)
{
    EXPECT_NEAR(a.chi2, b.chi2, tolerance*std::abs(b.chi2)) << "objective function";
    ASSERT_EQ(a.targets.size(), b.targets.size());
    for (const auto &t : b.targets)
    {
        auto at = a.targets.find(t.first);
        ASSERT_TRUE(at != a.targets.end());
        EXPECT_NEAR(at->second, t.second, tolerance*std::abs(t.second))
            << "target " << rmsName(t.first);
    }
    ASSERT_EQ(a.charges.size(), b.charges.size());
    for (size_t i = 0; i < b.charges.size(); i++)
    {
        EXPECT_NEAR(a.charges[i], b.charges[i], tolerance*(1 + std::abs(b.charges[i])))
            << "charge " << i;
    }
}

/*! \brief The kind of an optimized parameter
 * A kind is an interaction type with a parameter name, or a parameter
 * of the particle types.
 * \param[in] oi The optimization index of the parameter
 * \return The kind
 */
std::string kindOf(const OptimizationIndex &oi)
{
    return oi.particleType().empty() ?
        interactionTypeToString(oi.iType()) + ":" + oi.parameterType() :
        "particle:" + oi.parameterType();
}

/*! \brief Write a force field with other values of the optimized parameters
 * The force field is read from its file and written with the values, such
 * that reading it gives the parameters these values from the start,
 * independent of how training changes them.
 * \param[in] ffFile   The force field, relative to the test data or absolute
 * \param[in] optIndex The optimized parameters
 * \param[in] values   Their values
 * \param[in] fileName The file to write
 */
void writeValues(const std::string                    &ffFile,
                 const std::vector<OptimizationIndex> &optIndex,
                 const std::vector<double>            &values,
                 const std::string                    &fileName)
{
    ForceField pd;
    readForceField(ffFile[0] == '/' ? ffFile : gmx::test::TestFileManager::getInputFilePath(ffFile), &pd);
    for (size_t k = 0; k < optIndex.size(); k++)
    {
        OptimizationIndex mine = optIndex[k];
        mine.findForceFieldParameter(&pd);
        ASSERT_TRUE(mine.forceFieldParameter() != nullptr) << optIndex[k].name();
        mine.forceFieldParameter()->setValue(values[k]);
    }
    writeForceField(fileName, &pd, false);
}

/*! \brief Random values of the optimized parameters within their bounds, as train_ff -random_init
 * The force fields of the tests are not trained, and many parameters are
 * zero, for instance the van der Waals epsilon. With a combination rule
 * such as HogervorstEpsilon a pair then keeps a zero epsilon if one of its
 * particle types does, and changing one parameter changes nothing.
 * \param[in] sii The information shared by the individuals
 * \return The values
 */
std::vector<double> randomValues(StaticIndividualInfo *sii)
{
    std::mt19937 gen(1993);
    std::vector<double> values;
    for (size_t k = 0; k < sii->optIndex().size(); k++)
    {
        auto p = parameter(sii, k);
        std::uniform_real_distribution<double> dist(p->minimum(), p->maximum());
        values.push_back(p->maximum() > p->minimum() ? dist(gen) : p->value());
    }
    return values;
}

/*! \brief Change the parameters one at a time, compared with evaluating them from scratch
 * \param[in] fileManager For temporary files
 * \param[in] uc          The force field and molecules
 * \param[in] gpu         Whether to compute on the GPU
 * \param[in] tolerance   Relative tolerance
 */
void checkUpdates(gmx::test::TestFileManager *fileManager,
                  const UpdateCase           &uc,
                  bool                        gpu,
                  double                      tolerance)
{
    auto fit = trainableNames(uc.files.forceField);
    ASSERT_FALSE(fit.empty());
    // The shells start from the same positions each time, otherwise the
    // results depend on the evaluations before, see setKeepShells
    auto makeForceComputer = [gpu]()
    {
        auto fc = new ForceComputer();
        fc->setKeepShells(false);
        fc->setGpuResident(gpu);
        return fc;
    };
    // The force field with random values of the parameters, the start of both
    auto baseFiles = uc.files;
    {
        GaTestHelper original(1, fit, uc.targets, uc.files, makeForceComputer());
        ASSERT_FALSE(original.sii->optIndex().empty()) << "No parameters to train in " << uc.files.forceField;
        baseFiles.forceField = fileManager->getTemporaryFilePath("base.xml");
        writeValues(uc.files.forceField, original.sii->optIndex(), randomValues(original.sii),
                    baseFiles.forceField);
    }
    GaTestHelper incremental(1, fit, uc.targets, baseFiles, makeForceComputer());
    ASSERT_TRUE(incremental.msghandler->ok());
    auto *sii = incremental.sii;
    // The charges are only on the processor if it generates them
    const bool charges = !gpu;
    auto values = parameterValues(sii);
    auto before = evaluate(&incremental, values, {}, charges);
    EXPECT_GT(before.chi2, 0);
    // First each parameter is changed as in an MCMC step that is rejected,
    // which takes little time
    const auto         &optIndex = sii->optIndex();
    std::vector<double> probe(optIndex.size());
    std::vector<char>   probeChanged(optIndex.size(), 0);
    for (size_t k = 0; k < optIndex.size(); k++)
    {
        SCOPED_TRACE(gmx::formatString("parameter %zu %s", k, optIndex[k].name().c_str()));
        auto p = parameter(sii, k);
        ASSERT_TRUE(p != nullptr);
        probe[k]        = otherValue(*p);
        auto values1    = values;
        values1[k]      = probe[k];
        auto tried      = evaluate(&incremental, values1, { static_cast<int>(k) }, charges);
        probeChanged[k] = tried.chi2 != before.chi2;
        compare(evaluate(&incremental, values, { static_cast<int>(k) }, charges), before, tolerance);
    }
    // Then, for each kind of parameter, the first one that changed the
    // result, or the first one if none did, is compared with evaluating from
    // scratch, which takes longer. With ACT_TEST_ALL_PARAMETERS all are.
    // Since a mistake may also be that a change does not change the result,
    // whether a kind of parameter changes it is that of evaluating from scratch.
    static const bool             all = std::getenv("ACT_TEST_ALL_PARAMETERS") != nullptr;
    std::map<std::string, size_t> chosen;
    for (size_t k = 0; k < optIndex.size(); k++)
    {
        auto kind = kindOf(optIndex[k]);
        auto c    = chosen.find(kind);
        if (c == chosen.end())
        {
            chosen[kind] = k;
        }
        else if (!probeChanged[c->second] && probeChanged[k])
        {
            c->second = k;
        }
    }
    std::vector<size_t> toCompare;
    for (size_t k = 0; k < optIndex.size(); k++)
    {
        if (all || chosen[kindOf(optIndex[k])] == k)
        {
            toCompare.push_back(k);
        }
    }
    std::set<std::string> kinds, relevant;
    for (size_t k : toCompare)
    {
        const auto &oi = optIndex[k];
        SCOPED_TRACE(gmx::formatString("parameter %zu %s", k, oi.name().c_str()));
        kinds.insert(kindOf(oi));
        auto values1 = values;
        values1[k]   = probe[k];
        // All from scratch, with the changed force field read from a file
        auto files       = baseFiles;
        files.forceField = fileManager->getTemporaryFilePath(gmx::formatString("changed%zu.xml", k));
        writeValues(baseFiles.forceField, optIndex, values1, files.forceField);
        GaTestHelper fresh(1, fit, uc.targets, files, makeForceComputer());
        ASSERT_EQ(optIndex.size(), fresh.sii->optIndex().size());
        // The values as written in the file, of which only that of parameter k
        // differs, unless writing rounded others too
        auto changed = parameterValues(fresh.sii);
        std::set<int> which;
        for (size_t j = 0; j < changed.size(); j++)
        {
            if (changed[j] != values[j])
            {
                which.insert(j);
                EXPECT_TRUE(j == k || std::abs(changed[j] - values[j]) <= 1e-6*std::abs(values[j]))
                    << "parameter " << j << " was " << values[j] << " and is " << changed[j];
            }
        }
        ASSERT_EQ(1, which.count(k)) << "Parameter " << k << " did not change";
        auto reference = evaluate(&fresh, changed, {}, charges);
        if (reference.chi2 != before.chi2)
        {
            relevant.insert(kindOf(oi));
        }
        // As an MCMC step that is accepted
        compare(evaluate(&incremental, changed, which, charges), reference, tolerance);
        // The next step undoes it
        compare(evaluate(&incremental, values, which, charges), before, tolerance);
    }
    std::string irrelevant;
    for (const auto &kind : kinds)
    {
        if (relevant.count(kind) == 0)
        {
            irrelevant += " " + kind;
        }
    }
    printf("%s: %zu of %zu parameters compared with evaluating from scratch, %zu of %zu kinds change the result%s%s\n",
           uc.files.forceField.c_str(), toCompare.size(), optIndex.size(), relevant.size(), kinds.size(),
           irrelevant.empty() ? "" : ", these do not:", irrelevant.c_str());
    // Without changes of the result nothing is tested
    EXPECT_FALSE(relevant.empty());
}

//! The force fields and molecules
std::map<std::string, UpdateCase> updateCases()
{
    const std::vector<eRMS> inter = { eRMS::Interaction, eRMS::Electrostatics, eRMS::Exchange,
                                      eRMS::Dispersion, eRMS::Induction };
    std::map<std::string, UpdateCase> cases;
    {
        // Bonded and non-bonded parameters of monomers, without shells
        UpdateCase uc;
        uc.targets = { eRMS::EPOT, eRMS::MU };
        cases["alcohols"] = uc;
    }
    {
        // Polarizable, with Wang-Buckingham and Gaussian charges, dimers
        UpdateCase uc;
        uc.files.forceField = "forcefields/alc-pol-wbh.xml";
        uc.files.molprop    = "molprop/booze-inter.xml";
        uc.files.selection  = "molprop/booze-inter.dat";
        uc.files.oneH       = false;
        uc.targets          = inter;
        cases["alc-pol-wbh"] = uc;
    }
    {
        // Lennard-Jones 14-7 with virtual sites
        UpdateCase uc;
        uc.files.forceField = "forcefields/4S-lj14-7-pg-hvdw.xml";
        uc.files.molprop    = "molprop/hf-hf.xml";
        uc.files.selection  = "molprop/hf-hf.dat";
        uc.files.oneH       = false;
        uc.targets          = inter;
        cases["hf-lj14-7-vsites"] = uc;
    }
    return cases;
}

class ParameterUpdateTest : public ::testing::Test
{
protected:
    //! For the changed force fields
    gmx::test::TestFileManager fileManager_;
};

// Test: Changing one parameter of each kind in molecules that were
// evaluated before gives the result of evaluating them from scratch, on
// the processor. In double precision the results are the same up to
// rounding.
TEST_F(ParameterUpdateTest, OneAtATime)
{
    for (const auto &uc : updateCases())
    {
        SCOPED_TRACE(uc.first);
        checkUpdates(&fileManager_, uc.second, false, GMX_DOUBLE ? 1e-10 : 1e-5);
    }
}

// Test: The same with the molecules on the GPU, where only the rows of the
// parameter tables that may have changed are made again
TEST_F(ParameterUpdateTest, OneAtATimeGpu)
{
    if (!gmx::test::gpuTestRuns(gpu::useGpu() && !GMX_DOUBLE, "no GPU or double precision"))
    {
        return;
    }
    // Only the polarizable force field, the GPU computes no others
    auto cases = updateCases();
    SCOPED_TRACE("alc-pol-wbh");
    checkUpdates(&fileManager_, cases["alc-pol-wbh"], true, 1e-6);
}

} // namespace

} // namespace alexandria
