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
 * Helpers for the tests of the train_ff module, which run the module
 * with input derived from the benchmarks in ACTdata/Performance.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_TRAINFF_TEST_UTIL_H
#define ACT_TRAINFF_TEST_UTIL_H

#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include "act/alexandria/alex_modules.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace trainfftest
{

//! Wait for all MPI ranks
inline void barrier()
{
#if GMX_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/*! \brief Change the working directory for the lifetime of the object
 * The module writes some files, like the correlation plots of the
 * energy components, to the working directory.
 */
class ScopedWorkingDirectory
{
private:
    //! The original directory
    std::filesystem::path old_;
public:
    //! Change to \p dir
    explicit ScopedWorkingDirectory(const std::filesystem::path &dir)
        : old_(std::filesystem::current_path())
    {
        std::filesystem::current_path(dir);
    }
    //! Change back
    ~ScopedWorkingDirectory() { std::filesystem::current_path(old_); }
};

//! The output of a train_ff run
struct TrainOutput
{
    //! The return value of the module
    int                      status      = -1;
    //! The json output
    nlohmann::json           json;
    //! The chi2 of the training set at the start, if it is in the log
    double                   initialChi2 = -1;
    //! The lowest chi2 of the training set
    double                   finalChi2   = -1;
    //! The lines of the fitness file of the training set
    std::vector<std::string> fitness;
    //! The output force field
    std::string              forceField;
};

//! \return a number that is written as a number or a string in json
inline double jsonNumber(const nlohmann::json &value)
{
    return value.is_string() ? std::stod(value.get<std::string>()) : value.get<double>();
}

//! \return the lines of a text file
inline std::vector<std::string> readLines(const std::string &fileName)
{
    std::vector<std::string> lines;
    std::ifstream            in(fileName);
    std::string              line;
    while (std::getline(in, line))
    {
        lines.push_back(line);
    }
    return lines;
}

/*! \brief Find the chi2 of the training set in the log file
 * \param[in]  log     The lines of the log file
 * \param[out] output  The initial and final chi2
 */
inline void readChi2(const std::vector<std::string> &log, TrainOutput *output)
{
    bool initial = false;
    for (const auto &line : log)
    {
        if (line.find("initial parameter vector chi2 components") != std::string::npos)
        {
            initial = true;
        }
        else if (initial && line.rfind("TOT", 0) == 0 && line.find("Train") != std::string::npos)
        {
            auto pos = line.find("weighted:");
            if (pos != std::string::npos)
            {
                output->initialChi2 = std::stod(line.substr(pos + 9));
            }
            initial = false;
        }
        auto pos = line.find("Minimum chi2 for Train");
        if (pos != std::string::npos)
        {
            output->finalChi2 = std::stod(line.substr(pos + 22));
        }
    }
}

/*! \brief Run train_ff
 * The input files are in the trainff directory of the tests. Only the
 * master reads the output, since only it writes the output files.
 * \param[in] fileManager To find input and temporary files
 * \param[in] label       Unique label for the output files
 * \param[in] ffName      Force field file
 * \param[in] molName     Base name of the molprop and selection file
 * \param[in] options     Further options and their values, an empty
 *                        value for a boolean option
 * \return the output
 */
inline TrainOutput runTrainFF(gmx::test::TestFileManager                           *fileManager,
                              const std::string                                    &label,
                              const std::string                                    &ffName,
                              const std::string                                    &molName,
                              const std::vector<std::pair<std::string, std::string>> &options)
{
    gmx::test::CommandLine cmd;
    cmd.append("train_ff");
    cmd.addOption("-ff", fileManager->getInputFilePath("trainff/" + ffName));
    cmd.addOption("-mp", fileManager->getInputFilePath("trainff/" + molName + ".xml"));
    cmd.addOption("-sel", fileManager->getInputFilePath("trainff/" + molName + ".dat"));
    std::map<std::string, std::string> out = {
        { "-o",       label + "-ff.xml" },
        { "-otest",   label + "-fftest.xml" },
        { "-g",       label + ".log" },
        { "-json",    label + ".json" },
        { "-fitness", label + "-fitness.dat" },
        { "-gpout",   label + "-pool.csv" },
        { "-conv",    label + "-conv.xvg" },
        { "-chi2",    label + "-chi2.xvg" }
    };
    for (auto &o : out)
    {
        o.second = fileManager->getTemporaryFilePath(o.second);
        cmd.addOption(o.first.c_str(), o.second);
    }
    for (const auto &opt : options)
    {
        if (opt.second.empty())
        {
            cmd.append(opt.first);
        }
        else
        {
            cmd.addOption(opt.first.c_str(), opt.second);
        }
    }
    TrainOutput output;
    {
        std::filesystem::path dir = std::filesystem::path(out["-g"]).parent_path();
        ScopedWorkingDirectory cwd(dir);
        output.status = train_ff(cmd.argc(), cmd.argv());
    }
    if (gmx_node_rank() == 0)
    {
        std::ifstream in(out["-json"]);
        EXPECT_TRUE(in.good()) << "Cannot open " << out["-json"];
        if (in.good())
        {
            output.json = nlohmann::json::parse(in);
        }
        readChi2(readLines(out["-g"]), &output);
        // The fitness files are named after the data set
        auto fitness = out["-fitness"];
        output.fitness = readLines(fitness.substr(0, fitness.rfind(".")) + "_Train.dat");
        std::ifstream ff(out["-o"]);
        std::stringstream buffer;
        buffer << ff.rdbuf();
        output.forceField = buffer.str();
    }
    return output;
}

/*! \brief Check the output against reference data
 * The best parameters and the root mean square deviations of the
 * training set are checked.
 * \param[in] checker The checker
 * \param[in] output  The output
 */
inline void checkOutput(gmx::test::TestReferenceChecker *checker,
                        const TrainOutput               &output)
{
    const auto &tff = output.json["train_ff"];
    ASSERT_TRUE(tff.contains("parameters"));
    auto params = checker->checkCompound("Parameters", "Best");
    for (const auto &p : tff["parameters"])
    {
        auto name = p["class"].get<std::string>() + " " + p["name"].get<std::string>();
        params.checkReal(jsonNumber(p["best"]), name.c_str());
    }
    ASSERT_TRUE(tff.contains("statistics"));
    auto stats = checker->checkCompound("Statistics", "Train");
    for (const auto &s : tff["statistics"]["Train"].items())
    {
        stats.checkReal(jsonNumber(s.value()["RMSD"]), (s.key() + " RMSD").c_str());
    }
    checker->checkReal(output.finalChi2, "Train chi2");
}

} // namespace trainfftest

} // namespace alexandria

#endif
