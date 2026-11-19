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
 * Tests the description of the options of all alexandria commands as JSON,
 * written by "alexandria help -export json". The Python interface and the
 * graphical user interface are built from this description.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <gtest/gtest.h>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/utility/fileredirector.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textstream.h"
#include "testutils/cmdlinetest.h"

namespace alexandria
{

namespace
{

/*! \brief Output redirector that keeps what a program writes in memory.
 *
 * Both the standard output of the help module and the files it exports are
 * stored as strings, such that a test can inspect them without touching the
 * file system.
 */
class CapturingRedirector : public gmx::IFileOutputRedirector
{
public:
    //! \return the stream that collects the standard output
    gmx::TextOutputStream &standardOutput() override
    {
        return stdout_;
    }
    /*! \brief Open a file that is kept in memory
     * \param[in] filename The name of the file
     * \return the stream that collects the content of the file
     */
    gmx::TextOutputStreamPointer openTextOutputFile(const char *filename) override
    {
        auto stream = std::make_shared<gmx::StringOutputStream>();
        files_[filename] = stream;
        return stream;
    }
    //! \return what was written to the standard output
    const std::string &standardOutputText() const
    {
        return stdout_.toString();
    }
    /*! \brief Look up a file
     * \param[in] filename The name of the file
     * \return whether a file with this name was written
     */
    bool hasFile(const std::string &filename) const
    {
        return files_.find(filename) != files_.end();
    }
    /*! \brief Content of a file
     * \param[in] filename The name of the file
     * \return what was written to the file
     */
    const std::string &file(const std::string &filename) const
    {
        return files_.at(filename)->toString();
    }
private:
    //! Collects the standard output
    gmx::StringOutputStream                                         stdout_;
    //! Collects the files, by name
    std::map<std::string, std::shared_ptr<gmx::StringOutputStream>> files_;
};

/*! \brief Run alexandria with a command line and keep the output in memory
 * \param[in]  cmdline    The command line, starting with the program name
 * \param[out] redirector Receives standard output and exported files
 * \return the return code of the program
 */
int runAlexandria(const std::vector<const char *> &cmdline,
                  CapturingRedirector             *redirector)
{
    gmx::test::CommandLine         args(cmdline);
    gmx::CommandLineProgramContext context(args.argc(), args.argv());
    gmx::CommandLineModuleManager  manager("alexandria", &context);
    registerAlexandriaModules(&manager);
    manager.setQuiet(true);
    manager.setOutputRedirector(redirector);
    return manager.run(args.argc(), args.argv());
}

//! Name of the file written by "alexandria help -export json"
const char *const optionsFile = "alexandria-options.json";

//! All the commands of alexandria
const std::vector<std::string> allCommands = {
    "b2", "edit_ff", "edit_mp", "gen_ff", "gentop", "geometry_ff",
    "merge_ff", "min_complex", "minimize", "nma", "simulate", "train_ff"
};

/*! \brief Export the options of all commands and parse the result
 * \return the JSON document, or a null value if the export failed
 */
nlohmann::json exportOptions()
{
    CapturingRedirector redirector;
    int                 rc = -1;
    EXPECT_NO_THROW(rc = runAlexandria({ "alexandria", "help", "-export", "json" },
                                       &redirector));
    EXPECT_EQ(0, rc);
    EXPECT_TRUE(redirector.hasFile(optionsFile));
    if (!redirector.hasFile(optionsFile))
    {
        return nlohmann::json();
    }
    // Do not throw on invalid JSON, but return a discarded value
    return nlohmann::json::parse(redirector.file(optionsFile), nullptr, false);
}

/*! \brief Find an option of a command
 * \param[in] doc     The JSON document with all commands
 * \param[in] command The name of the command
 * \param[in] name    The name of the option, without dash
 * \return the JSON object of the option, or a null value if absent
 */
nlohmann::json findOption(const nlohmann::json &doc,
                          const std::string    &command,
                          const std::string    &name)
{
    if (!doc.is_object() || !doc.contains("commands") ||
        !doc["commands"].contains(command))
    {
        return nlohmann::json();
    }
    for (const auto &opt : doc["commands"][command]["options"])
    {
        if (opt.value("name", "") == name)
        {
            return opt;
        }
    }
    return nlohmann::json();
}

/*! \brief Extract the option names from the synopsis of a help text
 * \param[in] help The console help of a command
 * \return the names of the options, without dash and without [no]
 */
std::set<std::string> optionsInSynopsis(const std::string &help)
{
    std::set<std::string> names;
    auto                  start = help.find("SYNOPSIS");
    auto                  end   = help.find("DESCRIPTION");
    if (start == std::string::npos)
    {
        return names;
    }
    if (end == std::string::npos)
    {
        end = help.find("OPTIONS");
    }
    std::string synopsis = help.substr(start, end == std::string::npos ? end : end - start);
    // Options are written as -name or -[no]name, after a space or a bracket.
    const std::regex re("(?:^|[\\s\\[])-(?:\\[no\\])?([A-Za-z0-9_]+)");
    for (auto it = std::sregex_iterator(synopsis.begin(), synopsis.end(), re);
         it != std::sregex_iterator(); ++it)
    {
        names.insert((*it)[1].str());
    }
    return names;
}

TEST(OptionsJsonTest, ExportIsValidJsonWithAllCommands)
{
    auto doc = exportOptions();
    ASSERT_TRUE(doc.is_object());
    EXPECT_EQ("alexandria", doc.value("program", ""));
    EXPECT_EQ(1, doc.value("format", 0));
    ASSERT_TRUE(doc.contains("commands"));
    for (const auto &command : allCommands)
    {
        ASSERT_TRUE(doc["commands"].contains(command)) << command;
        const auto &cmd = doc["commands"][command];
        EXPECT_FALSE(cmd.value("description", "").empty()) << command;
        EXPECT_FALSE(cmd.value("help", "").empty()) << command;
        ASSERT_TRUE(cmd.contains("options")) << command;
        EXPECT_TRUE(cmd["options"].is_array()) << command;
        EXPECT_FALSE(cmd["options"].empty()) << command;
    }
    EXPECT_EQ(allCommands.size(), doc["commands"].size());
}

TEST(OptionsJsonTest, EveryOptionIsDescribed)
{
    auto                        doc = exportOptions();
    ASSERT_TRUE(doc.is_object());
    const std::set<std::string> knownTypes = {
        "bool", "int", "real", "string", "enum", "file"
    };
    for (const auto &command : allCommands)
    {
        ASSERT_TRUE(doc["commands"].contains(command)) << command;
        std::set<std::string> names;
        for (const auto &opt : doc["commands"][command]["options"])
        {
            const std::string name = opt.value("name", "");
            const std::string what = command + " -" + name;
            EXPECT_FALSE(name.empty()) << command;
            EXPECT_TRUE(names.insert(name).second) << "duplicate " << what;
            ASSERT_TRUE(opt.contains("type")) << what;
            const std::string type = opt["type"].get<std::string>();
            EXPECT_EQ(1U, knownTypes.count(type)) << what << " has type " << type;
            ASSERT_TRUE(opt.contains("description")) << what;
            EXPECT_TRUE(opt["description"].is_string()) << what;
            ASSERT_TRUE(opt.contains("hidden")) << what;
            EXPECT_TRUE(opt["hidden"].is_boolean()) << what;
            ASSERT_TRUE(opt.contains("required")) << what;
            EXPECT_TRUE(opt["required"].is_boolean()) << what;
            if (type == "enum")
            {
                ASSERT_TRUE(opt.contains("values")) << what;
                EXPECT_FALSE(opt["values"].empty()) << what;
            }
            if (type == "file")
            {
                ASSERT_TRUE(opt.contains("file")) << what;
                const auto &file = opt["file"];
                EXPECT_TRUE(file.value("input", false) || file.value("output", false)) << what;
                ASSERT_TRUE(file.contains("extensions")) << what;
                EXPECT_FALSE(file["extensions"].empty()) << what;
            }
            // Defaults have the JSON type that belongs to the option type
            if (opt.contains("default") && !opt["default"].is_null() &&
                !opt["default"].is_array())
            {
                const auto &def = opt["default"];
                if (type == "bool")
                {
                    EXPECT_TRUE(def.is_boolean()) << what;
                }
                else if (type == "int")
                {
                    EXPECT_TRUE(def.is_number_integer()) << what;
                }
                else if (type == "real")
                {
                    EXPECT_TRUE(def.is_number()) << what;
                }
                else
                {
                    EXPECT_TRUE(def.is_string()) << what;
                }
            }
        }
    }
}

TEST(OptionsJsonTest, ContainsEveryOptionOfHelp)
{
    auto doc = exportOptions();
    ASSERT_TRUE(doc.is_object());
    for (const auto &command : allCommands)
    {
        CapturingRedirector redirector;
        int                 rc = -1;
        ASSERT_NO_THROW(rc = runAlexandria({ "alexandria", command.c_str(), "-h", "-hidden" },
                                           &redirector));
        ASSERT_EQ(0, rc);
        const auto inHelp = optionsInSynopsis(redirector.standardOutputText());
        ASSERT_FALSE(inHelp.empty()) << command;
        ASSERT_TRUE(doc["commands"].contains(command)) << command;
        std::set<std::string> inJson;
        for (const auto &opt : doc["commands"][command]["options"])
        {
            inJson.insert(opt.value("name", ""));
        }
        for (const auto &name : inHelp)
        {
            EXPECT_EQ(1U, inJson.count(name))
            << command << " -" << name << " is in the help but not in the JSON";
        }
        for (const auto &name : inJson)
        {
            EXPECT_EQ(1U, inHelp.count(name))
            << command << " -" << name << " is in the JSON but not in the help";
        }
    }
}

TEST(OptionsJsonTest, HiddenOptionsAreMarked)
{
    auto doc = exportOptions();
    ASSERT_TRUE(doc.is_object());
    for (const auto &command : allCommands)
    {
        CapturingRedirector redirector;
        int                 rc = -1;
        ASSERT_NO_THROW(rc = runAlexandria({ "alexandria", command.c_str(), "-h" },
                                           &redirector));
        ASSERT_EQ(0, rc);
        const auto visible = optionsInSynopsis(redirector.standardOutputText());
        ASSERT_TRUE(doc["commands"].contains(command)) << command;
        for (const auto &opt : doc["commands"][command]["options"])
        {
            const std::string name = opt.value("name", "");
            EXPECT_EQ(visible.count(name) == 0, opt.value("hidden", false))
            << command << " -" << name;
        }
    }
}

TEST(OptionsJsonTest, KnownOptionsOfMinimize)
{
    const auto doc = exportOptions();
    ASSERT_TRUE(doc.is_object());

    auto ff = findOption(doc, "minimize", "ff");
    ASSERT_TRUE(ff.is_object());
    EXPECT_EQ("file", ff.value("type", ""));
    EXPECT_TRUE(ff["file"].value("input", false));
    EXPECT_FALSE(ff["file"].value("output", true));
    EXPECT_EQ(".xml", ff["file"]["extensions"][0].get<std::string>());

    auto log = findOption(doc, "minimize", "g");
    ASSERT_TRUE(log.is_object());
    EXPECT_TRUE(log["file"].value("output", false));
    EXPECT_FALSE(log["file"].value("input", true));

    auto maxiter = findOption(doc, "minimize", "maxiter");
    ASSERT_TRUE(maxiter.is_object());
    EXPECT_EQ("int", maxiter.value("type", ""));
    EXPECT_TRUE(maxiter["default"].is_number_integer());
    EXPECT_FALSE(maxiter.value("description", "").empty());

    auto toler = findOption(doc, "minimize", "toler");
    ASSERT_TRUE(toler.is_object());
    EXPECT_EQ("real", toler.value("type", ""));
    EXPECT_TRUE(toler["default"].is_number());

    auto minalg = findOption(doc, "minimize", "minalg");
    ASSERT_TRUE(minalg.is_object());
    EXPECT_EQ("enum", minalg.value("type", ""));
    ASSERT_TRUE(minalg["default"].is_string());
    bool found = false;
    for (const auto &value : minalg["values"])
    {
        found = found || (value == minalg["default"]);
    }
    EXPECT_TRUE(found) << "The default of -minalg is one of its values";

    auto json = findOption(doc, "minimize", "json");
    ASSERT_TRUE(json.is_object());
    EXPECT_EQ("bool", json.value("type", ""));
    EXPECT_TRUE(json["default"].is_boolean());

    auto f = findOption(doc, "minimize", "f");
    ASSERT_TRUE(f.is_object());
    EXPECT_EQ("string", f.value("type", ""));
    EXPECT_TRUE(f["default"].is_null());
}

TEST(OptionsJsonTest, OptionsWithSeveralValues)
{
    const auto doc = exportOptions();
    ASSERT_TRUE(doc.is_object());

    // A vector of three numbers
    auto box = findOption(doc, "gentop", "box");
    ASSERT_TRUE(box.is_object());
    EXPECT_EQ("real", box.value("type", ""));
    EXPECT_TRUE(box.value("multiple", false));
    EXPECT_EQ(3, box.value("max_values", 0));
    ASSERT_TRUE(box["default"].is_array());
    EXPECT_EQ(3U, box["default"].size());

    // Any number of force field files
    auto ff = findOption(doc, "train_ff", "ff");
    ASSERT_TRUE(ff.is_object());
    EXPECT_TRUE(ff.value("multiple", false));
    EXPECT_FALSE(ff.contains("max_values"));
    EXPECT_TRUE(ff["default"].is_array());

    // One force field file
    auto ff1 = findOption(doc, "minimize", "ff");
    ASSERT_TRUE(ff1.is_object());
    EXPECT_FALSE(ff1.value("multiple", true));
    EXPECT_TRUE(ff1["default"].is_string());
}

TEST(OptionsJsonTest, HelpTextHasNoMarkup)
{
    auto doc = exportOptions();
    ASSERT_TRUE(doc.is_object());
    for (const auto &command : allCommands)
    {
        ASSERT_TRUE(doc["commands"].contains(command)) << command;
        const auto &cmd  = doc["commands"][command];
        const std::string help = cmd.value("help", "");
        for (const char *markup : { "[TT]", "[tt]", "[PAR]", "[BR]", "[THISMODULE]" })
        {
            EXPECT_EQ(std::string::npos, help.find(markup)) << command << " help has " << markup;
            for (const auto &opt : cmd["options"])
            {
                EXPECT_EQ(std::string::npos, opt.value("description", "").find(markup))
                << command << " -" << opt.value("name", "") << " has " << markup;
            }
        }
    }
}

} // namespace

} // namespace alexandria
