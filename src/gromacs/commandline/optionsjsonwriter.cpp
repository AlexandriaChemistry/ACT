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
 * Implements gmx::OptionsJsonWriter.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_commandline
 */
#include "actpre.h"

#include "optionsjsonwriter.h"

#include <cstdlib>

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

namespace
{

//! The kinds of value an option can have in the JSON document
enum class JsonValueType
{
    Bool, Int, Real, String
};

/*! \brief
 * Converts the text of a default value to a JSON value.
 *
 * \param[in] text The default value as the help prints it.
 * \param[in] type The kind of value.
 * \returns   A JSON number, boolean or string. Text that is not a number
 *     gives null for the numeric types.
 */
nlohmann::json toJsonValue(const std::string &text, JsonValueType type)
{
    switch (type)
    {
        case JsonValueType::Bool:
            return text == "yes" || text == "true" || text == "1";
        case JsonValueType::Int:
        {
            char     *end   = nullptr;
            long long value = std::strtoll(text.c_str(), &end, 10);
            if (text.empty() || *end != '\0')
            {
                return nullptr;
            }
            return value;
        }
        case JsonValueType::Real:
        {
            char  *end   = nullptr;
            double value = std::strtod(text.c_str(), &end);
            if (text.empty() || *end != '\0')
            {
                return nullptr;
            }
            return value;
        }
        case JsonValueType::String:
            break;
    }
    return text;
}

/*! \internal \brief
 * Visitor that converts each option to a JSON object.
 *
 * \ingroup module_commandline
 */
class OptionsJsonVisitor : public OptionsVisitor
{
    public:
        /*! \brief
         * Initializes the visitor.
         *
         * \param[in] context Context that removes markup from descriptions.
         */
        explicit OptionsJsonVisitor(const HelpWriterContext &context)
            : context_(context), result_(nlohmann::json::array())
        {
        }

        //! Returns the array with one object per option.
        const nlohmann::json &result() const { return result_; }

        void visitSection(const OptionSectionInfo &section) override
        {
            OptionsIterator iterator(section);
            iterator.acceptSections(this);
            iterator.acceptOptions(this);
        }
        void visitOption(const OptionInfo &option) override;

    private:
        //! Removes the markup from a text, without wrapping the lines.
        std::string plainText(const std::string &text) const
        {
            return context_.substituteMarkupAndWrapToString(TextLineWrapperSettings(), text);
        }

        //! Context for the conversion of markup
        const HelpWriterContext &context_;
        //! The options converted so far
        nlohmann::json           result_;
};

void OptionsJsonVisitor::visitOption(const OptionInfo &option)
{
    nlohmann::json json;
    json["name"]        = option.name();
    json["description"] = plainText(option.formatDescription());
    json["hidden"]      = option.isHidden();
    json["required"]    = option.isRequired();

    // The type in the document and the kind of JSON value for defaults.
    std::string   type      = "string";
    JsonValueType valueType = JsonValueType::String;
    if (const auto *fileOption = option.toType<FileNameOptionInfo>())
    {
        type = "file";
        nlohmann::json file;
        file["input"]         = fileOption->isInputFile() || fileOption->isInputOutputFile();
        file["output"]        = fileOption->isOutputFile() || fileOption->isInputOutputFile();
        file["library"]       = fileOption->isLibraryFile();
        file["allow_missing"] = fileOption->allowMissing();
        file["directory"]     = fileOption->isDirectoryOption();
        file["extensions"]    = nlohmann::json::array();
        for (const char *extension : fileOption->extensions())
        {
            file["extensions"].push_back(extension);
        }
        json["file"] = file;
    }
    else if (option.isType<BooleanOptionInfo>())
    {
        type      = "bool";
        valueType = JsonValueType::Bool;
    }
    else if (option.isType<IntegerOptionInfo>() || option.isType<Int64OptionInfo>())
    {
        type      = "int";
        valueType = JsonValueType::Int;
    }
    else if (option.isType<FloatOptionInfo>() || option.isType<DoubleOptionInfo>())
    {
        type      = "real";
        valueType = JsonValueType::Real;
    }
    else if (const auto *enumOption = option.toType<EnumOptionInfo>())
    {
        type           = "enum";
        json["values"] = enumOption->allowedValues();
    }
    else if (const auto *stringOption = option.toType<StringOptionInfo>())
    {
        if (stringOption->isEnumerated())
        {
            type           = "enum";
            json["values"] = stringOption->allowedValues();
        }
    }
    json["type"] = type;

    // How many values the option takes. A boolean takes none on the command
    // line, the other types one unless stated otherwise.
    // A negative maximum means that there is no limit.
    const int  maxCount  = option.maxValueCount();
    const bool bMultiple = (valueType != JsonValueType::Bool && maxCount != 1);
    json["multiple"] = bMultiple;
    if (bMultiple && maxCount > 1)
    {
        json["max_values"] = maxCount;
    }

    nlohmann::json defaults = nlohmann::json::array();
    if (const auto *boolOption = option.toType<BooleanOptionInfo>())
    {
        defaults.push_back(boolOption->defaultValue());
    }
    else
    {
        for (const std::string &text : option.defaultValuesAsStrings())
        {
            defaults.push_back(toJsonValue(text, valueType));
        }
    }
    if (bMultiple)
    {
        json["default"] = defaults;
    }
    else if (!defaults.empty())
    {
        json["default"] = defaults[0];
    }
    else
    {
        json["default"] = nullptr;
    }
    result_.push_back(json);
}

}   // namespace

/*! \internal \brief
 * Private implementation class for OptionsJsonWriter.
 *
 * \ingroup module_commandline
 */
class OptionsJsonWriter::Impl
{
    public:
        //! Initializes the document for a program.
        explicit Impl(const std::string &binaryName)
            : discardWriter_(&discardStream_)
        {
            document_["program"]  = binaryName;
            document_["format"]   = OptionsJsonWriter::c_format;
            document_["commands"] = nlohmann::json::object();
        }

        //! Receives text written through the help context
        StringOutputStream discardStream_;
        //! Writer on the stream above
        TextWriter         discardWriter_;
        //! The document
        nlohmann::json     document_;
        //! Name of the command that is being exported
        std::string        currentCommand_;
};

OptionsJsonWriter::OptionsJsonWriter(const std::string &binaryName)
    : impl_(new Impl(binaryName))
{
}

OptionsJsonWriter::~OptionsJsonWriter()
{
}

TextWriter &OptionsJsonWriter::outputWriter()
{
    return impl_->discardWriter_;
}

void OptionsJsonWriter::startCommand(const std::string &name,
                                     const char        *description)
{
    impl_->currentCommand_ = name;
    nlohmann::json command;
    command["description"]             = description != nullptr ? description : "";
    impl_->document_["commands"][name] = command;
}

void OptionsJsonWriter::writeCommandOptions(const CommandLineHelpContext &context,
                                            const Options                &options,
                                            const std::string            &helpText)
{
    GMX_RELEASE_ASSERT(!impl_->currentCommand_.empty(),
                       "startCommand() must be called before the options are written");
    const HelpWriterContext &writerContext = context.writerContext();
    OptionsJsonVisitor       visitor(writerContext);
    visitor.visitSection(options.rootSection());

    nlohmann::json          &command = impl_->document_["commands"][impl_->currentCommand_];
    command["help"] = writerContext.substituteMarkupAndWrapToString(TextLineWrapperSettings(),
                                                                    helpText);
    command["options"] = visitor.result();
}

void OptionsJsonWriter::finishCommand()
{
    nlohmann::json &commands = impl_->document_["commands"];
    auto            command  = commands.find(impl_->currentCommand_);
    if (command != commands.end() && !command->contains("options"))
    {
        commands.erase(command);
    }
    impl_->currentCommand_.clear();
}

std::string OptionsJsonWriter::toString() const
{
    // Replace bytes that are not valid UTF-8 instead of throwing
    return impl_->document_.dump(1, ' ', false,
                                 nlohmann::json::error_handler_t::replace) + "\n";
}

} // namespace gmx
