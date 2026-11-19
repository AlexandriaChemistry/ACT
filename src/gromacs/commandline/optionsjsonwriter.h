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
 * Declares gmx::OptionsJsonWriter.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_OPTIONSJSONWRITER_H
#define GMX_COMMANDLINE_OPTIONSJSONWRITER_H

#include <string>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class CommandLineHelpContext;
class Options;
class TextWriter;

/*! \internal \brief
 * Collects the options of the commands of a program and writes them as JSON.
 *
 * The document is meant for programs that drive the command line tools, such
 * as the Python interface: it lists for each command its description, its
 * help text and for each option the name, the type, the default value, the
 * allowed values and, for files, whether the file is read or written. The
 * information comes from the same gmx::Options object that is used to print
 * the help and to parse the command line, so that it cannot differ from what
 * the program accepts.
 *
 * The layout of the document is
 * \code{.json}
 * { "program": "alexandria", "format": 1,
 *   "commands": {
 *     "minimize": {
 *       "description": "One line",
 *       "help": "Text without markup",
 *       "options": [
 *         { "name": "maxiter", "type": "int", "description": "...",
 *           "hidden": false, "required": false, "multiple": false,
 *           "default": 1000 } ] } } }
 * \endcode
 * The type is one of bool, int, real, string, enum or file. Options of type
 * enum have an array "values", options of type file an object "file" with
 * the fields input, output, library, allow_missing, directory and
 * extensions. Options with "multiple" set take more than one value and have
 * an array as default; if the number of values is limited, as for a vector,
 * the field "max_values" gives the limit. An option without default value
 * has null as default. A file option that is required has a default file
 * name, which the program uses when the option is not given.
 *
 * \ingroup module_commandline
 */
class OptionsJsonWriter
{
    public:
        //! Version of the layout of the document, increased when it changes
        static const int c_format = 1;

        /*! \brief
         * Initializes the writer.
         *
         * \param[in] binaryName Name of the program, e.g. alexandria.
         */
        explicit OptionsJsonWriter(const std::string &binaryName);
        ~OptionsJsonWriter();

        /*! \brief
         * Returns a writer for the help context.
         *
         * A help context needs a text writer. Nothing that matters is
         * written to it during the export; the text is discarded.
         */
        TextWriter &outputWriter();

        /*! \brief
         * Starts a new command.
         *
         * \param[in] name        Name of the command, e.g. minimize.
         * \param[in] description One-line description, may be nullptr.
         *
         * Must be called before the help of the command is written.
         */
        void startCommand(const std::string &name, const char *description);

        /*! \brief
         * Stores the options and the help text of the current command.
         *
         * \param[in] context  Help context, used to remove the markup from
         *                     the help text and the descriptions.
         * \param[in] options  The options of the command.
         * \param[in] helpText The help text of the command, with markup.
         *
         * Called by CommandLineHelpWriter::writeHelp() when the help
         * context is one for the export of options.
         */
        void writeCommandOptions(const CommandLineHelpContext &context,
                                 const Options                &options,
                                 const std::string            &helpText);

        /*! \brief
         * Finishes the current command.
         *
         * A command for which no options were written, such as the help
         * command itself, is removed from the document: there is nothing a
         * program could do with it.
         */
        void finishCommand();

        //! Returns the document as text, with a newline at the end.
        std::string toString() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
