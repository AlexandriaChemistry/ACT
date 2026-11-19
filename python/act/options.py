#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
"""!
@file
@brief The description of the commands of alexandria and their options.

The description is made by the program itself, with
`alexandria help -export json`, from the same tables that it uses to print
its help and to read its command line. The Python interface uses it to
check arguments before a program is started and to document the functions
in act.commands; a user interface can build its forms from it.
"""
from __future__ import annotations

import difflib
import json
import keyword
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .errors import ActError, ActOptionError

## The version of the layout of the description that this module understands
SUPPORTED_FORMAT = 1

## The types an option can have
OPTION_TYPES = ("bool", "int", "real", "string", "enum", "file")


def python_name(name: str) -> str:
    """!
    @brief The name of an option as a Python keyword argument.
    @param name The name of the option on the command line, without dash.
    @return The same name, with an underscore added if it is a reserved
    word of Python (del becomes del_) and with characters that cannot be
    in a Python name replaced by underscores.
    """
    result = "".join(c if (c.isalnum() or c == "_") else "_" for c in name)
    if not result or result[0].isdigit():
        result = "_" + result
    if keyword.iskeyword(result):
        result += "_"
    return result


@dataclass(frozen=True)
class OptionInfo:
    """!
    @brief One option of a command.
    """
    ## Name on the command line, without the dash
    name: str
    ## One of bool, int, real, string, enum or file
    type: str
    ## Help text
    description: str = ""
    ## Whether the option is left out of the normal help
    hidden: bool = False
    ## Whether the program always uses the option. For a file this means
    ## that the default file name is used when the option is not given.
    required: bool = False
    ## Whether the option takes more than one value
    multiple: bool = False
    ## The largest number of values, None if there is no limit
    max_values: Optional[int] = None
    ## The default value, a list for options with several values
    default: Any = None
    ## The allowed values of an enum
    values: tuple = ()
    ## For a file: whether the program reads it
    is_input: bool = False
    ## For a file: whether the program writes it
    is_output: bool = False
    ## For a file: whether it is looked for in the library of ACT as well
    library: bool = False
    ## For a file: whether the file may be absent
    allow_missing: bool = False
    ## For a file: the extensions it can have, the first is the default
    extensions: tuple = ()

    @property
    def python_name(self) -> str:
        """!
        @brief The name of the option as a Python keyword argument.
        @return The name, see python_name().
        """
        return python_name(self.name)

    @property
    def is_file(self) -> bool:
        """!
        @brief Whether the option is a file name.
        @return True for options of type file.
        """
        return self.type == "file"

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "OptionInfo":
        """!
        @brief Make an option from its description.
        @param data One element of the array "options" of a command.
        @return The option.
        @throws ActError if the description is not understood.
        """
        try:
            name = data["name"]
            otype = data["type"]
        except (KeyError, TypeError) as exc:
            raise ActError("An option in the description of alexandria "
                           "has no name or type") from exc
        if otype not in OPTION_TYPES:
            raise ActError("Option -%s has the unknown type '%s'" % (name, otype))
        fileinfo = data.get("file") or {}
        return OptionInfo(
            name=name,
            type=otype,
            description=data.get("description", ""),
            hidden=bool(data.get("hidden", False)),
            required=bool(data.get("required", False)),
            multiple=bool(data.get("multiple", False)),
            max_values=data.get("max_values"),
            default=data.get("default"),
            values=tuple(data.get("values", ())),
            is_input=bool(fileinfo.get("input", False)),
            is_output=bool(fileinfo.get("output", False)),
            library=bool(fileinfo.get("library", False)),
            allow_missing=bool(fileinfo.get("allow_missing", False)),
            extensions=tuple(fileinfo.get("extensions", ())),
        )


@dataclass(frozen=True)
class CommandInfo:
    """!
    @brief One command of alexandria, such as minimize, with its options.
    """
    ## Name of the command
    name: str
    ## One line that says what the command does
    description: str = ""
    ## The help text, without markup
    help: str = ""
    ## The options by their name on the command line
    options: Dict[str, OptionInfo] = field(default_factory=dict)

    def find_option(self, name: str) -> OptionInfo:
        """!
        @brief Look up an option by its name on the command line or by its
        name as a Python keyword argument.
        @param name The name, without dash.
        @return The option.
        @throws ActOptionError if the command has no such option. The
        message lists options with similar names.
        """
        if name in self.options:
            return self.options[name]
        for option in self.options.values():
            if option.python_name == name:
                return option
        close = difflib.get_close_matches(name, list(self.options), n=3)
        hint = (" Did you mean %s?" % " or ".join("'%s'" % c for c in close)) if close else ""
        raise ActOptionError("alexandria %s has no option '%s'.%s" % (self.name, name, hint))

    def has_option(self, name: str) -> bool:
        """!
        @brief Whether the command has an option.
        @param name The name on the command line, without dash.
        @return True if it has.
        """
        return name in self.options

    @staticmethod
    def from_json(name: str, data: Dict[str, Any]) -> "CommandInfo":
        """!
        @brief Make a command from its description.
        @param name The name of the command.
        @param data The description of the command.
        @return The command.
        """
        options = {}
        for item in data.get("options", []):
            option = OptionInfo.from_json(item)
            options[option.name] = option
        return CommandInfo(name=name,
                           description=data.get("description", ""),
                           help=data.get("help", ""),
                           options=options)


@dataclass(frozen=True)
class ProgramInfo:
    """!
    @brief All the commands of the alexandria program.
    """
    ## Name of the program
    program: str
    ## The commands by name
    commands: Dict[str, CommandInfo] = field(default_factory=dict)

    def find_command(self, name: str) -> CommandInfo:
        """!
        @brief Look up a command.
        @param name The name of the command, such as minimize.
        @return The command.
        @throws ActOptionError if there is no such command.
        """
        if name in self.commands:
            return self.commands[name]
        close = difflib.get_close_matches(name, list(self.commands), n=3)
        hint = (" Did you mean %s?" % " or ".join("'%s'" % c for c in close)) if close else ""
        raise ActOptionError("'%s' is not a command of %s.%s" % (name, self.program, hint))

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "ProgramInfo":
        """!
        @brief Make the description of the program from the JSON document.
        @param data The parsed document.
        @return The description.
        @throws ActError if the document has another layout than this
        module understands.
        """
        if not isinstance(data, dict) or "commands" not in data:
            raise ActError("The description of the options of alexandria is not understood")
        fmt = data.get("format")
        if fmt != SUPPORTED_FORMAT:
            raise ActError("The description of the options of alexandria has format %s, "
                           "but this version of the Python interface reads format %d. "
                           "Use the Python interface that belongs to the program."
                           % (fmt, SUPPORTED_FORMAT))
        commands = {name: CommandInfo.from_json(name, cmd)
                    for name, cmd in data["commands"].items()}
        return ProgramInfo(program=data.get("program", "alexandria"), commands=commands)

    @staticmethod
    def from_file(path: "os.PathLike[str] | str") -> "ProgramInfo":
        """!
        @brief Read the description from a file written by
        `alexandria help -export json`.
        @param path The file.
        @return The description.
        @throws ActError if the file is not valid JSON or not understood.
        """
        try:
            with open(path, "r", encoding="utf-8") as handle:
                data = json.load(handle)
        except (OSError, ValueError) as exc:
            raise ActError("Cannot read the description of the options of alexandria "
                           "from %s: %s" % (path, exc)) from exc
        return ProgramInfo.from_json(data)
