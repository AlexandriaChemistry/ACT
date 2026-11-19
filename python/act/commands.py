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
@brief One Python function for each command of alexandria.

The functions are made from the description that the program gives of its
own commands, the first time one of them is used. So
@code{.py}
import act
res = act.commands.minimize(ff="ACS-pg1.xml", f="acetone.sdf")
help(act.commands.minimize)
@endcode
runs `alexandria minimize -ff ACS-pg1.xml -f acetone.sdf`, and the help
lists the options of the program with their types, defaults and help texts.
All functions take the arguments of act.run() as well: workdir, processes,
timeout, check and env.
"""
from __future__ import annotations

import inspect
import textwrap
from typing import Any, Callable, Dict, List

from . import runner
from .options import CommandInfo, OptionInfo

## The functions that were made, by name of the command
_functions: Dict[str, Callable[..., runner.Result]] = {}

## The Python types that the options take, for the signatures
_ANNOTATIONS = {"bool": "bool", "int": "int", "real": "float",
                "string": "str", "enum": "str", "file": "path"}


def _describe(option: OptionInfo) -> str:
    """!
    @brief The lines of the help text for one option.
    @param option The option.
    @return The text, indented for the docstring of a function.
    """
    kind = _ANNOTATIONS[option.type]
    if option.is_file:
        direction = "input" if option.is_input and not option.is_output else \
            "output" if option.is_output and not option.is_input else "input and output"
        kind = "%s file %s" % (direction, "/".join(option.extensions))
    if option.type == "enum":
        kind = "one of " + ", ".join(option.values)
    if option.multiple:
        kind = "list of " + kind
    head = "%s (%s" % (option.python_name, kind)
    if option.default is not None and option.default != []:
        head += ", default %s" % (option.default,)
    head += ")"
    body = textwrap.fill(option.description or "No description.", width=72,
                         initial_indent="        ", subsequent_indent="        ")
    return "    %s\n%s" % (head, body)


def _documentation(command: CommandInfo) -> str:
    """!
    @brief The docstring of the function for a command.
    @param command The command.
    @return The description, the help text and the list of options.
    """
    parts = [command.description, ""]
    if command.help:
        for paragraph in command.help.split("\n\n"):
            parts.append(textwrap.fill(" ".join(paragraph.split()), width=76))
            parts.append("")
    visible = [o for o in command.options.values() if not o.hidden]
    hidden = [o for o in command.options.values() if o.hidden]
    parts.append("Options:")
    parts.extend(_describe(option) for option in visible)
    if hidden:
        parts.append("")
        parts.append("Options for experts:")
        parts.extend(_describe(option) for option in hidden)
    parts.append("")
    parts.append("Arguments of act.run(): workdir, processes, timeout, check, env.")
    parts.append("Returns an act.Result.")
    return "\n".join(parts)


def _make_function(command: CommandInfo) -> Callable[..., runner.Result]:
    """!
    @brief Make the function for a command.
    @param command The command.
    @return A function that runs the command with keyword arguments.
    """
    name = command.name

    def function(**options: Any) -> runner.Result:
        return runner.run(name, **options)

    function.__name__ = name
    function.__qualname__ = name
    function.__doc__ = _documentation(command)
    # A signature with the real options, for help() and for editors
    parameters: List[inspect.Parameter] = []
    for option in command.options.values():
        parameters.append(inspect.Parameter(option.python_name,
                                            inspect.Parameter.KEYWORD_ONLY,
                                            default=None))
    for extra in runner.RUN_ARGUMENTS:
        default = inspect.signature(runner.run).parameters[extra].default
        parameters.append(inspect.Parameter(extra, inspect.Parameter.KEYWORD_ONLY,
                                            default=default))
    function.__signature__ = inspect.Signature(  # type: ignore[attr-defined]
        parameters, return_annotation=runner.Result)
    return function


def available() -> List[str]:
    """!
    @brief The names of the commands.
    @return The sorted names, for which this module has functions.
    """
    return sorted(runner.program_info().commands)


def __getattr__(name: str) -> Callable[..., runner.Result]:
    """!
    @brief Make the function for a command when it is first asked for.
    @param name The name of the command.
    @return The function.
    @throws AttributeError if alexandria has no such command.
    """
    if name.startswith("_"):
        raise AttributeError(name)
    if name not in _functions:
        info = runner.program_info()
        if name not in info.commands:
            raise AttributeError("alexandria has no command '%s', the commands are: %s"
                                 % (name, ", ".join(sorted(info.commands))))
        _functions[name] = _make_function(info.commands[name])
    return _functions[name]


def __dir__() -> List[str]:
    """!
    @brief The names in this module, for completion at the Python prompt.
    @return The commands, if the program can be found, and available.
    """
    names = ["available"]
    try:
        names.extend(available())
    except Exception:  # pylint: disable=broad-except
        # Without the program there are no commands to list
        pass
    return names
