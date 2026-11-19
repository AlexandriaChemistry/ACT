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
@brief Python interface to the Alexandria Chemistry Toolkit (ACT).

The commands of the alexandria program can be run from Python, with their
options as keyword arguments and their results as Python objects:
@code{.py}
import act

res = act.commands.minimize(ff="ACS-pg1.xml", f="acetone.sdf")
print(res.data["rmsd"])           # {'value': 0.019, 'unit': 'nm'}
print(res.file("c"))              # the minimized structure

print(act.commands.available())   # the commands
help(act.commands.nma)            # the options of a command
@endcode
A program that fails raises act.ActError with the reason; a wrong option or
a missing input file raises an error before anything is started.
"""
from __future__ import annotations

import warnings
from typing import Any

from . import commands
from .errors import (ActCancelledError, ActError, ActNotFoundError, ActOptionError,
                     ActTimeoutError)
from .options import CommandInfo, OptionInfo, ProgramInfo
from .runner import Result, build_arguments, config, program_info, run

__all__ = [
    "ActCancelledError", "ActError", "ActNotFoundError", "ActOptionError", "ActTimeoutError",
    "CommandInfo", "OptionInfo", "ProgramInfo",
    "Result", "build_arguments", "commands", "config", "program_info", "run",
]


def __getattr__(name: str) -> Any:
    """!
    @brief Give the class ACT of the earlier module act.py, with a warning.
    @param name The name that was asked for.
    @return The class ACT, Alg or Target of act.legacy.
    @throws AttributeError for any other name.
    """
    if name in ("ACT", "Alg", "Target"):
        from . import legacy  # pylint: disable=import-outside-toplevel
        warnings.warn("act.%s has moved to act.legacy.%s and will be removed from act "
                      "in the next release" % (name, name), DeprecationWarning, stacklevel=2)
        return getattr(legacy, name)
    raise AttributeError("module 'act' has no attribute '%s'" % name)
