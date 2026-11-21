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
from .errors import (ActCancelledError, ActError, ActFileError, ActInputError,
                     ActInternalError, ActNotFoundError, ActOptionError, ActTimeoutError)
from .options import CommandInfo, OptionInfo, ProgramInfo
from .runner import Result, build_arguments, config, program_info, run

# The classes use the code of ACT directly, see the module _core. They are
# there when ACT was built with -DACT_PYTHON=ON, and the rest of the package
# works without them.
try:
    from . import _core
except ImportError:
    _core = None  # type: ignore[assignment]

## The classes that the module _core has
_CLASSES = ("ForceField", "Minimization", "Molecule", "NormalModes", "Parameter",
            "ParticleType", "Trajectory")


def classes_available() -> bool:
    """!
    @brief Whether this package has the classes of ACT, such as
    act.ForceField, which use the code of ACT directly instead of running
    the alexandria program.
    @return True if ACT was built with them.
    """
    return _core is not None


def data_directory() -> str:
    """!
    @brief The directory of the installation in which the classes of ACT
    look for the data files of ACT, such as the properties of the elements.
    The directory in the environment variable ACTDATA, if set, comes first.
    @return The directory, empty if this installation has none or was
    built without the classes.
    """
    return _core.data_directory() if _core is not None else ""


__all__ = [
    "ActCancelledError", "ActError", "ActFileError", "ActInputError", "ActInternalError",
    "ActNotFoundError", "ActOptionError", "ActTimeoutError",
    "classes_available", "data_directory",
    "CommandInfo", "OptionInfo", "ProgramInfo",
    "Result", "build_arguments", "commands", "config", "program_info", "run",
]


def __getattr__(name: str) -> Any:
    """!
    @brief Give a class of ACT, such as act.ForceField, or the class ACT of
    the earlier module act.py, with a warning.
    @param name The name that was asked for.
    @return The class of the module _core, or ACT, Alg or Target of act.legacy.
    @throws ImportError for a class of ACT when ACT was built without them.
    @throws AttributeError for any other name.
    """
    if name in _CLASSES:
        if _core is None:
            raise ImportError("act.%s needs the Python classes of ACT, which this installation "
                              "does not have: build ACT with -DACT_PYTHON=ON" % name)
        return getattr(_core, name)
    if name in ("ACT", "Alg", "Target"):
        from . import legacy  # pylint: disable=import-outside-toplevel
        warnings.warn("act.%s has moved to act.legacy.%s and will be removed from act "
                      "in the next release" % (name, name), DeprecationWarning, stacklevel=2)
        return getattr(legacy, name)
    raise AttributeError("module 'act' has no attribute '%s'" % name)
