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
@brief The class ACT of the earlier Python module act.py of ACT-utilities.

That module had the name that this package has now. Scripts that did
@code{.py}
from act import ACT
@endcode
keep working for one release, with a warning, and should change to
@code{.py}
from act.legacy import ACT
@endcode
or to the functions in act.commands. The class itself is still part of
ACT-utilities, where its module is called act_legacy now, since a module
act next to this package would hide one of the two.
"""
from __future__ import annotations

from types import ModuleType


def _load_old_module() -> ModuleType:
    """!
    @brief Load the module of ACT-utilities with the class ACT.
    @return The module.
    @throws ImportError if ACT-utilities is not installed, or a version
    of it from before the module got the name act_legacy.
    """
    try:
        import act_legacy  # pylint: disable=import-outside-toplevel
    except ImportError as error:
        raise ImportError(
            "The class ACT is part of ACT-utilities, in the module act_legacy, which "
            "cannot be imported (%s). Install or update ACT-utilities, for instance with "
            "'conda install -c conda-forge act-utilities'." % error) from error
    return act_legacy


## The module act_legacy of ACT-utilities
_old = _load_old_module()

## The class to run training with alexandria train_ff
ACT = _old.ACT
## The optimization algorithms of the class ACT
Alg = _old.Alg
## The training targets of the class ACT
Target = _old.Target
