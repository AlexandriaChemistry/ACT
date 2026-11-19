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
ACT-utilities; this module finds it there.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from types import ModuleType
from typing import Optional


def _find_old_module() -> Optional[str]:
    """!
    @brief Find the file act.py of ACT-utilities.
    @return The path of the file, None if ACT-utilities is not installed.
    ACT-utilities installs its modules in a directory that it adds to the
    path of Python; this package is installed in the same directory.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [os.path.join(here, "act.py")]
    for entry in sys.path:
        if entry and os.path.abspath(entry) != os.path.dirname(here):
            candidates.append(os.path.join(entry, "act.py"))
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    return None


def _load_old_module() -> ModuleType:
    """!
    @brief Load act.py of ACT-utilities under another name.
    @return The module.
    @throws ImportError if ACT-utilities is not installed.
    """
    path = _find_old_module()
    if path is None:
        raise ImportError("The class ACT is part of ACT-utilities, which is not installed. "
                          "Install it with 'conda install -c conda-forge act-utilities'.")
    spec = importlib.util.spec_from_file_location("act._legacy_act", path)
    if spec is None or spec.loader is None:
        raise ImportError("Cannot load %s" % path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


## The module act.py of ACT-utilities
_old = _load_old_module()

## The class to run training with alexandria train_ff
ACT = _old.ACT
## The optimization algorithms of the class ACT
Alg = _old.Alg
## The training targets of the class ACT
Target = _old.Target
