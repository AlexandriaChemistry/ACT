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
@brief The exceptions of the Python interface to ACT.
"""
from __future__ import annotations

from typing import List, Optional


class ActError(Exception):
    """!
    @brief Base class of all errors of the Python interface to ACT.

    When a program of ACT fails, the exception carries what is known about
    the run, such that a caller can show it or store it.
    """

    def __init__(self, message: str, *,
                 args: Optional[List[str]] = None,
                 returncode: Optional[int] = None,
                 stdout: str = "",
                 stderr: str = "",
                 log: str = ""):
        """!
        @brief Make the exception.
        @param message    What went wrong, for a person to read.
        @param args       The command line that was run, if a program was run.
        @param returncode The exit status of the program.
        @param stdout     What the program wrote to standard output.
        @param stderr     What the program wrote to standard error.
        @param log        The content of the log file of the program.
        """
        super().__init__(message)
        ## What went wrong
        self.message = message
        ## The command line that was run
        self.command_line = list(args) if args else []
        ## The exit status of the program
        self.returncode = returncode
        ## Standard output of the program
        self.stdout = stdout
        ## Standard error of the program
        self.stderr = stderr
        ## The content of the log file
        self.log = log


class ActOptionError(ActError, ValueError):
    """!
    @brief An option that the program does not have, or a value that the
    option does not accept. Raised before a program is started.
    """


class ActNotFoundError(ActError, FileNotFoundError):
    """!
    @brief The alexandria program, or a file it needs, cannot be found.
    """


class ActFileError(ActError, OSError):
    """!
    @brief A file that the classes of ACT cannot read or write: it is not
    there, it is not what it should be, or the place cannot be written.
    """


class ActInputError(ActError, ValueError):
    """!
    @brief Input that the classes of ACT do not accept, such as a force
    field with a value that ACT does not know.
    """


class ActInternalError(ActError, RuntimeError):
    """!
    @brief Something that should not happen inside ACT, which is a bug to
    report. The process is still alive.
    """


class ActCancelledError(ActError):
    """!
    @brief The program was stopped because the caller asked for it.
    """


class ActTimeoutError(ActError, TimeoutError):
    """!
    @brief The program ran longer than the time it was given and was stopped.
    """
