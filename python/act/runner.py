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
@brief Runs the commands of alexandria as child processes.

Every command is started with a list of arguments, never through a shell,
in a directory of its own. The arguments are checked against the
description of the options that the program gives of itself (see
act.options) before anything is started: an option that does not exist, a
value of the wrong type or an input file that is missing raises an
exception at once.
"""
from __future__ import annotations

import json
import os
import shutil
import signal
import stat
import subprocess
import tempfile
import time
import weakref
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from .errors import (ActCancelledError, ActError, ActNotFoundError, ActOptionError,
                     ActTimeoutError)
from .options import CommandInfo, OptionInfo, ProgramInfo

## A path as a string or a pathlib.Path
PathLike = Union[str, "os.PathLike[str]"]

## Keyword arguments of run() that are not options of the program
RUN_ARGUMENTS = ("workdir", "processes", "timeout", "check", "env", "on_progress", "stop")

## Name of the file that `alexandria help -export json` writes
_OPTIONS_FILE = "alexandria-options.json"

## Options that the program declares as strings, but that are names of input
## files. A relative path is made absolute for these, as for file options,
## because the program runs in another directory than the caller.
_FILE_NAME_STRINGS = ("f",)


class Config:
    """!
    @brief Settings of the Python interface: where the program is and how
    it is started on several processes.

    There is one object of this class, act.config. The settings can also be
    given with environment variables: ACT_ALEXANDRIA for the program and
    ACT_MPIEXEC for the launcher of MPI programs.
    """

    def __init__(self):
        """!
        @brief Make the settings with their defaults.
        """
        ## The alexandria program: a path, a name that is looked up in the
        ## PATH, or a list of words such as a container command followed by
        ## the name of the program. None means: ACT_ALEXANDRIA if it is
        ## set, otherwise the program that was installed together with this
        ## package, otherwise alexandria, or alexandria_d of a double
        ## precision build, from the PATH.
        self.executable: Union[None, PathLike, Sequence[str]] = None
        ## The launcher for several processes, a name or a list of words.
        ## None means: ACT_MPIEXEC if it is set, otherwise mpirun, mpiexec
        ## or srun from the PATH.
        self.mpi_launcher: Union[None, PathLike, Sequence[str]] = None
        ## The option of the launcher that sets the number of processes
        self.mpi_processes_flag: str = "-n"
        ## Seconds after which a program is stopped, None for no limit.
        ## Can be set for a single run with the argument timeout.
        self.timeout: Optional[float] = None

    def executable_words(self) -> List[str]:
        """!
        @brief The command that starts alexandria.
        @return A list of words, of which the first is an existing program.
        @throws ActNotFoundError if the program cannot be found.
        """
        choice = self.executable
        if choice is None:
            choice = (os.environ.get("ACT_ALEXANDRIA") or installed_program()
                      or _program_in_path())
        words = _as_words(choice)
        found = shutil.which(words[0])
        if found is None:
            raise ActNotFoundError(
                "Cannot find the program '%s'. Install ACT, for instance with "
                "'conda install -c conda-forge alexandria', or set "
                "act.config.executable or the environment variable "
                "ACT_ALEXANDRIA to the path of the alexandria program." % words[0])
        return [found] + words[1:]

    def launcher_words(self, processes: int) -> List[str]:
        """!
        @brief The words to put in front of the command to run it on
        several processes.
        @param processes The number of processes.
        @return A list of words, empty for one process.
        @throws ActNotFoundError if no launcher can be found.
        """
        if processes <= 1:
            return []
        choice = self.mpi_launcher
        if choice is None:
            choice = os.environ.get("ACT_MPIEXEC")
        if choice is None:
            for name in ("mpirun", "mpiexec", "srun"):
                if shutil.which(name):
                    choice = name
                    break
        if choice is None:
            raise ActNotFoundError(
                "Cannot find mpirun, mpiexec or srun to start %d processes. Set "
                "act.config.mpi_launcher or the environment variable ACT_MPIEXEC."
                % processes)
        words = _as_words(choice)
        found = shutil.which(words[0])
        if found is None:
            raise ActNotFoundError("Cannot find the launcher '%s'" % words[0])
        return [found] + words[1:] + [self.mpi_processes_flag, str(processes)]


## The names of the alexandria program: a double precision build installs
## alexandria_d, such that it can be next to a single precision build
PROGRAM_NAMES = ("alexandria", "alexandria_d")


def installed_program(package_dir: Optional[PathLike] = None) -> Optional[str]:
    """!
    @brief The alexandria program that was installed together with this
    package, which is the one whose options the package knows best.
    @param package_dir The directory of the package, by default that of
    this file.
    @return The path of the program. None if the package is used from the
    source code, was installed on its own, or if the program is gone.

    The installation of ACT writes the module _installed with the path of
    the program relative to the package, such that the installation can be
    moved as a whole.
    """
    try:
        from . import _installed  # pylint: disable=import-outside-toplevel
    except ImportError:
        return None
    relative = getattr(_installed, "EXECUTABLE", None)
    if not isinstance(relative, str) or not relative:
        return None
    if package_dir is None:
        package_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.normpath(os.path.join(os.fspath(package_dir), relative))
    if os.path.isfile(path) and os.access(path, os.X_OK):
        return path
    return None


def _program_in_path() -> str:
    """!
    @brief The name under which the alexandria program is in the PATH.
    @return The first of PROGRAM_NAMES that is found, alexandria if none is.
    """
    for name in PROGRAM_NAMES:
        if shutil.which(name):
            return name
    return PROGRAM_NAMES[0]


def _as_words(choice: Union[PathLike, Sequence[str]]) -> List[str]:
    """!
    @brief Turn a setting into a list of words.
    @param choice A path or a list of words.
    @return The list of words. A single path is one word, also when it has
    spaces in it: nothing is split and nothing is given to a shell.
    """
    if isinstance(choice, (str, os.PathLike)):
        return [os.fspath(choice)]
    words = [os.fspath(word) for word in choice]
    if not words:
        raise ActError("The command to start a program is empty")
    return words


## The settings of the Python interface
config = Config()

## The descriptions of the options that were read, by command line of the
## program, each with the state of the program files it was read from
_program_cache: Dict[Tuple[str, ...], Tuple[Tuple[Tuple[int, int], ...], ProgramInfo]] = {}


def _program_state(words: Sequence[str]) -> Tuple[Tuple[int, int], ...]:
    """!
    @brief What tells that a program was replaced by another one at the
    same path, as happens when ACT is built or installed again.
    @param words The command that starts the program.
    @return The time of the last change and the size of each word that is
    a file: the program, or an interpreter and its script.
    """
    state = []
    for word in words:
        try:
            status = os.stat(word)
        except OSError:
            continue
        if stat.S_ISREG(status.st_mode):
            state.append((status.st_mtime_ns, status.st_size))
    return tuple(state)


def program_info(refresh: bool = False) -> ProgramInfo:
    """!
    @brief The commands of alexandria and their options.

    The program is asked for the description the first time; after that the
    answer is kept for as long as Python runs and the program file stays
    the same. A program that is replaced by a newer one at the same path is
    asked again, such that a process that runs for a long time, like the
    server of a user interface, does not go on with the options of the
    program it started with.
    @param refresh Ask the program again in any case.
    @return The description.
    @throws ActNotFoundError if the program cannot be found.
    @throws ActError if the program does not give a description, which
    means that it is too old for this interface.
    """
    words = config.executable_words()
    key = tuple(words)
    state = _program_state(words)
    if not refresh and key in _program_cache and _program_cache[key][0] == state:
        return _program_cache[key][1]
    with tempfile.TemporaryDirectory(prefix="act-options-") as tmp:
        args = words + ["help", "-export", "json"]
        try:
            done = subprocess.run(args, cwd=tmp, stdin=subprocess.DEVNULL,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  universal_newlines=True, timeout=120, check=False)
        except (OSError, subprocess.SubprocessError) as exc:
            raise ActError("Cannot run %s: %s" % (words[0], exc), args=args) from exc
        path = os.path.join(tmp, _OPTIONS_FILE)
        if done.returncode != 0 or not os.path.exists(path):
            raise ActError(
                "%s does not describe its options. It is probably older than this "
                "Python interface, which needs 'alexandria help -export json'."
                % words[0], args=args, returncode=done.returncode,
                stdout=done.stdout, stderr=done.stderr)
        info = ProgramInfo.from_file(path)
    _program_cache[key] = (state, info)
    return info


def _is_number(text: str) -> bool:
    """!
    @brief Whether a text is a number.
    @param text The text.
    @return True if float() accepts it.
    """
    try:
        float(text)
    except ValueError:
        return False
    return True


def _check_text(option: OptionInfo, command: str, text: str) -> str:
    """!
    @brief Check a text that goes on the command line as the value of an option.

    A value that starts with a dash would be read by the program as the
    name of another option, such that a value could switch on something the
    caller did not ask for. Negative numbers are the exception.
    @param option  The option the value is for.
    @param command The name of the command, for the message.
    @param text    The value.
    @return The value.
    @throws ActOptionError if the value cannot be passed safely.
    """
    if "\0" in text:
        raise ActOptionError("The value of option '%s' of %s contains a null character"
                             % (option.name, command))
    if text.startswith("-") and not _is_number(text):
        raise ActOptionError("The value '%s' of option '%s' of %s starts with a dash, "
                             "which the program would read as another option"
                             % (text, option.name, command))
    return text


def _one_value(option: OptionInfo, command: str, value: Any, base: Path) -> str:
    """!
    @brief Convert one value of an option that is not a boolean or a file.
    @param option  The option.
    @param command The name of the command, for messages.
    @param value   The value given by the caller.
    @param base    The directory that relative paths are relative to.
    @return The text for the command line.
    @throws ActOptionError if the value has the wrong type or is not allowed.
    """
    where = "option '%s' of %s" % (option.name, command)
    if option.type == "int":
        # A bool is an int in Python, but True as a number of steps is a mistake
        if isinstance(value, bool) or not isinstance(value, int):
            raise ActOptionError("%s takes an integer, not %r" % (where.capitalize(), value))
        return str(value)
    if option.type == "real":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ActOptionError("%s takes a number, not %r" % (where.capitalize(), value))
        if value != value or value in (float("inf"), float("-inf")):
            raise ActOptionError("%s takes a finite number, not %r" % (where.capitalize(), value))
        return repr(float(value))
    if option.type == "enum":
        if not isinstance(value, str) or value not in option.values:
            raise ActOptionError("%s takes one of %s, not %r"
                                 % (where.capitalize(), ", ".join(option.values), value))
        return value
    # A string. A path object is made absolute, since the program runs in
    # another directory than the caller.
    if isinstance(value, os.PathLike):
        return _check_text(option, command, str(_absolute(value, base)))
    if not isinstance(value, str):
        raise ActOptionError("%s takes a string, not %r" % (where.capitalize(), value))
    if option.name in _FILE_NAME_STRINGS and value and (base / value).is_file():
        return _check_text(option, command, str(_absolute(value, base)))
    return _check_text(option, command, value)


def _absolute(path: PathLike, base: Path) -> Path:
    """!
    @brief Make a path absolute.
    @param path The path.
    @param base The directory a relative path is relative to.
    @return The absolute path, without resolving symbolic links.
    """
    path = Path(os.fspath(path)).expanduser()
    if not path.is_absolute():
        path = base / path
    return Path(os.path.normpath(str(path)))


def _as_list(option: OptionInfo, command: str, value: Any) -> List[Any]:
    """!
    @brief The values given for an option as a list.
    @param option  The option.
    @param command The name of the command, for messages.
    @param value   One value, or a list or tuple of values.
    @return The list of values.
    @throws ActOptionError if several values are given to an option that
    takes one, or more than it takes.
    """
    if isinstance(value, (list, tuple)):
        values = list(value)
        if not option.multiple and len(values) != 1:
            raise ActOptionError("Option '%s' of %s takes one value, not %d"
                                 % (option.name, command, len(values)))
        if option.max_values is not None and len(values) > option.max_values:
            raise ActOptionError("Option '%s' of %s takes at most %d values, not %d"
                                 % (option.name, command, option.max_values, len(values)))
        if not values:
            raise ActOptionError("Option '%s' of %s got an empty list" % (option.name, command))
        return values
    return [value]


def _default_names(option: OptionInfo) -> List[str]:
    """!
    @brief The default file names of a file option.
    @param option The option.
    @return The names, empty if the option has no default.
    """
    if option.default is None:
        return []
    if isinstance(option.default, list):
        return [str(name) for name in option.default]
    return [str(option.default)]


def build_arguments(command: CommandInfo,
                    options: Mapping[str, Any],
                    workdir: PathLike,
                    base: Optional[PathLike] = None) -> Tuple[List[str], Dict[str, List[Path]]]:
    """!
    @brief Turn keyword arguments into the command line of a command.

    The rules are:
    - A value None means that the option is not given.
    - A boolean option takes True or False and becomes -name or -noname.
    - Numbers, enums and strings are checked against the type of the option.
    - An input file must exist. Relative paths are relative to @p base and
      are made absolute, since the program runs in @p workdir.
    - An output file can be a path, or True for the default name in
      @p workdir, or False to leave an optional file out.
    - An option that takes several values takes a list.
    @param command The command.
    @param options The options by their name on the command line or their
    name as a Python keyword argument.
    @param workdir The directory the program will run in.
    @param base    The directory that relative paths are relative to, the
    current directory if None.
    @return The arguments after the name of the command, and for each
    output file option the files the program is expected to write. That
    includes the files with default names that the program always writes.
    @throws ActOptionError for an unknown option or a value it cannot take.
    @throws ActNotFoundError for an input file that does not exist.
    """
    workdir = Path(os.fspath(workdir))
    base = Path(os.fspath(base)) if base is not None else Path.cwd()
    args: List[str] = []
    outputs: Dict[str, List[Path]] = {}
    given = set()
    for key, value in options.items():
        option = command.find_option(key)
        if option.name in given:
            raise ActOptionError("Option '%s' of %s is given twice" % (option.name, command.name))
        given.add(option.name)
        if value is None:
            continue
        flag = "-" + option.name
        if option.type == "bool":
            if not isinstance(value, bool):
                raise ActOptionError("Option '%s' of %s takes True or False, not %r"
                                     % (option.name, command.name, value))
            args.append(flag if value else "-no" + option.name)
        elif option.is_file:
            paths = _file_values(option, command.name, value, workdir, base)
            if paths is None:
                continue
            args.append(flag)
            args.extend(str(path) for path in paths)
            if option.is_output:
                outputs[option.name] = paths
        else:
            values = _as_list(option, command.name, value)
            args.append(flag)
            args.extend(_one_value(option, command.name, v, base) for v in values)
    # Output files that the program writes without being asked
    for option in command.options.values():
        if option.is_file and option.is_output and option.required and option.name not in outputs:
            names = _default_names(option)
            if names:
                outputs[option.name] = [workdir / name for name in names]
    return args, outputs


def _file_values(option: OptionInfo, command: str, value: Any,
                 workdir: Path, base: Path) -> Optional[List[Path]]:
    """!
    @brief Convert the value of a file option to paths.
    @param option  The option.
    @param command The name of the command, for messages.
    @param value   What the caller gave: a path, a list of paths, or for an
    output file True or False.
    @param workdir The directory the program will run in.
    @param base    The directory that relative paths are relative to.
    @return The absolute paths, or None if the option is to be left out.
    @throws ActOptionError if the value cannot be used.
    @throws ActNotFoundError if an input file does not exist.
    """
    if isinstance(value, bool):
        if not option.is_output or option.is_input:
            raise ActOptionError("Option '%s' of %s takes the name of a file, not %r"
                                 % (option.name, command, value))
        if not value:
            if option.required:
                raise ActOptionError("The file of option '%s' of %s is always written, "
                                     "it cannot be switched off" % (option.name, command))
            return None
        names = _default_names(option)
        if not names:
            raise ActOptionError("Option '%s' of %s has no default file name, give a name"
                                 % (option.name, command))
        return [workdir / name for name in names]
    paths = []
    for item in _as_list(option, command, value):
        if not isinstance(item, (str, os.PathLike)):
            raise ActOptionError("Option '%s' of %s takes the name of a file, not %r"
                                 % (option.name, command, item))
        if not os.fspath(item):
            raise ActOptionError("Option '%s' of %s got an empty file name"
                                 % (option.name, command))
        path = _absolute(item, base)
        if "\0" in str(path):
            raise ActOptionError("The file name of option '%s' of %s contains a null character"
                                 % (option.name, command))
        if option.is_input and not option.is_output:
            if not path.exists():
                # A library file is looked up by the program in its own
                # data directory when it is given by name only.
                if option.library and os.sep not in os.fspath(item):
                    paths.append(Path(os.fspath(item)))
                    continue
                if not option.allow_missing:
                    raise ActNotFoundError("The file '%s' for option '%s' of %s does not exist"
                                           % (path, option.name, command))
        elif option.is_output and not path.parent.is_dir():
            raise ActNotFoundError("The directory '%s' for the output file of option '%s' "
                                   "of %s does not exist" % (path.parent, option.name, command))
        paths.append(path)
    return paths


class Result:
    """!
    @brief What a command of alexandria produced.

    If the command ran in a temporary directory, that directory is removed
    when the result is no longer used, or when close() is called. Copy the
    files that should be kept, or give a directory with the argument
    workdir of the command.
    """

    def __init__(self, command: str, args: List[str], returncode: int,
                 stdout: str, stderr: str, workdir: Path,
                 files: Dict[str, List[Path]], temporary: bool):
        """!
        @brief Store the outcome of a run.
        @param command    The name of the command.
        @param args       The complete command line.
        @param returncode The exit status.
        @param stdout     Standard output of the program.
        @param stderr     Standard error of the program.
        @param workdir    The directory the program ran in.
        @param files      The output files that exist, by option name.
        @param temporary  Whether the directory is removed with the result.
        """
        ## The name of the command, such as minimize
        self.command = command
        ## The complete command line that was run
        self.args = args
        ## The exit status of the program
        self.returncode = returncode
        ## Standard output of the program
        self.stdout = stdout
        ## Standard error of the program
        self.stderr = stderr
        ## The directory the program ran in
        self.workdir = workdir
        ## The output files that exist, by option name. Each is a list,
        ## since some options take several files.
        self.files = files
        self._json: Any = None
        self._json_read = False
        self._finalizer = None
        if temporary:
            self._finalizer = weakref.finalize(self, shutil.rmtree, str(workdir), True)

    @property
    def ok(self) -> bool:
        """!
        @brief Whether the program finished without error.
        @return True if the exit status is zero.
        """
        return self.returncode == 0

    def file(self, option: str) -> Optional[Path]:
        """!
        @brief The output file of an option.
        @param option The name of the option, such as c for the structure.
        @return The first file of the option, None if it was not written.
        """
        paths = self.files.get(option)
        return paths[0] if paths else None

    @property
    def json(self) -> Any:
        """!
        @brief The results the program wrote in JSON format.
        @return The parsed file of the option -json, None if the command
        has no such option or the file was not written.
        @throws ActError if the file is not valid JSON.
        """
        if not self._json_read:
            self._json_read = True
            path = self.file("json")
            if path is not None:
                try:
                    with open(path, "r", encoding="utf-8") as handle:
                        self._json = json.load(handle)
                except (OSError, ValueError) as exc:
                    raise ActError("Cannot read the results of alexandria %s from %s: %s"
                                   % (self.command, path, exc)) from exc
        return self._json

    @property
    def data(self) -> Any:
        """!
        @brief The results without the outer key.
        @return The JSON file of a command has one key, the name of the
        command, with everything below it. This is the content below that
        key, or None if there is no JSON file.
        """
        doc = self.json
        if isinstance(doc, dict) and len(doc) == 1:
            return next(iter(doc.values()))
        return doc

    @property
    def log(self) -> str:
        """!
        @brief The log file of the program.
        @return The text of the file of option -g, empty if there is none.
        """
        path = self.file("g")
        if path is None:
            return ""
        try:
            return path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return ""

    def close(self) -> None:
        """!
        @brief Remove the temporary directory of the run, if it has one.
        """
        if self._finalizer is not None:
            self._finalizer()

    def __enter__(self) -> "Result":
        """!
        @brief Use the result in a with statement.
        @return The result itself.
        """
        return self

    def __exit__(self, *exc: Any) -> None:
        """!
        @brief Remove the temporary directory at the end of a with statement.
        @param exc The exception that ended the block, if any.
        """
        self.close()

    def __repr__(self) -> str:
        """!
        @brief A short description for the Python prompt.
        @return Command, exit status and the options that have files.
        """
        return "<act.Result %s returncode=%d files=%s>" % (
            self.command, self.returncode, sorted(self.files))


def _error_message(command: str, returncode: int, stdout: str, stderr: str) -> str:
    """!
    @brief Extract the reason of a failure from the output of the program.
    @param command    The name of the command.
    @param returncode The exit status.
    @param stdout     Standard output of the program.
    @param stderr     Standard error of the program.
    @return A message for a person. The complete output is in the exception.
    """
    lines = stderr.splitlines()
    reason: List[str] = []
    # The block that the program prints for a fatal error: the reason is
    # between the line that names the kind of error and the footer.
    for start, line in enumerate(lines):
        text = line.strip()
        if text.endswith(":") and ("error" in text.lower()) and not text.startswith("Source"):
            for follow in lines[start:]:
                if follow.startswith("For more information") or follow.startswith("-----"):
                    break
                if follow.strip():
                    reason.append(follow.strip())
            break
    if not reason:
        # Otherwise the last lines that say something, without the notice
        # that MPI prints when a program stops
        noise = ("MPI_ABORT", "NOTE: invoking MPI_ABORT", "You may or may not",
                 "exactly when Open MPI", "with errorcode", "-----")
        for source in (stderr, stdout):
            useful = [l.strip() for l in source.splitlines()
                      if l.strip() and not l.strip().startswith(noise)]
            errors = [l for l in useful if "error" in l.lower() or "could not" in l.lower()
                      or "failed" in l.lower()]
            reason = errors[-3:] if errors else useful[-3:]
            if reason:
                break
    detail = " ".join(reason) if reason else "no message"
    return "alexandria %s failed with exit status %d: %s" % (command, returncode, detail)


class _ProgressReader:
    """!
    @brief Reads the lines that a program adds to its progress file.

    Each line of the file is a JSON object. A line is passed on when it is
    complete, such that a line that the program is still writing is not
    read half.
    """

    def __init__(self, path: Path):
        """!
        @brief Start following a file, which need not exist yet.
        @param path The progress file.
        """
        self._path = path
        self._offset = 0
        self._rest = b""

    def new_objects(self) -> List[Any]:
        """!
        @brief The objects on the lines that were added since the last call.
        @return The parsed lines, in the order of the file. Lines that are
        not valid JSON are skipped.
        """
        try:
            with open(self._path, "rb") as handle:
                handle.seek(self._offset)
                data = handle.read()
        except OSError:
            return []
        self._offset += len(data)
        lines = (self._rest + data).split(b"\n")
        self._rest = lines.pop()
        objects = []
        for line in lines:
            if not line.strip():
                continue
            try:
                objects.append(json.loads(line.decode("utf-8", errors="replace")))
            except ValueError:
                continue
        return objects


def _terminate(process: "subprocess.Popen[str]") -> None:
    """!
    @brief Stop a program and the processes it started.
    @param process The program, started in a session of its own.
    """
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except (OSError, AttributeError):
        process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except (OSError, AttributeError):
            process.kill()
        process.wait()


def run(command: str, *,
        workdir: Optional[PathLike] = None,
        processes: int = 1,
        timeout: Optional[float] = None,
        check: bool = True,
        env: Optional[Mapping[str, str]] = None,
        on_progress: Optional[Callable[[Any], None]] = None,
        stop: Any = None,
        **options: Any) -> Result:
    """!
    @brief Run a command of alexandria and return what it produced.

    For example
    @code{.py}
    import act
    res = act.run("minimize", ff="ACS-pg1.xml", f="acetone.sdf")
    print(res.data["rmsd"]["value"], res.file("c"))
    @endcode
    The functions in act.commands do the same with the name of the command
    as the name of the function: act.commands.minimize(ff=..., f=...).

    If the command can write its results in JSON format (option -json), it
    is asked to, unless json=False is given.
    @param command   The name of the command, such as minimize.
    @param workdir   The directory to run in, which is made if it does not
    exist and is kept. Without it the command runs in a temporary directory
    that is removed with the result.
    @param processes The number of processes; more than one starts the
    program with the MPI launcher.
    @param timeout   Seconds after which the program is stopped, None for
    act.config.timeout.
    @param check     Raise ActError if the program fails. With False the
    result is returned and its returncode tells what happened.
    @param env       Environment variables to set for the program, on top of
    those of Python.
    @param on_progress A function that is called while the program runs,
    for commands with the option -progress, such as simulate. It gets one
    argument, the object that the program wrote: for simulate a dict with
    the step, the time, the temperature and the energies. An exception in
    the function stops the program.
    @param stop      An object with a method is_set(), such as a
    threading.Event, with which another thread can stop the program: when
    is_set() returns True the program is stopped and ActCancelledError is
    raised.
    @param options   The options of the command, see build_arguments() for
    the rules and act.program_info() for the options of each command.
    @return The result.
    @throws ActOptionError   for an unknown command or option, or a wrong value.
    @throws ActNotFoundError if the program or an input file does not exist.
    @throws ActTimeoutError  if the program ran too long.
    @throws ActCancelledError if the program was stopped through @p stop.
    @throws ActError         if the program failed and check is True.
    """
    info = program_info().find_command(command)
    if isinstance(processes, bool) or not isinstance(processes, int) or processes < 1:
        raise ActOptionError("The number of processes must be a positive integer, not %r"
                             % (processes,))
    # Results in JSON format unless the caller decides otherwise
    options = dict(options)
    if info.has_option("json") and "json" not in options:
        options["json"] = True
    if on_progress is not None:
        if not info.has_option("progress"):
            raise ActOptionError("alexandria %s does not report its progress" % command)
        if not callable(on_progress):
            raise ActOptionError("on_progress must be a function, not %r" % (on_progress,))
        if options.get("progress") in (None, False):
            options["progress"] = True

    temporary = workdir is None
    if temporary:
        rundir = Path(tempfile.mkdtemp(prefix="act-%s-" % command))
    else:
        rundir = _absolute(workdir, Path.cwd())
        rundir.mkdir(parents=True, exist_ok=True)
    try:
        args, outputs = build_arguments(info, options, rundir)
        full = config.launcher_words(processes) + config.executable_words() + [command] + args
        environment = dict(os.environ)
        if env:
            environment.update({str(k): str(v) for k, v in env.items()})
        limit = timeout if timeout is not None else config.timeout
        # The output goes to files, not to pipes: a pipe that nobody reads
        # while the progress is followed would stop the program when it is full.
        with tempfile.TemporaryFile(mode="w+", encoding="utf-8", errors="replace") as out, \
                tempfile.TemporaryFile(mode="w+", encoding="utf-8", errors="replace") as err:
            try:
                process = subprocess.Popen(full, cwd=str(rundir), env=environment,
                                           stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                                           start_new_session=True)
            except (OSError, ValueError) as exc:
                raise ActError("Cannot start %s: %s" % (full[0], exc), args=full) from exc

            def output() -> Tuple[str, str]:
                """!
                @brief What the program wrote so far.
                @return Standard output and standard error.
                """
                out.seek(0)
                err.seek(0)
                return out.read(), err.read()

            reader = None
            if on_progress is not None:
                reader = _ProgressReader(outputs["progress"][0])
            deadline = None if limit is None else time.monotonic() + limit
            try:
                while True:
                    try:
                        # A short wait when there is progress to pass on
                        # or a request to stop to look for
                        wait = 0.05 if (reader is not None or stop is not None) else None
                        if deadline is not None:
                            left = max(0.0, deadline - time.monotonic())
                            wait = left if wait is None else min(wait, left)
                        process.wait(timeout=wait)
                        finished = True
                    except subprocess.TimeoutExpired:
                        finished = False
                    if reader is not None and on_progress is not None:
                        for item in reader.new_objects():
                            on_progress(item)
                    if finished:
                        break
                    if stop is not None and stop.is_set():
                        _terminate(process)
                        stdout, stderr = output()
                        raise ActCancelledError("alexandria %s was stopped" % command,
                                                args=full, returncode=process.returncode,
                                                stdout=stdout, stderr=stderr)
                    if deadline is not None and time.monotonic() >= deadline:
                        _terminate(process)
                        stdout, stderr = output()
                        raise ActTimeoutError("alexandria %s was stopped after %g seconds"
                                              % (command, limit), args=full,
                                              returncode=process.returncode,
                                              stdout=stdout, stderr=stderr)
            except BaseException:
                # For instance Ctrl-C or an error in on_progress: do not
                # leave the program running
                if process.poll() is None:
                    _terminate(process)
                raise
            stdout, stderr = output()
        files = {name: [p for p in paths if p.exists()] for name, paths in outputs.items()}
        files = {name: paths for name, paths in files.items() if paths}
        result = Result(command, full, process.returncode, stdout, stderr,
                        rundir, files, temporary)
    except BaseException:
        if temporary:
            shutil.rmtree(str(rundir), ignore_errors=True)
        raise
    if check and result.returncode != 0:
        error = ActError(_error_message(command, result.returncode, stdout, stderr),
                         args=full, returncode=result.returncode,
                         stdout=stdout, stderr=stderr, log=result.log)
        result.close()
        raise error
    return result
