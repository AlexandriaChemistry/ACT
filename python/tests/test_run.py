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
@brief Tests of running a program, with a stand-in for alexandria, such
that these tests need no build of ACT.
"""
import inspect
import os
import sys
import tempfile
import threading
import time
import unittest
import warnings
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import act  # noqa: E402  pylint: disable=wrong-import-position

## The stand-in program
FAKE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fake_alexandria.py")


class RunTest(unittest.TestCase):
    """!
    @brief Running commands and reading their results.
    """

    def setUp(self):
        """!
        @brief Use the stand-in program and make an input file.
        """
        self._executable = act.config.executable
        self._timeout = act.config.timeout
        act.config.executable = [sys.executable, FAKE]
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(os.path.realpath(self._tmp.name))
        self.ff = self.dir / "ff.xml"
        self.ff.write_text("x")

    def tearDown(self):
        """!
        @brief Restore the settings and remove the files.
        """
        act.config.executable = self._executable
        act.config.timeout = self._timeout
        self._tmp.cleanup()

    def test_description_is_read_from_the_program(self):
        """!
        @brief The program is asked for its commands.
        """
        self.assertEqual(["work"], act.commands.available())
        self.assertIn("maxiter", act.program_info().find_command("work").options)

    def test_run_gives_results(self):
        """!
        @brief A run returns the JSON results, the files and the output.
        """
        res = act.run("work", ff=self.ff, maxiter=5, oneH=True)
        self.assertTrue(res.ok)
        self.assertEqual(0, res.returncode)
        self.assertEqual("work", res.command)
        self.assertIn("done", res.stdout)
        # Results in JSON format without asking
        self.assertEqual(0.0192, res.data["rmsd"]["value"])
        self.assertEqual(res.json["work"], res.data)
        argv = res.data["argv"]
        self.assertEqual(["-maxiter", "5"], argv[argv.index("-maxiter"):argv.index("-maxiter")+2])
        self.assertIn("-oneH", argv)
        # The files that the program always writes are found
        self.assertEqual(res.workdir / "confout.pdb", res.file("c"))
        self.assertTrue(res.file("c").exists())
        self.assertEqual("log of work\n", res.log)
        self.assertIsNone(res.file("extra"))
        self.assertIn("work", repr(res))

    def test_program_runs_in_its_own_directory(self):
        """!
        @brief The program runs in a temporary directory that is removed.
        """
        res = act.run("work", ff=self.ff)
        workdir = res.workdir
        self.assertEqual(os.path.realpath(res.data["cwd"]), os.path.realpath(str(workdir)))
        self.assertNotEqual(os.path.realpath(os.getcwd()), os.path.realpath(str(workdir)))
        self.assertTrue(workdir.is_dir())
        res.close()
        self.assertFalse(workdir.exists())
        # Also in a with statement
        with act.run("work", ff=self.ff) as res2:
            workdir2 = res2.workdir
            self.assertTrue(workdir2.is_dir())
        self.assertFalse(workdir2.exists())

    def test_given_directory_is_kept(self):
        """!
        @brief A directory given by the caller is made and kept.
        """
        workdir = self.dir / "new" / "run1"
        res = act.run("work", ff=self.ff, workdir=workdir, extra=True)
        self.assertEqual(workdir, res.workdir)
        self.assertEqual(workdir / "extra.xvg", res.file("extra"))
        res.close()
        self.assertTrue((workdir / "confout.pdb").exists())
        self.assertTrue((workdir / "work.json").exists())

    def test_output_file_where_the_caller_wants_it(self):
        """!
        @brief An output file with a name is written there.
        """
        target = self.dir / "min.pdb"
        res = act.run("work", ff=self.ff, c=target)
        self.assertEqual(target, res.file("c"))
        res.close()
        self.assertTrue(target.exists())

    def test_json_can_be_switched_off(self):
        """!
        @brief Without the JSON file there are no data.
        """
        res = act.run("work", ff=self.ff, json=False)
        self.assertIsNone(res.json)
        self.assertIsNone(res.data)
        self.assertIsNone(res.file("json"))

    def test_environment(self):
        """!
        @brief Environment variables can be set for the program.
        """
        res = act.run("work", ff=self.ff, env={"ACT_TEST_VARIABLE": "hello"})
        self.assertEqual("hello", res.data["env"])

    def test_failure_raises_with_the_reason(self):
        """!
        @brief A program that fails raises an error with its message and output.
        """
        with self.assertRaises(act.ActError) as caught:
            act.run("work", ff=self.ff, fail=3)
        error = caught.exception
        self.assertEqual(3, error.returncode)
        self.assertIn("exit status 3", str(error))
        self.assertIn("The calculation went wrong.", str(error))
        self.assertNotIn("For more information", str(error))
        self.assertIn("Fatal error:", error.stderr)
        self.assertEqual("log of work\n", error.log)
        self.assertEqual("work", error.command_line[2])

    def test_failure_without_check(self):
        """!
        @brief Without check the result tells that the program failed.
        """
        res = act.run("work", ff=self.ff, fail=2, check=False)
        self.assertFalse(res.ok)
        self.assertEqual(2, res.returncode)
        self.assertIsNone(res.data)
        self.assertIsNone(res.file("c"))

    def test_nothing_runs_with_a_wrong_option(self):
        """!
        @brief Wrong arguments are refused before the program starts.
        """
        before = set(os.listdir(tempfile.gettempdir()))
        with self.assertRaises(act.ActOptionError):
            act.run("work", ff=self.ff, maxiter="many")
        with self.assertRaises(act.ActNotFoundError):
            act.run("work", ff=self.dir / "nothere.xml")
        with self.assertRaises(act.ActOptionError):
            act.run("nosuchcommand")
        with self.assertRaises(act.ActOptionError):
            act.run("work", ff=self.ff, processes=0)
        # No temporary directories are left behind
        left = [n for n in set(os.listdir(tempfile.gettempdir())) - before
                if n.startswith("act-work-")]
        self.assertEqual([], left)

    def test_timeout_stops_the_program(self):
        """!
        @brief A program that takes too long is stopped.
        """
        start = time.time()
        with self.assertRaises(act.ActTimeoutError) as caught:
            act.run("work", ff=self.ff, sleep=30.0, timeout=0.5)
        self.assertLess(time.time() - start, 15)
        self.assertIn("0.5 seconds", str(caught.exception))
        # The default for all runs
        act.config.timeout = 0.5
        with self.assertRaises(act.ActTimeoutError):
            act.run("work", ff=self.ff, sleep=30.0)

    def test_progress_is_passed_on_while_the_program_runs(self):
        """!
        @brief The function for the progress gets every line, the first
        well before the program is done.
        """
        seen = []
        start = time.time()

        def progress(item):
            seen.append((time.time() - start, item))

        res = act.run("work", ff=self.ff, on_progress=progress)
        total = time.time() - start
        self.assertEqual([0, 1, 2, 3, 4], [item["step"] for _, item in seen])
        self.assertEqual(5, seen[0][1]["nsteps"])
        # The program takes a second, the first line is there at once
        self.assertLess(seen[0][0], total - 0.5)
        self.assertTrue(res.file("progress").exists())

    def test_error_in_progress_function_stops_the_program(self):
        """!
        @brief An exception in the function ends the run at once.
        """
        def progress(item):
            raise RuntimeError("stop at step %d" % item["step"])

        start = time.time()
        with self.assertRaises(RuntimeError):
            act.run("work", ff=self.ff, on_progress=progress, sleep=30.0)
        self.assertLess(time.time() - start, 10)

    def test_progress_needs_a_function(self):
        """!
        @brief Something that cannot be called is refused.
        """
        with self.assertRaises(act.ActOptionError):
            act.run("work", ff=self.ff, on_progress="print")

    def test_another_thread_can_stop_the_program(self):
        """!
        @brief A program is stopped when the event is set.
        """
        stop = threading.Event()
        timer = threading.Timer(0.3, stop.set)
        timer.start()
        start = time.time()
        try:
            with self.assertRaises(act.ActCancelledError) as caught:
                act.run("work", ff=self.ff, sleep=30.0, stop=stop)
        finally:
            timer.cancel()
        self.assertLess(time.time() - start, 10)
        self.assertIn("stopped", str(caught.exception))
        # An event that is not set changes nothing
        res = act.run("work", ff=self.ff, stop=threading.Event())
        self.assertTrue(res.ok)

    def test_program_that_does_not_exist(self):
        """!
        @brief A missing program gives an error that says what to do.
        """
        act.config.executable = "no-such-alexandria-program"
        with self.assertRaises(act.ActNotFoundError) as caught:
            act.run("work", ff=self.ff)
        self.assertIn("ACT_ALEXANDRIA", str(caught.exception))

    def test_program_that_gives_no_description(self):
        """!
        @brief A program without the JSON export is reported as too old.
        """
        act.config.executable = [sys.executable, "-c", "print('old')"]
        with self.assertRaises(act.ActError) as caught:
            act.program_info()
        self.assertIn("older", str(caught.exception))

    def test_launcher_for_several_processes(self):
        """!
        @brief Several processes put the launcher in front of the command.
        """
        launcher = act.config.mpi_launcher
        try:
            # The stand-in for mpirun runs what follows the number of processes
            act.config.mpi_launcher = [sys.executable, "-c",
                                       "import subprocess, sys; "
                                       "sys.exit(subprocess.call(sys.argv[3:]))"]
            res = act.run("work", ff=self.ff, processes=4)
            self.assertEqual(["-n", "4"], res.args[3:5])
            self.assertTrue(res.ok)
            self.assertEqual(0.0192, res.data["rmsd"]["value"])
        finally:
            act.config.mpi_launcher = launcher


class CommandsTest(unittest.TestCase):
    """!
    @brief The functions in act.commands.
    """

    def setUp(self):
        """!
        @brief Use the stand-in program.
        """
        self._executable = act.config.executable
        act.config.executable = [sys.executable, FAKE]
        act.commands._functions.clear()  # pylint: disable=protected-access
        self._tmp = tempfile.TemporaryDirectory()
        self.ff = Path(os.path.realpath(self._tmp.name)) / "ff.xml"
        self.ff.write_text("x")

    def tearDown(self):
        """!
        @brief Restore the settings.
        """
        act.config.executable = self._executable
        act.commands._functions.clear()  # pylint: disable=protected-access
        self._tmp.cleanup()

    def test_function_runs_the_command(self):
        """!
        @brief A function runs its command.
        """
        res = act.commands.work(ff=self.ff, minalg="Steep")
        self.assertIn("Steep", res.data["argv"])

    def test_function_is_documented(self):
        """!
        @brief Name, help and signature come from the program.
        """
        work = act.commands.work
        self.assertEqual("work", work.__name__)
        self.assertIn("Pretend to do a calculation.", work.__doc__)
        self.assertIn("maxiter (int, default 100)", work.__doc__)
        self.assertIn("one of LBFGS, Newton, Steep", work.__doc__)
        self.assertIn("Options for experts", work.__doc__)
        parameters = inspect.signature(work).parameters
        for name in ("ff", "maxiter", "del_", "workdir", "processes", "timeout"):
            self.assertIn(name, parameters)

    def test_unknown_command(self):
        """!
        @brief A command that does not exist is an attribute error.
        """
        with self.assertRaises(AttributeError) as caught:
            act.commands.nosuchcommand  # pylint: disable=pointless-statement
        self.assertIn("work", str(caught.exception))
        self.assertIn("work", dir(act.commands))


class LegacyTest(unittest.TestCase):
    """!
    @brief The class ACT of the earlier module act.py.
    """

    def test_old_name_warns_or_explains(self):
        """!
        @brief act.ACT gives the old class with a warning when ACT-utilities
        is installed, and otherwise an error that says where the class is.
        """
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                old = act.ACT
            except ImportError as error:
                self.assertIn("ACT-utilities", str(error))
                return
        self.assertTrue(inspect.isclass(old))
        self.assertTrue(any(issubclass(w.category, DeprecationWarning) for w in caught))

    def test_other_names_do_not_exist(self):
        """!
        @brief Other names are not invented.
        """
        with self.assertRaises(AttributeError):
            act.NoSuchName  # pylint: disable=pointless-statement

    def test_class_comes_from_act_legacy(self):
        """!
        @brief act.legacy takes the class from the module act_legacy of
        ACT-utilities.
        """
        saved = {name: sys.modules.pop(name, None) for name in ("act.legacy", "act_legacy")}
        had_attribute = getattr(act, "legacy", None)
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "act_legacy.py").write_text(
                "class ACT:\n    marker = 'stand-in'\n"
                "class Alg:\n    pass\nclass Target:\n    pass\n")
            sys.path.insert(0, tmp)
            try:
                import importlib  # pylint: disable=import-outside-toplevel
                legacy = importlib.import_module("act.legacy")
                self.assertEqual("stand-in", legacy.ACT.marker)
                self.assertTrue(inspect.isclass(legacy.Alg))
                self.assertTrue(inspect.isclass(legacy.Target))
            finally:
                sys.path.remove(tmp)
                for name, module in saved.items():
                    sys.modules.pop(name, None)
                    if module is not None:
                        sys.modules[name] = module
                if had_attribute is None:
                    if "legacy" in vars(act):
                        delattr(act, "legacy")
                else:
                    setattr(act, "legacy", had_attribute)


class ProgramLookupTest(unittest.TestCase):
    """!
    @brief Which alexandria program is run when none is chosen.
    """

    def setUp(self):
        """!
        @brief No choice in the settings or in the environment, and a
        directory for stand-in programs.
        """
        self._executable = act.config.executable
        act.config.executable = None
        self._environment = {name: os.environ.get(name) for name in ("ACT_ALEXANDRIA", "PATH")}
        os.environ.pop("ACT_ALEXANDRIA", None)
        self._installed = sys.modules.pop("act._installed", None)
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(os.path.realpath(self._tmp.name))

    def tearDown(self):
        """!
        @brief Restore the settings and the environment.
        """
        act.config.executable = self._executable
        for name, value in self._environment.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
        sys.modules.pop("act._installed", None)
        if self._installed is not None:
            sys.modules["act._installed"] = self._installed
        if "_installed" in vars(act) and self._installed is None:
            delattr(act, "_installed")
        self._tmp.cleanup()

    def program(self, name: str) -> Path:
        """!
        @brief Make a stand-in program.
        @param name The path of the program relative to the directory.
        @return The path of the program.
        """
        path = self.dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("#!/bin/sh\nexit 0\n")
        path.chmod(0o755)
        return path

    def pretend_installed(self, relative: str):
        """!
        @brief Pretend that the installation wrote the module _installed.
        @param relative The path of the program relative to the package.
        """
        import types  # pylint: disable=import-outside-toplevel
        module = types.ModuleType("act._installed")
        module.EXECUTABLE = relative
        sys.modules["act._installed"] = module

    def test_double_precision_name(self):
        """!
        @brief alexandria_d is used if it is the only one in the PATH.
        """
        if self._installed is not None:
            self.skipTest("The package was installed with a program")
        program = self.program("bin/alexandria_d")
        os.environ["PATH"] = str(program.parent)
        self.assertEqual([str(program)], act.config.executable_words())
        # The single precision program goes first
        single = self.program("bin/alexandria")
        self.assertEqual([str(single)], act.config.executable_words())

    def test_program_of_the_installation(self):
        """!
        @brief The program that was installed with the package is found
        relative to the package, and only if it is there.
        """
        package = self.dir / "lib" / "site-packages" / "act"
        package.mkdir(parents=True)
        self.pretend_installed("../../../bin/alexandria")
        self.assertIsNone(act.runner.installed_program(package))
        program = self.program("bin/alexandria")
        self.assertEqual(str(program), act.runner.installed_program(package))

    def test_no_installation(self):
        """!
        @brief Without the module of the installation there is no such
        program, and an empty or wrong module is ignored.
        """
        if self._installed is None:
            self.assertIsNone(act.runner.installed_program(self.dir))
        self.pretend_installed("")
        self.assertIsNone(act.runner.installed_program(self.dir))
        self.pretend_installed("nothing-here")
        self.assertIsNone(act.runner.installed_program(self.dir))

    def test_environment_goes_first(self):
        """!
        @brief ACT_ALEXANDRIA goes before the PATH.
        """
        self.program("bin/alexandria")
        other = self.program("other/my-alexandria")
        os.environ["PATH"] = str(self.dir / "bin")
        os.environ["ACT_ALEXANDRIA"] = str(other)
        self.assertEqual([str(other)], act.config.executable_words())


if __name__ == "__main__":
    unittest.main()
