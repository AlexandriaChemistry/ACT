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
@brief Tests of the description of the options and of the conversion of
keyword arguments to a command line. These tests do not run a program.
"""
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import act  # noqa: E402  pylint: disable=wrong-import-position
from act.options import ProgramInfo, python_name  # noqa: E402  pylint: disable=wrong-import-position
from tests.fake_alexandria import DESCRIPTION  # noqa: E402  pylint: disable=wrong-import-position


class OptionsTest(unittest.TestCase):
    """!
    @brief The description of the commands and their options.
    """

    def setUp(self):
        """!
        @brief Read the description of the stand-in program.
        """
        self.info = ProgramInfo.from_json(DESCRIPTION)
        self.work = self.info.find_command("work")

    def test_commands_and_options(self):
        """!
        @brief Commands and options are found with their properties.
        """
        self.assertEqual(["work"], list(self.info.commands))
        self.assertEqual("Pretend to do a calculation.", self.work.description)
        ff = self.work.find_option("ff")
        self.assertTrue(ff.is_file and ff.is_input and not ff.is_output)
        self.assertEqual((".xml",), ff.extensions)
        self.assertTrue(self.work.find_option("c").is_output)
        self.assertEqual(("LBFGS", "Newton", "Steep"), self.work.find_option("minalg").values)
        self.assertEqual(3, self.work.find_option("box").max_values)
        self.assertTrue(self.work.find_option("fail").hidden)

    def test_unknown_command_suggests(self):
        """!
        @brief An unknown command gives an error with a suggestion.
        """
        with self.assertRaises(act.ActOptionError) as caught:
            self.info.find_command("wrok")
        self.assertIn("'work'", str(caught.exception))

    def test_unknown_option_suggests(self):
        """!
        @brief An unknown option gives an error with a suggestion.
        """
        with self.assertRaises(act.ActOptionError) as caught:
            self.work.find_option("maxitr")
        self.assertIn("'maxiter'", str(caught.exception))

    def test_python_names(self):
        """!
        @brief Reserved words and odd characters give usable argument names.
        """
        self.assertEqual("del_", python_name("del"))
        self.assertEqual("maxiter", python_name("maxiter"))
        self.assertEqual("a_b", python_name("a-b"))
        self.assertEqual("_2d", python_name("2d"))
        self.assertEqual("del", self.work.find_option("del_").name)

    def test_other_format_is_refused(self):
        """!
        @brief A description with another layout is refused with a clear message.
        """
        other = dict(DESCRIPTION)
        other["format"] = 2
        with self.assertRaises(act.ActError) as caught:
            ProgramInfo.from_json(other)
        self.assertIn("format 2", str(caught.exception))
        with self.assertRaises(act.ActError):
            ProgramInfo.from_json({"nothing": 1})

    def test_unknown_type_is_refused(self):
        """!
        @brief An option with a type this module does not know is refused.
        """
        bad = {"program": "alexandria", "format": 1,
               "commands": {"x": {"options": [{"name": "a", "type": "matrix"}]}}}
        with self.assertRaises(act.ActError):
            ProgramInfo.from_json(bad)

    def test_run_arguments_are_no_options(self):
        """!
        @brief The arguments of run() do not hide options of the program.
        """
        for name in act.runner.RUN_ARGUMENTS:
            self.assertFalse(self.work.has_option(name), name)


class ArgumentsTest(unittest.TestCase):
    """!
    @brief Keyword arguments become a command line.
    """

    def setUp(self):
        """!
        @brief Make a directory with input files and a directory to run in.
        """
        self.work = ProgramInfo.from_json(DESCRIPTION).find_command("work")
        self._tmp = tempfile.TemporaryDirectory()
        self.base = Path(os.path.realpath(self._tmp.name)) / "base"
        self.rundir = Path(os.path.realpath(self._tmp.name)) / "run"
        self.base.mkdir()
        self.rundir.mkdir()
        for name in ("ff.xml", "a.xml", "b.xml", "mol.sdf"):
            (self.base / name).write_text("x")

    def tearDown(self):
        """!
        @brief Remove the directories.
        """
        self._tmp.cleanup()

    def build(self, **options):
        """!
        @brief Convert keyword arguments.
        @param options The keyword arguments.
        @return The arguments and the expected output files.
        """
        return act.build_arguments(self.work, options, self.rundir, self.base)

    def test_booleans(self):
        """!
        @brief True gives -name, False gives -noname, other values are refused.
        """
        self.assertEqual(["-oneH"], self.build(oneH=True)[0])
        self.assertEqual(["-nooneH"], self.build(oneH=False)[0])
        self.assertEqual(["-del"], self.build(del_=True)[0])
        with self.assertRaises(act.ActOptionError):
            self.build(oneH=1)
        with self.assertRaises(act.ActOptionError):
            self.build(oneH="yes")

    def test_none_is_not_given(self):
        """!
        @brief An argument that is None leaves the option out.
        """
        self.assertEqual([], self.build(maxiter=None, oneH=None, ff=None)[0])

    def test_numbers(self):
        """!
        @brief Integers and reals are checked and written in full precision.
        """
        self.assertEqual(["-maxiter", "25"], self.build(maxiter=25)[0])
        self.assertEqual(["-toler", "1e-07"], self.build(toler=1e-7)[0])
        self.assertEqual(["-toler", "2.0"], self.build(toler=2)[0])
        self.assertEqual(["-toler", "0.1"], self.build(toler=0.1)[0])
        self.assertEqual(["-maxiter", "-3"], self.build(maxiter=-3)[0])
        for wrong in (2.5, "25", True, None.__class__):
            with self.assertRaises(act.ActOptionError):
                self.build(maxiter=wrong)
        for wrong in ("0.1", True, float("nan"), float("inf")):
            with self.assertRaises(act.ActOptionError):
                self.build(toler=wrong)

    def test_enum(self):
        """!
        @brief Only the values of an enum are accepted.
        """
        self.assertEqual(["-minalg", "Newton"], self.build(minalg="Newton")[0])
        with self.assertRaises(act.ActOptionError) as caught:
            self.build(minalg="newton")
        self.assertIn("LBFGS, Newton, Steep", str(caught.exception))

    def test_several_values(self):
        """!
        @brief Lists go to options that take several values, up to their maximum.
        """
        self.assertEqual(["-box", "1.0", "2.0", "3.5"], self.build(box=[1, 2, 3.5])[0])
        self.assertEqual(["-box", "4.0"], self.build(box=4)[0])
        with self.assertRaises(act.ActOptionError):
            self.build(box=[1, 2, 3, 4])
        with self.assertRaises(act.ActOptionError):
            self.build(box=[])
        with self.assertRaises(act.ActOptionError):
            self.build(maxiter=[1, 2])

    def test_strings(self):
        """!
        @brief Strings are passed as they are, but cannot pose as options.
        """
        self.assertEqual(["-name", "water dimer"], self.build(name="water dimer")[0])
        # A shell would expand these, but there is no shell
        self.assertEqual(["-name", "$(rm -rf x); `id` * > f"],
                         self.build(name="$(rm -rf x); `id` * > f")[0])
        self.assertEqual(["-name", "-0.5"], self.build(name="-0.5")[0])
        with self.assertRaises(act.ActOptionError) as caught:
            self.build(name="-ff")
        self.assertIn("dash", str(caught.exception))
        with self.assertRaises(act.ActOptionError):
            self.build(name="-0.8 0.4 0.4")
        with self.assertRaises(act.ActOptionError):
            self.build(name="a\0b")
        with self.assertRaises(act.ActOptionError):
            self.build(name=12)

    def test_molecule_file_name(self):
        """!
        @brief The molecule file, a string option, is found from the directory
        of the caller.
        """
        self.assertEqual(["-f", str(self.base / "mol.sdf")], self.build(f="mol.sdf")[0])
        self.assertEqual(["-f", str(self.base / "mol.sdf")], self.build(f=Path("mol.sdf"))[0])
        # A name that is not a file is left alone, the program reports it
        self.assertEqual(["-f", "nothere.sdf"], self.build(f="nothere.sdf")[0])
        # Other strings are never changed into paths
        self.assertEqual(["-name", "mol.sdf"], self.build(name="mol.sdf")[0])

    def test_input_files(self):
        """!
        @brief Input files must exist and are made absolute.
        """
        self.assertEqual(["-ff", str(self.base / "ff.xml")], self.build(ff="ff.xml")[0])
        self.assertEqual(["-ff", str(self.base / "ff.xml")],
                         self.build(ff=self.base / "ff.xml")[0])
        with self.assertRaises(act.ActNotFoundError) as caught:
            self.build(ff="nothere.xml")
        self.assertIn("nothere.xml", str(caught.exception))
        with self.assertRaises(act.ActOptionError):
            self.build(ff=True)
        with self.assertRaises(act.ActOptionError):
            self.build(ff="")
        with self.assertRaises(act.ActOptionError):
            self.build(ff=["ff.xml", "a.xml"])

    def test_several_input_files(self):
        """!
        @brief An option for several files takes a list.
        """
        args = self.build(many=["a.xml", "b.xml"])[0]
        self.assertEqual(["-many", str(self.base / "a.xml"), str(self.base / "b.xml")], args)
        self.assertEqual(["-many", str(self.base / "a.xml")], self.build(many="a.xml")[0])
        with self.assertRaises(act.ActNotFoundError):
            self.build(many=["a.xml", "nothere.xml"])

    def test_library_file(self):
        """!
        @brief A library file that is given by name only is left to the program.
        """
        self.assertEqual(["-lib", "atomprops.csv"], self.build(lib="atomprops.csv")[0])
        with self.assertRaises(act.ActNotFoundError):
            self.build(lib="sub/atomprops.csv")

    def test_output_files(self):
        """!
        @brief Output files can be a path, True or False.
        """
        args, outputs = self.build(extra=True)
        self.assertEqual(["-extra", str(self.rundir / "extra.xvg")], args)
        self.assertEqual([self.rundir / "extra.xvg"], outputs["extra"])
        args, outputs = self.build(extra="plot.xvg")
        self.assertEqual(["-extra", str(self.base / "plot.xvg")], args)
        self.assertEqual([self.base / "plot.xvg"], outputs["extra"])
        args, outputs = self.build(extra=False)
        self.assertEqual([], args)
        self.assertNotIn("extra", outputs)
        with self.assertRaises(act.ActNotFoundError):
            self.build(extra="no/such/dir/plot.xvg")
        # A file that is always written cannot be switched off
        with self.assertRaises(act.ActOptionError):
            self.build(c=False)

    def test_files_that_are_always_written(self):
        """!
        @brief Files with default names are expected in the run directory.
        """
        outputs = self.build()[1]
        self.assertEqual([self.rundir / "confout.pdb"], outputs["c"])
        self.assertEqual([self.rundir / "work.log"], outputs["g"])
        self.assertNotIn("extra", outputs)
        self.assertNotIn("ff", outputs)
        outputs = self.build(c="min.pdb")[1]
        self.assertEqual([self.base / "min.pdb"], outputs["c"])

    def test_unknown_and_double_options(self):
        """!
        @brief Unknown options and options given twice are refused.
        """
        with self.assertRaises(act.ActOptionError):
            self.build(nothere=1)
        with self.assertRaises(act.ActOptionError):
            act.build_arguments(self.work, {"del": True, "del_": True}, self.rundir, self.base)


if __name__ == "__main__":
    unittest.main()
