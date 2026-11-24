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
@brief Tests of the entries of the index of the manual that
docs/manual/act_index.py makes from the source code: the commands of
alexandria, their options and the variables of CMake that a text mentions.
"""
import importlib.util
import os
import unittest
from typing import Any

## The extension of the manual that makes the entries
EXTENSION = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "docs", "manual", "act_index.py")


def load_extension() -> Any:
    """!
    @brief Import docs/manual/act_index.py.
    @return The module.
    """
    spec = importlib.util.spec_from_file_location("act_index", EXTENSION)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ManualIndexTest(unittest.TestCase):
    """!
    @brief The names that the index finds in the text of the manual.
    """

    @classmethod
    def setUpClass(cls) -> None:
        """!
        @brief Read the commands and options from the source code.
        """
        cls.ext = load_extension()
        cls.commands = cls.ext.read_commands()
        cls.options = cls.ext.read_options()
        cls.terms = cls.ext.Terms(cls.commands, cls.options)

    def test_commands_from_source(self) -> None:
        """!
        @brief The commands are those that alexandria registers.
        """
        for name in ("train_ff", "geometry_ff", "simulate", "b2", "help"):
            self.assertIn(name, self.commands)
        self.assertNotIn("program", self.commands)

    def test_options_from_source(self) -> None:
        """!
        @brief The options are those of the code, without the dash.
        """
        for name in ("ff", "lockstep", "shell_keep", "sensitivity", "h"):
            self.assertIn(name, self.options)
        self.assertNotIn("bonded", self.options)

    def test_command_line(self) -> None:
        """!
        @brief A command line of alexandria gives the command and its options.
        """
        found = self.terms.find("alexandria train_ff -ff ff.xml -lockstep 8 -o out.xml")
        self.assertEqual(found, [("command", "train_ff"), ("option", "ff"),
                                 ("option", "lockstep"), ("option", "o")])

    def test_words_are_not_options(self) -> None:
        """!
        @brief Hyphens in words and unknown options are not indexed.
        """
        self.assertEqual(self.terms.find("non-bonded interactions of ACS-g with -1 kJ/mol"), [])
        self.assertEqual(self.terms.find("the option -nonexistent_option"), [])

    def test_boolean_option(self) -> None:
        """!
        @brief -noname is the option -name.
        """
        self.assertEqual(self.terms.find("-noshell_keep"), [("option", "shell_keep")])

    def test_short_options_need_alexandria(self) -> None:
        """!
        @brief Options of one letter of other programs are not indexed.
        """
        self.assertEqual(self.terms.find("conda create -n act -c conda-forge"), [])
        self.assertEqual(self.terms.find("alexandria minimize -c out.pdb"),
                         [("command", "minimize"), ("option", "c")])

    def test_ambiguous_commands(self) -> None:
        """!
        @brief A command that is an English word counts only on its own.
        """
        self.assertEqual(self.terms.find("minimize the energy"), [])
        self.assertEqual(self.terms.find("minimize", whole=True), [("command", "minimize")])
        self.assertEqual(self.terms.find("simulate -minimize", whole=True),
                         [("command", "simulate"), ("option", "minimize")])
        self.assertEqual(self.terms.find("written by geometry_ff"), [("command", "geometry_ff")])

    def test_cmake(self) -> None:
        """!
        @brief Variables of CMake have an underscore; -DOPTION=Value is an example.
        """
        self.assertEqual(self.terms.find("cmake -DACT_GPU=SYCL -DOPTION=Value ..", options=False),
                         [("cmake", "ACT_GPU")])

    def test_command_lines(self) -> None:
        """!
        @brief A line that starts with a dash continues the command line before it.
        """
        lines = self.ext.command_lines("alexandria b2 -ff TIP4P\n-T1 300\nconda list -c x")
        self.assertEqual(lines, [("alexandria b2 -ff TIP4P\n-T1 300", True), ("conda list -c x", False)])

    def test_entries(self) -> None:
        """!
        @brief The text of the entries.
        """
        self.assertEqual(self.ext.entry("command", "nma"), ("pair", "alexandria command; nma"))
        self.assertEqual(self.ext.entry("option", "ff"), ("single", "-ff (command-line option)"))
        self.assertEqual(self.ext.entry("cmake", "ACT_GPU"), ("single", "ACT_GPU (CMake variable)"))


if __name__ == "__main__":
    unittest.main()
