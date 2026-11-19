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
@brief Tests of the Python interface with the real alexandria program.

They run when the environment variable ACT_ALEXANDRIA points to the program
and ACT_TESTDATA to the directory src/act/testdata of the source of ACT, as
ctest sets them. Otherwise they are skipped.
"""
import math
import os
import sys
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import act  # noqa: E402  pylint: disable=wrong-import-position

## The program, if it is given
ALEXANDRIA = os.environ.get("ACT_ALEXANDRIA", "")
## The test data of ACT, if they are given
TESTDATA = os.environ.get("ACT_TESTDATA", "")
## All commands of alexandria
COMMANDS = ["b2", "edit_ff", "edit_mp", "gen_ff", "gentop", "geometry_ff",
            "merge_ff", "min_complex", "minimize", "nma", "simulate", "train_ff"]


@unittest.skipUnless(ALEXANDRIA and os.path.isdir(TESTDATA),
                     "ACT_ALEXANDRIA and ACT_TESTDATA are not set")
class AlexandriaTest(unittest.TestCase):
    """!
    @brief The commands of alexandria from Python.
    """

    @classmethod
    def setUpClass(cls):
        """!
        @brief Use the program that was given.
        """
        cls._executable = act.config.executable
        act.config.executable = ALEXANDRIA
        act.commands._functions.clear()  # pylint: disable=protected-access
        cls.data = Path(TESTDATA)
        cls.ff = cls.data / "forcefields" / "ACS-g1.xml"
        cls.polff = cls.data / "forcefields" / "ACS-pg1.xml"
        cls.acetone = cls.data / "molecules" / "acetone.sdf"

    @classmethod
    def tearDownClass(cls):
        """!
        @brief Restore the settings.
        """
        act.config.executable = cls._executable
        act.commands._functions.clear()  # pylint: disable=protected-access

    def test_all_commands_are_described(self):
        """!
        @brief Every command has a function with a help text.
        """
        self.assertEqual(COMMANDS, act.commands.available())
        for name in COMMANDS:
            function = getattr(act.commands, name)
            self.assertEqual(name, function.__name__)
            self.assertIn("Options:", function.__doc__)

    def test_run_arguments_do_not_hide_options(self):
        """!
        @brief No command has an option with the name of an argument of run(),
        and no two options of a command have the same Python name.
        """
        info = act.program_info()
        for command in info.commands.values():
            names = [option.python_name for option in command.options.values()]
            self.assertEqual(len(names), len(set(names)), command.name)
            for reserved in act.runner.RUN_ARGUMENTS:
                self.assertNotIn(reserved, names, command.name)

    def test_minimize(self):
        """!
        @brief Minimize acetone and read the results.
        """
        with act.commands.minimize(ff=self.ff, f=self.acetone) as res:
            data = res.data
            self.assertEqual("OK", data["status"])
            self.assertIs(True, data["converged"])
            self.assertEqual("nm", data["rmsd"]["unit"])
            self.assertTrue(0 < data["rmsd"]["value"] < 0.1)
            before = data["Energies before"]["EPOT"]["value"]
            after = data["Energies after"]["EPOT"]["value"]
            self.assertLess(after, before)
            molecule = data["molecule"]
            self.assertEqual("C3H6O", molecule["formula"])
            self.assertEqual(10, len(molecule["particles"]))
            self.assertEqual(9, len(molecule["bonds"]))
            self.assertEqual(10, len(data["coordinates"]["after"]))
            self.assertTrue(res.file("c").exists())
            self.assertIn("minimiz", res.log.lower())

    def test_nma(self):
        """!
        @brief Normal modes of acetone: frequencies and displacements.
        """
        with act.commands.nma(ff=self.ff, f=self.acetone, ir=True) as res:
            data = res.data
            modes = data["modes"]
            self.assertEqual(24, len(modes))
            for mode in modes:
                self.assertEqual(10, len(mode["displacements"]))
                norm = sum(x*x for d in mode["displacements"] for x in d)
                self.assertTrue(math.isclose(1.0, norm, abs_tol=1e-3))
                self.assertGreater(mode["reduced_mass"]["value"], 1.0)
            self.assertTrue(res.file("ir").exists())
            self.assertIn("Thermochemistry", data)

    def test_polarizable_model_has_shells(self):
        """!
        @brief With a polarizable force field the particles include shells.
        """
        with act.commands.minimize(ff=self.polff, f=self.acetone) as res:
            kinds = [p["particletype"] for p in res.data["molecule"]["particles"]]
            self.assertEqual(10, kinds.count("Atom"))
            self.assertEqual(10, kinds.count("Shell"))

    def test_molecule_from_smiles(self):
        """!
        @brief A molecule can be given as a SMILES string.
        """
        with act.commands.minimize(ff=self.ff, smiles="CC(=O)C") as res:
            self.assertEqual("C3H6O", res.data["molecule"]["formula"])
            self.assertIs(True, res.data["converged"])
        with self.assertRaises(act.ActError):
            act.commands.minimize(ff=self.ff, smiles="C(C")

    def test_simulation_reports_progress(self):
        """!
        @brief A simulation tells how far it is while it runs.
        """
        seen = []
        with act.commands.simulate(ff=self.ff, f=self.acetone, nsteps=200, nstener=20,
                                   temp=300.0, seed=1993, on_progress=seen.append) as res:
            self.assertEqual(list(range(0, 200, 20)), [item["step"] for item in seen])
            self.assertEqual(200, seen[0]["nsteps"])
            self.assertEqual("kJ/mol", seen[0]["energies"]["unit"])
            for item in seen:
                self.assertIn("EPOT", item["energies"])
                self.assertGreaterEqual(item["temperature"]["value"], 0)
            self.assertTrue(res.file("e").exists())
            self.assertTrue(res.file("o").exists())
        with self.assertRaises(act.ActOptionError):
            act.commands.minimize(ff=self.ff, f=self.acetone, on_progress=seen.append)

    def test_missing_molecule_file(self):
        """!
        @brief The program reports a molecule file that does not exist.
        """
        with self.assertRaises(act.ActError) as caught:
            act.commands.minimize(ff=self.ff, f=self.data / "molecules" / "nothere.sdf")
        self.assertEqual(1, caught.exception.returncode)
        self.assertIn("nothere.sdf", str(caught.exception))

    def test_missing_force_field_file(self):
        """!
        @brief A force field file that does not exist is caught before the run.
        """
        with self.assertRaises(act.ActNotFoundError):
            act.commands.minimize(ff=self.data / "nothere.xml", f=self.acetone)

    def test_wrong_value_reported_by_the_program(self):
        """!
        @brief An error of the program in reading its options gives its message.
        """
        # A file with an extension that the option does not take
        with self.assertRaises(act.ActError) as caught:
            act.commands.minimize(ff=self.acetone, f=self.acetone)
        self.assertIn("exit status", str(caught.exception))
        self.assertNotIn("MPI_ABORT", str(caught.exception))

    def test_option_cannot_be_smuggled_in(self):
        """!
        @brief A value cannot pose as another option.
        """
        with self.assertRaises(act.ActOptionError):
            act.commands.minimize(ff=self.ff, f=self.acetone, qalg="-maxiter")


if __name__ == "__main__":
    unittest.main()
