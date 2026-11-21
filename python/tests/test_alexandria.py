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
import time
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

    def test_force_field_as_json(self):
        """!
        @brief edit_ff describes a force field, and a change shows up in
        the description and in the file that is written.
        """
        with act.commands.edit_ff(ff=self.ff, json=True) as res:
            described = res.data["forcefield"]
            self.assertIsInstance(described["polarizable"], bool)
            kinds = {inter["type"]: inter for inter in described["interactions"]}
            self.assertIn("BONDS", kinds)
            first = kinds["BONDS"]["parameterlists"][0]
            self.assertIn("identifier", first)
            for parameter in first["parameters"]:
                self.assertLessEqual(parameter["minimum"], parameter["maximum"])
                self.assertIsInstance(parameter["nonnegative"], bool)
            self.assertTrue(res.data["validation"]["ok"])
            self.assertEqual([], res.data["validation"]["issues"])
            masses = [param["value"] for ptype in described["particletypes"]["types"]
                      for param in ptype["parameters"] if param["type"] == "mass"]
            self.assertGreater(len(masses), 5)
            self.assertGreater(min(masses), 0)

        def chi(data, identifier):
            """The electronegativity of a type in the description"""
            for inter in data["forcefield"]["interactions"]:
                if inter["type"] == "ELECTRONEGATIVITYEQUALIZATION":
                    for entry in inter["parameterlists"]:
                        if entry["identifier"] == identifier:
                            return [p for p in entry["parameters"] if p["type"] == "chi"][0]
            self.fail("No chi of %s in the description" % identifier)
            return None
        acm = self.data / "forcefields" / "ACM-g1.xml"
        with act.commands.edit_ff(ff=acm, json=True) as res:
            before = chi(res.data, "h_z")
        new = 0.5 * (before["value"] + before["maximum"])
        with act.commands.edit_ff(ff=acm, p="chi", a="h_z", val=new, o=True, json=True) as res:
            self.assertAlmostEqual(new, chi(res.data, "h_z")["value"], places=5)
            self.assertNotAlmostEqual(before["value"], chi(res.data, "h_z")["value"], places=3)
            # The file that was written has the same value
            written = res.file("o")
            self.assertTrue(written.is_file())
            with act.commands.edit_ff(ff=written, json=True) as again:
                self.assertAlmostEqual(new, chi(again.data, "h_z")["value"], places=5)

    def test_generated_force_field_as_json(self):
        """!
        @brief gen_ff describes the force field it makes, and geometry_ff
        the one it fills with the geometries of a molprop file.
        """
        atomtypes = self.data / "forcefields" / "genff-atomtypes.csv"
        with act.commands.gen_ff(f=atomtypes, o=True, json=True) as res:
            self.assertTrue(res.file("o").is_file())
            described = res.data["forcefield"]
            names = [ptype["identifier"] for ptype in described["particletypes"]["types"]]
            self.assertIn("hf", names)
            kinds = {inter["type"]: inter["function"] for inter in described["interactions"]}
            self.assertEqual("LJ12_6", kinds["VANDERWAALS"])
            self.assertTrue(res.data["validation"]["ok"])
        with act.commands.gen_ff(f=atomtypes, vdwfn="LJ14_7", o=True, json=True) as res:
            kinds = {inter["type"]: inter["function"]
                     for inter in res.data["forcefield"]["interactions"]}
            self.assertEqual("LJ14_7", kinds["VANDERWAALS"])
        with act.commands.geometry_ff(ff=self.data / "forcefields" / "ACS-g1.xml",
                                      mp=[self.data / "molprop" / "molprop.xml"],
                                      o=True, json=True) as res:
            self.assertEqual(3, res.data["compounds"])
            self.assertGreater(res.data["extracted"]["BONDS"], 3)
            bonds = [inter for inter in res.data["forcefield"]["interactions"]
                     if inter["type"] == "BONDS"][0]
            self.assertEqual(res.data["extracted"]["BONDS"], len(bonds["parameterlists"]))

    def test_molprop_file_as_json(self):
        """!
        @brief edit_mp describes the molecules of a molprop file and says
        whether they are consistent.
        """
        molprop = self.data / "molprop" / "ion-water.xml"
        with act.commands.edit_mp(mp=[molprop], o=True, validate=True, json=True) as res:
            summary = res.data["summary"]
            self.assertEqual(len(res.data["molecules"]), summary["molecules"])
            self.assertGreater(summary["dimers"], 0)
            self.assertEqual(summary["experiments"],
                             sum(mol["experiments"] for mol in res.data["molecules"]))
            first = res.data["molecules"][0]
            self.assertEqual(2, len(first["fragments"]))
            self.assertIn("InteractionEnergy", first["properties"])
            self.assertEqual({"checked": True, "ok": True, "problems": 0, "invalid": [],
                              "duplicates": 0, "with_duplicates": [], "messages": []},
                             res.data["validation"])
            self.assertTrue(res.file("o").is_file())

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

    def test_training_reports_progress(self):
        """!
        @brief Training a force field tells how far it is while it runs.
        """
        seen = []
        ffdir = self.data / "forcefields"
        mpdir = self.data / "molprop"
        with act.commands.train_ff(ff=ffdir / "halgren-PG.xml", mp=mpdir / "booze-mono.xml",
                                   sel=mpdir / "booze-mono.dat", optimizer="MCMC", pop_size=1,
                                   max_iter=3, seed=219463979, fit="De b c bondlength angle kt",
                                   fc_epot=1.0, on_progress=seen.append) as res:
            self.assertGreaterEqual(len(seen), 3)
            start = seen[0]
            self.assertEqual("start", start["event"])
            self.assertEqual("MCMC", start["optimizer"])
            self.assertEqual(1, start["ngenerations"])
            self.assertEqual(3, start["niterations"])
            self.assertGreater(start["nparameters"], 0)
            self.assertGreater(start["chi2"]["train"], 0)
            for item in seen[1:-1]:
                self.assertEqual("step", item["event"])
                self.assertLess(item["iteration"], 3)
            end = seen[-1]
            self.assertEqual("generation", end["event"])
            self.assertEqual(1, end["generation"])
            self.assertLessEqual(end["best"]["train"], start["chi2"]["train"])
            self.assertEqual("s", end["elapsed"]["unit"])
            self.assertTrue(res.file("o").exists())

    def test_training_stops_at_the_time_limit(self):
        """!
        @brief A training with more iterations than there is time for stops
        when the time of -maxh is up, and ends with a force field as a
        complete training does, without the sensitivity analysis.
        """
        seen = []
        ffdir = self.data / "forcefields"
        mpdir = self.data / "molprop"
        began = time.monotonic()
        # Two seconds, for iterations that would take days
        with act.commands.train_ff(ff=ffdir / "halgren-PG.xml", mp=mpdir / "booze-mono.xml",
                                   sel=mpdir / "booze-mono.dat", optimizer="MCMC", pop_size=1,
                                   max_iter=100000000, seed=219463979, fit="sigma epsilon",
                                   fc_epot=1.0, maxh=2.0/3600, on_progress=seen.append,
                                   timeout=300) as res:
            self.assertLess(time.monotonic() - began, 120)
            self.assertEqual("yes", res.data["time_limit"]["reached"])
            self.assertEqual("h", res.data["time_limit"]["limit"]["unit"])
            self.assertAlmostEqual(2.0/3600, res.data["time_limit"]["limit"]["value"], places=6)
            self.assertNotIn("sensitivity", res.data)
            self.assertIn("(-maxh) was reached", res.log)
            # It got somewhere, and wrote the force field of where it got
            steps = [item for item in seen if item["event"] == "step"]
            self.assertGreater(len(steps), 0)
            self.assertLess(steps[-1]["iteration"], 100000000)
            self.assertEqual("generation", seen[-1]["event"])
            self.assertLess(seen[-1]["best"]["train"], seen[0]["chi2"]["train"])
            self.assertTrue(res.file("o").exists())
            self.assertGreater(res.file("o").stat().st_size, 1000)

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
