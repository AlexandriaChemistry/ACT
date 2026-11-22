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
@brief Tests of the classes of ACT in Python, which use the code of ACT
directly instead of running the alexandria program. They need ACT built
with -DACT_PYTHON=ON and are skipped otherwise.
"""
import os
import signal
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from xml.etree import ElementTree

# The package of the build tree has the module with the classes; the tests
# of a build say where it is. Otherwise the package of the source code.
sys.path.insert(0, os.environ.get("ACT_PYTHON_PATH")
                or os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import act  # noqa: E402  pylint: disable=wrong-import-position

## The test data of ACT
TESTDATA = os.environ.get("ACT_TESTDATA", "")


def plain(value):
    """!
    @brief Lists instead of numpy arrays, also in dicts and lists, for comparisons.
    @param value What the classes returned.
    @return The same with lists.
    """
    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, dict):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    return value


class WithoutClassesTest(unittest.TestCase):
    """!
    @brief What a Python gets that has the package without its classes.
    """

    def test_tells_whether_the_classes_are_there(self):
        """!
        @brief The package says whether it was built with its classes, and
        works without them.
        """
        self.assertIsInstance(act.classes_available(), bool)
        self.assertTrue(callable(act.run))

    @unittest.skipUnless(os.environ.get("ACT_REQUIRE_CLASSES"), "The classes are not required")
    def test_classes_are_there_when_the_build_has_them(self):
        """!
        @brief A build with the classes tests them: nothing is skipped
        because the module could not be loaded.
        """
        self.assertTrue(act.classes_available())
        self.assertTrue(os.path.isdir(TESTDATA), "ACT_TESTDATA is not set")

    def test_data_directory(self):
        """!
        @brief The package tells where the classes look for the data files of
        ACT: a directory that exists, or nothing.
        """
        where = act.data_directory()
        self.assertIsInstance(where, str)
        if where:
            self.assertTrue(os.path.isdir(where), where)
        if not act.classes_available():
            self.assertEqual("", where)

    @unittest.skipIf(act.classes_available(), "ACT was built with its Python classes")
    def test_asking_for_a_class_says_how_to_get_it(self):
        """!
        @brief Without the classes, asking for one tells what is missing.
        """
        with self.assertRaises(ImportError) as caught:
            act.ForceField  # pylint: disable=pointless-statement
        self.assertIn("ACT_PYTHON", str(caught.exception))


@unittest.skipUnless(act.classes_available(), "ACT was built without its Python classes")
@unittest.skipUnless(os.path.isdir(TESTDATA), "ACT_TESTDATA is not set")
class ForceFieldTest(unittest.TestCase):
    """!
    @brief A force field in Python.
    """

    def setUp(self):
        """!
        @brief Find the test data and make a directory for files.
        """
        self.ffdir = Path(TESTDATA) / "forcefields"
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self._tmp.name)

    def tearDown(self):
        """!
        @brief Remove the files.
        """
        self._tmp.cleanup()

    def test_read_a_force_field(self):
        """!
        @brief A force field file is read, and tells what is in it.
        """
        ff = act.ForceField.read(self.ffdir / "ACS-pg1.xml")
        self.assertEqual(str(self.ffdir / "ACS-pg1.xml"), ff.filename)
        self.assertTrue(ff.polarizable)
        self.assertFalse(act.ForceField.read(str(self.ffdir / "ACS-g1.xml")).polarizable)
        # The particle types by name, atoms and shells
        types = ff.particle_types
        self.assertGreater(len(types), 20)
        self.assertEqual(len(types), len(set(types)))
        self.assertIn("ow", types)
        self.assertIn("c3", types)
        self.assertIn("c3_s", types)
        # The interactions it has, by the names of the file
        interactions = ff.interactions
        self.assertIn("BONDS", interactions)
        self.assertIn("VANDERWAALS", interactions)
        self.assertIn("POLARIZATION", interactions)
        self.assertEqual(len(interactions), len(set(interactions)))
        self.assertIn("ACS-pg1.xml", repr(ff))
        self.assertIn(str(len(types)), repr(ff))

    def test_file_that_cannot_be_read(self):
        """!
        @brief A file that is missing or broken raises an error that tells
        which file, and Python lives on.
        """
        with self.assertRaises(act.ActFileError) as caught:
            act.ForceField.read(self.dir / "no-such-force-field.xml")
        self.assertIn("no-such-force-field.xml", str(caught.exception))
        # An error of ACT, and an error of the operating system for those who catch that
        self.assertIsInstance(caught.exception, act.ActError)
        self.assertIsInstance(caught.exception, OSError)
        broken = self.dir / "broken.xml"
        broken.write_text("sigma 0.3\n")
        with self.assertRaises(act.ActFileError) as caught:
            act.ForceField.read(broken)
        self.assertIn("broken.xml", str(caught.exception))
        # Half of a force field
        whole = (self.ffdir / "ACM-g1.xml").read_text(errors="replace")
        half = self.dir / "half.xml"
        half.write_text(whole[:len(whole) // 2])
        with self.assertRaises(act.ActFileError):
            act.ForceField.read(half)
        # Still alive, and the next force field is read
        self.assertGreater(len(act.ForceField.read(self.ffdir / "ACM-g1.xml").particle_types), 0)

    def test_content_that_is_not_accepted(self):
        """!
        @brief A force field with a value that ACT does not know raises an
        error about the input, with what is wrong.
        """
        text = (self.ffdir / "ACM-g1.xml").read_text(errors="replace")
        self.assertIn('mutability="Bounded"', text)
        wrong = self.dir / "free.xml"
        wrong.write_text(text.replace('mutability="Bounded"', 'mutability="Free"', 1))
        with self.assertRaises(act.ActInputError) as caught:
            act.ForceField.read(wrong)
        self.assertIn("Free", str(caught.exception))
        self.assertIsInstance(caught.exception, ValueError)
        self.assertIsInstance(caught.exception, act.ActError)

    def test_write_a_force_field(self):
        """!
        @brief A force field is written to a file that reads back the same,
        and a place that cannot be written raises an error.
        """
        ff = act.ForceField.read(self.ffdir / "ACS-pg1.xml")
        out = self.dir / "written.xml"
        ff.write(out)
        self.assertGreater(out.stat().st_size, 1000)
        again = act.ForceField.read(out)
        self.assertEqual(ff.particle_types, again.particle_types)
        self.assertEqual(ff.interactions, again.interactions)
        # Compressed on request
        packed = self.dir / "packed.xml"
        ff.write(packed, compress=True)
        self.assertLess(packed.stat().st_size, out.stat().st_size)
        self.assertEqual(ff.particle_types, act.ForceField.read(packed).particle_types)
        with self.assertRaises(act.ActFileError) as caught:
            ff.write(self.dir / "no" / "such" / "directory" / "ff.xml")
        self.assertIn("ff.xml", str(caught.exception))

    def test_parameters_of_an_interaction(self):
        """!
        @brief The parameters of an interaction, each with its value, its
        bounds, its unit and whether training changes it.
        """
        ff = act.ForceField.read(self.ffdir / "ACS-pg1.xml")
        vdw = ff.parameters("VANDERWAALS")
        self.assertGreater(len(vdw), 100)
        self.assertEqual({"VANDERWAALS"}, {p.interaction for p in vdw})
        argon = [p for p in vdw if p.identifier == "Ar"]
        self.assertEqual(["epsilon", "gamma", "sigma"], [p.name for p in argon])
        sigma = ff.parameter("VANDERWAALS", "Ar", "sigma")
        self.assertEqual(("VANDERWAALS", "Ar", "sigma", 0.275, 0.15, 0.4, "nm", "Bounded"),
                         (sigma.interaction, sigma.identifier, sigma.name, sigma.value,
                          sigma.minimum, sigma.maximum, sigma.unit, sigma.mutability))
        self.assertIsInstance(sigma.ntrain, int)
        self.assertIsInstance(sigma.uncertainty, float)
        self.assertTrue(sigma.nonnegative)
        self.assertIn("sigma", repr(sigma))
        self.assertIn("0.275", repr(sigma))
        # A parameter is what the force field had when it was asked: it
        # cannot be changed by itself
        with self.assertRaises(AttributeError):
            sigma.value = 0.3
        # Bonds are between two types
        bonds = ff.parameters("BONDS")
        self.assertIn("br_b~br_b", {p.identifier for p in bonds})
        self.assertEqual(400, ff.parameter("BONDS", "br_b~br_b", "De").value)
        # What is not there
        for wrong, word in ((("NOSUCHTHING", "Ar", "sigma"), "NOSUCHTHING"),
                            (("VANDERWAALS", "Xx", "sigma"), "Xx"),
                            (("VANDERWAALS", "Ar", "kt"), "kt")):
            with self.assertRaises(act.ActInputError, msg=repr(wrong)) as caught:
                ff.parameter(*wrong)
            self.assertIn(word, str(caught.exception))
        with self.assertRaises(act.ActInputError):
            ff.parameters("NOSUCHTHING")
        # An interaction that ACT knows but this force field does not have
        plain = act.ForceField.read(self.ffdir / "ACS-g1.xml")
        self.assertNotIn("VSITE3OUT", plain.interactions)
        self.assertEqual([], plain.parameters("VSITE3OUT"))

    def test_particle_types(self):
        """!
        @brief A particle type with what it is, its options and its parameters.
        """
        ff = act.ForceField.read(self.ffdir / "ACS-pg1.xml")
        shell = ff.particle_type("c3_s")
        self.assertEqual(("c3_s", "Shell"), (shell.name, shell.kind))
        self.assertEqual("X", shell.options["element"])
        self.assertEqual(["charge", "mass", "ref_enthalpy", "zeta"],
                         [p.name for p in shell.parameters])
        charge = [p for p in shell.parameters if p.name == "charge"][0]
        self.assertEqual(("", "c3_s", -1.2, -2.0, -0.4, "e", "Bounded"),
                         (charge.interaction, charge.identifier, charge.value, charge.minimum,
                          charge.maximum, charge.unit, charge.mutability))
        self.assertEqual("Atom", ff.particle_type("c3").kind)
        self.assertEqual("C", ff.particle_type("c3").options["element"])
        self.assertIn("c3_s", repr(shell))
        with self.assertRaises(act.ActInputError) as caught:
            ff.particle_type("nosuchtype")
        self.assertIn("nosuchtype", str(caught.exception))

    def test_change_parameters(self):
        """!
        @brief A value, the bounds and the mutability of a parameter are
        changed, checked first, and are in the file that is written.
        """
        ff = act.ForceField.read(self.ffdir / "ACS-pg1.xml")
        ff.set_parameter("VANDERWAALS", "Ar", "sigma", value=0.3)
        self.assertEqual(0.3, ff.parameter("VANDERWAALS", "Ar", "sigma").value)
        ff.set_parameter("VANDERWAALS", "Ar", "sigma", minimum=0.2, maximum=0.5,
                         mutability="Fixed")
        sigma = ff.parameter("VANDERWAALS", "Ar", "sigma")
        self.assertEqual((0.3, 0.2, 0.5, "Fixed"),
                         (sigma.value, sigma.minimum, sigma.maximum, sigma.mutability))
        # A fixed parameter is changed when asked: fixed is for training
        ff.set_parameter("VANDERWAALS", "Ar", "sigma", value=0.25)
        self.assertEqual(0.25, ff.parameter("VANDERWAALS", "Ar", "sigma").value)
        # The value together with bounds that have it
        ff.set_parameter("VANDERWAALS", "Ar", "epsilon", value=2.0, maximum=3.0)
        epsilon = ff.parameter("VANDERWAALS", "Ar", "epsilon")
        self.assertEqual((2.0, 0.03, 3.0), (epsilon.value, epsilon.minimum, epsilon.maximum))
        # Of a particle type
        ff.set_particle_parameter("c3_s", "charge", value=-1.5, mutability="Fixed")
        charge = [p for p in ff.particle_type("c3_s").parameters if p.name == "charge"][0]
        self.assertEqual((-1.5, "Fixed"), (charge.value, charge.mutability))
        # The file has what was set, and nothing else changed
        out = self.dir / "changed.xml"
        ff.write(out)
        again = act.ForceField.read(out)
        sigma = again.parameter("VANDERWAALS", "Ar", "sigma")
        self.assertEqual((0.25, 0.2, 0.5, "Fixed"),
                         (sigma.value, sigma.minimum, sigma.maximum, sigma.mutability))
        self.assertEqual(2.0, again.parameter("VANDERWAALS", "Ar", "epsilon").value)
        self.assertEqual(-1.5, [p for p in again.particle_type("c3_s").parameters
                                if p.name == "charge"][0].value)
        before = act.ForceField.read(self.ffdir / "ACS-pg1.xml")
        self.assertEqual(before.parameter("VANDERWAALS", "Ar", "gamma").value,
                         again.parameter("VANDERWAALS", "Ar", "gamma").value)
        self.assertEqual(len(before.parameters("BONDS")), len(again.parameters("BONDS")))

    def test_changes_that_are_refused(self):
        """!
        @brief A change that would leave a parameter in a state that makes
        no sense is refused, and the parameter stays as it was.
        """
        ff = act.ForceField.read(self.ffdir / "ACS-pg1.xml")

        def sigma():
            p = ff.parameter("VANDERWAALS", "Ar", "sigma")
            return (p.value, p.minimum, p.maximum, p.mutability)
        before = sigma()
        wrong = [
            dict(value=0.5),                    # above the maximum
            dict(value=0.1),                    # below the minimum
            dict(minimum=0.3),                  # the value would be below it
            dict(maximum=0.2),                  # the value would be above it
            dict(minimum=0.5, maximum=0.4),     # the minimum above the maximum
            dict(value=-0.1, minimum=-0.2),     # negative, for a parameter that cannot be
            dict(mutability="Free"),            # not a mutability
            dict(value=float("nan")),
            dict(maximum=float("inf")),
            dict(),                             # nothing to change
        ]
        for change in wrong:
            with self.assertRaises(act.ActInputError, msg=repr(change)):
                ff.set_parameter("VANDERWAALS", "Ar", "sigma", **change)
            self.assertEqual(before, sigma(), change)
        for where in (("NOSUCHTHING", "Ar", "sigma"), ("VANDERWAALS", "Xx", "sigma"),
                      ("VANDERWAALS", "Ar", "kt")):
            with self.assertRaises(act.ActInputError, msg=repr(where)):
                ff.set_parameter(*where, value=0.3)
        with self.assertRaises(act.ActInputError):
            ff.set_particle_parameter("nosuchtype", "charge", value=1.0)
        with self.assertRaises(act.ActInputError):
            ff.set_particle_parameter("c3_s", "nosuchparameter", value=1.0)
        # The values by name only, and numbers only
        with self.assertRaises(TypeError):
            ff.set_parameter("VANDERWAALS", "Ar", "sigma", 0.3)  # pylint: disable=too-many-function-args
        with self.assertRaises(TypeError):
            ff.set_parameter("VANDERWAALS", "Ar", "sigma", value="large")

    def test_validate(self):
        """!
        @brief The problems with the parameters of a force field, as
        edit_ff -validate lists them.
        """
        ff = act.ForceField.read(self.ffdir / "ACS-pg1.xml")
        self.assertEqual([], ff.validate())
        # A value outside its bounds in the file is set to the nearest
        # bound when the file is read, which is told
        text = (self.ffdir / "ACM-g1.xml").read_text(errors="replace")
        marker = 'maximum="'
        at = text.index(marker, text.index('<parameter type="'))
        end = text.index('"', at + len(marker))
        outside = self.dir / "outside.xml"
        outside.write_text(text[:at + len(marker)] + "-1e9" + text[end:])
        try:
            issues = act.ForceField.read(outside).validate()
        except act.ActInputError:
            # A force field that cannot be read at all is not validated
            return
        self.assertGreater(len(issues), 0)
        issue = issues[0]
        self.assertEqual({"kind", "where", "identifier", "parameter", "text", "repaired"},
                         set(issue))
        self.assertIsInstance(issue["repaired"], bool)
        self.assertTrue(issue["text"])

    def test_wrong_arguments(self):
        """!
        @brief What is no file name is refused by Python, not by a crash.
        """
        for wrong in (None, 3, ["ff.xml"]):
            with self.assertRaises(TypeError):
                act.ForceField.read(wrong)
        with self.assertRaises(TypeError):
            act.ForceField.read()  # pylint: disable=no-value-for-parameter


@unittest.skipUnless(act.classes_available(), "ACT was built without its Python classes")
@unittest.skipUnless(os.path.isdir(TESTDATA), "ACT_TESTDATA is not set")
class MoleculeTest(unittest.TestCase):
    """!
    @brief A molecule in Python: what it is made of and its energy.
    """

    def setUp(self):
        """!
        @brief Read the force fields of the tests.
        """
        self.data = Path(TESTDATA)
        self.ffname = self.data / "forcefields" / "ACS-pg1.xml"
        self.ff = act.ForceField.read(self.ffname)
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self._tmp.name)

    def tearDown(self):
        """!
        @brief Remove the files.
        """
        self._tmp.cleanup()

    def test_molecule_from_smiles(self):
        """!
        @brief A molecule is made from a SMILES string, with a structure, a
        topology and charges for the force field.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        self.assertEqual("C3H6O", mol.formula)
        self.assertEqual(0, mol.charge)
        self.assertAlmostEqual(58.08, mol.mass, places=2)
        # Ten atoms, each with a shell in this polarizable force field
        atoms = mol.atoms
        self.assertEqual(20, len(atoms))
        self.assertEqual(10, len([a for a in atoms if a["particletype"] == "Atom"]))
        self.assertEqual(10, len([a for a in atoms if a["particletype"] == "Shell"]))
        self.assertEqual(["C", "C", "C", "H", "H", "H", "H", "H", "H", "O"],
                         sorted(a["element"] for a in atoms if a["particletype"] == "Atom"))
        self.assertTrue(all(a["type"] in self.ff.particle_types for a in atoms))
        # The charges are those of the particles, and add up to the charge
        charges = mol.charges
        self.assertEqual([a["charge"] for a in atoms], plain(charges))
        self.assertAlmostEqual(0.0, sum(charges), places=6)
        self.assertGreater(max(abs(q) for q in charges), 0.1)
        # Nine bonds, of which one is double, between atoms
        bonds = mol.bonds
        self.assertEqual(9, len(bonds))
        self.assertEqual([2], [b["order"] for b in bonds if b["order"] != 1])
        for bond in bonds:
            self.assertTrue(all(atoms[i]["particletype"] == "Atom" for i in bond["atoms"]))
        # A coordinate for each particle, in nm, of a molecule of that size
        coordinates = mol.coordinates
        self.assertEqual(20, len(coordinates))
        self.assertTrue(all(len(x) == 3 for x in coordinates))
        spread = max(x[0] for x in coordinates) - min(x[0] for x in coordinates)
        self.assertGreater(spread, 0.05)
        self.assertLess(spread, 1.0)
        self.assertIn("C3H6O", repr(mol))

    def test_energy_and_forces(self):
        """!
        @brief The energy by term in kJ/mol and the force on each particle.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        energy = mol.energy()
        self.assertIn("EPOT", energy)
        self.assertIn("BONDS", energy)
        self.assertIn("COULOMB", energy)
        self.assertTrue(all(isinstance(value, float) for value in energy.values()))
        # The terms add up to the potential energy
        parts = sum(value for name, value in energy.items()
                    if name not in ("EPOT", "ALLELEC", "EXCHIND", "EXCHDISP"))
        self.assertLess(energy["EPOT"], 0)
        self.assertAlmostEqual(energy["EPOT"], parts, delta=1e-3 * abs(energy["EPOT"]))
        forces = mol.forces()
        self.assertEqual(len(mol.atoms), len(forces))
        self.assertTrue(all(len(f) == 3 for f in forces))
        # No net force on the molecule
        for dim in range(3):
            self.assertAlmostEqual(0.0, sum(f[dim] for f in forces), delta=1e-2)
        # The same again: asking changes nothing
        self.assertEqual(energy, mol.energy())

    @unittest.skipUnless(os.environ.get("ACT_ALEXANDRIA"), "ACT_ALEXANDRIA is not set")
    def test_same_as_the_program(self):
        """!
        @brief The class gives what the alexandria program gives for the
        same molecule and force field: the particles with their charges,
        the structure and the energies.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        with act.commands.gentop(ff=self.ffname, smiles="CC(=O)C") as res:
            told = res.data["molecules"][0]
        self.assertEqual(told["molecule"]["formula"], mol.formula)
        self.assertEqual([p["type"] for p in told["molecule"]["particles"]],
                         [a["type"] for a in mol.atoms])
        for theirs, ours in zip(told["molecule"]["particles"], mol.charges):
            self.assertAlmostEqual(theirs["charge"], ours, places=5)
        self.assertEqual(told["molecule"]["bonds"], mol.bonds)
        for theirs, ours in zip(told["coordinates"]["positions"], mol.coordinates):
            for dim in range(3):
                self.assertAlmostEqual(theirs[dim], ours[dim], places=5)
        with act.commands.minimize(ff=self.ffname, smiles="CC(=O)C") as res:
            before = res.data["Energies before"]
        energy = mol.energy()
        self.assertEqual(set(before), set(energy))
        for name, value in before.items():
            self.assertAlmostEqual(value["value"], energy[name],
                                   delta=1e-3 + 1e-5 * abs(value["value"]), msg=name)

    def test_minimize(self):
        """!
        @brief The energy of a molecule is minimized: the molecule gets the
        structure of the minimum, and the result tells how it went.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        start = mol.coordinates
        before = mol.energy()["EPOT"]
        res = mol.minimize()
        self.assertTrue(res.converged)
        self.assertEqual("OK", res.status)
        self.assertAlmostEqual(before, res.energy_before["EPOT"], places=6)
        self.assertLess(res.energy_after["EPOT"], before)
        self.assertEqual(set(res.energy_before), set(res.energy_after))
        # It moved a little, and not far
        self.assertGreater(res.rmsd, 1e-5)
        self.assertLess(res.rmsd, 0.05)
        # The molecule has the new structure and its energy
        self.assertEqual(plain(res.coordinates), plain(mol.coordinates))
        self.assertNotEqual(plain(start), plain(mol.coordinates))
        self.assertAlmostEqual(res.energy_after["EPOT"], mol.energy()["EPOT"],
                               delta=1e-4 * abs(before))
        # The forces on the atoms are gone
        atoms = [i for i, a in enumerate(mol.atoms) if a["particletype"] == "Atom"]
        forces = mol.forces()
        self.assertLess(max(abs(forces[i][dim]) for i in atoms for dim in range(3)), 1.0)
        self.assertIn("converged", repr(res))
        # From the minimum there is nowhere to go
        again = mol.minimize()
        self.assertLess(again.rmsd, 1e-3)
        self.assertAlmostEqual(res.energy_after["EPOT"], again.energy_after["EPOT"],
                               delta=1e-4 * abs(before))

    def test_minimize_without_moving_the_molecule(self):
        """!
        @brief On request the molecule keeps its structure and only the
        result has the minimum.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        start = mol.coordinates
        res = mol.minimize(update=False)
        self.assertTrue(res.converged)
        self.assertEqual(plain(start), plain(mol.coordinates))
        self.assertNotEqual(plain(start), plain(res.coordinates))
        self.assertEqual(len(start), len(res.coordinates))
        self.assertAlmostEqual(res.energy_before["EPOT"], mol.energy()["EPOT"], places=6)

    @unittest.skipUnless(os.environ.get("ACT_ALEXANDRIA"), "ACT_ALEXANDRIA is not set")
    def test_minimize_same_as_the_program(self):
        """!
        @brief The minimization gives what alexandria minimize gives.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        res = mol.minimize()
        with act.commands.minimize(ff=self.ffname, smiles="CC(=O)C") as told:
            data = told.data
        self.assertEqual(data["converged"], res.converged)
        self.assertEqual(data["status"], res.status)
        self.assertAlmostEqual(data["rmsd"]["value"], res.rmsd, places=6)
        # The program may be another build than the module, as in the conda
        # packages. In single precision the two then find minima that differ
        # by the precision of the largest term, here the bonds, and the
        # smaller terms differ by as much.
        largest = max(abs(value["value"]) for value in data["Energies after"].values())
        for name, value in data["Energies after"].items():
            self.assertAlmostEqual(value["value"], res.energy_after[name],
                                   delta=1e-3 + 1e-5 * abs(value["value"]) + 2e-6 * largest,
                                   msg=name)

    def test_minimize_settings(self):
        """!
        @brief The algorithm, the number of iterations, the tolerance and
        atoms that stay where they are.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        # Two steps of steepest descent do not get there, which is told, not raised
        res = mol.minimize(algorithm="Steep", max_iter=2, update=False)
        self.assertFalse(res.converged)
        self.assertEqual("Too many steps", res.status)
        self.assertLessEqual(res.energy_after["EPOT"], res.energy_before["EPOT"])
        for algorithm in ("LBFGS", "Newton"):
            self.assertTrue(mol.minimize(algorithm=algorithm, update=False).converged, algorithm)
        # A looser tolerance is reached as well
        self.assertTrue(mol.minimize(tolerance=1e-3, update=False).converged)
        # Frozen atoms keep their place: the two that make the C=O bond
        double = [b for b in mol.bonds if b["order"] == 2][0]["atoms"]
        start = mol.coordinates
        res = mol.minimize(freeze=double, update=False)
        free = act.Molecule.from_smiles("CC(=O)C", self.ff).minimize(update=False)
        self.assertGreater(res.energy_after["EPOT"], free.energy_after["EPOT"] - 1e-6)

        def distance(coordinates, i, j):
            return sum((coordinates[i][d] - coordinates[j][d]) ** 2 for d in range(3)) ** 0.5
        self.assertAlmostEqual(distance(start, *double), distance(res.coordinates, *double),
                               places=7)
        self.assertNotAlmostEqual(distance(start, *double), distance(free.coordinates, *double),
                                  places=7)
        wrong = [dict(algorithm="Simplex"), dict(algorithm=""), dict(max_iter=0),
                 dict(max_iter=-5), dict(tolerance=0.0), dict(tolerance=-1e-6),
                 dict(tolerance=float("nan")), dict(freeze=[999]), dict(freeze=[-1]),
                 # A shell cannot be frozen, it follows its atom
                 dict(freeze=[i for i, a in enumerate(mol.atoms)
                              if a["particletype"] == "Shell"][:1])]
        for settings in wrong:
            with self.assertRaises(act.ActInputError, msg=repr(settings)):
                mol.minimize(**settings)
        self.assertEqual(plain(start), plain(mol.coordinates))
        with self.assertRaises(TypeError):
            mol.minimize("LBFGS")  # pylint: disable=too-many-function-args
        with self.assertRaises(TypeError):
            mol.minimize(max_iter="many")

    def test_settings_are_checked_as_the_programs_check_them(self):
        """!
        @brief Settings that a program refuses are refused with its message,
        which names its option; the rules of Python come on top of them.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        told = {
            "-toler": lambda: mol.minimize(tolerance=0.0),
            "Unknown minimization algorithm": lambda: mol.minimize(algorithm="Simplex"),
            "-maxiter": lambda: mol.minimize(max_iter=0),
            "There is no particle": lambda: mol.minimize(freeze=[999]),
            "-deltat": lambda: mol.simulate(steps=10, timestep=0.0),
            "-temp": lambda: mol.simulate(steps=10, temperature=-1.0),
            "-nstener": lambda: mol.simulate(steps=10, energy_interval=-1),
            "-nstxout": lambda: mol.simulate(steps=10, frame_interval=-1),
            "-nsteps": lambda: mol.simulate(steps=0),
            "-mDrude": lambda: mol.write_openmm(self.dir / "o.xml", drude_mass=-1.0),
            "-ntrain": lambda: mol.write_openmm(self.dir / "o.xml", min_training=-1),
        }
        for words, call in told.items():
            with self.assertRaises(act.ActInputError, msg=words) as caught:
                call()
            self.assertIn(words, str(caught.exception))
        with self.assertRaises(act.ActInputError) as caught:
            act.Molecule.from_smiles("O", self.ff).normal_modes(tolerance=-1.0)
        self.assertIn("-toler", str(caught.exception))

    def test_minimize_a_dimer(self):
        """!
        @brief The minimization of a dimer gives its interaction energy as
        well, that of the minimized structure; a molecule of one compound
        has none.
        """
        self.assertIsNone(act.Molecule.from_smiles("O", self.ff).minimize().interaction_energy)
        mol = act.Molecule.from_file(self.data / "molecules" / "water_dimer.pdb", self.ff)
        res = mol.minimize()
        self.assertTrue(res.converged)
        self.assertLess(res.interaction_energy["EPOT"], 0)
        # The same total as for the minimized structure that the molecule
        # has now. As alexandria minimize the minimization gives the
        # correction of the induction energy as a term of its own, which
        # interaction_energy adds to the induction energy, as simulate -sp.
        now = mol.interaction_energy()
        self.assertIn("INDUCTIONCORRECTION", res.interaction_energy)
        self.assertAlmostEqual(now["EPOT"], res.interaction_energy["EPOT"],
                               delta=1e-4 * max(1.0, abs(now["EPOT"])))

    @unittest.skipUnless(os.environ.get("ACT_ALEXANDRIA"), "ACT_ALEXANDRIA is not set")
    def test_minimize_a_dimer_same_as_the_program(self):
        """!
        @brief The interaction energy after a minimization is that of
        alexandria minimize.
        """
        dimer = self.data / "molecules" / "water_dimer.pdb"
        res = act.Molecule.from_file(dimer, self.ff).minimize()
        with act.commands.minimize(ff=self.ffname, f=dimer) as told:
            theirs = told.data["Interaction energies"]
        self.assertEqual(set(theirs), set(res.interaction_energy))
        for name, value in theirs.items():
            # The program writes four decimals; two builds may differ by
            # the precision of the minimum, see test_minimize_same_as_the_program
            self.assertAlmostEqual(float(value["value"]), res.interaction_energy[name],
                                   delta=2e-3 + 1e-4 * abs(float(value["value"])), msg=name)

    def test_normal_modes(self):
        """!
        @brief The vibrations of a molecule at its minimum: frequencies,
        intensities, the motion of the particles and the thermochemistry.
        """
        mol = act.Molecule.from_smiles("O", self.ff)
        start = mol.coordinates
        modes = mol.normal_modes()
        # Water has three vibrations: a bend and two stretches
        self.assertEqual(3, len(modes.frequencies))
        self.assertEqual(sorted(modes.frequencies), plain(modes.frequencies))
        self.assertGreater(modes.frequencies[0], 500)
        self.assertLess(modes.frequencies[0], 2500)
        self.assertGreater(modes.frequencies[1], 2500)
        self.assertLess(modes.frequencies[2], 5000)
        self.assertEqual(3, len(modes.intensities))
        self.assertTrue(all(value >= 0 for value in modes.intensities))
        # For each vibration how the particles move
        self.assertEqual(3, len(modes.modes))
        for k, mode in enumerate(modes.modes):
            self.assertEqual(modes.frequencies[k], mode["frequency"])
            self.assertEqual(modes.intensities[k], mode["intensity"])
            self.assertGreater(mode["reduced_mass"], 0.5)
            self.assertEqual(len(mol.atoms), len(mode["displacements"]))
            self.assertTrue(all(len(d) == 3 for d in mode["displacements"]))
        # The molecule was minimized first, and has that structure now
        self.assertTrue(modes.minimization.converged)
        self.assertNotEqual(plain(start), plain(mol.coordinates))
        thermo = modes.thermochemistry
        self.assertIn("Epot", thermo)
        self.assertGreater(thermo["T=0"]["Zero point energy"]["value"], 10)
        self.assertEqual("kJ/mol", thermo["T=0"]["Zero point energy"]["unit"])
        self.assertIn("3 modes", repr(modes))
        # At the minimum already: the same without minimizing again
        again = mol.normal_modes(minimize=False)
        self.assertIsNone(again.minimization)
        for first, second in zip(modes.frequencies, again.frequencies):
            self.assertAlmostEqual(first, second, delta=0.5)

    def test_spectrum(self):
        """!
        @brief The infrared spectrum: a line of a given width at each frequency.
        """
        modes = act.Molecule.from_smiles("O", self.ff).normal_modes()
        frequencies, intensities = modes.spectrum()
        self.assertEqual(len(frequencies), len(intensities))
        self.assertEqual(list(range(frequencies[0], frequencies[-1] + 1)), plain(frequencies))
        # From below the lowest to above the highest frequency
        self.assertLess(frequencies[0], modes.frequencies[0] - 100)
        self.assertGreater(frequencies[-1], modes.frequencies[-1] + 100)
        self.assertAlmostEqual(1.0, sum(intensities), places=6)
        # A maximum at the strongest vibration
        strongest = modes.frequencies[int(modes.intensities.argmax())]
        top = frequencies[int(intensities.argmax())]
        self.assertLess(abs(top - strongest), 2)
        # Narrow lines are higher
        narrow = modes.spectrum(linewidth=4.0)
        self.assertGreater(max(narrow[1]), max(intensities))
        for wrong in (0.0, -1.0, float("nan")):
            with self.assertRaises(act.ActInputError):
                modes.spectrum(linewidth=wrong)

    @unittest.skipUnless(os.environ.get("ACT_ALEXANDRIA"), "ACT_ALEXANDRIA is not set")
    def test_normal_modes_same_as_the_program(self):
        """!
        @brief The frequencies, intensities and thermochemistry are those of
        alexandria nma, and the spectrum is the one it writes.
        """
        modes = act.Molecule.from_smiles("CC(=O)C", self.ff).normal_modes()
        with act.commands.nma(ff=self.ffname, smiles="CC(=O)C", ir=True) as told:
            data = told.data
            spectrum = [line.split() for line in told.file("ir").read_text().splitlines()
                        if line and line[0] not in "#@&"]
        theirs = [mode["frequency"]["value"] for mode in data["modes"]]
        self.assertEqual(len(theirs), len(modes.frequencies))
        for first, second in zip(theirs, modes.frequencies):
            self.assertAlmostEqual(first, second, delta=0.05)
        for mode, intensity in zip(data["modes"], modes.intensities):
            self.assertAlmostEqual(mode["intensity"]["value"], intensity,
                                   delta=1e-6 + 1e-3 * abs(intensity))
        self.assertAlmostEqual(
            data["Thermochemistry"]["Alexandria"]["T=0"]["Zero point energy"]["value"],
            modes.thermochemistry["T=0"]["Zero point energy"]["value"], delta=0.01)
        frequencies, intensities = modes.spectrum()
        self.assertEqual([int(row[0]) for row in spectrum], plain(frequencies))
        for row, intensity in zip(spectrum, intensities):
            self.assertAlmostEqual(float(row[1]), intensity, delta=1e-7 + 1e-4 * intensity)

    def test_normal_modes_need_a_minimum(self):
        """!
        @brief Vibrations of a structure that is no minimum make no sense:
        a minimization that does not get there raises an error.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        with self.assertRaises(act.ActError) as caught:
            mol.normal_modes(max_iter=1, algorithm="Steep")
        self.assertIn("minim", str(caught.exception))
        # The molecule is still good for a minimization that does
        self.assertEqual(24, len(mol.normal_modes().frequencies))

    def test_simulate(self):
        """!
        @brief A molecular dynamics simulation gives the energies and the
        structure along the way, and leaves the molecule as it is.
        """
        mol = act.Molecule.from_smiles("O", self.ff)
        mol.minimize()
        start = mol.coordinates
        traj = mol.simulate(steps=200, temperature=300, seed=7, energy_interval=10,
                            frame_interval=20)
        # The energies every ten steps, the structure every twenty
        self.assertEqual(20, len(traj.energy_times))
        self.assertEqual(0.0, traj.energy_times[0])
        self.assertAlmostEqual(10 * 0.0002, traj.energy_times[1], places=9)
        self.assertEqual(list(range(0, 200, 10)), plain(traj.energy_steps))
        for name in ("EPOT", "EKIN", "ETOT", "BONDS"):
            self.assertEqual(20, len(traj.energies[name]), name)
        self.assertEqual(20, len(traj.temperature))
        # The velocities are drawn for the temperature that was asked for
        self.assertAlmostEqual(300, traj.temperature[0], delta=1)
        for k in range(20):
            self.assertAlmostEqual(traj.energies["ETOT"][k],
                                   traj.energies["EPOT"][k] + traj.energies["EKIN"][k], places=6)
        # The total energy stays, while the potential energy moves
        total, potential = traj.energies["ETOT"], traj.energies["EPOT"]
        self.assertGreater(max(potential) - min(potential), 0.1)
        self.assertLess(max(total) - min(total), 0.2 * (max(potential) - min(potential)))
        self.assertEqual(10, len(traj.frames))
        self.assertEqual(list(range(0, 200, 20)), plain(traj.frame_steps))
        self.assertAlmostEqual(20 * 0.0002, traj.times[1], places=9)
        for frame in traj.frames:
            self.assertEqual(len(mol.atoms), len(frame))
        # It starts where the molecule is, and moves
        for theirs, ours in zip(traj.frames[0], start):
            for dim in range(3):
                self.assertAlmostEqual(ours[dim], theirs[dim], places=6)
        self.assertNotEqual(plain(traj.frames[0]), plain(traj.frames[-1]))
        self.assertEqual(plain(start), plain(mol.coordinates))
        self.assertIn("200 steps", repr(traj))

    def test_simulation_is_reproducible(self):
        """!
        @brief The same seed gives the same simulation, another seed another one.
        """
        mol = act.Molecule.from_smiles("O", self.ff)
        mol.minimize()
        first = mol.simulate(steps=50, temperature=300, seed=7, energy_interval=5)
        second = mol.simulate(steps=50, temperature=300, seed=7, energy_interval=5)
        other = mol.simulate(steps=50, temperature=300, seed=8, energy_interval=5)
        self.assertEqual(plain(first.energies["EPOT"]), plain(second.energies["EPOT"]))
        self.assertNotEqual(plain(first.energies["EPOT"]), plain(other.energies["EPOT"]))
        # Without a temperature the atoms start at rest, at the minimum
        still = mol.simulate(steps=20, energy_interval=5)
        self.assertLess(max(still.energies["EKIN"]), 1e-3)

    def test_follow_and_stop_a_simulation(self):
        """!
        @brief A function is called with the energies while the simulation
        runs, and stops it by returning False.
        """
        mol = act.Molecule.from_smiles("O", self.ff)
        seen = []
        traj = mol.simulate(steps=100, temperature=300, seed=7, energy_interval=10,
                            on_energies=seen.append)
        self.assertEqual(10, len(seen))
        self.assertEqual(list(range(0, 100, 10)), [item["step"] for item in seen])
        self.assertEqual(plain(traj.energies["EPOT"]), [item["energies"]["EPOT"] for item in seen])
        self.assertEqual(plain(traj.temperature), [item["temperature"] for item in seen])
        self.assertAlmostEqual(seen[1]["time"], traj.energy_times[1], places=12)

        def three(item):
            return item["step"] < 20
        stopped = mol.simulate(steps=100000, temperature=300, seed=7, energy_interval=10,
                               on_energies=three)
        self.assertEqual([0, 10, 20], plain(stopped.energy_steps))
        # An error in the function ends the simulation and is raised

        def broken(item):
            raise KeyError("no such thing")
        with self.assertRaises(KeyError):
            mol.simulate(steps=100000, temperature=300, on_energies=broken)
        # The molecule is still good
        self.assertEqual(5, len(mol.simulate(steps=5, energy_interval=1).energy_steps))

    @unittest.skipUnless(signal.getsignal(signal.SIGINT) is signal.default_int_handler,
                         "Python does not handle the interrupt signal here")
    def test_interrupt_ends_a_simulation(self):
        """!
        @brief Ctrl-C ends a simulation soon, also for a molecule that was
        made from a SMILES string, which used to take the signal away from
        Python for good.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        mol.minimize()
        for follow in (None, lambda item: None):
            # Half a minute of simulation, if nothing ends it
            timer = threading.Timer(0.3, lambda: signal.raise_signal(signal.SIGINT))
            began = time.monotonic()
            timer.start()
            try:
                with self.assertRaises(KeyboardInterrupt):
                    mol.simulate(steps=1000000, temperature=100, seed=3, energy_interval=10,
                                 frame_interval=0, on_energies=follow)
            finally:
                timer.cancel()
                timer.join()
            self.assertLess(time.monotonic() - began, 10)
        # And the molecule is as good as before
        self.assertEqual(4, len(mol.simulate(steps=4).energy_steps))

    @unittest.skipUnless(os.environ.get("ACT_ALEXANDRIA"), "ACT_ALEXANDRIA is not set")
    def test_simulation_same_as_the_program(self):
        """!
        @brief The simulation is that of alexandria simulate with the same seed.
        """
        told = []
        with act.commands.simulate(ff=self.ffname, smiles="O", nsteps=60, temp=300, seed=1993,
                                   nstener=10, nstxout=0, on_progress=told.append):
            pass
        mol = act.Molecule.from_smiles("O", self.ff)
        traj = mol.simulate(steps=60, temperature=300, seed=1993, energy_interval=10,
                            frame_interval=0)
        self.assertEqual([item["step"] for item in told], plain(traj.energy_steps))
        self.assertEqual(0, len(traj.frames))
        for item, potential, temperature in zip(told, traj.energies["EPOT"], traj.temperature):
            self.assertAlmostEqual(item["energies"]["EPOT"], potential,
                                   delta=1e-6 * abs(potential) + 1e-6)
            self.assertAlmostEqual(item["temperature"]["value"], temperature, delta=1e-4)

    def test_simulation_settings(self):
        """!
        @brief Settings that make no sense are refused.
        """
        mol = act.Molecule.from_smiles("O", self.ff)
        wrong = [dict(steps=0), dict(steps=-10), dict(steps=10, timestep=0.0),
                 dict(steps=10, timestep=-0.001), dict(steps=10, timestep=float("nan")),
                 dict(steps=10, temperature=-1.0), dict(steps=10, energy_interval=-1),
                 dict(steps=10, frame_interval=-5), dict(steps=10, seed=-3)]
        for settings in wrong:
            with self.assertRaises(act.ActInputError, msg=repr(settings)):
                mol.simulate(**settings)
        with self.assertRaises(TypeError):
            mol.simulate()  # pylint: disable=missing-kwoa
        with self.assertRaises(TypeError):
            mol.simulate(100)  # pylint: disable=too-many-function-args,missing-kwoa
        with self.assertRaises(TypeError):
            mol.simulate(steps=10, on_energies="print")
        # A larger time step is taken as it is
        self.assertAlmostEqual(0.0005, mol.simulate(steps=3, timestep=0.0005).energy_times[1],
                               places=9)

    def test_write_openmm(self):
        """!
        @brief The force field file with which OpenMM simulates the molecule,
        and a structure file that goes with it.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        xml = self.dir / "acetone-openmm.xml"
        pdb = self.dir / "acetone.pdb"
        mol.write_openmm(xml, structure=pdb)
        text = xml.read_text()
        self.assertIn("<ForceField>", text)
        self.assertIn("<AtomTypes>", text)
        self.assertIn("<Residue ", text)
        atoms = [line for line in pdb.read_text().splitlines() if line.startswith("ATOM")]
        # The atoms without their shells, which OpenMM adds itself
        self.assertEqual(10, len(atoms))
        # The residue of the structure is the one of the force field file
        residue = atoms[0][17:20].strip()
        self.assertIn('<Residue name="%s"' % residue, text)
        # With the shells on request, and without a structure file
        shells = self.dir / "shells.pdb"
        mol.write_openmm(self.dir / "again.xml", structure=shells, shells=True)
        self.assertEqual(20, len([line for line in shells.read_text().splitlines()
                                  if line.startswith(("ATOM", "HETATM"))]))
        mol.write_openmm(self.dir / "only.xml")
        self.assertEqual(text, (self.dir / "only.xml").read_text())
        # The structure is the one the molecule has now
        mol.minimize()
        moved = self.dir / "minimized.pdb"
        mol.write_openmm(self.dir / "minimized.xml", structure=moved)
        self.assertNotEqual(pdb.read_text(), moved.read_text())
        for wrong in (dict(drude_mass=-0.1), dict(drude_mass=float("nan")),
                      dict(min_training=-1)):
            with self.assertRaises(act.ActInputError, msg=repr(wrong)):
                mol.write_openmm(self.dir / "wrong.xml", **wrong)
        with self.assertRaises(act.ActFileError):
            mol.write_openmm(self.dir / "no" / "such" / "directory" / "o.xml")
        with self.assertRaises(act.ActFileError):
            mol.write_openmm(self.dir / "o.xml",
                             structure=self.dir / "no" / "such" / "directory" / "c.pdb")
        with self.assertRaises(TypeError):
            mol.write_openmm()  # pylint: disable=no-value-for-parameter

    def test_write_openmm_residue_names(self):
        """!
        @brief The residue names are given as gentop -resnames takes them,
        and are the same in both files.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        xml = self.dir / "named.xml"
        pdb = self.dir / "named.pdb"
        mol.write_openmm(xml, structure=pdb, residue_names="ACT")
        atoms = [line for line in pdb.read_text().splitlines() if line.startswith("ATOM")]
        self.assertEqual({"ACT"}, {line[17:20].strip() for line in atoms})
        self.assertIn('<Residue name="ACT"', xml.read_text())
        # A name that OpenMM has already, and one that is too long
        for wrong in ("HOH", "TOOLONG"):
            with self.assertRaises(act.ActInputError, msg=wrong):
                mol.write_openmm(self.dir / "wrong.xml", residue_names=wrong)

    @unittest.skipUnless(os.environ.get("ACT_ALEXANDRIA"), "ACT_ALEXANDRIA is not set")
    def test_write_openmm_same_as_the_program(self):
        """!
        @brief The files are those that alexandria gentop -openmm writes.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        xml = self.dir / "o.xml"
        pdb = self.dir / "c.pdb"
        mol.write_openmm(xml, structure=pdb)
        with act.commands.gentop(ff=self.ffname, smiles="CC(=O)C", openmm=True, c=True) as told:
            theirs = ElementTree.parse(told.file("openmm")).getroot()
            self.assertEqual(told.file("c").read_text(), pdb.read_text())
        ours = ElementTree.parse(xml).getroot()

        def same(first, second, where):
            """Compare two elements: the same tags and attributes, numbers as numbers.

            The program has the mass of a shell in single precision, so that
            masses differ in the eighth digit.
            """
            self.assertEqual(first.tag, second.tag, where)
            self.assertEqual(sorted(first.attrib), sorted(second.attrib), where)
            for key, value in first.attrib.items():
                try:
                    self.assertAlmostEqual(float(value), float(second.attrib[key]),
                                           delta=1e-6 * max(1.0, abs(float(value))),
                                           msg="%s %s" % (where, key))
                except ValueError:
                    self.assertEqual(value, second.attrib[key], "%s %s" % (where, key))
            self.assertEqual(len(first), len(second), where)
            for one, other in zip(first, second):
                same(one, other, where + "/" + one.tag)
        same(theirs, ours, theirs.tag)

    def test_charges(self):
        """!
        @brief The charges are made again, for the structure the molecule has
        now, or set by hand.
        """
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        first = mol.charges
        # For the same structure the same charges
        mol.generate_charges()
        for before, after in zip(first, mol.charges):
            self.assertAlmostEqual(before, after, places=8)
        # The algorithm by name: the one of the force field, or its charges as they are
        mol.generate_charges("SQE")
        for before, after in zip(first, mol.charges):
            self.assertAlmostEqual(before, after, places=8)
        # The charges of the particle types, by the name the programs use
        fixed = act.Molecule.from_smiles("CC(=O)C", self.ff)
        fixed.generate_charges("None")
        self.assertEqual(len(first), len(fixed.charges))
        self.assertTrue(any(abs(a - b) > 1e-6 for a, b in zip(first, fixed.charges)))
        # In another structure the charges are other ones, with the same sum
        energy = mol.energy()["EPOT"]
        mol.minimize()
        mol.generate_charges()
        self.assertAlmostEqual(0.0, sum(mol.charges), places=6)
        self.assertGreater(max(abs(a - b) for a, b in zip(first, mol.charges)), 1e-6)
        self.assertNotEqual(energy, mol.energy()["EPOT"])
        # By hand: a charge for each particle
        mine = [0.0] * len(mol.atoms)
        mol.set_charges(mine)
        self.assertEqual(mine, plain(mol.charges))
        self.assertAlmostEqual(0.0, mol.energy()["COULOMB"], places=6)
        mine[0], mine[2] = 0.5, -0.5
        mol.set_charges(mine)
        self.assertEqual(mine, plain(mol.charges))
        self.assertNotAlmostEqual(0.0, mol.energy()["COULOMB"], places=3)
        self.assertEqual([a["charge"] for a in mol.atoms], plain(mol.charges))
        for wrong in ([0.0], mine + [0.1], [float("nan")] * len(mine)):
            with self.assertRaises(act.ActInputError):
                mol.set_charges(wrong)
        self.assertEqual(mine, plain(mol.charges))
        for wrong in ("Mulliken", "sqe", "", "ESP", "Read", "Custom"):
            with self.assertRaises(act.ActInputError, msg=wrong):
                mol.generate_charges(wrong)
        with self.assertRaises(TypeError):
            mol.set_charges("zero")
        # And back to those of the force field
        mol.generate_charges()
        self.assertGreater(max(abs(q) for q in mol.charges), 0.1)

    @unittest.skipUnless(os.environ.get("ACT_ALEXANDRIA"), "ACT_ALEXANDRIA is not set")
    def test_interaction_energy_same_as_the_program(self):
        """!
        @brief The interaction energy of a dimer is that of alexandria
        simulate -sp, also with shells.
        """
        dimer = str(self.data / "molecules" / "water_dimer.pdb")
        for ffname in ("ACS-g1.xml", "ACS-pg1.xml"):
            ffile = self.data / "forcefields" / ffname
            mol   = act.Molecule.from_file(dimer, act.ForceField.read(ffile))
            before = mol.energy()
            mine  = mol.interaction_energy()
            self.assertEqual(before, mol.energy(), ffname)
            with act.commands.simulate(ff=ffile, f=dimer, sp=True) as res:
                log = res.log
            theirs = {}
            for line in log.splitlines():
                if "Interaction energy" in line:
                    name, value = line.split("Interaction energy")[1].split(":")
                    theirs[name.strip()] = float(value)
            self.assertTrue(theirs, ffname)
            self.assertEqual(set(theirs), set(mine), ffname)
            for name, value in theirs.items():
                # The program writes six significant digits
                self.assertAlmostEqual(value, mine[name], delta=1e-5 * max(1.0, abs(value)), msg=ffname + " " + name)
        # The energy does not change the molecule, and asking again gives the same
        self.assertEqual(mine, mol.interaction_energy())

    def test_interaction_energy_needs_two_compounds(self):
        """!
        @brief A molecule of one compound has no interaction energy.
        """
        with self.assertRaises(act.ActInputError) as caught:
            act.Molecule.from_smiles("CC(=O)C", self.ff).interaction_energy()
        self.assertIn("two compounds", str(caught.exception))

    def test_numbers_are_numpy_arrays(self):
        """!
        @brief Coordinates, forces, charges, frequencies and trajectories are
        numpy arrays of the shapes one expects, and copies.
        """
        import numpy  # pylint: disable=import-outside-toplevel
        mol = act.Molecule.from_smiles("CC(=O)C", self.ff)
        n = len(mol.atoms)
        for name, value, shape in (("coordinates", mol.coordinates, (n, 3)),
                                   ("forces", mol.forces(), (n, 3)),
                                   ("charges", mol.charges, (n,))):
            self.assertIsInstance(value, numpy.ndarray, name)
            self.assertEqual(shape, value.shape, name)
            self.assertEqual(numpy.float64, value.dtype, name)
        # A copy: changing it changes nothing
        coords = mol.coordinates
        coords += 1.0
        self.assertFalse(numpy.allclose(coords, mol.coordinates))
        res = mol.minimize(update=False)
        self.assertEqual((n, 3), res.coordinates.shape)
        modes = mol.normal_modes()
        self.assertEqual(modes.frequencies.shape, modes.intensities.shape)
        self.assertEqual((n, 3), modes.modes[0]["displacements"].shape)
        frequencies, intensities = modes.spectrum()
        self.assertEqual(frequencies.shape, intensities.shape)
        traj = mol.simulate(steps=20, energy_interval=10, frame_interval=10)
        self.assertEqual((2, n, 3), traj.frames.shape)
        self.assertEqual((2,), traj.energies["EPOT"].shape)
        self.assertEqual(traj.times.shape, traj.frame_steps.shape)
        # Numbers in, as arrays as well
        mol.set_charges(numpy.zeros(n))
        self.assertEqual(0.0, float(numpy.abs(mol.charges).max()))

    def test_messages_go_to_the_logger(self):
        """!
        @brief The warnings and the information of ACT go to the logger act,
        which shows nothing unless the program asks for it.
        """
        import logging  # pylint: disable=import-outside-toplevel
        with self.assertLogs("act", level="INFO") as logged:
            act.Molecule.from_smiles("CC(=O)C", self.ff)
        self.assertTrue(any("propan-2-one" in line for line in logged.output), logged.output)
        self.assertTrue(all(line.split(":")[0] in ("INFO", "WARNING") for line in logged.output), logged.output)
        # The package does not print them by itself
        self.assertTrue(any(isinstance(h, logging.NullHandler) for h in logging.getLogger("act").handlers))
        # An error is raised, not logged
        with self.assertRaises(act.ActInputError):
            act.Molecule.from_smiles("C(C", self.ff)

    def test_molecule_from_file(self):
        """!
        @brief A molecule is read from a structure file.
        """
        mol = act.Molecule.from_file(self.data / "molecules" / "acetone.sdf", self.ff)
        self.assertEqual("C3H6O", mol.formula)
        self.assertEqual(20, len(mol.atoms))
        self.assertLess(mol.energy()["EPOT"], 0)
        # Another force field, without shells: another molecule
        plain = act.Molecule.from_file(str(self.data / "molecules" / "acetone.sdf"),
                                       act.ForceField.read(self.data / "forcefields" / "ACS-g1.xml"))
        self.assertEqual(10, len(plain.atoms))
        self.assertEqual(20, len(mol.atoms))

    def test_molecule_keeps_its_force_field(self):
        """!
        @brief A molecule works after the last name for its force field is gone.
        """
        import gc  # pylint: disable=import-outside-toplevel
        mol = act.Molecule.from_smiles("O", act.ForceField.read(self.ffname))
        gc.collect()
        self.assertEqual("H2O", mol.formula)
        self.assertIn("EPOT", mol.energy())
        self.assertIsInstance(mol.forcefield, act.ForceField)
        self.assertEqual(str(self.ffname), mol.forcefield.filename)

    def test_molecules_that_cannot_be_made(self):
        """!
        @brief What is no molecule, or one that the force field has no types
        for, raises an error that tells why, and Python lives on.
        """
        with self.assertRaises(act.ActInputError) as caught:
            act.Molecule.from_smiles("this is no molecule", self.ff)
        self.assertIn("SMILES", str(caught.exception))
        with self.assertRaises(act.ActError) as caught:
            act.Molecule.from_file(self.dir / "no-such-molecule.sdf", self.ff)
        self.assertIn("no-such-molecule.sdf", str(caught.exception))
        broken = self.dir / "broken.sdf"
        broken.write_text("not a structure\n")
        with self.assertRaises(act.ActInputError) as caught:
            act.Molecule.from_file(broken, self.ff)
        self.assertIn("broken.sdf", str(caught.exception))
        # A force field for water has no types for the atoms of ethanol
        water = act.ForceField.read(self.data / "forcefields" / "TIP5P.xml")
        with self.assertRaises(act.ActInputError) as caught:
            act.Molecule.from_smiles("CCO", water)
        self.assertTrue(str(caught.exception))
        # Still alive
        self.assertEqual("H2O", act.Molecule.from_smiles("O", self.ff).formula)
        for wrong in (None, 3):
            with self.assertRaises(TypeError):
                act.Molecule.from_smiles(wrong, self.ff)
        with self.assertRaises(TypeError):
            act.Molecule.from_smiles("O", "ACS-pg1.xml")
        with self.assertRaises(TypeError):
            act.Molecule.from_smiles("O")  # pylint: disable=no-value-for-parameter


@unittest.skipUnless(act.classes_available(), "ACT was built without its Python classes")
@unittest.skipUnless(os.path.isdir(TESTDATA), "ACT_TESTDATA is not set")
class MolPropsTest(unittest.TestCase):
    """!
    @brief The compounds of a molprop file with their reference data.
    """

    def setUp(self):
        """!
        @brief Find the test data and make a directory for files.
        """
        self.mpdir = Path(TESTDATA) / "molprop"
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self._tmp.name)

    def tearDown(self):
        """!
        @brief Remove the files.
        """
        self._tmp.cleanup()

    def test_select_and_exclude(self):
        """!
        @brief New sets with some of the compounds, by name or place.
        """
        mps = act.read_molprop(self.mpdir / "molprop.xml")
        names = mps.names
        chosen = mps.select([names[2], 0])
        self.assertEqual([names[2], names[0]], chosen.names)
        self.assertEqual(plain(mps.calculations(2)), plain(chosen.calculations(0)))
        rest = mps.exclude([names[1]])
        self.assertEqual([names[0], names[2]], rest.names)
        # The original is as it was
        self.assertEqual(names, mps.names)
        with self.assertRaises(KeyError):
            mps.select(["no such compound"])
        with self.assertRaises(IndexError):
            mps.exclude([7])

    @unittest.skipUnless(os.environ.get("ACT_ALEXANDRIA"), "ACT_ALEXANDRIA is not set")
    def test_merge_same_as_the_program(self):
        """!
        @brief Two sets are merged as alexandria edit_mp merges files.
        """
        first = self.mpdir / "booze-mono.xml"
        second = self.mpdir / "water-butanol.xml"
        merged = act.read_molprop(first).merge(act.read_molprop(second))
        mine = self.dir / "merged.xml"
        merged.write(mine)
        theirs = self.dir / "edit_mp.xml"
        act.run("edit_mp", mp=[first, second], o=theirs, molinfo=False)
        self.assertEqual(theirs.read_text(), mine.read_text())
        # A compound that both have is there once, with the calculations of both
        both = act.read_molprop(first).merge(act.read_molprop(first))
        self.assertEqual(len(act.read_molprop(first)), len(both))
        self.assertEqual(2 * len(act.read_molprop(first).calculations(0)), len(both.calculations(0)))

    def test_read_a_molprop_file(self):
        """!
        @brief The compounds of a file, what each is and what data it has.
        """
        mps = act.read_molprop(self.mpdir / "molprop.xml")
        self.assertEqual(3, len(mps))
        self.assertEqual(["1-butanol", "acetone", "water"], mps.names)
        self.assertEqual(str(self.mpdir / "molprop.xml"), mps.filename)
        self.assertIn("3 compounds", repr(mps))
        summary = mps.summary
        self.assertEqual((3, 3, 0), (summary["molecules"], summary["monomers"],
                                     summary["dimers"]))
        self.assertGreater(summary["properties"]["DeltaHform"], 0)
        # A compound by its place or by its name
        water = mps[2]
        self.assertEqual(water, mps["water"])
        self.assertEqual(water, mps[-1])
        self.assertEqual(("water", "H2O", 3, 0), (water["name"], water["formula"],
                                                  water["natom"], water["charge"]))
        self.assertEqual(1, len(water["fragments"]))
        self.assertGreater(water["experiments"], 0)
        self.assertEqual([c["name"] for c in mps], mps.names)
        for wrong in (3, -4):
            with self.assertRaises(IndexError):
                mps[wrong]  # pylint: disable=pointless-statement
        with self.assertRaises(KeyError):
            mps["no such compound"]  # pylint: disable=pointless-statement
        with self.assertRaises(TypeError):
            mps[1.5]  # pylint: disable=pointless-statement
        # Dimers
        dimers = act.read_molprop(str(self.mpdir / "mpdimer.xml"))
        self.assertEqual(3, dimers.summary["dimers"])
        self.assertEqual(2, len(dimers["1-propanol#ethane"]["fragments"]))

    def test_calculations_of_a_compound(self):
        """!
        @brief The calculations of a compound: how each was done, the
        structure and the values of the properties.
        """
        mps = act.read_molprop(self.mpdir / "molprop.xml")
        calculations = mps.calculations("water")
        self.assertEqual(mps["water"]["experiments"], len(calculations))
        self.assertEqual(plain(calculations), plain(mps.calculations(2)))
        first = calculations[0]
        for key in ("reference", "conformation", "method", "basisset", "program", "jobtype",
                    "datafile", "atoms", "coordinates", "properties"):
            self.assertIn(key, first)
        # The atoms of water with a position in nm for each
        self.assertEqual(["H", "H", "O"], sorted(atom["name"][0] for atom in first["atoms"]))
        self.assertEqual(3, len(first["coordinates"]))
        distances = [sum((first["coordinates"][0][d] - first["coordinates"][k][d]) ** 2
                         for d in range(3)) ** 0.5 for k in (1, 2)]
        self.assertTrue(all(0.05 < distance < 0.2 for distance in distances), distances)
        # Each property with what it is, its value and its unit
        properties = [p for c in calculations for p in c["properties"]]
        self.assertGreater(len(properties), 0)
        for item in properties:
            self.assertEqual({"observable", "type", "value", "error", "unit", "temperature",
                              "vector"}, set(item))
            # A number, or the components of a vector or tensor; None for
            # what a property is not, such as an electrostatic potential
            self.assertTrue(item["value"] is None or isinstance(item["value"], float))
            self.assertTrue(item["vector"] is None or isinstance(item["vector"], list))
        enthalpy = [item for item in properties if item["observable"] == "DeltaHform"][0]
        self.assertIsInstance(enthalpy["value"], float)
        self.assertTrue(enthalpy["unit"])
        dipoles = [item for item in properties if item["observable"] == "dipole"]
        if dipoles:
            self.assertEqual(3, len(dipoles[0]["vector"]))
        self.assertIn("DeltaHform", {item["observable"] for item in properties})
        # The number of values is that of the description
        counted = {}
        for item in properties:
            counted[item["observable"]] = counted.get(item["observable"], 0) + 1
        self.assertEqual(mps["water"]["properties"], counted)
        # Interaction energies of a dimer, one for each structure
        dimers = act.read_molprop(self.mpdir / "mpdimer.xml")
        energies = [p["value"] for c in dimers.calculations(0) for p in c["properties"]
                    if p["observable"] == "InteractionEnergy"]
        self.assertEqual(dimers[0]["properties"]["InteractionEnergy"], len(energies))
        self.assertGreater(len(set(energies)), 1)
        with self.assertRaises(KeyError):
            mps.calculations("no such compound")
        with self.assertRaises(IndexError):
            mps.calculations(17)

    def test_validate(self):
        """!
        @brief The check of edit_mp -validate, with what it found.
        """
        mps = act.read_molprop(self.mpdir / "molprop.xml")
        checked = mps.validate()
        self.assertEqual((True, True, 0, []), (checked["checked"], checked["ok"],
                                               checked["problems"], checked["invalid"]))
        self.assertIsInstance(checked["messages"], list)
        # Dimers, of which the compounds are compared with each other
        checked = act.read_molprop(self.mpdir / "booze-inter.xml").validate()
        self.assertEqual({"checked", "ok", "problems", "invalid", "duplicates", "with_duplicates",
                          "messages"}, set(checked))
        self.assertEqual(0, checked["duplicates"])
        self.assertTrue(checked["ok"])
        # A file with a problem is in the tests of edit_mp, which uses the same check

    def test_write_and_errors(self):
        """!
        @brief The compounds are written to a file that reads back the same,
        and a file that cannot be read or written raises an error.
        """
        mps = act.read_molprop(self.mpdir / "molprop.xml")
        out = self.dir / "written.xml"
        mps.write(out)
        again = act.read_molprop(out)
        self.assertEqual(mps.names, again.names)
        self.assertEqual(mps.summary, again.summary)
        self.assertEqual(plain(mps.calculations("water")), plain(again.calculations("water")))
        packed = self.dir / "packed.xml"
        mps.write(packed, compress=True)
        self.assertLess(packed.stat().st_size, out.stat().st_size)
        self.assertEqual(mps.names, act.read_molprop(packed).names)
        with self.assertRaises(act.ActFileError) as caught:
            act.read_molprop(self.dir / "no-such-file.xml")
        self.assertIn("no-such-file.xml", str(caught.exception))
        notxml = self.dir / "notxml.xml"
        notxml.write_text("water|Train\n")
        with self.assertRaises(act.ActFileError):
            act.read_molprop(notxml)
        with self.assertRaises(act.ActFileError):
            mps.write(self.dir / "no" / "such" / "directory" / "mp.xml")
        for wrong in (None, 3):
            with self.assertRaises(TypeError):
                act.read_molprop(wrong)
        # Still alive
        self.assertEqual(3, len(act.read_molprop(self.mpdir / "molprop.xml")))


if __name__ == "__main__":
    unittest.main()
