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
import sys
import tempfile
import unittest
from pathlib import Path

# The package of the build tree has the module with the classes; the tests
# of a build say where it is. Otherwise the package of the source code.
sys.path.insert(0, os.environ.get("ACT_PYTHON_PATH")
                or os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import act  # noqa: E402  pylint: disable=wrong-import-position

## The test data of ACT
TESTDATA = os.environ.get("ACT_TESTDATA", "")


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
        self.assertEqual([a["charge"] for a in atoms], charges)
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
        self.assertEqual(res.coordinates, mol.coordinates)
        self.assertNotEqual(start, mol.coordinates)
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
        self.assertEqual(start, mol.coordinates)
        self.assertNotEqual(start, res.coordinates)
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
        for name, value in data["Energies after"].items():
            self.assertAlmostEqual(value["value"], res.energy_after[name],
                                   delta=1e-3 + 1e-5 * abs(value["value"]), msg=name)

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
        self.assertEqual(start, mol.coordinates)
        with self.assertRaises(TypeError):
            mol.minimize("LBFGS")  # pylint: disable=too-many-function-args
        with self.assertRaises(TypeError):
            mol.minimize(max_iter="many")

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
        self.assertEqual(sorted(modes.frequencies), modes.frequencies)
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
        self.assertNotEqual(start, mol.coordinates)
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
        self.assertEqual(list(range(frequencies[0], frequencies[-1] + 1)), frequencies)
        # From below the lowest to above the highest frequency
        self.assertLess(frequencies[0], modes.frequencies[0] - 100)
        self.assertGreater(frequencies[-1], modes.frequencies[-1] + 100)
        self.assertAlmostEqual(1.0, sum(intensities), places=6)
        # A maximum at the strongest vibration
        strongest = modes.frequencies[modes.intensities.index(max(modes.intensities))]
        top = frequencies[intensities.index(max(intensities))]
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
        self.assertEqual([int(row[0]) for row in spectrum], frequencies)
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


if __name__ == "__main__":
    unittest.main()
