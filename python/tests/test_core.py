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


if __name__ == "__main__":
    unittest.main()
