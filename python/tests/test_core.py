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
