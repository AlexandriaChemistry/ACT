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
@brief Tests of the reference of the package act in the manual, which
docs/manual/act_api.py makes from the documentation in the code: every
public name is described, and every description has a text.
"""
import importlib.util
import os
import sys
import unittest
from typing import Any, Dict, List

## The directory of the package act of the source tree, or of a build
PACKAGE_DIR = os.environ.get("ACT_PYTHON_PATH") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

## The generator of the reference in the manual
GENERATOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "docs", "manual", "act_api.py")


def load_generator() -> Any:
    """!
    @brief Import docs/manual/act_api.py.
    @return The module.
    """
    spec = importlib.util.spec_from_file_location("act_api", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def descriptions(lines: List[str]) -> Dict[str, List[str]]:
    """!
    @brief The descriptions in a reference, by the name they describe.
    @param lines The reference in reStructuredText.
    @return For each name, with the class for a member, the lines of its text.
    """
    found: Dict[str, List[str]] = {}
    current = None
    cls = None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith(".. py:") and "::" in stripped:
            kind, _, signature = stripped[len(".. py:"):].partition("::")
            name = signature.strip().split("(")[0].strip()
            if kind in ("class", "exception"):
                cls = name
            elif kind in ("function", "data"):
                cls = None
            current = name if kind in ("class", "exception", "function", "data", "module") or cls is None \
                else cls + "." + name
            found[current] = []
        elif current and stripped and not stripped.startswith(":"):
            found[current].append(stripped)
    return found


class ApiReferenceTest(unittest.TestCase):
    """!
    @brief The reference of the package in the manual.
    """

    def setUp(self) -> None:
        """!
        @brief Make the reference of the package of the tests.
        """
        os.environ["ACT_PYTHON_PATH"] = PACKAGE_DIR
        self.generator = load_generator()
        self.act = self.generator.load_package()
        self.lines = self.generator.api_lines(self.act)
        self.found = descriptions(self.lines)

    def tearDown(self) -> None:
        """!
        @brief Forget the package that the generator imported.
        """
        for name in [m for m in sys.modules if m == "act" or m.startswith("act.")]:
            del sys.modules[name]

    def test_every_public_name_is_described(self) -> None:
        """!
        @brief The functions, classes and errors of act are in the reference.
        """
        names = list(self.act.__all__)
        if self.act.classes_available():
            names += list(self.act._CLASSES)  # pylint: disable=protected-access
        for name in names:
            self.assertIn(name, self.found, name)

    def test_every_description_has_a_text(self) -> None:
        """!
        @brief Nothing in the reference lacks its documentation.
        """
        for name, text in self.found.items():
            if name == "act":
                continue
            self.assertTrue(text, "%s has no documentation" % name)

    def test_without_classes_the_reference_says_so(self) -> None:
        """!
        @brief A package without its classes tells where they come from.
        """
        text = "\n".join(self.lines)
        self.assertEqual(not self.act.classes_available(), "built without the classes of ACT" in text)

    def test_doxygen_is_converted(self) -> None:
        """!
        @brief The fields and code of the docstrings become those of reStructuredText.
        """
        lines = self.generator.doxygen_to_rst("""!
        @brief Do something with *stars* and a_name.
        @code{.py}
        x = f(a_b)
        @endcode
        @param a   The first.
        @param b   The second, which
        goes on.
        @return    The result.
        @throws ValueError if wrong.
        """)
        text = "\n".join(lines)
        self.assertIn("Do something with \\*stars\\* and a\\_name.", text)
        self.assertIn(".. code-block:: python", text)
        self.assertIn("   x = f(a_b)", text)
        self.assertIn(":param a: The first.", text)
        self.assertIn(":param b: The second, which\n   goes on.", text)
        self.assertIn(":returns: The result.", text)
        self.assertIn(":raises ValueError: if wrong.", text)
        # The fields start after an empty line
        index = lines.index(":param a: The first.")
        self.assertEqual("", lines[index - 1])


if __name__ == "__main__":
    unittest.main()
