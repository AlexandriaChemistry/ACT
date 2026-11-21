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
@brief Test of the conda package alexandria-python in the environment it
was installed in, see recipe.yaml: the package act is that of the
environment, it has its classes, and it finds the alexandria program and the
data files of the environment. The test files are looked for below the
current directory.
"""
import math
import os
import sys

import act


def find(name: str) -> str:
    """!
    @brief Find a test file below the current directory.
    @param name The name of the file
    @return The path of the file
    """
    for top, _, files in os.walk("."):
        if name in files:
            return os.path.join(top, name)
    raise SystemExit("There is no test file %s below %s" % (name, os.getcwd()))


def inside(path: str) -> bool:
    """!
    @brief Tell whether a file or directory is in the environment of this Python.
    @param path The file or directory
    @return Whether it is
    """
    prefix = os.path.realpath(sys.prefix)
    return os.path.commonpath([prefix, os.path.realpath(path)]) == prefix


def main() -> None:
    """!
    @brief Check the installed package.
    """
    assert inside(act.__file__), act.__file__
    assert act.classes_available(), "The package act has no classes"
    # The data files and the program of the package alexandria
    data = act.data_directory()
    assert data and inside(data), data
    assert os.path.isfile(os.path.join(data, "atomprops.csv")), data
    words = act.config.executable_words()
    assert os.path.basename(words[0]) == "alexandria" and inside(words[0]), words
    assert "train_ff" in act.commands.available()
    # The classes compute what the program computes
    ffname = find("ACS-g1.xml")
    forcefield = act.ForceField.read(ffname)
    assert not forcefield.polarizable
    mol = act.Molecule.from_file(find("acetone.sdf"), forcefield)
    assert mol.formula == "C3H6O", mol.formula
    energy = mol.energy()["EPOT"]
    assert math.isfinite(energy), energy
    with act.commands.minimize(ff=ffname, f=find("acetone.sdf")) as told:
        theirs = told.data["Energies before"]["EPOT"]["value"]
    assert abs(theirs - energy) <= 1e-3 * abs(energy), (theirs, energy)
    minimum = mol.minimize()
    assert minimum.converged, minimum.status
    assert minimum.energy_after["EPOT"] <= energy
    print("act", act.__file__, "with", words[0], "and", data)


if __name__ == "__main__":
    main()
